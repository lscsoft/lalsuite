/*
 *  Copyright (C) 2013 Chris Pankow
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <lal/LALDetCharHveto.h>

static gint compare(gconstpointer a, gconstpointer b) {
        const SnglBurst *_a = a;
        const SnglBurst *_b = b;

        return XLALCompareSnglBurstByPeakTimeAndSNR(&_a, &_b);
}

// Program level functions, mostly utility

// Fake a series of triggers
void populate_trig_sequence( GSequence* trig_sequence );
// Read a series of triggers from XML
void populate_trig_sequence_from_file( GSequence* trig_sequence, const char* fname, double min_snr, GSequence* ignore_list );
// print a hash table of string / int entries
void print_hash_table( GHashTable* tbl );
// Read a sequence of channels to ignore
void get_ignore_list( const char* fname, GSequence* ignorelist );
// retrieve live time
void calculate_livetime( const char* fname, LALSegList* live_segs );

extern int lalDebugLevel;

//int main(int argc, char** argv){
int main(int argc, char** argv){

	GSequence* ignorelist = g_sequence_new(free);
	if( argc > 4 ){
		get_ignore_list( argv[4], ignorelist );
	}

	//lalDebugLevel=7;

	double sig_thresh = 3.0;
	// TODO: Significance or SNR option
	double min_de_snr = 30, min_aux_snr = 100;

	fprintf( stderr, "Name: %s\n", argv[1] );

	GSequence* trig_sequence = g_sequence_new((GDestroyNotify)XLALDestroySnglBurst);

	char refchan[1024];
	if( argc > 2 ){
		strcpy(refchan, "LSC-DARM_ERR");
		printf( "Reading file: %s\n", argv[2] );
		populate_trig_sequence_from_file( trig_sequence, argv[2], min_de_snr, ignorelist );
		fprintf( stderr, "list Length: %d\n", g_sequence_get_length( trig_sequence ) );
		printf( "Reading file: %s\n", argv[1] );
		populate_trig_sequence_from_file( trig_sequence, argv[1], min_aux_snr, ignorelist );
	} else {
		strcpy(refchan, "CHAN1");
		populate_trig_sequence( trig_sequence );
	}
	fprintf( stderr, "list Length: %d\n", g_sequence_get_length( trig_sequence ) );
	fprintf( stderr, "Reference channel: %s\n", refchan );

	GSequenceIter* igitr = g_sequence_get_begin_iter(ignorelist);
	while( !g_sequence_iter_is_end(igitr) ){
		printf( "Ignoring %s\n", (char *)g_sequence_get(igitr) );
		igitr = g_sequence_iter_next(igitr);
	}
	g_sequence_free(ignorelist);

	double livetime = 100;
	if( argc > 3 ){
		const char* livefname = argv[3];
		printf( "Livetime filename: %s\n", livefname );
		LALSegList live;
		XLALSegListInit( &live );
		calculate_livetime( argv[3], &live );
		size_t i;
		for( i=0; i<live.length; i++ ){
			LALSeg s = live.segs[i];
			livetime += XLALGPSDiff( &s.end, &s.start );
		}
		printf( "len: %d\n", g_sequence_get_length( trig_sequence ) );
		XLALDetCharPruneTrigs( trig_sequence, &live );
		printf( "len: %d\n", g_sequence_get_length( trig_sequence ) );
		XLALSegListClear( &live );
	}
	printf( "Livetime: %f\n", livetime );

	LALSegList vetoes;
	XLALSegListInit( &vetoes );

	GHashTable *chancount, *chanhist;
	//chanhist = g_hash_table_new_full( g_str_hash, g_str_equal, &g_free, &g_free );
	chanhist = g_hash_table_new( g_str_hash, g_str_equal );
	chancount = g_hash_table_new( g_str_hash, g_str_equal );

	int rnd = 1;
	double rnd_sig = 0;

	double wind = 0.8;
	double t_ratio = wind / livetime;
	int coinc_type = 1; // unique coincidences

	GList *channames = g_hash_table_get_keys( chancount );
	size_t nchans = g_list_length( channames ) - 1;
	char winner[1024];
	do {

		// Create our channel coincidence histogram
		// FIXME: Free memory?
		g_hash_table_remove_all( chanhist );
		// Create our channel count histogram
		// FIXME: Free memory?
		g_hash_table_remove_all( chancount );

		// TODO: Is there a way to avoid a rescan every round?
		XLALDetCharScanTrigs( chancount, chanhist, trig_sequence, refchan, wind, coinc_type );
		printf( "Trigger count:\n" );
		print_hash_table( chancount );
		printf( "Trigger coincidences with %s, window (%f):\n", refchan, wind );
		print_hash_table( chanhist );
		// TODO: Heuristic: check the highest # of coincidences first
		// (maybe normalize by number of triggers

		strcpy( winner, "" );
		// If there's no value for the target channel, there's no vetoes to do,
		// move on.
		printf( "Round %d\n", rnd );
		if( !g_hash_table_lookup(chancount, refchan) ){
			fprintf( stderr, "No triggers in target channel.\n" );
			break;
		} else {
			rnd_sig = XLALDetCharVetoRound( winner, chancount, chanhist, refchan, t_ratio );
		}
		printf( "Done, round %d, winner (%s) sig %g\n", rnd, winner, rnd_sig );
		//exit(0);

		if( rnd_sig > sig_thresh ){
			// TODO: Subtract vetoed livetime.
			//XLALSegListAppend( vetoes, wind );
			LALSeg veto;
			LIGOTimeGPS start;
			LIGOTimeGPS stop;
			XLALGPSSetREAL8( &start, -wind/2.0 );
			XLALGPSSetREAL8( &stop, wind/2.0 );
			XLALSegSet( &veto, &start, &stop, 0 );
			XLALDetCharRemoveTrigs( trig_sequence, veto, winner );

			// Remove the channel from consideration
			printf( "Removing %s from count\n", winner );

			// TODO: Reenable
			/*
			char* key;
			size_t *value;
			gpointer *cname = (gpointer*)&key, *n = (gpointer*)&value;
			g_hash_table_lookup_extended( chanhist, winner, cname, n );
			g_free(key);
			g_free(value);
			g_hash_table_remove( chanhist, winner );
			*/
		}
		rnd++;
		//if( rnd > 2 ) break;
	} while( rnd_sig > sig_thresh && (size_t)rnd < nchans );
	printf( "Last round did not pass significance threshold or all channels have been vetoed. Ending run.\n" );

	// TODO: Delete the key / value pairs still remaining
	g_hash_table_destroy( chanhist );
	g_hash_table_destroy( chancount );

	g_sequence_free( trig_sequence );

	XLALSegListClear( &vetoes );
	//LALCheckMemoryLeaks();
	return 0;
}

void populate_trig_sequence_from_file( GSequence* trig_sequence, const char* fname, double min_snr, GSequence* ignore_list ){
	SnglBurst* tbl = XLALSnglBurstTableFromLIGOLw( fname );
	if( !tbl ) return;

	//#pragma omp parallel
	//{
	//#pragma omp single
	//{
	//for( ; tbl; tbl=tbl->next ){
	do {
		//#pragma omp task
		//{
		gboolean ignore = FALSE;
		GSequenceIter* igitr;
		// FIXME: Support ignorelist == NULL
		igitr = g_sequence_get_begin_iter(ignore_list);
		while( !g_sequence_iter_is_end(igitr) ){
			/*
			 * Note this will cause incorrect behavior if the same channel name
			 * with different interferometers is included in the same run as
			 * the SB channel names do not include ifo by default.
			 */
			ignore = (strstr(g_sequence_get(igitr), tbl->channel) != NULL);
			if( !ignore ) {
				igitr = g_sequence_iter_next(igitr);
			} else {
				break;
			}
		}
		if( tbl->confidence > min_snr && !ignore ){
			g_sequence_insert_sorted( trig_sequence, tbl, (GCompareDataFunc)compare, NULL );
			tbl=tbl->next;
		} else {
			// Thought: move this node to a new list which you can delete later.
			// Question: Is that any faster since the delete is probably a big
			// portion of the CPU time
			SnglBurst* tmp = tbl;
			tbl=tbl->next;
			XLALDestroySnglBurst(tmp);
		}
		//} // end pragma task
	//}
	} while( tbl );
	//} // end pragma single
	//} // end pragma parallel
}

void populate_trig_sequence( GSequence* trig_sequence ){

	char const *channel_list[] = {
		"CHAN1",
		"CHAN2",
		"CHAN3",
		"CHAN4",
		"CHAN5",
		"CHAN6"
	};

	SnglBurst* t = NULL;
	int i = 0;
	//for(; i<1000; i++){
	for(; i<1000; i++){
		t = XLALCreateSnglBurst();
		XLALGPSSet(&t->peak_time, rand() % 100, rand() % XLAL_BILLION_INT8);
		t->snr = 1;
		t->central_freq = 100;
		t->event_id = i;
		strcpy( t->channel, channel_list[rand()%6] );
		strcpy( t->ifo, "H1" );
		g_sequence_insert_sorted( trig_sequence, t, (GCompareDataFunc)compare, NULL );
		fprintf( stderr, "%d %s %d.%d\n", i, t->channel, t->peak_time.gpsSeconds, t->peak_time.gpsNanoSeconds );
	}
	fprintf( stderr, "\nCreated %d events\n", i );
}

void print_hash_table( GHashTable* tbl ){
    GHashTableIter iter;
    gpointer key, val;
	g_hash_table_iter_init( &iter, tbl );
	while( g_hash_table_iter_next( &iter, &key, &val ) ){
		printf( "%s: %lu\n", (char *)key, *(size_t*)val );
	}
}

void get_ignore_list( const char* fname, GSequence* ignorel ){

	int cnt = 0;
	char* tmp;
	tmp = malloc( sizeof(char)*512 );
	FILE* lfile = fopen( fname, "r" );
	while(!feof(lfile)){
		cnt = fscanf( lfile, "%s", tmp );
		if( cnt == EOF ) break;
		g_sequence_append( ignorel, g_strdup(tmp) );
	}
	fclose(lfile);
}

void calculate_livetime( const char* fname, LALSegList* live_segs ){
	FILE* lfile = fopen( fname, "r" );
	//float st, end;
	int st, end;
	LALSeg *seg;
	int cnt = 0;
	while(!feof(lfile)){
		// FIXME: Having floats here causes a 32 second offset for both values
		//cnt = fscanf( lfile, "%f %f", &st, &end );
		cnt = fscanf( lfile, "%d %d", &st, &end );
		if( cnt == EOF ) break;
		LIGOTimeGPS start;
		LIGOTimeGPS stop;
		XLALGPSSetREAL8( &start, (double)st );
		XLALGPSSetREAL8( &stop, (double)end );
		seg = XLALSegCreate( &start, &stop, 0 );
		XLALSegListAppend( live_segs, seg );
		XLALFree( seg );
	}
	fclose(lfile);
}
