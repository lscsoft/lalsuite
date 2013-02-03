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
#include <lal/LIGOLwXML.h>

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
void print_hash_table( GHashTable* tbl, const char* file );
// print a hash table of string / float entries
void print_hash_table_float( GHashTable* tbl, const char* file );
// Read a sequence of channels to ignore
void get_ignore_list( const char* fname, GSequence* ignorelist );
// retrieve live time
void calculate_livetime( const char* fname, LALSegList* live_segs );
// Write trigger set
void write_triggers( GSequence* trig_sequence, const char* fname );

// Encode / Decode a subround winner
void encode_rnd_str(char* winner, double wind, double thresh);
void decode_rnd_str(char* winner_str, char* chan, double* wind, double* thresh);

extern int lalDebugLevel;

//int main(int argc, char** argv){
int main(int argc, char** argv){

	GSequence* ignorelist = g_sequence_new(free);
	if( argc > 4 ){
		get_ignore_list( argv[4], ignorelist );
	}

	//lalDebugLevel = 1 | 2 | 4 | 32;
	lalDebugLevel = 0;

	double sig_thresh = 3.0;
	// TODO: Significance or SNR option
	// Example 1
	double min_de_snr = 30, min_aux_snr = 50;
	int nthresh = 8;
	double aux_snr_thresh[8] = {3200.0, 1600.0, 800.0, 600.0, 400.0, 200.0, 100.0, 50.0};
	int nwinds = 5;
	double twins[5] = {0.1, 0.2, 0.4, 0.8, 1.0};

	// Example 2
	//double min_de_snr = 8, min_aux_snr = 10;
	/*
	int nthresh = 7, nwinds = 5;
	double aux_snr_thresh[7] = {300.0, 100.0, 40.0, 20.0, 15.0, 12.0, 10.0};
	double twins[5] = {0.1, 0.2, 0.4, 0.8, 1.0};
	//double twins[1] = {0.4};
	*/

	printf( "Name: %s\n", argv[1] );

	GSequence* trig_sequence = g_sequence_new((GDestroyNotify)XLALDestroySnglBurst);

	char refchan[1024];
	if( argc > 2 ){
		strcpy(refchan, "LSC-DARM_ERR");
		printf( "Reading file: %s\n", argv[2] );
		populate_trig_sequence_from_file( trig_sequence, argv[2], min_de_snr, ignorelist );
		printf( "list Length: %d\n", g_sequence_get_length( trig_sequence ) );
		printf( "Reading file: %s\n", argv[1] );
		populate_trig_sequence_from_file( trig_sequence, argv[1], min_aux_snr, ignorelist );
	} else {
		strcpy(refchan, "CHAN1");
		populate_trig_sequence( trig_sequence );
	}
	printf( "list Length: %d\n", g_sequence_get_length( trig_sequence ) );
	printf( "Reference channel: %s\n", refchan );

	GSequenceIter* igitr = g_sequence_get_begin_iter(ignorelist);
	while( !g_sequence_iter_is_end(igitr) ){
		printf( "Ignoring %s\n", (char *)g_sequence_get(igitr) );
		igitr = g_sequence_iter_next(igitr);
	}
	g_sequence_free(ignorelist);

	LALSegList live;
	XLALSegListInit( &live );
	double livetime = 100;
	if( argc > 3 ){
		const char* livefname = argv[3];
		printf( "Livetime filename: %s\n", livefname );
		calculate_livetime( argv[3], &live );
		size_t i;
		for( i=0; i<live.length; i++ ){
			LALSeg s = live.segs[i];
			livetime += XLALGPSDiff( &s.end, &s.start );
		}
		XLALSegListSort( &live );
		printf( "len: %d\n", g_sequence_get_length( trig_sequence ) );
		XLALDetCharPruneTrigs( trig_sequence, &live, min_de_snr, NULL );
		printf( "len: %d\n", g_sequence_get_length( trig_sequence ) );
	}
	printf( "Livetime: %f\n", livetime );

	LALSegList vetoes;
	XLALSegListInit( &vetoes );

	GHashTable *chancount, *chanhist, *subround_winners;
	//chanhist = g_hash_table_new_full( g_str_hash, g_str_equal, &g_free, &g_free );
	chanhist = g_hash_table_new( g_str_hash, g_str_equal );
	chancount = g_hash_table_new( g_str_hash, g_str_equal );
	subround_winners = g_hash_table_new( g_str_hash, g_str_equal );

	int rnd = 1;

	double wind = 0.8;
	int coinc_type = 1; // unique coincidences
	double rnd_sig = 0;

	char outpath[512];

	GList *channames = g_hash_table_get_keys( chancount );
	size_t nchans = g_list_length( channames ) - 1;

	/*
	 * Veto round loop.
	 */
	do {
		printf( "Round %d\n", rnd );
		
		// Set up the output path for this round
		sprintf( outpath, "mkdir -p round_%d", rnd );
		system( outpath );

		// Clear the winners of the previous subrounds
		// FIXME: Free memory?
		g_hash_table_remove_all( subround_winners );
		nwinds--;
		nthresh--;

		/*
		 * FIXME: We can do this more efficiently by passing this to the scan
		 * function, and only doing a single pass.
		 */
		char *winner;
		double *sig = malloc(sizeof(double));
		int nw = nwinds, nt = nthresh;
		for( nt=nthresh; nt >= 0; nt-- ){
			min_aux_snr = aux_snr_thresh[nt];
			printf( "len: %d\n", g_sequence_get_length( trig_sequence ) );
			XLALDetCharPruneTrigs( trig_sequence, &live, min_aux_snr, refchan );
			printf( "len: %d\n", g_sequence_get_length( trig_sequence ) );
			for( nw=nwinds; nw >= 0; nw-- ){
				printf( "Window: %d\n", nw );
				wind = twins[nw];
				printf( "wind: %f\n", wind );
				winner = malloc( sizeof(char)*1024 );
				printf( "Prev Winner str: %s\n", winner );

				// Create our channel coincidence histogram
				// FIXME: Free memory?
				printf( "Clearing coincidence table\n" );
				g_hash_table_remove_all( chanhist );
				// Create our channel count histogram
				// FIXME: Free memory?
				printf( "Clearing count table\n" );
				g_hash_table_remove_all( chancount );

				// TODO: Is there a way to avoid a rescan every round?
				XLALDetCharScanTrigs( chancount, chanhist, trig_sequence, refchan, wind, coinc_type );
				printf( "Trigger count:\n" );
				print_hash_table( chancount, NULL );
				sprintf( outpath, "round_%d_count.txt", rnd );
				print_hash_table( chancount, outpath );
				printf( "Trigger coincidences with %s, window (%g):\n", refchan, wind );
				print_hash_table( chanhist, NULL );
				sprintf( outpath, "round_%d_coinc.txt", rnd );
				print_hash_table( chanhist, outpath );

				strcpy( winner, "" );
				/* 
				 * If there's no value for the target channel, there's no 
				 * vetoes to do, move on.
				 */
				double t_ratio = wind / livetime;
				if( !g_hash_table_lookup(chancount, refchan) ){
					fprintf( stderr, "No triggers in target channel.\n" );
					break;
				} else {
					rnd_sig = XLALDetCharVetoRound( winner, chancount, chanhist, refchan, t_ratio );
				}
				printf( "Sub-round winner, window %2.2f, thresh %f: %s sig %g\n",  wind, min_aux_snr, winner, rnd_sig );

				sig = malloc(sizeof(double));
				*sig = rnd_sig;
				printf( "Sub-round winner, window %2.2f, thresh %f: %s sig %g\n",  wind, min_aux_snr, winner, *sig );
				encode_rnd_str( winner, wind, min_aux_snr), 
				g_hash_table_insert( subround_winners, winner, sig );

			}
		}
		print_hash_table_float( subround_winners, NULL );

		// Get the winner over all the subrounds
		GHashTableIter subrnd;
		g_hash_table_iter_init( &subrnd, subround_winners );
		char *chanwin;
		gpointer val, key;

		rnd_sig = 0;
		while( g_hash_table_iter_next( &subrnd, &key, &val ) ){
			if( *(double*)val > rnd_sig ){
				rnd_sig = *(double*)val;
				chanwin = (char*)key;
				printf( "rnd_sig: %g\n", rnd_sig );
				printf( "New winner: %s\n", chanwin );
			}
		}

		decode_rnd_str( chanwin, winner, &wind, &min_aux_snr );
		printf( "Done, round %d, winner: %s\n\tsignificance %g\n\twindow: %lf\n\tsnr thresh: %lf\n", rnd, winner, rnd_sig, wind, min_aux_snr );

		if( rnd_sig > sig_thresh ){
			// TODO: Subtract vetoed livetime.
			//XLALSegListAppend( vetoes, wind );
			LALSeg veto;
			LIGOTimeGPS start;
			LIGOTimeGPS stop;
			XLALGPSSetREAL8( &start, -wind/2.0 );
			XLALGPSSetREAL8( &stop, wind/2.0 );
			XLALSegSet( &veto, &start, &stop, 0 );
			// Remove the triggers veoted from the main list
			GSequence *vetoed_trigs = XLALDetCharRemoveTrigs( trig_sequence, veto, winner );
			sprintf( outpath, "./round_%d_vetoed_triggers.xml", rnd );
			// Write them for later use
			if( g_sequence_get_length(vetoed_trigs) > 0 ){
				write_triggers( vetoed_trigs, outpath );
			}

			// Remove livetime
			// livetime -= XLALDetCharHvetoExciseSegment( live, vetosegs );
			// TODO: Do this right
			livetime -= wind*g_sequence_get_length( vetoed_trigs );

			// Begone witcha!
			g_sequence_free( vetoed_trigs );

			/*
			 * TODO: This as supposed to remove the channel from the from list
			 * but we rescan anyway, so all of this is redundant -- regardless,
			 * this memory should be freed.
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
		//if( rnd > 1 ) break;
	} while( rnd_sig > sig_thresh && (size_t)rnd < nchans );
	printf( "Last round did not pass significance threshold or all channels have been vetoed. Ending run.\n" );

	// TODO: Delete the key / value pairs still remaining
	g_hash_table_destroy( chanhist );
	g_hash_table_destroy( chancount );

	g_sequence_free( trig_sequence );

	XLALSegListClear( &live );
	XLALSegListClear( &vetoes );
	//LALCheckMemoryLeaks();
	return 0;
}

void populate_trig_sequence_from_file( GSequence* trig_sequence, const char* fname, double min_snr, GSequence* ignore_list ){
	SnglBurst* tbl = XLALSnglBurstTableFromLIGOLw( fname );
	SnglBurst *begin = NULL, *deleteme = NULL;
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
		igitr = ignore_list ? g_sequence_get_begin_iter(ignore_list) : NULL;
		while( (igitr != NULL) & !g_sequence_iter_is_end(igitr) ){
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
			if( !deleteme ){
				begin = deleteme = tbl;
			} else {
				deleteme->next = tbl;
				deleteme = deleteme->next;
			}
			tbl=tbl->next;
		}
		//} // end pragma task
	//}
	} while( tbl );
	//} // end pragma single
	//} // end pragma parallel
	printf( "Deleting %d unused events.\n", XLALSnglBurstTableLength(begin) );
	deleteme = NULL;
	XLALDestroySnglBurstTable(begin);
	printf( "Done.\n" );
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

void print_hash_table( GHashTable* tbl, const char* file ){
    GHashTableIter iter;
    gpointer key, val;
	g_hash_table_iter_init( &iter, tbl );
	FILE* outfile = NULL;
	if( file ){
		outfile = fopen( file, "w" );
	}
	while( g_hash_table_iter_next( &iter, &key, &val ) ){
		if( file ){
			fprintf( outfile, "%s: %lu\n", (char *)key, *(size_t*)val );
		} else {
			printf( "%s: %lu\n", (char *)key, *(size_t*)val );
		}
	}
	if( outfile ){
		fclose( outfile );
	}
}
void print_hash_table_float( GHashTable* tbl, const char* file ){
    GHashTableIter iter;
    gpointer key, val;
	g_hash_table_iter_init( &iter, tbl );
	FILE* outfile = NULL;
	if( file ){
		outfile = fopen( file, "w" );
	}
	while( g_hash_table_iter_next( &iter, &key, &val ) ){
		if( file ){
			fprintf( outfile, "%s: %f\n", (char *)key, *(double*)val );
		} else {
			printf( "%s: %f\n", (char *)key, *(double*)val );
		}
	}
	if( outfile ){
		fclose( outfile );
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

void write_triggers( GSequence* trig_sequence, const char* fname ){

	// Convert a sequence back into SnglBurst
	GSequenceIter* sbit;
	SnglBurst *begin, *sb, *tmp;
	begin = sb = NULL;
	sbit = g_sequence_get_begin_iter( trig_sequence );
	while( !g_sequence_iter_is_end( sbit ) ){
		tmp = (SnglBurst*) g_sequence_get( sbit );
		sbit = g_sequence_iter_next( sbit );
		if( sb ){
			sb->next = tmp;
			sb = sb->next;
		} else {
			begin = sb = tmp;
		}
	}
	sb->next = NULL;
	
	LIGOLwXMLStream *str;
	str = XLALOpenLIGOLwXMLFile( fname );
	XLALWriteLIGOLwXMLSnglBurstTable( str, begin );
	XLALCloseLIGOLwXMLFile( str );
}

void encode_rnd_str(char* winner, double wind, double thresh){
	//winner = realloc( winner, sizeof(winner)+20*sizeof(char));
	sprintf( winner, "%s %1.5f %4.2f", winner, wind, thresh );
}
void decode_rnd_str(char* winner_str, char* chan, double* wind, double* thresh){
	sscanf( winner_str, "%s %lf %lf", chan, wind, thresh );
}
