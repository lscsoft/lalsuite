#include <lal/LALDetCharHveto.h>

size_t _count_subsequence( GSequence* subseq, SnglBurst* trig, float wind );
/*
 * This counts the number of coincidences between an input trigger; likely an
 * auxiliary trigger, and a sequence of triggers; likely h(t) triggers, within
 * a window length.
 * Additionally, since the input trigger is assumed to be from a sorted
 * sequence, the input trigger subsequence is pruned against triggers which
 * are guaranteed to no longer cause conincidences
 */
size_t _count_subsequence( GSequence* subseq, SnglBurst* trig, float wind ){
	GSequenceIter *tmp, *subitr = g_sequence_get_begin_iter(subseq);
	size_t counter = 0;
	while( !g_sequence_iter_is_end(subitr) ){
		SnglBurst *ht = (SnglBurst*)g_sequence_get(subitr);
		/*
		printf( "%d.%d %d.%d %f/%f ", trig->peak_time.gpsSeconds,
					trig->peak_time.gpsNanoSeconds,
					ht->peak_time.gpsSeconds,
					ht->peak_time.gpsNanoSeconds,
					diff, wind );
		*/
		// If |trig time - list| > counter increment
		if( fabs(XLALGPSDiff( &trig->peak_time, &ht->peak_time )) < wind ){
			counter++;
			subitr = g_sequence_iter_next(subitr);
		// If trig time > list + wind, trim list
		} else if( XLALGPSCmp( &trig->peak_time, &ht->peak_time ) == 1 ) {
			tmp = subitr;
			subitr = g_sequence_iter_next(subitr);
			g_sequence_remove(tmp);
		} else {
			//printf( "%d ", XLALGPSCmp( &trig->peak_time, &ht->peak_time ) );
			break;
			//subitr = g_sequence_iter_next(subitr);
		}
	}
	return counter;
}

/*
 * Scan through a list of triggers and count the instances of each trigger type
 * and count its coincidences with the target channel. The hash tables for the
 * channel count and coincidence count must already be initialized. The
 * trig_sequence parameter must contain a full list of the triggers to be
 * scanned. The chan parameter is the 'target' channel -- usually h(t). The
 * twind parameter controls the length of the coincidence window.
 *
 * coinc_type = 0: Allow all coincidences
 * coinc_type = 1: Allow only one coincidence between target channel and others
 * (e.g. if two triggers from the same channel are found in coincidence, only
 * one is recorded).
 */
void scan( GHashTable *chancount, GHashTable *chanhist, GSequence* trig_sequence, const char* chan, double twind, int coinctype ){

	/*
	 * pointer to the position in the list for the current auxiliary trigger
	 * being examined
	 */
	GSequenceIter *cur_trig = g_sequence_get_begin_iter( trig_sequence );

	LIGOTimeGPS start;
	LIGOTimeGPS stop;
	XLALGPSSetREAL8( &start, 0 );
	XLALGPSSetREAL8( &stop, 1 );

	// coincidence window
	LALSeg *wind;
	wind = XLALSegCreate( &start, &stop, 0 );

	// trigger pointer
	SnglBurst *sb_target;
	GHashTable *prevcoinc;

	// Iterate through trigger list -- it should be sorted by GPS time
	// TODO: Use foreach here instead?
	for( ; !g_sequence_iter_is_end(cur_trig) ; cur_trig = g_sequence_iter_next(cur_trig) ){

		// Current trigger
		sb_target = (SnglBurst*)g_sequence_get( cur_trig );

		/*
		 * Increment the trigger count for this channel.
		 * If it is not in the hash table, add it and initialize properly
		 */
		size_t *value;
		value = (size_t*)g_hash_table_lookup( chancount, sb_target->channel );
		if( value ){
			(*value)++;
			printf( "Count: Incrementing, %s value: %lu\n", sb_target->channel, *value );
			g_hash_table_insert( chancount, g_strdup(sb_target->channel), value );
		} else {
			value = g_new(size_t, 1);
			*value = 1;
			printf( "Count: Adding %s with time %d.%d\n", sb_target->channel, sb_target->peak_time.gpsSeconds, sb_target->peak_time.gpsNanoSeconds );
			g_hash_table_insert( chancount, g_strdup(sb_target->channel), value );
		}

		// Is it the channel we're looking at?
		// TODO: Consider doing this from the perspective of h(t), rather than
		// each aux channel, since we're checking coincidences of N to 1 rather
		// than 1 to N.
		// FIXME: check to make sure we don't have channels which are identical
		// up to a point in the string
		if( strstr( sb_target->channel, chan ) ){
			// Yes, create window segment
			start = stop = sb_target->peak_time;
			start = *XLALGPSAdd( &start, -twind/2.0 );
			stop = *XLALGPSAdd( &stop, twind/2.0 );
			printf( "creating segment from %s %d %d\n", sb_target->channel, start.gpsSeconds, stop.gpsSeconds );
			XLALSegSet( wind, &start, &stop, sb_target->event_id );
		} else { // No, go to the next
			continue;
		}

		// FIXME: Free memory
		prevcoinc = g_hash_table_new( g_str_hash, g_str_equal );

		printf( "Checking for event %d within %d %d\n", sb_target->peak_time.gpsSeconds, wind->start.gpsSeconds, wind->end.gpsSeconds );

		// This is our secondary pointer which counts out from the current
		// trigger
		GSequenceIter *trigp = cur_trig;
		SnglBurst* sb_aux;

		// Sweep backward, accumulate triggers in the window until we're outside
		// of it
		do {
			if( g_sequence_iter_is_begin(trigp) ) break;
			trigp = g_sequence_iter_prev(trigp);
			sb_aux = (SnglBurst*)g_sequence_get(trigp);

			// For now, don't use the target channel
			// FIXME: We may want this information in the future
			if( strstr( sb_aux->channel, chan ) ) continue;

			/*
			 * If we want unique coincidences, check to see if this channel has
			 * already been added.
			 */
			if( coinctype == 1 && g_hash_table_lookup( prevcoinc, sb_aux->channel ) ){
				 continue;
			} else {
				// FIXME: Use g_hash_table_add when compatible
				g_hash_table_insert( prevcoinc, &sb_aux->channel, &sb_aux->channel );
			}

			// TODO: Macroize?
			value = (size_t*)g_hash_table_lookup( chanhist, sb_aux->channel );
			// If we have an entry for this channel, use it, otherwise create a
			// new one
			if( value != NULL ){
				(*value)++;
				printf( "Coincidence: Incrementing, %s value: %lu\n", sb_aux->channel, *value );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			} else {
				value = g_new(size_t, 1);
				*value = 1;
				printf( "Coincidence: Adding %s with time %d.%d\n", sb_aux->channel, sb_aux->peak_time.gpsSeconds, sb_aux->peak_time.gpsNanoSeconds );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			}
		} while( XLALGPSInSeg( &sb_aux->peak_time, wind ) == 0 );

		// Sweep forward, accumulate triggers in the window until we're outside
		// of it
		do {
			trigp = g_sequence_iter_next(trigp);
			if( g_sequence_iter_is_end(trigp) ) break;
			sb_aux = (SnglBurst*)g_sequence_get(trigp);
			// Not the target channel?
			if( strstr( sb_aux->channel, chan ) ) continue;

			if( coinctype == 1 && g_hash_table_lookup( prevcoinc, sb_aux->channel ) ){
				 continue;
			} else {
				// FIXME: Use g_hash_table_add when compatible
				g_hash_table_insert( prevcoinc, &sb_aux->channel, &sb_aux->channel );
			}

			// TODO: Macroize?
			value = (size_t*)g_hash_table_lookup( chanhist, sb_aux->channel );
			// If we have an entry for this channel, use it, otherwise create a
			// new one
			if( value != NULL ){
				(*value)++;
				printf( "Coincidence: Incrementing, %s value: %lu\n", sb_aux->channel, *value );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			} else {
				value = g_new(size_t, 1);
				*value = 1;
				printf( "Coincidence: Adding %s with time %d.%d\n", sb_aux->channel, sb_aux->peak_time.gpsSeconds, sb_aux->peak_time.gpsNanoSeconds );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			}
		} while( XLALGPSInSeg( &sb_aux->peak_time, wind ) == 0 );
		g_hash_table_destroy(prevcoinc);
	}
	XLALFree( wind );
}

/*
 * Do a round of vetoes. In short, check the significance of each channel and
 * remove the highest signifiance channel, unless it falls below the
 * significance threshold. The winner is returned as the first parameter. The
 * chan parameter is the target channel. The t_ratio parameter is the ratio of
 * the veto window duration to the total examined livetime.
 */
double veto_round( char* winner, GHashTable* chancount, GHashTable* chanhist, const char* chan, double t_ratio ){
	double mu, sig, max_sig=-1;
	size_t *k;

	GHashTableIter iter;
	gpointer key, val;
	// Number of triggers in target channel
	size_t n_h = *(size_t *)g_hash_table_lookup(chancount, chan);

	g_hash_table_iter_init( &iter, chanhist );
	// Iterate over each channel and check the significance of the
	// accumulated trigger number
	while( g_hash_table_iter_next( &iter, &key, &val ) ){
		if( strstr( key, chan ) ){
			fprintf( stderr, "Skipping channel %s\n", (char *)key );
			continue;
		}
		// Number of triggers in auxillary channel
		size_t n_aux = *(size_t *)g_hash_table_lookup(chancount, key);
		mu = n_h*n_aux*t_ratio;

		k = (size_t *)val;

		fprintf( stderr, "Total coincidences for channel %s: %lu\n", (char *)key, *k );
		fprintf( stderr, "Mu for channel %s: %g\n", (char *)key, mu );
		sig = significance( mu, *k );
		fprintf( stderr, "Significance for this channel: %g\n", sig );
		if( sig > max_sig ){
				max_sig = sig;
				strcpy( winner, (char *)key );
		}
	}
	printf( "winner: %s\n", winner );
	return max_sig;
}

void prune_trigs( GSequence* trig_sequence, const LALSegList* onsource ){
	// FIXME: Actually, this should prune all triggers
	if( onsource->length == 0 ){
		return;
	}
	size_t i = 0;

	GSequenceIter *trigp, *tmp;
	trigp = g_sequence_get_begin_iter( trig_sequence );

	LALSeg onseg = onsource->segs[0];
	while( !g_sequence_iter_is_end( trigp ) ){
		SnglBurst* sb = g_sequence_get( trigp );
		int pos = XLALGPSInSeg( &sb->peak_time, &onseg );

		/*
		 * The case for this sentinel:
		 * Last trigger was within the segment previous to this
		 * This trigger is in the subsequent segment, thus we need to advance
		 * the pointer to account for it
		 */
		if( pos > 0 && i < onsource->length ){
			onseg = onsource->segs[i++];
			continue;
		/*
		 * We're beyond the extent of the onsource list
		 */
		} else if( pos > 0 && i >= onsource->length ) {
			tmp = trigp;
			trigp = g_sequence_iter_next( trigp );
			g_sequence_remove(tmp);
		}

		if( pos != 0 ){
			tmp = trigp;
			trigp = g_sequence_iter_next( trigp );
			g_sequence_remove(tmp);
		} else {
			trigp = g_sequence_iter_next( trigp );
		}
	}
}

/*
 * Remove triggers in window centered around peak times for triggers from vchan.
 * This is generally called after a veto round.
 *
 * TODO: Merge vetolist creation here
 * TODO: Can we also decrement the count / coincidences efficiently here?
 */
size_t remove_trigs( GSequence* trig_sequence, const LALSeg veto, const char* vchan ){

	size_t vetoed_events = 0;
	size_t nevents = g_sequence_get_length(trig_sequence);
	fprintf( stderr, "nevents: %lu\n", nevents );
	fprintf( stderr, "Channel to veto: %s\n", vchan );

	// Pointer to the current position in the trigger list
	GSequenceIter* trigp = g_sequence_get_begin_iter(trig_sequence);
	// Pointer to the trigger under examination
	SnglBurst* sb;

	// Reset to beginning for real delete pass
	trigp = g_sequence_get_begin_iter(trig_sequence);

	// Loop over our list, looking to remove trigs
	while( !g_sequence_iter_is_end(trigp) ){
		sb = (SnglBurst*)g_sequence_get(trigp);
		if( !sb ){
			fprintf( stderr, "Invalid pointer for top level iterator!\n" );
		}

		// Is it the channel to veto on?
		if( !strstr( sb->channel, vchan ) ){
			trigp = g_sequence_iter_next(trigp);
			continue; // no, move on
		}
		LALSeg *trig_veto = XLALSegCreate( &veto.start, &veto.end, (int)sb->event_id );
		// Center the window on the trigger
		XLALGPSAddGPS( &trig_veto->start, &sb->peak_time);
		XLALGPSAddGPS( &trig_veto->end, &sb->peak_time);

		GSequenceIter *st = trigp, *end = trigp;

		gboolean begin = g_sequence_iter_is_begin(st);
		if( !begin ){
			st = g_sequence_iter_prev(st);
			sb = (SnglBurst*)g_sequence_get(st);
		}

		// Backwards
		while( !begin & (XLALGPSInSeg( &sb->peak_time, trig_veto ) == 0) ){
			GSequenceIter *tmp;
			tmp = st;
			if( !g_sequence_iter_is_begin(st) ){
				st = g_sequence_iter_prev(st);
				sb = (SnglBurst*)g_sequence_get(st);
			}
			g_sequence_remove(tmp);
			vetoed_events++;
		}

		// Check to make sure that we're not at the end of the list
		if( g_sequence_iter_is_end(end) ){
			break;
		} else {
			sb = (SnglBurst*)g_sequence_get(end);
		}

		// Forwards
		while( XLALGPSInSeg( &sb->peak_time, trig_veto ) == 0 ){
			// don't invalidate the top level iterator
			GSequenceIter *tmp;
			tmp = end;
			trigp = end = g_sequence_iter_next(end);
			// Delete this trigger
			g_sequence_remove(tmp);
			vetoed_events++;
			if( g_sequence_iter_is_end(end) ){
				break;
			}
			sb = (SnglBurst*)g_sequence_get(end);
		}

		// FIXME: Add to veto list
		XLALFree( trig_veto );

		fprintf( stderr, "%lu events deleted so far.\n", vetoed_events );
		nevents = g_sequence_get_length(trig_sequence);
		fprintf( stderr, "%lu events remain\n", nevents );
	}
	fprintf( stderr, "Done, total events removed %lu\n", vetoed_events );

	return vetoed_events;
}

/*
 * Turn all the peak times for channel vchan into a segment list of vetoes.
 */
void trigs_to_vetolist( LALSegList* vetoes, GSequence* trig_sequence, const LALSeg veto, const char* vchan ){

	float wind = XLALGPSDiff(&veto.end, &veto.start);

	// Pointer to the current position in the trigger list
	GSequenceIter* trigp = g_sequence_get_begin_iter(trig_sequence);
	// Pointer to the trigger under examination
	SnglBurst* sb;

	while( !g_sequence_iter_is_end(trigp) ){
		sb = (SnglBurst*)g_sequence_get(trigp);
		if( !sb ){
			fprintf( stderr, "Invalid pointer for top level iterator!\n" );
		}

		if( strstr( sb->channel, vchan ) ){
            LALSeg vetotmp;
            LIGOTimeGPS start = sb->peak_time;
            LIGOTimeGPS stop = sb->peak_time;
            XLALGPSSetREAL8( &start, -wind/2.0 );
            XLALGPSSetREAL8( &stop, wind/2.0 );
            XLALSegSet( &vetotmp, &start, &stop, sb->event_id );
			XLALSegListAppend( vetoes, &vetotmp );
		}
		trigp = g_sequence_iter_next(trigp);
	}
}

/*
 * Calculate the signifiance of a set of triggers from the Poisson survival
 * function given the expected number of triggers.
 */
double significance( double mu, int k ){
	double sig = -log10( gsl_sf_gamma_inc_P(mu, k) );
	// FIXME: Arbitrary
	if( sig < 1e-15 ){
		return -k*log10(mu) + mu*log10(exp(1)) + gsl_sf_lngamma(k+1)/log(10);
	}
	return sig;
}
