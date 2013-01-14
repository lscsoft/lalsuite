#include <lal/LALDetCharHveto.h>

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
		double diff = fabs(XLALGPSDiff( &trig->peak_time, &ht->peak_time ));
		printf( "%d.%d %d.%d %f/%f ", trig->peak_time.gpsSeconds,
					trig->peak_time.gpsNanoSeconds,
					ht->peak_time.gpsSeconds,
					ht->peak_time.gpsNanoSeconds,
					diff, wind );
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
			printf( "%d ", XLALGPSCmp( &trig->peak_time, &ht->peak_time ) );
			printf( "%d\n", counter );
			break;
			//subitr = g_sequence_iter_next(subitr);
		}
		printf( "%d\n", counter );
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
 */
void scan( GHashTable *chancount, GHashTable *chanhist, GSequence* trig_sequence, const char* chan, double twind ){

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
	SnglBurst *sb;

	// Iterate through trigger list -- it should be sorted by GPS time
	// TODO: Use foreach here instead?
	for( ; !g_sequence_iter_is_end(cur_trig) ; cur_trig = g_sequence_iter_next(cur_trig) ){

		// Current trigger
		sb = (SnglBurst*)g_sequence_get( cur_trig );

		size_t *value;
		value = (size_t*)g_hash_table_lookup( chancount, sb->channel );
		/*
		 * Increment the trigger count for this channel.
		 * If it is not in the hash table, add it and initialize properly
		 */
		if( value ){
			(*value)++;
			printf( "Count: Incrementing, %s value: %lu\n", sb->channel, *value );
			g_hash_table_insert( chancount, g_strdup(sb->channel), value );
		} else {
			value = g_new(size_t, 1);
			*value = 1;
			printf( "Count: Adding %s with time %d.%d\n", sb->channel, sb->peak_time.gpsSeconds, sb->peak_time.gpsNanoSeconds );
			g_hash_table_insert( chancount, g_strdup(sb->channel), value );
		}

		// Is it the channel we're looking at?
		// TODO: Consider doing this from the perspective of h(t), rather than
		// each aux channel, since we're checking coincidences of N to 1 rather
		// than 1 to N.
		if( !strstr( sb->channel, chan ) ){
			// No, create window segment
			start = stop = sb->peak_time;
			start = *XLALGPSAdd( &start, -twind/2.0 );
			stop = *XLALGPSAdd( &stop, twind/2.0 );
			printf( "creating segment from %s %d %d\n", sb->channel, start.gpsSeconds, stop.gpsSeconds );
			XLALSegSet( wind, &start, &stop, sb->event_id );
		} else { // No, go to the next
			continue;
		}

		printf( "Checking for event %d within %d %d\n", sb->peak_time.gpsSeconds, wind->start.gpsSeconds, wind->end.gpsSeconds );

		// This is our secondary pointer which counts out from the current
		// trigger
		GSequenceIter *trigp = cur_trig;
		SnglBurst* sb_h;

		// Sweep backward, accumulate triggers in the window until we're outside
		// of it
		do {
			if( g_sequence_iter_is_begin(trigp) ) break;
			trigp = g_sequence_iter_prev(trigp);
			sb_h = (SnglBurst*)g_sequence_get(trigp);
			// Not the target channel?
			if( !strstr( sb_h->channel, chan ) ) continue;

			// TODO: Macroize?
			value = (size_t*)g_hash_table_lookup( chanhist, sb->channel );
			// If we have an entry for this channel, use it, otherwise create a
			// new one
			// TODO: Move this outside the do/while
			if( value != NULL ){
				(*value)++;
				printf( "Coincidence: Incrementing, %s value: %lu\n", sb->channel, *value );
				// FIXME: Necessary?
				g_hash_table_insert( chanhist, &sb->channel, value );
			} else {
				value = g_new(size_t, 1);
				*value = 1;
				printf( "Coincidence: Adding %s with time %d.%d\n", sb->channel, sb->peak_time.gpsSeconds, sb->peak_time.gpsNanoSeconds );
				g_hash_table_insert( chanhist, &sb->channel, value );
			}
		} while( XLALGPSInSeg( &sb_h->peak_time, wind ) == 0 );

		// Sweep forward, accumulate triggers in the window until we're outside
		// of it
		do {
			trigp = g_sequence_iter_next(trigp);
			if( g_sequence_iter_is_end(trigp) ) break;
			sb_h = (SnglBurst*)g_sequence_get(trigp);
			// Not the target channel?
			if( !strstr( sb_h->channel, chan ) ) continue;

			// TODO: Macroize?
			value = (size_t*)g_hash_table_lookup( chanhist, sb->channel );
			// If we have an entry for this channel, use it, otherwise create a
			// new one
			// TODO: Move this outside the do/while
			if( value != NULL ){
				(*value)++;
				printf( "Coincidence: Incrementing, %s value: %lu\n", sb->channel, *value );
				// FIXME: Necessary?
				g_hash_table_insert( chanhist, &sb->channel, value );
			} else {
				value = g_new(size_t, 1);
				*value = 1;
				printf( "Coincidence: Adding %s with time %d.%d\n", sb->channel, sb->peak_time.gpsSeconds, sb->peak_time.gpsNanoSeconds );
				g_hash_table_insert( chanhist, &sb->channel, value );
			}
		} while( XLALGPSInSeg( &sb_h->peak_time, wind ) == 0 );
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
 */
// TODO: Change to "veto_trigs" and return veto list
size_t remove_trigs( GSequence* trig_sequence, const LALSeg veto, GHashTable *chancount, GHashTable *chanhist, const char* vchan, const char* refchan ){

	size_t vetoed_events = 0;
	size_t nevents = g_sequence_get_length(trig_sequence);
	fprintf( stderr, "nevents: %lu\n", nevents );
	fprintf( stderr, "Channel to veto: %s\n", vchan );
	volatile float window = XLALGPSDiff(&veto.end, &veto.start);

	// Buffer of recent h(t) trigs for use to count how many coincidences
	// we lose when deleting triggers
	GSequence* h_trigs = g_sequence_new((GDestroyNotify)free);

	// Pointer to the current position in the trigger list
	GSequenceIter* trigp = g_sequence_get_begin_iter(trig_sequence);
	// Pointer to the trigger under examination
	SnglBurst* sb;

	/*
	 * TODO: We can probably avoid this by keeping track not only of the
	 * h(t) triggers we've examined but also those slightly ahead of us, but
	 * keeping all of them probably isn't too much burden since
	 * _count_subsequence loops fairly efficiently. It's just a memory drain.
	 */
	while( !g_sequence_iter_is_end(trigp) ){
		sb = (SnglBurst*)g_sequence_get(trigp);
		if( !sb ){
			fprintf( stderr, "Invalid pointer for top level iterator!\n" );
		}
		trigp = g_sequence_iter_next(trigp);

		/*
		 * For the record of recent h(t) trigs, we need to copy the memory
		 * or else we'll be reading garbage in _count_subsequence if the
		 * trigger is deleted by this function first.
		 *
		 * Note: _count_subsequence will free this memory as necessary
		 * and the rest will be deleted at the end of the function
		 */
		if( strstr( sb->channel, refchan ) ){
			SnglBurst *tcpy = malloc( sizeof(SnglBurst) );
			memcpy( tcpy, sb, sizeof(SnglBurst) );
			g_sequence_append( h_trigs, tcpy );
		}
	}

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

		volatile GSequenceIter *st = trigp, *end = trigp;
		volatile GSequenceIter *tmp;
		size_t *val;

		// Backwards
		while( XLALGPSInSeg( &sb->peak_time, trig_veto ) == 0 ){
			sb = (SnglBurst*)g_sequence_get(st);
			if( !sb ){
				fprintf( stderr, "BAD 2!\n" );
			}
			/*
			 * The trigger itself will be removed by the forwards sweep because
			 * it needs to increment the top level pointer as well. So we skip
			 * it on this pass.
			 */
			if( st != trigp ){
				val = g_hash_table_lookup( chancount, sb->channel );
				(*val)--;
				fprintf( stderr, "Left Count: Decrementing %s, value: %lu\n", sb->channel, *val );
				// TODO: check what it means for one aux trigger to be in
				// coincidence with several h(t) triggers
				val = g_hash_table_lookup( chanhist, sb->channel );
				if( val ){
					(*val) -= _count_subsequence( h_trigs, sb, window/2.0 );
					fprintf( stderr, "Left Coincidence: Decrementing %s, value: %lu\n", vchan, *val );
				}
			}

			// Check if we've hit the beginning of the list
			if( g_sequence_iter_is_begin(st) ){
				// If this is the first trigger in the list, don't remove it
				break;
			} else {
				// Increment the pointer so as not to invalidate it when we
				// remove the trigger
				tmp = st;
				st = g_sequence_iter_prev(st);
			}
			// If we're not saving it for the forward sweep, remove it.
			if( trigp != tmp ){
				g_sequence_remove(tmp);
				vetoed_events++;
			}
		}

		// Check to make sure that we're not at the end of the list
		if( g_sequence_iter_is_end(end) ){
			break;
		} else {
			sb = (SnglBurst*)g_sequence_get(end);
		}

		// Forwards
		while( XLALGPSInSeg( &sb->peak_time, trig_veto ) == 0 ){

			sb = (SnglBurst*)g_sequence_get(end);
			if( !sb ){
				fprintf( stderr, "BAD 4!\n" );
			}
			val = g_hash_table_lookup( chancount, sb->channel );
			(*val)--;
			fprintf( stderr, "Right Count: Decrementing %s, value: %lu\n", sb->channel, *val );

			val = g_hash_table_lookup( chanhist, sb->channel );
			// TODO: check what it means for one aux trigger to be in
			// coincidence with several h(t) triggers
			if( val ){
				(*val) -= _count_subsequence( h_trigs, sb, window/2.0 );
				fprintf( stderr, "Right Coincidence: Decrementing %s, value: %lu\n", sb->channel, *val );
			}

			vetoed_events++;
			// don't invalidate the top level iterator
			tmp = end;
			// Move the top level pointer over so that doesn't get invalidated
			trigp = end = g_sequence_iter_next(end);
			// Delete this trigger
			g_sequence_remove(tmp);
			if( g_sequence_iter_is_end(end) ){
				break;
			}
		}

		// FIXME: Add to veto list
		XLALFree( trig_veto );

		fprintf( stderr, "%lu events deleted so far.\n", vetoed_events );
		nevents = g_sequence_get_length(trig_sequence);
		fprintf( stderr, "%lu events remain\n", nevents );
	}
	fprintf( stderr, "Done, total events removed %lu\n", vetoed_events );
	g_sequence_free( h_trigs );

	return vetoed_events;
}

/*
 * Calculate the signifiance of a set of triggers from the Poisson survival
 * function given the expected number of triggers.
 */
double significance( double mu, int k ){
	return -log10( gsl_sf_gamma_inc_P(mu, k) );
}
