#include <lal/LALDetCharHveto.h>

/*
 * Okay, now we're playing a bit fast and loose. These are defined for use with
 * the SnglBurst's next pointer. Since we don't use the SnglBurst as a linked
 * list --- leaving that instead to the glib sequence --- the next pointer is
 * hijacked for use with "masking" trigs via the Prune function. Since several
 * rounds are usually required, when we "prune" a trig via its SNR, we don't
 * want to remove and then reload it later --- that takes a long time. Instead,
 * we "mask" it by assigning its next pointer one of the below values. Both
 * are near guaranteed to crash the program if someone attempts to use them
 * in the normal fashion, and are informative if looking at them in gdb.
 *
 * In short, don't use the SnglBurst next pointer for something in the hveto 
 * program. I hereby disclaim responsibility if by some magic, the program
 * actually allocates to 0xDEADBEEFDEADBEFF.
 */
#define TRIG_MASKED (SnglBurst*)0xDEADBEEFDEADBEEF
#define TRIG_NOT_MASKED (SnglBurst*)0x0

/*
static gint compare(gconstpointer a, gconstpointer b) {
        const SnglBurst *_a = a;
        const SnglBurst *_b = b;

        return XLALCompareSnglBurstByPeakTimeAndSNR(&_a, &_b);
}
*/

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
void XLALDetCharScanTrigs( GHashTable *chancount, GHashTable *chanhist, GSequence* trig_sequence, const char* chan, double twind, int coinctype ){

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
		if( sb_target->next == TRIG_MASKED ) continue;

		/*
		 * Increment the trigger count for this channel.
		 * If it is not in the hash table, add it and initialize properly
		 */
		size_t *value;
		value = (size_t*)g_hash_table_lookup( chancount, sb_target->channel );
		if( value ){
			(*value)++;
			XLALPrintInfo( "Count: Incrementing, %s value: %lu\n", sb_target->channel, *value );
			g_hash_table_insert( chancount, g_strdup(sb_target->channel), value );
		} else {
			value = g_new(size_t, 1);
			*value = 1;
			XLALPrintInfo( "Count: Adding %s with time %d.%d\n", sb_target->channel, sb_target->peak_time.gpsSeconds, sb_target->peak_time.gpsNanoSeconds );
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
			XLALPrintInfo( "creating segment from %s %d.%d %d.%d\n", sb_target->channel, start.gpsSeconds, start.gpsNanoSeconds, stop.gpsSeconds, stop.gpsNanoSeconds  );
			XLALSegSet( wind, &start, &stop, sb_target->event_id );
		} else { // No, go to the next
			continue;
		}

		// FIXME: Free memory
		prevcoinc = g_hash_table_new( g_str_hash, g_str_equal );

		XLALPrintInfo( "Checking for event %d within %d %d\n", sb_target->peak_time.gpsSeconds, wind->start.gpsSeconds, wind->end.gpsSeconds );

		// This is our secondary pointer which counts out from the current
		// trigger
		GSequenceIter *trigp = cur_trig;
		SnglBurst* sb_aux;
		gboolean begin = FALSE;
		if( g_sequence_iter_is_begin(trigp) ){
			begin = TRUE;
		} else {
			trigp = g_sequence_iter_prev(trigp);
			sb_aux = (SnglBurst*)g_sequence_get(trigp);
		}

		// Sweep backward, accumulate triggers in the window until we're outside
		// of it
		while( !begin && XLALGPSInSeg( &sb_aux->peak_time, wind ) == 0 ){

			/*
			 * If we want unique coincidences, check to see if this channel has
			 * already been added.
			 *
			 * For now, don't use the target channel
			 * FIXME: We may want this information in the future
			 */
			if( (coinctype == 1 && g_hash_table_lookup( prevcoinc, sb_aux->channel )) || strstr( sb_aux->channel, chan ) || sb_aux->next == TRIG_MASKED ){
				if( g_sequence_iter_is_begin(trigp) ) break;
				trigp = g_sequence_iter_prev(trigp);
				sb_aux = (SnglBurst*)g_sequence_get(trigp);
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
				XLALPrintInfo( "Left Coincidence: Incrementing, %s->%ld, time %d.%d value: %lu\n", sb_aux->channel, sb_target->event_id, sb_aux->peak_time.gpsSeconds, sb_aux->peak_time.gpsNanoSeconds, *value );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			} else {
				value = g_new(size_t, 1);
				*value = 1;
				XLALPrintInfo( "Left Coincidence: Adding %s->%ld with time %d.%d\n", sb_aux->channel, sb_target->event_id, sb_aux->peak_time.gpsSeconds, sb_aux->peak_time.gpsNanoSeconds );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			}
			if( g_sequence_iter_is_begin(trigp) ) break;
			trigp = g_sequence_iter_prev(trigp);
			sb_aux = (SnglBurst*)g_sequence_get(trigp);
		}

		trigp = cur_trig;
		trigp = g_sequence_iter_next(trigp);
		if( g_sequence_iter_is_end(trigp) ) break;
		sb_aux = (SnglBurst*)g_sequence_get(trigp);

		// Sweep forward, accumulate triggers in the window until we're outside
		// of it
		//do {
		while( XLALGPSInSeg( &sb_aux->peak_time, wind ) == 0 ){
			//trigp = g_sequence_iter_next(trigp);
			//if( g_sequence_iter_is_end(trigp) ) break;
			//sb_aux = (SnglBurst*)g_sequence_get(trigp);
			
			// Not the target channel?
			if( strstr( sb_aux->channel, chan ) || sb_aux->next == TRIG_MASKED ){ 
				trigp = g_sequence_iter_next(trigp);
				if( g_sequence_iter_is_end(trigp) ) break;
				sb_aux = (SnglBurst*)g_sequence_get(trigp);
				continue;
			}

			if( coinctype == 1 && g_hash_table_lookup( prevcoinc, sb_aux->channel ) ){
				trigp = g_sequence_iter_next(trigp);
				if( g_sequence_iter_is_end(trigp) ) break;
				sb_aux = (SnglBurst*)g_sequence_get(trigp);
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
				XLALPrintInfo( "Right Coincidence: Incrementing, %s->%ld, time %d.%d value: %lu\n", sb_aux->channel, sb_target->event_id, sb_aux->peak_time.gpsSeconds, sb_aux->peak_time.gpsNanoSeconds, *value );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			} else {
				value = g_new(size_t, 1);
				*value = 1;
				XLALPrintInfo( "Right Coincidence: Adding %s->%ld with time %d.%d\n", sb_aux->channel, sb_target->event_id, sb_aux->peak_time.gpsSeconds, sb_aux->peak_time.gpsNanoSeconds );
				g_hash_table_insert( chanhist, &sb_aux->channel, value );
			}
			trigp = g_sequence_iter_next(trigp);
			if( g_sequence_iter_is_end(trigp) ) break;
			sb_aux = (SnglBurst*)g_sequence_get(trigp);
		}
		//} while( XLALGPSInSeg( &sb_aux->peak_time, wind ) == 0 );
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
double XLALDetCharVetoRound( char* winner, GHashTable* chancount, GHashTable* chanhist, const char* chan, double t_ratio ){
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
			//XLALPrintInfo( stderr, "Skipping channel %s\n", (char *)key );
			continue;
		}
		// Number of triggers in auxillary channel
		size_t n_aux = *(size_t *)g_hash_table_lookup(chancount, key);
		mu = n_h*n_aux*t_ratio;

		k = (size_t *)val;

		printf( "Total coincidences for channel %s: %zu\n", (char *)key, *k );
		printf( "Mu for channel %s: %g\n", (char *)key, mu );
		sig = XLALDetCharHvetoSignificance( mu, *k );
		printf( "Significance for this channel: %g\n", sig );
		if( sig > max_sig && !strstr(chan, (char*)key) ){
				max_sig = sig;
				strcpy( winner, (char *)key );
		}
	}
	printf( "winner: %s\n", winner );
	return max_sig;
}

void XLALDetCharPruneTrigs( GSequence* trig_sequence, const LALSegList* onsource, double snr_thresh, const char* refchan ){
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
			i++;
			onseg = onsource->segs[i];
			continue;
		/*
		 * We're beyond the extent of the onsource list
		 */
		} else if( pos > 0 && i >= onsource->length ) {
			tmp = trigp;
			trigp = g_sequence_iter_next( trigp );
			g_sequence_remove(tmp);
			continue;
		}

		if( pos != 0 ){
			tmp = trigp;
			trigp = g_sequence_iter_next( trigp );
			g_sequence_remove(tmp);
			continue;
		}

		/*
		 * This is the reference channel, and we don't check it's snr to remove
		 * it. If no refchan is given, then remove anything according to 
		 * criteria.
		 */
		gboolean isrefchan = (refchan != NULL);
		if( isrefchan ) isrefchan = (strstr( refchan, sb->channel ) != NULL);

		/*
		 * Mark the trigger as unused, but don't delete it. This is a common
		 * where the snr threshold is raised, but we don't want to remove the
		 * trigger. See warnings and disclaimers at the top of the file.
		 */
		if( sb->snr < snr_thresh && !isrefchan ){
			sb->next = TRIG_MASKED;
		} else {
			sb->next = TRIG_NOT_MASKED;
		}
		trigp = g_sequence_iter_next( trigp );
	}
}

/*
 * Remove triggers in window centered around peak times for triggers from vchan.
 * This is generally called after a veto round.
 *
 * TODO: Merge vetolist creation here
 * TODO: Can we also decrement the count / coincidences efficiently here?
 */
GSequence* XLALDetCharRemoveTrigs( GSequence* trig_sequence, const LALSeg veto, const char* vchan, double snr_thresh ){

	char refchan[] = "LSC-DARM_ERR";
	size_t vetoed_events = 0;
	size_t de_vetoed_events = 0;
	size_t nevents = g_sequence_get_length(trig_sequence);
	XLALPrintInfo( "nevents: %lu\n", nevents );
	XLALPrintInfo( "Channel to veto: %s\n", vchan );

	// Pointer to the current position in the trigger list
	GSequenceIter* trigp = g_sequence_get_begin_iter(trig_sequence);
	// Pointer to the trigger under examination
	SnglBurst* sb;

	// Store the triggers to be deleted
	GSequence* tbd = g_sequence_new((GDestroyNotify)XLALDestroySnglBurst);
	GSequenceIter* tbdit = g_sequence_get_begin_iter(tbd);

	// Loop over our list, looking to remove trigs
	while( !g_sequence_iter_is_end(trigp) ){
		sb = (SnglBurst*)g_sequence_get(trigp);
		if( !sb ){
			XLALPrintError( "Invalid pointer for top level iterator!\n" );
		}

		// Is it the channel to veto on?
		if( (sb->snr < snr_thresh) | !strstr( sb->channel, vchan ) ){
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
			if( g_sequence_iter_is_begin(st) ){
				break;
			}
			st = g_sequence_iter_prev(st);
			g_sequence_move(tmp, tbdit);
			tbdit = g_sequence_iter_next(tbdit);
			vetoed_events++;
			XLALPrintInfo( "Backwards, deleting %s (%d.%d) id %ld\n", sb->channel, sb->peak_time.gpsSeconds, sb->peak_time.gpsNanoSeconds, sb->event_id );
			if( strstr(refchan, sb->channel) ){
				de_vetoed_events++;
			}
			sb = (SnglBurst*)g_sequence_get(st);
		}

		// Check to make sure that we're not at the end of the list
		if( g_sequence_iter_is_end(end) ){
			break;
		} else {
			sb = (SnglBurst*)g_sequence_get(end);
		}

		GSequenceIter *tmp;
		// Forwards
		while( XLALGPSInSeg( &sb->peak_time, trig_veto ) == 0 ){
			// don't invalidate the top level iterator
			tmp = end;
			end = g_sequence_iter_next(end);
			// Delete this trigger -- but don't delete anything we'd use later
			if( !strstr(sb->channel, vchan) ){
				g_sequence_move(tmp, tbdit);
				tbdit = g_sequence_iter_next(tbdit);
				vetoed_events++;
				XLALPrintInfo( "Forwards, deleting %s (%d.%d) id %ld\n", sb->channel, sb->peak_time.gpsSeconds, sb->peak_time.gpsNanoSeconds, sb->event_id );
				if( strstr(refchan, sb->channel) ){
					de_vetoed_events++;
				}
			}

			if( g_sequence_iter_is_end(end) ){
				break;
			}
			sb = (SnglBurst*)g_sequence_get(end);
		}
		tmp = trigp;
		sb = (SnglBurst*)g_sequence_get(tmp);
		trigp = g_sequence_iter_next(trigp);
		g_sequence_move( tmp, tbdit );
		tbdit = g_sequence_iter_next(tbdit);
		vetoed_events++;
		XLALPrintInfo( "Veto trig, deleting %s (%d.%d) id %ld\n", sb->channel, sb->peak_time.gpsSeconds, sb->peak_time.gpsNanoSeconds, sb->event_id );

		// FIXME: Add to veto list
		XLALFree( trig_veto );

		XLALPrintInfo( "%lu events deleted so far.\n", vetoed_events );
		nevents = g_sequence_get_length(trig_sequence);
		XLALPrintInfo( "%lu events remain\n", nevents );
	}
	XLALPrintInfo( "Done, total events removed %lu\n", vetoed_events );
	XLALPrintInfo( "Done, ref channel total events removed %lu\n", de_vetoed_events );

	return tbd;
}

/*
 * Turn all the peak times for channel vchan into a segment list of vetoes.
 */
void XLALDetCharTrigsToVetoList( LALSegList* vetoes, GSequence* trig_sequence, const LALSeg veto, const char* vchan ){

	float wind = XLALGPSDiff(&veto.end, &veto.start);

	// Pointer to the current position in the trigger list
	GSequenceIter* trigp = g_sequence_get_begin_iter(trig_sequence);
	// Pointer to the trigger under examination
	SnglBurst* sb;

	while( !g_sequence_iter_is_end(trigp) ){
		sb = (SnglBurst*)g_sequence_get(trigp);
		if( !sb ){
			XLALPrintError( "Invalid pointer for top level iterator!\n" );
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
 * Calculate the significance of a set of triggers from the Poisson survival
 * function given the expected number of triggers.
 */
double XLALDetCharHvetoSignificance( double mu, int k ){
	double sig = gsl_sf_gamma_inc_P(k, mu);
	if( sig != 0.0 ){
		return -log10(sig);
	} else {
		return -k*log10(mu) + mu*log10(exp(1)) + gsl_sf_lngamma(k+1)/log(10);
	}
}
