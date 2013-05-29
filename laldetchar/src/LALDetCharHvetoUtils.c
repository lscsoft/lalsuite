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

#include <lal/LALDetCharHvetoUtils.h>

/* 
 * SWIG refuses to let me typemap out double*, it thinks its a REAL8*, and it 
 * won't let me redefine that either... so... don't give this a key that would 
 * result in NULL
 */
double XLALGLibDoubleFromGHash( const char* key, GHashTable* tbl ){
	return *(double*)g_hash_table_lookup( tbl, key );
}
/* 
 * ditto
 */
size_t XLALGLibSizeTFromGHash( const char* key, GHashTable* tbl ){
	return *(size_t*)g_hash_table_lookup( tbl, key );
}

char* XLALGLibStringFromGHash( const char* key, GHashTable* tbl ){
	return (char*)g_hash_table_lookup( tbl, key );
}

/*
 * For some reason I can't return a char**, the SWIG bindings won't allow it
 */
char* XLALGLibKeyFromGHash( GHashTable* tbl, size_t indx ){
	GList* keys = g_hash_table_get_keys( tbl );
	size_t i = 0;
	while( keys ){
		if( i == indx ) return g_strdup( (char*)keys->data );
		i++;
		keys = keys->next;
	}
	return NULL;
}

void XLALGLibSeqDelete( GSequence* seq ){
	g_sequence_free( seq );
}

size_t XLALGLibSeqGetLength( GSequence* seq ){
	return g_sequence_get_length( seq );
}

size_t XLALGLibSeqGetUnmarkedLength( GSequence* seq ){
    GSequenceIter* itr = g_sequence_get_begin_iter(seq);

	size_t cntr = 0;
	while( !g_sequence_iter_is_end(itr) ){
		SnglBurst* sb = (SnglBurst*)g_sequence_get(itr);
		if( sb->next == NULL ){
			cntr++;
		}
		itr = g_sequence_iter_next(itr);
	}
	return cntr;
}

GHashTable* XLALGetChannelList( GSequence *triglist ){
    GSequenceIter* itr = g_sequence_get_begin_iter(triglist);

    GHashTable* channellist = g_hash_table_new( g_str_hash, g_str_equal );
    while( !g_sequence_iter_is_end(itr) ){
        SnglBurst *sb = (SnglBurst*)g_sequence_get(itr);
        g_hash_table_insert( channellist, g_strdup(sb->channel), g_strdup(sb->channel) );
        itr = g_sequence_iter_next(itr);
    }
    return channellist;
}

GHashTable* XLALCreateGHashTable( void ){
	GHashTable *tmp = g_hash_table_new( g_str_hash, g_str_equal );
	return tmp;
}

struct tagSnglBurst* XLALGLibSBListFromSeq( GSequence* in ){
    SnglBurst* outlist = NULL;
    SnglBurst* begin = NULL;
    GSequenceIter* iter = g_sequence_get_begin_iter(in);
    while( !g_sequence_iter_is_end(iter) ){
        SnglBurst *sb = (SnglBurst*)g_sequence_get(iter);
		if( !outlist ){
			begin = outlist = sb;
		} else {
			outlist = (outlist->next = sb);
		}
        iter = g_sequence_iter_next(iter);
    }
    //g_sequence_free(in);
	if( outlist ){ 
		outlist->next = NULL; 
	}
    return begin;
}

GSequence* XLALPopulateTrigSequenceFromFile( const char* fname, double min_snr, char* ignore_list, GSequence* trig_sequence ){
	if( !trig_sequence ){
    	trig_sequence = g_sequence_new(NULL);
	}
    GSequence* ignorel = g_sequence_new(free);

    SnglBurst* tbl = XLALSnglBurstTableFromLIGOLw( fname );
    SnglBurst *begin = NULL, *deleteme = NULL;
    if( !tbl ) {
        return trig_sequence;
    }

    int cnt = 0;
    char* tmp;
    tmp = malloc( sizeof(char)*512 );
    FILE* lfile = fopen( ignore_list, "r" );
    while(!feof(lfile)){
        cnt = fscanf( lfile, "%s", tmp );
        if( cnt == EOF ) break;
        g_sequence_append( ignorel, g_strdup(tmp) );
    }
    free(tmp);
    fclose(lfile);

    do {
        gboolean ignore = FALSE;
        GSequenceIter* igitr;
        // FIXME: Support ignorelist == NULL
        igitr = ignorel ? g_sequence_get_begin_iter(ignorel) : NULL;
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
        if( tbl->snr >= min_snr && !ignore ){
            //printf( "Adding event %p #%lu, channel: %s\n", tbl, tbl->event_id, tbl->channel );
            g_sequence_insert_sorted( trig_sequence, tbl, (GCompareDataFunc)compare, NULL );
            tbl=tbl->next;
        } else {
            //printf( "Ignoring event %p #%lu, channel: %s\n", tbl, tbl->event_id, tbl->channel );
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
    XLALPrintInfo( "Deleting %d unused events.\n", XLALSnglBurstTableLength(begin) );
	if( deleteme ){
    	deleteme->next = NULL; // Detach this from its sucessor in case that's used
    	deleteme = NULL;
	}
	if( begin ){
    	XLALDestroySnglBurstTable(begin);
	}	
    XLALPrintInfo( "Done.\n" );

    g_sequence_free( ignorel );

    return trig_sequence;
}

GSequence* XLALPopulateTrigSequenceFromTrigList( SnglBurst* tbl, GSequence* trig_sequence ){
	if( !trig_sequence ){
    	trig_sequence = g_sequence_new(NULL);
	}

    if( !tbl ) {
        return trig_sequence;
    }

    do {
        g_sequence_insert_sorted( trig_sequence, tbl, (GCompareDataFunc)compare, NULL );
        tbl=tbl->next;
    } while( tbl );

    return trig_sequence;
}
