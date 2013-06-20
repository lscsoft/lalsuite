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

#include <lal/LALString.h>
#include <lal/LALDetCharHvetoUtils.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

gint XLALGLibCompareSnglBurst(gconstpointer a, gconstpointer b, gpointer UNUSED user_data) {
  const SnglBurst *_a = a;
  const SnglBurst *_b = b;
  return XLALCompareSnglBurstByPeakTimeAndSNR(&_a, &_b);
}

GSnglBurstSeq* XLALCreateGSnglBurstSeq(void) {
  return g_sequence_new((GDestroyNotify)XLALDestroySnglBurst);
}

void XLALDestroyGSnglBurstSeq(GSnglBurstSeq* trig_sequence) {
  g_sequence_free(trig_sequence);
}

size_t XLALGetGSnglBurstSeqLength( GSnglBurstSeq* seq ){
  return g_sequence_get_length(seq);
}

size_t XLALGetGSnglBurstSeqUnmarkedLength( GSnglBurstSeq* trig_sequence ){
  GSnglBurstIter* itr = XLALGSnglBurstSeqBegin(trig_sequence);
  size_t cntr = 0;
  while( !g_sequence_iter_is_end(itr) ){
    SnglBurst* sb = XLALGSnglBurstIterGet(itr);
    if( sb->next == NULL ){
      cntr++;
    }
    itr = XLALGSnglBurstIterNext(itr);
  }
  return cntr;
}

GSnglBurstIter* XLALGSnglBurstSeqBegin(GSnglBurstSeq* trig_sequence) {
  return g_sequence_get_begin_iter(trig_sequence);
}

SnglBurst* XLALGSnglBurstIterGet(GSnglBurstIter* itr) {
  if(g_sequence_iter_is_end(itr)) {
    return NULL;
  }
  return (SnglBurst*)g_sequence_get(itr);
}

GSnglBurstIter* XLALGSnglBurstIterNext(GSnglBurstIter* itr) {
  return g_sequence_iter_is_end(itr) ? NULL : g_sequence_iter_next(itr);
}

GHashTable* XLALCreateGHashTable(void) {
  return g_hash_table_new_full(g_str_hash, g_str_equal, XLALFree, XLALFree);
}

void XLALDestroyGHashTable(GHashTable* tbl) {
  g_hash_table_destroy(tbl);
}

char* XLALGetStrGHashTableVal(GHashTable* tbl, const char* key) {
  return XLALStringDuplicate((char*)g_hash_table_lookup( tbl, key ));
}

double XLALGetDblGHashTableVal(GHashTable* tbl, const char* key) {
  return *(double*)g_hash_table_lookup( tbl, key );
}

size_t XLALGetIntGHashTableVal(GHashTable* tbl, const char* key) {
  return *(size_t*)g_hash_table_lookup( tbl, key );
}

char* XLALGetGHashTableKey(GHashTable* tbl, const size_t indx) {
  GList* keys = g_hash_table_get_keys( tbl );
  size_t i = 0;
  while( keys ){
    if( i == indx ) {
      return XLALStringDuplicate( (char*)keys->data );
    }
    i++;
    keys = keys->next;
  }
  return NULL;
}

GHashTable* XLALGetChannelList( GSequence *trig_sequence ){
  GSnglBurstIter* itr = XLALGSnglBurstSeqBegin(trig_sequence);
  GHashTable* channellist = XLALCreateGHashTable();
  while( !g_sequence_iter_is_end(itr) ){
    SnglBurst* sb = XLALGSnglBurstIterGet(itr);
    g_hash_table_insert( channellist, XLALStringDuplicate(sb->channel), XLALStringDuplicate(sb->channel) );
    itr = XLALGSnglBurstIterNext(itr);
  }
  return channellist;
}

GSnglBurstSeq* XLALPopulateTrigSequenceFromFile( GSnglBurstSeq* trig_sequence, const char* fname, double min_snr, char* ignore_list ){
  if( !trig_sequence ){
    trig_sequence = XLALCreateGSnglBurstSeq();
  }
  GSequence* ignorel = g_sequence_new(XLALFree);

  SnglBurst* tbl = XLALSnglBurstTableFromLIGOLw( fname );
  SnglBurst *begin = NULL, *deleteme = NULL;
  if( !tbl ) {
    return trig_sequence;
  }

  int cnt = 0;
  char* tmp;
  tmp = XLALMalloc( sizeof(char)*512 );
  FILE* lfile = fopen( ignore_list, "r" );
  while(!feof(lfile)){
    cnt = fscanf( lfile, "%s", tmp );
    if( cnt == EOF ) break;
    g_sequence_append( ignorel, XLALStringDuplicate(tmp) );
  }
  XLALFree(tmp);
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
      g_sequence_insert_sorted( trig_sequence, tbl, XLALGLibCompareSnglBurst, NULL );
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

GSnglBurstSeq* XLALPopulateTrigSequenceFromTrigList( GSnglBurstSeq* trig_sequence, SnglBurst* tbl ){
  if( !trig_sequence ){
    trig_sequence = g_sequence_new(NULL);
  }

  if( !tbl ) {
    return trig_sequence;
  }

  do {
    g_sequence_insert_sorted( trig_sequence, tbl, XLALGLibCompareSnglBurst, NULL );
    tbl=tbl->next;
  } while( tbl );

  return trig_sequence;
}
