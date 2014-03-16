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
#include <lal/LALString.h>
#include <lal/LALDetCharGlib.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

int XLALCountUnmarkedSnglBurst(LALGSequence* seq) {
  XLAL_CHECK(seq, XLAL_EFAULT);
  int unmarked = 0;
  LALGSequenceIter* itr = XLALGSequenceBegin(seq);
  while( itr ){
    SnglBurst* sb = XLALGetGSeqSnglBurst(itr);
    if( sb->next == NULL ){
      unmarked++;
    }
    if (!XLALGSequenceNext(itr)) {
      XLALDestroyGSequenceIter(itr);
      itr = NULL;
    }
  }
  return unmarked;
}

LALGHashTable* XLALGetChannelList(LALGSequence *trig_sequence) {
  LALGHashTable* channellist = XLALCreateGHashTable(LALGTYPE_STR);
  LALGSequenceIter* itr = XLALGSequenceBegin(trig_sequence);
  while( itr ){
    SnglBurst* sb = XLALGetGSeqSnglBurst(itr);
    XLALSetGHashTblStr(channellist, sb->channel, sb->channel);
    if (!XLALGSequenceNext(itr)) {
      XLALDestroyGSequenceIter(itr);
      itr = NULL;
    }
  }
  return channellist;
}

LALGSequence* XLALPopulateTrigSequenceFromFile( LALGSequence* trig_sequence, const char* fname, double min_snr, char* ignore_list ){
  if( !trig_sequence ){
    trig_sequence = XLALCreateGSequence(LALGTYPE_SNGL_BURST);
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
      XLALAddGSeqSnglBurst( trig_sequence, tbl );
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

  LALGSequenceIter* itr = XLALGSequenceBegin(trig_sequence);
  while( itr ){
    SnglBurst* sb = XLALGetGSeqSnglBurst(itr);
    sb->next = NULL;
    if (!XLALGSequenceNext(itr)) {
      XLALDestroyGSequenceIter(itr);
      itr = NULL;
    }
  }

  return trig_sequence;
}

LALGSequence* XLALPopulateTrigSequenceFromTrigList( LALGSequence* trig_sequence, SnglBurst* tbl ){
  if( !trig_sequence ){
    trig_sequence = XLALCreateGSequence(LALGTYPE_SNGL_BURST);
  }

  if( !tbl ) {
    return trig_sequence;
  }

  do {
    LALGSequenceIter* itr = XLALAddGSeqSnglBurst( trig_sequence, tbl );
    tbl=tbl->next;
    SnglBurst* sb = XLALGetGSeqSnglBurst(itr);
    sb->next = NULL;
    XLALDestroyGSequenceIter(itr);
  } while( tbl );

  return trig_sequence;
}
