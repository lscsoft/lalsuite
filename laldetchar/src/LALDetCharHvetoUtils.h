/*
 *  Copyright (C) 2013 Chris Pankow
 *
 *  This program is free software; ynu can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Fnundation; either version 2 of the License, or
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

#ifndef _LALDETCHARHVETOUTIL_H
#define _LALDETCHARHVETOUTIL_H

#include <lal/LALDetCharGlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/LIGOLwXMLBurstRead.h>

#ifdef  __cplusplus
extern "C" {
#endif

gint XLALGLibCompareSnglBurst(gconstpointer a, gconstpointer b, gpointer user_data);

#ifdef SWIG
SWIGLAL(EXTERNAL_STRUCT(GSequence, g_sequence_free));
#endif
typedef GSequence GSnglBurstSeq;
GSnglBurstSeq* XLALCreateGSnglBurstSeq(void);
void XLALDestroyGSnglBurstSeq(GSnglBurstSeq* trig_sequence);
size_t XLALGetGSnglBurstSeqLength( GSnglBurstSeq* seq );
size_t XLALGetGSnglBurstSeqUnmarkedLength( GSnglBurstSeq* seq );

typedef GSequenceIter GSnglBurstIter;
GSnglBurstIter* XLALGSnglBurstSeqBegin(GSnglBurstSeq* trig_sequence);
GSnglBurstIter* XLALGSnglBurstIterNext(GSnglBurstIter* iter);

#ifdef SWIG
SWIGLAL(RETURNS_PROPERTY(SnglBurst*, XLALGSnglBurstIterGet));
#endif
SnglBurst* XLALGSnglBurstIterGet(GSnglBurstIter* iter);


#ifdef SWIG
SWIGLAL(EXTERNAL_STRUCT(GHashTable, g_hash_table_destroy));
#endif
GHashTable* XLALCreateGHashTable(void);
void XLALDestroyGHashTable(GHashTable* tbl);
char* XLALGetStrGHashTableVal(GHashTable* tbl, const char* key);
double XLALGetDblGHashTableVal(GHashTable* tbl, const char* key);
size_t XLALGetIntGHashTableVal(GHashTable* tbl, const char* key);
char* XLALGetGHashTableKey(GHashTable* tbl, const size_t indx);

#if 0
typedef GHashTable GStrHashTable;
GStrHashTable* XLALCreateGStrHashTable(void);
void XLALDestroyGStrHashTable(GStrHashTable* tbl);
char* XLALGetGStrHashTableVal(GStrHashTable* tbl, const char* key);
char* XLALGetGStrHashTableKey(GStrHashTable* tbl, const size_t indx);

typedef GHashTable GIntHashTable;
GIntHashTable* XLALCreateGIntHashTable(void);
void XLALDestroyGIntHashTable(GIntHashTable* tbl);
size_t XLALGetGIntHashTableVal(GIntHashTable* tbl, const char* key);
char* XLALGetGIntHashTableKey(GIntHashTable* tbl, const size_t indx);

typedef GHashTable GDblHashTable;
GDblHashTable* XLALCreateGDblHashTable(void);
void XLALDestroyGDblHashTable(GDblHashTable* tbl);
double XLALGetGDblHashTableVal(GDblHashTable* tbl, const char* key);
char* XLALGetGDblHashTableKey(GDblHashTable* tbl, const size_t indx);
#endif

GHashTable* XLALGetChannelList( GSequence *trig_sequence );
GSnglBurstSeq* XLALPopulateTrigSequenceFromFile( GSnglBurstSeq* trig_sequence, const char* fname, double min_snr, char* ignore_list );
GSnglBurstSeq* XLALPopulateTrigSequenceFromTrigList( GSnglBurstSeq* trig_sequence, SnglBurst* tbl );

#ifdef  __cplusplus
}
#endif

#endif
