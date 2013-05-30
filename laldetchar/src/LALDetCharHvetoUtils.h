/*
 *  Copyright (C) 2013 Chris Pankow
 *
 e  This program is free software; ynu can redistribute it and/or modify
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

/*
 * These functions are mostly for the convenience of SWIG wrapping them and 
 * using them in python based applications.
 */

#ifndef _LALDETCHARHVETOUTIL_H
#define _LALDETCHARHVETOUTIL_H

#include <lal/LALDetCharGlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/LIGOLwXMLBurstRead.h>

static gint compare(gconstpointer a, gconstpointer b) {
        const SnglBurst *_a = (const SnglBurst *)a;
        const SnglBurst *_b = (const SnglBurst *)b;

        return XLALCompareSnglBurstByPeakTimeAndSNR(&_a, &_b);
}

//double get_double_from_ghash( const char* key, GHashTable* tbl );
double XLALGLibDoubleFromGHash( const char* key, GHashTable* tbl );
//size_t get_size_t_from_ghash( const char* key, GHashTable* tbl );
size_t XLALGLibSizeTFromGHash( const char* key, GHashTable* tbl );
//char* get_string_from_ghash( const char* key, GHashTable* tbl );
char* XLALGLibStringFromGHash( const char* key, GHashTable* tbl );
//char* get_key_from_ghash( GHashTable* tbl, size_t indx );
char* XLALGLibKeyFromGHash( GHashTable* tbl, size_t indx );
//size_t gseq_get_length( GSequence* seq );
size_t XLALGLibSeqGetLength( GSequence* seq );
//size_t gseq_get_unmarked_length( GSequence* seq );
size_t XLALGLibSeqGetUnmarkedLength( GSequence* seq );
//void gseq_delete( GSequence* seq );
void XLALGLibSeqDelete( GSequence* seq );
GHashTable* XLALGetChannelList( GSequence *triglist );
GHashTable* XLALCreateGHashTable( void );
// NOTE: This has to be the full struct name or else the attributes of the 
// SnglBurst won't be exposed properly
// TODO: Check if this is the case with the new lalmetaio bindings
//struct tagSnglBurst* list_from_g_sequence( GSequence* in );
struct tagSnglBurst* XLALGLibSBListFromSeq( GSequence* in );
GSequence* XLALPopulateTrigSequenceFromFile( GSequence* trig_sequence, const char* fname, double min_snr, char* ignore_list );
GSequence* XLALPopulateTrigSequenceFromTrigList( GSequence* trig_sequence, SnglBurst* tbl );

#endif
