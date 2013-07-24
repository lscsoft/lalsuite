/*
 *  Copyright (C) 2013 Chris Pankow, Karl Wette
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

#ifndef _LALDETCHARGLIBTYPES_H
#define _LALDETCHARGLIBTYPES_H

#include <stdbool.h>
#include <lal/LALDetCharGlib.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum {
  LALGTYPE_NONE = 0,
  LALGTYPE_INT,
  LALGTYPE_DBL,
  LALGTYPE_STR,
  LALGTYPE_SNGL_BURST,
  LALGTYOE_LAST
} LALGType;

typedef struct tagLALGSequence LALGSequence;
typedef struct tagLALGSequenceIter LALGSequenceIter;
typedef struct tagLALGHashTable LALGHashTable;

LALGSequence* XLALCreateGSequence(LALGType type);
void XLALDestroyGSequence(LALGSequence* seq);
LALGType XLALGetGSequenceType(LALGSequence* seq);
size_t XLALGetGSequenceLength(LALGSequence* seq);

void XLALDestroyGSequenceIter(LALGSequenceIter* itr);
LALGSequenceIter* XLALGSequenceBegin(LALGSequence* seq);
LALGSequenceIter* XLALGSequenceNext(LALGSequenceIter* itr);
#ifndef SWIG   // exclude from SWIG interface
GSequenceIter* XLALGSequenceBeginRaw(LALGSequence* seq, LALGType type);
#endif

#ifdef SWIG   // SWIG interface directives
SWIGLAL(RETURNS_PROPERTY(SnglBurst*, XLALGetGSeqSnglBurst));
#endif
SnglBurst* XLALGetGSeqSnglBurst(LALGSequenceIter* itr);
#ifdef SWIG   // SWIG interface directives
SWIGLAL(ACQUIRES_OWNERSHIP(SnglBurst*, sb));
#endif
void XLALAddGSeqSnglBurst(LALGSequence* seq, SnglBurst* sb);
#ifdef SWIG   // SWIG interface directives
SWIGLAL_CLEAR(ACQUIRES_OWNERSHIP(SnglBurst*, sb));
#endif

LALGHashTable* XLALCreateGHashTable(LALGType type);
void XLALDestroyGHashTable(LALGHashTable* tbl);
LALGType XLALGetGHashTableType(LALGHashTable* tbl);
bool XLALGHashTableKeyExists(LALGHashTable* tbl, const char* key);
char* XLALGetGHashTableKey(LALGHashTable* tbl, const size_t indx);
#ifndef SWIG   // exclude from SWIG interface
void XLALGHashTableBeginRaw(LALGHashTable* tbl, LALGType type, GHashTableIter* itr);
#endif

void XLALSetGHashTblInt(LALGHashTable* tbl, const char* key, const int val);
int XLALGetGHashTblInt(LALGHashTable* tbl, const char* key);
#ifndef SWIG   // exclude from SWIG interface
int* XLALGetGHashTblIntPtr(LALGHashTable* tbl, const char* key);
#endif

void XLALSetGHashTblDbl(LALGHashTable* tbl, const char* key, const double val);
double XLALGetGHashTblDbl(LALGHashTable* tbl, const char* key);
#ifndef SWIG   // exclude from SWIG interface
double* XLALGetGHashTblDblPtr(LALGHashTable* tbl, const char* key);
#endif

void XLALSetGHashTblStr(LALGHashTable* tbl, const char* key, const char* val);
char* XLALGetGHashTblStr(LALGHashTable* tbl, const char* key);

#ifdef  __cplusplus
}
#endif

#endif /* _LALDETCHARGLIBTYPES_H */
