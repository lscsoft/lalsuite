/*
*  Copyright (C) 2016 Jolien Creighton
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

#ifndef _LAL_DICT_H
#define _LAL_DICT_H

#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALValue.h>
#include <lal/LALList.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

#define LAL_KEYNAME_MAX 31

struct tagLALDictEntry;
typedef struct tagLALDictEntry LALDictEntry;

struct tagLALDict;
typedef struct tagLALDict LALDict;

struct tagLALDictIter {
	/* private data */
	struct tagLALDict *dict;
	struct tagLALDictEntry *next;
	size_t pos;
};
typedef struct tagLALDictIter LALDictIter;

void XLALDictEntryFree(LALDictEntry *list);
LALDictEntry * XLALDictEntryAlloc(size_t size);
LALDictEntry * XLALDictEntryRealloc(LALDictEntry *entry, size_t size);
LALDictEntry * XLALDictEntrySetKey(LALDictEntry *entry, const char *key);
LALDictEntry * XLALDictEntrySetValue(LALDictEntry *entry, const void *data, size_t size, LALTYPECODE type);

/* warning: shallow pointer */
const char * XLALDictEntryGetKey(const LALDictEntry *entry);
/* warning: shallow pointer */
const LALValue * XLALDictEntryGetValue(const LALDictEntry *entry);

void XLALDestroyDict(LALDict *dict);
LALDict * XLALCreateDict(void);

void XLALDictForeach(LALDict *dict, void (*func)(char *, LALValue *, void *), void *thunk);
LALDictEntry * XLALDictFind(LALDict *dict, int (*func)(const char *, const LALValue *, void *), void *thunk);
void XLALDictIterInit(LALDictIter *iter, LALDict *dict);
LALDictEntry * XLALDictIterNext(LALDictIter *iter);

LALList * XLALDictKeys(const LALDict *dict);
LALList * XLALDictValues(const LALDict *dict);

int XLALDictContains(const LALDict *dict, const char *key);
size_t XLALDictSize(const LALDict *dict);
int XLALDictRemove(LALDict *dict, const char *key);
int XLALDictInsert(LALDict *dict, const char *key, const void *data, size_t size, LALTYPECODE type);
int XLALDictInsertValue(LALDict *dict, const char *key, const LALValue *value);
int XLALDictInsertStringValue(LALDict *dict, const char *key, const char *value);
int XLALDictInsertCHARValue(LALDict *dict, const char *key, CHAR value);
int XLALDictInsertINT2Value(LALDict *dict, const char *key, INT2 value);
int XLALDictInsertINT4Value(LALDict *dict, const char *key, INT4 value);
int XLALDictInsertINT8Value(LALDict *dict, const char *key, INT8 value);
int XLALDictInsertUCHARValue(LALDict *dict, const char *key, UCHAR value);
int XLALDictInsertUINT2Value(LALDict *dict, const char *key, UINT2 value);
int XLALDictInsertUINT4Value(LALDict *dict, const char *key, UINT4 value);
int XLALDictInsertUINT8Value(LALDict *dict, const char *key, UINT8 value);
int XLALDictInsertREAL4Value(LALDict *dict, const char *key, REAL4 value);
int XLALDictInsertREAL8Value(LALDict *dict, const char *key, REAL8 value);
int XLALDictInsertCOMPLEX8Value(LALDict *dict, const char *key, COMPLEX8 value);
int XLALDictInsertCOMPLEX16Value(LALDict *dict, const char *key, COMPLEX16 value);

LALDictEntry *XLALDictLookup(LALDict *dict, const char *key);
/* warning: shallow pointer */
const char * XLALDictLookupStringValue(LALDict *dict, const char *key);
CHAR XLALDictLookupCHARValue(LALDict *dict, const char *key);
INT2 XLALDictLookupINT2Value(LALDict *dict, const char *key);
INT4 XLALDictLookupINT4Value(LALDict *dict, const char *key);
INT8 XLALDictLookupINT8Value(LALDict *dict, const char *key);
UCHAR XLALDictLookupUCHARValue(LALDict *dict, const char *key);
UINT2 XLALDictLookupUINT2Value(LALDict *dict, const char *key);
UINT4 XLALDictLookupUINT4Value(LALDict *dict, const char *key);
UINT8 XLALDictLookupUINT8Value(LALDict *dict, const char *key);
REAL4 XLALDictLookupREAL4Value(LALDict *dict, const char *key);
REAL8 XLALDictLookupREAL8Value(LALDict *dict, const char *key);
COMPLEX8 XLALDictLookupCOMPLEX8Value(LALDict *dict, const char *key);
COMPLEX16 XLALDictLookupCOMPLEX16Value(LALDict *dict, const char *key);

REAL8 XLALDictLookupValueAsREAL8(LALDict *dict, const char *key);

void XLALDictPrint(LALDict *dict, int fd);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LAL_DICT_H */
