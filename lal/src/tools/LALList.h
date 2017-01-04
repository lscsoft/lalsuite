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

#ifndef _LAL_LIST_H
#define _LAL_LIST_H

#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALValue.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

struct tagLALListItem;
typedef struct tagLALListItem LALListItem;

struct tagLALList;
typedef struct tagLALList LALList;

struct tagLALListIter {
	/* private data */
	struct tagLALListItem *next;
};
typedef struct tagLALListIter LALListIter;

LALListItem *XLALListItemAlloc(size_t size);
LALListItem *XLALListItemRealloc(LALListItem *item, size_t size);
LALListItem *XLALListItemSet(LALListItem *item, const void *data, size_t size, LALTYPECODE type);
LALListItem *XLALListItemSetValue(LALListItem *item, const LALValue *value);
LALListItem *XLALListItemDuplicate(const LALListItem *item);

const LALValue * XLALListItemGetValue(const LALListItem *item);
LALTYPECODE XLALListItemGetValueType(const LALListItem *item);
void * XLALListItemGetValueData(void * data, size_t size, LALTYPECODE type, const LALListItem *item);
const char * XLALListItemGetStringValue(const LALListItem *item);
CHAR XLALListItemGetCHARValue(const LALListItem *item);
INT2 XLALListItemGetINT2Value(const LALListItem *item);
INT4 XLALListItemGetINT4Value(const LALListItem *item);
INT8 XLALListItemGetINT8Value(const LALListItem *item);
UCHAR XLALListItemGetUCHARValue(const LALListItem *item);
UINT2 XLALListItemGetUINT2Value(const LALListItem *item);
UINT4 XLALListItemGetUINT4Value(const LALListItem *item);
UINT8 XLALListItemGetUINT8Value(const LALListItem *item);
REAL4 XLALListItemGetREAL4Value(const LALListItem *item);
REAL8 XLALListItemGetREAL8Value(const LALListItem *item);
COMPLEX8 XLALListItemGetCOMPLEX8Value(const LALListItem *item);
COMPLEX16 XLALListItemGetCOMPLEX16Value(const LALListItem *item);

REAL8 XLALListItemGetValueAsREAL8(const LALListItem *item);

void XLALDestroyList(LALList *list);
LALList * XLALCreateList(void);

LALList * XLALListDuplicate(const LALList *list);
int XLALListReverse(LALList *list);
int XLALListSort(LALList *list, int (*cmp)(const LALValue *, const LALValue *, void *), void *thunk);
size_t XLALListSize(const LALList *list);
void XLALListForeach(LALList *list, void (*func)(LALValue *, void *), void *thunk);
LALListItem * XLALListPop(LALList *list);
LALListItem * XLALListLast(LALList *list);
LALListItem * XLALListFind(LALList *list, int (*func)(const LALValue *, void *), void *thunk);
LALListItem * XLALListFindValue(LALList *list, const LALValue *value);
int XLALListReplace(LALList *list, int (*func)(const LALValue *, void *), void *thunk, const LALValue *replace);
int XLALListReplaceAll(LALList *list, int (*func)(const LALValue *, void *), void *thunk, const LALValue *replace);
int XLALListReplaceValue(LALList *list, const LALValue *value, const LALValue *replace);
int XLALListReplaceValueAll(LALList *list, const LALValue *value, const LALValue *replace);
int XLALListRemove(LALList *list, int (*func)(const LALValue *, void *), void *thunk);
int XLALListRemoveAll(LALList *list, int (*func)(const LALValue *, void *), void *thunk);
int XLALListRemoveValue(LALList *list, const LALValue *value);
int XLALListRemoveValueAll(LALList *list, const LALValue *value);
void XLALListIterInit(LALListIter *iter, LALList *list);
LALListItem * XLALListIterNext(LALListIter *iter);

int XLALListAdd(LALList *list, const void *data, size_t size, LALTYPECODE type);
int XLALListAddValue(LALList *list, const LALValue *value);
int XLALListAddStringValue(LALList *list, const char *value);
int XLALListAddCHARValue(LALList *list, CHAR value);
int XLALListAddINT2Value(LALList *list, INT2 value);
int XLALListAddINT4Value(LALList *list, INT4 value);
int XLALListAddINT8Value(LALList *list, INT8 value);
int XLALListAddUCHARValue(LALList *list, UCHAR value);
int XLALListAddUINT2Value(LALList *list, UINT2 value);
int XLALListAddUINT4Value(LALList *list, UINT4 value);
int XLALListAddUINT8Value(LALList *list, UINT8 value);
int XLALListAddREAL4Value(LALList *list, REAL4 value);
int XLALListAddREAL8Value(LALList *list, REAL8 value);
int XLALListAddCOMPLEX8Value(LALList *list, COMPLEX8 value);
int XLALListAddCOMPLEX16Value(LALList *list, COMPLEX16 value);

void XLALListPrint(LALList *list, int fd);


#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LAL_LIST_H */
