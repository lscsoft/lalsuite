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

#include <stdint.h>
#include <unistd.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALList.h>
#include "LALValue_private.h"

struct tagLALListItem {
	struct tagLALListItem *next;
	LALValue value;
};

struct tagLALList {
	struct tagLALListItem *head;
};

/* LIST ITEM ROUTINES */

LALListItem * XLALListItemAlloc(size_t size)
{
	LALListItem *item;
	item = XLALMalloc(sizeof(*item) + size);
	if (!item)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	item->value.size = size;
	return item;
}

LALListItem * XLALListItemRealloc(LALListItem *item, size_t size)
{
	if (item == NULL)
		return XLALListItemAlloc(size);
	item = XLALRealloc(item, sizeof(*item) + size);
	if (!item)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	item->value.size = size;
	return item;
}

LALListItem * XLALListItemSet(LALListItem *item, const void *data, size_t size, LALTYPECODE type)
{
	if (XLALValueSet(&item->value, data, size, type) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return item;
}

LALListItem * XLALListItemSetValue(LALListItem *item, const LALValue *value)
{
	if (XLALValueCopy(&item->value, value) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return item;
}

LALListItem * XLALListItemDuplicate(const LALListItem *item)
{
	size_t size = sizeof(LALListItem) + item->value.size;
	LALListItem *copy = XLALMalloc(size);
	if (!copy)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	return memcpy(copy, item, size);
}

const LALValue * XLALListItemGetValue(const LALListItem *item)
{
	return &item->value;
}

LALTYPECODE XLALListItemGetValueType(const LALListItem *item)
{
	const LALValue *value = XLALListItemGetValue(item);
	if (value == NULL)
		XLAL_ERROR(XLAL_EFUNC);
	return XLALValueGetType(value);
}

void * XLALListItemGetValueData(void * data, size_t size, LALTYPECODE type, const LALListItem *item)
{
	const LALValue *value = XLALListItemGetValue(item);
	if (value == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return XLALValueGetData(data, size, type, value);
}

/* warning: shallow pointer */
const char * XLALListItemGetStringValue(const LALListItem *item)
{
	const LALValue *value = XLALListItemGetValue(item);
	if (value == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return XLALValueGetString(value);
}

#define DEFINE_GET_FUNC(TYPE, FAILVAL) \
	TYPE XLALListItemGet ## TYPE ## Value(const LALListItem *item) \
	{ \
		const LALValue *value = XLALListItemGetValue(item); \
		if (value == NULL) \
			XLAL_ERROR_VAL(FAILVAL, XLAL_EFUNC); \
		return XLALValueGet ## TYPE (value); \
	}

DEFINE_GET_FUNC(CHAR, XLAL_FAILURE)
DEFINE_GET_FUNC(INT2, XLAL_FAILURE)
DEFINE_GET_FUNC(INT4, XLAL_FAILURE)
DEFINE_GET_FUNC(INT8, XLAL_FAILURE)
DEFINE_GET_FUNC(UCHAR, XLAL_FAILURE)
DEFINE_GET_FUNC(UINT2, XLAL_FAILURE)
DEFINE_GET_FUNC(UINT4, XLAL_FAILURE)
DEFINE_GET_FUNC(UINT8, XLAL_FAILURE)
DEFINE_GET_FUNC(REAL4, XLAL_REAL4_FAIL_NAN)
DEFINE_GET_FUNC(REAL8, XLAL_REAL8_FAIL_NAN)
DEFINE_GET_FUNC(COMPLEX8, XLAL_REAL4_FAIL_NAN)
DEFINE_GET_FUNC(COMPLEX16, XLAL_REAL8_FAIL_NAN)

#undef DEFINE_GET_FUNC

REAL8 XLALListItemGetValueAsREAL8(const LALListItem *item)
{
	const LALValue *value = XLALListItemGetValue(item);
	if (value == NULL)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	return XLALValueGetAsREAL8(value);
}

/* LIST ROUTINES */

void XLALDestroyList(LALList *list)
{
	if (list != NULL) {
		LALListItem *item = list->head;
		while (item) {
			LALListItem *next = item->next;
			LALFree(item);
			item = next;
		}
		LALFree(list);
	}
	return;
}

LALList * XLALCreateList(void)
{
	return XLALCalloc(1, sizeof(LALList));
}

LALList * XLALListDuplicate(const LALList *list)
{
	LALList *newlist;
	const LALListItem *item = list->head;
	LALListItem *head = NULL;
	LALListItem *prev = NULL;

	while (item) {
		LALListItem *copy = NULL;
		copy = XLALListItemDuplicate(item);
		copy->next = NULL;
		if (prev == NULL)
			head = copy;
		else
			prev->next = copy;
		prev = copy;
		item = item->next;
	}

	newlist = XLALCreateList();
	newlist->head = head;
	return newlist;
}

int XLALListReverse(LALList *list)
{
	if (list != NULL) {
		LALListItem *item = list->head;
		LALListItem *tsil = NULL; // reversed list
		while (item) {
			LALListItem *next = item->next;
			item->next = tsil;
			tsil = item;
			item = next;
		}
		list->head = tsil;
	}
	return 0;
}

static LALListItem *XLALListSortMerge(LALListItem *a, LALListItem *b, int (*cmp)(const LALValue *, const LALValue *, void *), void *thunk)
{
	LALListItem *result = NULL;
	if (a == NULL)
		return b;
	if (b == NULL)
		return a;

	if (cmp(&a->value, &b->value, thunk) <= 0) {
		result = a;
		result->next = XLALListSortMerge(a->next, b, cmp, thunk);
	} else {
		result = b;
		result->next = XLALListSortMerge(a, b->next, cmp, thunk);
	}

	return result;
}

int XLALListSort(LALList *list, int (*cmp)(const LALValue *, const LALValue *, void *), void *thunk)
{
	LALListItem *slow;
	LALListItem *fast;
	LALList head;
	LALList tail;

	if (list->head == NULL || list->head->next == NULL)
		return 0;

	/* find midpoint */
	slow = list->head;
	fast = slow->next;
	while (fast != NULL) {
		fast = fast->next;
		if (fast != NULL) {
			slow = slow->next;
			fast = fast->next;
		}
	}

	/* split at midpoint to get back half */
	head.head = list->head;
	tail.head = slow->next;
	slow->next = NULL;

	/* recursively sort front and back */
	XLALListSort(&head, cmp, thunk);
	XLALListSort(&tail, cmp, thunk);

	/* merge sort front and back */
	list->head = XLALListSortMerge(head.head, tail.head, cmp, thunk);
	return 0;
}

size_t XLALListSize(const LALList *list)
{
	size_t size = 0;
	if (list != NULL) {
		const LALListItem *item = list->head;
		while (item) {
			item = item->next;
			++size;
		}
	}
	return size;
}

void XLALListForeach(LALList *list, void (*func)(LALValue *, void *), void *thunk)
{
	if (list != NULL) {
		LALListItem *item = list->head;
		while (item) {
			func(&item->value, thunk);
			item = item->next;
		}
	}
	return;
}

LALListItem * XLALListPop(LALList *list)
{
	LALListItem *item = NULL;
	if (list != NULL && list->head != NULL) {
		item = list->head;
		list->head = item->next;
	}
	return item;
}

LALListItem * XLALListLast(LALList *list)
{
	LALListItem *item = NULL;
	if (list != NULL && list->head != NULL) {
		item = list->head;
		while (item->next != NULL)
			item = item->next;
	}
	return item;
}

LALListItem * XLALListFind(LALList *list, int (*func)(const LALValue *, void *), void *thunk)
{
	LALListItem *item = NULL;
	if (list != NULL) {
		item = list->head;
		while (item) {
			if (func(&item->value, thunk))
				break;
			item = item->next;
		}
	}
	return item;
}

static int XLALListFindValueFunc(const LALValue *value, void *thunk)
{
	const LALValue *target = (const LALValue *)thunk;
	return XLALValueEqual(value, target);
}

LALListItem * XLALListFindValue(LALList *list, const LALValue *value)
{
	void *thunk = (void *)(uintptr_t)value; /* discard const qualifier */
	return XLALListFind(list, XLALListFindValueFunc, thunk);
}

int XLALListReplace(LALList *list, int (*func)(const LALValue *, void *), void *thunk, const LALValue *replace)
{
	size_t size = XLALValueGetSize(replace);
	if (list != NULL) {
		LALListItem *item = list->head;
		LALListItem *prev = NULL;
		while (item) {
			if (func(&item->value, thunk)) { /* found it! */
				LALListItem *orig = item;

				if (XLALValueGetSize(&item->value) != size)
					item = XLALListItemRealloc(orig, size);

				XLALListItemSetValue(item, replace);

				/* repair links if necessary */
				if (orig != item) {
					if (prev == NULL) /* head is changed */
						list->head = item;
					else
						prev->next = item;
				}

				return 0;
			}
			prev = item;
			item = item->next;
		}
	}
	return -1; /* not found */
}

int XLALListReplaceAll(LALList *list, int (*func)(const LALValue *, void *), void *thunk, const LALValue *replace)
{
	int replaced = 0;
	size_t size = XLALValueGetSize(replace);
	if (list != NULL) {
		LALListItem *item = list->head;
		LALListItem *prev = NULL;
		while (item) {
			if (func(&item->value, thunk)) { /* found it! */
				LALListItem *orig = item;

				if (XLALValueGetSize(&item->value) != size)
					item = XLALListItemRealloc(orig, size);

				XLALListItemSetValue(item, replace);

				/* repair links if necessary */
				if (orig != item) {
					if (prev == NULL) /* head is changed */
						list->head = item;
					else
						prev->next = item;
				}

				++replaced;
			}
			prev = item;
			item = item->next;
		}
	}
	return replaced;
}

int XLALListReplaceValue(LALList *list, const LALValue *value, const LALValue *replace)
{
	void *thunk = (void *)(uintptr_t)value; /* discard const qualifier */
	return XLALListReplace(list, XLALListFindValueFunc, thunk, replace);
}

int XLALListReplaceValueAll(LALList *list, const LALValue *value, const LALValue *replace)
{
	void *thunk = (void *)(uintptr_t)value; /* discard const qualifier */
	return XLALListReplaceAll(list, XLALListFindValueFunc, thunk, replace);
}

int XLALListRemove(LALList *list, int (*func)(const LALValue *, void *), void *thunk)
{
	if (list != NULL) {
		LALListItem *item = list->head;
		LALListItem *prev = NULL;
		while (item) {
			if (func(&item->value, thunk)) { /* found it! */
				if (prev == NULL) /* head is removed */
					list->head = item->next;
				else
					prev->next = item->next;
				LALFree(item);
				return 0;
			}
			prev = item;
			item = item->next;
		}
	}
	return -1; /* not found */
}

int XLALListRemoveAll(LALList *list, int (*func)(const LALValue *, void *), void *thunk)
{
	int removed = 0;
	if (list != NULL) {
		LALListItem *item = list->head;
		LALListItem *prev = NULL;
		while (item) {
			if (func(&item->value, thunk)) { /* found it! */
				LALListItem *next = item->next;
				if (prev == NULL) /* head is removed */
					list->head = next;
				else
					prev->next = next;
				LALFree(item);
				item = next;
				++removed;
			} else {
				prev = item;
				item = item->next;
			}
		}
	}
	return removed;
}

int XLALListRemoveValue(LALList *list, const LALValue *value)
{
	void *thunk = (void *)(uintptr_t)value; /* discard const qualifier */
	return XLALListRemove(list, XLALListFindValueFunc, thunk);
}

int XLALListRemoveValueAll(LALList *list, const LALValue *value)
{
	void *thunk = (void *)(uintptr_t)value; /* discard const qualifier */
	return XLALListRemoveAll(list, XLALListFindValueFunc, thunk);
}

void XLALListIterInit(LALListIter *iter, LALList *list)
{
	iter->next = list == NULL ? NULL : list->head;
}

LALListItem * XLALListIterNext(LALListIter *iter)
{
	if (iter->next == NULL) /* iteration terminated */
		return NULL;
	LALListItem *next = iter->next;
	iter->next = next->next;
	return next;
}

int XLALListAdd(LALList *list, const void *data, size_t size, LALTYPECODE type)
{
	LALListItem *item;

	XLAL_CHECK(list != NULL, XLAL_EFAULT);

	item = XLALListItemAlloc(size);
	if (item == NULL)
                XLAL_ERROR(XLAL_EFUNC);

	if (XLALListItemSet(item, data, size, type) == NULL) {
		LALFree(item);
                XLAL_ERROR(XLAL_EFUNC);
	}

	item->next = list->head;
	list->head = item;
	return 0;
}

int XLALListAddValue(LALList *list, const LALValue *value)
{
	LALTYPECODE type = XLALValueGetType(value);
	size_t size = XLALValueGetSize(value);
	const void * data = XLALValueGetDataPtr(value);
	return XLALListAdd(list, data, size, type);
}

int XLALListAddStringValue(LALList *list, const char *string)
{
	size_t size = strlen(string) + 1;
	return XLALListAdd(list, string, size, LAL_CHAR_TYPE_CODE);
}

#define DEFINE_ADD_FUNC(TYPE, TCODE) \
	int XLALListAdd ## TYPE ## Value (LALList *list, TYPE value) \
	{ return XLALListAdd(list, &value, sizeof(value), TCODE); } \

DEFINE_ADD_FUNC(CHAR, LAL_CHAR_TYPE_CODE)
DEFINE_ADD_FUNC(INT2, LAL_I2_TYPE_CODE)
DEFINE_ADD_FUNC(INT4, LAL_I4_TYPE_CODE)
DEFINE_ADD_FUNC(INT8, LAL_I8_TYPE_CODE)
DEFINE_ADD_FUNC(UCHAR, LAL_UCHAR_TYPE_CODE)
DEFINE_ADD_FUNC(UINT2, LAL_U2_TYPE_CODE)
DEFINE_ADD_FUNC(UINT4, LAL_U4_TYPE_CODE)
DEFINE_ADD_FUNC(UINT8, LAL_U8_TYPE_CODE)
DEFINE_ADD_FUNC(REAL4, LAL_S_TYPE_CODE)
DEFINE_ADD_FUNC(REAL8, LAL_D_TYPE_CODE)
DEFINE_ADD_FUNC(COMPLEX8, LAL_C_TYPE_CODE)
DEFINE_ADD_FUNC(COMPLEX16, LAL_Z_TYPE_CODE)

#undef DEFINE_ADD_FUNC
static void XLALListValuePrintFunc(LALValue *value, void *thunk)
{
	int fd = *(int *)(thunk);
	int fd2 = dup(fd);
	FILE *fp = fdopen(fd2, "w");
	XLALValuePrint(value, fd);
	fprintf(fp, "\n");
	fclose(fp);
	return;
}

void XLALListPrint(LALList *list, int fd)
{
	XLALListForeach(list, XLALListValuePrintFunc, &fd);
	return;
}
