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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/*
 * Dictionary is implemented as a hash search with algorithm adopted
 * from "The C Programming Language" by Kernighan and Ritchie, 2nd ed.
 * section 6.6.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/LALDict.h>
#include "LALValue_private.h"
#include "config.h"

#define LAL_DICT_HASHSIZE 101

struct tagLALDictEntry {
	struct tagLALDictEntry *next;
	char *key;
	LALValue value;
};

struct tagLALDict {
	size_t size;
	struct tagLALDictEntry *hashes[];
};

static size_t hash(const char *s)
{
	size_t hashval;
	for (hashval = 0; *s != '\0'; ++s)
		hashval = *s + 31 * hashval;
	return hashval;
}

/* DICT ENTRY ROUTINES */

void XLALDictEntryFree(LALDictEntry *list)
{
	while (list) {
		LALDictEntry *next = list->next;
		if (list->key)
			LALFree(list->key);
		LALFree(list);
		list = next;
	}
	return;
}

LALDictEntry * XLALDictEntryAlloc(size_t size)
{
	LALDictEntry *entry;
	entry = XLALMalloc(sizeof(*entry) + size);
	if (!entry)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	entry->key = NULL;
	entry->value.size = size;
	return entry;
}

LALDictEntry * XLALDictEntryRealloc(LALDictEntry *entry, size_t size)
{
	if (entry == NULL)
		return XLALDictEntryAlloc(size);
	if (entry->value.size == size)
		return entry;
	entry = XLALRealloc(entry, sizeof(*entry) + size);
	if (!entry)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	entry->value.size = size;
	return entry;
}

LALDictEntry * XLALDictEntrySetKey(LALDictEntry *entry, const char *key)
{
	if (entry->key)
		LALFree(entry->key);
	if ((entry->key = XLALStringDuplicate(key)) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return entry;
}

LALDictEntry * XLALDictEntrySetValue(LALDictEntry *entry, const void *data, size_t size, LALTYPECODE type)
{
	if (XLALValueSet(&entry->value, data, size, type) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return entry;
}

/* warning: shallow pointer */
const char * XLALDictEntryGetKey(const LALDictEntry *entry)
{
	return entry->key;
}

/* warning: shallow pointer */
const LALValue * XLALDictEntryGetValue(const LALDictEntry *entry)
{
	return &entry->value;
}

/* DICT ROUTINES */

void XLALClearDict(LALDict *dict)
{
	if (dict) {
		for (size_t i = 0; i < dict->size; ++i)
			XLALDictEntryFree(dict->hashes[i]);
		memset(dict, 0, sizeof(dict) + LAL_DICT_HASHSIZE * sizeof(*dict->hashes));
		dict->size = LAL_DICT_HASHSIZE;
	}
}

void XLALDestroyDict(LALDict *dict)
{
	if (dict) {
		size_t i;
		for (i = 0; i < dict->size; ++i)
			XLALDictEntryFree(dict->hashes[i]);
		LALFree(dict);
	}
	return;
}

LALDict * XLALCreateDict(void)
{
	LALDict *dict;
	dict = XLALCalloc(1, sizeof(dict) + LAL_DICT_HASHSIZE * sizeof(*dict->hashes));
	if (!dict)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	dict->size = LAL_DICT_HASHSIZE;
	return dict;
}

void XLALDictForeach(LALDict *dict, void (*func)(char *, LALValue *, void *), void *thunk)
{
	size_t i;
	for (i = 0; i < dict->size; ++i) {
		LALDictEntry *entry;
		for (entry = dict->hashes[i]; entry != NULL; entry = entry->next)
			func(entry->key, &entry->value, thunk);
	}
	return;
}

LALDictEntry * XLALDictFind(LALDict *dict, int (*func)(const char *, const LALValue *, void *), void *thunk)
{
	size_t i;
	for (i = 0; i < dict->size; ++i) {
		LALDictEntry *entry;
		for (entry = dict->hashes[i]; entry != NULL; entry = entry->next)
			if (func(entry->key, &entry->value, thunk))
				return entry;
	}
	return NULL;
}

void XLALDictIterInit(LALDictIter *iter, LALDict *dict)
{
	iter->dict = dict;
	iter->pos = 0;
	iter->next = NULL;
	return;
}

LALDictEntry * XLALDictIterNext(LALDictIter *iter)
{
	while (1) {

		if (iter->next) {
			LALDictEntry *entry = iter->next;
			iter->next = entry->next;
			return entry;
		}

		/* check end of iteration */
		if (iter->pos >= iter->dict->size)
			return NULL;

		iter->next = iter->dict->hashes[iter->pos++];
	}
	return NULL;
}

int XLALDictUpdate(LALDict *dst, const LALDict *src)
{
	size_t i;
	XLAL_CHECK(dst, XLAL_EFAULT);
	XLAL_CHECK(src, XLAL_EFAULT);
	for (i = 0; i < src->size; ++i) {
		const LALDictEntry *entry;
		for (entry = src->hashes[i]; entry != NULL; entry = entry->next) {
			const char *key = XLALDictEntryGetKey(entry);
			const LALValue *value = XLALDictEntryGetValue(entry);
			if (XLALDictInsertValue(dst, key, value) < 0)
				XLAL_ERROR(XLAL_EFUNC);
		}
	}
	return XLAL_SUCCESS;
}

LALDict * XLALDictMerge(const LALDict *dict1, const LALDict *dict2)
{
	LALDict *new = XLALCreateDict();
	XLAL_CHECK_NULL(new, XLAL_EFUNC);
	if (dict1)
		XLAL_CHECK_FAIL(XLALDictUpdate(new, dict1) == XLAL_SUCCESS, XLAL_EFUNC);
	if (dict2)
		XLAL_CHECK_FAIL(XLALDictUpdate(new, dict2) == XLAL_SUCCESS, XLAL_EFUNC);
	return new;
XLAL_FAIL:
	XLALDestroyDict(new);
	return NULL;
}

LALDict * XLALDictDuplicate(const LALDict *orig)
{
	return XLALDictMerge(orig, NULL);
}

LALList * XLALDictKeys(const LALDict *dict)
{
	LALList *list;
	size_t i;
	list = XLALCreateList();
	if (!list)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	for (i = 0; i < dict->size; ++i) {
		const LALDictEntry *entry;
		for (entry = dict->hashes[i]; entry != NULL; entry = entry->next) {
			const char *key = XLALDictEntryGetKey(entry);
			if (XLALListAddStringValue(list, key) < 0) {
				XLALDestroyList(list);
				XLAL_ERROR_NULL(XLAL_EFUNC);
			}
		}
	}
	return list;
}

LALList * XLALDictValues(const LALDict *dict)
{
	LALList *list;
	size_t i;
	list = XLALCreateList();
	if (!list)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	for (i = 0; i < dict->size; ++i) {
		const LALDictEntry *entry;
		for (entry = dict->hashes[i]; entry != NULL; entry = entry->next) {
			const LALValue *value = XLALDictEntryGetValue(entry);
			if (XLALListAddValue(list, value) < 0) {
				XLALDestroyList(list);
				XLAL_ERROR_NULL(XLAL_EFUNC);
			}
		}
	}
	return list;
}

int XLALDictContains(const LALDict *dict, const char *key)
{
	const LALDictEntry *entry;
	for (entry = dict->hashes[hash(key) % dict->size]; entry != NULL; entry = entry->next)
		if (strcmp(key, entry->key) == 0)
			return 1;
	return 0;
}

size_t XLALDictSize(const LALDict *dict)
{
	size_t size = 0;
	size_t i;
	for (i = 0; i < dict->size; ++i) {
		const LALDictEntry *entry;
		for (entry = dict->hashes[i]; entry != NULL; entry = entry->next)
			++size;
	}
	return size;
}

LALDictEntry *XLALDictLookup(const LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	for (entry = dict->hashes[hash(key) % dict->size]; entry != NULL; entry = entry->next)
		if (strcmp(key, entry->key) == 0)
			return entry;
	return NULL;
}

LALDictEntry *XLALDictPop(LALDict *dict, const char *key)
{
	size_t hashidx = hash(key) % dict->size;
	LALDictEntry *this = dict->hashes[hashidx];
	LALDictEntry *prev = this;
	while (this) {
		if (strcmp(this->key, key) == 0) { /* found it! */
			if (prev == this) /* head is removed */
				dict->hashes[hashidx] = this->next;
			else
				prev->next = this->next;
			this->next = NULL;
			return this;
		}
		prev = this;
		this = this->next;
	}
	/* not found */
	XLAL_ERROR_NULL(XLAL_ENAME, "Key `%s' not found", key);
}

int XLALDictRemove(LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	entry = XLALDictPop(dict, key);
	XLAL_CHECK(entry, XLAL_EFUNC); /* not found */
	XLALDictEntryFree(entry);
	return XLAL_SUCCESS;
}

int XLALDictInsert(LALDict *dict, const char *key, const void *data, size_t size, LALTYPECODE type)
{
	size_t hashidx = hash(key) % dict->size;
	LALDictEntry *this = dict->hashes[hashidx];
	LALDictEntry *prev = NULL;
	LALDictEntry *entry;

	/* see if entry already exists */
	while (this) {
		if (strcmp(this->key, key) == 0) { /* found it! */
			entry = XLALDictEntryRealloc(this, size);
			if (entry == NULL)
				XLAL_ERROR(XLAL_EFUNC);
			if (entry != this) { /* relink */
				if (prev == NULL) /* head is moved */
					dict->hashes[hashidx] = entry;
				else
					prev->next = entry;
			}
			entry = XLALDictEntrySetValue(entry, data, size, type);

			if (entry == NULL)
				XLAL_ERROR(XLAL_EFUNC);

			return 0;
		}
		prev = this;
		this = this->next;
	}

	/* not found: create new entry */
	entry = XLALDictEntryAlloc(size);
	if (entry == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	entry = XLALDictEntrySetKey(entry, key);
	if (entry == NULL) {
		LALFree(entry);
		XLAL_ERROR(XLAL_EFUNC);
	}

	entry = XLALDictEntrySetValue(entry, data, size, type);
	if (entry == NULL) {
		LALFree(entry);
		XLAL_ERROR(XLAL_EFUNC);
	}

	entry->next = dict->hashes[hashidx];
	dict->hashes[hashidx] = entry;
	return 0;
}

int XLALDictInsertValue(LALDict *dict, const char *key, const LALValue *value)
{
	LALTYPECODE type = XLALValueGetType(value);
	size_t size = XLALValueGetSize(value);
	const void * data = XLALValueGetDataPtr(value);
	return XLALDictInsert(dict, key, data, size, type);
}

int XLALDictInsertBLOBValue(LALDict *dict, const char *key, const void *blob, size_t size)
{
	if (XLALDictInsert(dict, key, blob, size, LAL_UCHAR_TYPE_CODE) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
}

int XLALDictInsertStringValue(LALDict *dict, const char *key, const char *string)
{
	size_t size = strlen(string) + 1;
	if (XLALDictInsert(dict, key, string, size, LAL_CHAR_TYPE_CODE) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
}

#define DEFINE_INSERT_FUNC(TYPE, TCODE) \
	int XLALDictInsert ## TYPE ## Value(LALDict *dict, const char *key, TYPE value) \
	{ \
		if (XLALDictInsert(dict, key, &value, sizeof(value), TCODE) < 0) \
			XLAL_ERROR(XLAL_EFUNC); \
		return 0; \
	}

DEFINE_INSERT_FUNC(CHAR, LAL_CHAR_TYPE_CODE)
DEFINE_INSERT_FUNC(INT2, LAL_I2_TYPE_CODE)
DEFINE_INSERT_FUNC(INT4, LAL_I4_TYPE_CODE)
DEFINE_INSERT_FUNC(INT8, LAL_I8_TYPE_CODE)
DEFINE_INSERT_FUNC(UCHAR, LAL_UCHAR_TYPE_CODE)
DEFINE_INSERT_FUNC(UINT2, LAL_U2_TYPE_CODE)
DEFINE_INSERT_FUNC(UINT4, LAL_U4_TYPE_CODE)
DEFINE_INSERT_FUNC(UINT8, LAL_U8_TYPE_CODE)
DEFINE_INSERT_FUNC(REAL4, LAL_S_TYPE_CODE)
DEFINE_INSERT_FUNC(REAL8, LAL_D_TYPE_CODE)
DEFINE_INSERT_FUNC(COMPLEX8, LAL_C_TYPE_CODE)
DEFINE_INSERT_FUNC(COMPLEX16, LAL_Z_TYPE_CODE)

#undef DEFINE_INSERT_FUNC

void * XLALDictLookupBLOBValue(const LALDict *dict, const char *key)
{
	LALDictEntry *entry = XLALDictLookup(dict, key);
	const LALValue *value;
	if (entry == NULL)
		XLAL_ERROR_NULL(XLAL_ENAME, "Key `%s' not found", key);
	value = XLALDictEntryGetValue(entry);
	if (value == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return XLALValueGetBLOB(value);
}

/* warning: shallow pointer */
const char * XLALDictLookupStringValue(const LALDict *dict, const char *key)
{
	LALDictEntry *entry = XLALDictLookup(dict, key);
	const LALValue *value;
	if (entry == NULL)
		XLAL_ERROR_NULL(XLAL_ENAME, "Key `%s' not found", key);
	value = XLALDictEntryGetValue(entry);
	if (value == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return XLALValueGetString(value);
}

#define DEFINE_LOOKUP_FUNC(TYPE, FAILVAL) \
	TYPE XLALDictLookup ## TYPE ## Value(const LALDict *dict, const char *key) \
	{ \
		LALDictEntry *entry; \
		const LALValue *value; \
		entry = XLALDictLookup(dict, key); \
		if (entry == NULL) \
			XLAL_ERROR_VAL(FAILVAL, XLAL_ENAME, "Key `%s' not found", key); \
		value = XLALDictEntryGetValue(entry); \
		if (value == NULL) \
			XLAL_ERROR_VAL(FAILVAL, XLAL_EFUNC); \
		return XLALValueGet ## TYPE (value); \
	}

DEFINE_LOOKUP_FUNC(CHAR, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(INT2, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(INT4, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(INT8, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(UCHAR, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(UINT2, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(UINT4, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(UINT8, XLAL_FAILURE)
DEFINE_LOOKUP_FUNC(REAL4, XLAL_REAL4_FAIL_NAN)
DEFINE_LOOKUP_FUNC(REAL8, XLAL_REAL8_FAIL_NAN)
DEFINE_LOOKUP_FUNC(COMPLEX8, XLAL_REAL4_FAIL_NAN)
DEFINE_LOOKUP_FUNC(COMPLEX16, XLAL_REAL8_FAIL_NAN)

REAL8 XLALDictLookupValueAsREAL8(const LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	const LALValue *value;
	entry = XLALDictLookup(dict, key);
	if (entry == NULL)
		XLAL_ERROR_REAL8(XLAL_ENAME, "Key `%s' not found", key);
	value = XLALDictEntryGetValue(entry);
	if (value == NULL)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	return XLALValueGetAsREAL8(value);
}

LALValue * XLALDictPopValue(LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	LALValue *value;
	entry = XLALDictPop(dict, key);
	XLAL_CHECK_NULL(entry, XLAL_EFUNC);
	value = XLALValueDuplicate(XLALDictEntryGetValue(entry));
	XLALDictEntryFree(entry);
	return value;
}

void * XLALDictPopBLOBValue(LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	void *value;
	entry = XLALDictPop(dict, key);
	XLAL_CHECK_NULL(entry, XLAL_EFUNC);
	value = XLALValueGetBLOB(XLALDictEntryGetValue(entry));
	XLALDictEntryFree(entry);
	return value;
}

char * XLALDictPopStringValue(LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	char *value;
	entry = XLALDictPop(dict, key);
	XLAL_CHECK_NULL(entry, XLAL_EFUNC);
	value = XLALStringDuplicate(XLALValueGetString(XLALDictEntryGetValue(entry)));
	XLALDictEntryFree(entry);
	XLAL_CHECK_NULL(value, XLAL_EFUNC);
	return value;
}

#define DEFINE_POP_FUNC(TYPE, FAILVAL) \
	TYPE XLALDictPop ## TYPE ## Value(LALDict *dict, const char *key) \
	{ \
		LALDictEntry *entry; \
		TYPE value; \
		entry = XLALDictPop(dict, key); \
		XLAL_CHECK_VAL(FAILVAL, entry, XLAL_EFUNC); \
		value = XLALValueGet ## TYPE(XLALDictEntryGetValue(entry)); \
		XLALDictEntryFree(entry); \
		return value; \
	}

DEFINE_POP_FUNC(CHAR, XLAL_FAILURE)
DEFINE_POP_FUNC(INT2, XLAL_FAILURE)
DEFINE_POP_FUNC(INT4, XLAL_FAILURE)
DEFINE_POP_FUNC(INT8, XLAL_FAILURE)
DEFINE_POP_FUNC(UCHAR, XLAL_FAILURE)
DEFINE_POP_FUNC(UINT2, XLAL_FAILURE)
DEFINE_POP_FUNC(UINT4, XLAL_FAILURE)
DEFINE_POP_FUNC(UINT8, XLAL_FAILURE)
DEFINE_POP_FUNC(REAL4, XLAL_REAL4_FAIL_NAN)
DEFINE_POP_FUNC(REAL8, XLAL_REAL8_FAIL_NAN)
DEFINE_POP_FUNC(COMPLEX8, XLAL_REAL4_FAIL_NAN)
DEFINE_POP_FUNC(COMPLEX16, XLAL_REAL8_FAIL_NAN)

REAL8 XLALDictPopValueAsREAL8(LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	REAL8 value;
	entry = XLALDictPop(dict, key);
	XLAL_CHECK_REAL8(entry, XLAL_EFUNC);
	value = XLALValueGetAsREAL8(XLALDictEntryGetValue(entry));
	XLALDictEntryFree(entry);
	return value;
}

struct LALDictAsStringAppendValueFuncParams {char *s; int first;};

static void XLALDictAsStringAppendValueFunc(char *key, LALValue *value, void *thunk)
{
	struct LALDictAsStringAppendValueFuncParams *p = (struct LALDictAsStringAppendValueFuncParams *)(thunk);
	if (p->first)
		p->first = 0;
	else
		p->s = XLALStringAppend(p->s, ", ");
	p->s = XLALStringAppendFmt(p->s, "\"%s\": ", key);
	p->s = XLALValueAsStringAppend(p->s, value);
	return;
}

char * XLALDictAsStringAppend(char *s, const LALDict *dict)
{
	LALDict *d = (LALDict *)(uintptr_t)dict; // discarding const qual is harmless
	struct LALDictAsStringAppendValueFuncParams p = {s, 1};
	p.s = XLALStringAppend(p.s, "{");
	XLALDictForeach(d, XLALDictAsStringAppendValueFunc, &p);
	p.s = XLALStringAppend(p.s, "}");
	return p.s;
}

void XLALDictPrint(const LALDict *dict, int fd)
{
	char *s = NULL;
	s = XLALDictAsStringAppend(s, dict);
	XLAL_CHECK_VOID(s, XLAL_EFUNC);
#if HAVE_DPRINTF
	dprintf(fd, "%s", s);
#else
	/* hack... */
	switch (fd) {
	case 1:
		fprintf(stdout, "%s", s);
		break;
	case 2:
		fprintf(stderr, "%s", s);
		break;
	default:
		LALFree(s);
		XLAL_ERROR_VOID(XLAL_EIO, "Don't know what to do with file descriptor %d", fd);
	}
#endif
	LALFree(s);
	return;
}
