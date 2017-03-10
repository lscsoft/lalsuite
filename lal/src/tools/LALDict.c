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
#include <lal/LALDict.h>
#include "LALValue_private.h"

#define LAL_DICT_HASHSIZE 101

struct tagLALDictEntry {
        struct tagLALDictEntry *next;
        char key[LAL_KEYNAME_MAX + 1];
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
	if ((size_t)snprintf(entry->key, sizeof(entry->key), "%s", key) >= sizeof(entry->key))
		XLAL_ERROR_NULL(XLAL_ENAME, "Key name `%s' too long (max %d characters)", key, LAL_KEYNAME_MAX);
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

/* DICT ROUTINES */

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

LALDictEntry *XLALDictLookup(LALDict *dict, const char *key)
{
	LALDictEntry *entry;
	for (entry = dict->hashes[hash(key) % dict->size]; entry != NULL; entry = entry->next)
		if (strcmp(key, entry->key) == 0)
			return entry;
	return NULL;
}

int XLALDictRemove(LALDict *dict, const char *key)
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
			LALFree(this);
			return 0;
		}
		prev = this;
		this = this->next;
	}
	return -1; /* not found */
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

int XLALDictInsertStringValue(LALDict *dict, const char *key, const char *value)
{
	size_t size = strlen(value) + 1;
	if (XLALDictInsert(dict, key, value, size, LAL_CHAR_TYPE_CODE) < 0)
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

/* warning: shallow pointer */
const char * XLALDictLookupStringValue(LALDict *dict, const char *key)
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
	TYPE XLALDictLookup ## TYPE ## Value(LALDict *dict, const char *key) \
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

REAL8 XLALDictLookupValueAsREAL8(LALDict *dict, const char *key)
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

static void XLALDictEntryPrintFunc(char *key, LALValue *value, void *thunk)
{
	int fd = *(int *)(thunk);
	int fd2 = dup(fd);
	FILE *fp = fdopen(fd2, "w");
	fprintf(fp, "\"%s\": ", key);
	XLALValuePrint(value, fd);
	fprintf(fp, "\n");
	fclose(fp);
	return;
}

void XLALDictPrint(LALDict *dict, int fd)
{
	XLALDictForeach(dict, XLALDictEntryPrintFunc, &fd);
	return;
}
