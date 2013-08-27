/*
 *  Copyright (C) 2013 Chris Pankow, Karl Wette
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

#include <lal/LALDetCharGlibTypes.h>
#include <lal/XLALError.h>
#include <lal/LALString.h>
#include <lal/LIGOMetadataBurstUtils.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

struct tagLALGSequence {
  LALGType type;
  GSequence* s;
};

struct tagLALGSequenceIter {
  LALGType type;
  GSequenceIter* i;
};

struct tagLALGHashTable {
  LALGType type;
  GHashTable* t;
};

static gint g_compare_SnglBurst(gconstpointer a, gconstpointer b, gpointer UNUSED user_data) {
  const SnglBurst *_a = a;
  const SnglBurst *_b = b;
  return XLALCompareSnglBurstByPeakTimeAndSNR(&_a, &_b);
}

LALGSequence* XLALCreateGSequence(LALGType type) {
  LALGSequence* seq = (LALGSequence*)XLALMalloc(sizeof(LALGSequence*));
  XLAL_CHECK_NULL(seq, XLAL_ENOMEM);
  seq->type = type;
  switch (type) {
  case LALGTYPE_SNGL_BURST:
    seq->s = g_sequence_new((GDestroyNotify)XLALDestroySnglBurst);
    break;
  default:
    XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid LALGType=%i", type);
  }
  XLAL_CHECK_NULL(seq->s, XLAL_ENOMEM);
  return seq;
}

void XLALDestroyGSequence(LALGSequence* seq) {
  if (seq) {
    g_sequence_free(seq->s);
    XLALFree(seq);
  }
}

LALGType XLALGetGSequenceType(LALGSequence* seq) {
  XLAL_CHECK_VAL(LALGTYPE_NONE, seq, XLAL_EFAULT);
  return seq->type;
}

size_t XLALGetGSequenceLength(LALGSequence* seq){
  XLAL_CHECK_VAL(0, seq, XLAL_EFAULT);
  return g_sequence_get_length(seq->s);
}

void XLALDestroyGSequenceIter(LALGSequenceIter* itr) {
  if (itr) {
    XLALFree(itr);
  }
}

LALGSequenceIter* XLALGSequenceBegin(LALGSequence* seq) {
  XLAL_CHECK_NULL(seq, XLAL_EFAULT);
  GSequenceIter* i = g_sequence_get_begin_iter(seq->s);
  if (g_sequence_iter_is_end(i)) {
    return NULL;
  } else {
    LALGSequenceIter* itr = (LALGSequenceIter*)XLALMalloc(sizeof(LALGSequenceIter*));
    XLAL_CHECK_NULL(itr, XLAL_ENOMEM);
    itr->type = seq->type;
    itr->i = i;
    return itr;
  }
}

LALGSequenceIter* XLALGSequenceNext(LALGSequenceIter* itr) {
  XLAL_CHECK_NULL(itr, XLAL_EFAULT);
  itr->i = g_sequence_iter_next(itr->i);
  if (g_sequence_iter_is_end(itr->i)) {
    XLALDestroyGSequenceIter(itr);
    return NULL;
  }
  return itr;
}

GSequenceIter* XLALGSequenceBeginRaw(LALGSequence* seq, LALGType type) {
  XLAL_CHECK_NULL(seq, XLAL_EFAULT);
  XLAL_CHECK_NULL(seq->type == type, XLAL_EINVAL);
  return g_sequence_get_begin_iter(seq->s);
}

SnglBurst* XLALGetGSeqSnglBurst(LALGSequenceIter* itr) {
  XLAL_CHECK_NULL(itr, XLAL_EFAULT);
  XLAL_CHECK_NULL(itr->type == LALGTYPE_SNGL_BURST, XLAL_EINVAL);
  return (SnglBurst*)g_sequence_get(itr->i);
}

LALGSequenceIter* XLALAddGSeqSnglBurst(LALGSequence* seq, SnglBurst* sb) {
  XLAL_CHECK_NULL(seq, XLAL_EFAULT);
  XLAL_CHECK_NULL(seq->type == LALGTYPE_SNGL_BURST, XLAL_EINVAL);
  XLAL_CHECK_NULL(sb, XLAL_EFAULT);
  GSequenceIter* i = g_sequence_insert_sorted(seq->s, sb, g_compare_SnglBurst, NULL);
  if (g_sequence_iter_is_end(i)) {
    return NULL;
  } else {
    LALGSequenceIter* itr = (LALGSequenceIter*)XLALMalloc(sizeof(LALGSequenceIter*));
    XLAL_CHECK_NULL(itr, XLAL_ENOMEM);
    itr->type = seq->type;
    itr->i = i;
    return itr;
  }
}

LALGHashTable* XLALCreateGHashTable(LALGType type) {
  LALGHashTable* tbl = (LALGHashTable*)XLALMalloc(sizeof(LALGHashTable*));
  XLAL_CHECK_NULL(tbl, XLAL_ENOMEM);
  tbl->type = type;
  switch (type) {
  case LALGTYPE_INT:
  case LALGTYPE_DBL:
  case LALGTYPE_STR:
    tbl->t = g_hash_table_new_full(g_str_hash, g_str_equal, XLALFree, XLALFree);
    break;
  default:
    XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid LALGType=%i", type);
  }
  XLAL_CHECK_NULL(tbl->t, XLAL_ENOMEM);
  return tbl;
}

void XLALDestroyGHashTable(LALGHashTable* tbl) {
  if (tbl) {
    g_hash_table_destroy(tbl->t);
    XLALFree(tbl);
  }
}

LALGType XLALGetGHashTableType(LALGHashTable* tbl) {
  XLAL_CHECK_VAL(LALGTYPE_NONE, tbl, XLAL_EFAULT);
  return tbl->type;
}

bool XLALGHashTableKeyExists(LALGHashTable* tbl, const char* key) {
  XLAL_CHECK_VAL(false, tbl, XLAL_EFAULT);
  XLAL_CHECK_VAL(false, key, XLAL_EFAULT);
  return g_hash_table_lookup(tbl->t, key) != NULL;
}

char* XLALGetGHashTableKey(LALGHashTable* tbl, const size_t indx) {
  XLAL_CHECK_NULL(tbl, XLAL_EFAULT);
  GList* keys = g_hash_table_get_keys(tbl->t);
  size_t i = 0;
  while (keys){
    if (i == indx) {
      return XLALStringDuplicate((char*)keys->data);
    }
    i++;
    keys = keys->next;
  }
  return NULL;
}

void XLALGHashTableBeginRaw(LALGHashTable* tbl, LALGType type, GHashTableIter* itr) {
  XLAL_CHECK_VOID(tbl, XLAL_EFAULT);
  XLAL_CHECK_VOID(tbl->type == type, XLAL_EINVAL);
  XLAL_CHECK_VOID(itr, XLAL_EFAULT);
  g_hash_table_iter_init(itr, tbl->t);
}

void XLALSetGHashTblInt(LALGHashTable* tbl, const char* key, const int val) {
  XLAL_CHECK_VOID(tbl, XLAL_EFAULT);
  XLAL_CHECK_VOID(tbl->type == LALGTYPE_INT, XLAL_EINVAL);
  XLAL_CHECK_VOID(key, XLAL_EFAULT);
  int *newval = XLALMalloc(sizeof(int));
  XLAL_CHECK_VOID(newval, XLAL_ENOMEM);
  *newval = val;
  g_hash_table_insert(tbl->t, XLALStringDuplicate(key), newval);
}

int XLALGetGHashTblInt(LALGHashTable* tbl, const char* key) {
  XLAL_CHECK(tbl, XLAL_EFAULT);
  XLAL_CHECK(tbl->type == LALGTYPE_INT, XLAL_EINVAL);
  XLAL_CHECK(key, XLAL_EFAULT);
  int *val = (int*)g_hash_table_lookup(tbl->t, key);
  XLAL_CHECK(val, XLAL_EDOM);
  return *val;
}

int* XLALGetGHashTblIntPtr(LALGHashTable* tbl, const char* key) {
  XLAL_CHECK_NULL(tbl, XLAL_EFAULT);
  XLAL_CHECK_NULL(tbl->type == LALGTYPE_INT, XLAL_EINVAL);
  XLAL_CHECK_NULL(key, XLAL_EFAULT);
  return (int*)g_hash_table_lookup(tbl->t, key);
}

void XLALSetGHashTblDbl(LALGHashTable* tbl, const char* key, const double val) {
  XLAL_CHECK_VOID(tbl, XLAL_EFAULT);
  XLAL_CHECK_VOID(tbl->type == LALGTYPE_DBL, XLAL_EINVAL);
  XLAL_CHECK_VOID(key, XLAL_EFAULT);
  double *newval = XLALMalloc(sizeof(double));
  XLAL_CHECK_VOID(newval, XLAL_ENOMEM);
  *newval = val;
  g_hash_table_insert(tbl->t, XLALStringDuplicate(key), newval);
}

double XLALGetGHashTblDbl(LALGHashTable* tbl, const char* key) {
  XLAL_CHECK_REAL8(tbl, XLAL_EFAULT);
  XLAL_CHECK_REAL8(tbl->type == LALGTYPE_DBL, XLAL_EINVAL);
  XLAL_CHECK_REAL8(key, XLAL_EFAULT);
  double *val = (double*)g_hash_table_lookup(tbl->t, key);
  XLAL_CHECK_REAL8(val, XLAL_EDOM);
  return *val;
}

double* XLALGetGHashTblDblPtr(LALGHashTable* tbl, const char* key) {
  XLAL_CHECK_NULL(tbl, XLAL_EFAULT);
  XLAL_CHECK_NULL(tbl->type == LALGTYPE_DBL, XLAL_EINVAL);
  XLAL_CHECK_NULL(key, XLAL_EFAULT);
  return (double*)g_hash_table_lookup(tbl->t, key);
}

void XLALSetGHashTblStr(LALGHashTable* tbl, const char* key, const char* val) {
  XLAL_CHECK_VOID(tbl, XLAL_EFAULT);
  XLAL_CHECK_VOID(tbl->type == LALGTYPE_STR, XLAL_EINVAL);
  XLAL_CHECK_VOID(key, XLAL_EFAULT);
  XLAL_CHECK_VOID(val, XLAL_EFAULT);
  g_hash_table_insert(tbl->t, XLALStringDuplicate(key), XLALStringDuplicate(val));
}

char* XLALGetGHashTblStr(LALGHashTable* tbl, const char* key) {
  XLAL_CHECK_NULL(tbl, XLAL_EFAULT);
  XLAL_CHECK_NULL(tbl->type == LALGTYPE_STR, XLAL_EINVAL);
  XLAL_CHECK_NULL(key, XLAL_EFAULT);
  char *val = (char*)g_hash_table_lookup(tbl->t, key);
  XLAL_CHECK_NULL(val, XLAL_EDOM);
  return XLALStringDuplicate(val);
}
