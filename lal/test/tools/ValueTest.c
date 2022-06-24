#include <lal/LALConfig.h>

#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALValue.h>
#include <lal/LALDict.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define CHAR_VALUE LAL_INT8_C(-111)
#define INT2_VALUE LAL_INT8_C(-11111)
#define INT4_VALUE LAL_INT8_C(-1111111111)
#define INT8_VALUE LAL_INT8_C(-111111111111111111)
#define UCHAR_VALUE LAL_UINT8_C(222)
#define UINT2_VALUE LAL_UINT8_C(44444)
#define UINT4_VALUE LAL_UINT8_C(3333333333)
#define UINT8_VALUE LAL_UINT8_C(11111111111111111111)
#define REAL4_VALUE 100.0
#define REAL8_VALUE 1e100
#define COMPLEX8_VALUE 3.0 + 4.0 * I
#define COMPLEX16_VALUE 3e100 + 4e100 * I
char BLOB_VALUE[5] = "\x68\x65\x6c\x6c\x6f";
char String_VALUE[] = "world";

static LALDict * create_dict(void)
{
    LALDict *dict;
    int err = 0;
    dict = XLALCreateDict();
    if (!dict)
        return NULL;
    err |= XLALDictInsertCHARValue(dict, "CHAR", CHAR_VALUE);
    err |= XLALDictInsertINT2Value(dict, "INT2", INT2_VALUE);
    err |= XLALDictInsertINT4Value(dict, "INT4", INT4_VALUE);
    err |= XLALDictInsertINT8Value(dict, "INT8", INT8_VALUE);
    err |= XLALDictInsertUCHARValue(dict, "UCHAR", UCHAR_VALUE);
    err |= XLALDictInsertUINT2Value(dict, "UINT2", UINT2_VALUE);
    err |= XLALDictInsertUINT4Value(dict, "UINT4", UINT4_VALUE);
    err |= XLALDictInsertUINT8Value(dict, "UINT8", UINT8_VALUE);
    err |= XLALDictInsertREAL4Value(dict, "REAL4", REAL4_VALUE);
    err |= XLALDictInsertREAL8Value(dict, "REAL8", REAL8_VALUE);
    err |= XLALDictInsertCOMPLEX8Value(dict, "COMPLEX8", COMPLEX8_VALUE);
    err |= XLALDictInsertCOMPLEX16Value(dict, "COMPLEX16", COMPLEX16_VALUE);
    err |= XLALDictInsertBLOBValue(dict, "BLOB", BLOB_VALUE, sizeof(BLOB_VALUE));
    err |= XLALDictInsertStringValue(dict, "String", String_VALUE);
    return err ? NULL : dict;
}

static LALList * create_list(void)
{
    LALList *list;
    int err = 0;
    list = XLALCreateList();
    if (!list)
        return NULL;
    /* put these in alphabetical order */
    err |= XLALListAddStringValue(list, "BLOB");
    err |= XLALListAddStringValue(list, "CHAR");
    err |= XLALListAddStringValue(list, "COMPLEX16");
    err |= XLALListAddStringValue(list, "COMPLEX8");
    err |= XLALListAddStringValue(list, "INT2");
    err |= XLALListAddStringValue(list, "INT4");
    err |= XLALListAddStringValue(list, "INT8");
    err |= XLALListAddStringValue(list, "REAL4");
    err |= XLALListAddStringValue(list, "REAL8");
    err |= XLALListAddStringValue(list, "String");
    err |= XLALListAddStringValue(list, "UCHAR");
    err |= XLALListAddStringValue(list, "UINT2");
    err |= XLALListAddStringValue(list, "UINT4");
    err |= XLALListAddStringValue(list, "UINT8");
    /* inserting them in this order puts them in reverse order... */
    XLALListReverse(list);
    return err ? NULL : list;
}

static int lists_are_equal(LALList *list1, LALList *list2)
{
    LALListIter iter1;
    LALListIter iter2;
    LALListItem *item1;
    LALListItem *item2;
    if (XLALListSize(list1) != XLALListSize(list2))
        return 0;
    XLALListIterInit(&iter1, list1);
    XLALListIterInit(&iter2, list2);
    while ((item1 = XLALListIterNext(&iter1)) && (item2 = XLALListIterNext(&iter2))) {
        const LALValue *value1 = XLALListItemGetValue(item1);
        const LALValue *value2 = XLALListItemGetValue(item2);
        if (!XLALValueEqual(value1, value2))
            return 0;
    }
    return 1;
}

static int string_value_cmp(const LALValue *value1, const LALValue *value2, void UNUSED *thunk)
{
    return strcmp(XLALValueGetString(value1), XLALValueGetString(value2));
}

#define COMPARE(v, TYPE) (v == TYPE ## _VALUE)

#define TEST(TYPE) \
    fprintf(stderr, "Testing %s...", #TYPE); \
    entry = XLALDictLookup(dict, #TYPE); \
    orig = XLALDictEntryGetValue(entry); \
    size = XLALValueGetSize(orig); \
    copy = XLALValueRealloc(copy, size); \
    copy = XLALValueCopy(copy, orig); \
    if (!XLALValueEqual(copy, orig)) \
        return 1; \
    if (!COMPARE(XLALValueGet ## TYPE(copy), TYPE)) \
        return 1; \
    XLALDictRemove(dict, #TYPE); \
    fprintf(stderr, " passed\n");

int main(void)
{
    LALDict *dict;
    LALList *list;
    LALList *keys;

    LALDictEntry *entry;
    const LALValue *orig;
    LALValue *copy = NULL;
    size_t size;
    void *blob;

    /* make sure that the keys in the dict are what they should be */
    list = create_list();
    dict = create_dict();
    keys = XLALDictKeys(dict);
    XLALListSort(keys, string_value_cmp, NULL);
    if (!lists_are_equal(list, keys))
        return 1;

    /* make sure the values in the dict are what they should be */
    TEST(CHAR)
    TEST(INT2)
    TEST(INT4)
    TEST(INT8)
    TEST(UCHAR)
    TEST(UINT2)
    TEST(UINT4)
    TEST(UINT8)
    TEST(REAL4)
    TEST(REAL8)
    TEST(COMPLEX8)
    TEST(COMPLEX16)

#undef COMPARE
#define COMPARE(v, TYPE) (!memcmp(blob = v, BLOB_VALUE, sizeof(BLOB_VALUE)))
    TEST(BLOB)
    LALFree(blob);

#undef COMPARE
#define COMPARE(v, TYPE) (!strcmp(v, String_VALUE))
    TEST(String)

    /* dict should now be empty */
    if (XLALDictSize(dict) != 0)
        return 1;

    XLALDestroyDict(dict);
    XLALDestroyList(keys);
    XLALDestroyList(list);
    XLALDestroyValue(copy);

    LALCheckMemoryLeaks();
    return 0;
}
