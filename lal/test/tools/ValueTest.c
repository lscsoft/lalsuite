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

#define CHAR_VALUE LAL_INT8_C(+111)
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
char BLOB_VALUE[5] = {'\x68', '\x65', '\x6c', '\x6c', '\x6f'};
char String_VALUE[] = "world";

static LALDict * create_dict(void)
{
    LALDict *dict;
    LALDict *intdict;
    LALDict *floatdict;
    LALDict *complexdict;
    LALDict *blobdict;
    LALDict *strdict;
    int err = 0;
    intdict = XLALCreateDict();
    floatdict = XLALCreateDict();
    complexdict = XLALCreateDict();
    blobdict = XLALCreateDict();
    strdict = XLALCreateDict();
    if (!intdict || !floatdict || !complexdict || !blobdict || !strdict)
        return NULL;
    err |= XLALDictInsertCHARValue(intdict, "CHAR", CHAR_VALUE);
    err |= XLALDictInsertINT2Value(intdict, "INT2", INT2_VALUE);
    err |= XLALDictInsertINT4Value(intdict, "INT4", INT4_VALUE);
    err |= XLALDictInsertINT8Value(intdict, "INT8", INT8_VALUE);
    err |= XLALDictInsertUCHARValue(intdict, "UCHAR", UCHAR_VALUE);
    err |= XLALDictInsertUINT2Value(intdict, "UINT2", UINT2_VALUE);
    err |= XLALDictInsertUINT4Value(intdict, "UINT4", UINT4_VALUE);
    err |= XLALDictInsertUINT8Value(intdict, "UINT8", UINT8_VALUE);
    err |= XLALDictInsertREAL4Value(floatdict, "REAL4", REAL4_VALUE);
    err |= XLALDictInsertREAL8Value(floatdict, "REAL8", REAL8_VALUE);
    err |= XLALDictInsertCOMPLEX8Value(complexdict, "COMPLEX8", COMPLEX8_VALUE);
    err |= XLALDictInsertCOMPLEX16Value(complexdict, "COMPLEX16", COMPLEX16_VALUE);
    err |= XLALDictInsertBLOBValue(blobdict, "BLOB", BLOB_VALUE, sizeof(BLOB_VALUE));
    err |= XLALDictInsertStringValue(strdict, "String", String_VALUE);
    dict = XLALDictMerge(intdict, floatdict);
    XLALDictUpdate(dict, complexdict);
    XLALDictUpdate(dict, blobdict);
    XLALDictUpdate(dict, strdict);
    XLALDestroyDict(strdict);
    XLALDestroyDict(blobdict);
    XLALDestroyDict(complexdict);
    XLALDestroyDict(floatdict);
    XLALDestroyDict(intdict);
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
    fprintf(stderr, "Testing %s... ", #TYPE); \
    entry = XLALDictLookup(dict, #TYPE); \
    orig = XLALDictEntryGetValue(entry); \
    size = XLALValueGetSize(orig); \
    copy = XLALValueRealloc(copy, size); \
    copy = XLALValueCopy(copy, orig); \
    if (!XLALValueEqual(copy, orig)) { \
        fprintf(stderr, "failed:"); \
        XLALValuePrint(copy, 2); \
        fprintf(stderr, " != "); \
        XLALValuePrint(orig, 2); \
        fprintf(stderr, "\n"); \
        return 1; \
    } \
    if (!COMPARE(XLALValueGet ## TYPE(copy), TYPE)) { \
        fprintf(stderr, "failed: incorrect value\n"); \
        return 1; \
    } \
    XLALDictRemove(dict, #TYPE); \
    fprintf(stderr, "passed\n");

#define TEST2(TYPE) \
    fprintf(stderr, "Testing %s Pop... ", #TYPE); \
    if (!COMPARE(XLALDictPop ## TYPE ## Value(dict2, #TYPE), TYPE)) { \
        fprintf(stderr, "failed: incorrect value\n"); \
        return 1; \
    } \
    fprintf(stderr, "passed\n");

int main(void)
{
    LALDict *dict;
    LALDict *dict2;
    LALList *list;
    LALList *keys;

    LALDictEntry *entry;
    const LALValue *orig;
    LALValue *copy = NULL;
    size_t size;
    void *blob;
    char *str;

    /* make sure that the keys in the dict are what they should be */
    list = create_list();
    dict = create_dict();
    keys = XLALDictKeys(dict);
    XLALListSort(keys, string_value_cmp, NULL);
    if (!lists_are_equal(list, keys))
        return 1;

    dict2 = XLALDictDuplicate(dict);

    /* make sure the values in the dict are what they should be */
    TEST(CHAR)
    TEST2(CHAR)
    TEST(INT2)
    TEST2(INT2)
    TEST(INT4)
    TEST2(INT4)
    TEST(INT8)
    TEST2(INT8)
    TEST(UCHAR)
    TEST2(UCHAR)
    TEST(UINT2)
    TEST2(UINT2)
    TEST(UINT4)
    TEST2(UINT4)
    TEST(UINT8)
    TEST2(UINT8)
    TEST(REAL4)
    TEST2(REAL4)
    TEST(REAL8)
    TEST2(REAL8)
    TEST(COMPLEX8)
    TEST2(COMPLEX8)
    TEST(COMPLEX16)
    TEST2(COMPLEX16)

#undef COMPARE
#define COMPARE(v, TYPE) (!memcmp(blob = v, BLOB_VALUE, sizeof(BLOB_VALUE)))
    TEST(BLOB)
    LALFree(blob);
    TEST2(BLOB)
    LALFree(blob);

#undef COMPARE
#define COMPARE(v, TYPE) (!strcmp(v, String_VALUE))
    TEST(String)
#undef COMPARE
#define COMPARE(v, TYPE) (!strcmp(str = v, String_VALUE))
    TEST2(String)
    LALFree(str);

    /* dict should now be empty */
    if (XLALDictSize(dict) != 0) {
        fprintf(stderr, "Error: Dictionary has %zu leftover elements\n", XLALDictSize(dict));
        return 1;
    }

    /* dict2 should now be empty */
    if (XLALDictSize(dict2) != 0) {
        fprintf(stderr, "Error: Dictionary has %zu leftover elements\n", XLALDictSize(dict2));
        return 1;
    }

    XLALDestroyDict(dict2);
    XLALDestroyDict(dict);
    XLALDestroyList(keys);
    XLALDestroyList(list);
    XLALDestroyValue(copy);

    LALCheckMemoryLeaks();
    return 0;
}
