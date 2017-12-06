/*
*  Copyright (C) 2016 Karl Wette
*  Copyright (C) 2007 Jolien Creighton, Josh Willis
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#include <config.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdio.h>
#include <lal/LALError.h>

/* global variables */
size_t lalMallocTotal = 0;	/**< current amount of memory allocated by process */
size_t lalMallocTotalPeak = 0;	/**< peak amount of memory allocated so far */

/*
 *
 * XLAL Routines.
 *
 */

#define XLAL_TEST_POINTER( ptr, size )                                    \
    if ( ! (ptr) && (size) )                                              \
       XLAL_ERROR_NULL( XLAL_ENOMEM );                                    \
    else (void)(0)
#define XLAL_TEST_POINTER_LONG( ptr, size, file, line )                   \
    if ( ! (ptr) && (size) )                                              \
    {                                                                     \
       XLALPrintError( "XLALError - %s in %s:%d", __func__, file, line ); \
       XLAL_ERROR_NULL( XLAL_ENOMEM );                                    \
    }                                                                     \
    else (void)(0)
#define XLAL_TEST_POINTER_ALIGNED( ptr, size, retval )                    \
    if ( ! (ptr) && (size) && (retval) )                                  \
       XLAL_ERROR_NULL( XLAL_ENOMEM );                                    \
    else (void)(0)
#define XLAL_TEST_POINTER_ALIGNED_LONG( ptr, size, retval, file, line )   \
    if ( ! (ptr) && (size) && (retval) )                                  \
    {                                                                     \
       XLALPrintError( "XLALError - %s in %s:%d", __func__, file, line ); \
       XLAL_ERROR_NULL( XLAL_ENOMEM );                                    \
    }                                                                     \
    else (void)(0)

void *(XLALMalloc) (size_t n) {
    void *p;
    p = LALMallocShort(n);
    XLAL_TEST_POINTER(p, n);
    return p;
}

void *XLALMallocLong(size_t n, const char *file, int line)
{
    void *p;
    p = LALMallocLong(n, file, line);
    XLAL_TEST_POINTER_LONG(p, n, file, line);
    return p;
}

void *(XLALCalloc) (size_t m, size_t n) {
    void *p;
    p = LALCallocShort(m, n);
    XLAL_TEST_POINTER(p, m * n);
    return p;
}

void *XLALCallocLong(size_t m, size_t n, const char *file, int line)
{
    void *p;
    p = LALCallocLong(m, n, file, line);
    XLAL_TEST_POINTER_LONG(p, m * n, file, line);
    return p;
}

void *(XLALRealloc) (void *p, size_t n) {
    p = LALReallocShort(p, n);
    XLAL_TEST_POINTER(p, n);
    return p;
}

void *XLALReallocLong(void *p, size_t n, const char *file, int line)
{
    p = LALReallocLong(p, n, file, line);
    XLAL_TEST_POINTER_LONG(p, n, file, line);
    return p;
}

void XLALFree(void *p)
{
    if (p)
        LALFree(p);
    return;
}

/*
 * Aligned memory routines.
 */

#if LAL_FFTW3_MEMALIGN_ENABLED

#ifndef HAVE_POSIX_MEMALIGN
#error no posix_memalign available
#endif

int XLALIsMemoryAligned(void *ptr)
{
	return LAL_IS_MEMORY_ALIGNED(ptr);
}

void *XLALMallocAlignedLong(size_t size, const char *file, int line)
{
	void *p=NULL;
	int retval;
	retval = posix_memalign(&p, LAL_MEM_ALIGNMENT, size);
	XLAL_TEST_POINTER_ALIGNED_LONG(p, size, retval, file, line);
	return p;
}

void *(XLALMallocAligned)(size_t size)
{
	void *p=NULL;
	int retval;
	retval = posix_memalign(&p, LAL_MEM_ALIGNMENT, size);
	XLAL_TEST_POINTER_ALIGNED(p, size, retval);
	return p;
}

void *XLALCallocAlignedLong(size_t nelem, size_t elsize, const char *file, int line)
{
	size_t size = nelem * elsize;
	void *p = XLALMallocAlignedLong(size, file, line);
	XLAL_TEST_POINTER_LONG(p, size, file, line);
	memset(p, 0, size);
	return p;
}

void *(XLALCallocAligned)(size_t nelem, size_t elsize)
{
	size_t size = nelem * elsize;
	void *p = (XLALMallocAligned)(size);
	XLAL_TEST_POINTER(p, size);
	memset(p, 0, size);
	return p;
}

void *XLALReallocAlignedLong(void *ptr, size_t size, const char *file, int line)
{
	void *p;
	if (ptr == NULL)
		return XLALMallocAlignedLong(size, file, line);
	if (size == 0) {
		XLALFreeAligned(ptr);
		return NULL;
	}
	p = realloc(ptr, size); /* use ordinary realloc */
	if (XLALIsMemoryAligned(p))
		return p;
	/* need to do a new allocation and a memcpy, inefficient... */
	ptr = XLALMallocAlignedLong(size, file, line);
	memcpy(ptr, p, size);
	XLALFree(p);
	return ptr;
}

void *(XLALReallocAligned)(void *ptr, size_t size)
{
	void *p;
	if (ptr == NULL)
		return XLALMallocAligned(size);
	if (size == 0) {
		XLALFreeAligned(ptr);
		return NULL;
	}
	p = realloc(ptr, size); /* use ordinary realloc */
	if (XLALIsMemoryAligned(p))
		return p;
	/* need to do a new allocation and a memcpy, inefficient... */
	ptr = XLALMallocAligned(size);
	memcpy(ptr, p, size);
	XLALFree(p);
	return ptr;
}

void XLALFreeAligned(void *ptr)
{
	free(ptr); /* use ordinary free */
}

#endif /* LAL_FFTW3_MEMALIGN_ENABLED */

/*
 *
 * LAL Routines... only if compiled with debugging enabled.
 * (otherwise the LALMalloc-family reverts to the standard malloc-family).
 *
 */


#if ! defined NDEBUG

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#else
#define pthread_mutex_lock( pmut )
#define pthread_mutex_unlock( pmut )
#endif

#include <lal/LALStdlib.h>

/* global variables to assist in memory debugging */
/* watch the value of these variables to find a particular alloc/free */
char *lalMemDbgArgPtr = NULL;   /* set to ptr arg in free or realloc */
char *lalMemDbgRetPtr = NULL;   /* set to ptr returned in alloc functions */
char *lalMemDbgPtr = NULL;      /* set in both cases */
char *lalMemDbgUsrPtr = NULL;   /* avaliable global memory pointer for user */
void **lalMemDbgUsrHndl = NULL; /* avaliable global memory handle for user */
int lalIsMemDbgArgPtr;  /* ( lalMemDbgUsrPtr == lalMemDbgArgPtr ) */
int lalIsMemDbgRetPtr;  /* ( lalMemDbgUsrPtr == lalMemDbgRetPtr ) */
int lalIsMemDbgPtr;     /* ( lalMemDbgUsrPtr == lalMemDbgPtr ) */


enum { nprefix = 2 };
static const size_t prefix = nprefix * sizeof(size_t);
static const size_t padFactor = 2;
static const size_t padding = 0xDeadBeef;
static const size_t repadding = 0xBeefDead;
static const size_t magic = 0xABadCafe;

#define allocsz(n) ((lalDebugLevel & LALMEMPADBIT) ? (padFactor * (n) + prefix) : (n))

/* Hash table implementation taken from src/utilities/LALHashTbl.c */

static struct allocNode {
    void *addr;
    size_t size;
    const char *file;
    int line;
} **alloc_data = NULL;		/* Allocation hash table with open addressing and linear probing */
static int alloc_data_len = 0;	/* Size of the memory block 'alloc_data', in number of elements */
static int alloc_n = 0;		/* Number of valid elements in the hash */
static int alloc_q = 0;		/* Number of non-NULL elements in the hash */

/* Special allocation hash table element value to indicate elements that have been deleted */
static const void *hash_del = 0;
#define DEL   ((struct allocNode*) &hash_del)

/* Evaluates to the hash value of x, restricted to the length of the allocation hash table */
#define HASHIDX(x)   ((int)( ((intptr_t)( (x)->addr )) % alloc_data_len ))

/* Increment the next hash index, restricted to the length of the allocation hash table */
#define INCRIDX(i)   do { if (++(i) == alloc_data_len) { (i) = 0; } } while(0)

/* Evaluates true if the elements x and y are equal */
#define EQUAL(x, y)   ((x)->addr == (y)->addr)

/* need this to turn off gcc warnings about unused functions */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Resize and rebuild the allocation allocation hash table */
UNUSED static int AllocHashTblResize(void)
{
    struct allocNode **old_data = alloc_data;
    int old_data_len = alloc_data_len;
    alloc_data_len = 2;
    while (alloc_data_len < 3*alloc_n) {
        alloc_data_len *= 2;
    }
    alloc_data = calloc(alloc_data_len, sizeof(alloc_data[0]));
    if (alloc_data == NULL) {
        return 0;
    }
    alloc_q = alloc_n;
    for (int k = 0; k < old_data_len; ++k) {
        if (old_data[k] != NULL && old_data[k] != DEL) {
            int i = HASHIDX(old_data[k]);
            while (alloc_data[i] != NULL) {
                INCRIDX(i);
            }
            alloc_data[i] = old_data[k];
        }
    }
    free(old_data);
    return 1;
}

/* Find node in allocation hash table */
UNUSED static struct allocNode *AllocHashTblFind(struct allocNode *x)
{
    struct allocNode *y = NULL;
    if (alloc_data_len > 0) {
        int i = HASHIDX(x);
        while (alloc_data[i] != NULL) {
            y = alloc_data[i];
            if (y != DEL && EQUAL(x, y)) {
                return y;
            }
            INCRIDX(i);
        }
    }
    return NULL;
}

/* Add node to allocation hash table */
UNUSED static int AllocHashTblAdd(struct allocNode *x)
{
    if (2*(alloc_q + 1) > alloc_data_len) {
        /* Resize allocation hash table to preserve maximum 50% occupancy */
        if (!AllocHashTblResize()) {
            return 0;
        }
    }
    int i = HASHIDX(x);
    while (alloc_data[i] != NULL && alloc_data[i] != DEL) {
        INCRIDX(i);
    }
    if (alloc_data[i] == NULL) {
        ++alloc_q;
    }
    ++alloc_n;
    alloc_data[i] = x;
    return 1;
}

/* Extract node from allocation hash table */
UNUSED static struct allocNode *AllocHashTblExtract(struct allocNode *x)
{
    if (alloc_data_len > 0) {
        int i = HASHIDX(x);
        while (alloc_data[i] != NULL) {
            struct allocNode *y = alloc_data[i];
            if (y != DEL && EQUAL(x, y)) {
                alloc_data[i] = DEL;
                --alloc_n;
                if (alloc_n == 0) {
                    /* Free all hash table memory */
                    free(alloc_data);
                    alloc_data = NULL;
                    alloc_data_len = 0;
                    alloc_q = 0;
                } else if (8*alloc_n < alloc_data_len) {
                    /* Resize hash table to preserve minimum 50% occupancy */
                    if (!AllocHashTblResize()) {
                        return NULL;
                    }
                }
                return y;
            }
            INCRIDX(i);
        }
    }
    return NULL;
}


/* Useful function for debugging */
/* Checks to make sure alloc list is OK */
/* Returns 0 if list is corrupted; 1 if list is OK */
UNUSED static int CheckAllocList(void)
{
    int count = 0;
    size_t total = 0;
    for (int k = 0; k < alloc_data_len; ++k) {
        if (alloc_data[k] != NULL && alloc_data[k] != DEL) {
            ++count;
            total += alloc_data[k]->size;
        }
    }
    return count == alloc_n && total == lalMallocTotal;
}

/* Useful function for debugging */
/* Finds the node of the alloc list for the desired alloc */
/* Returns NULL if not found  */
UNUSED static struct allocNode *FindAlloc(void *p)
{
    struct allocNode key = { .addr = p };
    return AllocHashTblFind(&key);
}


static void *PadAlloc(size_t * p, size_t n, int keep, const char *func)
{
    size_t i;

    if (!(lalDebugLevel & LALMEMPADBIT)) {
        return p;
    }

    if (!p) {
        return NULL;
    }

    if (lalDebugLevel & LALMEMINFOBIT) {
        XLALPrintError("%s meminfo: allocating %zu bytes at address %p\n",
                      func, n, p + nprefix);
    }

    /* store the size in a known position */
    p[0] = n;
    p[1] = magic;

    /* pad the memory */
    for (i = keep ? n : 0; i < padFactor * n; ++i) {
        ((char *) p)[i + prefix] = (char) (i ^ padding);
    }

    pthread_mutex_lock(&mut);
    lalMallocTotal += n;
    lalMallocTotalPeak = (lalMallocTotalPeak > lalMallocTotal) ? lalMallocTotalPeak : lalMallocTotal;
    pthread_mutex_unlock(&mut);

    return (void *) (((char *) p) + prefix);
}


static void *UnPadAlloc(void *p, int keep, const char *func)
{
    size_t n;
    size_t i;
    size_t *q;
    char *s;

    if (!(lalDebugLevel & LALMEMPADBIT)) {
        return p;
    }

    if (!p || !(q = ((size_t *) p) - nprefix)) {
        lalRaiseHook(SIGSEGV, "%s error: tried to free NULL pointer\n",
                     func);
        return NULL;
    }

    n = q[0];
    s = (char *) q;

    if (lalDebugLevel & LALMEMINFOBIT) {
        XLALPrintError("%s meminfo: freeing %zu bytes at address %p\n",
                      func, n, p);
    }

    if (n == (size_t)(-1)) {
        lalRaiseHook(SIGSEGV,
                     "%s error: tried to free a freed pointer at address %p\n",
                     func, p);
        return NULL;
    }

    if (q[1] != magic) {
        lalRaiseHook(SIGSEGV,
                     "%s error: wrong magic for pointer at address %p\n",
                     func, p);
        return NULL;
    }

    if (((long) n) < 0) {
        lalRaiseHook(SIGSEGV,
                     "%s error: corrupt size descriptor for pointer at address %p\n",
                     func, p);
        return NULL;
    }

    /* check for writing past end of array: */
    for (i = n; i < padFactor * n; ++i) {
        if (s[i + prefix] != (char) (i ^ padding)) {
            lalRaiseHook(SIGSEGV, "%s error: array bounds overwritten\n"
                         "Byte %ld past end of array has changed\n"
                         "Corrupted address: %p\nArray address: %p\n",
                         func, i - n + 1, s + i + prefix, s + prefix);
            return NULL;
        }
    }

    /* see if there is enough allocated memory to be freed */
    if (lalMallocTotal < n) {
        lalRaiseHook(SIGSEGV, "%s error: lalMallocTotal too small\n",
                     func);
        return NULL;
    }

    /* repad the memory */
    for (i = keep ? n : 0; i < padFactor * n; ++i) {
        s[i + prefix] = (char) (i ^ repadding);
    }

    q[0] = -1;  /* set negative to detect duplicate frees */
    q[1] = ~magic;

    pthread_mutex_lock(&mut);
    lalMallocTotal -= n;
    pthread_mutex_unlock(&mut);

    return q;
}


static void *PushAlloc(void *p, size_t n, const char *file, int line)
{
    struct allocNode *newnode;
    if (!(lalDebugLevel & LALMEMTRKBIT)) {
        return p;
    }
    if (!p) {
        return NULL;
    }
    if (!(newnode = malloc(sizeof(*newnode)))) {
        return NULL;
    }
    pthread_mutex_lock(&mut);
    newnode->addr = p;
    newnode->size = n;
    newnode->file = file;
    newnode->line = line;
    if (!AllocHashTblAdd(newnode)) {
        free(newnode);
        return NULL;
    }
    pthread_mutex_unlock(&mut);
    return p;
}


static void *PopAlloc(void *p, const char *func)
{
    if (!(lalDebugLevel & LALMEMTRKBIT)) {
        return p;
    }
    if (!p) {
        return NULL;
    }
    pthread_mutex_lock(&mut);
    struct allocNode key = { .addr = p };
    struct allocNode *node = AllocHashTblExtract(&key);
    if (node == NULL) {
        pthread_mutex_unlock(&mut);
        lalRaiseHook(SIGSEGV, "%s error: alloc %p not found\n", func, p);
        return NULL;
    }
    free(node);
    pthread_mutex_unlock(&mut);
    return p;
}


static void *ModAlloc(void *p, void *q, size_t n, const char *func,
                      const char *file, int line)
{
    if (!(lalDebugLevel & LALMEMTRKBIT)) {
        return q;
    }
    if (!p || !q) {
        return NULL;
    }
    pthread_mutex_lock(&mut);
    struct allocNode key = { .addr = p };
    struct allocNode *node = AllocHashTblExtract(&key);
    if (node == NULL) {
        pthread_mutex_unlock(&mut);
        lalRaiseHook(SIGSEGV, "%s error: alloc %p not found\n", func, p);
        return NULL;
    }
    node->addr = q;
    node->size = n;
    node->file = file;
    node->line = line;
    if (!AllocHashTblAdd(node)) {
        free(node);
        return NULL;
    }
    pthread_mutex_unlock(&mut);
    return q;
}



void *LALMallocShort(size_t n)
{
    return (lalDebugLevel & LALMEMDBGBIT) ? LALMallocLong(n, "unknown", -1) : malloc(n);
}



void *LALMallocLong(size_t n, const char *file, int line)
{
    void *p;
    void *q;

    if (!(lalDebugLevel & LALMEMDBGBIT)) {
        return malloc(n);
    }

    p = malloc(allocsz(n));
    q = PushAlloc(PadAlloc(p, n, 0, "LALMalloc"), n, file, line);
    lalMemDbgPtr = lalMemDbgRetPtr = q;
    lalIsMemDbgPtr = lalIsMemDbgRetPtr = (lalMemDbgRetPtr == lalMemDbgUsrPtr);
    if (!q) {
        XLALPrintError("LALMalloc: failed to allocate %zd bytes of memory\n", n);
        XLALPrintError("LALMalloc: %zd bytes of memory already allocated\n", lalMallocTotal);
        if (lalDebugLevel & LALMEMINFOBIT) {
            XLALPrintError("LALMalloc meminfo: out of memory\n");
        }
        if (p) {
            free(p);
        }
    }
    return q;
}



void *LALCallocShort(size_t m, size_t n)
{
    return (lalDebugLevel & LALMEMDBGBIT) ? LALCallocLong(m, n, "unknown", -1) :
 calloc(m, n);
}



void *LALCallocLong(size_t m, size_t n, const char *file, int line)
{
    size_t sz;
    void *p;
    void *q;

    if (!(lalDebugLevel & LALMEMDBGBIT)) {
        return calloc(m, n);
    }

    sz = m * n;
    p = malloc(allocsz(sz));
    q = PushAlloc(PadAlloc(p, sz, 1, "LALCalloc"), sz, file, line);
    lalMemDbgPtr = lalMemDbgRetPtr = q;
    lalIsMemDbgPtr = lalIsMemDbgRetPtr = (lalMemDbgRetPtr == lalMemDbgUsrPtr);
    if (!q) {
        XLALPrintError("LALMalloc: failed to allocate %zd bytes of memory\n", n);
        XLALPrintError("LALMalloc: %zd bytes of memory already allocated\n", lalMallocTotal);
        if (lalDebugLevel & LALMEMINFOBIT) {
            XLALPrintError("LALCalloc meminfo: out of memory\n");
        }
        if (p) {
            free(p);
        }
    }
    return q ? memset(q, 0, sz) : NULL;
}



void *LALReallocShort(void *p, size_t n)
{
    return (lalDebugLevel & LALMEMDBGBIT) ? LALReallocLong(p, n, "unknown", -1): realloc(p, n);
}



void *LALReallocLong(void *q, size_t n, const char *file, const int line)
{
    void *p;
    if (!(lalDebugLevel & LALMEMDBGBIT)) {
        return realloc(q, n);
    }

    lalMemDbgPtr = lalMemDbgArgPtr = q;
    lalIsMemDbgPtr = lalIsMemDbgArgPtr = (lalMemDbgArgPtr == lalMemDbgUsrPtr);
    if (!q) {
        p = malloc(allocsz(n));
        q = PushAlloc(PadAlloc(p, n, 0, "LALRealloc"), n, file, line);
        if (!q) {
            XLALPrintError("LALMalloc: failed to allocate %zd bytes of memory\n", n);
            XLALPrintError("LALMalloc: %zd bytes of memory already allocated\n", lalMallocTotal);
            if (lalDebugLevel & LALMEMINFOBIT) {
                XLALPrintError("LALRealloc meminfo: out of memory\n");
            }
            if (p) {
                free(p);
            }
        }
        return q;
    }

    if (!n) {
        p = UnPadAlloc(PopAlloc(q, "LALRealloc"), 0, "LALRealloc");
        if (p) {
            free(p);
        }
        return NULL;
    }

    p = UnPadAlloc(q, 1, "LALRealloc");
    if (!p) {
        return NULL;
    }

    q = ModAlloc(q, PadAlloc(realloc(p, allocsz(n)), n, 1, "LALRealloc"), n, "LALRealloc", file, line);
    lalMemDbgPtr = lalMemDbgRetPtr = q;
    lalIsMemDbgPtr = lalIsMemDbgRetPtr = (lalMemDbgRetPtr == lalMemDbgUsrPtr);

    return q;
}



void LALFree(void *q)
{
    void *p;
    if (q == NULL)
        return;
    if (!(lalDebugLevel & LALMEMDBGBIT)) {
        free(q);
        return;
    }
    lalMemDbgPtr = lalMemDbgArgPtr = q;
    lalIsMemDbgPtr = lalIsMemDbgArgPtr = (lalMemDbgArgPtr == lalMemDbgUsrPtr);
    p = UnPadAlloc(PopAlloc(q, "LALFree"), 0, "LALFree");
    if (p) {
        free(p);
    }
    return;
}



void LALCheckMemoryLeaks(void)
{
    int leak = 0;
    if (!(lalDebugLevel & LALMEMDBGBIT)) {
        return;
    }

    /* alloc_data_len should be zero */
    if ((lalDebugLevel & LALMEMTRKBIT) && alloc_data_len > 0) {
        XLALPrintError("LALCheckMemoryLeaks: allocation list\n");
        for (int k = 0; k < alloc_data_len; ++k) {
            if (alloc_data[k] != NULL && alloc_data[k] != DEL) {
                XLALPrintError("%p: %zu bytes (%s:%d)\n", alloc_data[k]->addr,
                               alloc_data[k]->size, alloc_data[k]->file,
                               alloc_data[k]->line);
            }
        }
        leak = 1;
    }

    /* lalMallocTotal and alloc_n should be zero */
    if ((lalDebugLevel & LALMEMPADBIT) && (lalMallocTotal || alloc_n)) {
        XLALPrintError("LALCheckMemoryLeaks: %d allocs, %zd bytes\n", alloc_n, lalMallocTotal);
        leak = 1;
    }

    if (leak) {
        lalRaiseHook(SIGSEGV, "LALCheckMemoryLeaks: memory leak\n");
    } else if (lalDebugLevel & LALMEMINFOBIT) {
        XLALPrintError
            ("LALCheckMemoryLeaks meminfo: no memory leaks detected\n");
    }

    return;
}

#else

void (LALCheckMemoryLeaks)(void) { return; }

#endif /* ! defined NDEBUG */
