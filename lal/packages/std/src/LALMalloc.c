/*
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

#include <lal/LALConfig.h>
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
 *
 */

/* Only define these if configured to do so */

#ifdef WHATEVER_CONFIG_FFTALIGNED

#define XLAL_TEST_POINTER_ALIGNED( ptr, size, retval )	  \
  if ( ! (ptr) && (size) && (retval) )			  \
       XLAL_ERROR_NULL( XLAL_ENOMEM );                    \
  else (void)(0)


void *XLALAlignedMalloc(size_t n){
  void *p;
  int retval;

  retval = posix_memalign(&p,LAL_MEM_ALIGNMENT,n);
  XLAL_TEST_POINTER_ALIGNED(p,n,retval);
  return p;
}

void *XLALAlignedCalloc(size_t m, size_t n){
  void *p;
  int retval;

  retval = posix_memalign(&p,LAL_MEM_ALIGNMENT,m*n);
  XLAL_TEST_POINTER_ALIGNED(p,m*n,retval);
  p = memset(p,0,m*n);
  return p;
}

void *XLALAlignedRealloc(void *p, size_t n){
  /* We don't really expect reallocing aligned memory to
     be efficient, but we'd like it not to break.  Hence
     we just save a copy and see if we either didn't change
     the pointer, or if we did that the newptr is still
     a multiple of LAL_MEM_ALIGNMENT.  Only if both of these
     fail do we call posix_memalign and memcopy */
  void *saveptr, *newptr;
  int retval;

  saveptr = p;
  newptr = realloc(p,n);
  XLAL_TEST_POINTER(newptr,n);
  /* If n = 0 realloc should behave like free */
  if (n) {
    if (saveptr == newptr) {
      /* If the pointer didn't change we're good */
      return newptr;
    } else if (! (int) (((uintptr_t) newptr) % LAL_MEM_ALIGNMENT)) {
      /* If it did but is still a multiple of LAL_MEM_ALIGNMENT we're also good.
         Note that realloc will already have freed p == saveptr */
      return newptr;
    } else {
      /* If it changed but is NOT a multiple of LAL_MEM_ALIGNMENT, we've got to
	 try a new allocation and copy over.  We reuse saveptr since it will
	 have been freed by realloc when realloc returned a new location */
      retval = posix_memalign(&saveptr,LAL_MEM_ALIGNMENT,n);
      XLAL_TEST_POINTER_ALIGNED(p,n,retval);
      saveptr = memcpy(saveptr,newptr,n);
      /* Don't need newptr anymore */
      free(newptr);
      return saveptr;
    }
  } else {
    /* Since n was zero, return our original pointer */
    return p;
  }

}

#endif /* WHATEVER_CONFIG_FFTALIGNED */


/*
 *
 * LAL Routines... only if compiled with debugging enabled.
 * (otherwise the LALMalloc-family reverts to the standard malloc-family).
 *
 */


#if ! defined NDEBUG && ! defined LAL_NDEBUG

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

static struct allocNode {
    void *addr;
    size_t size;
    const char *file;
    int line;
    struct allocNode *next;
} *allocList = NULL;
static int lalMallocCount = 0;

/* need this to turn off gcc warnings about unused functions */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Useful function for debugging */
/* Checks to make sure alloc list is OK */
/* Returns 0 if list is corrupted; 1 if list is OK */
UNUSED static int CheckAllocList(void)
{
    int count = 0;
    size_t total = 0;
    struct allocNode *node = allocList;
    while (node) {
        ++count;
        total += node->size;
        node = node->next;
    }
    return count == lalMallocCount && total == lalMallocTotal;
}

/* Useful function for debugging */
/* Finds the node of the alloc list previous to the desired alloc */
/* Returns NULL if not found or if the alloc node is the head */
UNUSED static struct allocNode *FindPrevAlloc(void *p)
{
    struct allocNode *node = allocList;
    if (p == node->addr)        /* top of list */
        return NULL;
    /* scroll through list to find node before the alloc */
    while (node->next)
        if (node->next->addr == p)
            return node;
        else
            node = node->next;
    return NULL;
}

/* Useful function for debugging */
/* Finds the node of the alloc list for the desired alloc */
/* Returns NULL if not found  */
UNUSED static struct allocNode *FindAlloc(void *p)
{
    struct allocNode *node = allocList;
    while (node)
        if (p == node->addr)
            return node;
        else
            node = node->next;
    return NULL;
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
        LALPrintError("%s meminfo: allocating %ld bytes at address %p\n",
                      func, n, p + nprefix);
    }

    if (lalDebugLevel & LALWARNING && n == 0) {
        LALPrintError("%s warning: zero size allocation at address %p\n",
                      func, p + nprefix);
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
    ++lalMallocCount;
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
        LALPrintError("%s meminfo: freeing %ld bytes at address %p\n",
                      func, n, p);
    }

    if (lalDebugLevel & LALWARNING && n == 0) {
        LALPrintError
            ("%s warning: tried to free a freed pointer at address %p\n",
             func, p);
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
    --lalMallocCount;
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
    newnode->next = allocList;
    allocList = newnode;
    pthread_mutex_unlock(&mut);
    return p;
}


static void *PopAlloc(void *p, const char *func)
{
    struct allocNode *node;
    if (!(lalDebugLevel & LALMEMTRKBIT)) {
        return p;
    }
    if (!p) {
        return NULL;
    }
    pthread_mutex_lock(&mut);
    if (!(node = allocList)) {  /* empty allocation list */
        pthread_mutex_unlock(&mut);
        lalRaiseHook(SIGSEGV, "%s error: alloc %p not found\n", func, p);
        return NULL;
    }
    if (p == node->addr) {      /* free the top of the list */
        allocList = node->next;
        free(node);
    } else {    /* free somewhere within the list */

        while (node->next && p != node->next->addr) {
            node = node->next;
        }
        if (!node->next) {      /* bottom of list reached */
            pthread_mutex_unlock(&mut);
            lalRaiseHook(SIGSEGV, "%s error: alloc %p not found\n", func,
                         p);
            return NULL;
        } else {        /* found the alloc */

            struct allocNode *tmp = node->next;
            node->next = node->next->next;
            free(tmp);
        }
    }
    pthread_mutex_unlock(&mut);
    return p;
}


static void *ModAlloc(void *p, void *q, size_t n, const char *func,
                      const char *file, int line)
{
    struct allocNode *node;
    if (!(lalDebugLevel & LALMEMTRKBIT)) {
        return q;
    }
    if (!p || !q) {
        return NULL;
    }
    pthread_mutex_lock(&mut);
    if (!(node = allocList)) {  /* empty allocation list */
        pthread_mutex_unlock(&mut);
        lalRaiseHook(SIGSEGV, "%s error: alloc %p not found\n", func, p);
        return NULL;
    }
    while (p != node->addr) {
        if (!(node = node->next)) {
            pthread_mutex_unlock(&mut);
            lalRaiseHook(SIGSEGV, "%s error: alloc %p not found\n", func,
                         p);
            return NULL;
        }
    }
    node->addr = q;
    node->size = n;
    node->file = file;
    node->line = line;
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
            LALPrintError("LALMalloc meminfo: out of memory\n");
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
            LALPrintError("LALCalloc meminfo: out of memory\n");
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
                LALPrintError("LALRealloc meminfo: out of memory\n");
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

    /* allocList should be NULL */
    if ((lalDebugLevel & LALMEMTRKBIT) && allocList) {
        struct allocNode *node = allocList;
        LALPrintError("LALCheckMemoryLeaks: allocation list\n");
        while (node) {
            LALPrintError("%p: %lu bytes (%s:%d)\n", node->addr,
                          (unsigned long) node->size, node->file,
                          node->line);
            node = node->next;
        }
        leak = 1;
    }

    /* lalMallocTotal and lalMallocCount should be zero */
    if ((lalDebugLevel & LALMEMPADBIT)
        && (lalMallocTotal || lalMallocCount)) {
        LALPrintError("LALCheckMemoryLeaks: lalMallocCount = %d allocs, "
                      "lalMallocTotal = %ld bytes\n", lalMallocCount,
                      (long) lalMallocTotal);
        leak = 1;
    }

    if (leak) {
        lalRaiseHook(SIGSEGV, "LALCheckMemoryLeaks: memory leak\n");
    } else if (lalDebugLevel & LALMEMINFOBIT) {
        LALPrintError
            ("LALCheckMemoryLeaks meminfo: no memory leaks detected\n");
    }

    return;
}

#else

void (LALCheckMemoryLeaks)(void) { return; }

#endif /* ! defined NDEBUG && ! defined LAL_NDEBUG */
