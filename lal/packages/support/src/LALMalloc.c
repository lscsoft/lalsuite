/*----------------------------------------------------------------------- 
 * 
 * File Name: LALMalloc.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _STDLIB_H
#include <stdlib.h>
#ifndef _STDLIB_H
#define _STDLIB_H
#endif
#endif

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _SIGNAL_H
#include <signal.h>
#ifndef _SIGNAL_H
#define _SIGNAL_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _LALERROR_H
#include "LALError.h"
#ifndef _LALERROR_H
#define _LALERROR_H
#endif
#endif

#ifndef _LALMALLOC_H
#include "LALMalloc.h"
#ifndef _LALMALLOC_H
#define _LALMALLOC_H
#endif
#endif

NRCSID (LALMALLOCC, "$Id$");

static const size_t padFactor = 2;
static const size_t padding   = 0xDeadBeef;
static const size_t repadding = 0xBeefDead;
static const size_t magic     = 0xABadCafe;

static size_t lalMallocTotal = 0;
static int    lalMallocCount = 0;
extern int    debuglevel;


void *
LALMalloc (size_t n)
{
  if (debuglevel == 0)
  {
    return malloc (n);
  }
  else
  {
    size_t  prefix = 2*sizeof(size_t);
    char   *p;
    int     i;

    if (debuglevel > 2)
    {
      LALPrintError ("LALMalloc: allocating %ld bytes\n", n);
    }

    /*
    if (n == 0)
    {
      LALPrintError ("LALMalloc: zero size allocation\n");
      raise (SIGSEGV);
      return NULL;
    }
    */
  
    p = (char *) malloc (padFactor*n + prefix);
    if (!p)
    {
      LALPrintError ("LALMalloc: out of memory\n");
      raise (SIGSEGV);
      return NULL;
    }
  
    /* store the size in a known position */
    ((size_t *) p)[0] = n;
    ((size_t *) p)[1] = magic;
    for (i = 0; i < padFactor*n; ++i)
    {
      p[i + prefix] = (char) (i ^ padding);
    }

    if (debuglevel > 2)
    {
      LALPrintError ("LALMalloc: successful allocation\n");
    }

    lalMallocTotal += n;
    ++lalMallocCount;

    /* skip the size we stored previously */
    return (void *) (p + prefix);
  }
}


void
LALFree (void *p)
{
  if (debuglevel == 0)
  {
    free (p);
    return;
  }
  else
  {
    size_t  prefix = 2*sizeof(size_t);
    char   *q      = ((char *) p) - prefix;

    if (!p)
    {
      LALPrintError ("LALFree: tried to free NULL pointer\n");
      raise (SIGSEGV);
      return;
    }

    if (!q)
    {
      LALPrintError ("LALFree: tried to free NULL+TWOINTS pointer\n");
      raise (SIGSEGV);
      return;
    }

    {
      size_t n       = ((size_t *) q)[0];
      size_t myMagic = ((size_t *) q)[1];
      size_t i;
          
      if (debuglevel > 2)
      {
        LALPrintError ("LALFree: freeing %ld bytes\n", n);
      }
          
      /*
      if (n == 0)
      {
        LALPrintError ("LALFree: tried to free a freed pointer\n");
        raise (SIGSEGV);
        return;
      }
      */
          
      if (myMagic != magic)
      {
        LALPrintError ("LALFree: wrong magic\n");
        raise (SIGSEGV);
        return;
      }

      if (n < 0)
      {
        LALPrintError ("LALFree: corrupt size descriptor\n");
        raise (SIGSEGV);
        return;
      }
          
      /* check for writing past end of array: */
      for (i = n; i < padFactor*n; ++i)
      {
        if (q[i + prefix] != (char) (i ^ padding))
        {
          LALPrintError ("LALFree: array bounds overwritten\n");
          LALPrintError ("Byte %ld past end of array has changed\n",
                   i - n + 1);
          raise (SIGSEGV);
          return;
        }
      }
          
      /* see if there is enough allocated memory to be freed */
      if (lalMallocTotal < n)
      {
        LALPrintError ("LALFree: lalMallocTotal too small\n");
        raise (SIGSEGV);
        return;
      }

      /* repad the memory */
      for (i = 0; i < padFactor*n; ++i)
      {
        q[i + prefix] = (char) (i ^ repadding);
      }

      if (debuglevel > 2)
      {
        LALPrintError ("LALFree: successful freeing\n");
      }
          
      *((size_t *) q) = 0; /* set to zero to detect duplicate frees */
      ((size_t *) q)[1] = ~magic;

      lalMallocTotal -= n;
      --lalMallocCount;
      free (q);
    }
  }
}


void *
LALCalloc (size_t m, size_t n)
{
  if (debuglevel == 0)
  {
    return calloc (m, n);
  }
  else
  {
    void *result = LALMalloc (m*n);

    if (result)
    {
      memset (result, 0, m*n);
    }

    return result;
  }
}


void *
LALRealloc (void *p, size_t n)
{
  if (debuglevel == 0)
  {
    return realloc (p, n);
  }
  else if (!p)
  {
    return LALMalloc (n);
  }
  else if (!n)
  {
    LALFree (p);
    return NULL;
  }
  else
  {
    size_t  prefix = 2*sizeof(size_t);
    char   *q = ((char *) p) - prefix;
    size_t  m = ((size_t *) q)[0];      /* size of old array */

    if (m == n) /* no resizing necessary! */
    {
      return p;
    }
    else
    { /* create new vector, copy old to new, and free old */
      void *new = LALMalloc (n);

      if (new)
      {
        memcpy (new, p, m < n ? m : n);
        LALFree (p);
      }

      return new;
    }
  }
}


void
LALCheckMemoryLeaks (void)
{
  if (debuglevel == 0)
  {
    return;
  }

  if (lalMallocTotal || lalMallocCount)
  {
    LALPrintError ("LALCheckMemoryLeaks: memory leak\n");
    LALPrintError ("lalMallocCount = %d allocs\n", lalMallocCount);
    LALPrintError ("lalMallocTotal = %ld bytes\n", (long)lalMallocTotal);
    raise (SIGSEGV);
    return;
  }

  if (debuglevel > 1)
  {
    LALPrintError ("LALCheckMemoryLeaks: no memory leaks detected\n");
  }

  return;
}
