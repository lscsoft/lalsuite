/************************************ <lalVerbatim file="LALMallocCV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALMalloc.c}}
\label{ss:LALMalloc.c}

LAL memory allocation routines.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALMallocCP}
\index{\verb&LALMalloc()&}
\index{\verb&LALCalloc()&}
\index{\verb&LALRealloc()&}
\index{\verb&LALFree()&}
\index{\verb&LALCheckMemoryLeaks()&}

\subsubsection*{Description}

These functions are the LAL replacements for \verb+malloc+, \verb+calloc+,
\verb+realloc+, and \verb+free+.  When \verb+debuglevel+ is zero, they are
the essentially the same as the standard functions but when \verb+debuglevel+
is greater then memory leak detection code is invoked.  Any time an object is
freed, \verb+LALFree+ checks to make sure that the memory bounds were not
over-written.  The function \verb+LALCheckMemoryLeaks+ is to be called at the
end of a program when all the allocated memory is expected to have been freed.
If there is memory that has been allocated but not freed then this routine
reports an error.  Whenever a memory leak is detected, the routines raise a
segmentation violation signal.

Use of these routines with \verb+debuglevel+ greater than zero can cause the
performance of the code to suffer.  Also, these routines use static memory and
are not thread-safe when \verb+debuglevel+ is greater than zero.  Production
code should therefore be run with \verb+debuglevel+ equal to zero.

\subsubsection*{Algorithm}

When \verb+debuglevel+ is greater than zero, \verb+LALMalloc+ allocates, in
addition to the requested memory, storage at the beginning of the object where
a magic number and the size of the object is recorded, and padding at the end
of the object.  The number of allocations and the total size of allocated
memory is stored in static memory.  When \verb+LALFree+ is executed, the
padding at the end of the object and the magic number are examined to see if
the bounds of the object were over-written.  The total number of allocations
and the total memory allocated are decreased.  \verb+LALCheckMemoryLeaks+
is called when all memory should have been freed.  If the number of allocations
or the total memory allocated is not zero, this routine reports an error.

The code only accesses the heap using \verb+malloc+ and \verb+free+.  This
means that the routine \verb+LALRealloc+ is particularly clumsy when run with
\verb+debuglevel+ greater than zero.

However, when \verb+debuglevel+ is set to zero, these functions revert to their
standard system versions, and \verb+LALCheckMemoryLeaks+ does nothing.

\subsubsection*{Uses}

\begin{verbatim}
debuglevel
LALPrintError()
\end{verbatim}

\subsubsection*{Notes}

Use of these routines with \verb+debuglevel+ greater than zero can cause the
performance of the code to suffer.  Also, these routines use static memory and
are not thread-safe when \verb+debuglevel+ is greater than zero.  When memory
leaks are detected with \verb+debuglevel+ greater than zero, a segmentation
violation signal is raised.  Production code should therefore be run with
\verb+debuglevel+ equal to zero.

\vfill{\footnotesize\input{LALMallocCV}}

</lalLaTeX> */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "LALStdlib.h"
#include "LALError.h"
#include "LALMalloc.h"

NRCSID( LALMALLOCC, "$Id$" );

static const size_t padFactor = 2;
static const size_t padding   = 0xDeadBeef;
static const size_t repadding = 0xBeefDead;
static const size_t magic     = 0xABadCafe;

static size_t lalMallocTotal = 0;
static int    lalMallocCount = 0;
extern int    debuglevel;


/* <lalVerbatim file="LALMallocCP"> */
void *
LALMalloc( size_t n )
{ /* </lalVerbatim> */
  if ( debuglevel == 0 )
  {
    return malloc( n );
  }
  else
  {
    size_t  prefix = 2*sizeof( size_t );
    char   *p;
    int     i;

    if ( debuglevel > 2 )
    {
      LALPrintError( "LALMalloc: allocating %ld bytes\n", n );
    }

    /* This should be a warning, not an error.
    if (n == 0)
    {
      LALPrintError( "LALMalloc: zero size allocation\n" );
      raise( SIGSEGV );
      return NULL;
    }
    */
  
    p = (char *) malloc( padFactor*n + prefix );
    if ( !p )
    {
      LALPrintError( "LALMalloc: out of memory\n" );
      raise( SIGSEGV );
      return NULL;
    }
  
    /* store the size in a known position */
    ((size_t *) p)[0] = n;
    ((size_t *) p)[1] = magic;
    for ( i = 0; i < padFactor*n; ++i )
    {
      p[i + prefix] = (char) (i ^ padding);
    }

    if ( debuglevel > 2 )
    {
      LALPrintError( "LALMalloc: successful allocation\n" );
    }

    lalMallocTotal += n;
    ++lalMallocCount;

    /* skip the size we stored previously */
    return (void *) (p + prefix);
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void
LALFree( void *p )
{ /* </lalVerbatim> */
  if ( debuglevel == 0 )
  {
    free( p );
    return;
  }
  else
  {
    size_t  prefix = 2*sizeof( size_t );
    char   *q      = ((char *) p) - prefix;

    if ( !p )
    {
      LALPrintError( "LALFree: tried to free NULL pointer\n" );
      raise( SIGSEGV );
      return;
    }

    if ( !q )
    {
      LALPrintError( "LALFree: tried to free NULL+TWOINTS pointer\n" );
      raise( SIGSEGV );
      return;
    }

    {
      size_t n       = ((size_t *) q)[0];
      size_t myMagic = ((size_t *) q)[1];
      size_t i;
          
      if ( debuglevel > 2 )
      {
        LALPrintError( "LALFree: freeing %ld bytes\n", n );
      }
          
      /* This should be a warning, not an error.
      if (n == 0)
      {
        LALPrintError( "LALFree: tried to free a freed pointer\n" );
        raise( SIGSEGV );
        return;
      }
      */
          
      if ( myMagic != magic )
      {
        LALPrintError( "LALFree: wrong magic\n" );
        raise( SIGSEGV );
        return;
      }

      if ( n < 0 )
      {
        LALPrintError( "LALFree: corrupt size descriptor\n" );
        raise( SIGSEGV );
        return;
      }
          
      /* check for writing past end of array: */
      for ( i = n; i < padFactor*n; ++i )
      {
        if ( q[i + prefix] != (char) (i ^ padding) )
        {
          LALPrintError( "LALFree: array bounds overwritten\n" );
          LALPrintError( "Byte %ld past end of array has changed\n",
                         i - n + 1 );
          raise( SIGSEGV );
          return;
        }
      }
          
      /* see if there is enough allocated memory to be freed */
      if ( lalMallocTotal < n )
      {
        LALPrintError( "LALFree: lalMallocTotal too small\n" );
        raise( SIGSEGV );
        return;
      }

      /* repad the memory */
      for ( i = 0; i < padFactor*n; ++i )
      {
        q[i + prefix] = (char) (i ^ repadding);
      }

      if ( debuglevel > 2 )
      {
        LALPrintError( "LALFree: successful freeing\n" );
      }
          
      *((size_t *) q) = 0; /* set to zero to detect duplicate frees */
      ((size_t *) q)[1] = ~magic;

      lalMallocTotal -= n;
      --lalMallocCount;
      free( q );
    }
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALCalloc( size_t m, size_t n )
{ /* </lalVerbatim> */
  if ( debuglevel == 0 )
  {
    return calloc( m, n );
  }
  else
  {
    void *result = LALMalloc( m*n );

    if ( result )
    {
      memset( result, 0, m*n );
    }

    return result;
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALRealloc( void *p, size_t n )
{ /* </lalVerbatim> */
  if ( debuglevel == 0 )
  {
    return realloc( p, n );
  }
  else if ( !p )
  {
    return LALMalloc( n );
  }
  else if ( !n )
  {
    LALFree( p );
    return NULL;
  }
  else
  {
    size_t  prefix = 2*sizeof( size_t );
    char   *q = ((char *) p) - prefix;
    size_t  m = ((size_t *) q)[0];      /* size of old array */

    if ( m == n ) /* no resizing necessary! */
    {
      return p;
    }
    else
    { /* create new vector, copy old to new, and free old */
      void *new = LALMalloc( n );

      if ( new )
      {
        memcpy( new, p, m < n ? m : n );
        LALFree( p );
      }

      return new;
    }
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void
LALCheckMemoryLeaks( void )
{ /* </lalVerbatim> */
  if ( debuglevel == 0 )
  {
    return;
  }

  if ( lalMallocTotal || lalMallocCount )
  {
    LALPrintError( "LALCheckMemoryLeaks: memory leak\n" );
    LALPrintError( "lalMallocCount = %d allocs\n", lalMallocCount );
    LALPrintError( "lalMallocTotal = %ld bytes\n", (long)lalMallocTotal );
    raise( SIGSEGV );
    return;
  }

  if ( debuglevel > 1 )
  {
    LALPrintError( "LALCheckMemoryLeaks: no memory leaks detected\n" );
  }

  return;
}
