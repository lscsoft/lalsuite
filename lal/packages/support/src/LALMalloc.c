/************************************ <lalVerbatim file="LALMallocCV">
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{LALMalloc.c}}
\label{ss:LALMalloc.c}

LAL memory allocation routines.

\subsubsection*{Prototypes}
\input{LALMallocCP}
\index{\texttt{LALMalloc()}}
\index{\texttt{LALCalloc()}}
\index{\texttt{LALRealloc()}}
\index{\texttt{LALFree()}}
\index{\texttt{LALCheckMemoryLeaks()}}

\subsubsection*{Description}

These functions are the LAL replacements for \verb+malloc()+,
\verb+calloc()+, \verb+realloc()+, and \verb+free()+, with extra
functionality to check for memory leaks (i.e.\ unfreed memory and
segmentation violations).  Any time an object is freed,
\verb+LALFree()+ checks to make sure that the memory bounds were not
over-written.  The function \verb+LALCheckMemoryLeaks()+ is to be
called at the end of a program when all the allocated memory is
expected to have been freed.  If there is memory that has been
allocated but not freed then this routine reports an error.  Whenever
a memory leak is detected, the routines raise a segmentation violation
signal \verb+SIGSEGV+.

Memory leak detection adds significant computational overhead to a
program.  It also requires the use of static memory, making the code
non-thread-safe.  Production code should suppress memory leak
detection at runtime by setting the global \verb+lalDebugLevel+ equal to
zero, or at compile time by compiling all modules with the
\verb+NDEBUG+ flag set.  This causes \verb+LALCheckMemoryLeaks()+ to
do nothing, and the other functions to revert to their standard C
counterparts.

\subsubsection*{Algorithm}

When memory leak detection is active, \verb+LALMalloc()+ allocates, in
addition to the requested memory, storage at the beginning of the
object where a magic number and the size of the object is recorded,
and padding at the end of the object.  The number of allocations and
the total size of allocated memory are stored in static memory.  When
\verb+LALFree()+ is executed, the padding at the end of the object and
the magic number are examined to see if the bounds of the object were
over-written.  The total number of allocations and the total memory
allocated are decreased.  \verb+LALCheckMemoryLeaks()+ is called when
all memory should have been freed.  If the number of allocations or
the total memory allocated is not zero, this routine reports an error.

Also, with leak detection active the code can only access the heap
using \verb+malloc()+ and \verb+free()+.  This makes the routine
\verb+LALRealloc()+ particularly clumsy.

When any of these routines encounter an error, they will issue an
error message using \verb+LALPrintError()+ and will raise a
\verb+SIGSEGV+ signal, which will normally cause execution to
terminate.

These routines also issue status messages indicating how much memory
is being allocated or freed with each function call.  These memory
information messages are considered a distinct class of status
message, and can be activated or suppressed independently of other
status messages.  See the discussion in \verb+LALStatusMacros.h+.

When \verb+lalDebugLevel+ is set to zero, or when compiled with the
\verb+NDEBUG+ flag set, these functions revert to their standard
system versions, and \verb+LALCheckMemoryLeaks()+ does nothing.

\subsubsection*{Uses}

\begin{verbatim}
lalDebugLevel
LALPrintError()
\end{verbatim}

\subsubsection*{Notes}

Memory leak detection only occurs when \verb+lalDebugLevel+$\neq0$.  To
turn on leak detection independent of error reporting, simply switch
on the most-significant bit of \verb+lalDebugLevel+, which is reserved
not to be associated with any type of status message.  See the
discussion in \verb+LALStatusMacros.h+ for more information about
\verb+lalDebugLevel+.

It is assumed that pointers of type \verb+size_t *+ have the most restrictive
alignment.  If this is not true, then this code may not work except in
non-debugging mode.  (It will probably produce bus errors.)

\vfill{\footnotesize\input{LALMallocCV}}

******************************************************* </lalLaTeX> */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <pthread.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALMalloc.h>
#undef LALMalloc
#undef LALFree
#undef LALCalloc
#undef LALRealloc
#undef LALCheckMemoryLeaks

NRCSID( LALMALLOCC, "$Id$" );

static const size_t padFactor = 2;
static const size_t padding   = 0xDeadBeef;
static const size_t repadding = 0xBeefDead;
static const size_t magic     = 0xABadCafe;

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
static size_t lalMallocTotal = 0;
static int    lalMallocCount = 0;
extern int    lalDebugLevel;


/* <lalVerbatim file="LALMallocCP"> */
void *
LALMalloc( size_t n )
{ /* </lalVerbatim> */
  if ( lalDebugLevel == 0 )
  {
    return malloc( n );
  }
  else
  {
    size_t  i;
    size_t *p       = NULL;
    char   *pc      = NULL;
    size_t  nprefix = 2;
    size_t  prefix  = nprefix*sizeof( size_t );
    int     newline = 0; /* need a new line */

    if ( lalDebugLevel & LALMEMINFO )
    {
      newline = 1;
      LALPrintError( "LALMalloc meminfo: allocating %ld bytes", n );
    }

    if ( lalDebugLevel & LALWARNING && n == 0 )
    {
      newline = newline ? LALPrintError( "\n" ),0 : 0;
      LALPrintError( "LALMalloc warning: zero size allocation\n" );
    }
  
    p  = malloc( padFactor*n + prefix );
    pc = (char *) p;
    if ( !p )
    {
      if ( lalDebugLevel & LALERROR )
      {
        newline = newline ? LALPrintError( "\n" ),0 : 0;
        LALPrintError( "LALMalloc error: out of memory\n" );
      }
      raise( SIGSEGV );
      return NULL;
    }
  
    /* store the size in a known position */
    p[0] = n;
    p[1] = magic;
    for ( i = 0; i < padFactor*n; ++i )
    {
      pc[i + prefix] = (char) (i ^ padding);
    }

    if ( lalDebugLevel & LALMEMINFO )
    {
      if ( newline )
      {
        LALPrintError( "... successful allocation\n" );
      }
      else
      {
        LALPrintError( "LALMalloc meminfo: successful allocation\n" );
      }
    }

    
    pthread_mutex_lock( &mut );
    lalMallocTotal += n;
    ++lalMallocCount;
    pthread_mutex_unlock( &mut );

    /* skip the size we stored previously */
    return (void *) (pc + prefix);
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void
LALFree( void *p )
{ /* </lalVerbatim> */
  if ( lalDebugLevel == 0 )
  {
    free( p );
    return;
  }
  else
  {
    size_t  nprefix = 2;
    size_t  prefix  = nprefix*sizeof( size_t );
    size_t *q       = ((size_t *) p) - nprefix;
    int     newline = 0; /* need a new line */

    if ( !p )
    {
      if ( lalDebugLevel & LALERROR )
      {
        LALPrintError( "LALFree error: tried to free NULL pointer\n" );
      }
      raise( SIGSEGV );
      return;
    }

    if ( !q )
    {
      if ( lalDebugLevel & LALERROR )
      {
        LALPrintError( "LALFree error: tried to free NULL+TWOINTS pointer\n" );
      }
      raise( SIGSEGV );
      return;
    }

    {
      size_t  n       = q[0];
      size_t  myMagic = q[1];
      char   *qc      = (char *) q;
      size_t  i;
          
      if ( lalDebugLevel & LALMEMINFO )
      {
        newline = 1;
        LALPrintError( "LALFree meminfo: freeing %ld bytes", n );
      }
          
      if ( lalDebugLevel & LALWARNING && n == 0)
      {
        newline = newline ? LALPrintError( "\n" ),0 : 0;
        LALPrintError( "LALFree warning: tried to free a freed pointer\n" );
      }
          
      if ( myMagic != magic )
      {
        if ( lalDebugLevel & LALERROR )
        {
          newline = newline ? LALPrintError( "\n" ),0 : 0;
          LALPrintError( "LALFree error: wrong magic\n" );
        }
        raise( SIGSEGV );
        return;
      }

      if ( ( (int) n ) < 0 )
      {
        if ( lalDebugLevel & LALERROR )
        {
          newline = newline ? LALPrintError( "\n" ),0 : 0;
          LALPrintError( "LALFree error: corrupt size descriptor\n" );
        }
        raise( SIGSEGV );
        return;
      }
          
      /* check for writing past end of array: */
      for ( i = n; i < padFactor*n; ++i )
      {
        if ( qc[i + prefix] != (char) (i ^ padding) )
        {
          if ( lalDebugLevel & LALERROR )
          {
            newline = newline ? LALPrintError( "\n" ),0 : 0;
            LALPrintError( "LALFree error: array bounds overwritten\n" );
            LALPrintError( "Byte %ld past end of array has changed\n",
                           i - n + 1 );
          }
          raise( SIGSEGV );
          return;
        }
      }
          
      /* see if there is enough allocated memory to be freed */
      if ( lalMallocTotal < n )
      {
        if ( lalDebugLevel & LALERROR )
        {
          newline = newline ? LALPrintError( "\n" ),0 : 0;
          LALPrintError( "LALFree error: lalMallocTotal too small\n" );
        }
        raise( SIGSEGV );
        return;
      }

      /* repad the memory */
      for ( i = 0; i < padFactor*n; ++i )
      {
        qc[i + prefix] = (char) (i ^ repadding);
      }

      if ( lalDebugLevel & LALMEMINFO )
      {
        if ( newline )
        {
          LALPrintError( "... successful freeing\n" );
        }
        else
        {
          LALPrintError( "LALFree meminfo: successful freeing\n" );
        }
      }
          
      q[0] = 0; /* set to zero to detect duplicate frees */
      q[1] = ~magic;

      pthread_mutex_lock( &mut );
      lalMallocTotal -= n;
      --lalMallocCount;
      pthread_mutex_unlock( &mut );

      free( q );
    }
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALCalloc( size_t m, size_t n )
{ /* </lalVerbatim> */
  if ( lalDebugLevel == 0 )
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
  if ( lalDebugLevel == 0 )
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
    size_t  nprefix = 2;
    size_t *q = ((size_t *) p) - nprefix;
    size_t  m = q[0];      /* size of old array */

    if ( m == n ) /* no resizing necessary! */
    {
      return p;
    }
    else
    { /* create new vector, copy old to new, and free old */
      void *pnew = LALMalloc( n );

      if ( pnew )
      {
        memcpy( pnew, p, m < n ? m : n );
        LALFree( p );
      }

      return pnew;
    }
  }
}


/* <lalVerbatim file="LALMallocCP"> */
void
LALCheckMemoryLeaks( void )
{ /* </lalVerbatim> */
  if ( lalDebugLevel == 0 )
  {
    return;
  }

  if ( lalMallocTotal || lalMallocCount )
  {
    if ( lalDebugLevel & LALERROR || lalDebugLevel & LALWARNING )
    {
      LALPrintError( "LALCheckMemoryLeaks: memory leak\n" );
      LALPrintError( "lalMallocCount = %d allocs\n", lalMallocCount );
      LALPrintError( "lalMallocTotal = %ld bytes\n", (long)lalMallocTotal );
    }
    raise( SIGSEGV );
    return;
  }

  if ( lalDebugLevel & LALMEMINFO )
  {
    LALPrintError( "LALCheckMemoryLeaks meminfo: no memory leaks detected\n" );
  }

  return;
}
