/************************************ <lalVerbatim file="LALMallocCV">
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{LALMalloc.c}}
\label{ss:LALMalloc.c}

LAL memory allocation routines.

\subsubsection*{Macros}
\begin{verbatim}
#if defined NDEBUG || defined LAL_NDEBUG

#define LALMalloc                           malloc
#define LALMallocShort                      malloc
#define LALMallocLong( n, file, line )      malloc( n )
#define LALCalloc                           calloc
#define LALCallocShort                      calloc
#define LALCallocLong( m, n, file, line )   calloc( m, n )
#define LALRealloc                          realloc
#define LALReallocShort                     realloc
#define LALReallocLong( p, n, file, line )  realloc( p, n )
#define LALFree                             free
#define LALCheckMemoryLeaks()

#else

#define LALMalloc( n )      LALMallocLong( n, __FILE__, __LINE__ )
#define LALCalloc( m, n )   LALCallocLong( m, n, __FILE__, __LINE__ )
#define LALRealloc( p, n )  LALReallocLong( p, n, __FILE__, __LINE__ )

#endif
\end{verbatim}
\idx[Macro]{LALMalloc()}
\idx[Macro]{LALCalloc()}
\idx[Macro]{LALRealloc()}

\subsubsection*{Prototypes}
\input{LALMallocCP}
\idx{LALMallocLong()}
\idx{LALCallocLong()}
\idx{LALReallocLong()}
\idx{LALMallocShort()}
\idx{LALCallocShort()}
\idx{LALReallocShort()}
\idx{LALFree()}
\idx{LALCheckMemoryLeaks()}

\subsubsection*{Description}

These functions are the LAL replacements for \verb+malloc()+, \verb+calloc()+,
\verb+realloc()+, and \verb+free()+, with extra functionality to check for
memory leaks (i.e.\ unfreed memory and segmentation violations).
The \verb+LALMallocLong()+, \verb+LALCallocLong()+, and \verb+LALReallocLong()+
functions have two extra arguments giving the file name and line number of the
calling statement; \verb+LALMallocShort()+, \verb+LALCallocShort()+, and
\verb+LALReallocShort()+ do not have these extra arguments, and are merely
call the corresponding long alloc functions with a file name of
\verb+"unknown"+ and a line number of \verb+-1+ (they are useful if you want
to replace hooks to \verb+malloc()+, \verb+calloc()+, and \verb+realloc()+ of
an external package that provides suitable hooks).  \verb+LALMalloc()+,
\verb+LALCalloc()+, and \verb+LALRealloc()+ are actually macros which call the
functions \verb+LALMallocLong()+, \verb+LALCallocLong()+, and
\verb+LALReallocLong+ with the appropriate file name and line number
information.  In practice, it is almost sufficient to use \verb+LALMalloc()+,
\verb+LALCalloc()+, and \verb+LALRealloc()+ as you would \verb+malloc()+,
\verb+calloc()+, and \verb+realloc()+.

Any time an object is freed, \verb+LALFree()+ checks to make sure that the
memory bounds were not over-written, and that the memory address is valid.  The
function \verb+LALCheckMemoryLeaks()+ is to be called at the end of a program
when all the allocated memory is expected to have been freed.  If there is
memory that has been allocated but not freed then this routine reports an
error.  Whenever a memory leak is detected, the routines raise a segmentation
violation signal \verb+SIGSEGV+.  (The signal is raised using the signal
raising hook \verb+lalRaiseHook+, which can be reset to a different handler if
desired.)

Memory leak detection adds significant computational overhead to a
program.  It also requires the use of static memory, making the code
non-thread-safe (but it can be made posix-thread-safe using the
\verb+--enable-pthread-lock+ configure option).  Production code should
suppress memory leak detection at runtime by setting the global
\verb+lalDebugLevel+ equal to zero or by setting the \verb+LALNMEMDBG+ bit of
\verb+lalDebugLevel+, or at compile time by compiling all modules with the
\verb+NDEBUG+ flag set or by using the \verb+--disable-debug+ configure option.
This causes \verb+LALCheckMemoryLeaks()+ to do nothing, and the other functions
to revert to their standard C counterparts.  In addition, you can turn off
individual components of the memory debugging tools.  Setting the
\verb+LALNMEMPAD+ bit of \verb+lalDebugLevel+ prevents the allocation routines
from ``padding out'' the arrays in an effort to detect buffer overflows.
Setting the \verb+LALNMEMTRK+  bit of \verb+lalDebugLevel+ prevents tracking
the allocations/frees.  Setting the \verb+LALMEMINFO+ bit of
\verb+lalDebugLevel+ produces copious output describing each memory allocation
and deallocation.

\subsubsection*{Algorithm}

When buffer overflow detection is active, \verb+LALMalloc()+ allocates, in
addition to the requested memory, storage at the beginning of the object where
a magic number and the size of the object is recorded, and padding at the end
of the object.  The number of allocations and the total size of allocated
memory are stored in static memory.  When \verb+LALFree()+ is executed, the
padding at the end of the object and the magic number are examined to see if
the bounds of the object were over-written.  The total number of allocations
and the total memory allocated are decreased.  \verb+LALCheckMemoryLeaks()+ is
called when all memory should have been freed.  If the number of allocations or
the total memory allocated is not zero, this routine reports an error.

When memory tracking is active, \verb+LALMalloc()+ keeps a linked list
containing information about each allocation: the memory address, the size of
the allocation, and the file name and line number of the calling statement.
Subsequent calls to \verb+LALFree()+ make sure that the address to be freed was
correctly allocated.  In addition, in the case of a memory leak in which some
memory that was allocated was not freed, \verb+LALCheckMemoryLeaks()+ prints a
list of all allocations and the information about the allocations.

When any of these routines encounter an error, they will issue an error message
using \verb+LALPrintError()+ and will raise a \verb+SIGSEGV+ signal, which will
normally cause execution to terminate.  The signal is raised using the hook
\verb+lalRaiseHook+, which can be set to perform a different action if desired.

These routines also issue status messages indicating how much memory is being
allocated or freed with each function call.  These memory information messages
are considered a distinct class of status message, and can be activated or
suppressed independently of other status messages.  See the discussion in
\verb+LALStatusMacros.h+.

When \verb+lalDebugLevel+ is set to zero or the \verb+LALNMEMDBG+ bit is set,
or when compiled with the \verb+NDEBUG+ flag set, these functions revert to
their standard system versions, and \verb+LALCheckMemoryLeaks()+ does nothing.

\subsubsection*{Uses}

\begin{verbatim}
lalDebugLevel
lalRaiseHook
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

#include <lal/LALConfig.h>

#if ! defined NDEBUG && ! defined LAL_NDEBUG

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#else
#define pthread_mutex_lock( pmut )
#define pthread_mutex_unlock( pmut )
#endif

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALMalloc.h>

NRCSID( LALMALLOCC, "$Id$" );

enum              { nprefix   = 2 };
static const size_t prefix    = nprefix * sizeof( size_t );
static const size_t padFactor = 2;
static const size_t padding   = 0xDeadBeef;
static const size_t repadding = 0xBeefDead;
static const size_t magic     = 0xABadCafe;

#define allocsz( n ) ( lalDebugLevel & LALNMEMPAD ? n : padFactor * n + prefix )

static struct allocNode
{
  void             *addr;
  size_t            size;
  const char       *file;
  int               line;
  struct allocNode *next;
} *allocList = NULL;
static size_t lalMallocTotal = 0;
static int    lalMallocCount = 0;
extern int    lalDebugLevel;

static int CheckAllocList( void )
{
  int    count = 0;
  size_t total = 0;
  struct allocNode *node = allocList;
  while ( node )
  {
    ++count;
    total += node->size;
    node = node->next;
  }
  return count == lalMallocCount && total == lalMallocTotal;
}

static void *PadAlloc( size_t *p, size_t n, int keep, const char *func )
{
  size_t i;

  if ( lalDebugLevel & LALNMEMPAD )
  {
    return p;
  }

  if ( ! p )
  {
    return NULL;
  }

  if ( lalDebugLevel & LALMEMINFO )
  {
    LALPrintError( "%s meminfo: allocating %ld bytes\n", func, n );
  }

  if ( lalDebugLevel & LALWARNING && n == 0 )
  {
    LALPrintError( "%s warning: zero size allocation\n", func );
  }
  
  /* store the size in a known position */
  p[0] = n;
  p[1] = magic;

  /* pad the memory */
  for ( i = keep ? n : 0; i < padFactor * n; ++i )
  {
    ((char *)p)[i + prefix] = (char) (i ^ padding);
  }

  pthread_mutex_lock( &mut );
  lalMallocTotal += n;
  ++lalMallocCount;
  pthread_mutex_unlock( &mut );

  return (void *)(((char *)p) + prefix);
}


static void *UnPadAlloc( void *p, int keep, const char *func )
{
  size_t  n;
  size_t  i;
  size_t *q;
  char   *s;

  if ( lalDebugLevel & LALNMEMPAD )
  {
    return p;
  }

  if ( ! p || ! ( q = ((size_t *) p ) - nprefix ) )
  {
    lalRaiseHook( SIGSEGV, "%s error: tried to free NULL pointer\n", func );
    return NULL;
  }

  n = q[0];
  s = (char *)q;

  if ( lalDebugLevel & LALMEMINFO )
  {
    LALPrintError( "%s meminfo: freeing %ld bytes\n", func, n );
  }

  if ( lalDebugLevel & LALWARNING && n == 0 )
  {
    LALPrintError( "%s warning: tried to free a freed pointer\n", func );
  }
          
  if ( q[1] != magic )
  {
    lalRaiseHook( SIGSEGV, "%s error: wrong magic\n", func );
    return NULL;
  }

  if ( ( (int) n ) < 0 )
  {
    lalRaiseHook( SIGSEGV, "%s error: corrupt size descriptor\n", func );
    return NULL;
  }
          
  /* check for writing past end of array: */
  for ( i = n; i < padFactor * n; ++i )
  {
    if ( s[i + prefix] != (char) (i ^ padding) )
    {
      lalRaiseHook( SIGSEGV, "%s error: array bounds overwritten\n"
          "Byte %ld past end of array has changed\n", func, i - n + 1 );
      return NULL;
    }
  }
          
  /* see if there is enough allocated memory to be freed */
  if ( lalMallocTotal < n )
  {
    lalRaiseHook( SIGSEGV, "%s error: lalMallocTotal too small\n", func );
    return NULL;
  }

  /* repad the memory */
  for ( i = keep ? n : 0; i < padFactor * n; ++i )
  {
    s[i + prefix] = (char) (i ^ repadding);
  }

  q[0] = -1; /* set negative to detect duplicate frees */
  q[1] = ~magic;

  pthread_mutex_lock( &mut );
  lalMallocTotal -= n;
  --lalMallocCount;
  pthread_mutex_unlock( &mut );

  return q;
}


static void *PushAlloc( void *p, size_t n, const char *file, int line )
{
  struct allocNode *new;
  if ( lalDebugLevel & LALNMEMTRK )
  {
    return p;
  }
  if ( ! p )
  {
    return NULL;
  }
  if ( ! ( new = malloc( sizeof( *new ) ) ) )
  {
    return NULL;
  }
  pthread_mutex_lock( &mut );
  new->addr = p;
  new->size = n;
  new->file = file;
  new->line = line;
  new->next = allocList;
  allocList = new;
  pthread_mutex_unlock( &mut );
  return p;
}


static void *PopAlloc( void *p, const char *func )
{
  struct allocNode *node;
  if ( lalDebugLevel & LALNMEMTRK )
  {
    return p;
  }
  if ( ! p )
  {
    return NULL;
  }
  pthread_mutex_lock( &mut );
  if ( ! ( node = allocList ) ) /* empty allocation list */
  {
    lalRaiseHook( SIGSEGV, "%s error: alloc not found\n", func );
    return NULL;
  }
  if ( p == node->addr ) /* free the top of the list */
  {
    allocList = node->next;
    free( node );
  }
  else /* free somewhere within the list */
  {
    while ( node->next && p != node->next->addr )
    {
      node = node->next;
    }
    if ( ! node->next ) /* bottom of list reached */
    {
      lalRaiseHook( SIGSEGV, "%s error: alloc not found\n", func );
      return NULL;
    }
    else /* found the alloc */
    {
      struct allocNode *tmp = node->next;
      node->next = node->next->next;
      free( tmp );
    }
  }
  pthread_mutex_unlock( &mut );
  return p;
}


static void *ModAlloc( void *p, void *q, size_t n, const char *func,
    const char *file, int line )
{
  struct allocNode *node;
  if ( lalDebugLevel & LALNMEMTRK )
  {
    return q;
  }
  if ( ! p || ! q )
  {
    return NULL;
  }
  pthread_mutex_lock( &mut );
  if ( ! ( node = allocList ) ) /* empty allocation list */
  {
    lalRaiseHook( SIGSEGV, "%s error: alloc not found\n", func );
    return NULL;
  }
  while ( p != node->addr )
  {
    if ( ! ( node = node->next ) )
    {
      lalRaiseHook( SIGSEGV, "%s error: alloc not found\n", func );
      return NULL;
    }
  }
  node->addr = q;
  node->size = n;
  node->file = file;
  node->line = line;
  pthread_mutex_unlock( &mut );
  return q;
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALMallocShort( size_t n )
{ /* </lalVerbatim> */
  return ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG ) ? malloc( n ) :
    LALMallocLong( n, "unknown", -1 );
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALMallocLong( size_t n, const char *file, int line )
{ /* </lalVerbatim> */
  void *p;
  void *q;

  if ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG )
  {
    return malloc( n );
  }

  p = malloc( allocsz( n ) );
  q = PushAlloc( PadAlloc( p, n, 0, "LALMalloc" ), n, file, line );
  if ( ! q )
  {
    if ( lalDebugLevel & LALMEMINFO )
    {
      LALPrintError( "LALMalloc meminfo: out of memory\n" );
    }
    if ( p )
    {
      free( p );
    }
  }    
  return q;
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALCallocShort( size_t m, size_t n )
{ /* </lalVerbatim> */
  return ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG ) ? calloc( m, n ) :
    LALCallocLong( m, n, "unknown", -1 );
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALCallocLong( size_t m, size_t n, const char *file, int line )
{ /* </lalVerbatim> */
  size_t sz;
  void *p;
  void *q;

  if ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG )
  {
    return calloc( m, n );
  }

  sz = m * n;
  p  = malloc( allocsz( sz ) );
  q  = PushAlloc( PadAlloc( p, sz, 1, "LALCalloc" ), sz, file, line );
  if ( ! q )
  {
    if ( lalDebugLevel & LALMEMINFO )
    {
      LALPrintError( "LALCalloc meminfo: out of memory\n" );
    }
    if ( p )
    {
      free( p );
    }
  }    
  return q ? memset( q, 0, sz ) : NULL;
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALReallocShort( void *p, size_t n )
{ /* </lalVerbatim> */
  return ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG ) ? realloc( p, n ) :
    LALReallocLong( p, n, "unknown", -1 );
}


/* <lalVerbatim file="LALMallocCP"> */
void *
LALReallocLong( void *q, size_t n, const char *file, const int line )
{ /* </lalVerbatim> */
  void *p;
  if ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG )
  {
    return realloc( q, n );
  }

  if ( ! q )
  {
    p = malloc( allocsz( n ) );
    q = PushAlloc( PadAlloc( p, n, 0, "LALRealloc" ), n, file, line );
    if ( ! q )
    {
      if ( lalDebugLevel & LALMEMINFO )
      {
        LALPrintError( "LALRealloc meminfo: out of memory\n" );
      }
      if ( p )
      {
        free( p );
      }
    }    
    return q;
  }

  if ( ! n )
  {
    p = UnPadAlloc( PopAlloc( q, "LALRealloc" ), 0, "LALRealloc" );
    if ( p )
    {
      free( p );
    }
    return NULL;
  }

  p = UnPadAlloc( q, 1, "LALRealloc" );
  if ( ! p )
  {
    return NULL;
  }

  q = ModAlloc( q, PadAlloc( realloc( p, allocsz( n ) ), n, 1, "LALRealloc" ),
      n, "LALRealloc", file, line );

  return q;
}


/* <lalVerbatim file="LALMallocCP"> */
void
LALFree( void *q )
{ /* </lalVerbatim> */
  void *p;
  if ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG )
  {
    free( q );
    return;
  }
  p = UnPadAlloc( PopAlloc( q, "LALFree" ), 0, "LALFree" );
  if ( p )
  {
    free( p );
  }
  return;
}


/* <lalVerbatim file="LALMallocCP"> */
void
LALCheckMemoryLeaks( void )
{ /* </lalVerbatim> */
  int leak = 0;
  if ( ! lalDebugLevel || lalDebugLevel & LALNMEMDBG )
  {
    return;
  }

  /* allocList should be NULL */
  if ( ! ( lalDebugLevel & LALNMEMTRK ) && allocList )
  {
    struct allocNode *node = allocList;
    LALPrintError( "LALCheckMemoryLeaks: allocation list\n" );
    while ( node )
    {
      LALPrintError( "%p: %lu bytes (%s:%d)\n", node->addr,
          (unsigned long)node->size, node->file, node->line );
      node = node->next;
    }
    leak = 1;
  }

  /* lalMallocTotal and lalMallocCount should be zero */
  if ( !( lalDebugLevel & LALNMEMPAD ) && ( lalMallocTotal || lalMallocCount ) )
  {
    LALPrintError( "LALCheckMemoryLeaks: lalMallocCount = %d allocs, "
        "lalMallocTotal = %ld bytes\n", lalMallocCount, (long)lalMallocTotal );
    leak = 1;
  }

  if ( leak )
  {
    lalRaiseHook( SIGSEGV, "LALCheckMemoryLeaks: memory leak\n" );
  }
  else if ( lalDebugLevel & LALMEMINFO )
  {
    LALPrintError( "LALCheckMemoryLeaks meminfo: no memory leaks detected\n" );
  }

  return;
}

#endif /* ! defined NDEBUG && ! defined LAL_NDEBUG */
