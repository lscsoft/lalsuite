/*
*  Copyright (C) 2007 Jolien Creighton
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


\subsubsection*{Debugging memory leak tips}

Programs should end by calling \verb+LALCheckMemoryLeaks()+.  This will
ensure that all memory that has been allocated has been freed.  Making
sure that all memory allocated is freed is a good idea in order to
make sure (i) that memory isn't being ``lost'' (which may mean that
the computer will run out of memory when the program is run under more
extensive use), (ii) that array bounds are not being exceeded (since this
will usually overwrite the pad area at the end of the array, and this
overwrite is detected when the array is freed).  \verb+LALCheckMemoryLeaks()+
should pass silently---if it doesn't, then there is probably some memory
that has not been freed; \verb+LALCheckMemoryLeaks()+ will give information
about where this memory was allocated.

The most common problem (after forgetting to free some memory) is overwriting
of array bounds.  When this is detected, \verb+LALFree()+ reports the memory
address that was overwritten, as well as the address of the array that
\verb+LALFree()+ attempted to free.  In order to find out where the overwrite
occurs, run the program in the debugger and stop the execution as soon as the
array that is being overwritten has been allocated.  The \verb+LALMalloc+
module has some secret memory debugging tools (for use in debugging only!).
One is the global variable \verb+lalMemDbgUsrPtr+, which is of type
\verb+char *+.  Set this variable to be equal to the memory address
where the overwrite occurs.  Then watch the contents of the variable
to find out where the overwrite occurs.  This is done in \verb+gdb+ using
the commands:
\begin{verbatim}
set var lalMemDbgUsrPtr=0x20f530
watch *lalMemDgbUsrPtr
cont
\end{verbatim}
where \verb+0x20f530+ is the corrupted memory address.  The program
will run until the value of this address is changed, thereby allowing
you to find out where in the program the overwrite occurs.

If you don't know where the memory was allocated, you can locate this
too.  To do so, set \verb+lalMemDbgUsrPtr+ to be the address of the array.
Then, every time \verb+LALMalloc()+ is called, it sets the value of the
global variable \verb+lalIsMemDbgRetPtr+ to be one zero if the array
address produced by \verb+LALMalloc()+ is not the address in
\verb+lalMemDbgUsrPtr+, and one if it is.  Then you can watch the value
of \verb+lalIsMemDbgRetPtr+ in a debugger until it changes to one, which stops
execution at that point.  (Note: it is possible that a given address is
allocated, then freed, the allocated again---you may need to watch
\verb+lalIsMemDbgRetPtr+ for a while.)

Here's an example debugging session: first we run the program, identify
the address of the array whose bounds are being overwritten, and find
out where that array is allocated.
\begin{verbatim}
(gdb) run
LALFree error: array bounds overwritten
Byte 4 past end of array has changed
Corrupted address: 0x1cf530
Array address: 0x1cf528

Program received signal SIGSEGV, Segmentation fault.
0x9001b46c in kill ()
(gdb) list 1,11
1       #include <lal/LALStdlib.h>
2       int main( void )
3       {
4         char *s;
5         lalDebugLevel = 1;
6         s = LALMalloc( 5 );
7         s[8] = 'x';
8         LALFree( s );
9         LALCheckMemoryLeaks();
10        return 0;
11      }
(gdb) break 5
Breakpoint 1 at 0x1b60: file bad.c, line 5.
(gdb) run

Breakpoint 1, main () at bad.c:5
5         lalDebugLevel = 1;
(gdb) set var lalMemDbgUsrPtr = 0x1cf528
(gdb) watch lalIsMemDbgRetPtr
Hardware watchpoint 2: lalIsMemDbgRetPtr
(gdb) cont
Continuing.
Hardware watchpoint 2: lalIsMemDbgRetPtr

Old value = 0
New value = 1
0x0088d63c in LALMallocLong (n=5, file=0x1ff8 "bad.c", line=6) at LALMalloc.c:575
575       lalIsMemDbgPtr = lalIsMemDbgRetPtr = ( lalMemDbgRetPtr == lalMemDbgUsrPtr );
(gdb) up
#1  0x00001b84 in main () at bad.c:6
6         s = LALMalloc( 5 );
\end{verbatim}
So here is where the memory is allocated.  We want to find out where the
memory is being corrupted.
\begin{verbatim}
(gdb) set var lalMemDbgUsrPtr = 0x1cf530
(gdb) watch *lalMemDbgUsrPtr
Hardware watchpoint 3: *lalMemDbgUsrPtr
(gdb) cont
Continuing.
Hardware watchpoint 3: *lalMemDbgUsrPtr

Old value = -25
New value = 120 'x'
main () at bad.c:8
8         LALFree( s );
(gdb) list
3       {
4         char *s;
5         lalDebugLevel = 1;
6         s = LALMalloc( 5 );
7         s[8] = 'x';
8         LALFree( s );
9         LALCheckMemoryLeaks();
10        return 0;
11      }
\end{verbatim}
Notice that the program has stopped just \emph{after} the line in which
the array bounds were overwritten.


\vfill{\footnotesize\input{LALMallocCV}}

******************************************************* </lalLaTeX> */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#include <lal/LALConfig.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdio.h>
#include <lal/LALError.h>

/*
 *
 * XLAL Routines.
 *
 */

#define XLAL_TEST_POINTER( ptr, size, func )                               \
    if ( ! (ptr) && (size) )                                               \
       XLAL_ERROR_NULL( func, XLAL_ENOMEM );                               \
    else (void)(0)
#define XLAL_TEST_POINTER_LONG( ptr, size, func, file, line )              \
    if ( ! (ptr) && (size) )                                               \
    {                                                                      \
       char msg[64];                                                       \
       snprintf( msg, sizeof( msg ), "%s in %s:%d", func, file, line ); \
       XLAL_ERROR_NULL( msg, XLAL_ENOMEM );                                \
    }                                                                      \
    else (void)(0)


void * (XLALMalloc)( size_t n )
{
  static const char *func = "XLALMalloc";
  void *p;
  p = LALMallocShort( n );
  XLAL_TEST_POINTER( p, n, func );
  return p;
}

void * XLALMallocLong( size_t n, const char *file, int line )
{
  static const char *func = "XLALMallocLong";
  void *p;
  p = LALMallocLong( n, file, line );
  XLAL_TEST_POINTER_LONG( p, n, func, file, line );
  return p;
}

void * (XLALCalloc)( size_t m, size_t n )
{
  static const char *func = "XLALCalloc";
  void *p;
  p = LALCallocShort( m, n );
  XLAL_TEST_POINTER( p, m * n, func );
  return p;
}

void * XLALCallocLong( size_t m, size_t n, const char *file, int line )
{
  static const char *func = "XLALCallocLong";
  void *p;
  p = LALCallocLong( m, n, file, line );
  XLAL_TEST_POINTER_LONG( p, m * n, func, file, line );
  return p;
}

void * (XLALRealloc)( void *p, size_t n )
{
  static const char *func = "XLALRealloc";
  p = LALReallocShort( p, n );
  XLAL_TEST_POINTER( p, n, func );
  return p;
}

void * XLALReallocLong( void *p, size_t n, const char *file, int line )
{
  static const char *func = "XLALReallocLong";
  p = LALReallocLong( p, n, file, line );
  XLAL_TEST_POINTER_LONG( p, n, func, file, line );
  return p;
}

void XLALFree( void *p )
{
  if ( p ) LALFree( p );
  return;
}


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

NRCSID( LALMALLOCC, "$Id$" );

/* global variables to assist in memory debugging */
/* watch the value of these variables to find a particular alloc/free */
char  *lalMemDbgArgPtr  = NULL; /* set to ptr arg in free or realloc */
char  *lalMemDbgRetPtr  = NULL; /* set to ptr returned in alloc functions */
char  *lalMemDbgPtr     = NULL; /* set in both cases */
char  *lalMemDbgUsrPtr  = NULL; /* avaliable global memory pointer for user */
void **lalMemDbgUsrHndl = NULL; /* avaliable global memory handle for user */
int lalIsMemDbgArgPtr; /* ( lalMemDbgUsrPtr == lalMemDbgArgPtr ) */
int lalIsMemDbgRetPtr; /* ( lalMemDbgUsrPtr == lalMemDbgRetPtr ) */
int lalIsMemDbgPtr;    /* ( lalMemDbgUsrPtr == lalMemDbgPtr ) */


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

/* need this to turn off gcc warnings about unused functions */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Useful function for debugging */
/* Checks to make sure alloc list is OK */
/* Returns 0 if list is corrupted; 1 if list is OK */
UNUSED static int CheckAllocList( void )
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

/* Useful function for debugging */
/* Finds the node of the alloc list previous to the desired alloc */
/* Returns NULL if not found or if the alloc node is the head */
UNUSED static struct allocNode *FindPrevAlloc( void *p )
{
  struct allocNode *node = allocList;
  if ( p == node->addr ) /* top of list */
    return NULL;
  /* scroll through list to find node before the alloc */
  while ( node->next )
    if ( node->next->addr == p )
      return node;
    else
      node = node->next;
  return NULL;
}

/* Useful function for debugging */
/* Finds the node of the alloc list for the desired alloc */
/* Returns NULL if not found  */
UNUSED static struct allocNode *FindAlloc( void *p )
{
  struct allocNode *node = allocList;
  while ( node )
    if ( p == node->addr )
      return node;
    else
      node = node->next;
  return NULL;
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
    LALPrintError( "%s meminfo: allocating %ld bytes at address %p\n",
        func, n, p + nprefix );
  }

  if ( lalDebugLevel & LALWARNING && n == 0 )
  {
    LALPrintError( "%s warning: zero size allocation at address %p\n",
        func, p + nprefix );
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
    LALPrintError( "%s meminfo: freeing %ld bytes at address %p\n",
        func, n, p );
  }

  if ( lalDebugLevel & LALWARNING && n == 0 )
  {
    LALPrintError( "%s warning: tried to free a freed pointer at address %p\n",
        func, p );
  }

  if ( q[1] != magic )
  {
    lalRaiseHook( SIGSEGV, "%s error: wrong magic for pointer at address %p\n",
        func, p );
    return NULL;
  }

  if ( ( (long) n ) < 0 )
  {
    lalRaiseHook( SIGSEGV, "%s error: corrupt size descriptor for pointer at address %p\n",
        func, p );
    return NULL;
  }

  /* check for writing past end of array: */
  for ( i = n; i < padFactor * n; ++i )
  {
    if ( s[i + prefix] != (char) (i ^ padding) )
    {
      lalRaiseHook( SIGSEGV, "%s error: array bounds overwritten\n"
          "Byte %ld past end of array has changed\n"
          "Corrupted address: %p\nArray address: %p\n",
          func, i - n + 1, s + i + prefix, s + prefix );
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
  struct allocNode *newnode;
  if ( lalDebugLevel & LALNMEMTRK )
  {
    return p;
  }
  if ( ! p )
  {
    return NULL;
  }
  if ( ! ( newnode = malloc( sizeof( *newnode ) ) ) )
  {
    return NULL;
  }
  pthread_mutex_lock( &mut );
  newnode->addr = p;
  newnode->size = n;
  newnode->file = file;
  newnode->line = line;
  newnode->next = allocList;
  allocList = newnode;
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
    pthread_mutex_unlock( &mut );
    lalRaiseHook( SIGSEGV, "%s error: alloc %p not found\n", func, p );
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
      pthread_mutex_unlock( &mut );
      lalRaiseHook( SIGSEGV, "%s error: alloc %p not found\n", func, p );
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
    pthread_mutex_unlock( &mut );
    lalRaiseHook( SIGSEGV, "%s error: alloc %p not found\n", func, p );
    return NULL;
  }
  while ( p != node->addr )
  {
    if ( ! ( node = node->next ) )
    {
      pthread_mutex_unlock( &mut );
      lalRaiseHook( SIGSEGV, "%s error: alloc %p not found\n", func, p );
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
  lalMemDbgPtr = lalMemDbgRetPtr = q;
  lalIsMemDbgPtr = lalIsMemDbgRetPtr = ( lalMemDbgRetPtr == lalMemDbgUsrPtr );
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
  lalMemDbgPtr = lalMemDbgRetPtr = q;
  lalIsMemDbgPtr = lalIsMemDbgRetPtr = ( lalMemDbgRetPtr == lalMemDbgUsrPtr );
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

  lalMemDbgPtr = lalMemDbgArgPtr = q;
  lalIsMemDbgPtr = lalIsMemDbgArgPtr = ( lalMemDbgArgPtr == lalMemDbgUsrPtr );
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
  lalMemDbgPtr = lalMemDbgRetPtr = q;
  lalIsMemDbgPtr = lalIsMemDbgRetPtr = ( lalMemDbgRetPtr == lalMemDbgUsrPtr );

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
  lalMemDbgPtr = lalMemDbgArgPtr = q;
  lalIsMemDbgPtr = lalIsMemDbgArgPtr = ( lalMemDbgArgPtr == lalMemDbgUsrPtr );
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
