/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

/*************** <lalVerbatim file="LALMallocTestCV"> *************
$Id$
**************** </lalVerbatim> ***********************************/

/* <lalLaTeX>

\subsection{Program \texttt{LALMallocTest.c}}
\label{s:LALMallocTest.c}

Tests the routines in \verb@LALMalloc.h@.

\subsubsection*{Usage}
\begin{verbatim}
LALMallocTest
\end{verbatim}

\subsubsection*{Description}

This program has ugly code for testing the LAL memory allocation and freeing
routines.

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success.                      \\
\tt 1 & Failure.                      \\
\hline
\end{tabular}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()
LALCalloc()
LALRealloc()
LALFree()
LALCheckMemoryLeaks()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALMallocTestCV}}

</lalLaTeX> */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <setjmp.h>
#include <signal.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>

NRCSID (LALMALLOCTESTC,"$Id$");

char caughtMessage[1024];
jmp_buf jump;
FILE *mystderr;

/* replacement for LALRaise */
static int TestRaise( int sig, const char *fmt, ... )
{
  va_list ap;
  va_start( ap, fmt );
  vsnprintf( caughtMessage, sizeof( caughtMessage ), fmt, ap );
  va_end( ap );
  longjmp( jump, sig );
  return -1;
}

#define STR( a ) #a
#define XSTR( a ) STR( a )
#define LINE ":" XSTR( __LINE__ ) ")\n"
#define trial( func, sig, msg ) \
do { \
  int val; \
  if ( ! ( val = setjmp( jump ) ) ) \
  { \
    func; \
    if ( sig ) \
    { \
      fprintf( mystderr, "Error: no signal raised! (" #func LINE ); \
      return 1; \
    } \
  } \
  else \
  { \
    if ( val != sig ) \
    { \
      fprintf( mystderr, "Error: wrong signal raised! (" #func LINE ); \
      fprintf( mystderr, "Received: %d %s", val, caughtMessage ); \
      fprintf( mystderr, "Expected: %d %s\n", sig, msg ); \
      return 1; \
    } \
    if ( ! strstr( caughtMessage, msg ) ) \
    { \
      fprintf( mystderr, "Error: wrong message! (" #func LINE ); \
      fprintf( mystderr, "Received: %d %s", val, caughtMessage ); \
      fprintf( mystderr, "Expected: %d %s\n", sig, msg ); \
      return 1; \
    } \
  } \
} \
while ( 0 )

#define die( msg ) ( fputs( "Error: " #msg "\n", mystderr ), exit( 1 ), 1 )


int lalDebugLevel = LALMEMDBG;

/* make these global so they don't get clobbered by longjmp */
size_t   i;
size_t   j;
size_t   n;
size_t  *p;
size_t  *q;
size_t  *r;
size_t  *s;
size_t **v;

/* do a bunch of allocations/deallocations that are OK */
static int testOK( void )
{
  int keep = lalDebugLevel;

  lalDebugLevel &= ~( LALNMEMDBG | LALNMEMPAD | LALNMEMTRK );
  trial( p = LALMalloc( 1024 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) p[i] = i;
  trial( q = LALCalloc( 1024, sizeof( *q ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( q[i] ) die( memory not blanked );
  trial( p = LALRealloc( p, 4096 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( p[i] != i ) die( memory not copied );
  trial( q = LALRealloc( q, 0 ), 0, "" );
  if ( q ) die( memory not freed );
  trial( LALCheckMemoryLeaks(), SIGSEGV, "LALCheckMemoryLeaks: memory leak\n" );
  if ( *( p - 1 ) != (size_t)0xABadCafe ) die( wrong magic );
  if ( *( p - 2 ) != 4096 * sizeof( *p ) ) die( wrong size );
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  lalDebugLevel |= LALNMEMPAD;
  trial( p = LALMalloc( 1024 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) p[i] = i;
  trial( q = LALCalloc( 1024, sizeof( *q ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( q[i] ) die( memory not blanked );
  trial( p = LALRealloc( p, 4096 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( p[i] != i ) die( memory not copied );
  trial( q = LALRealloc( q, 0 ), 0, "" );
  if ( q ) die( memory not freed );
  trial( LALCheckMemoryLeaks(), SIGSEGV, "LALCheckMemoryLeaks: memory leak\n" );
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  lalDebugLevel &= ~LALNMEMPAD;
  lalDebugLevel |= LALNMEMTRK;
  trial( p = LALMalloc( 1024 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) p[i] = i;
  trial( q = LALCalloc( 1024, sizeof( *q ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( q[i] ) die( memory not blanked );
  trial( p = LALRealloc( p, 4096 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( p[i] != i ) die( memory not copied );
  trial( q = LALRealloc( q, 0 ), 0, "" );
  if ( q ) die( memory not freed );
  trial( LALCheckMemoryLeaks(), SIGSEGV, "LALCheckMemoryLeaks: memory leak\n" );
  if ( *( p - 1 ) != (size_t)0xABadCafe ) die( wrong magic );
  if ( *( p - 2 ) != 4096 * sizeof( *p ) ) die( wrong size );
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  lalDebugLevel |= LALNMEMDBG;
  trial( p = LALMalloc( 1024 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) p[i] = i;
  trial( q = LALCalloc( 1024, sizeof( *q ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( q[i] ) die( memory not blanked );
  trial( p = LALRealloc( p, 4096 * sizeof( *p ) ), 0, "" );
  for ( i = 0; i < 1024; ++i ) if ( p[i] != i ) die( memory not copied );
  trial( q = LALRealloc( q, 0 ), 0, "" );
  /* if ( q ) die( memory not freed ); */
  trial( LALCheckMemoryLeaks(), 0, "" );
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  lalDebugLevel = keep;
  return 0;
}


/* test to make sure padding does what it's supposed to do */
static int testPadding( void )
{
  int keep = lalDebugLevel;

  lalDebugLevel |= LALNMEMTRK;
  lalDebugLevel &= ~( LALNMEMDBG | LALNMEMPAD );

  /* try to free NULL pointer */
  trial( LALFree( NULL ), SIGSEGV, "error: tried to free NULL pointer" );

  /* wrong magic */
  trial( p = LALMalloc( 2 * sizeof( *p ) ), 0, "" );
  p[-1] = 4;
  trial( LALFree( p ), SIGSEGV, "error: wrong magic" );
  p[-1] = 0xABadCafe;
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  /* corrupt size */
  trial( p = LALMalloc( 4 * sizeof( *p ) ), 0, "");
  n = p[-2];
  p[-2] = -1;
  trial( LALFree( p ), SIGSEGV, "error: corrupt size descriptor" );
  p[-2] = n;
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  /* overwritten array bounds */
  trial( p = LALMalloc( 8 * sizeof( *p ) ), 0, "" );
  n = p[8];
  p[8] = 0;
  trial( LALFree( p ), SIGSEGV, "error: array bounds overwritten" );
  p[8] = n;
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  /* free too much memory */
  q = malloc( 4 * sizeof( *p ) );
  trial( p = LALMalloc( sizeof( *p ) ), 0, "" );
  memcpy( q, p - 2, 4 * sizeof( *p ) );
  trial( LALFree( p ), 0, "" );
  trial( LALFree( q + 2 ), SIGSEGV, "error: lalMallocTotal too small" );
  free( q );
  trial( LALCheckMemoryLeaks(), 0, "" );

  lalDebugLevel = keep;
  return 0;
}

/* test to make sure alloc list does what it's supposed to do */
static int testAllocList( void )
{
  int keep = lalDebugLevel;

  s = malloc( sizeof( *s ) );

  lalDebugLevel |= LALNMEMPAD;
  lalDebugLevel &= ~( LALNMEMDBG | LALNMEMTRK );

  /* empty allocation list */
  trial( LALCheckMemoryLeaks(), 0, "" );
  trial( LALFree( s ), SIGSEGV, "not found" );

  /* can't find allocation in PopAlloc */
  trial( p = LALMalloc( 2 * sizeof( *p ) ), 0, "" );
  trial( q = LALMalloc( 4 * sizeof( *q ) ), 0, "" );
  trial( r = LALMalloc( 8 * sizeof( *r ) ), 0, "" );
  trial( LALFree( s ), SIGSEGV, "not found" );
  trial( LALFree( p ), 0, "" );
  trial( LALFree( r ), 0, "" );
  trial( LALCheckMemoryLeaks(), SIGSEGV, "memory leak" );
  trial( LALFree( q ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  /* can't fine allocation in ModAlloc */
  trial( s = LALRealloc( s, 1024 ), SIGSEGV, "not found" );
  trial( p = LALRealloc( NULL, 2 * sizeof( *p ) ), 0, "" );
  /* trial( s = LALRealloc( s, 1024 ), SIGSEGV, "not found" ); */
  trial( LALFree( p ), 0, "" );
  trial( LALCheckMemoryLeaks(), 0, "" );

  free( s );
  lalDebugLevel = keep;
  return 0;
}

/* stress test the realloc routine */
static int stressTestRealloc( void )
{
  const size_t nmax = 256;
  int keep = lalDebugLevel;

  v = NULL;

  lalDebugLevel &= ~( LALMEMINFO | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK );

  /* ascending */
  for ( n = 1; n <= nmax; ++n )
  {
    size_t *u;
    trial( v = LALRealloc( v, n * sizeof( *v ) ), 0, "" );
    trial( u = v[n - 1] = LALRealloc( NULL, n * sizeof( **v ) ), 0, "" );
    for ( i = 0; i < n; ++i ) u[i] = n - 1;
    for ( i = 0; i < n; ++i )
    {
      trial( u = v[i] = LALRealloc( v[i], n * sizeof( *u ) ), 0, "" );
      for ( j = 0; j < n - 1; ++j )
        if ( u[j] != n - 1 ) die( wrong contents );
      for ( j = 0; j < n; ++j )
        u[j] = n;
    }
  }

  for ( n = 0; n < nmax; ++n )
  {
    trial( v[n] = LALRealloc( v[n], 0 ), 0, "" );
  }
  trial( v = LALRealloc( v, 0 ), 0, "" );

  trial( LALCheckMemoryLeaks(), 0, "" );
  lalDebugLevel = keep;
  return 0;
}


int main( void )
{
#if defined(NDEBUG) || defined(LAL_NDEBUG) /* debugging is turned off */
  return 77; /* don't do any testing */
#else
  /* get rid of annoying messages from elsewhere */
  setvbuf( mystderr = stdout, NULL, _IONBF, 0 );
  freopen( "/dev/null", "w", stderr );

  lalRaiseHook = TestRaise;

  if ( lalNoDebug ) /* library was not compiled with debugging */
    return 77; /* don't do any testing */

  if ( testOK() ) return 1;
  if ( testPadding() ) return 1;
  if ( testAllocList() ) return 1;
  if ( stressTestRealloc() ) return 1;

  trial( LALCheckMemoryLeaks(), 0, "" );

  return 0;
#endif
}
