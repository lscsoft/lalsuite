/************************************ <lalVerbatim file="LALMallocHV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALMalloc.h}}
\label{s:LALMalloc.h}

Provides standard LAL memory allocation/deallocation routines.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALMalloc.h>
\end{verbatim}

\noindent This header covers routines that replace the standard
\verb+malloc()+, \verb+calloc()+, \verb+realloc()+, and \verb+free()+.
All memory allocation and deallocation in LAL should use these
replacement functions.  If the \verb+NDEBUG+ flag is set at compile
time, the LAL routines are \verb+#define+d to be the same as the
standard C routines.

\vfill{\footnotesize\input{LALMallocHV}}
\newpage\input{LALMallocC}
\newpage\input{LALMallocTestC}

</lalLaTeX> */



#ifndef _LALMALLOC_H
#define _LALMALLOC_H

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALMALLOCH, "$Id$" );

void *XLALMalloc( size_t n );
void *XLALMallocLong( size_t n, const char *file, int line );
void *XLALCalloc( size_t m, size_t n );
void *XLALCallocLong( size_t m, size_t n, const char *file, int line );
void *XLALRealloc( void *p, size_t n );
void *XLALReallocLong( void *p, size_t n, const char *file, int line );
void  XLALFree( void *p );
#define XLALMalloc( n )        XLALMallocLong( n, __FILE__, __LINE__ )
#define XLALCalloc( m, n )     XLALCallocLong( m, n, __FILE__, __LINE__ )
#define XLALRealloc( p, n )    XLALReallocLong( p, n, __FILE__, __LINE__ )

#if defined NDEBUG || defined LAL_NDEBUG

#define LALMalloc                          malloc     
#define LALMallocShort                     malloc     
#define LALMallocLong( n, file, line )     malloc( n )
#define LALCalloc                          calloc        
#define LALCallocShort                     calloc        
#define LALCallocLong( m, n, file, line )  calloc( m, n )
#define LALRealloc                         realloc
#define LALReallocShort                    realloc
#define LALReallocLong( p, n, file, line ) realloc( p, n )
#define LALFree                            free      
#define LALCheckMemoryLeaks()

#else

#define LALMalloc( n )        LALMallocLong( n, __FILE__, __LINE__ )
#define LALCalloc( m, n )     LALCallocLong( m, n, __FILE__, __LINE__ )
#define LALRealloc( p, n )    LALReallocLong( p, n, __FILE__, __LINE__ )

/* global variables to assist in memory debugging */
/* watch the value of these variables to find a particular alloc/free */
extern char  *lalMemDbgArgPtr;   /* set to ptr arg in free or realloc */
extern char  *lalMemDbgRetPtr;   /* set to ptr returned in alloc functions */
extern char  *lalMemDbgPtr;      /* set in both cases */
extern char  *lalMemDbgUsrPtr;   /* avaliable global memory pointer for user */
extern void **lalMemDbgUsrHndl;  /* avaliable global memory handle for user */
extern int    lalIsMemDbgArgPtr; /* ( lalMemDbgUsrPtr == lalMemDbgArgPtr ) */
extern int    lalIsMemDbgRetPtr; /* ( lalMemDbgUsrPtr == lalMemDbgRetPtr ) */
extern int    lalIsMemDbgPtr;    /* ( lalMemDbgUsrPtr == lalMemDbgPtr ) */


void *
LALMallocShort( size_t n );

void *
LALMallocLong( size_t n, const char *file, int line );

void
LALFree( void *p );

void *
LALCallocShort( size_t m, size_t n );

void *
LALCallocLong( size_t m, size_t n, const char *file, int line );

void *
LALReallocShort( void *p, size_t n );

void *
LALReallocLong( void *p, size_t n, const char *file, int line );

void
LALCheckMemoryLeaks( void );

#endif /* NDEBUG || LAL_NDEBUG */

#ifdef  __cplusplus
}
#endif

#endif /* _LALMALLOC_H */
