/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
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

#ifndef _LALMALLOC_H
#define _LALMALLOC_H

#include <stddef.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/** \addtogroup LALMalloc_h *//*@{ */
void *XLALMalloc(size_t n);
void *XLALMallocLong(size_t n, const char *file, int line);
void *XLALCalloc(size_t m, size_t n);
void *XLALCallocLong(size_t m, size_t n, const char *file, int line);
void *XLALRealloc(void *p, size_t n);
void *XLALReallocLong(void *p, size_t n, const char *file, int line);
void XLALFree(void *p);
#ifndef SWIG    /* exclude from SWIG interface */
#define XLALMalloc( n )        XLALMallocLong( n, __FILE__, __LINE__ )
#define XLALCalloc( m, n )     XLALCallocLong( m, n, __FILE__, __LINE__ )
#define XLALRealloc( p, n )    XLALReallocLong( p, n, __FILE__, __LINE__ )
#endif /* SWIG */
/*@}*/


#if defined NDEBUG || defined LAL_NDEBUG

#ifndef SWIG    /* exclude from SWIG interface */
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
#endif /* SWIG */

#else

#ifndef SWIG    /* exclude from SWIG interface */
#define LALMalloc( n )        LALMallocLong( n, __FILE__, __LINE__ )
#define LALCalloc( m, n )     LALCallocLong( m, n, __FILE__, __LINE__ )
#define LALRealloc( p, n )    LALReallocLong( p, n, __FILE__, __LINE__ )
#endif /* SWIG */

/* global variables to assist in memory debugging */
/* watch the value of these variables to find a particular alloc/free */
#ifndef SWIG    /* exclude from SWIG interface */
extern char *lalMemDbgArgPtr;   /* set to ptr arg in free or realloc */
extern char *lalMemDbgRetPtr;   /* set to ptr returned in alloc functions */
extern char *lalMemDbgPtr;      /* set in both cases */
extern char *lalMemDbgUsrPtr;   /* avaliable global memory pointer for user */
extern void **lalMemDbgUsrHndl; /* avaliable global memory handle for user */
extern int lalIsMemDbgArgPtr;   /* ( lalMemDbgUsrPtr == lalMemDbgArgPtr ) */
extern int lalIsMemDbgRetPtr;   /* ( lalMemDbgUsrPtr == lalMemDbgRetPtr ) */
extern int lalIsMemDbgPtr;      /* ( lalMemDbgUsrPtr == lalMemDbgPtr ) */
#endif /* SWIG */


/** \addtogroup LALMalloc_h *//*@{ */
void *LALMallocShort(size_t n);
void *LALMallocLong(size_t n, const char *file, int line);
void *LALCallocShort(size_t m, size_t n);
void LALFree(void *p);
void *LALCallocLong(size_t m, size_t n, const char *file, int line);
void *LALReallocShort(void *p, size_t n);
void *LALReallocLong(void *p, size_t n, const char *file, int line);
/*@}*/

#endif /* NDEBUG || LAL_NDEBUG */

void (LALCheckMemoryLeaks) (void);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALMALLOC_H */
