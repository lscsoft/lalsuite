/*----------------------------------------------------------------------- 
 * 
 * File Name: LALMalloc.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALMALLOC_H
#define _LALMALLOC_H

#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (LALMALLOCH, "$Id$");

void *
LALMalloc (size_t n);

void
LALFree (void *p);

void *
LALCalloc (size_t m, size_t n);

void *
LALRealloc (void *p, size_t n);

void
LALCheckMemoryLeaks (void);


#ifdef  __cplusplus
}
#endif

#endif /* _LALMALLOC_H */
