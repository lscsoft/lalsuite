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

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
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

#endif
