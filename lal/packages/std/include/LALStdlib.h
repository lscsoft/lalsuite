/*-----------------------------------------------------------------------
 *
 * File Name: LALstdlib.h
 *
 * Author: Finn, L. S.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * LALStdlib.h
 *
 * SYNOPSIS
 * #include "LALStdlib.h"
 *
 * DESCRIPTION
 * Defines LIGO/LSC Analysis Library support functions for file I/O, 
 * memory allocation, etc. Provides prototypes for support functions 
 * used to manipulate status structures.
 *
 * DIAGNOSTICS
 *
 *----------------------------------------------------------------------
 */

#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H

#include <stdio.h>
#include "LALRCSID.h"
#include "LALMalloc.h"
#include "LALDatatypes.h"
#include "LALStatusMacros.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALSTDLIBH, "$Id$");

#define LALFopen  fopen    /* will be replaced by custom unit */
#define LALFclose fclose   /* will be replaced by custom unit */

#ifdef  __cplusplus
}
#endif

#endif /* _LALSTDLIB_H */
