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

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
#endif

NRCSID (LALSTDLIBH, "$Id$");

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _LALMALLOC_H
#include "LALMalloc.h"
#ifndef _LALMALLOC_H
#define _LALMALLOC_H
#endif
#endif

#define LALFopen  fopen    /* will be replaced by custom unit */
#define LALFclose fclose   /* will be replaced by custom unit */

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

#ifndef _LALSTATUSMACROS_H
#include "LALStatusMacros.h"
#ifndef _LALSTATUSMACROS_H
#define _LALSTATUSMACROS_H
#endif
#endif

#endif
