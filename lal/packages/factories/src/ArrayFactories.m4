/*----------------------------------------------------------------------- 
 * 
 * File Name: ArrayFactories.c
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * SYNOPSIS 
 * void TYPECODECreateArray( Status *, TYPEArray **array, UINT4Vector *dimlength );
 * void TYPECODEDestroyArray( Status *, TYPEArray **array);
 * 
 * DESCRIPTION 
 * Create/destroy a TYPEArray object. 
 * 
 * DIAGNOSTICS 
 * Illegal length, array == NULL, *array != NULL, malloc failure
 *
 * CALLS
 * LALMalloc
 * LALFree
 * 
 *-----------------------------------------------------------------------
 */

#include <string.h>
#include "LALStdlib.h"
#include "AVFactories.h"

NRCSID( ARRAYFACTORIESC, "$Id$" );

define(`TYPECODE',`Z')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`C')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`D')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`S')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I2')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I4')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I8')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U2')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U4')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U8')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`')
include(`CreateArray.m4')
include(`DestroyArray.m4')
