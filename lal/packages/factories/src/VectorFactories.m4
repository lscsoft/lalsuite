/*----------------------------------------------------------------------- 
 * 
 * File Name: VectorFactories.c
 * 
 * Author: Finn, L. S.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * SYNOPSIS 
 * void TYPECODECreateVector( Status *, TYPEVector **vector, UINT4 length );
 * void TYPECODEDestroyVector( Status *, TYPEVector **vector);
 * 
 * DESCRIPTION 
 * Create/destroy a TYPEVector object. 
 * 
 * DIAGNOSTICS 
 * Illegal length, vector == NULL, *vector != NULL, malloc failure
 *
 * CALLS
 * LALMalloc
 * 
 *-----------------------------------------------------------------------
 */

#include "LALStdlib.h"
#include "AVFactories.h"

NRCSID( VECTORFACTORIESC, "$Id$" );

define(`TYPECODE',`Z')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`C')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`D')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`S')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I2')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I4')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I8')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U2')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U4')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U8')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`CHAR')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`')
include(`CreateVector.m4')
include(`DestroyVector.m4')
