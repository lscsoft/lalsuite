/*----------------------------------------------------------------------- 
 * 
 * File Name: VectorSequenceFactories.c
 * 
 * Author: Finn, L. S.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * SYNOPSIS 
 * void TYPECODECreateVectorSequence(
 *     Status                  *status,
 *     TYPEVectorSequence     **vseq,
 *     CreateVectorSequenceIn  *in
 *     );
 * void TYPECODEDestroyVectorSequence(
 *     Status              *status,
 *     TYPEVectorSequence **vector
 *     );
 * 
 * DESCRIPTION 
 * Create/destroy a TYPEVectorSequence object. 
 * 
 * DIAGNOSTICS 
 * Illegal sequence length, illegal vectorLength, vseq == NULL, 
 * *vseq != NULL, malloc failure 
 *
 * CALLS
 * LALMalloc
 * LALFree
 * 
 *-----------------------------------------------------------------------
 */

#include "LALStdlib.h"
#include "SeqFactories.h"

NRCSID( VECTORSEQUENCEFACTORIESC, "$Id$" );

define(`TYPECODE',`Z')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`C')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`D')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`S')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I2')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I4')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I8')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U2')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U4')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U8')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`CHAR')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')
