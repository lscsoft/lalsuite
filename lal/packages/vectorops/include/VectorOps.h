/*-----------------------------------------------------------------------
 *
 * File Name: VectorOps.h
 *
 * Authors: Creighton, J. D. E., Creighton, T. D., Sintes, A. M.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * VectorOps.h
 *
 * SYNOPSIS
 * #include "VectorOps.h"
 *
 * DESCRIPTION
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _VECTOROPS_H
#define _VECTOROPS_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (VECTOROPSH, "$Id$");

#define VECTOROPS_ENULL 1
#define VECTOROPS_ESIZE 2
#define VECTOROPS_ESZMM 4
#define VECTOROPS_ESAME 8

#define VECTOROPS_MSGENULL "Null pointer"
#define VECTOROPS_MSGESIZE "Invalid input size"
#define VECTOROPS_MSGESZMM "Size mismatch"
#define VECTOROPS_MSGESAME "Input/Output data vectors are the same"

void
LALCCVectorMultiply (
    LALStatus               *,
    COMPLEX8Vector       *, 
    const COMPLEX8Vector *,
    const COMPLEX8Vector *
    );

void
LALCCVectorMultiplyConjugate (
    LALStatus               *,
    COMPLEX8Vector       *, 
    const COMPLEX8Vector *,
    const COMPLEX8Vector *
    );

void
LALCCVectorDivide (
    LALStatus               *,
    COMPLEX8Vector       *, 
    const COMPLEX8Vector *,
    const COMPLEX8Vector *
    );

void
LALCVectorAbs (
    LALStatus               *,
    REAL4Vector          *,
    const COMPLEX8Vector *
    );

void
LALCVectorAngle (
    LALStatus               *,
    REAL4Vector          *,
    const COMPLEX8Vector *
    );

void
LALUnwrapREAL4Angle (
    LALStatus               *,
    REAL4Vector          *,
    const REAL4Vector    *
    );

void
LALZZVectorMultiply (
    LALStatus                *,
    COMPLEX16Vector       *, 
    const COMPLEX16Vector *,
    const COMPLEX16Vector *
    );

void
LALZZVectorMultiplyConjugate (
    LALStatus                *,
    COMPLEX16Vector       *, 
    const COMPLEX16Vector *,
    const COMPLEX16Vector *
    );

void
LALZZVectorDivide (
    LALStatus                *,
    COMPLEX16Vector       *, 
    const COMPLEX16Vector *,
    const COMPLEX16Vector *
    );

void
LALZVectorAbs (
    LALStatus                *,
    REAL8Vector           *,
    const COMPLEX16Vector *
    );

void
LALZVectorAngle (
    LALStatus                *,
    REAL8Vector           *,
    const COMPLEX16Vector *
    );

void
LALUnwrapREAL8Angle (
    LALStatus               *,
    REAL8Vector          *,
    const REAL8Vector    *
    );

void
LALSCVectorMultiply(
    LALStatus               *,
    COMPLEX8Vector       *,
    const REAL4Vector    *,
    const COMPLEX8Vector *
    );

void
LALSSVectorMultiply(
    LALStatus               *,
    REAL4Vector          *,
    const REAL4Vector    *,
    const REAL4Vector    *
    );

void
LALDZVectorMultiply(
    LALStatus                *,
    COMPLEX16Vector       *,
    const REAL8Vector     *,
    const COMPLEX16Vector *
    );

void
LALDDVectorMultiply(
    LALStatus               *,
    REAL8Vector          *,
    const REAL8Vector    *,
    const REAL8Vector    *
    );


#ifdef  __cplusplus
}
#endif

#endif /* _VECTOROPS_H */
