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
CCVectorMultiply (
    Status               *,
    COMPLEX8Vector       *, 
    const COMPLEX8Vector *,
    const COMPLEX8Vector *
    );

void
CCVectorMultiplyConjugate (
    Status               *,
    COMPLEX8Vector       *, 
    const COMPLEX8Vector *,
    const COMPLEX8Vector *
    );

void
CCVectorDivide (
    Status               *,
    COMPLEX8Vector       *, 
    const COMPLEX8Vector *,
    const COMPLEX8Vector *
    );

void
CVectorAbs (
    Status               *,
    REAL4Vector          *,
    const COMPLEX8Vector *
    );

void
CVectorAngle (
    Status               *,
    REAL4Vector          *,
    const COMPLEX8Vector *
    );

void
UnwrapREAL4Angle (
    Status               *,
    REAL4Vector          *,
    const REAL4Vector    *
    );

void
ZZVectorMultiply (
    Status                *,
    COMPLEX16Vector       *, 
    const COMPLEX16Vector *,
    const COMPLEX16Vector *
    );

void
ZZVectorMultiplyConjugate (
    Status                *,
    COMPLEX16Vector       *, 
    const COMPLEX16Vector *,
    const COMPLEX16Vector *
    );

void
ZZVectorDivide (
    Status                *,
    COMPLEX16Vector       *, 
    const COMPLEX16Vector *,
    const COMPLEX16Vector *
    );

void
ZVectorAbs (
    Status                *,
    REAL8Vector           *,
    const COMPLEX16Vector *
    );

void
ZVectorAngle (
    Status                *,
    REAL8Vector           *,
    const COMPLEX16Vector *
    );

void
UnwrapREAL8Angle (
    Status               *,
    REAL8Vector          *,
    const REAL8Vector    *
    );

void
SCVectorMultiply(
    Status               *,
    COMPLEX8Vector       *,
    const REAL4Vector    *,
    const COMPLEX8Vector *
    );

void
SSVectorMultiply(
    Status               *,
    REAL4Vector          *,
    const REAL4Vector    *,
    const REAL4Vector    *
    );

void
DZVectorMultiply(
    Status                *,
    COMPLEX16Vector       *,
    const REAL8Vector     *,
    const COMPLEX16Vector *
    );

void
DDVectorMultiply(
    Status               *,
    REAL8Vector          *,
    const REAL8Vector    *,
    const REAL8Vector    *
    );


#ifdef  __cplusplus
}
#endif

#endif /* _VECTOROPS_H */
