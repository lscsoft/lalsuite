/*-----------------------------------------------------------------------
 *
 * File Name: VectorOps.h
 *
 * Authors: Creighton, J. D. E., Creighton, T. D.
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

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
#endif

NRCSID (VECTOROPSH, "$Id$");

#define VECTOROPS_ENULL 1
#define VECTOROPS_ESIZE 2
#define VECTOROPS_ESZMM 4

#define VECTOROPS_MSGENULL "Null pointer"
#define VECTOROPS_MSGESIZE "Invalid input size"
#define VECTOROPS_MSGESZMM "Size mismatch"

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

#endif
