/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**** <lalVerbatim file="VectorOpsHV">
 * Author: J. D. E. Creighton, T. D. Creighton, A. M. Sintes
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{VectorOps.h}}
 *
 * Basic vector manipulation operations.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/VectorOps.h>
 * \end{verbatim}
 *
 * \subsection*{Error conditions}
 * \input{VectorOpsHE}
 *
 * \vfill{\footnotesize\input{VectorOpsHV}}
 * \newpage\input{VectorMultiplyC}
 * \newpage\input{VectorPolarC}
 * \newpage\input{VectorOpsTestC}
 *
 **** </lalLaTeX> */

#ifndef _VECTOROPS_H
#define _VECTOROPS_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (VECTOROPSH, "$Id$");

/**** <lalErrTable file="VectorOpsHE"> */
#define VECTOROPSH_ENULL 1
#define VECTOROPSH_ESIZE 2
#define VECTOROPSH_ESZMM 4
#define VECTOROPSH_ESAME 8

#define VECTOROPSH_MSGENULL "Null pointer"
#define VECTOROPSH_MSGESIZE "Invalid input size"
#define VECTOROPSH_MSGESZMM "Size mismatch"
#define VECTOROPSH_MSGESAME "Input/Output data vectors are the same"
/**** </lalErrTable> */

/*
 *
 * XLAL Routines.
 *
 */

/* sigle precision */
COMPLEX8Vector * XLALCCVectorDivide( COMPLEX8Vector *out, const COMPLEX8Vector *in1, const COMPLEX8Vector *in2 );
COMPLEX8Vector * XLALCCVectorMultiply( COMPLEX8Vector *out, const COMPLEX8Vector *in1, const COMPLEX8Vector *in2 );
COMPLEX8Vector * XLALCCVectorMultiplyConjugate( COMPLEX8Vector *out, const COMPLEX8Vector *in1, const COMPLEX8Vector *in2 );
COMPLEX8Vector * XLALSCVectorMultiply( COMPLEX8Vector *out, const REAL4Vector *in1, const COMPLEX8Vector *in2 );
REAL4Vector * XLALSSVectorMultiply( REAL4Vector *out, const REAL4Vector *in1, const REAL4Vector *in2 );

int XLALCOMPLEX8VectorAbs( REAL4Vector *out, const COMPLEX8Vector *in );
int XLALCOMPLEX8VectorArg( REAL4Vector *out, const COMPLEX8Vector *in );
int XLALREAL4VectorUnwrapAngle( REAL4Vector *out, const REAL4Vector *in );


/* double precision */
COMPLEX16Vector * XLALZZVectorDivide( COMPLEX16Vector *out, const COMPLEX16Vector *in1, const COMPLEX16Vector *in2 );
COMPLEX16Vector * XLALZZVectorMultiply( COMPLEX16Vector *out, const COMPLEX16Vector *in1, const COMPLEX16Vector *in2 );
COMPLEX16Vector * XLALZZVectorMultiplyConjugate( COMPLEX16Vector *out, const COMPLEX16Vector *in1, const COMPLEX16Vector *in2 );
COMPLEX16Vector * XLALDZVectorMultiply( COMPLEX16Vector *out, const REAL8Vector *in1, const COMPLEX16Vector *in2 );
REAL8Vector * XLALDDVectorMultiply( REAL8Vector *out, const REAL8Vector *in1, const REAL8Vector *in2 );

int XLALCOMPLEX16VectorAbs( REAL8Vector *out, const COMPLEX16Vector *in );
int XLALCOMPLEX16VectorArg( REAL8Vector *out, const COMPLEX16Vector *in );
int XLALREAL8VectorUnwrapAngle( REAL8Vector *out, const REAL8Vector *in );

/*
 *
 * LAL Routines.
 *
 */

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
#pragma {
}
#endif

#endif /* _VECTOROPS_H */
