/******** <lalVerbatim file="MatrixHV"> ********
Author: Matthew M. Tibbits
$Id$
********* </lalVerbatim> ********/

/* <lalLaTeX>
 * 
 * \section{Header \texttt{Matrix.h}}
 *
 * Matlab Routines to handle Matrices \& Vectors.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Matrix.h>
 * \end{verbatim}
 * 
 * \subsection*{Error conditions}
 * \input{MatrixHE}
 * 
 * \vfill{\footnotesize\input{MatrixHV}}
 * \newpage\input{MatrixMultiplyC}
 * \newpage\input{MatrixDivideC}
 * \newpage\input{MatrixPowerC}
 * \newpage\input{MiscMatlabC}
 * 
 * </lalLaTeX> */

#include <math.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/AVFactories.h>

#ifndef __MATLAB_MATRIX_H__

#define __MATLAB_MATRIX_H__

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( MATLABMATRIXH, "$Id$");

/**************************** <lalErrTable file="MatrixHE"> */
#define	MATLABMATRIXH_EARG 1
#define	MATLABMATRIXH_ECHK 2
#define	MATLABMATRIXH_EFLS 3
#define	MATLABMATRIXH_EUSE 4
#define	MATLABMATRIXH_ENULL 5
#define	MATLABMATRIXH_EALOC 6
#define	MATLABMATRIXH_EFPMS 7
#define MATLABMATRIXH_ENUMZ 8
#define	MATLABMATRIXH_ELNTH 9
#define MATLABMATRIXH_ENNUL 10
#define	MATLABMATRIXH_MSGEARG "Error parsing command-line arguments"
#define	MATLABMATRIXH_MSGECHK "Error checking failed to catch bad data"
#define MATLABMATRIXH_MSGEFLS "Incorrect answer for valid data"
#define	MATLABMATRIXH_MSGEUSE "Bad user-entered data"
#define	MATLABMATRIXH_MSGENULL "Null Pointer."
#define MATLABMATRIXH_MSGEALOC "Memory Allocation Error"
#define MATLABMATRIXH_MSGEFPMS "Filter Parameter Structure Error"
#define MATLABMATRIXH_MSGENUMZ "Incorrect number of command line arguments"
#define MATLABMATRIXH_MSGELNTH "Vector/Array of Improper Length"
#define MATLABMATRIXH_MSGENNUL "Non-Null Pointer that should be NULL"
/***************************** </lalErrTable> */

/* Multiply */

void LALDDotStarDVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			B,
        REAL8Vector		*A
);

void LALDVectorDotStarDVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        REAL8Vector		*A
);

void LALDDotStarDArray (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        REAL8Array		*B
);

void LALDArrayDotStarDArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL8Array		*B
);

void LALDDotStarSVector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL8			B,
        REAL4Vector		*A
);

void LALDVectorDotStarSVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        REAL4Vector		*A
);

void LALDDotStarSArray (
        LALStatus               *status,
        REAL4Array		**result,
        REAL8			A,
        REAL4Array		*B
);

void LALDArrayDotStarSArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL4Array		*B
);

void LALDDotStarI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        REAL8			B,
        INT2Vector		*A
);

void LALDVectorDotStarI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        INT2Vector		*A
);

void LALDDotStarI2Array (
        LALStatus               *status,
        INT2Array		**result,
        REAL8			A,
        INT2Array		*B
);

void LALDArrayDotStarI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT2Array		*B
);

void LALDDotStarI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        REAL8			B,
        INT4Vector		*A
);

void LALDVectorDotStarI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        INT4Vector		*A
);

void LALDDotStarI4Array (
        LALStatus               *status,
        INT4Array		**result,
        REAL8			A,
        INT4Array		*B
);

void LALDArrayDotStarI4Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT4Array		*B
);

void LALDDotStarI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        REAL8			B,
        INT8Vector		*A
);

void LALDVectorDotStarI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        INT8Vector		*A
);

void LALDDotStarI8Array (
        LALStatus               *status,
        INT8Array		**result,
        REAL8			A,
        INT8Array		*B
);

void LALDArrayDotStarI8Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT8Array		*B
);

void LALSDotStarSVector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4			B,
        REAL4Vector		*A
);

void LALSVectorDotStarSVector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*B,
        REAL4Vector		*A
);

void LALSDotStarSArray (
        LALStatus               *status,
        REAL4Array		**result,
        REAL4			A,
        REAL4Array		*B
);

void LALSArrayDotStarSArray (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        REAL4Array		*B
);

void LALSDotStarI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        REAL4			B,
        INT2Vector		*A
);

void LALSVectorDotStarI2Vector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*B,
        INT2Vector		*A
);

void LALSDotStarI2Array (
        LALStatus               *status,
        INT2Array		**result,
        REAL4			A,
        INT2Array		*B
);

void LALSArrayDotStarI2Array (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT2Array		*B
);

void LALSDotStarI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        REAL4			B,
        INT4Vector		*A
);

void LALSVectorDotStarI4Vector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*B,
        INT4Vector		*A
);

void LALSDotStarI4Array (
        LALStatus               *status,
        INT4Array		**result,
        REAL4			A,
        INT4Array		*B
);

void LALSArrayDotStarI4Array (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT4Array		*B
);

void LALSDotStarI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        REAL4			B,
        INT8Vector		*A
);

void LALSVectorDotStarI8Vector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*B,
        INT8Vector		*A
);

void LALSDotStarI8Array (
        LALStatus               *status,
        INT8Array		**result,
        REAL4			A,
        INT8Array		*B
);

void LALSArrayDotStarI8Array (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT8Array		*B
);

void LALI2DotStarI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT2			B,
        INT2Vector		*A
);

void LALI2VectorDotStarI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT2Vector		*B,
        INT2Vector		*A
);

void LALI2DotStarI2Array (
        LALStatus               *status,
        INT2Array		**result,
        INT2			A,
        INT2Array		*B
);

void LALI2ArrayDotStarI2Array (
        LALStatus		*status,
        INT2Array		**result,
        INT2Array		*A,
        INT2Array		*B
);

void LALI4DotStarI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT4			B,
        INT2Vector		*A
);

void LALI4VectorDotStarI2Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT4Vector		*B,
        INT2Vector		*A
);

void LALI4DotStarI2Array (
        LALStatus               *status,
        INT2Array		**result,
        INT4			A,
        INT2Array		*B
);

void LALI4ArrayDotStarI2Array (
        LALStatus		*status,
        INT4Array		**result,
        INT4Array		*A,
        INT2Array		*B
);

void LALI4DotStarI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT4			B,
        INT4Vector		*A
);

void LALI4VectorDotStarI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT4Vector		*B,
        INT4Vector		*A
);

void LALI4DotStarI4Array (
        LALStatus               *status,
        INT4Array		**result,
        INT4			A,
        INT4Array		*B
);

void LALI4ArrayDotStarI4Array (
        LALStatus		*status,
        INT4Array		**result,
        INT4Array		*A,
        INT4Array		*B
);

void LALI8DotStarI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT8			B,
        INT2Vector		*A
);

void LALI8VectorDotStarI2Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8Vector		*B,
        INT2Vector		*A
);

void LALI8DotStarI2Array (
        LALStatus               *status,
        INT2Array		**result,
        INT8			A,
        INT2Array		*B
);

void LALI8ArrayDotStarI2Array (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT2Array		*B
);

void LALI8DotStarI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT8			B,
        INT4Vector		*A
);

void LALI8VectorDotStarI4Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8Vector		*B,
        INT4Vector		*A
);

void LALI8DotStarI4Array (
        LALStatus               *status,
        INT4Array		**result,
        INT8			A,
        INT4Array		*B
);

void LALI8ArrayDotStarI4Array (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT4Array		*B
);

void LALI8DotStarI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8			B,
        INT8Vector		*A
);

void LALI8VectorDotStarI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8Vector		*B,
        INT8Vector		*A
);

void LALI8DotStarI8Array (
        LALStatus               *status,
        INT8Array		**result,
        INT8			A,
        INT8Array		*B
);

void LALI8ArrayDotStarI8Array (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT8Array		*B
);

/* Matrix Division */

void LALDDotSlashDVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			A,
        REAL8Vector		*B
);

void LALDVectorDotSlashD (
        LALStatus               *status,
        REAL8Vector          **result,
        REAL8Vector          *A,
        REAL8                    B
);

void LALDVectorDotSlashDVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        REAL8Vector		*B
);

void LALDDotSlashDArray (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        REAL8Array		*B
);

void LALDArrayDotSlashD (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL8			B
);

void LALDArrayDotSlashDArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL8Array		*B
);

void LALDDotSlashSVector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL8			A,
        REAL4Vector		*B
);

void LALDVectorDotSlashS (
        LALStatus               *status,
        REAL4Vector          **result,
        REAL4Vector          *A,
        REAL8                    B
);

void LALDVectorDotSlashSVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        REAL4Vector		*B
);

void LALDDotSlashSArray (
        LALStatus               *status,
        REAL4Array		**result,
        REAL8			A,
        REAL4Array		*B
);

void LALDArrayDotSlashS (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL4			B
);

void LALDArrayDotSlashSArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL4Array		*B
);

void LALDDotSlashI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        REAL8			A,
        INT2Vector		*B
);

void LALDVectorDotSlashI2 (
        LALStatus               *status,
        INT2Vector          **result,
        INT2Vector          *A,
        REAL8                    B
);

void LALDVectorDotSlashI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        INT2Vector		*B
);

void LALDDotSlashI2Array (
        LALStatus               *status,
        INT2Array		**result,
        REAL8			A,
        INT2Array		*B
);

void LALDArrayDotSlashI2 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT2			B
);

void LALDArrayDotSlashI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT2Array		*B
);

void LALDDotSlashI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        REAL8			A,
        INT4Vector		*B
);

void LALDVectorDotSlashI4 (
        LALStatus               *status,
        INT4Vector          **result,
        INT4Vector          *A,
        REAL8                    B
);

void LALDVectorDotSlashI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        INT4Vector		*B
);

void LALDDotSlashI4Array (
        LALStatus               *status,
        INT4Array		**result,
        REAL8			A,
        INT4Array		*B
);

void LALDArrayDotSlashI4 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT4			B
);

void LALDArrayDotSlashI4Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT4Array		*B
);

void LALDDotSlashI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        REAL8			A,
        INT8Vector		*B
);

void LALDVectorDotSlashI8 (
        LALStatus               *status,
        INT8Vector          **result,
        INT8Vector          *A,
        REAL8                    B
);

void LALDVectorDotSlashI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        INT8Vector		*B
);

void LALDDotSlashI8Array (
        LALStatus               *status,
        INT8Array		**result,
        REAL8			A,
        INT8Array		*B
);

void LALDArrayDotSlashI8 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT8			B
);

void LALDArrayDotSlashI8Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT8Array		*B
);

void LALSDotSlashSVector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4			A,
        REAL4Vector		*B
);

void LALSVectorDotSlashS (
        LALStatus               *status,
        REAL4Vector          **result,
        REAL4Vector          *A,
        REAL4                    B
);

void LALSVectorDotSlashSVector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*A,
        REAL4Vector		*B
);

void LALSDotSlashSArray (
        LALStatus               *status,
        REAL4Array		**result,
        REAL4			A,
        REAL4Array		*B
);

void LALSArrayDotSlashS (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        REAL4			B
);

void LALSArrayDotSlashSArray (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        REAL4Array		*B
);

void LALSDotSlashI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        REAL4			A,
        INT2Vector		*B
);

void LALSVectorDotSlashI2 (
        LALStatus               *status,
        INT2Vector          **result,
        INT2Vector          *A,
        REAL4                    B
);

void LALSVectorDotSlashI2Vector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*A,
        INT2Vector		*B
);

void LALSDotSlashI2Array (
        LALStatus               *status,
        INT2Array		**result,
        REAL4			A,
        INT2Array		*B
);

void LALSArrayDotSlashI2 (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT2			B
);

void LALSArrayDotSlashI2Array (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT2Array		*B
);

void LALSDotSlashI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        REAL4			A,
        INT4Vector		*B
);

void LALSVectorDotSlashI4 (
        LALStatus               *status,
        INT4Vector          **result,
        INT4Vector          *A,
        REAL4                    B
);

void LALSVectorDotSlashI4Vector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*A,
        INT4Vector		*B
);

void LALSDotSlashI4Array (
        LALStatus               *status,
        INT4Array		**result,
        REAL4			A,
        INT4Array		*B
);

void LALSArrayDotSlashI4 (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT4			B
);

void LALSArrayDotSlashI4Array (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT4Array		*B
);

void LALSDotSlashI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        REAL4			A,
        INT8Vector		*B
);

void LALSVectorDotSlashI8 (
        LALStatus               *status,
        INT8Vector          **result,
        INT8Vector          *A,
        REAL4                    B
);

void LALSVectorDotSlashI8Vector (
        LALStatus		*status,
        REAL4Vector		**result,
        REAL4Vector		*A,
        INT8Vector		*B
);

void LALSDotSlashI8Array (
        LALStatus               *status,
        INT8Array		**result,
        REAL4			A,
        INT8Array		*B
);

void LALSArrayDotSlashI8 (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT8			B
);

void LALSArrayDotSlashI8Array (
        LALStatus		*status,
        REAL4Array		**result,
        REAL4Array		*A,
        INT8Array		*B
);

void LALI2DotSlashI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT2			A,
        INT2Vector		*B
);

void LALI2VectorDotSlashI2 (
        LALStatus               *status,
        INT2Vector          **result,
        INT2Vector          *A,
        INT2                    B
);

void LALI2VectorDotSlashI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT2Vector		*A,
        INT2Vector		*B
);

void LALI2DotSlashI2Array (
        LALStatus               *status,
        INT2Array		**result,
        INT2			A,
        INT2Array		*B
);

void LALI2ArrayDotSlashI2 (
        LALStatus		*status,
        INT2Array		**result,
        INT2Array		*A,
        INT2			B
);

void LALI2ArrayDotSlashI2Array (
        LALStatus		*status,
        INT2Array		**result,
        INT2Array		*A,
        INT2Array		*B
);

void LALI4DotSlashI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT4			A,
        INT2Vector		*B
);

void LALI4VectorDotSlashI2 (
        LALStatus               *status,
        INT2Vector          **result,
        INT2Vector          *A,
        INT4                    B
);

void LALI4VectorDotSlashI2Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT4Vector		*A,
        INT2Vector		*B
);

void LALI4DotSlashI2Array (
        LALStatus               *status,
        INT2Array		**result,
        INT4			A,
        INT2Array		*B
);

void LALI4ArrayDotSlashI2 (
        LALStatus		*status,
        INT4Array		**result,
        INT4Array		*A,
        INT2			B
);

void LALI4ArrayDotSlashI2Array (
        LALStatus		*status,
        INT4Array		**result,
        INT4Array		*A,
        INT2Array		*B
);

void LALI4DotSlashI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT4			A,
        INT4Vector		*B
);

void LALI4VectorDotSlashI4 (
        LALStatus               *status,
        INT4Vector          **result,
        INT4Vector          *A,
        INT4                    B
);

void LALI4VectorDotSlashI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT4Vector		*A,
        INT4Vector		*B
);

void LALI4DotSlashI4Array (
        LALStatus               *status,
        INT4Array		**result,
        INT4			A,
        INT4Array		*B
);

void LALI4ArrayDotSlashI4 (
        LALStatus		*status,
        INT4Array		**result,
        INT4Array		*A,
        INT4			B
);

void LALI4ArrayDotSlashI4Array (
        LALStatus		*status,
        INT4Array		**result,
        INT4Array		*A,
        INT4Array		*B
);

void LALI8DotSlashI2Vector (
        LALStatus		*status,
        INT2Vector		**result,
        INT8			A,
        INT2Vector		*B
);

void LALI8VectorDotSlashI2 (
        LALStatus               *status,
        INT2Vector          **result,
        INT2Vector          *A,
        INT8                    B
);

void LALI8VectorDotSlashI2Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8Vector		*A,
        INT2Vector		*B
);

void LALI8DotSlashI2Array (
        LALStatus               *status,
        INT2Array		**result,
        INT8			A,
        INT2Array		*B
);

void LALI8ArrayDotSlashI2 (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT2			B
);

void LALI8ArrayDotSlashI2Array (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT2Array		*B
);

void LALI8DotSlashI4Vector (
        LALStatus		*status,
        INT4Vector		**result,
        INT8			A,
        INT4Vector		*B
);

void LALI8VectorDotSlashI4 (
        LALStatus               *status,
        INT4Vector          **result,
        INT4Vector          *A,
        INT8                    B
);

void LALI8VectorDotSlashI4Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8Vector		*A,
        INT4Vector		*B
);

void LALI8DotSlashI4Array (
        LALStatus               *status,
        INT4Array		**result,
        INT8			A,
        INT4Array		*B
);

void LALI8ArrayDotSlashI4 (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT4			B
);

void LALI8ArrayDotSlashI4Array (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT4Array		*B
);

void LALI8DotSlashI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8			A,
        INT8Vector		*B
);

void LALI8VectorDotSlashI8 (
        LALStatus               *status,
        INT8Vector          **result,
        INT8Vector          *A,
        INT8                    B
);

void LALI8VectorDotSlashI8Vector (
        LALStatus		*status,
        INT8Vector		**result,
        INT8Vector		*A,
        INT8Vector		*B
);

void LALI8DotSlashI8Array (
        LALStatus               *status,
        INT8Array		**result,
        INT8			A,
        INT8Array		*B
);

void LALI8ArrayDotSlashI8 (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT8			B
);

void LALI8ArrayDotSlashI8Array (
        LALStatus		*status,
        INT8Array		**result,
        INT8Array		*A,
        INT8Array		*B
);

/* Matrix Power */

void LALDDotPowerDVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			A,
        REAL8Vector		*B
);

void LALDVectorDotPowerD (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        REAL8			B
);

void LALDVectorDotPowerDVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        REAL8Vector		*A
);

void LALDDotPowerDArray (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        REAL8Array		*B
);

void LALDArrayDotPowerD (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL8			B
);

void LALDArrayDotPowerDArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL8Array		*B
);

void LALDDotPowerSVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			A,
        REAL4Vector		*B
);

void LALDVectorDotPowerS (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        REAL4			B
);

void LALDVectorDotPowerSVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        REAL4Vector		*A
);

void LALDDotPowerSArray (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        REAL4Array		*B
);

void LALDArrayDotPowerS (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL4			B
);

void LALDArrayDotPowerSArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        REAL4Array		*B
);

void LALDDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			A,
        INT2Vector		*B
);

void LALDVectorDotPowerI2 (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        INT2			B
);

void LALDVectorDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        INT2Vector		*A
);

void LALDDotPowerI2Array (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        INT2Array		*B
);

void LALDArrayDotPowerI2 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT2			B
);

void LALDArrayDotPowerI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT2Array		*B
);

void LALDDotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			A,
        INT4Vector		*B
);

void LALDVectorDotPowerI4 (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        INT4			B
);

void LALDVectorDotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        INT4Vector		*A
);

void LALDDotPowerI4Array (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        INT4Array		*B
);

void LALDArrayDotPowerI4 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT4			B
);

void LALDArrayDotPowerI4Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT4Array		*B
);

void LALDDotPowerI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8			A,
        INT8Vector		*B
);

void LALDVectorDotPowerI8 (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL8Vector		*A,
        INT8			B
);

void LALDVectorDotPowerI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL8Vector		*B,
        INT8Vector		*A
);

void LALDDotPowerI8Array (
        LALStatus               *status,
        REAL8Array		**result,
        REAL8			A,
        INT8Array		*B
);

void LALDArrayDotPowerI8 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT8			B
);

void LALDArrayDotPowerI8Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL8Array		*A,
        INT8Array		*B
);

void LALSDotPowerSVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4			A,
        REAL4Vector		*B
);

void LALSVectorDotPowerS (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL4Vector		*A,
        REAL4			B
);

void LALSVectorDotPowerSVector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4Vector		*B,
        REAL4Vector		*A
);

void LALSDotPowerSArray (
        LALStatus               *status,
        REAL8Array		**result,
        REAL4			A,
        REAL4Array		*B
);

void LALSArrayDotPowerS (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        REAL4			B
);

void LALSArrayDotPowerSArray (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        REAL4Array		*B
);

void LALSDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4			A,
        INT2Vector		*B
);

void LALSVectorDotPowerI2 (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL4Vector		*A,
        INT2			B
);

void LALSVectorDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4Vector		*B,
        INT2Vector		*A
);

void LALSDotPowerI2Array (
        LALStatus               *status,
        REAL8Array		**result,
        REAL4			A,
        INT2Array		*B
);

void LALSArrayDotPowerI2 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        INT2			B
);

void LALSArrayDotPowerI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        INT2Array		*B
);

void LALSDotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4			A,
        INT4Vector		*B
);

void LALSVectorDotPowerI4 (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL4Vector		*A,
        INT4			B
);

void LALSVectorDotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4Vector		*B,
        INT4Vector		*A
);

void LALSDotPowerI4Array (
        LALStatus               *status,
        REAL8Array		**result,
        REAL4			A,
        INT4Array		*B
);

void LALSArrayDotPowerI4 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        INT4			B
);

void LALSArrayDotPowerI4Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        INT4Array		*B
);

void LALSDotPowerI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4			A,
        INT8Vector		*B
);

void LALSVectorDotPowerI8 (
        LALStatus               *status,
        REAL8Vector		**result,
        REAL4Vector		*A,
        INT8			B
);

void LALSVectorDotPowerI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        REAL4Vector		*B,
        INT8Vector		*A
);

void LALSDotPowerI8Array (
        LALStatus               *status,
        REAL8Array		**result,
        REAL4			A,
        INT8Array		*B
);

void LALSArrayDotPowerI8 (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        INT8			B
);

void LALSArrayDotPowerI8Array (
        LALStatus		*status,
        REAL8Array		**result,
        REAL4Array		*A,
        INT8Array		*B
);

void LALI2DotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT2			A,
        INT2Vector		*B
);

void LALI2VectorDotPowerI2 (
        LALStatus               *status,
        REAL8Vector		**result,
        INT2Vector		*A,
        INT2			B
);

void LALI2VectorDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT2Vector		*B,
        INT2Vector		*A
);

void LALI2DotPowerI2Array (
        LALStatus               *status,
        REAL8Array		**result,
        INT2			A,
        INT2Array		*B
);

void LALI2ArrayDotPowerI2 (
        LALStatus		*status,
        REAL8Array		**result,
        INT2Array		*A,
        INT2			B
);

void LALI2ArrayDotPowerI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        INT2Array		*A,
        INT2Array		*B
);

void LALI4DotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT4			A,
        INT2Vector		*B
);

void LALI4VectorDotPowerI2 (
        LALStatus               *status,
        REAL8Vector		**result,
        INT4Vector		*A,
        INT2			B
);

void LALI4VectorDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT4Vector		*B,
        INT2Vector		*A
);

void LALI4DotPowerI2Array (
        LALStatus               *status,
        REAL8Array		**result,
        INT4			A,
        INT2Array		*B
);

void LALI4ArrayDotPowerI2 (
        LALStatus		*status,
        REAL8Array		**result,
        INT4Array		*A,
        INT2			B
);

void LALI4ArrayDotPowerI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        INT4Array		*A,
        INT2Array		*B
);

void LALI4DotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT4			A,
        INT4Vector		*B
);

void LALI4VectorDotPowerI4 (
        LALStatus               *status,
        REAL8Vector		**result,
        INT4Vector		*A,
        INT4			B
);

void LALI4VectorDotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT4Vector		*B,
        INT4Vector		*A
);

void LALI4DotPowerI4Array (
        LALStatus               *status,
        REAL8Array		**result,
        INT4			A,
        INT4Array		*B
);

void LALI4ArrayDotPowerI4 (
        LALStatus		*status,
        REAL8Array		**result,
        INT4Array		*A,
        INT4			B
);

void LALI4ArrayDotPowerI4Array (
        LALStatus		*status,
        REAL8Array		**result,
        INT4Array		*A,
        INT4Array		*B
);

void LALI8DotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT8			A,
        INT2Vector		*B
);

void LALI8VectorDotPowerI2 (
        LALStatus               *status,
        REAL8Vector		**result,
        INT8Vector		*A,
        INT2			B
);

void LALI8VectorDotPowerI2Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT8Vector		*B,
        INT2Vector		*A
);

void LALI8DotPowerI2Array (
        LALStatus               *status,
        REAL8Array		**result,
        INT8			A,
        INT2Array		*B
);

void LALI8ArrayDotPowerI2 (
        LALStatus		*status,
        REAL8Array		**result,
        INT8Array		*A,
        INT2			B
);

void LALI8ArrayDotPowerI2Array (
        LALStatus		*status,
        REAL8Array		**result,
        INT8Array		*A,
        INT2Array		*B
);

void LALI8DotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT8			A,
        INT4Vector		*B
);

void LALI8VectorDotPowerI4 (
        LALStatus               *status,
        REAL8Vector		**result,
        INT8Vector		*A,
        INT4			B
);

void LALI8VectorDotPowerI4Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT8Vector		*B,
        INT4Vector		*A
);

void LALI8DotPowerI4Array (
        LALStatus               *status,
        REAL8Array		**result,
        INT8			A,
        INT4Array		*B
);

void LALI8ArrayDotPowerI4 (
        LALStatus		*status,
        REAL8Array		**result,
        INT8Array		*A,
        INT4			B
);

void LALI8ArrayDotPowerI4Array (
        LALStatus		*status,
        REAL8Array		**result,
        INT8Array		*A,
        INT4Array		*B
);

void LALI8DotPowerI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT8			A,
        INT8Vector		*B
);

void LALI8VectorDotPowerI8 (
        LALStatus               *status,
        REAL8Vector		**result,
        INT8Vector		*A,
        INT8			B
);

void LALI8VectorDotPowerI8Vector (
        LALStatus		*status,
        REAL8Vector		**result,
        INT8Vector		*B,
        INT8Vector		*A
);

void LALI8DotPowerI8Array (
        LALStatus               *status,
        REAL8Array		**result,
        INT8			A,
        INT8Array		*B
);

void LALI8ArrayDotPowerI8 (
        LALStatus		*status,
        REAL8Array		**result,
        INT8Array		*A,
        INT8			B
);

void LALI8ArrayDotPowerI8Array (
        LALStatus		*status,
        REAL8Array		**result,
        INT8Array		*A,
        INT8Array		*B
);

/* Miscellaneous Routines */

void LALDCumSum (
	LALStatus		*status,
	REAL8Vector			**result,
	REAL8Vector		*data
);

void LALDSum (
        LALStatus		*status,
        REAL8			*result,
        REAL8Vector		*data
);

void LALDMax (
        LALStatus		*status,
        REAL8			*result,
        REAL8Vector		*data,
	INT4			*myindex
);

void LALDFlipVector (
        LALStatus		*status,
        REAL8Vector			**result,
        REAL8Vector           *data
);

void LALSCumSum (
	LALStatus		*status,
	REAL4Vector			**result,
	REAL4Vector		*data
);

void LALSSum (
        LALStatus		*status,
        REAL4			*result,
        REAL4Vector		*data
);

void LALSMax (
        LALStatus		*status,
        REAL4			*result,
        REAL4Vector		*data,
	INT4			*myindex
);

void LALSFlipVector (
        LALStatus		*status,
        REAL4Vector			**result,
        REAL4Vector           *data
);

void LALI2CumSum (
	LALStatus		*status,
	INT2Vector			**result,
	INT2Vector		*data
);

void LALI2Sum (
        LALStatus		*status,
        INT2			*result,
        INT2Vector		*data
);

void LALI2Max (
        LALStatus		*status,
        INT2			*result,
        INT2Vector		*data,
	INT4			*myindex
);

void LALI2FlipVector (
        LALStatus		*status,
        INT2Vector			**result,
        INT2Vector           *data
);

void LALI4CumSum (
	LALStatus		*status,
	INT4Vector			**result,
	INT4Vector		*data
);

void LALI4Sum (
        LALStatus		*status,
        INT4			*result,
        INT4Vector		*data
);

void LALI4Max (
        LALStatus		*status,
        INT4			*result,
        INT4Vector		*data,
	INT4			*myindex
);

void LALI4FlipVector (
        LALStatus		*status,
        INT4Vector			**result,
        INT4Vector           *data
);

void LALI8CumSum (
	LALStatus		*status,
	INT8Vector			**result,
	INT8Vector		*data
);

void LALI8Sum (
        LALStatus		*status,
        INT8			*result,
        INT8Vector		*data
);

void LALI8Max (
        LALStatus		*status,
        INT8			*result,
        INT8Vector		*data,
	INT4			*myindex
);

void LALI8FlipVector (
        LALStatus		*status,
        INT8Vector			**result,
        INT8Vector           *data
);

#ifdef  __cplusplus
}
#endif /* C++ protection */

#endif
