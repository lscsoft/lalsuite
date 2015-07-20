/*
*  Copyright (C) 2007 Teviet Creighton
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

#ifndef _MATRIXUTILS_H
#define _MATRIXUTILS_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup MatrixUtils_h Header MatrixUtils.h
 * \ingroup lal_utilities
 * \author Creighton, T. D.
 *
 * \brief Provides routines to solve linear systems.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/MatrixUtils.h>
 * \endcode
 *
 * This header covers routines to solve linear systems of equations and
 * eigensystems, using algorithms adapted from Chapters~2 and~11 of
 * Numerical Recipes \cite ptvf1992 .  The only routines at present are
 * for computing eigenvalues and eigenvectors of real symmetric matrices.
 * Routines for inverting or computing the determinant of arbitrary
 * square matrices will likely follow.
 *
 * ### Notation ###
 *
 * A \e matrix is represented in LAL by a <tt>\<datatype\>Array</tt>
 * structure with a <tt>dimLength-\>length</tt> field of 2; the
 * <tt>dimLength-\>data</tt> field gives the dimensions \f$[M,N]\f$ of the
 * matrix.  Using the place-index notation common in tensor calculus, a
 * matrix is a two-index tensor:
 * \f{equation}{
 * \mathsf{A}^a{}_b = \left[\begin{array}{cccc}
 * A^1{}_1 & A^1{}_2 & \cdots & A^1{}_N \\
 * A^2{}_1 & A^2{}_2 & \cdots & A^2{}_N \\
 * \vdots  & \vdots  & \ddots & \vdots  \\
 * A^M{}_1 & A^M{}_2 & \cdots & A^M{}_N
 * \end{array}\right] \;,
 * \f}
 * that is, the first (raised) index represents the row number and the
 * second (lowered) index the column number.  The standard C array
 * structure declares this object as, say, <tt>float a[M][N]</tt>, where
 * the element \f$A^i{}_j\f$ is stored in <tt>a[i][j]</tt>.  The LAL array
 * structure <tt>REAL4Array a</tt> stores data in a flattened block of
 * memory, where the element \f$A^i{}_j\f$ is stored in <tt>a.data[i*M+j]</tt>.
 *
 * A <em>row vector</em> is a matrix with only one row (\f$M=1\f$).  In the
 * above place-index notation, it is represented with a single lowered
 * index:
 * \f{equation}{
 * \mathsf{r}_a = \left[\begin{array}{cccc} r_1 & r_2 & \cdots & r_N
 * \end{array}\right] \;.
 * \f}
 * A <em>column vector</em> is a matrix with only one column (\f$N=1\f$).  In
 * the above place-index notation, it is represented with a single raised
 * index:
 * \f{equation}{
 * \mathsf{c}^a = \left[\begin{array}{c} c^1 \\ c^2 \\ \vdots \\ c^M
 * \end{array}\right] \;.
 * \f}
 * In LAL, both of these objects are conventionally represented as a LAL
 * vector structure.  Whether the object is to be used as a row or column
 * vector must be determined from context; it is not specified by the
 * datatype.
 *
 * ### Properties ###
 *
 * The basic matrix operations are addition, scalar multiplication, and
 * vector multiplication.  We assume the reader is familiar with these.
 * In addition, we will refer to the following unary operations on
 * \e square matrices:
 *
 * \e Inversion: The inverse \f$(\mathsf{A}^{-1}){}^a{}_b\f$ of a
 * matrix \f$\mathsf{A}^a{}_b\f$ is one such that their matrix product is the
 * identity matrix \f$\delta^a{}_b\f$ (whose elements \f$\delta^i{}_j\f$ are just
 * the Kronecker delta function).
 *
 * \e Transposition: The transpose \f$(\mathsf{A}^T){}^a{}_b\f$ of a
 * matrix \f$\mathsf{A}^a{}_b\f$ is given by interchanging the indecies on
 * each component: \f$(A^T){}^i{}_j=A^j{}_i\f$.
 *
 * \e Conjugation: The Hermitian conjugate (adjoint)
 * \f$(\mathsf{A}^\dag){}^a{}_b\f$ of a matrix \f$\mathsf{A}^a{}_b\f$ is given by
 * interchanging the indecies and taking the complex conjugate of each
 * component: \f$(A^\dag){}^i{}_j={A^j{}_i}^*\f$.
 *
 * A matrix that is identical to its own transpose is called
 * \e symmetric.  A matrix whose transpose is identical to the
 * original matrix's inverse is called \e orthogonal.  A matrix that
 * is identical to its own Hermitian conjugate is called \e Hermitian
 * (or <em>self-adjoint</em>.  A matrix whose Hermitian conjugate is
 * identical to the original matrix's inverse is called \e unitary.
 *
 * At present, the routines under this header only deal with \e real
 * matrices (i.e.\ matrices, vectors, and scalars whose components are
 * all real).  In this case, symmetric is equivalent to Hermitian, and
 * orthogonal is equivalent to unitary.
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define MATRIXUTILSH_ENUL  1		/**< Unexpected null pointer in arguments */
#define MATRIXUTILSH_EDIM  2		/**< Bad matrix dimensions */
#define MATRIXUTILSH_EITER 3		/**< Did not converge after maximum iterations */
#define MATRIXUTILSH_ESING 4		/**< Singular matrix */
#define MATRIXUTILSH_EMEM  5		/**< Memory allocation error */
/*@}*/

/*@}*/

#define MATRIXUTILSH_MSGENUL  "Unexpected null pointer in arguments"
#define MATRIXUTILSH_MSGEDIM  "Bad matrix dimensions"
#define MATRIXUTILSH_MSGEITER "Did not converge after maximum iterations"
#define MATRIXUTILSH_MSGESING "Singular matrix"
#define MATRIXUTILSH_MSGEMEM  "Memory allocation error"


/* ---------- Function prototypes. ---------- */
/** \addtogroup MatrixOps_c */
/*@{*/
void
LALDMatrixMultiply( LALStatus *, REAL8Array *out, REAL8Array *in1, REAL8Array *in2 );

void
LALDMatrixTranspose( LALStatus *, REAL8Array *out, REAL8Array *in1 );
/*@}*/


void
LALSMatrixInverse( LALStatus *, REAL4 *det, REAL4Array *matrix, REAL4Array *inverse );

void
LALDMatrixDeterminant( LALStatus *, REAL8 *det, REAL8Array *matrix );

void
LALDMatrixInverse( LALStatus *, REAL8 *det, REAL8Array *matrix, REAL8Array *inverse );

void
LALDMatrixDeterminantErr( LALStatus *, REAL8 det[2], REAL8Array *matrix, REAL8Array *matrixErr );


void
LALSSymmetricEigenVectors( LALStatus *, REAL4Vector *values, REAL4Array *matrix );

void
LALDSymmetricEigenVectors( LALStatus *, REAL8Vector *values, REAL8Array *matrix );


#ifdef  __cplusplus
}
#endif

#endif /* _MATRIXUTILS_H */
