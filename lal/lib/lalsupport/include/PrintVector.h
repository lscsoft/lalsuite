/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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

#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H

#ifndef _LALSTDLIB_H
#include <lal/LALStdlib.h>
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifdef  __cplusplus
extern "C" {
#endif


/**
 * \addtogroup PrintVector_h
 * \author Allen, B.; generalized by J. T. Whelan <jtwhelan@loyno.edu>
 *
 * \brief This is a simple utility to print vectors into a file.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/PrintVector.h>
 * \endcode
 *
 * Contains the prototypes for the LAL <tt>\<datatype\>PrintVector</tt> functions
 *
 * \heading{Description}
 *
 * Each member of this family of functions prints the elements of
 * \f$\langle\mbox{datatype}\rangle\f$\c Vector into a file.  Note: the
 * file names are <tt>\<datatype\>PrintVector.000</tt>,
 * <tt>\<datatype\>PrintVector.001</tt>, and so on.
 * (<tt>\<datatype\></tt> is the abbreviation for the datatype,
 * included in the function names above.) The file numbers are
 * incremented with each additional call.  This function is for debugging
 * use only: it uses a static internal variable to keep track of the file
 * number so it should not be used in any real analysis codes.
 *
 * \heading{Notes}
 *
 * This function uses an internal static variable to keep track of file
 * numbers.  For this reason it should only be used for debugging
 * purposes in test functions, not in any production code.
 *
 * Additionally, since printf cannot handle INT8 as integers, the
 * functions <tt>LALI8PrintVector()</tt> and <tt>LALU8PrintVector()</tt> use
 * a typecast to REAL8 and are thus only valid for numbers between around
 * \f$-10^{15}\f$ and \f$10^{15}\f$.
 *
 * The output format is two or three space-separated columns: the first
 * column is the index of the element; the second is the element itself
 * for real and integer vectors and the real part of the element for
 * complex vectors; complex vectors have a third column containing the
 * imaginary part of the element.
 *
 */
/*@{*/
void LALCHARPrintVector( CHARVector *vector );
void LALI2PrintVector( INT2Vector *vector );
void LALI4PrintVector( INT4Vector *vector );
void LALI8PrintVector( INT8Vector *vector );
void LALU2PrintVector( UINT2Vector *vector );
void LALU4PrintVector( UINT4Vector *vector );
void LALU8PrintVector( UINT8Vector *vector );
void LALPrintVector( REAL4Vector *vector );
void LALSPrintVector( REAL4Vector *vector );
void LALDPrintVector( REAL8Vector *vector );
void LALCPrintVector( COMPLEX8Vector *vector );
void LALZPrintVector( COMPLEX16Vector *vector );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _PRINTVECTOR_H */
