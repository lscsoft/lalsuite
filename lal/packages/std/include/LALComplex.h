/** \cond DONT_DOXYGEN */
/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

#ifndef _LALCOMPLEX_H
#define _LALCOMPLEX_H

/* include only if using old LAL complex struct types */
#ifdef LAL_USE_OLD_COMPLEX_STRUCTS

#include <lal/LALAtomicDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

COMPLEX16 XLALCOMPLEX16Add (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16AddImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16AddReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccos (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArccosReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccosh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArccoshReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccot (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arccoth (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arccsc (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArccscReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccsch (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arcsec (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArcsecReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arcsech (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arcsin (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArcsinReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arcsinh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arctan (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arctanh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArctanhReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Conjugate (COMPLEX16 z);
COMPLEX16 XLALCOMPLEX16Cos (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Cosh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Cot (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Coth (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Csc (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Csch (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Div (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16DivImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16DivReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16Exp (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Inverse (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Log (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Log10 (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16LogB (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16Mul (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16MulImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16MulReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16Negative (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Polar (REAL8 x, REAL8 y);
COMPLEX16 XLALCOMPLEX16Pow (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16PowReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16Rect (REAL8 x, REAL8 y);
COMPLEX16 XLALCOMPLEX16Sec (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Sech (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Sin (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Sinh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Sqrt (COMPLEX16 z);
COMPLEX16 XLALCOMPLEX16SqrtReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Sub (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16SubImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16SubReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16Tan (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Tanh (COMPLEX16 a);
COMPLEX8 XLALCOMPLEX8Add (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8AddImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8AddReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccos (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArccosReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccosh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArccoshReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccot (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arccoth (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arccsc (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArccscReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccsch (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arcsec (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArcsecReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arcsech (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arcsin (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArcsinReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arcsinh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arctan (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arctanh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArctanhReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Conjugate (COMPLEX8 z);
COMPLEX8 XLALCOMPLEX8Cos (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Cosh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Cot (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Coth (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Csc (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Csch (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Div (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8DivImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8DivReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8Exp (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Inverse (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Log (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Log10 (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8LogB (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8Mul (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8MulImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8MulReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8Negative (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Polar (REAL4 x, REAL4 y);
COMPLEX8 XLALCOMPLEX8Pow (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8PowReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8Rect (REAL4 x, REAL4 y);
COMPLEX8 XLALCOMPLEX8Sec (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Sech (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Sin (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Sinh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Sqrt (COMPLEX8 z);
COMPLEX8 XLALCOMPLEX8SqrtReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Sub (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8SubImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8SubReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8Tan (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Tanh (COMPLEX8 a);
REAL4 XLALCOMPLEX8Abs (COMPLEX8 z);
REAL4 XLALCOMPLEX8Abs2 (COMPLEX8 z);
REAL4 XLALCOMPLEX8Arg (COMPLEX8 z);
REAL4 XLALCOMPLEX8LogAbs (COMPLEX8 z);
REAL8 XLALCOMPLEX16Abs (COMPLEX16 z);
REAL8 XLALCOMPLEX16Abs2 (COMPLEX16 z);
REAL8 XLALCOMPLEX16Arg (COMPLEX16 z);
REAL8 XLALCOMPLEX16LogAbs (COMPLEX16 z);

// keep these macros as a reminder of what is the equivalent C99 function
/* #define cabs2(z) (XLALCOMPLEX16Abs2(z)) */
/* #define cabs2f(z) (XLALCOMPLEX8Abs2(z)) */
/* #define cacoshr(x) (XLALCOMPLEX16ArccoshReal(x)) */
/* #define cacoshrf(x) (XLALCOMPLEX8ArccoshReal(x)) */
/* #define cacosr(x) (XLALCOMPLEX16ArccosReal(x)) */
/* #define cacosrf(x) (XLALCOMPLEX8ArccosReal(x)) */
/* #define cacot(a) (XLALCOMPLEX16Arccot(a)) */
/* #define cacotf(a) (XLALCOMPLEX8Arccot(a)) */
/* #define cacoth(a) (XLALCOMPLEX16Arccoth(a)) */
/* #define cacothf(a) (XLALCOMPLEX8Arccoth(a)) */
/* #define cacsc(a) (XLALCOMPLEX16Arccsc(a)) */
/* #define cacscf(a) (XLALCOMPLEX8Arccsc(a)) */
/* #define cacsch(a) (XLALCOMPLEX16Arccsch(a)) */
/* #define cacschf(a) (XLALCOMPLEX8Arccsch(a)) */
/* #define cacscr(x) (XLALCOMPLEX16ArccscReal(x)) */
/* #define cacscrf(x) (XLALCOMPLEX8ArccscReal(x)) */
/* #define cadd(a,b) (XLALCOMPLEX16Add((a),(b))) */
/* #define caddf(a,b) (XLALCOMPLEX8Add((a),(b))) */
/* #define caddi(a,y) (XLALCOMPLEX16AddImag((a),(y))) */
/* #define caddif(a,y) (XLALCOMPLEX8AddImag((a),(y))) */
/* #define caddr(a,x) (XLALCOMPLEX16AddReal((a),(x))) */
/* #define caddrf(a,x) (XLALCOMPLEX8AddReal((a),(x))) */
/* #define casec(a) (XLALCOMPLEX16Arcsec(a)) */
/* #define casecf(a) (XLALCOMPLEX8Arcsec(a)) */
/* #define casech(a) (XLALCOMPLEX16Arcsech(a)) */
/* #define casechf(a) (XLALCOMPLEX8Arcsech(a)) */
/* #define casecr(x) (XLALCOMPLEX16ArcsecReal(x)) */
/* #define casecrf(x) (XLALCOMPLEX8ArcsecReal(x)) */
/* #define casinr(x) (XLALCOMPLEX16ArcsinReal(x)) */
/* #define casinrf(x) (XLALCOMPLEX8ArcsinReal(x)) */
/* #define catanhr(x) (XLALCOMPLEX16ArctanhReal(x)) */
/* #define catanhrf(x) (XLALCOMPLEX8ArctanhReal(x)) */
/* #define ccot(a) (XLALCOMPLEX16Cot(a)) */
/* #define ccotf(a) (XLALCOMPLEX8Cot(a)) */
/* #define ccoth(a) (XLALCOMPLEX16Coth(a)) */
/* #define ccothf(a) (XLALCOMPLEX8Coth(a)) */
/* #define ccsc(a) (XLALCOMPLEX16Csc(a)) */
/* #define ccscf(a) (XLALCOMPLEX8Csc(a)) */
/* #define ccsch(a) (XLALCOMPLEX16Csch(a)) */
/* #define ccschf(a) (XLALCOMPLEX8Csch(a)) */
/* #define cdiv(a,b) (XLALCOMPLEX16Div((a),(b))) */
/* #define cdivf(a,b) (XLALCOMPLEX8Div((a),(b))) */
/* #define cdivi(a,y) (XLALCOMPLEX16DivImag((a),(y))) */
/* #define cdivif(a,y) (XLALCOMPLEX8DivImag((a),(y))) */
/* #define cdivr(a,x) (XLALCOMPLEX16DivReal((a),(x))) */
/* #define cdivrf(a,x) (XLALCOMPLEX8DivReal((a),(x))) */
/* #define cinv(a) (XLALCOMPLEX16Inverse(a)) */
/* #define cinvf(a) (XLALCOMPLEX8Inverse(a)) */
/* #define clogabs(z) (XLALCOMPLEX16LogAbs(z)) */
/* #define clogabsf(z) (XLALCOMPLEX8LogAbs(z)) */
/* #define clogb(a,b) (XLALCOMPLEX16LogB((a),(b))) */
/* #define clogbf(a,b) (XLALCOMPLEX8LogB((a),(b))) */
/* #define cmul(a,b) (XLALCOMPLEX16Mul((a),(b))) */
/* #define cmulf(a,b) (XLALCOMPLEX8Mul((a),(b))) */
/* #define cmuli(a,y) (XLALCOMPLEX16MulImag((a),(y))) */
/* #define cmulif(a,y) (XLALCOMPLEX8MulImag((a),(y))) */
/* #define cmulr(a,x) (XLALCOMPLEX16MulReal((a),(x))) */
/* #define cmulrf(a,x) (XLALCOMPLEX8MulReal((a),(x))) */
/* #define cneg(a) (XLALCOMPLEX16Negative(a)) */
/* #define cnegf(a) (XLALCOMPLEX8Negative(a)) */
/* #define cpowr(a,x) (XLALCOMPLEX16PowReal((a),(x))) */
/* #define cpowrf(a,x) (XLALCOMPLEX8PowReal((a),(x))) */
/* #define csec(a) (XLALCOMPLEX16Sec(a)) */
/* #define csecf(a) (XLALCOMPLEX8Sec(a)) */
/* #define csech(a) (XLALCOMPLEX16Sech(a)) */
/* #define csechf(a) (XLALCOMPLEX8Sech(a)) */
/* #define csqrtr(x) (XLALCOMPLEX16SqrtReal(x)) */
/* #define csqrtrf(x) (XLALCOMPLEX8SqrtReal(x)) */
/* #define csub(a,b) (XLALCOMPLEX16Sub((a),(b))) */
/* #define csubf(a,b) (XLALCOMPLEX8Sub((a),(b))) */
/* #define csubi(a,y) (XLALCOMPLEX16SubImag((a),(y))) */
/* #define csubif(a,y) (XLALCOMPLEX8SubImag((a),(y))) */
/* #define csubr(a,x) (XLALCOMPLEX16SubReal((a),(x))) */
/* #define csubrf(a,x) (XLALCOMPLEX8SubReal((a),(x))) */

#ifdef __cplusplus
}
#endif

#endif /* LAL_USE_OLD_COMPLEX_STRUCTS */

#endif /* _LALCOMPLEX8_H */

/** \endcond */
