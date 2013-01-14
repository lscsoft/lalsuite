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

#ifdef __cplusplus
extern "C" {
#endif

#include <lal/LALAtomicDatatypes.h>

#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#if defined __COMPLEX__ || defined _COMPLEX_H
#error "don't use both <complex.h> and LAL_USE_SHORT_MACROS"
#endif
#endif

#ifdef LAL_NO_COMPLEX_MACROS
#undef LAL_USE_COMPLEX_MACROS
#undef LAL_USE_COMPLEX_SHORT_MACROS
#else
#ifndef LAL_USE_COMPLEX_MACROS
#define LAL_USE_COMPLEX_MACROS /* default is to have these macros */
#endif
#endif


#define LAL_REAL(z) ((z).re)
#define LAL_IMAG(z) ((z).im)
#define LAL_COMPLEX_EQ(z1,z2) (((z1).re == (z2).re) && ((z1).im == (z2).im))

#define LAL_SET_COMPLEX(zp,x,y) do {(zp)->re=(x); (zp)->im=(y);} while(0)
#define LAL_SET_REAL(zp,x) do {(zp)->re=(x);} while(0)
#define LAL_SET_IMAG(zp,y) do {(zp)->im=(y);} while(0)

#define LAL_COMPLEX16_ONE (XLALCOMPLEX16Rect(1.0,0.0))
#define LAL_COMPLEX16_ZERO (XLALCOMPLEX16Rect(0.0,0.0))
#define LAL_COMPLEX16_NEGONE (XLALCOMPLEX16Rect(-1.0,0.0))
#define LAL_COMPLEX16_I (XLALCOMPLEX16Rect(0.0,1.0))

#ifdef LAL_USE_COMPLEX_MACROS
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define crect(x,y) (XLALCOMPLEX16Rect((x),(y)))
#define csetr(x) (XLALCOMPLEX16Rect((x),0.0))
#define cseti(y) (XLALCOMPLEX16Rect(0.0,(y)))
#define cpolar(x,y) (XLALCOMPLEX16Polar((x),(y)))
#define cisequal(z1,z2) (((z1).re == (z2).re) && ((z1).im == (z2).im))

#define cunit (XLALCOMPLEX16Rect(1.0,0.0))
#define czero (XLALCOMPLEX16Rect(0.0,0.0))
#define cnegone (XLALCOMPLEX16Rect(-1.0,0.0))
#endif
#endif

COMPLEX16 XLALCOMPLEX16Rect (REAL8 x, REAL8 y);
COMPLEX16 XLALCOMPLEX16Polar (REAL8 x, REAL8 y);

REAL8 XLALCOMPLEX16Arg (COMPLEX16 z);
REAL8 XLALCOMPLEX16Abs (COMPLEX16 z);
REAL8 XLALCOMPLEX16Abs2 (COMPLEX16 z);
REAL8 XLALCOMPLEX16LogAbs (COMPLEX16 z);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CARG(z) (XLALCOMPLEX16Arg(z))
#define LAL_CABS(z) (XLALCOMPLEX16Abs(z))
#define LAL_CABS2(z) (XLALCOMPLEX16Abs2(z))
#define LAL_CLOGABS(z) (XLALCOMPLEX16LogAbs(z))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define cabs2(z) (XLALCOMPLEX16Abs2(z))
#define clogabs(z) (XLALCOMPLEX16LogAbs(z))
#endif
#endif


COMPLEX16 XLALCOMPLEX16Add (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16Sub (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16Mul (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16Div (COMPLEX16 a, COMPLEX16 b);

COMPLEX16 XLALCOMPLEX16AddReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16SubReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16MulReal (COMPLEX16 a, REAL8 x);
COMPLEX16 XLALCOMPLEX16DivReal (COMPLEX16 a, REAL8 x);

COMPLEX16 XLALCOMPLEX16AddImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16SubImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16MulImag (COMPLEX16 a, REAL8 y);
COMPLEX16 XLALCOMPLEX16DivImag (COMPLEX16 a, REAL8 y);

COMPLEX16 XLALCOMPLEX16Conjugate (COMPLEX16 z);
COMPLEX16 XLALCOMPLEX16Inverse (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Negative (COMPLEX16 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CADD(a,b) (XLALCOMPLEX16Add((a),(b)))
#define LAL_CSUB(a,b) (XLALCOMPLEX16Sub((a),(b)))
#define LAL_CMUL(a,b) (XLALCOMPLEX16Mul((a),(b)))
#define LAL_CDIV(a,b) (XLALCOMPLEX16Div((a),(b)))

#define LAL_CADD_REAL(a,x) (XLALCOMPLEX16AddReal((a),(x)))
#define LAL_CSUB_REAL(a,x) (XLALCOMPLEX16SubReal((a),(x)))
#define LAL_CMUL_REAL(a,x) (XLALCOMPLEX16MulReal((a),(x)))
#define LAL_CDIV_REAL(a,x) (XLALCOMPLEX16DivReal((a),(x)))

#define LAL_CADD_IMAG(a,y) (XLALCOMPLEX16AddImag((a),(y)))
#define LAL_CSUB_IMAG(a,y) (XLALCOMPLEX16SubImag((a),(y)))
#define LAL_CMUL_IMAG(a,y) (XLALCOMPLEX16MulImag((a),(y)))
#define LAL_CDIV_IMAG(a,y) (XLALCOMPLEX16DivImag((a),(y)))

#define LAL_CONJ(z) (XLALCOMPLEX16Conjugate(z))
#define LAL_CINV(a) (XLALCOMPLEX16Inverse(a))
#define LAL_CNEG(a) (XLALCOMPLEX16Negative(a))

#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define cadd(a,b) (XLALCOMPLEX16Add((a),(b)))
#define csub(a,b) (XLALCOMPLEX16Sub((a),(b)))
#define cmul(a,b) (XLALCOMPLEX16Mul((a),(b)))
#define cdiv(a,b) (XLALCOMPLEX16Div((a),(b)))

#define caddr(a,x) (XLALCOMPLEX16AddReal((a),(x)))
#define csubr(a,x) (XLALCOMPLEX16SubReal((a),(x)))
#define cmulr(a,x) (XLALCOMPLEX16MulReal((a),(x)))
#define cdivr(a,x) (XLALCOMPLEX16DivReal((a),(x)))

#define caddi(a,y) (XLALCOMPLEX16AddImag((a),(y)))
#define csubi(a,y) (XLALCOMPLEX16SubImag((a),(y)))
#define cmuli(a,y) (XLALCOMPLEX16MulImag((a),(y)))
#define cdivi(a,y) (XLALCOMPLEX16DivImag((a),(y)))

#define cinv(a) (XLALCOMPLEX16Inverse(a))
#define cneg(a) (XLALCOMPLEX16Negative(a))
#endif
#endif


COMPLEX16 XLALCOMPLEX16Sqrt (COMPLEX16 z);
COMPLEX16 XLALCOMPLEX16SqrtReal (REAL8 x);

COMPLEX16 XLALCOMPLEX16Pow (COMPLEX16 a, COMPLEX16 b);
COMPLEX16 XLALCOMPLEX16PowReal (COMPLEX16 a, REAL8 x);

COMPLEX16 XLALCOMPLEX16Exp (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Log (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Log10 (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16LogB (COMPLEX16 a, COMPLEX16 b);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CSQRT(z) (XLALCOMPLEX16Sqrt(z))
#define LAL_CSQRT_REAL(x) (XLALCOMPLEX16SqrtReal(x))

#define LAL_CPOW(a,b) (XLALCOMPLEX16Pow((a),(b)))
#define LAL_CPOW_REAL(a,x) (XLALCOMPLEX16PowReal((a),(x)))

#define LAL_CEXP(a) (XLALCOMPLEX16Exp(a))
#define LAL_CLOG(a) (XLALCOMPLEX16Log(a))
#define LAL_CLOG10(a) (XLALCOMPLEX16Log10(a))
#define LAL_CLOGB(a,b) (XLALCOMPLEX16LogB((a),(b)))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define csqrtr(x) (XLALCOMPLEX16SqrtReal(x))

#define cpowr(a,x) (XLALCOMPLEX16PowReal((a),(x)))

#define clogb(a,b) (XLALCOMPLEX16LogB((a),(b)))
#endif
#endif


COMPLEX16 XLALCOMPLEX16Sin (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Cos (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Sec (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Csc (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Tan (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Cot (COMPLEX16 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CSIN(a) (XLALCOMPLEX16Sin(a))
#define LAL_CCOS(a) (XLALCOMPLEX16Cos(a))
#define LAL_CSEC(a) (XLALCOMPLEX16Sec(a))
#define LAL_CCSC(a) (XLALCOMPLEX16Csc(a))
#define LAL_CTAN(a) (XLALCOMPLEX16Tan(a))
#define LAL_CCOT(a) (XLALCOMPLEX16Cot(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define csec(a) (XLALCOMPLEX16Sec(a))
#define ccsc(a) (XLALCOMPLEX16Csc(a))
#define ccot(a) (XLALCOMPLEX16Cot(a))
#endif
#endif



COMPLEX16 XLALCOMPLEX16Arcsin (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArcsinReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccos (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArccosReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arcsec (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArcsecReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccsc (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArccscReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arctan (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arccot (COMPLEX16 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CASIN(a) (XLALCOMPLEX16Arcsin(a))
#define LAL_CASIN_REAL(x) (XLALCOMPLEX16ArcsinReal(x))
#define LAL_CACOS(a) (XLALCOMPLEX16Arccos(a))
#define LAL_CACOS_REAL(x) (XLALCOMPLEX16ArccosReal(x))
#define LAL_CASEC(a) (XLALCOMPLEX16Arcsec(a))
#define LAL_CASEC_REAL(x) (XLALCOMPLEX16ArcsecReal(x))
#define LAL_CACSC(a) (XLALCOMPLEX16Arccsc(a))
#define LAL_CACSC_REAL(x) (XLALCOMPLEX16ArccscReal(x))
#define LAL_CATAN(a) (XLALCOMPLEX16Arctan(a))
#define LAL_CACOT(a) (XLALCOMPLEX16Arccot(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define casinr(x) (XLALCOMPLEX16ArcsinReal(x))
#define cacosr(x) (XLALCOMPLEX16ArccosReal(x))
#define casec(a) (XLALCOMPLEX16Arcsec(a))
#define casecr(x) (XLALCOMPLEX16ArcsecReal(x))
#define cacsc(a) (XLALCOMPLEX16Arccsc(a))
#define cacscr(x) (XLALCOMPLEX16ArccscReal(x))
#define cacot(a) (XLALCOMPLEX16Arccot(a))
#endif
#endif



COMPLEX16 XLALCOMPLEX16Sinh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Cosh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Sech (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Csch (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Tanh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Coth (COMPLEX16 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CSINH(a) (XLALCOMPLEX16Sinh(a))
#define LAL_CCOSH(a) (XLALCOMPLEX16Cosh(a))
#define LAL_CSECH(a) (XLALCOMPLEX16Sech(a))
#define LAL_CCSCH(a) (XLALCOMPLEX16Csch(a))
#define LAL_CTANH(a) (XLALCOMPLEX16Tanh(a))
#define LAL_CCOTH(a) (XLALCOMPLEX16Coth(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define csech(a) (XLALCOMPLEX16Sech(a))
#define ccsch(a) (XLALCOMPLEX16Csch(a))
#define ccoth(a) (XLALCOMPLEX16Coth(a))
#endif
#endif



COMPLEX16 XLALCOMPLEX16Arcsinh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arccosh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArccoshReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arcsech (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arccsch (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16Arctanh (COMPLEX16 a);
COMPLEX16 XLALCOMPLEX16ArctanhReal (REAL8 x);
COMPLEX16 XLALCOMPLEX16Arccoth (COMPLEX16 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CASINH(a) (XLALCOMPLEX16Arcsinh(a))
#define LAL_CACOSH(a) (XLALCOMPLEX16Arccosh(a))
#define LAL_CACOSH_REAL(x) (XLALCOMPLEX16ArccoshReal(x))
#define LAL_CASECH(a) (XLALCOMPLEX16Arcsech(a))
#define LAL_CACSCH(a) (XLALCOMPLEX16Arccsch(a))
#define LAL_CATANH(a) (XLALCOMPLEX16Arctanh(a))
#define LAL_CATANH_REAL(x) (XLALCOMPLEX16ArctanhReal(x))
#define LAL_CACOTH(a) (XLALCOMPLEX16Arccoth(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define cacoshr(x) (XLALCOMPLEX16ArccoshReal(x))
#define casech(a) (XLALCOMPLEX16Arcsech(a))
#define cacsch(a) (XLALCOMPLEX16Arccsch(a))
#define catanhr(x) (XLALCOMPLEX16ArctanhReal(x))
#define cacoth(a) (XLALCOMPLEX16Arccoth(a))
#endif
#endif




#define LAL_COMPLEX8_ONE (XLALCOMPLEX8Rect(1.0,0.0))
#define LAL_COMPLEX8_ZERO (XLALCOMPLEX8Rect(0.0,0.0))
#define LAL_COMPLEX8_NEGONE (XLALCOMPLEX8Rect(-1.0,0.0))
#define LAL_COMPLEX8_I (XLALCOMPLEX8Rect(0.0,1.0))

#ifdef LAL_USE_COMPLEX_MACROS
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define crectf(x,y) (XLALCOMPLEX8Rect((x),(y)))
#define csetrf(x) (XLALCOMPLEX8Rect((x),0.0))
#define csetif(y) (XLALCOMPLEX8Rect(0.0,(y)))
#define cpolarf(x,y) (XLALCOMPLEX8Polar((x),(y)))

#define cunitf (XLALCOMPLEX8Rect(1.0,0.0))
#define czerof (XLALCOMPLEX8Rect(0.0,0.0))
#define cnegonef (XLALCOMPLEX8Rect(-1.0,0.0))
#define If (XLALCOMPLEX8Rect(0.0,1.0))
#endif
#endif

COMPLEX8 XLALCOMPLEX8Rect (REAL4 x, REAL4 y);
COMPLEX8 XLALCOMPLEX8Polar (REAL4 x, REAL4 y);

REAL4 XLALCOMPLEX8Arg (COMPLEX8 z);
REAL4 XLALCOMPLEX8Abs (COMPLEX8 z);
REAL4 XLALCOMPLEX8Abs2 (COMPLEX8 z);
REAL4 XLALCOMPLEX8LogAbs (COMPLEX8 z);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CARGF(z) (XLALCOMPLEX8Arg(z))
#define LAL_CABSF(z) (XLALCOMPLEX8Abs(z))
#define LAL_CABS2F(z) (XLALCOMPLEX8Abs2(z))
#define LAL_CLOGABSF(z) (XLALCOMPLEX8LogAbs(z))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define cabs2f(z) (XLALCOMPLEX8Abs2(z))
#define clogabsf(z) (XLALCOMPLEX8LogAbs(z))
#endif
#endif


COMPLEX8 XLALCOMPLEX8Add (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8Sub (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8Mul (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8Div (COMPLEX8 a, COMPLEX8 b);

COMPLEX8 XLALCOMPLEX8AddReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8SubReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8MulReal (COMPLEX8 a, REAL4 x);
COMPLEX8 XLALCOMPLEX8DivReal (COMPLEX8 a, REAL4 x);

COMPLEX8 XLALCOMPLEX8AddImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8SubImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8MulImag (COMPLEX8 a, REAL4 y);
COMPLEX8 XLALCOMPLEX8DivImag (COMPLEX8 a, REAL4 y);

COMPLEX8 XLALCOMPLEX8Conjugate (COMPLEX8 z);
COMPLEX8 XLALCOMPLEX8Inverse (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Negative (COMPLEX8 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CADDF(a,b) (XLALCOMPLEX8Add((a),(b)))
#define LAL_CSUBF(a,b) (XLALCOMPLEX8Sub((a),(b)))
#define LAL_CMULF(a,b) (XLALCOMPLEX8Mul((a),(b)))
#define LAL_CDIVF(a,b) (XLALCOMPLEX8Div((a),(b)))

#define LAL_CADD_REALF(a,x) (XLALCOMPLEX8AddReal((a),(x)))
#define LAL_CSUB_REALF(a,x) (XLALCOMPLEX8SubReal((a),(x)))
#define LAL_CMUL_REALF(a,x) (XLALCOMPLEX8MulReal((a),(x)))
#define LAL_CDIV_REALF(a,x) (XLALCOMPLEX8DivReal((a),(x)))

#define LAL_CADD_IMAGF(a,y) (XLALCOMPLEX8AddImag((a),(y)))
#define LAL_CSUB_IMAGF(a,y) (XLALCOMPLEX8SubImag((a),(y)))
#define LAL_CMUL_IMAGF(a,y) (XLALCOMPLEX8MulImag((a),(y)))
#define LAL_CDIV_IMAGF(a,y) (XLALCOMPLEX8DivImag((a),(y)))

#define LAL_CONJF(z) (XLALCOMPLEX8Conjugate(z))
#define LAL_CINVF(a) (XLALCOMPLEX8Inverse(a))
#define LAL_CNEGF(a) (XLALCOMPLEX8Negative(a))

#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define caddf(a,b) (XLALCOMPLEX8Add((a),(b)))
#define csubf(a,b) (XLALCOMPLEX8Sub((a),(b)))
#define cmulf(a,b) (XLALCOMPLEX8Mul((a),(b)))
#define cdivf(a,b) (XLALCOMPLEX8Div((a),(b)))

#define caddrf(a,x) (XLALCOMPLEX8AddReal((a),(x)))
#define csubrf(a,x) (XLALCOMPLEX8SubReal((a),(x)))
#define cmulrf(a,x) (XLALCOMPLEX8MulReal((a),(x)))
#define cdivrf(a,x) (XLALCOMPLEX8DivReal((a),(x)))

#define caddif(a,y) (XLALCOMPLEX8AddImag((a),(y)))
#define csubif(a,y) (XLALCOMPLEX8SubImag((a),(y)))
#define cmulif(a,y) (XLALCOMPLEX8MulImag((a),(y)))
#define cdivif(a,y) (XLALCOMPLEX8DivImag((a),(y)))

#define cinvf(a) (XLALCOMPLEX8Inverse(a))
#define cnegf(a) (XLALCOMPLEX8Negative(a))
#endif
#endif


COMPLEX8 XLALCOMPLEX8Sqrt (COMPLEX8 z);
COMPLEX8 XLALCOMPLEX8SqrtReal (REAL4 x);

COMPLEX8 XLALCOMPLEX8Pow (COMPLEX8 a, COMPLEX8 b);
COMPLEX8 XLALCOMPLEX8PowReal (COMPLEX8 a, REAL4 x);

COMPLEX8 XLALCOMPLEX8Exp (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Log (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Log10 (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8LogB (COMPLEX8 a, COMPLEX8 b);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CSQRTF(z) (XLALCOMPLEX8Sqrt(z))
#define LAL_CSQRT_REALF(x) (XLALCOMPLEX8SqrtReal(x))

#define LAL_CPOWF(a,b) (XLALCOMPLEX8Pow((a),(b)))
#define LAL_CPOW_REALF(a,x) (XLALCOMPLEX8PowReal((a),(x)))

#define LAL_CEXPF(a) (XLALCOMPLEX8Exp(a))
#define LAL_CLOGF(a) (XLALCOMPLEX8Log(a))
#define LAL_CLOG10F(a) (XLALCOMPLEX8Log10(a))
#define LAL_CLOGBF(a,b) (XLALCOMPLEX8LogB((a),(b)))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define csqrtrf(x) (XLALCOMPLEX8SqrtReal(x))

#define cpowrf(a,x) (XLALCOMPLEX8PowReal((a),(x)))

#define clogbf(a,b) (XLALCOMPLEX8LogB((a),(b)))
#endif
#endif


COMPLEX8 XLALCOMPLEX8Sin (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Cos (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Sec (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Csc (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Tan (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Cot (COMPLEX8 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CSINF(a) (XLALCOMPLEX8Sin(a))
#define LAL_CCOSF(a) (XLALCOMPLEX8Cos(a))
#define LAL_CSECF(a) (XLALCOMPLEX8Sec(a))
#define LAL_CCSCF(a) (XLALCOMPLEX8Csc(a))
#define LAL_CTANF(a) (XLALCOMPLEX8Tan(a))
#define LAL_CCOTF(a) (XLALCOMPLEX8Cot(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define csecf(a) (XLALCOMPLEX8Sec(a))
#define ccscf(a) (XLALCOMPLEX8Csc(a))
#define ccotf(a) (XLALCOMPLEX8Cot(a))
#endif
#endif



COMPLEX8 XLALCOMPLEX8Arcsin (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArcsinReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccos (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArccosReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arcsec (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArcsecReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccsc (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArccscReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arctan (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arccot (COMPLEX8 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CASINF(a) (XLALCOMPLEX8Arcsin(a))
#define LAL_CASIN_REALF(x) (XLALCOMPLEX8ArcsinReal(x))
#define LAL_CACOSF(a) (XLALCOMPLEX8Arccos(a))
#define LAL_CACOS_REALF(x) (XLALCOMPLEX8ArccosReal(x))
#define LAL_CASECF(a) (XLALCOMPLEX8Arcsec(a))
#define LAL_CASEC_REALF(x) (XLALCOMPLEX8ArcsecReal(x))
#define LAL_CACSCF(a) (XLALCOMPLEX8Arccsc(a))
#define LAL_CACSC_REALF(x) (XLALCOMPLEX8ArccscReal(x))
#define LAL_CATANF(a) (XLALCOMPLEX8Arctan(a))
#define LAL_CACOTF(a) (XLALCOMPLEX8Arccot(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define casinrf(x) (XLALCOMPLEX8ArcsinReal(x))
#define cacosrf(x) (XLALCOMPLEX8ArccosReal(x))
#define casecf(a) (XLALCOMPLEX8Arcsec(a))
#define casecrf(x) (XLALCOMPLEX8ArcsecReal(x))
#define cacscf(a) (XLALCOMPLEX8Arccsc(a))
#define cacscrf(x) (XLALCOMPLEX8ArccscReal(x))
#define cacotf(a) (XLALCOMPLEX8Arccot(a))
#endif
#endif



COMPLEX8 XLALCOMPLEX8Sinh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Cosh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Sech (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Csch (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Tanh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Coth (COMPLEX8 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CSINHF(a) (XLALCOMPLEX8Sinh(a))
#define LAL_CCOSHF(a) (XLALCOMPLEX8Cosh(a))
#define LAL_CSECHF(a) (XLALCOMPLEX8Sech(a))
#define LAL_CCSCHF(a) (XLALCOMPLEX8Csch(a))
#define LAL_CTANHF(a) (XLALCOMPLEX8Tanh(a))
#define LAL_CCOTHF(a) (XLALCOMPLEX8Coth(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define csechf(a) (XLALCOMPLEX8Sech(a))
#define ccschf(a) (XLALCOMPLEX8Csch(a))
#define ccothf(a) (XLALCOMPLEX8Coth(a))
#endif
#endif



COMPLEX8 XLALCOMPLEX8Arcsinh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arccosh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArccoshReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arcsech (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arccsch (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8Arctanh (COMPLEX8 a);
COMPLEX8 XLALCOMPLEX8ArctanhReal (REAL4 x);
COMPLEX8 XLALCOMPLEX8Arccoth (COMPLEX8 a);
#ifdef LAL_USE_COMPLEX_MACROS
#define LAL_CASINHF(a) (XLALCOMPLEX8Arcsinh(a))
#define LAL_CACOSHF(a) (XLALCOMPLEX8Arccosh(a))
#define LAL_CACOSH_REALF(x) (XLALCOMPLEX8ArccoshReal(x))
#define LAL_CASECHF(a) (XLALCOMPLEX8Arcsech(a))
#define LAL_CACSCHF(a) (XLALCOMPLEX8Arccsch(a))
#define LAL_CATANHF(a) (XLALCOMPLEX8Arctanh(a))
#define LAL_CATANH_REALF(x) (XLALCOMPLEX8ArctanhReal(x))
#define LAL_CACOTHF(a) (XLALCOMPLEX8Arccoth(a))
#ifdef LAL_USE_COMPLEX_SHORT_MACROS
#define cacoshrf(x) (XLALCOMPLEX8ArccoshReal(x))
#define casechf(a) (XLALCOMPLEX8Arcsech(a))
#define cacschf(a) (XLALCOMPLEX8Arccsch(a))
#define catanhrf(x) (XLALCOMPLEX8ArctanhReal(x))
#define cacothf(a) (XLALCOMPLEX8Arccoth(a))
#endif
#endif

#ifdef __cplusplus
}
#endif

#endif /* LAL_USE_OLD_COMPLEX_STRUCTS */

#endif /* _LALCOMPLEX8_H */
