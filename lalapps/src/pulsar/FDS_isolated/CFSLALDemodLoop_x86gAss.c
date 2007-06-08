/*
*  Copyright (C) 2007 Bernd Machenschalk
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

          {
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#ifdef __APPLE__ /* specials for Apples assembler */

#define AD_FLOAT   ".single "
#define AD_ASCII   ".ascii "
#define AD_ALIGN16 ".align 4"
#define AD_ALIGN64 ".align 6"
#else /* x86 gas */
#define AD_FLOAT   ".float "
#define AD_ASCII   ".string "
#define AD_ALIGN16 ".align 16"
#define AD_ALIGN64 ".align 64"
#endif
	    
            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;

	    /* prepare values for SSE */
	    REAL4 tempFreqX = tempFreq1; /* REAL4 because of SSE */
	    COMPLEX8 *Xalpha_kX = Xalpha_k; /* -> SSE values */

	    REAL4 XRes, XIms;   /* sums of Xa.re and Xa.im */
	    COMPLEX8 XSumsX;    /* sums of Xa.re and Xa.im for SSE */

	    /* shift values for FPU calculation */
	    COMPLEX8 *Xalpha_kF = Xalpha_k + 16;
	    REAL8    tempFreqF  = tempFreq1 - 15.0;

	    /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
	    /* one SSE register will consist 4 REAL4 values */
	    /* 4 REAL4 vaules = 2 ReIm pairs */

	    /* we want to do this SSE code 14 times */
#define ADD4SSE(a,b) \
	       "movlps "#a"(%[Xalpha_kX]),%%xmm4 \n\t"\
	       "movhps "#b"(%[Xalpha_kX]),%%xmm4 \n\t"\
	       "rcpps %%xmm0,%%xmm1       \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */\
	       "subps %%xmm5,%%xmm0       \n\t" /* XMM2: f-3 f-3 f-2 f-2 */\
	       "mulps %%xmm4,%%xmm1       \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
	       "addps %%xmm1,%%xmm2       \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */

	    __asm __volatile
	      (
	       /* constants */
	       "jmp cntcode              \n"
	       AD_ALIGN16 "              \n\t"
	       "C_V0011:                 \n\t"
	       AD_FLOAT "0.0             \n\t"
	       AD_FLOAT "0.0             \n\t"
	       "C_F1:                    \n\t"
	       AD_FLOAT "1.0             \n\t"
	       AD_FLOAT "1.0             \n"
	       "C_F2:                    \n\t"
	       "C_V2222:                 \n\t"
	       AD_FLOAT "2.0             \n\t"
	       AD_FLOAT "2.0             \n\t"
	       AD_FLOAT "2.0             \n\t"
	       AD_FLOAT "2.0             \n"

	       AD_ASCII "\"$Id$\"\n"
	       AD_ALIGN16 "              \n\t"
	       "cntcode:                 \n\t"

	       /* SSE prelude */
	       "movss %[tempFreqX],%%xmm0\n\t"
	       "shufps $0,%%xmm0,%%xmm0  \n\t" /* XMM0: f   f   f   f */
	       "subps  C_V0011,%%xmm0    \n\t" /* XMM0: f-1 f-1 f   f */
	       "xorps  %%xmm2,%%xmm2     \n\t" /* XMM2 will collect the low-precision values */
	       "movups C_V2222,%%xmm5    \n\t"
	       /* add two complex elements at a time */
	       ADD4SSE(0,8)
	       ADD4SSE(16,24)
	       ADD4SSE(32,40)

	       "FLDL %[tempFreq1]        \n\t" /* FLD D [TempFreqMinus15]   ;A15 */
	       "FLDS -16(%[Xalpha_k])    \n\t" /* FLD D [EBP-10h]   ;X14 A15 */
	       "FMUL %%ST(1),%%ST        \n\t" /* FMUL ST,ST1       ;X14A15 A15 */

	       ADD4SSE(48,56)

	       "FLDS -12(%[Xalpha_k])    \n\t" /* FLD D [EBP-0Ch]   ;Y14 X14A15 A15 */
	       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;Y14A15 X14A15 A15 */
	       "FXCH %%ST(2)             \n\t" /* FXCH ST2          ;A15 X14A15 Y14A15 */
	       "FLDS C_F1                \n\t" /* FLD D [_ONE_]     ;1 A15 X14A15 Y14A15 */
	       "FADD %%ST(1),%%ST        \n\t" /* FADD ST,ST1       ;A14 A15 X14A15 Y14A15 */
	       
	       ADD4SSE(64,72)

	       "FLDS -8(%[Xalpha_k])     \n\t" /* FLD D [EBP-08h]   ;X15 A14 A15 X14A15 Y14A15 */
	       "FMUL %%ST(1),%%ST        \n\t" /* FMUL ST,ST1       ;X15A14 A14 A15 X14A15 Y14A15 */
	       "FADDP %%ST,%%ST(3)       \n\t" /* FADDP ST3,ST      ;A14 A15 X' Y14A15 */
	       "FMUL %%ST,%%ST(1)        \n\t" /* FMUL ST1,ST       ;A14 Q145 X' Y14A15 */
	       "FLDS -4(%[Xalpha_k])     \n\t" /* FLD D [EBP-04h]   ;Y15 A14 Q145 X' Y14A15 */

	       ADD4SSE(80,88)

	       "FMUL %%ST(1),%%ST        \n\t" /* FMUL ST,ST1       ;Y15A14 A14 Q145 X' Y14A15 */
	       "FADDP %%ST,%%ST(4)       \n\t" /* FADDP ST4,ST      ;A14 Q145 X' Y' */
	       "FSUBS C_F2               \n\t" /* FSUB D [_TWO_]    ;A16 Q145 X' Y' */
	       "FMUL %%ST,%%ST(2)        \n\t" /* FMUL ST2,ST       ;A16 Q145 X'A16 Y' */
	       "FMUL %%ST,%%ST(3)        \n\t" /* FMUL ST3,ST       ;A16 Q145 X'A16 Y'A16 */
	       
	       ADD4SSE(96,104) 

	       /* SSE: skip FPU calculated values */
	       "subps %%xmm5,%%xmm0      \n\t"
	       "addl  $144,%[Xalpha_kX]  \n\t" /* Xalpha_kX = Xalpha_kX + 4; */
	       "subps %%xmm5,%%xmm0      \n\t"

	       "FLDS 0(%[Xalpha_k])      \n\t" /* FLD D [EBP+00h]   ;X16 A16 Q145 X'A16 Y'A16 */
	       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;X16Q145 A16 Q145 X'A16 Y'A16 */
	       "FADDP %%ST,%%ST(3)       \n\t" /* FADDP ST3,ST      ;A16 Q145 X" Y'A16 */
	       "FLDS 4(%[Xalpha_k])      \n\t" /* FLD D [EBP+04h]   ;Y16 A16 Q145 X" Y'A16 */
	       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;Y16Q145 A16 Q145 X" Y'A16 */

	       ADD4SSE(0,8)

	       "FADDP %%ST,%%ST(4)       \n\t" /* FADDP ST4,ST      ;A16 Q145 X" Y" */
	       "FMUL %%ST,%%ST(1)        \n\t" /* FMUL ST1,ST       ;A16 Q146 X" Y" */
	       "FSUBS C_F1               \n\t" /* FSUB D [_ONE_]    ;A17 Q146 X" Y" */
	       "FMUL %%ST,%%ST(2)        \n\t" /* FMUL ST2,ST       ;A17 Q146 X"A17 Y" */
	       "FMUL %%ST,%%ST(3)        \n\t" /* FMUL ST3,ST       ;A17 Q146 X"A17 Y"A17 */
	       
	       ADD4SSE(16,24)

	       "FLDS 8(%[Xalpha_k])      \n\t" /* FLD D [EBP+08h]   ;X17 A17 Q146 X"A17 Y"A17 */
	       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;X17Q146 A17 Q146 X"A17 Y"A17 */
	       "FADDP %%ST,%%ST(3)       \n\t" /* FADDP ST3,ST      ;A17 Q146 X! Y"A17 */

	       ADD4SSE(32,40)

	       "FLDS 12(%[Xalpha_k])     \n\t" /* FLD D [EBP+0Ch]   ;Y17 A17 Q146 X! Y"A17 */
	       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;Y17Q146 A17 Q146 X! Y"A17 */
	       "FADDP %%ST,%%ST(4)       \n\t" /* FADDP ST4,ST      ;A17 Q146 X! Y! */
	       "FMULP %%ST,%%ST(1)       \n\t" /* FMULP ST1,ST      ;Q147 X! Y! */
	       
	       ADD4SSE(48,56)

	       "FDIVRS C_F1              \n\t" /* FDIVR D [_ONE_]   ;1/Q x y */
	       "FMUL %%ST,%%ST(1)        \n\t" /* FMUL ST1,ST       ;1/Q xq y */
	       "FMULP %%ST,%%ST(2)       \n\t" /* FMULP ST2,ST      ;xq yq */
	       
	       ADD4SSE(64,72)
	       ADD4SSE(80,88)

	       "FSTPS %[XRes]            \n\t"
	       "FSTPS %[XIms]            \n\t"

	       ADD4SSE(96,104) 

	       /* add two complex elements at a time */

	       /* add the two calculated parts of the real and imaginary part */
	       "movhlps %%xmm2,%%xmm3    \n\t" /* XMM3: ? ? C_ImH C_ReH */
	       "addps %%xmm3,%%xmm2      \n\t" /* XMM2: - - C_Im C_Re */

	       /* store the result */
	       "movlps %%xmm2, %[XSumsX] \n\t"

	       : /* output  (here: to memory)*/
	       [XSumsX]  "=m" (XSumsX),
	       [XRes]  "=m" (XRes),
	       [XIms]  "=m" (XIms),
 	         /* input */
	       [Xalpha_kX] "+r" (Xalpha_kX) /* is changed by the code, so put into output section */
	       :
	       [Xalpha_k]  "r" (Xalpha_kF),
	       [tempFreq1] "m" (tempFreqF),
	       [tempFreqX] "m" (tempFreqX)

	       : /* clobbered registers */
#ifndef IGNORE_XMM_REGISTERS
	       "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5",
#endif
	       "st","st(1)","st(2)","st(3)","st(4)"
	       );	    

	    /* And last, we add the single and double precision values */
	    XRes = XRes + XSumsX.re;
	    XIms = XIms + XSumsX.im;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          }
