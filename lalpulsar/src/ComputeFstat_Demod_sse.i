//
// Copyright (C) 2007, 2008, 2009, 2010, 2012 Bernd Machenschalk, Reinhard Prix, Fekete Akos
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

/** SSE version from Akos */
{
  {

    COMPLEX8 XSums __attribute__ ((aligned (16))); /* sums of Xa.re and Xa.im for SSE */
    REAL4 kappa_m = kappa_max; /* single precision version of kappa_max */

    static REAL4 *scd               =  &(sincosLUTdiff[0]);
    static REAL4 *scb               =  &(sincosLUTbase[0]);
    static REAL4 M1 = -1.0f;
    static REAL8 sincos_adds = 402653184.0;
    REAL8 tmp;
    REAL8 _lambda_alpha = lambda_alpha;
    /* vector constants */
    /* having these not aligned will crash the assembler code */
    static REAL4 V0011[4] __attribute__ ((aligned (16))) = { 0,0,1,1 };
    static REAL4 V2222[4] __attribute__ ((aligned (16))) = { 2,2,2,2 };

    /* hand-coded SSE version from Akos */

    /* one loop iteration as a macro */

  #define VEC_LOOP_AV(a)\
       "MOVUPS " #a "(%[Xa]),%%xmm3       \n\t" \
       "SUBPS   %%xmm4,%%xmm2           \n\t" \
       "MULPS   %%xmm0,%%xmm3           \n\t" \
       "MULPS   %%xmm2,%%xmm1           \n\t" \
       "MULPS   %%xmm2,%%xmm0           \n\t" \
       "ADDPS   %%xmm3,%%xmm1           \n\t"


  #ifdef EAH_HOTLOOP_INTERLEAVED

  #define LIN_SIN_COS_TRIM_P0A(alpha) \
          "fldl %[" #alpha "] \n\t" \
          "fistpll %[tmp] \n\t" \
          "fld1 \n\t"     \
          "fildll %[tmp] \n\t"

  #define LIN_SIN_COS_TRIM_P0B(alpha)\
          "fsubrp %%st,%%st(1) \n\t" \
          "faddl %[" #alpha "] \n\t" \
          "faddl  %[sincos_adds]  \n\t" \
          "fstpl  %[tmp]    \n\t"

  #define LIN_SIN_COS_P0(alpha) \
          "fldl %[" #alpha "] \n\t" \
          "faddl  %[sincos_adds]  \n\t" \
          "fstpl  %[tmp]    \n\t"

  #define LIN_SIN_COS_P1 \
          "mov  %[tmp],%%eax \n\t"  \
          "mov  %%eax,%%edx  \n\t" \
          "and  $0x3fff,%%eax \n\t"
  #define LIN_SIN_COS_P2 \
          "mov  %%eax,%[tmp] \n\t"   \
          "mov  %[scd], %%eax \n\t"\
          "and  $0xffffff,%%edx \n\t" \
          "fildl %[tmp]\n\t"
  #define LIN_SIN_COS_P3 \
          "sar $0xe,%%edx \n\t" \
          "fld %%st  \n\t"   \
          "fmuls (%%eax,%%edx,4)   \n\t" \
          "mov  %[scb], %%edi \n\t"
  #define LIN_SIN_COS_P4(sin)\
          "fadds (%%edi,%%edx,4)   \n\t" \
          "add $0x100,%%edx \n\t"   \
          "fstps %[" #sin "] \n\t" \
          "fmuls (%%eax,%%edx,4)   \n\t"
  #define LIN_SIN_COS_P5(cos) \
          "fadds (%%edi,%%edx,4)   \n\t" \
          "fstps %[" #cos "] \n\t"


  #else
  #define LIN_SIN_COS_TRIM_P0A(alpha) ""
  #define LIN_SIN_COS_TRIM_P0B(alpha) ""
  #define LIN_SIN_COS_P0(alpha) ""
  #define LIN_SIN_COS_P1 ""
  #define LIN_SIN_COS_P2 ""
  #define LIN_SIN_COS_P3 ""
  #define LIN_SIN_COS_P4(sin) ""
  #define LIN_SIN_COS_P5(cos) ""

    SINCOS_2PI_TRIMMED ( &s_alpha, &c_alpha, kappa_star );

    SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
  #endif
    __asm __volatile
      (
       /* -------------------------------------------------------------------; */
       /* Prepare common divisor method for 4 values ( two Re-Im pair ) */
       /*  Im1, Re1, Im0, Re0 */
       "MOVAPS  %[V0011],%%xmm6         \n\t"
       "MOVSS   %[kappa_m],%%xmm2       \n\t"   /* X2:  -   -   -   C */
       "MOVUPS  (%[Xa]),%%xmm1          \n\t"   /* X1: Y01 X01 Y00 X00 */
       "SHUFPS  $0,%%xmm2,%%xmm2        \n\t"   /* X2:  C   C   C   C */
       "MOVAPS  %[V2222],%%xmm4         \n\t"   /* X7:  2   2   2   2 */
       "SUBPS   %%xmm6,%%xmm2           \n\t"   /* X2: C-1 C-1  C   C */
       /* -------------------------------------------------------------------; */
       "MOVAPS  %%xmm2,%%xmm0           \n\t"   /* X0: C-1 C-1  C   C */
       /* -------------------------------------------------------------------; */
       /* xmm0: collected denumerators -> a new element will multiply by this */
       /* xmm1: collected numerators -> we will divide it by the denumerator last */
       /* xmm2: current denumerator ( counter type ) */
       /* xmm3: current numerator ( current Re,Im elements ) */
       /* -------------------------------------------------------------------; */

       /* seven "loop iterations" (unrolled) */
  LIN_SIN_COS_P0(kappa_star)
       VEC_LOOP_AV(16)
  LIN_SIN_COS_P1
       VEC_LOOP_AV(32)
  LIN_SIN_COS_P2
       VEC_LOOP_AV(48)
  LIN_SIN_COS_P3
       VEC_LOOP_AV(64)
  LIN_SIN_COS_P4(sin)
       VEC_LOOP_AV(80)
  LIN_SIN_COS_P5(cos)
       VEC_LOOP_AV(96)

  LIN_SIN_COS_TRIM_P0A(_lambda_alpha)


  #if EAH_HOTLOOP_DIVS == EAH_HOTLOOP_DIVS_RECIPROCAL
       "movhlps %%xmm6,%%xmm6     \n\t"
  #endif
       "movss %[M1] , %%xmm5 \n\t"
       VEC_LOOP_AV(112)


  #if EAH_HOTLOOP_DIVS == EAH_HOTLOOP_DIVS_RECIPROCAL
       "movhlps %%xmm0,%%xmm2     \n\t"
       "mulss %%xmm0,%%xmm2     \n\t"

  #ifdef EAH_HOTLOOP_RENR
       "RCPSS %%xmm2,%%xmm6     \n\t"
       "MULSS %%xmm6,%%xmm2     \n\t"
       "MULSS %%xmm6,%%xmm2     \n\t"
  LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
       "ADDSS %%xmm6,%%xmm6     \n\t"
       "SUBSS %%xmm2,%%xmm6     \n\t"
  #else
       "divss %%xmm2,%%xmm6     \n\t"
  LIN_SIN_COS_TRIM_P0B(_lambda_alpha)

  #endif

       "shufps $78,%%xmm0,%%xmm0  \n\t"
       "mulps  %%xmm1,%%xmm0      \n\t"
       "movhlps %%xmm0,%%xmm4   \n\t"
       "addps %%xmm0,%%xmm4     \n\t"

       "shufps $160,%%xmm6,%%xmm6   \n\t"
       "mulps %%xmm6,%%xmm4     \n\t"
  #else
       /* -------------------------------------------------------------------; */
       /* Four divisions at once ( two for real parts and two for imaginary parts ) */
  #ifdef EAH_HOTLOOP_RENR
       "RCPPS %%xmm0,%%xmm6     \n\t"
       "MULPS %%xmm6,%%xmm0     \n\t"
       "MULPS %%xmm6,%%xmm0     \n\t"
  LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
       "ADDPS %%xmm6,%%xmm6     \n\t"
       "SUBPS %%xmm0,%%xmm6     \n\t"
       "MULPS %%xmm6,%%xmm1     \n\t"
  #else
       "DIVPS   %%xmm0,%%xmm1           \n\t"   /* X1: Y0G X0G Y1F X1F */
  LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
  #endif
       /* -------------------------------------------------------------------; */
       /* So we have to add the two real and two imaginary parts */
       "MOVHLPS   %%xmm1,%%xmm4           \n\t" /* X4:  -   -  Y0G X0G */
       "ADDPS   %%xmm1,%%xmm4           \n\t"   /* X4:  -   -  YOK XOK */
       /* -------------------------------------------------------------------; */
  #endif

  /*
  c_alpha-=1.0f;
    realXP = s_alpha * XSums.re - c_alpha * XSums.im;
    imagXP = c_alpha * XSums.re + s_alpha * XSums.im;
  */

          "movaps %%xmm4,%%xmm3 \n\t"
          "shufps $1,%%xmm3,%%xmm3 \n\t"
          "movss %[cos],%%xmm2 \n\t"
  LIN_SIN_COS_P1
          "movss %[sin],%%xmm1 \n\t"
          "addss %%xmm5,%%xmm2\n\t"
          "movss %%xmm2,%%xmm6 \n\t"
  LIN_SIN_COS_P2
          "movss %%xmm1,%%xmm5  \n\t"
          "mulss %%xmm4,%%xmm1 \n\t"
          "mulss %%xmm4,%%xmm2 \n\t"
  LIN_SIN_COS_P3
          "mulss %%xmm3,%%xmm5 \n\t"
          "mulss %%xmm3,%%xmm6 \n\t"
  LIN_SIN_COS_P4(Qimag)
          "addss %%xmm5,%%xmm2 \n\t"
          "subss %%xmm6,%%xmm1 \n\t"
  LIN_SIN_COS_P5(Qreal)
          "MOVss        %%xmm2,%[XPimag]        \n\t"   /*  */
          "MOVss        %%xmm1,%[XPreal]        \n\t"   /*  */

       /* interface */
       :
       /* output  (here: to memory)*/
       [XPreal]      "=m" (realXP),
       [XPimag]      "=m" (imagXP),
       [Qreal]      "=m" (realQ),
       [Qimag]      "=m" (imagQ),
       [tmp]        "=m" (tmp),
       [sin]        "=m" (s_alpha),
       [cos]        "=m" (c_alpha)

       :
       /* input */
       [Xa]          "r"  (Xalpha_l),
       [kappa_m]     "m"  (kappa_m),
       [kappa_star]  "m"  (kappa_star),
       [_lambda_alpha] "m" (_lambda_alpha),
       [scd]         "m"  (scd),
       [scb]         "m"  (scb),
       [sincos_adds]       "m"  (sincos_adds),
       [M1]         "m" (M1),


       /* vector constants */
       [V0011]       "m"  (V0011[0]),
       [V2222]       "m"  (V2222[0])

  #ifndef IGNORE_XMM_REGISTERS
       :
       /* clobbered registers */
       "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","st","st(1)","st(2)","eax","edx","edi","cc"
  #endif
       );

    /* moved the sin/cos call down here to avoid the store/forward stall of Core2s */

    /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
     * the trig-functions need to be calculated only once!
     * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
     * closest to zero and will pose no numerical difficulties !
     */

  }

  #ifndef EAH_HOTLOOP_INTERLEAVED
  {
    REAL8 _lambda_alpha = lambda_alpha;
    SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
    SINCOS_2PI_TRIMMED( &imagQ, &realQ, _lambda_alpha );
  }
  #endif
}
