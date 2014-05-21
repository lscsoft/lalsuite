//
// Copyright (C) 2007, 2008, 2009, 2010, 2012 Bernd Machenschalk, Reinhard Prix
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

#include "xmmintrin.h"

{
  {
    __declspec(align(16)) static struct { REAL4 a,b,c,d; } v0011 = {0.0, 0.0, 1.0, 1.0};
    __declspec(align(16)) static struct { REAL4 a,b,c,d; } v2222 = {2.0, 2.0, 2.0, 2.0};
    __declspec(align(16)) COMPLEX8 STn;

    REAL4 kappa_max = kappa_star + 1.0f * DTERMS - 1.0f;
    REAL4 kappa_m = kappa_max; /* single precision version of kappa_max */

    /* prelude */
    __asm {
        mov      esi , Xalpha_l                 /* Xal = Xalpha_l         */
        movss    xmm2, kappa_m                  /* pn[0] = kappa_max      */
        movlps   xmm1, MMWORD PTR [esi]         /* STnV = Xal ...         */
        movhps   xmm1, MMWORD PTR [esi+8]       /* ... continued          */
        shufps   xmm2, xmm2, 0                  /* pn[3]=pn[2]=pn[1]=pn[0]*/
        movaps   xmm4, XMMWORD PTR v2222        /* xmm4 = V2222           */
        subps    xmm2, XMMWORD PTR v0011        /* pn[2]-=1.0; pn[3]-=1.0 */
        movaps   xmm0, xmm2                     /* qn = pn                */
        };

    /* one loop iteration as a macro */
  #define VEC_LOOP_AV(a,b)\
    { \
        __asm movlps   xmm3, MMWORD PTR [esi+a] /* Xai = Xal[a]  ...*/\
        __asm movhps   xmm3, MMWORD PTR [esi+b] /* ... continued    */\
        __asm subps    xmm2, xmm4               /* pn   -= V2222    */\
        __asm mulps    xmm3, xmm0               /* Xai  *= qn       */\
        __asm mulps    xmm1, xmm2               /* STnV *= pn       */\
        __asm mulps    xmm0, xmm2               /* qn   *= pn       */\
        __asm addps    xmm1, xmm3               /* STnV += Xai      */\
        }

    /* seven macro calls i.e. loop iterations */
    VEC_LOOP_AV(16,24);
    VEC_LOOP_AV(32,40);
    VEC_LOOP_AV(48,56);
    VEC_LOOP_AV(64,72);
    VEC_LOOP_AV(80,88);
    VEC_LOOP_AV(96,104);
    VEC_LOOP_AV(112,120);

    /* four divisions and summing in SSE, then write out the result */
    __asm {
        divps    xmm1, xmm0                     /* STnV      /= qn       */
        movhlps  xmm4, xmm1                     /* / STnV[0] += STnV[2] \ */
        addps    xmm4, xmm1                     /* \ STnV[1] += STnV[3] / */
        movlps   STn, xmm4                      /* STn = STnV */
        };

    /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
     * the trig-functions need to be calculated only once!
     * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
     * closest to zero and will pose no numerical difficulties !
     */
    REAL4 s_alpha, c_alpha;   /* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
    SINCOS_2PI_TRIMMED ( &s_alpha, &c_alpha, kappa_star );
    c_alpha -= 1.0f;

    realXP = s_alpha * STn.re - c_alpha * STn.im;
    imagXP = c_alpha * STn.re + s_alpha * STn.im;

  }

  {
    REAL8 _lambda_alpha = lambda_alpha;
    SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
    SINCOS_2PI_TRIMMED( &imagQ, &realQ, _lambda_alpha );
  }
}
