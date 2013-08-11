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

#include <altivec.h>

{
  {
    SINCOS_2PI_TRIMMED ( &s_alpha, &c_alpha, kappa_star );
    c_alpha -= 1.0f;
    {
      REAL4 *Xalpha_kR4 = (REAL4*)(Xalpha_l);
      REAL4 kappa_max = kappa_star + 1.0f * DTERMS - 1.0f;

      float STn[4] __attribute__ ((aligned (16)));       /* aligned for vector output */
      /* the vectors actually become registers in the AVUnit */
      vector unsigned char perm;   /* permutation pattern for unaligned memory access */
      vector float load0, load1, load2; /* temp registers for unaligned memory access */
      vector float XaiV   /* xmm3 */;                  /* SFT data loaded from memory */
      vector float STnV   /* xmm1 */;                         /* sums up the dividend */
      vector float V0000             = {0,0,0,0};             /* zero vector constant */
      vector float V2222  /* xmm4 */ = {2,2,2,2};                  /* vector constant */
      vector float pnV    /* xmm2 */ = {((float)(kappa_max)),
                                        ((float)(kappa_max)),
                                        ((float)(kappa_max - 1)),
                                        ((float)(kappa_max - 1)) };
      vector float qnV    /* xmm0 */ = pnV;   /* common divisor, initally = 1.0 * pnV */
      /*    this column above (^) lists the corresponding register in the SSE version */

      vector float tV;          /* temporary vector used for Newton-Rhapson iterarion */

      /* init the memory access (load0,load1) */
      load0   = vec_ld  (0,(Xalpha_kR4));
      perm    = vec_lvsl(0,(Xalpha_kR4));
      load1   = vec_ld  (0,(Xalpha_kR4+4));

      /* first "iteration" & initialization */
      XaiV    = vec_perm(load0,load1,perm);
      qnV     = vec_re(pnV);
      STnV    = vec_madd(XaiV, qnV, V0000);

      /* use a reciprocal estimate as a replacement for a division.
         in our case this is only valid for the "outer" elements of the kernel loop */
#define VEC_LOOP_RE(n,a,b)\
      pnV     = vec_sub(pnV,V2222);\
      perm    = vec_lvsl(0,(Xalpha_kR4+(n)));\
      load##b = vec_ld(0,(Xalpha_kR4+(n)+4));\
      XaiV    = vec_perm(load##a,load##b,perm);\
      qnV     = vec_re(pnV);\
      STnV    = vec_madd(XaiV, qnV, STnV);  /* STnV = XaiV * qnV + STnV */

      /* refine the reciprocal estimate to by a Newton-Rhapson iteration.
         re1(x) = re0(x) * (2 - x * re0(x))
         (see http://en.wikipedia.org/wiki/Division_(digital)#Newton-Raphson_division)
         this should give as much precision as a normal float division */
#define VEC_LOOP_RE_NR(n,a,b)\
      pnV     = vec_sub(pnV,V2222);\
      perm    = vec_lvsl(0,(Xalpha_kR4+(n)));\
      load##b = vec_ld(0,(Xalpha_kR4+(n)+4));\
      XaiV    = vec_perm(load##a,load##b,perm);\
      qnV     = vec_re(pnV);\
      tV      = vec_madd(qnV,pnV,V0000);\
      tV      = vec_sub(V2222,tV);\
      qnV     = vec_madd(qnV,tV,V0000);\
      STnV    = vec_madd(XaiV, qnV, STnV);

      /* actual "hot loop" (unrolled) */
      VEC_LOOP_RE     (4,1,2);
      VEC_LOOP_RE     (8,2,0);
      VEC_LOOP_RE_NR (12,0,1);
      VEC_LOOP_RE_NR (16,1,2);
      VEC_LOOP_RE_NR (20,2,0);
      VEC_LOOP_RE    (24,0,1);
      VEC_LOOP_RE    (28,1,0);

      /* output the vector */
      vec_st(STnV,0,STn);

      /* combine the sums */
      {
        REAL4 U_alpha = STn[0] + STn[2];
        REAL4 V_alpha = STn[1] + STn[3];

        realXP = s_alpha * U_alpha - c_alpha * V_alpha;
        imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
      }
    }
  }

  {
    REAL8 _lambda_alpha = lambda_alpha;
    SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
    SINCOS_2PI_TRIMMED( &imagQ, &realQ, _lambda_alpha );
  }
}
