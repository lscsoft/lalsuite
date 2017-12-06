//
// Copyright (C) 2005, 2006  Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
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

{
  /* old 'vanilla' (pre-Akos) LALDemod hotloop algorithm, unrestricted
   * Dterms: based on version 5b0343e65a5a820d3e21a2afd9ba72123b05309c of
   * XLALComputeFaFb().
   */

  /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ] = sin [ 2pi kappa_star ],
   * therefore the trig-functions need to be calculated only once!
   * We choose the value sin[ 2pi kappa_star ] because it is the
   * closest to zero and will pose no numerical difficulties !
   * As kappa in [0, 1) we can skip the trimming step.
   */
  REAL4 s_alpha, c_alpha;   /* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
  XLALSinCos2PiLUTtrimmed ( &s_alpha, &c_alpha, kappa_star);
  c_alpha -= 1.0f;

  REAL8 kappa_max = kappa_star + 1.0f * Dterms - 1.0f;
  REAL8 x0 = kappa_max;

  realXP = 0;
  imagXP = 0;
  /* count down 2*Dterms values */
  for ( UINT4 l = 2 * Dterms; l > 0; l -- )
    {
      REAL4 realP, imagP;	/* real and imaginary parts of Dirichlet-kernel P_alpha_k */
      COMPLEX8 Xa = *Xalpha_l;
      REAL8 xinv = 1.0 / x0;

      /* calculate P_alpha_k */
      realP = s_alpha * xinv;
      imagP = c_alpha * xinv;

      /* calculate P_alpha_k * X_alpha_k */
      realXP += realP * crealf(Xa) - imagP * cimagf(Xa);
      imagXP += imagP * crealf(Xa) + realP * cimagf(Xa);

      Xalpha_l ++;	/* point to next frequency-bin */
      x0 -= 1.0 ;	/* x0-value for next iteration */

    } /* for k=kstar-Dterms to kstar+Dterms */

  /* real- and imaginary part of e^{i 2 pi lambda_alpha } */
  XLALSinCos2PiLUT ( &imagQ, &realQ, lambda_alpha );
}
