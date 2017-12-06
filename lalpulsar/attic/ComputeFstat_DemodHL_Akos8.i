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

{
  {
    /* improved hotloop algorithm by Fekete Akos:
     * take out repeated divisions into a single common denominator,
     * plus use extra cleverness to compute the nominator efficiently...
     */
    REAL4 kappa_max = kappa_star + 1.0f * DTERMS - 1.0f;
    REAL4 Sn = crealf(*Xalpha_l);
    REAL4 Tn = cimagf(*Xalpha_l);
    REAL4 pn = kappa_max;
    REAL4 qn = pn;
    REAL4 U_alpha, V_alpha;

    /* 2*DTERMS iterations */
    UINT4 l;
    for ( l = 1; l < 2*DTERMS; l ++ )
      {
        Xalpha_l ++;
        pn = pn - 1.0f;                         /* p_(n+1) */
        Sn = pn * Sn + qn * crealf(*Xalpha_l);  /* S_(n+1) */
        Tn = pn * Tn + qn * cimagf(*Xalpha_l);  /* T_(n+1) */
        qn *= pn;                               /* q_(n+1) */
      } /* for l < 2*DTERMS */

    { /* only one division left */
      REAL4 r_qn = 1.0 / qn;
      U_alpha = Sn * r_qn;
      V_alpha = Tn * r_qn;
    }

    /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
     * the trig-functions need to be calculated only once!
     * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
     * closest to zero and will pose no numerical difficulties !
     */
    REAL4 s_alpha, c_alpha;   /* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
    XLALSinCos2PiLUTtrimmed ( &s_alpha, &c_alpha, kappa_star );
    c_alpha -= 1.0f;

    realXP = s_alpha * U_alpha - c_alpha * V_alpha;
    imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
  }

  XLALSinCos2PiLUT ( &imagQ, &realQ, lambda_alpha );
}
