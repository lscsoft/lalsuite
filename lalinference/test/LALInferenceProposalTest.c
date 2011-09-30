/* 
 *  LALInferenceProposalTest.c: Testing the jump propsals in LALInferenceProposal.c
 *
 *  Copyright (C) 2011 Will M. Farr
 *
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

#include<lal/LALInferenceProposal.h>
#include<lal/LALInference.h>
#include<lal/XLALError.h>
#include<math.h>

/** Efficient integer power computation. */
REAL8 pow_int(const REAL8, const INT4);
REAL8
pow_int(const REAL8 x, const INT4 n) {
  if (n < 0) {
    return 1.0/pow_int(x, -n);
  } else if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  } else if (n % 2 == 0) {
    REAL8 sqrt_x = pow_int(x, n/2);
    return sqrt_x * sqrt_x;
  } else {
    return x*pow_int(x, n-1);
  }
}

/** Cumulative distribution function for KS statistic.  Algorithm from
    Numerical Recipes, Third Edition by Press, Teukolsky, Vetterling
    and Flannery.  Cambridge University Press, 2007. Section
    6.14.12 */
REAL8 PKS(const REAL8);
REAL8
PKS(const REAL8 z) {
  if (z < 0.0) {
    XLALError("PKS", __FILE__, __LINE__, XLAL_FAILURE);
    return -1.0;
  } else if (z < 0.042) {
    return 0.0;
  } else if (z < 1.18) {
    REAL8 x = exp(-1.23370055013616983/(z*z));
    return 2.25675833419102515*sqrt(-log(x))*(x + pow_int(x,9) + pow_int(x, 25) + pow_int(x, 49));
  } else {
    REAL8 x = exp(-2.0*z*z);
    return 1.0 - 2.0*(x - pow_int(x,4) + pow_int(x,9));
  }
}

/** Compliment of PKS(). */
REAL8 QKS(const REAL8);
REAL8
QKS(const REAL8 z) {
  if (z < 0.0) {
    XLALError("QKS", __FILE__, __LINE__, XLAL_FAILURE);
    return -1.0;
  } else if (z == 0.0) {
    return 1.0;
  } else if (z < 1.18) {
    return 1.0 - PKS(z);
  } else {
    REAL8 x = exp(-2.0*z*z);
    return 2.0*(x - pow_int(x, 4) + pow_int(x, 9));
  }
}

/** Computes the p of the KS-statistic comparing the cumulative
    distribution given by the discrete points in \a points with the
    corresponding analytic cumulative distribution values in \a
    cumValues (assumed to be evaluated at the locations in points).
    The input array \a points must be sorted. */
REAL8 KSPValue(const REAL8Vector *, const REAL8Vector *);
REAL8
KSPValue(const REAL8Vector *points, const REAL8Vector *cumValues) {
  UINT4 i;
  REAL8 maxD = 0.0;

  if (points->length != cumValues->length) {
    XLALError("KSPValue", __FILE__, __LINE__, XLAL_FAILURE);
    return -1.0;
  }

  /* Check for sorted points. */
  for (i = 0; i < points->length-1; i++) {
    if (points->data[i+1] < points->data[i]) {
      XLALError("KSPValue", __FILE__, __LINE__, XLAL_FAILURE);
      return -1.0;
    }
  }

  maxD = 0.0;
  for (i = 0; i < points->length; i++) {
    REAL8 D = fabs((REAL8)i/((REAL8) points->length) - cumValues->data[i]);
    maxD = (D > maxD ? D : maxD);
  }

  return QKS((sqrt(points->length) + 0.12 + 0.11/sqrt(points->length))*maxD);
}

int main(void) {
  return 0;
}
