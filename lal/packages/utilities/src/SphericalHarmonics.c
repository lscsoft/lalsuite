/*
 * Copyright (C) 2007 S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson,
 * C. Pankow
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

#include <lal/SphericalHarmonics.h>
#include <lal/LALError.h>
#include <lal/XLALGSL.h>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

/**
 * Computes the (s)Y(l,m) spin-weighted spherical harmonic.
 *
 * From somewhere ....
 *
 * See also:
 * Implements Equations (II.9)-(II.13) of
 * D. A. Brown, S. Fairhurst, B. Krishnan, R. A. Mercer, R. K. Kopparapu,
 * L. Santamaria, and J. T. Whelan,
 * "Data formats for numerical relativity waves",
 * arXiv:0709.0093v1 (2007).
 *
 * Currently only supports s=-2, l=2,3,4,5,6,7,8 modes.
 */
COMPLEX16 XLALSpinWeightedSphericalHarmonic(
                                   REAL8 theta,  /**< polar angle (rad) */
                                   REAL8 phi,    /**< azimuthal angle (rad) */
                                   int s,        /**< spin weight */
                                   int l,        /**< mode number l */
                                   int m         /**< mode number m */
    )
{
  REAL8 fac;
  COMPLEX16 ans;

  /* sanity checks ... */
  if ( l < abs(s) ) 
  {
    XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", __func__, s, l, m );
    XLAL_ERROR_VAL(0, XLAL_EINVAL);
  }
  if ( l < abs(m) ) 
  {
    XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
    XLAL_ERROR_VAL(0, XLAL_EINVAL);
  }

  if ( s == -2 ) 
  {
    if ( l == 2 ) 
    {
      switch ( m ) 
      {
        case -2:
          fac = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
          break;
        case -1:
          fac = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
          break;

        case 0:
          fac = sqrt( 15.0 / ( 32.0 * LAL_PI ) ) * sin( theta )*sin( theta );
          break;

        case 1:
          fac = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
          break;

        case 2:
          fac = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      } /*  switch (m) */
    }  /* l==2*/
    else if ( l == 3 ) 
    {
      switch ( m ) 
      {
        case -3:
          fac = sqrt(21.0/(2.0*LAL_PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
          break;
        case -2:
          fac = sqrt(7.0/(4.0*LAL_PI))*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
          break;
        case -1:
          fac = sqrt(35.0/(2.0*LAL_PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
          break;
        case 0:
          fac = (sqrt(105.0/(2.0*LAL_PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
          break;
        case 1:
          fac = -sqrt(35.0/(2.0*LAL_PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
          break;

        case 2:
          fac = sqrt(7.0/LAL_PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
          break;

        case 3:
          fac = -sqrt(21.0/(2.0*LAL_PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
          break;

        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      }
    }   /* l==3 */
    else if ( l == 4 ) 
    {
      switch ( m ) 
      {
        case -4:
          fac = 3.0*sqrt(7.0/LAL_PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
          break;
        case -3:
          fac = 3.0*sqrt(7.0/(2.0*LAL_PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
          break;

        case -2:
          fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(LAL_PI));
          break;
        case -1:
          fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*LAL_PI));
          break;
        case 0:
          fac = (3.0*sqrt(5.0/(2.0*LAL_PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
          break;
        case 1:
          fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*LAL_PI));
          break;
        case 2:
          fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(LAL_PI));
          break;
        case 3:
          fac = -3.0*sqrt(7.0/(2.0*LAL_PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
          break;
        case 4:
          fac = 3.0*sqrt(7.0/LAL_PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      }
    }    /* l==4 */
    else if ( l == 5 ) 
    {
      switch ( m ) 
      {
        case -5:
          fac = sqrt(330.0/LAL_PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
          break;
        case -4:
          fac = sqrt(33.0/LAL_PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
          break;
        case -3:
          fac = (sqrt(33.0/(2.0*LAL_PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
          break;
        case -2:
          fac = (sqrt(11.0/LAL_PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
          break;
        case -1:
          fac = (sqrt(77.0/LAL_PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
          break;
        case 0:
          fac = (sqrt(1155.0/(2.0*LAL_PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
          break;
        case 1:
          fac = sqrt(77.0/LAL_PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
          break;
        case 2:
          fac = sqrt(11.0/LAL_PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
          break;
        case 3:
          fac = -sqrt(33.0/(2.0*LAL_PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
          break;
        case 4:
          fac = sqrt(33.0/LAL_PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
          break;
        case 5:
          fac = -sqrt(330.0/LAL_PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      }
    }  /* l==5 */
    else if ( l == 6 )
    {
      switch ( m )
      {
        case -6:
          fac = (3.*sqrt(715./LAL_PI)*pow(cos(theta/2.0),4)*pow(sin(theta/2.0),8))/2.0;
          break;
        case -5:
          fac = (sqrt(2145./LAL_PI)*pow(cos(theta/2.0),3)*(1. + 3.*cos(theta))*pow(sin(theta/2.0),7))/2.0;
          break;
        case -4:
          fac = (sqrt(195./(2.0*LAL_PI))*pow(cos(theta/2.0),2)*(35. + 44.*cos(theta) 
          + 33.*cos(2.*theta))*pow(sin(theta/2.0),6))/8.0;
          break;
        case -3:
          fac = (3.*sqrt(13./LAL_PI)*cos(theta/2.0)*(98. + 185.*cos(theta) + 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*pow(sin(theta/2.0),5))/32.0;
          break;
        case -2:
          fac = (sqrt(13./LAL_PI)*(1709. + 3096.*cos(theta) + 2340.*cos(2.*theta) + 1320.*cos(3.*theta) 
          + 495.*cos(4.*theta))*pow(sin(theta/2.0),4))/256.0;
          break;
        case -1:
          fac = (sqrt(65./(2.0*LAL_PI))*cos(theta/2.0)*(161. + 252.*cos(theta) + 252.*cos(2.*theta) 
          + 132.*cos(3.*theta) + 99.*cos(4.*theta))*pow(sin(theta/2.0),3))/64.0;
          break;
        case 0:
          fac = (sqrt(1365./LAL_PI)*(35. + 60.*cos(2.*theta) + 33.*cos(4.*theta))*pow(sin(theta),2))/512.0;
          break;
        case 1:
          fac = (sqrt(65./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(161. - 252.*cos(theta) + 252.*cos(2.*theta) 
          - 132.*cos(3.*theta) + 99.*cos(4.*theta))*sin(theta/2.0))/64.0;
          break;
        case 2:
          fac = (sqrt(13./LAL_PI)*pow(cos(theta/2.0),4)*(1709. - 3096.*cos(theta) + 2340.*cos(2.*theta) 
          - 1320*cos(3*theta) + 495*cos(4*theta)))/256.0;
          break;
        case 3:
          fac = (-3.*sqrt(13./LAL_PI)*pow(cos(theta/2.0),5)*(-98. + 185.*cos(theta) - 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*sin(theta/2.0))/32.0;
          break;
        case 4:
          fac = (sqrt(195./(2.0*LAL_PI))*pow(cos(theta/2.0),6)*(35. - 44.*cos(theta) 
          + 33.*cos(2*theta))*pow(sin(theta/2.0),2))/8.0;
          break;
        case 5:
          fac = -(sqrt(2145./LAL_PI)*pow(cos(theta/2.0),7)*(-1. + 3.*cos(theta))*pow(sin(theta/2.0),3))/2.0;
          break;
        case 6:
          fac = (3.*sqrt(715./LAL_PI)*pow(cos(theta/2.0),8)*pow(sin(theta/2.0),4))/2.0;
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      }
    } /* l==6 */
    else if ( l == 7 )
    {
      switch ( m )
      {
        case -7:
          fac = sqrt(15015./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*pow(sin(theta/2.0),9);
          break;
        case -6:
          fac = (sqrt(2145./LAL_PI)*pow(cos(theta/2.0),4)*(2. + 7.*cos(theta))*pow(sin(theta/2.0),8))/2.0;
          break;
        case -5:
          fac = (sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(93. + 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),7))/8.0;
          break;
        case -4:
          fac = (sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),2)*(140. + 285.*cos(theta) 
          + 156.*cos(2.*theta) + 91.*cos(3.*theta))*pow(sin(theta/2.0),6))/16.0;
          break;
        case -3:
          fac = (sqrt(15./(2.0*LAL_PI))*cos(theta/2.0)*(3115. + 5456.*cos(theta) + 4268.*cos(2.*theta) 
          + 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*pow(sin(theta/2.0),5))/128.0;
          break;
        case -2:
          fac = (sqrt(15./LAL_PI)*(5220. + 9810.*cos(theta) + 7920.*cos(2.*theta) + 5445.*cos(3.*theta) 
          + 2860.*cos(4.*theta) + 1001.*cos(5.*theta))*pow(sin(theta/2.0),4))/512.0;
          break;
        case -1:
          fac = (3.*sqrt(5./(2.0*LAL_PI))*cos(theta/2.0)*(1890. + 4130.*cos(theta) + 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) + 1430.*cos(4.*theta) + 1001.*cos(5*theta))*pow(sin(theta/2.0),3))/512.0;
          break;
        case 0:
          fac = (3.*sqrt(35./LAL_PI)*cos(theta)*(109. + 132.*cos(2.*theta) 
          + 143.*cos(4.*theta))*pow(sin(theta),2))/512.0;
          break;
        case 1:
          fac = (3.*sqrt(5./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(-1890. + 4130.*cos(theta) - 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) - 1430.*cos(4.*theta) + 1001.*cos(5.*theta))*sin(theta/2.0))/512.0;
          break;
        case 2:
          fac = (sqrt(15./LAL_PI)*pow(cos(theta/2.0),4)*(-5220. + 9810.*cos(theta) - 7920.*cos(2.*theta) 
          + 5445.*cos(3.*theta) - 2860.*cos(4.*theta) + 1001.*cos(5.*theta)))/512.0;
          break;
        case 3:
          fac = -(sqrt(15./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(3115. - 5456.*cos(theta) + 4268.*cos(2.*theta) 
          - 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*sin(theta/2.0))/128.0;
          break;  
        case 4:
          fac = (sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),6)*(-140. + 285.*cos(theta) - 156.*cos(2*theta) 
          + 91.*cos(3.*theta))*pow(sin(theta/2.0),2))/16.0;
          break;
        case 5:
          fac = -(sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),7)*(93. - 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),3))/8.0;
          break;
        case 6:
          fac = (sqrt(2145./LAL_PI)*pow(cos(theta/2.0),8)*(-2. + 7.*cos(theta))*pow(sin(theta/2.0),4))/2.0;
          break;
        case 7:
          fac = -(sqrt(15015./(2.0*LAL_PI))*pow(cos(theta/2.0),9)*pow(sin(theta/2.0),5));
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      }
    } /* l==7 */
    else if ( l == 8 )
    {
      switch ( m )
      {
        case -8:
          fac = sqrt(34034./LAL_PI)*pow(cos(theta/2.0),6)*pow(sin(theta/2.0),10);
          break;
        case -7:
          fac = sqrt(17017./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(1. + 4.*cos(theta))*pow(sin(theta/2.0),9);
          break;
        case -6:
          fac = sqrt(255255./LAL_PI)*pow(cos(theta/2.0),4)*(1. + 2.*cos(theta))
          *sin(LAL_PI/4.0 - theta/2.0)*sin(LAL_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),8);
          break;
        case -5:
          fac = (sqrt(12155./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(19. + 42.*cos(theta) 
          + 21.*cos(2.*theta) + 14.*cos(3.*theta))*pow(sin(theta/2.0),7))/8.0;
          break;
        case -4:
          fac = (sqrt(935./(2.0*LAL_PI))*pow(cos(theta/2.0),2)*(265. + 442.*cos(theta) + 364.*cos(2.*theta) 
          + 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),6))/32.0;
          break;
        case -3:
          fac = (sqrt(561./(2.0*LAL_PI))*cos(theta/2.0)*(869. + 1660.*cos(theta) + 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) + 455.*cos(4.*theta) + 182.*cos(5.*theta))*pow(sin(theta/2.0),5))/128.0;
          break;
        case -2:
          fac = (sqrt(17./LAL_PI)*(7626. + 14454.*cos(theta) + 12375.*cos(2.*theta) + 9295.*cos(3.*theta) 
          + 6006.*cos(4.*theta) + 3003.*cos(5.*theta) + 1001.*cos(6.*theta))*pow(sin(theta/2.0),4))/512.0;
          break;
        case -1:
          fac = (sqrt(595./(2.0*LAL_PI))*cos(theta/2.0)*(798. + 1386.*cos(theta) + 1386.*cos(2.*theta) 
          + 1001.*cos(3.*theta) + 858.*cos(4.*theta) + 429.*cos(5.*theta) + 286.*cos(6.*theta))*pow(sin(theta/2.0),3))/512.0;
          break;
        case 0:
          fac = (3.*sqrt(595./LAL_PI)*(210. + 385.*cos(2.*theta) + 286.*cos(4.*theta) 
          + 143.*cos(6.*theta))*pow(sin(theta),2))/4096.0;
          break;
        case 1:
          fac = (sqrt(595./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(798. - 1386.*cos(theta) + 1386.*cos(2.*theta) 
          - 1001.*cos(3.*theta) + 858.*cos(4.*theta) - 429.*cos(5.*theta) + 286.*cos(6.*theta))*sin(theta/2.0))/512.0;
          break;
        case 2:
          fac = (sqrt(17./LAL_PI)*pow(cos(theta/2.0),4)*(7626. - 14454.*cos(theta) + 12375.*cos(2.*theta) 
          - 9295.*cos(3.*theta) + 6006.*cos(4.*theta) - 3003.*cos(5.*theta) + 1001.*cos(6.*theta)))/512.0;
          break;
        case 3:
          fac = -(sqrt(561./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(-869. + 1660.*cos(theta) - 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) - 455.*cos(4.*theta) + 182.*cos(5.*theta))*sin(theta/2.0))/128.0;
          break;
        case 4:
          fac = (sqrt(935./(2.0*LAL_PI))*pow(cos(theta/2.0),6)*(265. - 442.*cos(theta) + 364.*cos(2.*theta) 
          - 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),2))/32.0;
          break;
        case 5:
          fac = -(sqrt(12155./(2.0*LAL_PI))*pow(cos(theta/2.0),7)*(-19. + 42.*cos(theta) - 21.*cos(2.*theta) 
          + 14.*cos(3.*theta))*pow(sin(theta/2.0),3))/8.0;
          break;
        case 6:
          fac = sqrt(255255./LAL_PI)*pow(cos(theta/2.0),8)*(-1. + 2.*cos(theta))*sin(LAL_PI/4.0 - theta/2.0)
          *sin(LAL_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),4);
          break;
        case 7:
          fac = -(sqrt(17017./(2.0*LAL_PI))*pow(cos(theta/2.0),9)*(-1. + 4.*cos(theta))*pow(sin(theta/2.0),5));
          break;
        case 8:
          fac = sqrt(34034./LAL_PI)*pow(cos(theta/2.0),10)*pow(sin(theta/2.0),6);
          break;
        default:
          XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", __func__, s, l, m );
          XLAL_ERROR_VAL(0, XLAL_EINVAL);
          break;
      }
    } /* l==8 */
    else 
    {
      XLALPrintError("XLAL Error - %s: Unsupported mode l=%d (only l in [2,8] implemented)\n", __func__, l);
      XLAL_ERROR_VAL(0, XLAL_EINVAL);
    }
  }
  else 
  {
    XLALPrintError("XLAL Error - %s: Unsupported mode s=%d (only s=-2 implemented)\n", __func__, s);
    XLAL_ERROR_VAL(0, XLAL_EINVAL);
  }
  if (m)
    ans = cpolar(1.0, m*phi) * fac;
  else
    ans = fac;
  return ans;
}


/**
 * Computes the scalar spherical harmonic \f$ Y_{lm}(\theta, \phi) \f$.
 */
int
XLALScalarSphericalHarmonic(
                         COMPLEX16 *y, /**< output */
                         UINT4 l,      /**< value of l */
                         INT4  m,      /**< value of m */
                         REAL8 theta,  /**< angle theta */
                         REAL8 phi     /**< angle phi */
                         )
{

  int   gslStatus;
  gsl_sf_result pLm;

  INT4 absM = abs( m );

  if ( absM > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* For some reason GSL will not take negative m */
  /* We will have to use the relation between sph harmonics of +ve and -ve m */
  XLAL_CALLGSL( gslStatus = gsl_sf_legendre_sphPlm_e((INT4)l, absM, cos(theta), &pLm ) );
  if (gslStatus != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Compute the values for the spherical harmonic */
  *y = cpolar(pLm.val, m * phi);

  /* If m is negative, perform some jiggery-pokery */
  if ( m < 0 && absM % 2  == 1 )
  {
    *y = - *y;
  }

  return XLAL_SUCCESS;
}

/**
 * Computes the spin 2 weighted spherical harmonic. This function is now
 * deprecated and will be removed soon. All calls should be replaced with
 * calls to XLALSpinWeightedSphericalHarmonic().
 */
INT4 XLALSphHarm ( COMPLEX16 *out, /**< output */
                   UINT4   L,      /**< value of L */
                   INT4 M,         /**< value of M */
                   REAL4 theta,    /**< angle with respect to the z axis */
                   REAL4   phi     /**< angle with respect to the x axis */
                   )
{

  XLALPrintDeprecationWarning( "XLALSphHarm", "XLALSpinWeightedSphericalHarmonic" );

  *out = XLALSpinWeightedSphericalHarmonic( theta, phi, -2, L, M );
  if ( xlalErrno )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  return XLAL_SUCCESS;
}

/**
 * Computes the n-th Jacobi polynomial for polynomial weights alpha and beta.
 * The implementation here is only valid for real x -- enforced by the argument
 * type. An extension to complex values would require evaluation of several 
 * gamma functions.
 *
 * See http://en.wikipedia.org/wiki/Jacobi_polynomials
 */
double XLALJacobiPolynomial(int n, int alpha, int beta, double x){
	double f1 = (x-1)/2.0, f2 = (x+1)/2.0;
	int s=0;
	double sum=0, val=0;
	if( n == 0 ) return 1.0;
	for( s=0; n-s >= 0; s++ ){
		val=1.0;
		val *= gsl_sf_choose( n+alpha, s );
		val *= gsl_sf_choose( n+beta, n-s );
		if( n-s != 0 ) val *= pow( f1, n-s );
		if( s != 0 ) val*= pow( f2, s );

		sum += val;
	}
	return sum;
}

/**
 * Computes the 'little' d Wigner matrix for the Euler angle beta. Single angle 
 * small d transform with major index 'l' and minor index transition from m to 
 * mp. 
 *
 * Uses a slightly unconventional method since the intuitive version by Wigner 
 * is less suitable to algorthmic development. 
 *
 * See http://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_.28small.29_d-matrix
 */
#define MIN(a,b) ((a) < (b) ? (a) : (b))
double XLALWignerdMatrix(
                                   int l,        /**< mode number l */
                                   int mp,        /**< mode number m' */
                                   int m,        /**< mode number m */
                                   double beta  /**< euler angle (rad) */
    )
{

	int k = MIN( l+m, MIN( l-m, MIN( l+mp, l-mp )));
	double a=0, lam=0;
	if(k == l+m){
		a = mp-m;
		lam = mp-m;
	} else if(k == l-m) {
		a = m-mp;
		lam = 0;
	} else if(k == l+mp) {
		a = m-mp;
		lam = 0;
	} else if(k == l-mp) {
		a = mp-m;
		lam = mp-m;
	}

	int b = 2*l-2*k-a;
	double pref = pow(-1, lam) * sqrt(gsl_sf_choose( 2*l-k, k+a )) / sqrt(gsl_sf_choose( k+b, b ));

	return pref * pow(sin(beta/2.0), a) * pow( cos(beta/2.0), b) * XLALJacobiPolynomial(k, a, b, cos(beta));

}

/**
 * Computes the full Wigner D matrix for the Euler angle alpha, beta, and gamma
 * with major index 'l' and minor index transition from m to mp. 
 *
 * Uses a slightly unconventional method since the intuitive version by Wigner 
 * is less suitable to algorthmic development. 
 *
 * See http://en.wikipedia.org/wiki/Wigner_D-matrix
 *
 * Currently only supports the modes which are implemented for the spin
 * weighted spherical harmonics.
 */
COMPLEX16 XLALWignerDMatrix(
                                   int l,        /**< mode number l */
                                   int mp,        /**< mode number m' */
                                   int m,        /**< mode number m */
                                   double alpha,  /**< euler angle (rad) */
                                   double beta, /**< euler angle (rad) */
                                   double gam  /**< euler angle (rad) */
    )
{
	 return cexp( -(1.0I)*mp*alpha ) *
			XLALWignerdMatrix( l, mp, m, beta ) * 
			cexp( -(1.0I)*m*gam );
}
