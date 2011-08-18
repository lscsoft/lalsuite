/*
 * Copyright (C) 2007 S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson
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


/** 
 * \defgroup SphericalHarmonics Spin-weighted Spherical Harmonics
 * \ingroup support
 * \author S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson
 *
 * \brief Library of spherical harmonic functions
 *

 *
 */

#include <lal/SphericalHarmonics.h>
#include <lal/LALError.h>
#include <lal/XLALGSL.h>

#include <gsl/gsl_sf_legendre.h>

/** 
  * Spin 2 weighted spherical Harmonic
  */
INT4 XLALSpinWeightedSphHarm ( COMPLEX16 *out, /**< output */
                               UINT4   L,      /**< value of L */
                               INT4    M,      /**< value of M */
                               REAL4   theta,  /**< angle with respect to the z axis */
                               REAL4   phi     /**< angle with respect to the x axis */)

{
  REAL4      deptheta; /* dependency on theta */

  /* check arguments are sensible */
  if ( !out ) {
    LALPrintError ("\nOutput pointer is NULL !\n\n");
    XLAL_ERROR ( __func__, XLAL_EINVAL);
  }

  if (L == 2)
  {
    switch ( M ) {
      case -2:
        deptheta = (sqrt(5./LAL_PI)*pow(sin(theta/2.0),4))/2.0;
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = sqrt(5./LAL_PI)*cos(theta/2.0)*pow(sin(theta/2.0),3);
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (sqrt(15./(2.0*LAL_PI))*pow(sin(theta),2))/4.0;
        out->re = deptheta;
        out->im = 0.;
        break;

      case 1:
        deptheta = sqrt(5./LAL_PI)*pow(cos(theta/2.0),3)*sin(theta/2.0);
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (sqrt(5./LAL_PI)*pow(cos(theta/2.0),4))/2.0;
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 2 */

  else if (L == 3)
  {
    switch ( M ) {
      case -3:
        deptheta = sqrt(21./(2.0*LAL_PI))*cos(theta/2.0)*pow(sin(theta/2.0),5);
        out->re = deptheta * cos( -3.*phi );
        out->im = deptheta * sin( -3.*phi );
        break;

      case -2:
        deptheta = (sqrt(7./LAL_PI)*(2. + 3.*cos(theta))*pow(sin(theta/2.0),4))/2.0;
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = (sqrt(35./(2.0*LAL_PI))*cos(theta/2.0)*(1. + 3.*cos(theta))*pow(sin(theta/2.0),3))/2.0;
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (sqrt(105./(2.0*LAL_PI))*cos(theta)*pow(sin(theta),2))/4.0;
        out->re = deptheta;
        out->im = 0.;
        break;

      case 1:
        deptheta = (sqrt(35./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(-1. + 3.*cos(theta))*sin(theta/2.0))/2.0;
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (sqrt(7./LAL_PI)*pow(cos(theta/2.0),4)*(-2. + 3.*cos(theta)))/2.0;
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      case 3:
        deptheta = -(sqrt(21./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*sin(theta/2.0));
        out->re = deptheta * cos( 3.*phi );
        out->im = deptheta * sin( 3.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 3 */

  else if (L == 4)
  {
    switch ( M ) {
      case -4:
        deptheta = 3.*sqrt(7./LAL_PI)*pow(cos(theta/2.0),2)*pow(sin(theta/2.0),6);
        out->re = deptheta * cos( -4.*phi );
        out->im = deptheta * sin( -4.*phi );
        break;

      case -3:
        deptheta = 3.*sqrt(7./(2.0*LAL_PI))*cos(theta/2.0)*(1. + 2.*cos(theta))*pow(sin(theta/2.0),5);
        out->re = deptheta * cos( -3.*phi );
        out->im = deptheta * sin( -3.*phi );
        break;

      case -2:
        deptheta = (3.*(9. + 14.*cos(theta) + 7.*cos(2.*theta))*pow(sin(theta/2.0),4))/(4.0*sqrt(LAL_PI));
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = (3.*cos(theta/2.0)*(6. + 7.*cos(theta) + 7.*cos(2.*theta))*pow(sin(theta/2.0),3))/(2.0*sqrt(2.*LAL_PI));
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (3.*sqrt(5./(2.0*LAL_PI))*(5. + 7.*cos(2.*theta))*pow(sin(theta),2))/16.0;
        out->re = deptheta;
        out->im = 0.;
        break;

      case 1:
        deptheta = (3.*pow(cos(theta/2.0),3)*(6. - 7.*cos(theta) + 7.*cos(2.*theta))*sin(theta/2.0))/(2.0*sqrt(2.*LAL_PI));
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (3.*pow(cos(theta/2.0),4)*(9. - 14.*cos(theta) + 7.*cos(2.*theta)))/(4.0*sqrt(LAL_PI));
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      case 3:
        deptheta = -3.*sqrt(7./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(-1. + 2.*cos(theta))*sin(theta/2.0);
        out->re = deptheta * cos( 3.*phi );
        out->im = deptheta * sin( 3.*phi );
        break;

      case 4:
        deptheta = 3.*sqrt(7./LAL_PI)*pow(cos(theta/2.0),6)*pow(sin(theta/2.0),2);
        out->re = deptheta * cos( 4.*phi );
        out->im = deptheta * sin( 4.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 4 */

  else if (L == 5)
  {
    switch ( M ) {
      case -5:
        deptheta = sqrt(330./LAL_PI)*pow(cos(theta/2.0),3)*pow(sin(theta/2.0),7);
        out->re = deptheta * cos( -5.*phi );
        out->im = deptheta * sin( -5.*phi );
        break;

      case -4:
        deptheta = sqrt(33./LAL_PI)*pow(cos(theta/2.0),2)*(2. + 5.*cos(theta))*pow(sin(theta/2.0),6);
        out->re = deptheta * cos( -4.*phi );
        out->im = deptheta * sin( -4.*phi );
        break;

      case -3:
        deptheta = (sqrt(33./(2.0*LAL_PI))*cos(theta/2.0)*(17. + 24*cos(theta) + 15.*cos(2.*theta))*pow(sin(theta/2.0),5))/4.0;
        out->re = deptheta * cos( -3.*phi );
        out->im = deptheta * sin( -3.*phi );
        break;

      case -2:
        deptheta = (sqrt(11./LAL_PI)*(32. + 57.*cos(theta) + 36.*cos(2.*theta) + 15.*cos(3.*theta))*pow(sin(theta/2.0),4))/8.0;
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = (sqrt(77./LAL_PI)*cos(theta/2.0)*(14. + 33.*cos(theta) + 18.*cos(2.*theta) + 15.*cos(3.*theta))*pow(sin(theta/2.0),3))/16.0;
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (sqrt(1155./(2.0*LAL_PI))*cos(theta)*(1. + 3.*cos(2.*theta))*pow(sin(theta),2))/16.0;
        out->re = deptheta;
        out->im = 0;
        break;

      case 1:
        deptheta = (sqrt(77./LAL_PI)*pow(cos(theta/2.0),3)*(-14. + 33.*cos(theta) - 18.*cos(2*theta) + 15.*cos(3.*theta))*sin(theta/2.0))/16.0;
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (sqrt(11./LAL_PI)*pow(cos(theta/2.0),4)*(-32. + 57.*cos(theta) - 36.*cos(2.*theta) + 15.*cos(3.*theta)))/8.0;
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      case 3:
        deptheta = -(sqrt(33./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(17. - 24.*cos(theta) + 15.*cos(2.*theta))*sin(theta/2.0))/4.0;
        out->re = deptheta * cos( 3.*phi );
        out->im = deptheta * sin( 3.*phi );
        break;

      case 4:
        deptheta = sqrt(33./LAL_PI)*pow(cos(theta/2.0),6)*(-2. + 5.*cos(theta))*pow(sin(theta/2.0),2);
        out->re = deptheta * cos( 4.*phi );
        out->im = deptheta * sin( 4.*phi );
        break;

      case 5:
        deptheta = -(sqrt(330./LAL_PI)*pow(cos(theta/2.0),7)*pow(sin(theta/2.0),3));
        out->re = deptheta * cos( 5.*phi );
        out->im = deptheta * sin( 5.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 5 */

  else if (L == 6)
  {
    switch ( M ) {
      case -6:
        deptheta = (3.*sqrt(715./LAL_PI)*pow(cos(theta/2.0),4)*pow(sin(theta/2.0),8))/2.0;
        out->re = deptheta * cos( -6.*phi );
        out->im = deptheta * sin( -6.*phi );
        break;

      case -5:
        deptheta = (sqrt(2145./LAL_PI)*pow(cos(theta/2.0),3)*(1. + 3.*cos(theta))*pow(sin(theta/2.0),7))/2.0;
        out->re = deptheta * cos( -5.*phi );
        out->im = deptheta * sin( -5.*phi );
        break;

      case -4:
        deptheta = (sqrt(195./(2.0*LAL_PI))*pow(cos(theta/2.0),2)*(35. + 44.*cos(theta) 
          + 33.*cos(2.*theta))*pow(sin(theta/2.0),6))/8.0;
        out->re = deptheta * cos( -4.*phi );
        out->im = deptheta * sin( -4.*phi );
        break;

      case -3:
        deptheta = (3.*sqrt(13./LAL_PI)*cos(theta/2.0)*(98. + 185.*cos(theta) + 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*pow(sin(theta/2.0),5))/32.0;
        out->re = deptheta * cos( -3.*phi );
        out->im = deptheta * sin( -3.*phi );
        break;

      case -2:
        deptheta = (sqrt(13./LAL_PI)*(1709. + 3096.*cos(theta) + 2340.*cos(2.*theta) + 1320.*cos(3.*theta) 
          + 495.*cos(4.*theta))*pow(sin(theta/2.0),4))/256.0;
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = (sqrt(65./(2.0*LAL_PI))*cos(theta/2.0)*(161. + 252.*cos(theta) + 252.*cos(2.*theta) 
          + 132.*cos(3.*theta) + 99.*cos(4.*theta))*pow(sin(theta/2.0),3))/64.0;
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (sqrt(1365./LAL_PI)*(35. + 60.*cos(2.*theta) + 33.*cos(4.*theta))*pow(sin(theta),2))/512.0;
        out->re = deptheta;
        out->im = 0;
        break;

      case 1:
        deptheta = (sqrt(65./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(161. - 252.*cos(theta) + 252.*cos(2.*theta) 
          - 132.*cos(3.*theta) + 99.*cos(4.*theta))*sin(theta/2.0))/64.0;
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (sqrt(13./LAL_PI)*pow(cos(theta/2.0),4)*(1709. - 3096.*cos(theta) + 2340.*cos(2.*theta) 
          - 1320*cos(3*theta) + 495*cos(4*theta)))/256.0;
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      case 3:
        deptheta = (-3.*sqrt(13./LAL_PI)*pow(cos(theta/2.0),5)*(-98. + 185.*cos(theta) - 110.*cos(2*theta) 
          + 55.*cos(3.*theta))*sin(theta/2.0))/32.0;
        out->re = deptheta * cos( 3.*phi );
        out->im = deptheta * sin( 3.*phi );
        break;

      case 4:
        deptheta = (sqrt(195./(2.0*LAL_PI))*pow(cos(theta/2.0),6)*(35. - 44.*cos(theta) 
          + 33.*cos(2*theta))*pow(sin(theta/2.0),2))/8.0;
        out->re = deptheta * cos( 4.*phi );
        out->im = deptheta * sin( 4.*phi );
        break;

      case 5:
        deptheta = -(sqrt(2145./LAL_PI)*pow(cos(theta/2.0),7)*(-1. + 3.*cos(theta))*pow(sin(theta/2.0),3))/2.0;
        out->re = deptheta * cos( 5.*phi );
        out->im = deptheta * sin( 5.*phi );
        break;

      case 6:
        deptheta = (3.*sqrt(715./LAL_PI)*pow(cos(theta/2.0),8)*pow(sin(theta/2.0),4))/2.0;
        out->re = deptheta * cos( 6.*phi );
        out->im = deptheta * sin( 6.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 6 */

  else if (L == 7)
  {
    switch ( M ) {
      case -7:
        deptheta = sqrt(15015./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*pow(sin(theta/2.0),9);
        out->re = deptheta * cos( -7.*phi );
        out->im = deptheta * sin( -7.*phi );
        break;

      case -6:
        deptheta = (sqrt(2145./LAL_PI)*pow(cos(theta/2.0),4)*(2. + 7.*cos(theta))*pow(sin(theta/2.0),8))/2.0;
        out->re = deptheta * cos( -6.*phi );
        out->im = deptheta * sin( -6.*phi );
        break;

      case -5:
        deptheta = (sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(93. + 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),7))/8.0;
        out->re = deptheta * cos( -5.*phi );
        out->im = deptheta * sin( -5.*phi );
        break;

      case -4:
        deptheta = (sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),2)*(140. + 285.*cos(theta) 
          + 156.*cos(2.*theta) + 91.*cos(3.*theta))*pow(sin(theta/2.0),6))/16.0;
        out->re = deptheta * cos( -4.*phi );
        out->im = deptheta * sin( -4.*phi );
        break;

      case -3:
        deptheta = (sqrt(15./(2.0*LAL_PI))*cos(theta/2.0)*(3115. + 5456.*cos(theta) + 4268.*cos(2.*theta) 
          + 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*pow(sin(theta/2.0),5))/128.0;
        out->re = deptheta * cos( -3.*phi );
        out->im = deptheta * sin( -3.*phi );
        break;

      case -2:
        deptheta = (sqrt(15./LAL_PI)*(5220. + 9810.*cos(theta) + 7920.*cos(2.*theta) + 5445.*cos(3.*theta) 
          + 2860.*cos(4.*theta) + 1001.*cos(5.*theta))*pow(sin(theta/2.0),4))/512.0;
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = (3.*sqrt(5./(2.0*LAL_PI))*cos(theta/2.0)*(1890. + 4130.*cos(theta) + 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) + 1430.*cos(4.*theta) + 1001.*cos(5*theta))*pow(sin(theta/2.0),3))/512.0;
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (3.*sqrt(35./LAL_PI)*cos(theta)*(109. + 132.*cos(2.*theta) 
          + 143.*cos(4.*theta))*pow(sin(theta),2))/512.0;
        out->re = deptheta;
        out->im = 0.;
        break;

      case 1:
        deptheta = (3.*sqrt(5./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(-1890. + 4130.*cos(theta) - 3080.*cos(2.*theta) 
          + 2805.*cos(3.*theta) - 1430.*cos(4.*theta) + 1001.*cos(5.*theta))*sin(theta/2.0))/512.0;
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (sqrt(15./LAL_PI)*pow(cos(theta/2.0),4)*(-5220. + 9810.*cos(theta) - 7920.*cos(2.*theta) 
          + 5445.*cos(3.*theta) - 2860.*cos(4.*theta) + 1001.*cos(5.*theta)))/512.0;
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      case 3:
        deptheta = -(sqrt(15./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(3115. - 5456.*cos(theta) + 4268.*cos(2.*theta) 
          - 2288.*cos(3.*theta) + 1001.*cos(4.*theta))*sin(theta/2.0))/128.0;
        out->re = deptheta * cos( 3.*phi );
        out->im = deptheta * sin( 3.*phi );
        break;

      case 4:
        deptheta = (sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),6)*(-140. + 285.*cos(theta) - 156.*cos(2*theta) 
          + 91.*cos(3.*theta))*pow(sin(theta/2.0),2))/16.0;
        out->re = deptheta * cos( 4.*phi );
        out->im = deptheta * sin( 4.*phi );
        break;

      case 5:
        deptheta = -(sqrt(165./(2.0*LAL_PI))*pow(cos(theta/2.0),7)*(93. - 104.*cos(theta) 
          + 91.*cos(2.*theta))*pow(sin(theta/2.0),3))/8.0;
        out->re = deptheta * cos( 5.*phi );
        out->im = deptheta * sin( 5.*phi );
        break;

      case 6:
        deptheta = (sqrt(2145./LAL_PI)*pow(cos(theta/2.0),8)*(-2. + 7.*cos(theta))*pow(sin(theta/2.0),4))/2.0;
        out->re = deptheta * cos( 6.*phi );
        out->im = deptheta * sin( 6.*phi );
        break;

      case 7:
        deptheta = -(sqrt(15015./(2.0*LAL_PI))*pow(cos(theta/2.0),9)*pow(sin(theta/2.0),5));
        out->re = deptheta * cos( 7.*phi );
        out->im = deptheta * sin( 7.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 7 */

  else if (L == 8)
  {
    switch ( M ) {
      case -8:
        deptheta = sqrt(34034./LAL_PI)*pow(cos(theta/2.0),6)*pow(sin(theta/2.0),10);
        out->re = deptheta * cos( -8.*phi );
        out->im = deptheta * sin( -8.*phi );
        break;

      case -7:
        deptheta = sqrt(17017./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(1. + 4.*cos(theta))*pow(sin(theta/2.0),9);
        out->re = deptheta * cos( -7.*phi );
        out->im = deptheta * sin( -7.*phi );
        break;

      case -6:
        deptheta = sqrt(255255./LAL_PI)*pow(cos(theta/2.0),4)*(1. + 2.*cos(theta))
          *sin(LAL_PI/4.0 - theta/2.0)*sin(LAL_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),8);
        out->re = deptheta * cos( -6.*phi );
        out->im = deptheta * sin( -6.*phi );
        break;

      case -5:
        deptheta = (sqrt(12155./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(19. + 42.*cos(theta) 
          + 21.*cos(2.*theta) + 14.*cos(3.*theta))*pow(sin(theta/2.0),7))/8.0;
        out->re = deptheta * cos( -5.*phi );
        out->im = deptheta * sin( -5.*phi );
        break;

      case -4:
        deptheta = (sqrt(935./(2.0*LAL_PI))*pow(cos(theta/2.0),2)*(265. + 442.*cos(theta) + 364.*cos(2.*theta) 
          + 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),6))/32.0;
        out->re = deptheta * cos( -4.*phi );
        out->im = deptheta * sin( -4.*phi );
        break;

      case -3:
        deptheta = (sqrt(561./(2.0*LAL_PI))*cos(theta/2.0)*(869. + 1660.*cos(theta) + 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) + 455.*cos(4.*theta) + 182.*cos(5.*theta))*pow(sin(theta/2.0),5))/128.0;
        out->re = deptheta * cos( -3.*phi );
        out->im = deptheta * sin( -3.*phi );
        break;

      case -2:
        deptheta = (sqrt(17./LAL_PI)*(7626. + 14454.*cos(theta) + 12375.*cos(2.*theta) + 9295.*cos(3.*theta) 
          + 6006.*cos(4.*theta) + 3003.*cos(5.*theta) + 1001.*cos(6.*theta))*pow(sin(theta/2.0),4))/512.0;
        out->re = deptheta * cos( -2.*phi );
        out->im = deptheta * sin( -2.*phi );
        break;

      case -1:
        deptheta = (sqrt(595./(2.0*LAL_PI))*cos(theta/2.0)*(798. + 1386.*cos(theta) + 1386.*cos(2.*theta) 
          + 1001.*cos(3.*theta) + 858.*cos(4.*theta) + 429.*cos(5.*theta) + 286.*cos(6.*theta))*pow(sin(theta/2.0),3))/512.0;
        out->re = deptheta * cos( -phi );
        out->im = deptheta * sin( -phi );
        break;

      case 0:
        deptheta = (3.*sqrt(595./LAL_PI)*(210. + 385.*cos(2.*theta) + 286.*cos(4.*theta) 
          + 143.*cos(6.*theta))*pow(sin(theta),2))/4096.0;
        out->re = deptheta;
        out->im = 0.;
        break;

      case 1:
        deptheta = (sqrt(595./(2.0*LAL_PI))*pow(cos(theta/2.0),3)*(798. - 1386.*cos(theta) + 1386.*cos(2.*theta) 
          - 1001.*cos(3.*theta) + 858.*cos(4.*theta) - 429.*cos(5.*theta) + 286.*cos(6.*theta))*sin(theta/2.0))/512.0;
        out->re = deptheta * cos( phi );
        out->im = deptheta * sin( phi );
        break;

      case 2:
        deptheta = (sqrt(17./LAL_PI)*pow(cos(theta/2.0),4)*(7626. - 14454.*cos(theta) + 12375.*cos(2.*theta) 
          - 9295.*cos(3.*theta) + 6006.*cos(4.*theta) - 3003.*cos(5.*theta) + 1001.*cos(6.*theta)))/512.0;
        out->re = deptheta * cos( 2.*phi );
        out->im = deptheta * sin( 2.*phi );
        break;

      case 3:
        deptheta = -(sqrt(561./(2.0*LAL_PI))*pow(cos(theta/2.0),5)*(-869. + 1660.*cos(theta) - 1300.*cos(2.*theta) 
          + 910.*cos(3.*theta) - 455.*cos(4.*theta) + 182.*cos(5.*theta))*sin(theta/2.0))/128.0;
        out->re = deptheta * cos( 3.*phi );
        out->im = deptheta * sin( 3.*phi );
        break;

      case 4:
        deptheta = (sqrt(935./(2.0*LAL_PI))*pow(cos(theta/2.0),6)*(265. - 442.*cos(theta) + 364.*cos(2.*theta) 
          - 182.*cos(3.*theta) + 91.*cos(4.*theta))*pow(sin(theta/2.0),2))/32.0;
        out->re = deptheta * cos( 4.*phi );
        out->im = deptheta * sin( 4.*phi );
        break;

      case 5:
        deptheta = -(sqrt(12155./(2.0*LAL_PI))*pow(cos(theta/2.0),7)*(-19. + 42.*cos(theta) - 21.*cos(2.*theta) 
          + 14.*cos(3.*theta))*pow(sin(theta/2.0),3))/8.0;
        out->re = deptheta * cos( 5.*phi );
        out->im = deptheta * sin( 5.*phi );
        break;

      case 6:
        deptheta = sqrt(255255./LAL_PI)*pow(cos(theta/2.0),8)*(-1. + 2.*cos(theta))*sin(LAL_PI/4.0 - theta/2.0)
          *sin(LAL_PI/4.0 + theta/2.0)*pow(sin(theta/2.0),4);
        out->re = deptheta * cos( 6.*phi );
        out->im = deptheta * sin( 6.*phi );
        break;

      case 7:
        deptheta = -(sqrt(17017./(2.0*LAL_PI))*pow(cos(theta/2.0),9)*(-1. + 4.*cos(theta))*pow(sin(theta/2.0),5));
        out->re = deptheta * cos( 7.*phi );
        out->im = deptheta * sin( 7.*phi );
        break;

      case 8:
        deptheta = sqrt(34034./LAL_PI)*pow(cos(theta/2.0),10)*pow(sin(theta/2.0),6);
        out->re = deptheta * cos( 8.*phi );
        out->im = deptheta * sin( 8.*phi );
        break;

      default:
        /*Error message informing that the chosen M is incompatible with L*/
        LALPrintError ("\n Inconsistent (L, M) values \n\n");
        XLAL_ERROR(__func__, XLAL_EINVAL);
        break;
    } /* switch (M) */
  }   /* L== 8 */

  else
  {
    /* Error message informing that L!=2 is not yet implemented*/
    LALPrintError ("\n These (L, M) values not implemented yet \n\n");
    XLAL_ERROR ( __func__, XLAL_EINVAL);
  }

  return 0;
}


/**
 * This function computes the regular scalar spherical harmonic.
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
    XLAL_ERROR( __func__, XLAL_EINVAL );
  }

  /* For some reason GSL will not take negative m */
  /* We will have to use the relation between sph harmonics of +ve and -ve m */
  XLAL_CALLGSL( gslStatus = gsl_sf_legendre_sphPlm_e((INT4)l, absM, cos(theta), &pLm ) );
  if (gslStatus != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( __func__, XLAL_EFUNC );
  }

  /* Compute the values for the spherical harmonic */
  y->re = pLm.val * cos(m * phi);
  y->im = pLm.val * sin(m * phi);

  /* If m is negative, perform some jiggery-pokery */
  if ( m < 0 && absM % 2  == 1 )
  {
    y->re = - y->re;
    y->im = - y->im;
  }

  return XLAL_SUCCESS;
}

/**
 * Computes the spin 2 weighted spherical harmonic. This function is now
 * deprecated and will be removed soon. All calls should be replaced with
 * calls to XLALSpinWeightedSphHarm.
 */
INT4 XLALSphHarm ( COMPLEX16 *out, /**< output */
                   UINT4   L,      /**< value of L */
                   INT4 M,         /**< value of M */
                   REAL4 theta,    /**< angle with respect to the z axis */
                   REAL4   phi     /**< angle with respect to the x axis */
                   )
{

  XLALPrintDeprecationWarning( "XLALSphHarm", "XLALSpinWeightedSphHarm" );

  if ( XLALSpinWeightedSphHarm( out, L, M, theta, phi ) == XLAL_FAILURE )
  {
    XLAL_ERROR( __func__, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;
}
