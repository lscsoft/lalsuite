/*
* Copyright (C) 2016 Lionel London
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


/* .................... */
/* HEADER SECTION       */
/* .................... */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/FileIO.h>

#include <lal/XLALError.h>

#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LALSimSphHarmSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/Window.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include "LALSimRingdownMMRDNS.h"
#include "LALSimRingdownData.h"
#include <lal/LALSimInspiralTestGRParams.h>
#include <lal/SphericalHarmonics.h>

/* note: use double-precision variables, but demand single-precision accuracy */
#define EPS LAL_REAL4_EPS
#define TINY LAL_REAL4_MIN

/* Mode included within the model */
UINT4 XLALMMRDNS_NUM_MODES = 18;
UINT4 XLALMMRDNS_NUM_MULTIPOLES = 14;

/* Define list of Spherical Multipoles to be used */
INT4 XLALMMRDNS_MULTIPOLES[14][2] = { {2,2}, {2,-2},
                                      {2,1}, {2,-1},
                                      {3,3}, {3,-3},
                                      {3,2}, {3,-2},
                                      {4,4}, {4,-4},
                                      {4,3}, {4,-3},
                                      {5,5}, {5,-5}
                                    };
/* Define list of QuasiNormal Modes to be used */
INT4 XLALMMRDNS_MODES[18][3] = { {2,2,0}, {2,-2,0},
                                 {2,2,1}, {2,-2,1},
                                 {3,3,0}, {3,-3,0},
                                 {3,3,1}, {3,-3,1},
                                 {4,4,0}, {4,-4,0},
                                 {5,5,0}, {5,-5,0},
                                 {2,1,0}, {2,-1,0},
                                 {3,2,0}, {3,-2,0},
                                 {4,3,0}, {4,-3,0}
                               };


/* Fit for Spherical-Spheroidal inner-products */
COMPLEX16 XLALQNM_YSPROD( REAL8 jf, UINT4 ll, INT4 mm, UINT4 l, INT4 m, UINT4 n ){

  //
  if ( mm<0 ) return pow(-1,l+ll) * conj( XLALQNM_YSPROD( jf, ll, -mm, l, -m, n ) );

  //
  COMPLEX16 ans = 0;
  REAL8 kap = XLALKAPPA( jf, ll, mm );
  REAL8 kap2 = kap*kap;
  REAL8 kap3 = kap2*kap;

  //
  if ( mm != m ){
    ans = 0;
  } else {

    if ( ll==2 && mm==2 && l==2 && n==0 ) {

      /* Spherical (ll,mm) = (2,2), Spheroidal (l,m,n) = (2,2,0) */
      ans = 0.99733430*cexp(6.28131482*I)  +
            7.53361604e-03 * (  14.59237582*cexp(5.06005232*I)*kap +
            28.76128310*cexp(1.62899784*I)*kap2 +
            14.51135075*cexp(4.63622189*I)*kap3 +
            1.96240925*cexp(3.01126131*I) ) /
            ( 1.0 +  0.88674277*cexp(6.22025328*I)*(kap) + 1.00197923*cexp(3.27370110*I)*kap2 + 0.08214826*cexp(2.49527907*I)*kap3 );

    } else if ( ll==2 && mm==2 && l==2 && n==1 ) {

      /* Spherical (ll,mm) = (2,2), Spheroidal (l,m,n) = (2,2,1) */
      ans = 0.99682510*cexp(6.27824650*I)  +
            2.07578052e-02 * (  15.07688773*cexp(4.83228133*I)*kap +
            31.13883493*cexp(1.58504497*I)*kap2 +
            15.44864065*cexp(4.67273252*I)*kap3 +
            0.71897499*cexp(2.80840474*I) ) / ( 1.0 +  0.80592193*cexp(0.25792353*I)*kap +
            0.69501985*cexp(3.68434414*I)*kap2 + 0.35613451*cexp(2.81292594*I)*kap3 );

    } else if ( ll==2 && mm==1 && l==2 && n==0 ) {

      /* Spherical (ll,mm) = (2,1), Spheroidal (l,m,n) = (2,1,0) */
      ans = 0.99715771*cexp(6.28150883*I)  +  6.35419212e-03 * (
            143454.06162814*cexp(4.50608133*I)*kap + 354688.47508323*cexp(1.73266939*I)*kap2 +
            240378.28541255*cexp(5.16291028*I)*kap3 + 6026.92526290*cexp(1.88805840*I) ) / (
            1.0 +  73780.21592675*cexp(1.41292540*I)*kap + 97493.61722607*cexp(4.53958149*I)*kap2 +
            34814.97776969*cexp(1.42067340*I)*kap3 );

    } else if ( ll==3 && mm==3 && l==3 && n==0 ) {

      /* Spherical (ll,mm) = (3,3), Spheroidal (l,m,n) = (3,3,0) */
      ans = 0.99569351*cexp(6.27848234*I)  +  1.45455997e-02 * (  7.21123775*cexp(0.62810602*I)*kap +
            6.53812261*cexp(4.62156760*I)*kap2 + 4.45103475*cexp(2.92276566*I)*kap3 +
            1.71127741*cexp(2.95271076*I) ) / ( 1.0 +  1.49740070*cexp(4.81030326*I)*kap +
            1.52875732*cexp(2.24686417*I)*kap2 + 0.52113989*cexp(5.68866116*I)*kap3 );

    } else if ( ll==3 && mm==3 && l==3 && n==1 ) {

      /* Spherical (ll,mm) = (3,3), Spheroidal (l,m,n) = (3,3,1) */
      ans = 0.99477991*cexp(6.26878602*I)  +  4.04784842e-02 * (  4.41127340*cexp(1.25014631*I)*kap +
            11.58762382*cexp(0.27959121*I)*kap2 + 17.32181411*cexp(3.79036872*I)*kap3 +
            0.67723741*cexp(2.57968532*I) ) / ( 1.0 +  3.87820757*cexp(5.42802120*I)*kap +
            3.49129316*cexp(2.52389826*I)*kap2 + 1.03677334*cexp(6.04982333*I)*kap3 );

    } else if ( ll==3 && mm==2 && l==3 && n==0 ) {

      /* Spherical (ll,mm) = (3,2), Spheroidal (l,m,n) = (3,2,0) */
      ans = 0.99009488*cexp(6.28042487*I)  +  2.36900822e-02 * (
            71893.02846884*cexp(1.23951995*I)*kap + 170547.75350114*cexp(5.03706013*I)*kap2 +
            129473.06324008*cexp(2.35898347*I)*kap3 + 1935.50631290*cexp(4.66800539*I) ) / (
            1.0 +  38206.44397720*cexp(4.36699213*I)*kap + 35811.01565015*cexp(0.82023766*I)*kap2 +
            8378.33659199*cexp(3.25885566*I)*kap3 );

    } else if ( ll==3 && mm==2 && l==2 && n==0 ) {

      /* Spherical (ll,mm) = (3,2), Spheroidal (l,m,n) = (2,2,0) */
      ans = 0.02059787*cexp(0.04742953*I)  +  6.91901622e-02 * (  2.76568066*cexp(2.13298172*I)*kap +
            3.95621590*cexp(4.65301637*I)*kap2 + 2.33636985*cexp(2.64436467*I)*kap3 +
            2.39898037*cexp(6.27665812*I) ) / ( 1.0 +  1.05953986*cexp(1.64486125*I)*kap +
            0.91307743*cexp(6.02859766*I)*kap2 + 0.69468169*cexp(3.33278829*I)*kap3 );

    } else if ( ll==3 && mm==2 && l==2 && n==1 ) {

      /* Spherical (ll,mm) = (3,2), Spheroidal (l,m,n) = (2,2,1) */
      ans = 0.02203038*cexp(0.16452240*I)  +  7.32333443e-02 * (  24.93221135*cexp(1.01806613*I)*kap +
            30.19732202*cexp(4.40469119*I)*kap2 + 11.27412168*cexp(2.98101589*I)*kap3 +
            2.43736370*cexp(6.19589588*I) ) / ( 1.0 +  11.39671828*cexp(0.85375374*I)*kap +
            10.91505457*cexp(2.66095722*I)*kap2 + 7.21957683*cexp(4.95919117*I)*kap3 );

    } else if ( ll==4 && mm==4 && l==4 && n==0 ) {

      /* Spherical (ll,mm) = (4,4), Spheroidal (l,m,n) = (4,4,0) */
      ans = 0.99478347*cexp(6.27759288*I)  +  2.47909254e-02 * (  6.51717929*cexp(0.79834752*I)*kap +
            7.77481969*cexp(4.24851207*I)*kap2 + 1.15767997*cexp(1.59054524*I)*kap3 +
            1.24344271*cexp(2.96155263*I) ) / ( 1.0 +  0.44548033*cexp(4.39121823*I)*kap +
            0.59436720*cexp(2.53165593*I)*kap2 + 0.24743370*cexp(5.97075782*I)*kap3 );

    } else if ( ll==4 && mm==3 && l==4 && n==0 ) {

      /* Spherical (ll,mm) = (4,3), Spheroidal (l,m,n) = (4,3,0) */
      ans = 0.98735236*cexp(6.27946510*I)  +  3.30278035e-02 * (
            700838.90444419*cexp(1.10668663*I)*kap + 1843013.12162912*cexp(4.88083492*I)*kap2 +
            1436658.36121341*cexp(2.14120510*I)*kap3 + 13844.41115020*cexp(4.56010157*I) ) / (
            1.0 +  356669.13441405*cexp(4.15650528*I)*kap +
            327401.16374453*cexp(0.63297714*I)*kap2 + 88620.86197318*cexp(3.24254228*I)*kap3 );

    } else if ( ll==4 && mm==3 && l==3 && n==0 ) {

      /* Spherical (ll,mm) = (4,3), Spheroidal (l,m,n) = (3,3,0) */
      ans = 0.02811216*cexp(0.04848786*I)  +  8.63825699e-02 * (  12.08728952*cexp(0.47221221*I)*kap +
            30.62626264*cexp(3.32811584*I)*kap2 + 16.32819913*cexp(6.17853751*I)*kap3 +
            2.36032677*cexp(6.26623516*I) ) / ( 1.0 +  4.96381416*cexp(0.45147018*I)*kap +
            6.25524720*cexp(3.05848542*I)*kap2 + 1.45381517*cexp(5.69553887*I)*kap3 );

    } else if ( ll==5 && mm==5 && l==5 && n==0 ) {

      /* Spherical (ll,mm) = (5,5), Spheroidal (l,m,n) = (5,5,0) */
      ans = 0.99433654*cexp(6.27734803*I)  +  3.12598758e-02 * (  6.55080118*cexp(0.93398061*I)*kap +
            8.05578622*cexp(4.28811876*I)*kap2 + 0.92971171*cexp(1.04364364*I)*kap3 +
            1.09036011*cexp(2.97119894*I) ) / ( 1.0 +  0.23127696*cexp(4.90814345*I)*kap +
            0.54957916*cexp(2.77620813*I)*kap2 + 0.21299995*cexp(6.15083758*I)*kap3 );

    } else {

      ans = 0;

    }

  } /* End of IF m==m */

  //
  return ans;

}

/* Interpolate tabulated data for QNM frequency */
COMPLEX16 XLALQNM_CW( REAL8 jf, UINT4 l, INT4 input_m, UINT4 n ){

  /* Setup GSL Splines: acceleration */
  gsl_interp_accel *acc_real = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_imag = gsl_interp_accel_alloc();

  /* Setup GSL Splines: spline "objects" */
  gsl_spline *CWREAL = gsl_spline_alloc(gsl_interp_cspline, QNMDataLength);
  gsl_spline *CWIMAG = gsl_spline_alloc(gsl_interp_cspline, QNMDataLength);

  /**/
  REAL8 cwreal=0, cwimag=0;
  INT4 m = abs(input_m);

  /**/
  if ( 2==l && 2==m && 0==n  ){
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW220REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW220REAL_DATA, CW220IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 2==l && 2==m && 1==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW221REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW221REAL_DATA, CW221IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 2==l && 1==m && 0==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW210REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW210REAL_DATA, CW210IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 3==l && 3==m && 0==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW330REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW330REAL_DATA, CW330IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 3==l && 3==m && 1==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW331REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW331REAL_DATA, CW331IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 3==l && 2==m && 0==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW320REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW320REAL_DATA, CW320IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 4==l && 4==m && 0==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW440REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW440REAL_DATA, CW440IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 4==l && 3==m && 0==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW430REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW430REAL_DATA, CW430IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  } else if ( 5==l && 5==m && 0==n ) {
    /**/
    gsl_spline_init( CWREAL, JF_DATA, CW550REAL_DATA, QNMDataLength);
    cwreal = gsl_spline_eval( CWREAL, jf, acc_real );
    /**/
    gsl_spline_init( CWIMAG, CW550REAL_DATA, CW550IMAG_DATA, QNMDataLength);
    cwimag = gsl_spline_eval( CWIMAG, cwreal, acc_imag );

  }

  //
  gsl_spline_free(CWREAL); gsl_spline_free(CWIMAG);
  gsl_interp_accel_free(acc_real); gsl_interp_accel_free(acc_imag);

  if ( input_m > 0 ) {
    return  cwreal + I*cwimag;
  } else {
    return -cwreal + I*cwimag;
  }

}


/* Interpolate tabulated data for QNM separation constant */
COMPLEX16 XLALQNM_SC( REAL8 jf, UINT4 l, INT4 input_m, UINT4 n ){

  /* Setup GSL Splines: acceleration */
  gsl_interp_accel *acc_real = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_imag = gsl_interp_accel_alloc();

  /* Setup GSL Splines: spline "objects" */
  gsl_spline *SCREAL = gsl_spline_alloc(gsl_interp_cspline, QNMDataLength);
  gsl_spline *SCIMAG = gsl_spline_alloc(gsl_interp_cspline, QNMDataLength);

  /**/
  REAL8 screal=0, scimag=0;
  INT4 m = abs(input_m);

  /**/
  if ( 2==l && 2==m && 0==n  ){
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC220REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC220IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 2==l && 2==m && 1==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC221REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC221IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 2==l && 1==m && 0==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC210REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC210IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 3==l && 3==m && 0==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC330REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC330IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 3==l && 3==m && 1==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC331REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC331IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 3==l && 2==m && 0==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC320REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC320IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 4==l && 4==m && 0==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC440REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC440IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 4==l && 3==m && 0==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC430REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC430IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  } else if ( 5==l && 5==m && 0==n ) {
    /**/
    gsl_spline_init( SCREAL, JF_DATA, SC550REAL_DATA, QNMDataLength);
    screal = gsl_spline_eval( SCREAL, jf, acc_real );
    /**/
    gsl_spline_init( SCIMAG, JF_DATA, SC550IMAG_DATA, QNMDataLength);
    scimag = gsl_spline_eval( SCIMAG, jf, acc_imag );

  }

  //
  gsl_spline_free(SCREAL); gsl_spline_free(SCIMAG);
  gsl_interp_accel_free(acc_real); gsl_interp_accel_free(acc_imag);

  /**/
  if ( input_m > 0 ) {
    return  screal + I*scimag;
  } else {
    return  screal - I*scimag;
  }

}

REAL8 XLALE_rad_nonspinning_UIB2016(REAL8 eta){
    REAL8 a1, a2, a3, a4;
    REAL8 eta2, eta3, eta4;

    /********************************************************************************************************************
     * Formula taken from eq. (21) and table VII of Jimenez, Keitel, Husa et al.  arXiv:1611.00332                      *
     * Note that the formula given there corresponds to E_rad = M - M_f = 1 - M_f (in units of the initial total mass). *
     * Then M = M_f + E_rad * M ==> M = M_f/(1-E_rad)                                                                   *
     ********************************************************************************************************************/

    a1   =  (1-(2.*sqrt(2.))/3.); /* Std error: 0 (Th) ; Rel error: 0   % */
    a2   =  0.5610;               /* Std error: 0.0026 ; Rel error: 0.5 % */
    a3   = -0.847;                /* Std error: 0.0270 ; Rel error: 3.2 % */
    a4   =  3.145;                /* Std error: 0.0690 ; Rel error: 2.2 % */

    eta2 = eta  * eta;
    eta3 = eta  * eta * eta;
    eta4 = eta2 * eta2;

    return a4*eta4 + a3*eta3 + a2*eta2 + a1*eta;

}

REAL8 XLALMf_to_M_nonspinning_UIB2016(REAL8 eta, REAL8 M_f){
    REAL8 M;

    M = M_f/(1-XLALE_rad_nonspinning_UIB2016(eta));

    return M;
}


/*
* Domain mapping for dimnesionless BH spin
*/
REAL8 XLALKAPPA( double jf, int l, int m ){
  /* */
  /* if ( jf > 1.0 ) XLAL_ERROR(XLAL_EDOM, "Spin (dimensionless Kerr parameter) must not be greater than 1.0\n"); */
  /**/
  double alpha = log( 2.0 - jf ) * 0.91023922662; // 1/log(3)=1.09861228867
  double beta  = 1.0 / ( 2.0 + l-abs(m) );
  return pow( alpha , beta );
}

/*
* Dimensionless QNM Frequencies: Note that name encodes date of writing
*/
COMPLEX16 XLALcomplexOmega( double kappa,  /* Domain mapping for  remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue*/
                          int n ) {     /* Overtone Number*/

  /* Predefine powers to increase efficiency*/
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken*/
  int m = abs(input_m);

  /**/
  complex double j = _Complex_I;

  /* Initialize the answer*/
  double complex ans = 0.0;

  /* Use If-Else ladder to determine which mode function to evaluate*/
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit.*/
    ans = 1.0 + kappa * (  1.557847  *cexp(2.903124*j) +
                           1.95097051*cexp(5.920970*j)*kappa +
                           2.09971716*cexp(2.760585*j)*kappa2 +
                           1.41094660*cexp(5.914340*j)*kappa3 +
                           0.41063923*cexp(2.795235*j)*kappa4  );

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.0 + kappa * (  1.870939*cexp(2.511247*j) +
                             2.71924916*cexp(5.424999*j)*kappa +
                             3.05648030*cexp(2.285698*j)*kappa2 +
                             2.05309677*cexp(5.486202*j)*kappa3 +
                             0.59549897*cexp(2.422525*j)*kappa4  );

    } else if ( 3==l && 2==m && 0==n ) {

      /* Define extra powers as needed*/
      double kappa5 = kappa4 * kappa;
      double kappa6 = kappa5 * kappa;

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.022464*cexp(0.004870*j) +
            0.24731213*cexp(0.665292*j)*kappa +
            1.70468239*cexp(3.138283*j)*kappa2 +
            0.94604882*cexp(0.163247*j)*kappa3 +
            1.53189884*cexp(5.703573*j)*kappa4 +
            2.28052668*cexp(2.685231*j)*kappa5 +
            0.92150314*cexp(5.841704*j)*kappa6;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 2.0 + kappa * (  2.658908*cexp(3.002787*j) +
                             2.97825567*cexp(6.050955*j)*kappa +
                             3.21842350*cexp(2.877514*j)*kappa2 +
                             2.12764967*cexp(5.989669*j)*kappa3 +
                             0.60338186*cexp(2.830031*j)*kappa4  );

    } else if ( 2==l && 1==m && 0==n ) {

      /* Define extra powers as needed*/
      double kappa5 = kappa4 * kappa;
      double kappa6 = kappa5 * kappa;

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.589113*cexp(0.043525*j) +
            0.18896353*cexp(2.289868*j)*kappa +
            1.15012965*cexp(5.810057*j)*kappa2 +
            6.04585476*cexp(2.741967*j)*kappa3 +
            11.12627777*cexp(5.844130*j)*kappa4 +
            9.34711461*cexp(2.669372*j)*kappa5 +
            3.03838318*cexp(5.791518*j)*kappa6;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.5 + kappa*(  2.095657*cexp(2.964973*j) +
                           2.46964352*cexp(5.996734*j)*kappa +
                           2.66552551*cexp(2.817591*j)*kappa2 +
                           1.75836443*cexp(5.932693*j)*kappa3 +
                           0.49905688*cexp(2.781658*j)*kappa4  );

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.5 + kappa*(  2.339070*cexp(2.649692*j) +
                           3.13988786*cexp(5.552467*j)*kappa +
                           3.59156756*cexp(2.347192*j)*kappa2 +
                           2.44895997*cexp(5.443504*j)*kappa3 +
                           0.70040804*cexp(2.283046*j)*kappa4  );

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.5 + kappa*(  0.205046*cexp(0.595328*j) +
                           3.10333396*cexp(3.016200*j)*kappa +
                           4.23612166*cexp(6.038842*j)*kappa2 +
                           3.02890198*cexp(2.826239*j)*kappa3 +
                           0.90843949*cexp(5.915164*j)*kappa4  );

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 2.5 + kappa*(  3.240455*cexp(3.027869*j) +
                           3.49056455*cexp(6.088814*j)*kappa +
                           3.74704093*cexp(2.921153*j)*kappa2 +
                           2.47252790*cexp(6.036510*j)*kappa3 +
                           0.69936568*cexp(2.876564*j)*kappa4  );

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  /* If m<0, then take the *Negative* conjugate */
  if ( input_m < 0 ) {
    /**/
    ans = -conj( ans );
  }

  return ans;

} /* END of complexOmega */


/*
* -------------------------------------------------------------------------------- *
* Low level models: QNM Frequencies, Separation Constants and Spheroidal Harmonics
* -------------------------------------------------------------------------------- *
*/

/*
* QNM Separation Constants: Note that name encodes date of writing
*/
COMPLEX16 XLALseparationConstant( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue */
                          int n ) {     /* Overtone Number */

  /* Predefine powers to increase efficiency */
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;
  double kappa5 = kappa4 * kappa;
  double kappa6 = kappa5 * kappa;
  double kappa7 = kappa6 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken */
  int m = abs(input_m);

  /**/
  complex double j = _Complex_I;

  /* Initialize the answer */
  double complex ans = 0.0;

  /* Use If-Else ladder to determine which mode function to evaluate */
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit. */
    ans = 0.55262405 + 6.54272463*cexp(0.24443847*j)*kappa + 5.94664565*cexp(3.88409012*j)*kappa2 + 5.39298183*cexp(1.01651284*j)*kappa3 + 3.58701474*cexp(4.53395559*j)*kappa4 + 1.36858235*cexp(1.57079633*j)*kappa5 + 0.18520700*cexp(4.71238898*j)*kappa6 ;

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 0.55229247 + 7.94074969*cexp(0.64081239*j)*kappa + 12.55567057*cexp(4.41980669*j)*kappa2 + 13.68518711*cexp(1.48039237*j)*kappa3 + 10.43884041*cexp(4.72599435*j)*kappa4 + 4.20731453*cexp(1.57079633*j)*kappa5 + 0.76232588*cexp(4.71238898*j)*kappa6;

    } else if ( 3==l && 2==m && 0==n ) {

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.  */
      ans = 8.18542769*cexp(6.27603422*j) + 1.55192720*cexp(1.79088081*j)*kappa + 8.94654695*cexp(5.18681710*j)*kappa2 + 28.66050158*cexp(1.63658858*j)*kappa3 + 60.77789497*cexp(4.72114050*j)*kappa4 + 72.13239907*cexp(1.57079633*j)*kappa5 + 45.38115278*cexp(4.71238898*j)*kappa6 + 11.84706755*cexp(1.57079633*j)*kappa7;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 13.05294185 + 9.23462388*cexp(0.14179514*j)*kappa + 7.09045393*cexp(3.69184561*j)*(kappa2) + 6.46711175*cexp(0.89254551*j)*(kappa3) + 4.96905278*cexp(4.43853588*j)*(kappa4) + 2.62299932*cexp(1.57079633*j)*(kappa5) + 0.58168681*cexp(4.71238898*j)*(kappa6);

    } else if ( 2==l && 1==m && 0==n ) {

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit. */
      ans = 3.10089518*cexp(6.25822093*j) + 2.69208437*cexp(1.95853947*j)*kappa + 16.58575360*cexp(4.98423605*j)*kappa2 + 57.84090876*cexp(1.63720921*j)*kappa3 + 118.21761290*cexp(4.72674943*j)*kappa4 + 135.93985738*cexp(1.57079633*j)*kappa5 + 82.81742189*cexp(4.71238898*j)*kappa6 + 20.85173245*cexp(1.57079633*j)*kappa7;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 5.70465254 + 7.94433155*cexp(0.18039136*j)*kappa + 6.55099749*cexp(3.77926384*j)*kappa2 + 6.31422768*cexp(0.93863733*j)*kappa3 + 4.81214531*cexp(4.46906976*j)*kappa4 + 2.38927043*cexp(1.57079633*j)*kappa5 + 0.48077965*cexp(4.71238898*j)*kappa6;

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 5.70318420 + 8.94926548*cexp(0.49834140*j)*kappa + 12.70528736*cexp(4.31772419*j)*kappa2 + 15.63533560*cexp(1.39390017*j)*kappa3 + 14.19057659*cexp(4.66913674*j)*kappa4 + 7.33238119*cexp(1.57079633*j)*kappa5 + 1.53701758*cexp(4.71238898*j)*kappa6;

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 15.28866348 + 0.75297352*cexp(0.22048290*j)*kappa + 3.64936150*cexp(0.61644055*j)*kappa2 + 8.02530641*cexp(4.82756576*j)*kappa3 + 12.47205664*cexp(1.67334685*j)*kappa4 + 10.30282199*cexp(4.71238898*j)*kappa5 + 3.52885679*cexp(1.57079633*j)*kappa6;

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 22.52292196 + 10.44137664*cexp(0.11607502*j)*kappa + 7.79707643*cexp(3.61247422*j)*kappa2 + 6.59989026*cexp(0.83792606*j)*kappa3 + 4.90367451*cexp(4.40545635*j)*kappa4 + 2.59913853*cexp(1.57079633*j)*kappa5 + 0.58985077*cexp(4.71238898*j)*kappa6;

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  /* If m<0, then take the conjugate */
  if ( input_m < 0 ) {
    /**/
    ans = conj( ans );
  }

  return ans;

} /* END of separationConstant */


/*
* Spheroidal Harmonic Normalization Constants: Note that name encodes date of writing
*/
double XLALspheroidalHarmonicNormalization( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless)*/
                          int l,        /* Polar eigenvalue*/
                          int input_m,  /* Azimuthal eigenvalue*/
                          int n ) {     /* Overtone Number*/

  /* Predefine powers to increase efficiency*/
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;
  double kappa5 = kappa4 * kappa;
  double kappa6 = kappa5 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken*/
  int m = abs(input_m);

  /* Initialize the answer*/
  double complex ans = 0.0;

  /* Use If-Else ladder to determine which mode function to evaluate */
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit.*/
    ans = 7.86366171 - 3.61447483*kappa + 3.48996689*kappa2 - 2.29347705*kappa3 + 0.74425069*kappa4 ;

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 7.86298703 - 3.59872285*kappa + 2.88459437*kappa2 - 0.92740734*kappa3 - 0.04445478*kappa4;

    } else if ( 3==l && 2==m && 0==n ) {

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.74845717 - 0.08157463*kappa + 1.03748092*kappa2 - 3.27926931*kappa3 + 7.24584503*kappa4 - 7.41316799*kappa5 + 3.06056035*kappa6;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.75389888 + 1.00111258*kappa + 1.55498487*kappa2 - 1.22344804*kappa3 + 1.64621074*kappa4;

    } else if ( 2==l && 1==m && 0==n ) {

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.04393302 - 0.06877527*kappa + 0.87671129*kappa2 - 3.92206769*kappa3 + 8.59631959*kappa4 - 8.52199526*kappa5 + 3.31150324*kappa6;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.51631915 + 0.16499714*kappa + 1.30114387*kappa2 - 0.83622153*kappa3 + 0.82020713*kappa4;

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.51530809 + 0.19285707*kappa + 0.96814190*kappa2 - 0.00547882*kappa3 + 0.24982172*kappa4;

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.39542385 - 0.09918352*kappa + 1.52850262*kappa2 - 5.09932727*kappa3 + 10.95647104*kappa4 - 10.99914124*kappa5 + 4.52212985*kappa6;

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 0.91349889 + 0.89568178*kappa + 2.54404526*kappa2 - 2.82437113*kappa3 + 3.28143852*kappa4;

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  return ans;

} /* END of XLALspheroidalHarmonicNormalization */


/*
* Final Spin (Dimensionless)
*/
UNUSED static double finalSpinFit( double eta );
UNUSED static double finalSpinFit( double eta ) {
  /* Implement Final Spin Fit from arXiv:1404.3197 */
  return eta*( 3.4339 - eta*( 3.7988 + eta*(5.7733 - 6.3780*eta) ) );
}

/*
* Final Mass (Dimensionless: Relative to M=1 initially)
*/
UNUSED static double finalMassFit( double eta );
UNUSED static double finalMassFit( double eta ) {
  /* Implement Final Mass Fit from arXiv:1404.3197 */
  return 1.0 - eta*(0.046297 - eta*(0.71006 + eta*(1.5028 - eta*(4.0124 - eta*0.28448))));
}


/* ------------------------------------------------
          Angular parameter functions
 ------------------------------------------------ */
// double XLALK1( int m, int s );
double XLALK1( int m, int s ){
  return 0.5*abs(m-s);
}
// double XLALK2( int m, int s );
double XLALK2( int m, int s ){
  return 0.5*abs(m+s);
}
// complex double XLALALPHA_RD( int m, int s, int p );
COMPLEX16 XLALALPHA_RD( int m, int s, int p ){
  /**/
  double k1 = XLALK1(m,s);
  return -2.0*(p+1.0)*(p+2.0*k1+1.0);
}
// COMPLEX16 XLALBETA_RD( int m, int s, int p, COMPLEX16 aw, COMPLEX16 A_lm );
COMPLEX16 XLALBETA_RD( int m, int s, int p, COMPLEX16 aw, COMPLEX16 A_lm ){
  /**/
  double k1 = XLALK1(m,s);
  double k2 = XLALK2(m,s);
  return  p*(p-1.0)+2.0*p*(k1+k2+1.0-2.0*aw) - ( 2.0*aw*(2.0*k1+s+1.0)-(k1+k2) * (k1+k2+1.0) ) - ( aw*aw + s*(s+1.0) + A_lm);
}
// COMPLEX16 XLALGAMMA_RD( int m, int s, int p, COMPLEX16 aw );
COMPLEX16 XLALGAMMA_RD( int m, int s, int p, COMPLEX16 aw ){
  /**/
  double k1 = XLALK1(m,s);
  double k2 = XLALK2(m,s);
  return 2.0*aw*(p+k1+k2+s);
}

/*
* QNM Ampliutde models for MMRDNS
* NOTE that the terms here differ from 1404.3197v3 for accurate relative phases
*/
COMPLEX16 XLALMMRDNSAmplitudeOverOmegaSquared( REAL8 eta, UINT4 l, INT4 m, INT4 n ){

  /* Handle m<0 through symmetry */
  if ( m<0 ) return pow(-1,l) * conj( XLALMMRDNSAmplitudeOverOmegaSquared(eta,l,-m,n) );

  /* Initialize the answer */
  COMPLEX16 ans = 0.0;

  /**/
  REAL8 eta2 = eta*eta;
  REAL8 eta3 = eta2*eta;
  REAL8 eta4 = eta3*eta;
  REAL8 eta5 = eta4*eta;

  /*  Evaluate QNM amplitude models for input l m n */
  if ( l==2 && m==2 && n==0 ) {

    /*A2201*/
    ans = 0.95846504*cexp(2.99318408*I)*eta + 0.47588079*cexp(0.82658128*I)*(eta2) + 1.23853419*cexp(2.30528861*I)*(eta3);

  } else if ( l==2 && m==2 && n==1 ) {

    /*A2211*/
    ans = 0.12750415*cexp(0.05809736*I)*eta + 1.18823931*cexp(1.51798243*I)*(eta2) + 8.27086561*cexp(4.42014780*I)*(eta3) + 26.23294960*cexp(1.16782950*I)*(eta4);

  } else if ( l==2 && m==1 && n==0 ) {

    /*A2101*/
    ans = sqrt(1.0-4.0*eta)  * (  0.47952344*cexp(5.96556090*I)*eta + 1.17357614*cexp(3.97472217*I)*(eta2) + 1.23033028*cexp(2.17322465*I)*(eta3)  );

  } else if ( l==3 && m==3 && n==0 ) {

    /*A3301*/
    ans = sqrt(1.0-4.0*eta)  * (  0.42472339*cexp(4.54734400*I)*eta + 1.47423728*cexp(2.70187807*I)*(eta2) + 4.31385024*cexp(5.12815819*I)*(eta3) + 15.72642073*cexp(2.25473854*I)*(eta4)  );

  } else if ( l==3 && m==3 && n==1 ) {

    /*A3311*/
    ans = sqrt(1.0-4.0*eta)  * (  0.14797161*cexp(2.03957081*I)*eta + 1.48738894*cexp(5.89538621*I)*(eta2) + 10.16366839*cexp(3.28354928*I)*(eta3) + 29.47859659*cexp(0.81061521*I)*(eta4)  );

  } else if ( l==3 && m==2 && n==0 ) {

    /*A3201*/
    ans = 0.19573228*cexp(0.54325509*I)*eta + 1.58299638*cexp(4.24509590*I)*(eta2) + 5.03380859*cexp(1.71003281*I)*(eta3) + 3.73662711*cexp(5.14735754*I)*(eta4);

  } else if ( l==4 && m==4 && n==0 ) {

    /*A4401*/
    ans = 0.25309908*cexp(5.16320109*I)*eta + 2.40404787*cexp(2.46899414*I)*(eta2) + 14.72733952*cexp(5.56235208*I)*(eta3) + 67.36237809*cexp(2.19824119*I)*(eta4) + 126.58579931*cexp(5.41735031*I)*(eta5);

  } else if ( l==4 && m==3 && n==0 ) {

    /*A4301*/
    ans = sqrt(1.0-4.0*eta)  * (  0.09383417*cexp(2.30765661*I)*eta + 0.82734483*cexp(6.10053234*I)*(eta2) + 3.33846327*cexp(3.87329126*I)*(eta3) + 4.66386840*cexp(1.75165690*I)*(eta4)  );

  } else if ( l==5 && m==5 && n==0 ) {

    /*A5501*/
    ans = sqrt(1.0-4.0*eta)  * (  0.15477314*cexp(1.06752431*I)*eta + 1.50914172*cexp(4.54983062*I)*(eta2) + 8.93331690*cexp(1.28981042*I)*(eta3) + 42.34309620*cexp(4.10035598*I)*(eta4) + 89.19466498*cexp(1.02508947*I)*(eta5)  );

  }

  /*NOTE that the MATLAB code used to perform the fitting uses a different convention when handling the real and imaginary parts of psi4 than we will use here. The conjugation below makes the output of MMRDNS consistent with nrutils, which injects no manual minus signs when handling psi4, but enforces a phase convention: m>0 has frequencies >0 (non-precessing). NOTE that this may change in the future if significantly precessing systems are found to not sufficiently obey this property. See https://github.com/llondon6/nrutils_dev/blob/master/nrutils/core/nrsc.py#L1714-L1728 for more details.*/
  ans = conj( ans );

  return ans;

}

/*
* Spheroical Harmonic Functions (Leaver's Formulation circa 1986/85)
*/
COMPLEX16 XLALSpinWeightedSpheroidalHarmonic( REAL8 jf,           /* Spin of remnant */
                   int l, int m, int n, /* QNM indices */
                   REAL8 theta,        /* polar angle */
                   REAL8 phi          /* azimuthal angle */
                 ) {

  if (m<0) return (1.0-2.0*((l+m)%2))*conj(XLALSpinWeightedSpheroidalHarmonic(jf, l, -m, n, LAL_PI-theta, LAL_PI+phi));

  /* Set spin weight */
  const REAL8 s = -2.0;

  /* Use tabulated cw and sc values from the core package*/
  COMPLEX16 cw, sc, aw;
  REAL8 kappa = XLALKAPPA(jf,l,m);

  /*cw = conj( XLALcomplexOmega( kappa, l, m, n ) );
  sc = XLALseparationConstant( kappa, l, m, n );*/
  cw = XLALQNM_CW(jf,l,m,n);
  sc = XLALQNM_SC(jf,l,m,n);

  /* Define dimensionless deformation parameter */
  aw = jf*cw;

  /* ------------------------------------------------
      Calculate the angular eigenfunction
   ------------------------------------------------ */

  /* Variable map for theta */
  REAL8 u = cos(theta);

  /* the non-sum part of eq 18 */
  COMPLEX16 X = 1.0;
  X = X * cexp(aw*u) * pow(1.0+u, XLALK1(m,s) );
  X = X * pow(1.0-u,XLALK2(m,s));

  /* NOTE that here we apply the normalization constant */
  X = X / XLALspheroidalHarmonicNormalization(kappa,l,m,n);
  // X = X / sqrt( cabs(XLALSpheroidalSpheroidalInnerProduct( jf, l,m,n, l,m,n, 512 )) );

  /* initial series values */
  COMPLEX16 a0 = 1.0;

  COMPLEX16 a1 = -XLALBETA_RD( m, s, 0, aw, sc )/XLALALPHA_RD( m, s, 0 );

  /* the sum part */
  COMPLEX16 Y = a0;
  Y = Y + a1*(1.0+u);
  COMPLEX16 S = 1.0;
  COMPLEX16 dY = 0;
  REAL8 xx = 0;
  COMPLEX16 a2 = 0;
  int done = 0;
  REAL8 k = 0, j = 0;
  k = 1.0;
  int kmax = 5e3;
  REAL8 et = 1e-8;
  REAL8 lastxx = 0;
  while ( ! done ) {
    k = k + 1.0;
    j = k - 1.0;
    a2 = -1.0*( XLALBETA_RD( m, s, j, aw, sc )*a1 + XLALGAMMA_RD(m,s,j,aw)*a0 ) / XLALALPHA_RD(m,s,j);
    dY = a2 * pow(1.0+u,k);
    Y = Y + dY;
    xx = cabs( dY );

    done = (k>=l) && ( (xx<et && k>30) || k>kmax );
    done = done || xx<et;

    /* Jump out of the WHILE-LOOP is the convergence is non-monotonic */
    if ( (k>3) && (lastxx<xx) ){
      done = true;
    }
    lastxx = xx;

    a0 = a1;
    a1 = a2;
  }
  // if (k > kmax) XLAL_ERROR(XLAL_EDOM, "sum did not converge\n");

  /* together now */
  S = X*Y*cexp( _Complex_I * m * phi );

  /* Use the same sign convention as spherical harmonics
  e.g http://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics#Calculating */
  REAL8 C = 1.0;
  C = C * pow(-1,fmax(-m,-s)) * pow(-1,l);
  S = C * S;

  /**/
  // return Y;
  return  S;

} /* END of Spheroical Harmonic function */

/* SpheroidalHarmonic Plus and Cross analogous to ones for SphericalHarmonics in LALSimBlackHoleRingdownTiger.c */

REAL8 XLALSimSpheroidalHarmonicPlus(REAL8 jf, UINT4 l, INT4 m, UINT4 n, REAL8 iota, REAL8 phi){
        REAL8 Yplus = 0.0;
        Yplus = creal(XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, phi) + (1.0 - 2.0*(l % 2))*XLALSpinWeightedSpheroidalHarmonic(jf, l, -m, n, iota, phi));
        return Yplus;
}

REAL8 XLALSimSpheroidalHarmonicCross(REAL8 jf, UINT4 l, INT4 m, UINT4 n, REAL8 iota, REAL8 phi){
        REAL8 Ycross = 0.0;
        Ycross = creal(XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, phi) - (1.0 - 2.0*(l % 2))*XLALSpinWeightedSpheroidalHarmonic(jf, l, -m, n, iota, phi));
        return Ycross;
}


/* Function to compute spherical and spheroidal inner-products */
COMPLEX16 XLALSphericalSpheroidalInnerProduct( REAL8 jf, UINT4 ll, INT4 mm, UINT4 l, INT4 m, UINT4 n, UINT4 N ){
  /**/
  REAL8 dtheta = LAL_PI/(N-1);
  /**/
  COMPLEX16 th = 0;
  COMPLEX16 ans = 0;
  /* NOTE that the k=0 and k=N terms of the trapezoidal rule do not contribute for integrands proportional to sin(th). */
  for ( UINT4 k=1; k<N; k++ ) {
    /**/
    th   = k*dtheta;
    /**/
    ans += sin(th) * XLALSpinWeightedSpheroidalHarmonic( jf, l,m,n, th, 0 ) * conj( XLALSpinWeightedSphericalHarmonic( th, 0, -2, ll, mm  ) );
  }
  /**/
  ans *= dtheta * LAL_TWOPI;
  /**/
  return ans;
}

/* Function to compute SPHEROIDAL and spheroidal inner-products */
COMPLEX16 XLALSpheroidalSpheroidalInnerProduct( REAL8 jf, UINT4 ll, INT4 mm, UINT4 nn, UINT4 l, INT4 m, UINT4 n, UINT4 N ){
  /**/
  REAL8 dtheta = LAL_PI/(N-1);
  /**/
  COMPLEX16 th = 0;
  COMPLEX16 ans = 0;
  /* NOTE that the k=0 and k=N terms of the trapezoidal rule do not contribute for integrands proportional to sin(th). */
  for ( UINT4 k=1; k<N; k++ )
  {
    /**/
    th   = k*dtheta;
    /**/
    ans += sin(th) * XLALSpinWeightedSpheroidalHarmonic( jf, l,m,n, th, 0 ) * conj( XLALSpinWeightedSpheroidalHarmonic( jf, ll,mm,nn, th, 0 ) );
  }
  /**/
  ans *= dtheta * LAL_TWOPI;
  /**/
  return ans;
}


/* XLALSimRingdownMMRDNSFD: Frequency domain waveformgenerator for all QNM with angular dependence */
int XLALSimRingdownMMRDNSFD(
        COMPLEX16FrequencySeries **hptilde,          /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde,          /**< OUTPUT FD h_x polarization */
        const REAL8 deltaF,                          /**< Frequency resolution (Hz) */
        const REAL8 fStart,                          /**< Start GW frequency (Hz) */
        const REAL8 fEnd,                            /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams       /**< testing GR parameters */
    ){
        /* Perform some initial checks */
        if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
        if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
        if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
        if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
        if (r <= 0) XLAL_ERROR(XLAL_EDOM);

        /* Declarations */
        UINT4 jStart, jMax;
        REAL8 f_max;
        LIGOTimeGPS tC = {0,0};
        //XLALGPSAdd(&tC, -1 / deltaF);

        REAL8 dfreq220=0., dfreq221=0., dfreq330=0., dfreq331=0., dfreq440=0., dfreq550=0., dfreq210=0., dfreq320=0., dfreq430=0.;
        REAL8 dtau220=0.,  dtau221=0.,  dtau330=0.,  dtau331=0.,  dtau440=0.,  dtau550=0.,  dtau210=0.,  dtau320=0.,  dtau430=0.;

        COMPLEX16FrequencySeries *hp220, *hp221, *hp330, *hp331, *hp440, *hp550, *hp210, *hp320, *hp430;
        COMPLEX16FrequencySeries *hc220, *hc221, *hc330, *hc331, *hc440, *hc550, *hc210, *hc320, *hc430;

        /* Get nonGRparams */
        char *nonGRParamName = malloc(512*sizeof(char));

        sprintf(nonGRParamName,"dfreq220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;

        /* allocate hptilde and hctilde */
        if ( fEnd == 0. ){
          f_max = 2048.0;
        } else {
          f_max = fEnd;
        }
        jMax = (UINT4)(f_max/deltaF + 1);

        XLALSimRingdownGenerateSingleModeFD(&hp220, &hc220, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 2, 2, 0, r, dfreq220, dtau220);
        XLALSimRingdownGenerateSingleModeFD(&hp221, &hc221, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 2, 2, 1, r, dfreq221, dtau221);
        XLALSimRingdownGenerateSingleModeFD(&hp330, &hc330, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 3, 3, 0, r, dfreq330, dtau330);
        XLALSimRingdownGenerateSingleModeFD(&hp331, &hc331, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 3, 3, 1, r, dfreq331, dtau331);
        XLALSimRingdownGenerateSingleModeFD(&hp440, &hc440, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 4, 4, 0, r, dfreq440, dtau440);
        XLALSimRingdownGenerateSingleModeFD(&hp550, &hc550, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 5, 5, 0, r, dfreq550, dtau550);
        XLALSimRingdownGenerateSingleModeFD(&hp210, &hc210, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 2, 1, 0, r, dfreq210, dtau210);
        XLALSimRingdownGenerateSingleModeFD(&hp320, &hc320, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 3, 2, 0, r, dfreq320, dtau320);
        XLALSimRingdownGenerateSingleModeFD(&hp430, &hc430, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 4, 3, 0, r, dfreq430, dtau430);

        *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hptilde)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hptilde)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hctilde)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hctilde)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

        jStart   = (UINT4) ceil(fStart / deltaF);

        for ( UINT4 j=jStart ; j<jMax ; j++ ) {
          (*hptilde)->data->data[j] = hp220->data->data[j];
          (*hptilde)->data->data[j] += hp221->data->data[j];
          (*hptilde)->data->data[j] += hp330->data->data[j];
          (*hptilde)->data->data[j] += hp331->data->data[j];
          (*hptilde)->data->data[j] += hp440->data->data[j];
          (*hptilde)->data->data[j] += hp550->data->data[j];
          (*hptilde)->data->data[j] += hp210->data->data[j];
          (*hptilde)->data->data[j] += hp320->data->data[j];
          (*hptilde)->data->data[j] += hp430->data->data[j];

          (*hctilde)->data->data[j] = hc220->data->data[j];
          (*hctilde)->data->data[j] += hc221->data->data[j];
          (*hctilde)->data->data[j] += hc330->data->data[j];
          (*hctilde)->data->data[j] += hc331->data->data[j];
          (*hctilde)->data->data[j] += hc440->data->data[j];
          (*hctilde)->data->data[j] += hc550->data->data[j];
          (*hctilde)->data->data[j] += hc210->data->data[j];
          (*hctilde)->data->data[j] += hc320->data->data[j];
          (*hctilde)->data->data[j] += hc430->data->data[j];
        }

        /* Destroy the COMPLEX16FrequencySeries object */
        if (hp220) XLALDestroyCOMPLEX16FrequencySeries(hp220);
        if (hp221) XLALDestroyCOMPLEX16FrequencySeries(hp221);
        if (hp330) XLALDestroyCOMPLEX16FrequencySeries(hp330);
        if (hp331) XLALDestroyCOMPLEX16FrequencySeries(hp331);
        if (hp440) XLALDestroyCOMPLEX16FrequencySeries(hp440);
        if (hp550) XLALDestroyCOMPLEX16FrequencySeries(hp550);
        if (hp210) XLALDestroyCOMPLEX16FrequencySeries(hp210);
        if (hp320) XLALDestroyCOMPLEX16FrequencySeries(hp320);
        if (hp430) XLALDestroyCOMPLEX16FrequencySeries(hp430);

        if (hc220) XLALDestroyCOMPLEX16FrequencySeries(hc220);
        if (hc221) XLALDestroyCOMPLEX16FrequencySeries(hc221);
        if (hc330) XLALDestroyCOMPLEX16FrequencySeries(hc330);
        if (hc331) XLALDestroyCOMPLEX16FrequencySeries(hc331);
        if (hc440) XLALDestroyCOMPLEX16FrequencySeries(hc440);
        if (hc550) XLALDestroyCOMPLEX16FrequencySeries(hc550);
        if (hc210) XLALDestroyCOMPLEX16FrequencySeries(hc210);
        if (hc320) XLALDestroyCOMPLEX16FrequencySeries(hc320);
        if (hc430) XLALDestroyCOMPLEX16FrequencySeries(hc430);

        /* Cleanup nonGRParamName */
        free(nonGRParamName);

  return XLAL_SUCCESS;
}

/* XLALSimRingdownGenerateSingleModeFD: Frequency domain waveformgenerator for single QNM with angular dependence */
int XLALSimRingdownGenerateSingleModeFD(
        COMPLEX16FrequencySeries **hptilde_lmn,      /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde_lmn,      /**< OUTPUT FD h_x polarization */
        const REAL8 deltaF,                          /**< Frequency resolution (Hz) */
        const REAL8 fStart,                          /**< Start GW frequency (Hz) */
        const REAL8 fEnd,                            /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset (in rad) */
        UINT4 l,                                     /**< Polar eigenvalue */
        UINT4 m,                                     /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau                                   /**< relative shift in the damping time parameter */
        ){

        /* Perform some initial checks */
        if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
        if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
        if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
        if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
        if (r <= 0) XLAL_ERROR(XLAL_EDOM);

        /* Declarations */
        UINT4 jStart, jMax;
        COMPLEX16 h_f, hconj_mf;
        REAL8 f, f_max;
        LIGOTimeGPS tC = {0,0};
        //XLALGPSAdd(&tC, -1 / deltaF);

        REAL8 kappa   = XLALKAPPA( jf, l, m );
        REAL8 Mf_sec  = Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
        REAL8 r_sec   = r/LAL_C_SI;

        COMPLEX16 A_lmn, S_lmn, Omega_lmn, Prefactor;

        /* Mode Component Calculation*/
        Omega_lmn = XLALcomplexOmega(kappa, l, m, n)/Mf_sec;
        A_lmn = XLALMMRDNSAmplitudeOverOmegaSquared(eta, l, m, n);
        Omega_lmn = creal(Omega_lmn)*(1.+dfreq) + I * cimag(Omega_lmn)/(1.+dtau);
        S_lmn = XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, 0.0);
        Prefactor = cexp(I*phi_offset)*(Mf_sec/r_sec)*(A_lmn*S_lmn)*(-I);

        /* allocate htilde_p and htilde_c */
        /* The COMPLEX16FrequencySeries has to be destroyed by whom created it. */
        if ( fEnd == 0. ){
          f_max = 2048.0;
        } else {
          f_max = fEnd;
        }

        jMax = (UINT4)(f_max/deltaF + 1);

        *hptilde_lmn = XLALCreateCOMPLEX16FrequencySeries("hptilde_lmn: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hptilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hptilde_lmn)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hptilde_lmn)->sampleUnits, &(*hptilde_lmn)->sampleUnits, &lalSecondUnit);
        *hctilde_lmn = XLALCreateCOMPLEX16FrequencySeries("hctilde_lmn: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hctilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hctilde_lmn)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hctilde_lmn)->sampleUnits, &(*hctilde_lmn)->sampleUnits, &lalSecondUnit);

        jStart   = (UINT4) ceil(fStart / deltaF);
        f        = jStart*deltaF;
        h_f      = 0.0;
        hconj_mf = 0.0;

        for ( UINT4 j=jStart ; j<jMax ; j++ ) {
        h_f      =      Prefactor/(Omega_lmn-LAL_TWOPI*f)*cexp(I*Omega_lmn*10.0*Mf_sec-I*LAL_TWOPI*f*10.0*Mf_sec);
        hconj_mf = conj(Prefactor/(Omega_lmn+LAL_TWOPI*f)*cexp(I*Omega_lmn*10.0*Mf_sec+I*LAL_TWOPI*f*10.0*Mf_sec));

        (*hptilde_lmn)->data->data[j] = 0.5 * (h_f + hconj_mf);
        (*hctilde_lmn)->data->data[j] = 0.5 * I * (h_f - hconj_mf);

        h_f      = 0.0;
        hconj_mf = 0.0;
        f+=deltaF;

        }

  return XLAL_SUCCESS;
}

/* Functions analogous to functions used in LALSimBlackHoleRingdownTiger.c for Kamaretsos fits */
int XLALSimRingdownMMRDNS_time(
        REAL8TimeSeries **hplus,                     /**< plus-polarization waveform [returned] */
        REAL8TimeSeries **hcross,                    /**< cross-polarization waveform [returned] */
        const LIGOTimeGPS *T0,                       /**< start time of ringdown => NEEDS TO BE CHECKED! */
        REAL8 deltaT,                                /**< sampling interval (s) */
        UINT4 Nsamples,                              /**< Number of samples (effective T_End) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams       /**< testing GR parameters */
    ){
        /* Perform some initial checks */
        if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
        if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
        if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
        if (r <= 0) XLAL_ERROR(XLAL_EDOM);

        /* Declarations */
        REAL8 dfreq220=0., dfreq221=0., dfreq330=0., dfreq331=0., dfreq440=0., dfreq550=0., dfreq210=0., dfreq320=0., dfreq430=0.;
        REAL8 dtau220=0.,  dtau221=0.,  dtau330=0.,  dtau331=0.,  dtau440=0.,  dtau550=0.,  dtau210=0.,  dtau320=0.,  dtau430=0.;
        COMPLEX16TimeSeries *h220,  *h221,  *h330,  *h331,  *h440,  *h550,  *h210,  *h320,  *h430;
        COMPLEX16TimeSeries *h2m20, *h2m21, *h3m30, *h3m31, *h4m40, *h5m50, *h2m10, *h3m20, *h4m30;

        /* Get nonGRparams */
        char *nonGRParamName = malloc(512*sizeof(char));

        sprintf(nonGRParamName,"dfreq220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;

        /* Use total mass to agree with NR conventions. Inside the template t=0 corresponds to 10M after the merger. 
           Tend is computed inside LALInferenceTemplate.c to give the same number of samples to the window and to the generator.
           Inside this routine Tend is effectively represented by Nsamples passed as input.*/
        REAL8 M        = XLALMf_to_M_nonspinning_UIB2016(eta, Mf);
        REAL8 Tstart   = 0.0*M*LAL_MTSUN_SI/LAL_MSUN_SI;


        /* Compute the modes seperately */
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h220, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2, 2, 0, r, dfreq220, dtau220, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h221, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2, 2, 1, r, dfreq221, dtau221, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h330, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3, 3, 0, r, dfreq330, dtau330, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h331, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3, 3, 1, r, dfreq331, dtau331, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h440, T0, deltaT, Mf, jf, eta, iota, phi_offset, 4, 4, 0, r, dfreq440, dtau440, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h550, T0, deltaT, Mf, jf, eta, iota, phi_offset, 5, 5, 0, r, dfreq550, dtau550, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h210, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2, 1, 0, r, dfreq210, dtau210, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h320, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3, 2, 0, r, dfreq320, dtau320, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h430, T0, deltaT, Mf, jf, eta, iota, phi_offset, 4, 3, 0, r, dfreq430, dtau430, Nsamples, Tstart );

        /* Compute the modes seperately */
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h2m20, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2,-2, 0, r, dfreq220, dtau220, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h2m21, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2,-2, 1, r, dfreq221, dtau221, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h3m30, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3,-3, 0, r, dfreq330, dtau330, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h3m31, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3,-3, 1, r, dfreq331, dtau331, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h4m40, T0, deltaT, Mf, jf, eta, iota, phi_offset, 4,-4, 0, r, dfreq440, dtau440, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h5m50, T0, deltaT, Mf, jf, eta, iota, phi_offset, 5,-5, 0, r, dfreq550, dtau550, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h2m10, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2,-1, 0, r, dfreq210, dtau210, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h3m20, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3,-2, 0, r, dfreq320, dtau320, Nsamples, Tstart );
        XLALSimRingdownGenerateSingleModeMMRDNS_time( &h4m30, T0, deltaT, Mf, jf, eta, iota, phi_offset, 4,-3, 0, r, dfreq430, dtau430, Nsamples, Tstart );

        /* Add the modes to get the final waveform and  get cross and plus polarization */
        *hplus = XLALCreateREAL8TimeSeries( "hplus: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
        if (!(*hplus)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hplus)->data->data, 0, Nsamples * sizeof(REAL8));
        *hcross = XLALCreateREAL8TimeSeries( "hcross: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
        if (!(*hcross)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hcross)->data->data, 0, Nsamples * sizeof(REAL8));


        COMPLEX16 h_val = 0.0;
        for ( UINT4 i=0 ; i<Nsamples ; i++ )
	      {

          h_val = 0;

          h_val += h220->data->data[i];
          h_val += h221->data->data[i];
          h_val += h330->data->data[i];
          h_val += h331->data->data[i];
          h_val += h440->data->data[i];
          h_val += h550->data->data[i];
          h_val += h210->data->data[i];
          h_val += h320->data->data[i];
          h_val += h430->data->data[i];

          h_val += h2m20->data->data[i];
          h_val += h2m21->data->data[i];
          h_val += h3m30->data->data[i];
          h_val += h3m31->data->data[i];
          h_val += h4m40->data->data[i];
          h_val += h5m50->data->data[i];
          h_val += h2m10->data->data[i];
          h_val += h3m20->data->data[i];
          h_val += h4m30->data->data[i];

          (*hplus)->data->data[i]  =        creal(h_val);
          (*hcross)->data->data[i] = -1.0 * cimag(h_val);

        }

        /* Destroy the COMPLEX16TimeSeries object */
        if (h220) XLALDestroyCOMPLEX16TimeSeries(h220);
        if (h221) XLALDestroyCOMPLEX16TimeSeries(h221);
        if (h330) XLALDestroyCOMPLEX16TimeSeries(h330);
        if (h331) XLALDestroyCOMPLEX16TimeSeries(h331);
        if (h440) XLALDestroyCOMPLEX16TimeSeries(h440);
        if (h550) XLALDestroyCOMPLEX16TimeSeries(h550);
        if (h210) XLALDestroyCOMPLEX16TimeSeries(h210);
        if (h320) XLALDestroyCOMPLEX16TimeSeries(h320);
        if (h430) XLALDestroyCOMPLEX16TimeSeries(h430);
        /**/
        if (h2m20) XLALDestroyCOMPLEX16TimeSeries(h2m20);
        if (h2m21) XLALDestroyCOMPLEX16TimeSeries(h2m21);
        if (h3m30) XLALDestroyCOMPLEX16TimeSeries(h3m30);
        if (h3m31) XLALDestroyCOMPLEX16TimeSeries(h3m31);
        if (h4m40) XLALDestroyCOMPLEX16TimeSeries(h4m40);
        if (h5m50) XLALDestroyCOMPLEX16TimeSeries(h5m50);
        if (h2m10) XLALDestroyCOMPLEX16TimeSeries(h2m10);
        if (h3m20) XLALDestroyCOMPLEX16TimeSeries(h3m20);
        if (h4m30) XLALDestroyCOMPLEX16TimeSeries(h4m30);

        free(nonGRParamName);
        return 0;

}


int XLALSimRingdownGenerateSingleModeMMRDNS_time(
        COMPLEX16TimeSeries **htilde_lmn,            /**< OUTPUT TD waveform mode lmn */
        const LIGOTimeGPS *T0,                       /**< start time of ringdown */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic ORBITAL phase offset (in rad) */
        UINT4 l,                                     /**< Polar eigenvalue */
        INT4 m,                                     /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        UINT4 Nsamples,                              /**< waveform length */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */
        ){

        /*  */
        COMPLEX16 S = XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, phi_offset);

        /* allocate htilde_lmn */
        if ( ! *htilde_lmn ){
          *htilde_lmn = XLALCreateCOMPLEX16TimeSeries("htilde_lmn: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
          if (!(*htilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
          memset((*htilde_lmn)->data->data, 0, Nsamples * sizeof(COMPLEX16));
        }

        XLALSimRingdownGenerateSingleBareModeMMRDNS_time( htilde_lmn, T0, deltaT, Mf, jf, eta, l, m, n, r, dfreq, dtau, Nsamples, Tstart );

        COMPLEX16 hk = 0;
        for ( UINT4 k=0; k<Nsamples; k++ )
        {
          /**/
          hk = S * ((*htilde_lmn)->data->data[k]);
          (*htilde_lmn)->data->data[k] = hk;
        }

  return XLAL_SUCCESS;

}

/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
/* Generate a single TD QNM without angular dependence */
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int XLALSimRingdownGenerateSingleBareModeMMRDNS_time(
        COMPLEX16TimeSeries **htilde_lmn,            /**< OUTPUT TD waveform mode lmn */
        UNUSED const LIGOTimeGPS *T0,                       /**< start time of ringdown */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        UINT4 l,                                     /**< Polar eigenvalue */
        INT4 m,                                     /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        UINT4 Nsamples,                              /**< waveform length */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */
        ){

        /* Declarations: M is the initial total mass, Mf is the remnant mass */
        COMPLEX16 h_lmn     = 0.0;
        REAL8 M             = XLALMf_to_M_nonspinning_UIB2016(eta, Mf);
        REAL8 M_sec         = M*LAL_MTSUN_SI/LAL_MSUN_SI;
        REAL8 Mf_sec        = Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
        REAL8 r_sec         = r/LAL_C_SI;
        REAL8 t             = 0.0;
        COMPLEX16 A_lmn     = 0.0;
        COMPLEX16 A         = 0.0;
        COMPLEX16 CW        = 0.0;
        COMPLEX16 Omega_lmn = 0.0;


        // if (m<0) return (1.0-2.0*((l+m)%2))*conj(XLALSimRingdownGenerateSingleBareModeMMRDNS_time(htilde_lmn, T0, deltaT, Mf, jf, eta, l, m, n, r, dfreq, dtau, Nsamples, Tstart));

        /* Mode Component Calculation*/
        A_lmn 	     = XLALMMRDNSAmplitudeOverOmegaSquared(eta, l, m, n);
        A            = A_lmn*M_sec/(r_sec);

        /* NOTE: Whether a fit for interopolation is used here makes not appreciable difference. But let's use a fit for consistency with the spheroidal function. */
        CW        = conj( XLALQNM_CW(jf,l,m,n) ) / Mf_sec ;
        Omega_lmn = creal(CW)*(1.0+dfreq) + I*cimag(CW)/(1.0+dtau);

        /* allocate htilde_lmn */
        *htilde_lmn = XLALCreateCOMPLEX16TimeSeries("htilde_lmn: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
        if (!(*htilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
        memset((*htilde_lmn)->data->data, 0, Nsamples * sizeof(COMPLEX16));

        /* fill waveform */
        for ( UINT4 i=0 ; i<Nsamples ; i++ )
	      {
            t     = Tstart+i*deltaT;

            h_lmn = A * cexp( I * t * Omega_lmn );

            (*htilde_lmn)->data->data[i] = h_lmn;
            h_lmn = 0.0;
            t     = 0.0;
        }

  return XLAL_SUCCESS;

}



/* Full waveform generator that uses the Spherical basis. */
int XLALSimRingdownGenerateFullSphericalWaveform_time
(
    REAL8TimeSeries **hplus,                     /**< OUTPUT TD waveform */
    REAL8TimeSeries **hcross,                    /**< OUTPUT TD waveform */
    const LIGOTimeGPS *T0,                       /**< start time of ringdown => NEEDS TO BE CHECKED! */
    REAL8 deltaT,                                /**< sampling interval (s) */
    UINT4 Nsamples,                              /**< Number of samples (effective T_End) */
    REAL8 Mf,                                    /**< Final BH Mass (kg) */
    REAL8 jf,                                    /**< Final BH dimensionaless spin */
    REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
    REAL8 iota,                                  /**< inclination angle (in rad) */
    REAL8 phi_offset,                            /**< intrinsic phase offset */
    REAL8 r,                                     /**< distance of source (m) */
    LALSimInspiralTestGRParam *nonGRparams
)
{

    /* Perform some initial checks */
    if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
    if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
    if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);
    
    /* allocate htilde */
    COMPLEX16TimeSeries *htilde = XLALCreateCOMPLEX16TimeSeries("htilde: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
    if (! htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, Nsamples * sizeof(COMPLEX16));

    /* initialize outputs */
    *hplus = XLALCreateREAL8TimeSeries( "hplus: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
    if (!(*hplus)) XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, Nsamples * sizeof(REAL8));
    *hcross = XLALCreateREAL8TimeSeries( "hcross: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
    if (!(*hcross)) XLAL_ERROR(XLAL_EFUNC);
    memset((*hcross)->data->data, 0, Nsamples * sizeof(REAL8));

    /* Declare helper variables */
    COMPLEX16 Prefactor=0;
    UINT4 ll=0; INT4 mm=0;
    
    /* Use total mass to agree with NR conventions. Inside the template t=0 corresponds to 10M after the psi_22 peak luminosity.
       Tend is computed inside LALInferenceTemplate.c to give the same number of samples to the window and to the generator.
       Inside this routine Tend is effectively represented by Nsamples passed as input.*/
    REAL8 M        = XLALMf_to_M_nonspinning_UIB2016(eta, Mf);
    REAL8 Tstart   = 0.0*M*LAL_MTSUN_SI/LAL_MSUN_SI;
    
    
  /* FOR all MULTIPOLES in MMRDNS */
    for ( UINT4 k=0; k<XLALMMRDNS_NUM_MULTIPOLES; k++ )
    {
        /* Extract the mode indeces from the master list */
        ll = XLALMMRDNS_MULTIPOLES[k][0]; mm = XLALMMRDNS_MULTIPOLES[k][1];
        /* Here the prefactor is the Spin weighted spherical harmonic */
        Prefactor = XLALSpinWeightedSphericalHarmonic( iota, phi_offset, -2, ll, mm  );
        /* Add the current mode to the intermediate timeseries */
        XLALSimRingdownAddSphericalMultipoleTD( &htilde, T0, deltaT, Mf, jf, eta, ll, mm, r, nonGRparams, Nsamples, Tstart, Prefactor );
    } /* END of FOR all MULTIPOLES */

    /* Extract hplus and hcross */
    COMPLEX16 h_val = 0.0;
    for ( UINT4 i=0 ; i<Nsamples ; i++ )
    {
        h_val = htilde->data->data[i];
        (*hplus)->data->data[i]  =        creal(h_val);
        (*hcross)->data->data[i] = -1.0 * cimag(h_val);
    } /* END of FOR all MULTIPOLES */

    /* Destroy the COMPLEX16TimeSeries hlmn object */
    if (htilde) XLALDestroyCOMPLEX16TimeSeries(htilde);

    return XLAL_SUCCESS;

}

/* ADD a QNM timeseries to an existing waveform time series with a prefactor (e.g. harmonic for inner-product) */
int XLALSimRingdownAddSphericalMultipoleTD(
     COMPLEX16TimeSeries **htilde,                /**< OUTPUT TD waveform mode lmn */
     const LIGOTimeGPS *T0,                       /**< start time of ringdown */
     REAL8 deltaT,                                /**< sampling interval (s) */
     REAL8 Mf,                                    /**< Final BH Mass (kg) */
     REAL8 jf,                                    /**< Final BH dimensionaless spin */
     REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
     UINT4 ll,                                     /**< Polar eigenvalue */
     INT4 mm,                                      /**< Azimuthal eigenvalue */
     REAL8 r,                                     /**< distance of source (m) */
     LALSimInspiralTestGRParam *nonGRparams,      /**< Testing GR params: fractional dfreq and dtau */
     UINT4 Nsamples,                              /**< waveform length */
     REAL8 Tstart,                                /**< starting time of waveform (10M at zero) */
     COMPLEX16 Prefactor                         /* The mode time series will be sclaed by this before being added to the htilde input. This can be a harmonic function or inner-product value. */
        ){

  /* Mode to be added to htilde */
  COMPLEX16TimeSeries *hllmm = XLALCreateCOMPLEX16TimeSeries("htilde: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
  if (! hllmm) XLAL_ERROR(XLAL_EFUNC);
  memset(hllmm->data->data, 0, Nsamples * sizeof(COMPLEX16));
  /* Compute the mode to be added */
  XLALSimRingdownGenerateSphericalMultipoleMMRDNS_time( &hllmm, T0, deltaT, Mf, jf, eta, ll, mm, r, nonGRparams, Nsamples, Tstart );

  /* Add the multipole to the reference timeseries including the prefactor */
  for ( UINT4 k=0 ; k<Nsamples ; k++ ) (*htilde)->data->data[k] += Prefactor * (hllmm->data->data[k]);

  /* Destroy the COMPLEX16TimeSeries hllmm object */
  if (hllmm) XLALDestroyCOMPLEX16TimeSeries(hllmm);

  return XLAL_SUCCESS;

}




/* Generate a single spherical harmonic multipole for MMRDNS */
int XLALSimRingdownGenerateSphericalMultipoleMMRDNS_time(
        COMPLEX16TimeSeries **htildeLM,            /**< OUTPUT TD waveform mode lmn */
        const LIGOTimeGPS *T0,                       /**< start time of ringdown */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        UINT4 ll,                                    /**< SPHERICAL Polar eigenvalue */
        INT4 mm,                                     /**< SPHERICAL Azimuthal eigenvalue */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams,      /**< Testing GR params: fractional dfreq and dtau */
        UINT4 Nsamples,                              /**< waveform length */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */
        ){

        /* allocate htildeLM */
        if ( ! *htildeLM ){
          *htildeLM = XLALCreateCOMPLEX16TimeSeries("htildeLM: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
          if (!(*htildeLM)) XLAL_ERROR(XLAL_EFUNC);
          memset((*htildeLM)->data->data, 0, Nsamples * sizeof(COMPLEX16));
        }

        /* Declare helper variables */
        COMPLEX16 Prefactor=0;
        UINT4 l=0,n=0; INT4 m=0;

        /* FOR all MODES in MMRDNS */
        for ( UINT4 k=0; k<XLALMMRDNS_NUM_MODES; k++ ) {
          /* Extract the mode indeces from the master list */
          l = XLALMMRDNS_MODES[k][0]; m = XLALMMRDNS_MODES[k][1]; n = XLALMMRDNS_MODES[k][2];
          /* NOTE that only QNMs with m==mm can appear in the ll,mm spherical multipole moment */
          if ( m == mm ) {
            /* In this setting the prefactor is the mixing coefficient between spherical and spheroidal harmonics */
            /* NOTE that XLALSphericalSpheroidalInnerProduct produces a result by numerical integration. This means many calls the the Spheroidal WHILE-LOOP, and many calls to the QNM tables. This results in very slow C code! (but the same thing is fast in python. Why?). XLALQNM_YSPROD uses a fit. */
            // Prefactor = XLALSphericalSpheroidalInnerProduct( jf, ll, mm, l, m, n, 256 );
            Prefactor = XLALQNM_YSPROD( jf, ll, mm, l, m, n );
            /* Add the current mode to the output timeseries */
            XLALSimRingdownAddSpheroidalModeTD( htildeLM, T0, deltaT, Mf, jf, eta, l, m, n, r, nonGRparams, Nsamples, Tstart, Prefactor );
          } /* END of IF m==mm */
        } /* END of FOR all MODES */

  return XLAL_SUCCESS;
}

/* ADD a QNM timeseries to an existing waveform time series with a prefactor (e.g. harmonic for inner-product) */
int XLALSimRingdownAddSpheroidalModeTD(
          COMPLEX16TimeSeries **htilde,                /**< OUTPUT TD waveform mode lmn */
          const LIGOTimeGPS *T0,                       /**< start time of ringdown */
          REAL8 deltaT,                                /**< sampling interval (s) */
          REAL8 Mf,                                    /**< Final BH Mass (kg) */
          REAL8 jf,                                    /**< Final BH dimensionaless spin */
          REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
          UINT4 l,                                     /**< Polar eigenvalue */
          INT4 m,                                      /**< Azimuthal eigenvalue */
          UINT4 n,                                     /**< Overtone Number */
          REAL8 r,                                     /**< distance of source (m) */
          LALSimInspiralTestGRParam *nonGRparams,      /**< Testing GR params: fractional dfreq and dtau */
          UINT4 Nsamples,                              /**< waveform length */
          REAL8 Tstart,                                /**< starting time of waveform (10M at zero) */
          COMPLEX16 Prefactor                         /* The mode time series will be sclaed by this before being added to the htilde input. This can be a harmonic function or inner-product value. */
        ){

  /* Mode to be added to htilde */
  COMPLEX16TimeSeries *hlmn;

  /* Extract NonGR Params */
  // NOTE: Recall that dfreq and dtau represent fractional changes, not shifts. This means this affects the handling of positive and negative m modes below.
  REAL8 dfreq=0,dtau=0;
  char *nonGRParamName = malloc(512*sizeof(char));
  sprintf(nonGRParamName,"dfreq%d%d%d",l,abs(m),n) ;
  if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
    dfreq = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
  sprintf(nonGRParamName,"dtau%d%d%d",l,abs(m),n) ;
  if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
    dtau = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;

  /* Compute the mode to be added */
  XLALSimRingdownGenerateSingleBareModeMMRDNS_time( &hlmn, T0, deltaT, Mf, jf, eta, l, m, n, r, dfreq, dtau, Nsamples, Tstart );

  /* Add the multipole to the reference timeseries including the prefactor */
  for ( UINT4 k=0 ; k<Nsamples ; k++ ) (*htilde)->data->data[k] += Prefactor * hlmn->data->data[k];

  /* Destroy the COMPLEX16TimeSeries hlmn object */
  if (hlmn) XLALDestroyCOMPLEX16TimeSeries(hlmn);
  /* Cleanup nonGRParamName */
  free(nonGRParamName);

  return XLAL_SUCCESS;

}







/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
/* Generate a full FD signal with angular dependence */
/* NOTE that this evaulates at an input frequency    */
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int XLALSimRingdownGenerateMMRDNS_freq(
        COMPLEX16FrequencySeries **hptilde,          /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde,          /**< OUTPUT FD h_x polarization */
        REAL8 df,
        REAL8 fMin,                                  /**< min frequency at which to evaluate */
        UINT4 Nsamples,
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< Inclination */
        REAL8 phi,                                   /**< Orbital phase */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */ ){

  /**/
  COMPLEX16FrequencySeries *hp220, *hp221, *hp330, *hp331, *hp440, *hp550, *hp210, *hp320, *hp430;
  COMPLEX16FrequencySeries *hc220, *hc221, *hc330, *hc331, *hc440, *hc550, *hc210, *hc320, *hc430;
  /**/
  COMPLEX16FrequencySeries *hp2m20, *hp2m21, *hp3m30, *hp3m31, *hp4m40, *hp5m50, *hp2m10, *hp3m20, *hp4m30;
  COMPLEX16FrequencySeries *hc2m20, *hc2m21, *hc3m30, *hc3m31, *hc4m40, *hc5m50, *hc2m10, *hc3m20, *hc4m30;
  /**/
  COMPLEX16 hp=0, hc=0;

  /**/
  LIGOTimeGPS tC = {0,0};
  /**/
  *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, df, &lalStrainUnit, Nsamples);
  if (!(*hptilde)) XLAL_ERROR(XLAL_EFUNC);
  memset((*hptilde)->data->data, 0, Nsamples * sizeof(COMPLEX16));
  XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
  /**/
  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, df, &lalStrainUnit, Nsamples);
  if (!(*hctilde)) XLAL_ERROR(XLAL_EFUNC);
  memset((*hctilde)->data->data, 0, Nsamples * sizeof(COMPLEX16));
  XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  //
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp220, &hc220, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 2,2,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp221, &hc221, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 2,2,1, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp210, &hc210, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 2,1,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp330, &hc330, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 3,3,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp331, &hc331, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 3,3,1, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp320, &hc320, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 3,2,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp440, &hc440, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 4,4,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp430, &hc430, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 4,3,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp550, &hc550, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 5,5,0, r, dfreq, dtau, Tstart );
  //
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp2m20, &hc2m20, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 2,-2,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp2m21, &hc2m21, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 2,-2,1, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp2m10, &hc2m10, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 2,-1,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp3m30, &hc3m30, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 3,-3,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp3m31, &hc3m31, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 3,-3,1, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp3m20, &hc3m20, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 3,-2,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp4m40, &hc4m40, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 4,-4,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp4m30, &hc4m30, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 4,-3,0, r, dfreq, dtau, Tstart );
  XLALSimRingdownGenerateSingleModeMMRDNS_freq( &hp5m50, &hc5m50, df, fMin, Nsamples, Mf, jf, eta, iota, phi, 5,-5,0, r, dfreq, dtau, Tstart );

  for ( UINT4 k=0; k<Nsamples; k++ ){
    //
    hp = 0;
    hp += hp220->data->data[k];  hp += hp2m20->data->data[k];
    hp += hp221->data->data[k];  hp += hp2m21->data->data[k];
    hp += hp210->data->data[k];  hp += hp2m10->data->data[k];
    hp += hp330->data->data[k];  hp += hp3m30->data->data[k];
    hp += hp331->data->data[k];  hp += hp3m31->data->data[k];
    hp += hp320->data->data[k];  hp += hp3m20->data->data[k];
    hp += hp440->data->data[k];  hp += hp4m40->data->data[k];
    hp += hp430->data->data[k];  hp += hp4m30->data->data[k];
    hp += hp550->data->data[k];  hp += hp5m50->data->data[k];
    //
    hc = 0;
    hc += hc220->data->data[k];  hc += hc2m20->data->data[k];
    hc += hc221->data->data[k];  hc += hc2m21->data->data[k];
    hc += hc210->data->data[k];  hc += hc2m10->data->data[k];
    hc += hc330->data->data[k];  hc += hc3m30->data->data[k];
    hc += hc331->data->data[k];  hc += hc3m31->data->data[k];
    hc += hc320->data->data[k];  hc += hc3m20->data->data[k];
    hc += hc440->data->data[k];  hc += hc4m40->data->data[k];
    hc += hc430->data->data[k];  hc += hc4m30->data->data[k];
    hc += hc550->data->data[k];  hc += hc5m50->data->data[k];
    //
    (*hptilde)->data->data[k]  =  hp;
    (*hctilde)->data->data[k]  =  hc;
  }

  /* Destroy the COMPLEX16TimeSeries objects */
  if (hp220) { XLALDestroyCOMPLEX16FrequencySeries(hp220); };   if (hp2m20) { XLALDestroyCOMPLEX16FrequencySeries(hp2m20); };
  if (hp221) { XLALDestroyCOMPLEX16FrequencySeries(hp221); };   if (hp2m21) { XLALDestroyCOMPLEX16FrequencySeries(hp2m21); };
  if (hp330) { XLALDestroyCOMPLEX16FrequencySeries(hp330); };   if (hp3m30) { XLALDestroyCOMPLEX16FrequencySeries(hp3m30); };
  if (hp331) { XLALDestroyCOMPLEX16FrequencySeries(hp331); };   if (hp3m31) { XLALDestroyCOMPLEX16FrequencySeries(hp3m31); };
  if (hp440) { XLALDestroyCOMPLEX16FrequencySeries(hp440); };   if (hp4m40) { XLALDestroyCOMPLEX16FrequencySeries(hp4m40); };
  if (hp550) { XLALDestroyCOMPLEX16FrequencySeries(hp550); };   if (hp5m50) { XLALDestroyCOMPLEX16FrequencySeries(hp5m50); };
  if (hp210) { XLALDestroyCOMPLEX16FrequencySeries(hp210); };   if (hp2m10) { XLALDestroyCOMPLEX16FrequencySeries(hp2m10); };
  if (hp320) { XLALDestroyCOMPLEX16FrequencySeries(hp320); };   if (hp3m20) { XLALDestroyCOMPLEX16FrequencySeries(hp3m20); };
  if (hp430) { XLALDestroyCOMPLEX16FrequencySeries(hp430); };   if (hp4m30) { XLALDestroyCOMPLEX16FrequencySeries(hp4m30); };
  /**/
  if (hc220) { XLALDestroyCOMPLEX16FrequencySeries(hc220); };   if (hc2m20) { XLALDestroyCOMPLEX16FrequencySeries(hc2m20); };
  if (hc221) { XLALDestroyCOMPLEX16FrequencySeries(hc221); };   if (hc2m21) { XLALDestroyCOMPLEX16FrequencySeries(hc2m21); };
  if (hc330) { XLALDestroyCOMPLEX16FrequencySeries(hc330); };   if (hc3m30) { XLALDestroyCOMPLEX16FrequencySeries(hc3m30); };
  if (hc331) { XLALDestroyCOMPLEX16FrequencySeries(hc331); };   if (hc3m31) { XLALDestroyCOMPLEX16FrequencySeries(hc3m31); };
  if (hc440) { XLALDestroyCOMPLEX16FrequencySeries(hc440); };   if (hc4m40) { XLALDestroyCOMPLEX16FrequencySeries(hc4m40); };
  if (hc550) { XLALDestroyCOMPLEX16FrequencySeries(hc550); };   if (hc5m50) { XLALDestroyCOMPLEX16FrequencySeries(hc5m50); };
  if (hc210) { XLALDestroyCOMPLEX16FrequencySeries(hc210); };   if (hc2m10) { XLALDestroyCOMPLEX16FrequencySeries(hc2m10); };
  if (hc320) { XLALDestroyCOMPLEX16FrequencySeries(hc320); };   if (hc3m20) { XLALDestroyCOMPLEX16FrequencySeries(hc3m20); };
  if (hc430) { XLALDestroyCOMPLEX16FrequencySeries(hc430); };   if (hc4m30) { XLALDestroyCOMPLEX16FrequencySeries(hc4m30); };
  /**/
  return XLAL_SUCCESS;

}


/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
/* Generate a single FD QNM WITH angular dependence */
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int XLALSimRingdownGenerateSingleModeMMRDNS_freq(
        COMPLEX16FrequencySeries **hptilde_lmn,      /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde_lmn,      /**< OUTPUT FD h_x polarization */
        REAL8 df,                                    /**< Frequency resolution (Hz) */
        REAL8 fMin,                                  /**< Start GW frequency (Hz) */
        UINT4 Nsamples,                              /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< Inclination */
        REAL8 phi,                                   /**< Orbital phase */
        UINT4 l,                                     /**< Polar eigenvalue */
        INT4 m,                                      /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */){

  /**/
  COMPLEX16 S = XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, phi);

  /**/
  XLALSimRingdownGenerateSingleBareModeMMRDNS_freq( hptilde_lmn, hctilde_lmn, df, fMin, Nsamples, Mf, jf, eta, l, m, n, r, dfreq, dtau, Tstart, S );

  // /* Evualte the frequency domain QNM at this frequency */
  // for ( UINT4 k=0; k<Nsamples; k++ ) {
  //   /**/
  //   (*hptilde_lmn)->data->data[k] =  (*hptilde_lmn)->data->data[k];
  //   (*hctilde_lmn)->data->data[k] =  (*hctilde_lmn)->data->data[k];
  //   // /**/
  //   // (*hptilde_lmn)->data->data[k] =  XLALSimSpheroidalHarmonicPlus(jf, l, m, n, iota, phi) * (*hptilde_lmn)->data->data[k];
  //   // (*hctilde_lmn)->data->data[k] = XLALSimSpheroidalHarmonicCross(jf, l, m, n, iota, phi) * (*hctilde_lmn)->data->data[k];
  // }

  /**/
  return XLAL_SUCCESS;

}


/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
/* Generate a single FD QNM without angular dependence */
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int XLALSimRingdownGenerateSingleBareModeMMRDNS_freq(
        COMPLEX16FrequencySeries **hptilde,          /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde,          /**< OUTPUT FD h_x polarization */
        REAL8 df,                                    /**< Frequency resolution (Hz) */
        REAL8 fMin,                                  /**< Start GW frequency (Hz) */
        UINT4 Nsamples,                              /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        UINT4 l,                                     /**< Polar eigenvalue */
        INT4 m,                                      /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        REAL8 Tstart,                                 /**< starting time of waveform (10M at zero) */
        COMPLEX16 Slmn){

  /* Declarations: M is the initial total mass, Mf is the remnant mass */
  REAL8 M             = XLALMf_to_M_nonspinning_UIB2016(eta, Mf);
  REAL8 M_sec         = M*LAL_MTSUN_SI/LAL_MSUN_SI;
  REAL8 Mf_sec        = Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
  REAL8 r_sec         = r/LAL_C_SI;
  COMPLEX16 A_lmn     = 0.0;
  COMPLEX16 A         = 0.0;
  //COMPLEX16 conj_A    = 0.0;
  COMPLEX16 CW        = 0.0;
  COMPLEX16 Omega_lmn = 0.0;
  COMPLEX16 chr = 0, chl = 0;//, conj_Omega_lmn = 0;
  REAL8 w, wmin, dw;

  /* Convert frequency input to angular freq from Hz */
  wmin = LAL_TWOPI * fMin;
  dw   = LAL_TWOPI * df;

  /* Evaluate QNM Amplitude model */
  A_lmn 	     = XLALMMRDNSAmplitudeOverOmegaSquared(eta, l, m, n);
  A            = A_lmn*M_sec/(r_sec);
  //conj_A       = conj(A);

  /* Get QNM frequency and allow for deviations */
  CW        = conj( XLALQNM_CW(jf,l,m,n) ) / Mf_sec ;
  Omega_lmn = creal(CW)*(1.0+dfreq) + I*cimag(CW)/(1.0+dtau);
  //conj_Omega_lmn = creal(CW)*(1.0+dfreq) - I*cimag(CW)/(1.0+dtau);

  /**/
  LIGOTimeGPS tC = {0,0};
  /**/
  *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde_lmn: FD waveform", &tC, 0.0, df, &lalStrainUnit, Nsamples);
  if (!(*hptilde)) XLAL_ERROR(XLAL_EFUNC);
  memset((*hptilde)->data->data, 0, Nsamples * sizeof(COMPLEX16));
  XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
  /**/
  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde_lmn: FD waveform", &tC, 0.0, df, &lalStrainUnit, Nsamples);
  if (!(*hctilde)) XLAL_ERROR(XLAL_EFUNC);
  memset((*hctilde)->data->data, 0, Nsamples * sizeof(COMPLEX16));
  XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  /* Evaluate the frequency domain QNM */
  for ( UINT4 k=0; k<Nsamples; k++ ) {
    w = wmin + k*dw;
    chr =       cexp( I * Omega_lmn*Tstart) * A / ( Omega_lmn - w );
    chl = conj( cexp( I * Omega_lmn*Tstart) * A / ( Omega_lmn + w ) );
    /**/
    (*hctilde)->data->data[k] = Slmn * ( chr + chl ) / 2.0 ;
    (*hptilde)->data->data[k] = Slmn * - I * ( chr - chl ) / 2.0 ;
  }

  return XLAL_SUCCESS;

}



/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
/* Generate a single FD QNM without angular dependence: NOT separated into + and x */
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int XLALSimRingdownEvalSinlgeModeMMRDNS_freq(
        COMPLEX16FrequencySeries **htilde,          /**< OUTPUT FD h_+ polarization */
        REAL8 df,                                    /**< Frequency resolution (Hz) */
        REAL8 fMin,                                  /**< Start GW frequency (Hz) */
        UINT4 Nsamples,                              /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        UINT4 l,                                     /**< Polar eigenvalue */
        INT4 m,                                      /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        REAL8 Tstart,                                 /**< starting time of waveform (10M at zero) */
        COMPLEX16 Slmn ){

  /* Declarations: M is the initial total mass, Mf is the remnant mass */
  REAL8 M             = XLALMf_to_M_nonspinning_UIB2016(eta, Mf);
  REAL8 M_sec         = M*LAL_MTSUN_SI/LAL_MSUN_SI;
  REAL8 Mf_sec        = Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
  REAL8 r_sec         = r/LAL_C_SI;
  COMPLEX16 A_lmn     = 0.0;
  COMPLEX16 A         = 0.0;
  //COMPLEX16 conj_A    = 0.0;
  COMPLEX16 CW        = 0.0;
  COMPLEX16 Omega_lmn = 0.0;
  REAL8 w, wmin, dw;

  /* Convert frequency input to angular freq from Hz */
  wmin = LAL_TWOPI * fMin;
  dw   = LAL_TWOPI * df;

  /* Evaluate QNM Amplitude model */
  A_lmn 	     = XLALMMRDNSAmplitudeOverOmegaSquared(eta, l, m, n);
  A            = A_lmn*M_sec/(r_sec);
  //conj_A       = conj(A);

  /* Get QNM frequency and allow for deviations */
  CW        = conj( XLALQNM_CW(jf,l,m,n) ) / Mf_sec ;
  Omega_lmn = creal(CW)*(1.0+dfreq) + I*cimag(CW)/(1.0+dtau);

  /**/
  LIGOTimeGPS tC = {0,0};
  /**/
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde_lmn: FD waveform", &tC, 0.0, df, &lalStrainUnit, Nsamples);
  if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);
  memset((*htilde)->data->data, 0, Nsamples * sizeof(COMPLEX16));
  XLALUnitMultiply(&(*htilde)->sampleUnits, &(*htilde)->sampleUnits, &lalSecondUnit);

  /* Evaluate the frequency domain QNM */
  for ( UINT4 k=0; k<Nsamples; k++ ) {
    w = wmin + k*dw;
    (*htilde)->data->data[k] = Slmn * cexp( I * Omega_lmn*Tstart) * A / ( Omega_lmn - w );
  }

  return XLAL_SUCCESS;

}
