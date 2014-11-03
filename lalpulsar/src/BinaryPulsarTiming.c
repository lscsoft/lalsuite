/*
*  Copyright (C) 2007 Jolien Creighton, Matt Pitkin
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
 * \author Matt Pitkin
 * \date 2006
 * \file
 * \ingroup pulsarTODO
 * \brief Functions to calculate binary system time delays and read TEMPO pulsar parameter files
 *
 * Functions for calculating the timing delay to a signal from a pulsar in a
 * binary system and reading pulsar parameters from TEMPO .par
 * files.
 * Models are taken from Taylor and Weisberg (1989) and use the
 * naming conventions therein and used by TEMPO .
 *
 * ### Prototypes ###
 *
 *
 * ### Description ###
 *
 * The main function computes the time delay of a signal from a pulsar in a
 * binary system due to doppler shifts and relativistic delays,
 * \f{equation}{
 * \Delta{}t = t_\textrm{Roemer} + t_\textrm{Shapiro} + t_\textrm{Einstein} + t_\textrm{
 * Abberation},
 * \f}
 * where \f$t_\textrm{Roemer}\f$ is the light travel time, \f$t_\textrm{Shapiro}\f$ is the
 * General relativistic time delay, \f$t_\textrm{Einstein}\f$ is the special
 * relativistic time delay, and \f$t_\textrm{Abberation}\f$ is the delay caused by the
 * pulsars' rotation. There are several models of the binary systems, described
 * in \cite TaylorWeisberg1989, of which the four most common are so far
 * implemented. The four models are the Blandford-Teukolsky model (BT)
 * \cite BlandfordTeukolsky1976, the low ellipticity model (ELL1)
 * \cite ChLangeetal2001, Damour-Deruelle model (DD) \cite DamourDeruelle1985,
 * and the main sequence system model (MSS) \cite Wex1998.
 * These four models all use the five main binary parameters: the longitude of
 * periastron \f$\omega_0\f$, the eccentricity of the orbit \f$e\f$, the orbital period
 * \f$P\f$, the time of periastron/or the time of ascension of the first node
 * \f$T_0\f$/\f$T_{\textrm{asc}}\f$, and the projected semi-major axis \f$a\sin{}i\f$. The are
 * also many other model dependent parameters. These routines closely follow
 * those used in the radio astronomy package TEMPO. A further model from TEMPO2
 * called T2 is also implemented in a basic form. The model is generally based
 * on the DD model, but will convert to ELL1 if the \c eps parameters are set.
 * At the moment this (T2) does not include multiple companions in the orbit,
 * but does encompass the DDS model. It also can include Kopeikin terms that
 * take account of the effect of the binary orbit on the parallax.
 *
 * Radio astronomers fit pulsar parameters using TEMPO which will output
 * the parameters in a <tt>.par</tt> file. The values allowed in this file can be
 * found in the TEMPO documentation. A function is included to extract these
 * parameters from the <tt>.par</tt> files and put them into a
 * \c BinaryPulsarParams structure, it will set any unused parameters to
 * zero or \c NULL. All parameters are in the units used by TEMPO with any
 * conversion to SI units occuring within the binary timing routines. A function
 * is also included which converts a string containing the right ascension or
 * declination in the format <tt>ddd/hh:mm:ss.s</tt> or <tt>ddd/hhmmss.s</tt>
 * (as is given in the <tt>.par</tt> file) into a \c REAL8 value in
 * radians.
 *
 * ### Notes ###
 *
 */

/* Matt Pitkin 29/04/04 */

#include <lal/BinaryPulsarTiming.h>

#include <string.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/SSBtimes.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define AULTSC 499.00478364 /* number of light seconds in AU (from tempo2.h) */

/**
 * XLAL function to compute the eccentric anomaly iteratively from Kelper's
 * equation.
 */
void XLALComputeEccentricAnomaly( REAL8 phase, REAL8 ecc, REAL8 *u){
  REAL8 du;

  *u = phase + ecc*sin(phase) / sqrt(1.0 - 2.*ecc*cos(phase) + ecc*ecc);
  do {
    du = (phase - (*u - ecc*sin(*u))) / (1.0 - ecc*cos(*u));
    (*u) += du;
  } while ( fabs(du) > 1.e-14 );
}


/**
 * XLAL function to compute Kopeikin terms that include the effect of
 * binary orbital parameters of parallax
 */
void XLALComputeKopeikinTerms( KopeikinTerms *kop,
                               BinaryPulsarParams *params,
                               BinaryPulsarInput *in ){
  REAL8 sini, cosi, tani;
  REAL8 sin_delta, cos_delta, sin_alpha, cos_alpha;
  REAL8 delta_i0, delta_j0;
  REAL8 sin_omega, cos_omega;
  REAL8 ca, sa, cd, sd;
  REAL8 delt;
  REAL8 posPulsar[3], velPulsar[3], psrPos[3];
  REAL8 tt0;
  REAL8 dpara; /* parallax */
  REAL8 x;

  sini = sin(params->kin);
  cosi = cos(params->kin);
  tani = sini / cosi;

  sin_omega = sin(params->kom);
  cos_omega = cos(params->kom);

  /* ki_dot is set in tempo2 function, but not used */
  /* ki_dot = -params.pmra * sin_omega + params.pmdec*cos_omega; */
  /* Equation 8 in Kopeikin 1996 */
  //  (*x) += ((*x)*ki_dot/tani)*tt0;
  /* Equation 9 in Kopeikin 1996 */
  //(*omz) += (pmra*cos_omega+pmdec*sin_omega)/sini*tt0;

  /* Now modify x and omega due to the annual-orbital parallax term
   * as described in Kopeikin 1995
   *
   * Require knowledge of the barycentric earth position vector - earth_ssb
   */

  /* get pulsar vector */
  ca = cos(params->ra);
  sa = sin(params->ra);
  cd = cos(params->dec);
  sd = sin(params->dec);

  posPulsar[0] = ca*cd;
  posPulsar[1] = sa*cd;
  posPulsar[2] = sd;

  velPulsar[0] = -params->pmra/cos(params->dec)*sa*cd - params->pmdec*ca*sd;
  velPulsar[1] = params->pmra/cos(params->dec)*ca*cd - params->pmdec*sa*sd;
  velPulsar[2] = params->pmdec*cd;

  delt = in->tb - params->posepoch;
  /* add proper motion onto the pulsar position */
  for( UINT4 i = 0; i < 3; i++ ) psrPos[i] = posPulsar[i] + delt*velPulsar[i];

  /* Obtain vector pointing at the pulsar */
  sin_delta = psrPos[2];
  cos_delta = cos(asin(sin_delta));
  sin_alpha = psrPos[1] / cos_delta;
  cos_alpha = psrPos[0] / cos_delta;

  /* Equation 15 in Kopeikin 1995 */
  delta_i0 = -in->earth.posNow[0]/AULTSC*sin_alpha +
    in->earth.posNow[1]/AULTSC*cos_alpha;
  /* Equation 16 in Kopeikin 1995 */
  delta_j0 = -in->earth.posNow[0]/AULTSC * sin_delta*cos_alpha -
    in->earth.posNow[1]/AULTSC * sin_delta*sin_alpha +
    in->earth.posNow[2]/AULTSC * cos_delta;

  dpara = params->px;
  x = params->x;

  /* xpr and ypr are set in tempo2 function, but not used */
  /* xpr = delta_i0*sin_omega - delta_j0*cos_omega;
  ypr = delta_i0*cos_omega + delta_j0*sin_omega; */

  /* Equations 18 and 19 in Kopeikin 1995 */
  if( params->daopset ){
    REAL8 daop = params->daop;

    kop->DK011 = - x / daop / sini*delta_i0*sin_omega;
    kop->DK012 = - x / daop / sini*delta_j0*cos_omega;
    kop->DK013 = - x / daop / sini*delta_i0*cos_omega;
    kop->DK014 = x / daop / sini*delta_j0*sin_omega;

    kop->DK021 = x / daop / tani*delta_i0*cos_omega;
    kop->DK022 = -x / daop / tani*delta_j0*sin_omega;
    kop->DK023 = x / daop / tani*delta_i0*sin_omega;
    kop->DK024 = x / daop / tani*delta_j0*cos_omega;
  }
  else{
    kop->DK011 = -x * dpara / sini*delta_i0*sin_omega;
    kop->DK012 = -x * dpara / sini*delta_j0*cos_omega;
    kop->DK013 = -x * dpara / sini*delta_i0*cos_omega;
    kop->DK014 = x * dpara / sini*delta_j0*sin_omega;

    kop->DK021 = x * dpara / tani*delta_i0*cos_omega;
    kop->DK022 = -x * dpara / tani*delta_j0*sin_omega;
    kop->DK023 = x * dpara / tani*delta_i0*sin_omega;
    kop->DK024 = x * dpara / tani*delta_j0*cos_omega;
  }

    if( params->T0 != 0. ) tt0 = in->tb - params->T0;
    else if( params->Tasc != 0. ) tt0 = in->tb - params->Tasc;
    else{
      XLALPrintError("%s: Neither T0 or Tasc is defined!\n", __func__);
      XLAL_ERROR_VOID( XLAL_EINVAL );
    }

    kop->DK031 = x * tt0 / sini*params->pmra*sin_omega;
    kop->DK032 = x * tt0 / sini*params->pmdec*cos_omega;
    kop->DK033 = x * tt0 / sini*params->pmra*cos_omega;
    kop->DK034 = -x * tt0 / sini*params->pmdec*sin_omega;

    kop->DK041 = x * tt0 / tani*params->pmra*cos_omega;
    kop->DK042 = -x * tt0 / tani*params->pmdec*sin_omega;
    kop->DK043 = -x * tt0 / tani*params->pmra*sin_omega;
    kop->DK044 = -x * tt0 / tani*params->pmdec*cos_omega;
}


/**
 * Calculate the binary system time delay using the pulsar parameters in
 * \c params
 */
void
LALBinaryPulsarDeltaT( LALStatus            *status,
                       BinaryPulsarOutput   *output,
                       BinaryPulsarInput    *input,
                       BinaryPulsarParams   *params ){
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* Check input arguments */
  ASSERT(input != (BinaryPulsarInput *)NULL, status,
  BINARYPULSARTIMINGH_ENULLINPUT, BINARYPULSARTIMINGH_MSGENULLINPUT);

  ASSERT(output != (BinaryPulsarOutput *)NULL, status,
  BINARYPULSARTIMINGH_ENULLOUTPUT, BINARYPULSARTIMINGH_MSGENULLOUTPUT);

  ASSERT(params != (BinaryPulsarParams *)NULL, status,
  BINARYPULSARTIMINGH_ENULLPARAMS, BINARYPULSARTIMINGH_MSGENULLPARAMS);

  ASSERT((!strcmp(params->model, "BT")) ||
         (!strcmp(params->model, "BT1P")) ||
         (!strcmp(params->model, "BT2P")) ||
         (!strcmp(params->model, "BTX")) ||
         (!strcmp(params->model, "ELL1")) ||
         (!strcmp(params->model, "DD")) ||
         (!strcmp(params->model, "DDS")) ||
         (!strcmp(params->model, "MSS")) ||
         (!strcmp(params->model, "T2")), status,
         BINARYPULSARTIMINGH_ENULLBINARYMODEL,
         BINARYPULSARTIMINGH_MSGNULLBINARYMODEL);

  XLALBinaryPulsarDeltaT( output, input, params );

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/**
 * XLAL function to compute the binary time delay
 */
void
XLALBinaryPulsarDeltaT( BinaryPulsarOutput   *output,
                        BinaryPulsarInput    *input,
                        BinaryPulsarParams   *params )
{
  REAL8 dt=0.; /* binary pulsar deltaT */
  REAL8 x, xdot;	/* x = asini/c */
  REAL8 w=0;  /* longitude of periastron */
  REAL8 e, edot;  /* eccentricity */
  REAL8 eps1, eps2;
  REAL8 eps1dot, eps2dot;
  REAL8 w0, wdot;
  REAL8 Pb, pbdot;
  REAL8 xpbdot;
  REAL8 T0, Tasc, tb=0.; /* time parameters */

  REAL8 s, r; /* Shapiro shape and range params */
  REAL8 lal_gamma; /* time dilation and grav redshift */
  REAL8 dr, dth;
  REAL8 shapmax; /* Shapiro max parameter for DDS model */
  REAL8 a0, b0;	/* abberation parameters */

  REAL8 m2;
  const REAL8 c3 = (REAL8)LAL_C_SI*(REAL8)LAL_C_SI*(REAL8)LAL_C_SI;

  CHAR *model = params->model;

  /* Check input arguments */
  if( input == (BinaryPulsarInput *)NULL ){
    XLAL_ERROR_VOID( BINARYPULSARTIMINGH_ENULLINPUT );
  }

  if( output == (BinaryPulsarOutput *)NULL ){
    XLAL_ERROR_VOID( BINARYPULSARTIMINGH_ENULLOUTPUT );
  }

  if( params == (BinaryPulsarParams *)NULL ){
    XLAL_ERROR_VOID( BINARYPULSARTIMINGH_ENULLPARAMS );
  }

  if((!strcmp(params->model, "BT")) &&
     (!strcmp(params->model, "BT1P")) &&
     (!strcmp(params->model, "BT2P")) &&
     (!strcmp(params->model, "BTX")) &&
     (!strcmp(params->model, "ELL1")) &&
     (!strcmp(params->model, "DD")) &&
     (!strcmp(params->model, "DDS")) &&
     (!strcmp(params->model, "MSS")) &&
     (!strcmp(params->model, "T2"))){
    XLAL_ERROR_VOID( BINARYPULSARTIMINGH_ENULLBINARYMODEL );
  }

  /* convert certain params to SI units */
  w0 = params->w0;
  wdot = params->wdot; /* wdot in rads/s */

  Pb = params->Pb; /* period in secs */
  pbdot = params->Pbdot;

  T0 = params->T0; /* these should be in TDB in seconds */
  Tasc = params->Tasc;

  e = params->e;
  edot = params->edot;
  eps1 = params->eps1;
  eps2 = params->eps2;
  eps1dot = params->eps1dot;
  eps2dot = params->eps2dot;

  x = params->x;
  xdot = params->xdot;
  xpbdot = params->xpbdot;

  lal_gamma = params->gamma;
  s = params->s; /* sin i */
  dr = params->dr;
  dth = params->dth;
  shapmax = params->shapmax;

  a0 = params->a0;
  b0 = params->b0;

  m2 = params->m2;

  /* Shapiro range parameter r defined as Gm2/c^3 (secs) */
  r = LAL_G_SI*m2/c3;

  /* if T0 is not defined, but Tasc is */
  if(T0 == 0.0 && Tasc != 0.0 && eps1 == 0.0 && eps2 == 0.0){
    REAL8 fe, uasc, Dt; /* see TEMPO tasc2t0.f */

    fe = sqrt((1.0 - e)/(1.0 + e));
    uasc = 2.0*atan(fe*tan(w0/2.0));
    Dt = (Pb/LAL_TWOPI)*(uasc-e*sin(uasc));

    T0 = Tasc + Dt;
  }

  /* set time at which to calculate the binary time delay */
  tb = input->tb;

  /* for BT, BT1P and BT2P models (and BTX model, but only for one orbit) */
  if(strstr(model, "BT") != NULL){
    REAL8 tt0;
    REAL8 orbits=0.;
    INT4 norbits=0.;
    REAL8 phase; /* same as mean anomaly */
    REAL8 u = 0.0; /* eccentric anomaly */

    INT4 nplanets=1; /* number of orbiting bodies in system */
    INT4 i=1, j=1;
    REAL8 fac=1.; /* factor in front of fb coefficients */

    REAL8 su = 0., cu = 0.;
    REAL8 sw = 0., cw = 0.;

    /* work out number of orbits i.e. have we got a BT1P or BT2P model */
    if( !strcmp(model, "BT1P") ) nplanets = 2;
    if( !strcmp(model, "BT2P") ) nplanets = 3;

    for ( i=1 ; i < nplanets+1 ; i++){

      /* set some vars for bnrybt.f (TEMPO) method */
      /*REAL8 tt;
      REAL8 som;
      REAL8 com;
      REAL8 alpha, beta;*/
      /*REAL8 q, r, s;*/

      /*fprintf(stderr, "You are using the Blandford-Teukolsky (BT) binary
        model.\n");*/

      if(i==2){
        T0 = params->T02;
        w0 = params->w02;
        x = params->x2;
        e = params->e2;
        Pb = params->Pb2;
      }
      else if(i==3){
        T0 = params->T03;
        w0 = params->w03;
        x = params->x3;
        e = params->e3;
        Pb = params->Pb3;
      }

      tt0 = tb - T0;

      /* only do relativistic corrections for first orbit */
      if(i==1){
        x = x + xdot*tt0;
        e = e + edot*tt0;
        w = w0 + wdot*tt0; /* calculate w */

        if( !strcmp(model, "BTX") ){
          fac = 1.;
          for ( j=1 ; j < params->nfb + 1; j++){
            fac /= (REAL8)j;
            orbits += fac*params->fb[j-1]*pow(tt0,j);
          }
        }
        else{
          orbits = tt0/Pb - 0.5*(pbdot+xpbdot)*(tt0/Pb)*(tt0/Pb);
        }
      }
      else{
        orbits = tt0/Pb;
      }

      norbits = (INT4)orbits;

      if(orbits < 0.) norbits--;

      phase = LAL_TWOPI*(orbits - (REAL8)norbits); /* called phase in TEMPO */

      /* compute eccentric anomaly */
      XLALComputeEccentricAnomaly( phase, e, &u );

      su = sin(u);
      cu = cos(u);

      /*fprintf(stderr, "Eccentric anomaly = %f, phase = %f.\n", u, phase);*/
      sw = sin(w);
      cw = cos(w);

      /* see eq 5 of Taylor and Weisberg (1989) */
      /**********************************************************/
      if( !strcmp(model, "BTX") ){
        /* dt += (x*sin(w)*(cos(u)-e) + (x*cos(w)*sqrt(1.0-e*e) +
          lal_gamma)*sin(u))*(1.0 - params->fb[0]*(x*cos(w)*sqrt(1.0 -
          e*e)*cos(u) - x*sin(w)*sin(u))/(1.0 - e*cos(u))); */
        dt += (x*sw*(cu-e) + (x*cw*sqrt(1.0-e*e) +
          lal_gamma)*su)*(1.0 - LAL_TWOPI*params->fb[0]*(x*cw*sqrt(1.0 -
          e*e)*cu - x*sw*su)/(1.0 - e*cu));
      }
      else{
        /* dt += (x*sin(w)*(cos(u)-e) + (x*cos(w)*sqrt(1.0-e*e) +
          lal_gamma)*sin(u))*(1.0 - (LAL_TWOPI/Pb)*(x*cos(w)*sqrt(1.0 -
          e*e)*cos(u) - x*sin(w)*sin(u))/(1.0 - e*cos(u))); */
        dt += (x*sw*(cu-e) + (x*cw*sqrt(1.0-e*e) +
          lal_gamma)*su)*(1.0 - (LAL_TWOPI/Pb)*(x*cw*sqrt(1.0 -
          e*e)*cu - x*sw*su)/(1.0 - e*cu));
      }
    /**********************************************************/
    }

    /* use method from Taylor etal 1976 ApJ Lett and used in bnrybt.f */
    /**********************************************************/
    /*tt = 1.0-e*e;
    som = sin(w);
    com = cos(w);
    alpha = x*som;
    beta = x*com*sqrt(tt);
    q = alpha*(cos(u)-e) + (beta+lal_gamma)*sin(u);
    r = -alpha*sin(u) + beta*cos(u);
    s = 1.0/(1.0-e*cos(u));
    dt = -(-q+(LAL_TWOPI/Pb)*q*r*s);*/
    /**********************************************************/
    /* There appears to be NO difference between either method */

    output->deltaT = -dt;
  }

  /* for ELL1 model (low eccentricity orbits so use eps1 and eps2) */
  /* see Appendix A, Ch. Lange etal, MNRAS (2001) (also accept T2 model if
   eps values are set - this will include Kopeikin terms if necessary) */
  if( !strcmp(model, "ELL1") || (!strcmp(model, "T2") && eps1 != 0. ) ){
    REAL8 nb = LAL_TWOPI/Pb;
    REAL8 tt0;
    REAL8 w_int; /* omega internal to this model */
    REAL8 orbits, phase;
    INT4 norbits;
    REAL8 e1, e2, ecc;
    REAL8 DRE, DREp, DREpp; /* Roemer and Einstein delays (cf DD) */
    REAL8 dlogbr;
    REAL8 DS, DA; /* Shapiro delay and Abberation delay terms */
    REAL8 Dbb;
    REAL8 DAOP, DSR; /* Kopeikin delay terms */

    KopeikinTerms kt;
    REAL8 Ck, Sk;

    REAL8 sp = 0., cp = 0., s2p = 0., c2p = 0.;

    /* fprintf(stderr, "You are using the ELL1 low eccentricity orbit model.\n");*/

    /*********************************************************/
    /* CORRECT CODE (as in TEMPO bnryell1.f) FROM HERE       */

    ecc = sqrt(eps1*eps1 + eps2*eps2);

    /* if Tasc is not defined convert T0 */
    if(Tasc == 0.0 && T0 != 0.0){
      REAL8 fe, uasc, Dt; /* see TEMPO tasc2t0.f */

      fe = sqrt((1.0 - ecc)/(1.0 + ecc));
      uasc = 2.0*atan(fe*tan(w0/2.0));
      Dt = (Pb/LAL_TWOPI)*(uasc-ecc*sin(uasc));

      /* rearrange from what's in tasc2t0.f */
      Tasc = T0 - Dt;
    }

    tt0 = tb - Tasc;

    orbits = tt0/Pb - 0.5*(pbdot+xpbdot)*(tt0/Pb)*(tt0/Pb);
    norbits = (INT4)orbits;
    if(orbits < 0.0) norbits--;

    phase=LAL_TWOPI*(orbits - (REAL8)norbits);

    x = x + xdot*tt0;

    /* depending on whether we have eps derivs or w time derivs calculate e1 and e2 accordingly */
    if(params->nEll == 0){
      e1 = eps1 + eps1dot*tt0;
      e2 = eps2 + eps2dot*tt0;
    }
    else{
      ecc = sqrt(eps1*eps1 + eps2*eps2);
      ecc += edot*tt0;
      w_int = atan2(eps1, eps2);
      w_int = w_int + wdot*tt0;

      e1 = ecc*sin(w_int);
      e2 = ecc*cos(w_int);
    }

    //sin_cos_LUT(&sp, &cp, phase);
    //sin_cos_LUT(&s2p, &c2p, 2.*phase);
    sp = sin(phase);
    cp = cos(phase);
    s2p = sin(2.*phase);
    c2p = cos(2.*phase);

    /* this timing delay (Roemer + Einstein) should be most important in most cases */
    /* DRE = x*(sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase)));
    DREp = x*cos(phase);
    DREpp = -x*sin(phase); */
    DRE = x*(sp-0.5*(e1*c2p-e2*s2p));
    DREp = x*cp;
    DREpp = -x*sp;

    /* these params will normally be negligable */
    dlogbr = log(1.0-s*sp);
    DS = -2.0*r*dlogbr;
    DA = a0*sp + b0*cp;

    /* compute Kopeikin terms */
    if( params->kinset && params->komset && ( params->pmra != 0. ||
        params->pmdec != 0. ) ){
      XLALComputeKopeikinTerms( &kt, params, input );

      Ck = sp - 0.5*(e1*c2p - e2*s2p);
      Sk = cp + 0.5*(e2*c2p + e1*s2p);

      DAOP = (kt.DK011 + kt.DK012)*Ck - (kt.DK021 + kt.DK022)*Sk;
      DSR = (kt.DK031 + kt.DK032)*Ck + (kt.DK041 + kt.DK042)*Sk;
    }
    else{
      DAOP = 0.;
      DSR = 0.;
    }

    Dbb = DRE*(1.0-nb*DREp+(nb*DREp)*(nb*DREp) + 0.5*nb*nb*DRE*DREpp) + DS + DA
      + DAOP + DSR;

    output->deltaT = -Dbb;
    /********************************************************/
  }

  /* for DD model - code partly adapted from TEMPO bnrydd.f */
  /* also used for MSS model (Wex 1998) - main sequence star orbit - this only has two lines
different than DD model - TEMPO bnrymss.f */
  /* also DDS model and (partial) T2 model (if EPS params not set) from TEMPO2
T2model.C */
  if( !strcmp(model, "DD") || !strcmp(model, "MSS") || !strcmp(model, "DDS") ||
      (!strcmp(model, "T2") && eps1 == 0.) ){
    REAL8 u;        /* new eccentric anomaly */
    REAL8 Ae;       /* eccentricity parameter */
    REAL8 DRE;      /* Roemer delay + Einstein delay */
    REAL8 DREp, DREpp; /* see DD eqs 48 - 50 */
    REAL8 DS;       /* Shapiro delay */
    REAL8 DA;       /* aberation caused by pulsar rotation delay */
    REAL8 DAOP, DSR; /* Kopeikin term delays */

    REAL8 tt0;
    /* various variable use during calculation */
    REAL8 er, eth, an, k;
    REAL8 orbits, phase;
    INT4 norbits;
    REAL8 onemecu, cae, sae;
    REAL8 alpha, beta, bg;
    REAL8 anhat, sqr1me2, cume, brace, dlogbr;
    REAL8 Dbb;    /* Delta barbar in DD eq 52 */

    REAL8 xi; /* parameter for MSS model - the only other one needed */
    REAL8 sdds = 0.; /* parameter for DDS model */

    REAL8 su = 0., cu = 0.;
    REAL8 sw = 0., cw = 0.;

    KopeikinTerms kt;
    REAL8 Ck, Sk;

    /* fprintf(stderr, "You are using the Damour-Deruelle (DD) binary model.\n");*/

    /* part of code adapted from TEMPO bnrydd.f */
    an = LAL_TWOPI/Pb;
    k = wdot/an;
    xi = xdot/an; /* MSS parameter */

    tt0 = tb - T0;
    /* x = x + xdot*tt0; */
    e = e + edot*tt0;
    er = e*(1.0+dr);
    eth = e*(1.0+dth);

    orbits = (tt0/Pb) - 0.5*(pbdot+xpbdot)*(tt0/Pb)*(tt0/Pb);
    norbits = (INT4)orbits;

    if(orbits < 0.0) norbits--;

    phase = LAL_TWOPI*(orbits - (REAL8)norbits);

    /* compute eccentric anomaly */
    XLALComputeEccentricAnomaly( phase, e, &u );

    su = sin(u);
    cu = cos(u);

    /* compute Ae as in TEMPO bnrydd.f */
    onemecu = 1.0 - e*cu;
    cae = (cu - e)/onemecu;
    sae = sqrt(1.0 - e*e)*su/onemecu;

    Ae = atan2(sae,cae);

    if(Ae < 0.0)
      Ae = Ae + LAL_TWOPI;

    Ae = LAL_TWOPI*orbits + Ae - phase;

    w = w0 + k*Ae; /* add corrections to omega */ /* MSS also uses (om2dot, but not defined) */

    /* small difference between MSS and DD */
    if( !strcmp(model, "MSS") ){
      x = x + xi*Ae; /* in bnrymss.f they also include a second time derivative of x (x2dot), but
this isn't defined for either of the two pulsars currently using this model */
    }
    else
      x = x + xdot*tt0;

    /* now compute time delays as in DD eqs 46 - 52 */

    /* calculate Einstein and Roemer delay */
    sw = sin(w);
    cw = cos(w);
    alpha = x*sw;
    beta = x*sqrt(1.0-eth*eth)*cw;
    bg = beta + lal_gamma;
    DRE = alpha*(cu-er)+bg*su;
    DREp = -alpha*su + bg*cu;
    DREpp = -alpha*cu - bg*su;
    anhat = an/onemecu;

    /* calculate Shapiro and abberation delays DD eqs 26, 27 */
    sqr1me2 = sqrt(1.0-e*e);
    cume = cu-e;
    if( !strcmp(model, "DDS") ){
      sdds = 1. - exp(-1.*shapmax);
      brace = onemecu-sdds*(sw*cume + sqr1me2*cw*su);
    }
    else brace = onemecu-s*(sw*cume + sqr1me2*cw*su);
    dlogbr = log(brace);
    DS = -2.0*r*dlogbr;

    /* this abberation delay is prob fairly small */
    DA = a0*(sin(w+Ae)+e*sw) + b0*(cos(w+Ae)+e*cw);

    /* compute Kopeikin terms */
    if( params->kinset && params->komset && ( params->pmra != 0. ||
        params->pmdec != 0. ) ){
      XLALComputeKopeikinTerms( &kt, params, input );

      Ck = cw*(cu-er) - sqrt(1.-eth*eth)*sw*su;
      Sk = sw*(cu-er) + sqrt(1.-eth*eth)*cw*su;

      DAOP = (kt.DK011 + kt.DK012)*Ck - (kt.DK021 + kt.DK022)*Sk;
      DSR = (kt.DK031 + kt.DK032)*Ck + (kt.DK041 + kt.DK042)*Sk;
    }
    else{
      DAOP = 0.;
      DSR = 0.;
    }

    /* timing difference */
    Dbb = DRE*(1.0 - anhat*DREp+anhat*anhat*DREp*DREp + 0.5*anhat*anhat*DRE*DREpp -
          0.5*e*su*anhat*anhat*DRE*DREp/onemecu) + DS + DA + DAOP + DSR;

    output->deltaT = -Dbb;
  }

  /* for DDGR model */

  /* for Epstein-Haugan (EH) model - see Haugan, ApJ (1985) eqs 69 and 71 */

  /* check that the returned value is not a NaN */
  if( isnan(output->deltaT) ){
    XLAL_ERROR_VOID( BINARYPULSARTIMINGH_ENAN );
  }
}
