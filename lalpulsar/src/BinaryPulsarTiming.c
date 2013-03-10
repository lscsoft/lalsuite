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
   Functions for calculating the timing delay to a signal from a pulsar in a
   binary system and reading pulsar parameters from TEMPO .par
   files.
   Models are taken from Taylor and Weisberg (1989) and use the
   naming conventions therein and used by TEMPO .

   \heading{Prototypes}



   \heading{Description}

   The main function computes the time delay of a signal from a pulsar in a
   binary system due to doppler shifts and relativistic delays,
   \f{equation}{
   \Delta{}t = t_\textrm{Roemer} + t_\textrm{Shapiro} + t_\textrm{Einstein} + t_\textrm{
   Abberation},
   \f}
   where \f$t_\textrm{Roemer}\f$ is the light travel time, \f$t_\textrm{Shapiro}\f$ is the
   General relativistic time delay, \f$t_\textrm{Einstein}\f$ is the special
   relativistic time delay, and \f$t_\textrm{Abberation}\f$ is the delay caused by the
   pulsars' rotation. There are several models of the binary systems, described
   in [\ref TaylorWeisberg1989], of which the four most common are so far
   implemented. The four models are the Blandford-Teukolsky model (BT)
   [\ref BlandfordTeukolsky1976], the low ellipticity model (ELL1)
   [\ref ChLangeetal2001], Damour-Deruelle model (DD) [\ref DamourDeruelle1985],
   and the main sequence system model (MSS) [\ref Wex1998].
   These four models all use the five main binary parameters: the longitude of
   periastron \f$\omega_0\f$, the eccentricity of the orbit \f$e\f$, the orbital period
   \f$P\f$, the time of periastron/or the time of ascension of the first node
   \f$T_0\f$/\f$T_{\textrm{asc}}\f$, and the projected semi-major axis \f$a\sin{}i\f$. The are
   also many other model dependent parameters. These routines closely follow
   those used in the radio astronomy package TEMPO. A further model from TEMPO2
   called T2 is also implemented in a basic form. The model is generally based
   on the DD model, but will convert to ELL1 if the \c eps parameters are set.
   At the moment this (T2) does not include multiple companions in the orbit,
   but does encompass the DDS model. It also can include Kopeikin terms that
   take account of the effect of the binary orbit on the parallax.

   Radio astronomers fit pulsar parameters using TEMPO which will output
   the parameters in a <tt>.par</tt> file. The values allowed in this file can be
   found in the TEMPO documentation. A function is included to extract these
   parameters from the <tt>.par</tt> files and put them into a
   \c BinaryPulsarParams structure, it will set any unused parameters to
   zero or \c NULL. All parameters are in the units used by TEMPO with any
   conversion to SI units occuring within the binary timing routines. A function
   is also included which converts a string containing the right ascension or
   declination in the format <tt>ddd/hh:mm:ss.s</tt> or <tt>ddd/hhmmss.s</tt>
   (as is given in the <tt>.par</tt> file) into a \c REAL8 value in
   radians.

   \heading{Notes}

*/

/* Matt Pitkin 29/04/04 */

#include <lal/BinaryPulsarTiming.h>

#include <string.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/ComputeFstat.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define DAYSTOSECS 86400.0 /* number of seconds in a day */
#define AULTSC 499.00478364 /* number of light seconds in AU (from tempo2.h) */

/** XLAL function to compute the eccentric anomaly iteratively from Kelper's
 * equation. */
void XLALComputeEccentricAnomaly( REAL8 phase, REAL8 ecc, REAL8 *u){
  REAL8 du;

  *u = phase + ecc*sin(phase) / sqrt(1.0 - 2.*ecc*cos(phase) + ecc*ecc);
  do {
    du = (phase - (*u - ecc*sin(*u))) / (1.0 - ecc*cos(*u));
    (*u) += du;
  } while ( fabs(du) > 1.e-14 );
}

/** XLAL function to compute Kopeikin terms that include the effect of
  * binary orbital parameters of parallax */
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
      XLAL_ERROR_VOID( XLAL_EFUNC );
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


/** Calculate the binary system time delay using the pulsar parameters in
 *  \c params
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


/** XLAL function to compute the binary time delay
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
  REAL8 c3 = (REAL8)LAL_C_SI*(REAL8)LAL_C_SI*(REAL8)LAL_C_SI;

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
  w0 = params->w0*LAL_PI_180; /* convert w to rads from degs */
  wdot = params->wdot*LAL_PI_180/(365.25*DAYSTOSECS); /* convert wdot to rads/s from degs/yr */

  Pb = params->Pb*DAYSTOSECS; /* covert period from days to secs */
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

  m2 = params->m2*LAL_MSUN_SI;

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
        w0 = params->w02*LAL_PI_180;
        x = params->x2;
        e = params->e2;
        Pb = params->Pb2*DAYSTOSECS;
      }
      else if(i==3){
        T0 = params->T03;
        w0 = params->w03*LAL_PI_180;
        x = params->x3;
        e = params->e3;
        Pb = params->Pb3*DAYSTOSECS;
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


void
LALReadTEMPOParFile(  LALStatus *status,
                      BinaryPulsarParams *output,
                      CHAR      *pulsarAndPath )
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(output != (BinaryPulsarParams *)NULL, status,
  BINARYPULSARTIMINGH_ENULLOUTPUT, BINARYPULSARTIMINGH_MSGENULLOUTPUT);

  XLALReadTEMPOParFile( output, pulsarAndPath );

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* NOTE: Convert this function to be more like readParfile.C in TEMPO2 - read
 * in a line at a time using fgets and make each parameter a structure */
void
XLALReadTEMPOParFile( BinaryPulsarParams *output,
                      CHAR      *pulsarAndPath )
{
  FILE *fp=NULL;
  CHAR val[500][40]; /* string array to hold all the read in values
                        500 strings of max 40 characters is enough */
  INT4 i=0, j=1, k;
  int UNUSED c;

  if( output == (BinaryPulsarParams *)NULL ){
    XLAL_ERROR_VOID( XLAL_EFAULT );
  }

  output->name = NULL;
  output->jname = NULL;
  output->bname = NULL;

  output->model = NULL; /* set binary model to null - in case not a binary */

  /* set all output params to zero*/
  output->e=0.0;      /* orbital eccentricity */
  output->Pb=0.0;     /* orbital period (days) */
  output->w0=0.0;     /* longitude of periastron (deg) */
  output->x=0.0;      /* projected semi-major axis/speed of light (light secs) */
  output->T0=0.0;     /* time of orbital periastron as measured in TDB (MJD) */

  output->e2=0.0;
  output->Pb2=0.0;
  output->w02=0.0;
  output->x2=0.0;
  output->T02=0.0;

  output->e3=0.0;
  output->Pb3=0.0;
  output->w03=0.0;
  output->x3=0.0;
  output->T03=0.0;

  output->xpbdot=0.0;  /* (10^-12) */

  output->eps1=0.0;       /* e*sin(w) */
  output->eps2=0.0;       /* e*cos(w) */
  output->eps1dot=0.0;
  output->eps2dot=0.0;
  output->Tasc=0.0;   /* time of the ascending node (used rather than T0) */

  output->fb = NULL;
  output->fbErr = NULL;
  output->nfb=0;

  output->wdot=0.0;   /* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
  output->gamma=0.0;  /* gravitational redshift and time dilation parameter (s)*/
  output->Pbdot=0.0;  /* rate of change of Pb (dimensionless 10^-12) */
  output->xdot=0.0;   /* rate of change of x(=asini/c) - optional (10^-12)*/
  output->edot=0.0;   /* rate of change of e (10^-12)*/

  output->s=0.0;      /* Shapiro 'shape' parameter sin i */
  output->sstr=NULL;

  output->shapmax=0.;

  /*output.r=0.0; Shapiro 'range' parameter */
  output->dr=0.0;
  output->dth=0.0;    /* (10^-6) */
  output->a0=0.0;
  output->b0=0.0; /* abberation delay parameters */

  output->M=0.0;     /* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
  output->m2=0.0;    /* companion mass */

  output->f0=0.0;
  output->f1=0.0;
  output->f2=0.0;
  output->f3=0.0;
  output->f4=0.0;
  output->f5=0.0;
  output->f6=0.0;
  output->f7=0.0;
  output->f8=0.0;
  output->f9=0.0;

  output->waveSin = NULL;
  output->waveCos = NULL;
  output->wave_om = 0.0;
  output->waveepoch = 0.0;
  output->nwaves = 0;

  output->ra=0.0;
  output->dec=0.0;
  output->pmra=0.0;
  output->pmdec=0.0;

  output->px=0.;    /* parallax (mas) */
  output->dist=0.;  /* distance (kpc) */

  output->DM=0.;    /* dispersion measure */
  output->DM1=0.;   /* first derivative of dispersion measure */

  output->daop=0.;
  output->daopset=0;
  output->kin=0.;
  output->kinset=0;
  output->kom=0.;
  output->komset=0;

  /* set all errors on params to zero */
  output->raErr=0.0;
  output->decErr=0.0;
  output->pmraErr=0.0;
  output->pmdecErr=0.0;

  output->posepoch=0.0;
  output->pepoch=0.0;

  output->posepochErr=0.0;
  output->pepochErr=0.0;

  output->xpbdotErr=0.0;  /* (10^-12) */

  output->eps1Err=0.0;        /* e*sin(w) */
  output->eps2Err=0.0;        /* e*cos(w) */
  output->eps1dotErr=0.0;
  output->eps2dotErr=0.0;
  output->TascErr=0.0;    /* time of the ascending node (used rather than T0) */

  output->wdotErr=0.0;   /* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
  output->gammaErr=0.0;  /* gravitational redshift and time dilation parameter (s)*/
  output->PbdotErr=0.0;  /* rate of change of Pb (dimensionless 10^-12) */
  output->xdotErr=0.0;   /* rate of change of x(=asini/c) - optional (10^-12)*/
  output->edotErr=0.0;   /* rate of change of e (10^-12)*/

  output->sErr=0.0;     /* Shapiro 'shape' parameter sin i */
  output->shapmaxErr=0.;

  /*output->rErr=0.0;  Shapiro 'range' parameter */
  output->drErr=0.0;
  output->dthErr=0.0;   /* (10^-6) */
  output->a0Err=0.0;
  output->b0Err=0.0;    /* abberation delay parameters */

  output->MErr=0.0;     /* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
  output->m2Err=0.0;    /* companion mass */

  output->f0Err=0.0;
  output->f1Err=0.0;
  output->f2Err=0.0;
  output->f3Err=0.0;
  output->f4Err=0.0;
  output->f5Err=0.0;
  output->f6Err=0.0;
  output->f7Err=0.0;
  output->f8Err=0.0;
  output->f9Err=0.0;

  output->eErr =0.0;
  output->w0Err=0.0;
  output->PbErr=0.0;
  output->xErr=0.0;
  output->T0Err=0.0;

  output->e2Err =0.0;
  output->w02Err=0.0;
  output->Pb2Err=0.0;
  output->x2Err=0.0;
  output->T02Err=0.0;

  output->e3Err =0.0;
  output->w03Err=0.0;
  output->Pb3Err=0.0;
  output->x3Err=0.0;
  output->T03Err=0.0;

  output->pxErr=0.;
  output->distErr=0.;

  output->DMErr=0.;
  output->DM1Err=0.;

  output->h0=0.;
  output->cosiota=0.;
  output->psi=0.;
  output->phi0=0.;
  output->Aplus=0.;
  output->Across=0.;
  output->I21=0.;
  output->I31=0.;
  output->r=0.;
  output->lambda=0.;
  output->costheta=0.;

  output->h0Err=0.;
  output->cosiotaErr=0.;
  output->psiErr=0.;
  output->phi0Err=0.;
  output->AplusErr=0.;
  output->AcrossErr=0.;
  output->I21Err=0.;
  output->I31Err=0.;
  output->rErr=0.;
  output->lambdaErr=0.;
  output->costhetaErr=0.;

  output->wave_omErr = 0.0;

  output->units = NULL;
  output->ephem = NULL;

  if((fp = fopen(pulsarAndPath, "r")) == NULL){
    XLALPrintError("Error... Cannot open .par file %s\n", pulsarAndPath);
    XLAL_ERROR_VOID( XLAL_EIO );
  }

  /* read all the pulsar data into the string array */
  while(!feof(fp)){
    /* make sure val[i] is clear first */
    sprintf(val[i], "%s", "");

    c = fscanf(fp, "%s", val[i]);

    /* if line starts with a '#' then skip to end of line */
    if( val[i][0] == '#' ){
       /* skip to the end of the line */
      c = fscanf(fp, "%*[^\n]");
      if ( feof(fp) ) break;
      continue;
    }

    i++;
  }

  k=i; /* k is the end number */
  i=0; /* reset i */

  /* set pulsar values for output */
  /* in .par files first column will param name, second will be param value,
     if third is defined it will be an integer to tell TEMPO whether to fit
     the param or not (don't need this), fourth will be the error on the
     param (in same units as the param) */

  /* convert all epochs given in MJD in .par files to secs in TDB  */
  while(1){
    j=i;
    if(!strcmp(val[i], "NAME") || !strcmp(val[i], "name")){
      output->name = XLALStringDuplicate(val[i+1]);
      j++;
    }
    else if(!strcmp(val[i], "PSRJ") || !strcmp(val[i], "psrj") ){
      output->jname = XLALStringDuplicate(val[i+1]);
      j++;
    }
    else if(!strcmp(val[i], "PSRB") || !strcmp(val[i], "psrb") ){
      output->bname = XLALStringDuplicate(val[i+1]);
      j++;
    }
    else if(!strcmp(val[i],"ra") || !strcmp(val[i],"RA") || !strcmp(val[i],"RAJ")){
      /* this can be in form hh:mm:ss.ss or hhmmss.ss */
      output->ra = XLALhmsToRads(val[i+1]);
      j++;

      /* only try to get error if one exists */
      if(atoi(val[i+2])==1 && i+2<k){
        /* assuming at the moment that error is in arcsec */
        output->raErr = LAL_TWOPI*atof(val[i+3])/(24.0*60.0*60.0);
        j+=2;
      }
    }
    else if(!strcmp(val[i],"dec") || !strcmp(val[i],"DEC") || !strcmp(val[i],"DECJ")) {
      output->dec = XLALdmsToRads(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* assuming at the moment that error is in arcsec */
        output->decErr = LAL_TWOPI*atof(val[i+3])/(360.0*60.0*60.0);
        j+=2;
      }
    }
    else if(!strcmp(val[i],"pmra") || !strcmp(val[i],"PMRA")) {
      /* convert pmra from mas/year to rads/sec */
      output->pmra = LAL_PI_180*atof(val[i+1])/(60.0*60.0*1000.*365.25*86400.);
      j++;
      if(atoi(val[i+2])==1 && i+2<k){
        output->pmraErr =
          LAL_PI_180*atof(val[i+3])/(60.0*60.0*1000.*365.25*86400.);
        j+=2;
      }
    }
    else if(!strcmp(val[i],"pmdec") || !strcmp(val[i],"PMDEC")) {
      /* convert pmdec from mas/year to rads/sec */
      output->pmdec = LAL_PI_180*atof(val[j+1])/(60.0*60.0*1000.*365.25*86400.);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->pmdecErr =
          LAL_PI_180*atof(val[i+3])/(60.0*60.0*1000.*365.25*86400.);
        j+=2;
      }
    }
    else if(!strcmp(val[i],"pepoch") || !strcmp(val[i],"PEPOCH")) {
      output->pepoch = XLALTTMJDtoGPS(atof(val[i+1])); /* convert all epochs to
        from MJD to GPS seconds in TDB */
      j++;

    }
    else if( !strcmp(val[i],"posepoch") || !strcmp(val[i],"POSEPOCH")){
      output->posepoch = XLALTTMJDtoGPS(atof(val[i+1]));
      j++;
      /* position epoch in GPS seconds TDB */
    }
    else if( !strcmp(val[i],"f0") || !strcmp(val[i],"F0")) {
      /* in .par files exponents sometimes shown as D/d rather than e/E
         need way to check this as atof will not convert D (but will
         work for e/E (if a d/D is present atof will convert the number
         before the d/D but not the exponent */
      CHAR *loc;

      output->f0 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f0Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f0Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f1") || !strcmp(val[i],"F1")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f1 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f1 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f1Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f1Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f2") || !strcmp(val[i],"F2")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f2 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f2 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f2Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f2Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f3") || !strcmp(val[i],"F3")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f3 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f3 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f3Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f3Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f4") || !strcmp(val[i],"F4")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f4 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f4 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f4Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f4Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f5") || !strcmp(val[i],"F5")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f5 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f5 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f5Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f5Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f6") || !strcmp(val[i],"F6")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f6 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f6 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f6Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f6Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f7") || !strcmp(val[i],"F7")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f7 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f7 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f7Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f7Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f8") || !strcmp(val[i],"F8")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f8 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f8 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f8Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f8Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"f9") || !strcmp(val[i],"F9")) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->f9 = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->f9 = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->f9Err = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->f9Err = atof(val[i+3]);
        }
        j+=2;
      }
    }
    else if( !strcmp(val[i],"WAVE_OM") || !strcmp(val[i],"wave_om") ) {
      output->wave_om = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->wave_omErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "WAVEEPOCH") || !strcmp(val[i], "waveepoch") ){
      output->waveepoch = XLALTTMJDtoGPS( atof(val[i+1]) );
      j++;
    }
    else if( strstr(val[i],"WAVE") != NULL || strstr(val[i],"wave") != NULL ) {
      INT4 wnum = 0;

      if( sscanf(val[i]+4, "%d", &wnum) != 1 ){
        fprintf(stderr, "Error reading WAVE number from par file\n");
        exit(1);
      }

      if ( wnum > output->nwaves ){
        output->nwaves = wnum;
        output->waveSin = XLALRealloc(output->waveSin, wnum*sizeof(REAL8));
        output->waveCos = XLALRealloc(output->waveCos, wnum*sizeof(REAL8));
      }

      output->waveSin[wnum-1] = atof(val[i+1]);
      output->waveCos[wnum-1] = atof(val[i+2]);

      j++;
    }
    else if( !strcmp(val[i],"binary") || !strcmp(val[i],"BINARY")) {
      output->model = XLALStringDuplicate(val[i+1]);
      j++;
    }
    else if( !strcmp(val[i],"units") || !strcmp(val[i],"UNITS")){
      output->units = XLALStringDuplicate(val[i+1]);
      j++;
    }
    else if( !strcmp(val[i],"ephem") || !strcmp(val[i],"EPHEM")){
      output->ephem = XLALStringDuplicate(val[i+1]);
      j++;
    }
    else if( !strcmp(val[i],"a1") || !strcmp(val[i],"A1")) {
      output->x = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->xErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"e") || !strcmp(val[i],"E") || !strcmp(val[i],"ECC") ||
      !strcmp(val[i],"ecc")) {
      output->e = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"pb") || !strcmp(val[i],"PB")) {
      output->Pb = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->PbErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"om") || !strcmp(val[i],"OM")) {
      output->w0 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->w0Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "T0")){
      output->T0 = XLALTTMJDtoGPS(atof(val[i+1]));
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->T0Err = atof(val[i+3])*DAYSTOSECS; /* convert to seconds */
        j+=2;
      }
    }
    else if( !strcmp(val[i], "Tasc") || !strcmp(val[i], "TASC")){
      output->Tasc = XLALTTMJDtoGPS(atof(val[i+1]));
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->TascErr = atof(val[i+3])*DAYSTOSECS; /* convert to seconds; */
        j+=2;
      }
    }
    else if( !strcmp(val[i], "eps1") || !strcmp(val[i], "EPS1")){
      output->eps1 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eps1Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "eps2") || !strcmp(val[i], "EPS2")){
      output->eps2 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eps2Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "eps1dot") || !strcmp(val[i], "EPS1DOT")){
      output->eps1dot = atof(val[i+1]);
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if( fabs( output->eps1dot ) > 1e-7  ) output->eps1dot *= 1.e-12;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eps1dotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "eps2dot") || !strcmp(val[i], "EPS2DOT")){
      output->eps2dot = atof(val[i+1]);
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if( fabs( output->eps2dot ) > 1e-7 ) output->eps2dot *= 1.e-12;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eps2dotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "xpbdot") || !strcmp(val[i], "XPBDOT")){
      output->xpbdot = atof(val[i+1]);
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if( fabs( output->xpbdot ) > 1e-7 ) output->xpbdot *= 1.e-12;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->xpbdotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "omdot") || !strcmp(val[i], "OMDOT")){
      output->wdot = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->wdotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "pbdot") || !strcmp(val[i], "PBDOT")){
      output->Pbdot = atof(val[i+1]);
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if( fabs( output->Pbdot ) > 1e-7 ) output->Pbdot *= 1.e-12;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->PbdotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "xdot") || !strcmp(val[i], "XDOT")){
      output->xdot = atof(val[i+1]);
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if( fabs( output->xdot ) > 1e-7 ) output->xdot *= 1.e-12;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->xdotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "edot") || !strcmp(val[i], "EDOT")){
      output->edot = atof(val[i+1]);
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if( fabs( output->edot ) > 1e-7 ) output->edot *= 1.e-12;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->edotErr = atof(val[i+3]);
        j+=2;
      }

      /* some of the parameter files in the ATNF catalogue have values
         of EDOT that are stupidly large e.g. O(1e33). These can cause
         the time delay routines to fail, so if values of EDOT are
         greater than 10000 ignore them and set it to zero */
      if( output->edot > 10000 ){
        output->edot = 0.;
        output->edotErr = 0.;
      }
    }
    else if( !strcmp(val[i], "gamma") || !strcmp(val[i], "GAMMA")){
      output->gamma = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->gammaErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "sini") || !strcmp(val[i], "SINI")){
      output->sstr = XLALStringDuplicate(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->sErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "mtot") || !strcmp(val[i], "MTOT")){
      output->M = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->MErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "m2") || !strcmp(val[i], "M2")){
      output->m2 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->m2Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "a0") || !strcmp(val[i], "A0")){
      output->a0 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->a0Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "b0") || !strcmp(val[i], "B0")){
      output->b0 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->b0Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "dr") || !strcmp(val[i], "DR")){
      output->dr = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->drErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "dtheta") || !strcmp(val[i], "DTHETA")){
      output->dth = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->dthErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "shapmax") || !strcmp(val[i], "SHAPMAX") ){
      output->shapmax = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->shapmaxErr = atof(val[i+3]);
        j+=2;
      }
    }

    /* parameters for Kopeikin terms */
    else if( !strcmp(val[i],"D_AOP") || !strcmp(val[i],"d_aop") ){
      /* convert into 1/rads (factor from T2model.C in TEMPO2 */
      output->daop = atof(val[i+1]) * 3600.0 / LAL_PI_180;
      output->daopset = 1;
      j++;
    }
    else if( !strcmp(val[i], "KIN") || !strcmp(val[i], "kin") ){
      output->kin = atof(val[i+1]) * LAL_PI_180; /* convert degs to rads */
      output->kinset = 1;
      j++;
    }
    else if( !strcmp(val[i], "KOM") || !strcmp(val[i], "kom") ){
      output->kom = atof(val[i+1]) * LAL_PI_180; /* convert degs to rads */
      output->komset = 1;
      j++;
    }

    /* parameters for distance */
    else if( !strcmp(val[i],"px") || !strcmp(val[i],"PX") ) { /* parallax */
      /* convert from mas to rads (factor from T2model.C in TEMPO2) */
      output->px = atof(val[i+1]) * LAL_PI_180 / 3600.0e3;
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->pxErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"dist") || !strcmp(val[i],"DIST") ) { /* distance */
      output->dist = atof(val[i+1]); /* in kpc */
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->distErr = atof(val[i+3]);
        j+=2;
      }
    }

    /* dispersion measure parameters */
    else if( !strcmp(val[i],"dm") || !strcmp(val[i],"DM") ) {
      output->DM = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->DMErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"dm1") || !strcmp(val[i],"DM1") ) {
      output->DM1 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->DM1Err = atof(val[i+3]);
        j+=2;
      }
    }

    /* add parameters extra orbital parameters for the BT1P and BT2P models */
    else if( !strcmp(val[i],"a1_2") || !strcmp(val[i],"A1_2")) {
      output->x2 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->x2Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"e_2") || !strcmp(val[i],"E_2") ||
      !strcmp(val[i],"ECC_2")      || !strcmp(val[i],"ecc_2")) {
      output->e2 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->e2Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"pb_2") || !strcmp(val[i],"PB_2")) {
      output->Pb2 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->Pb2Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"om_2") || !strcmp(val[i],"OM_2")) {
      output->w02 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->w02Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "T0_2")){
      output->T02 = XLALTTMJDtoGPS(atof(val[i+1]));
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->T02Err = atof(val[i+3])*DAYSTOSECS; /* convert to seconds */
        j+=2;
      }
    }
    else if( !strcmp(val[i],"a1_3") || !strcmp(val[i],"A1_3")) {
      output->x3 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->x3Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"e_3") || !strcmp(val[i],"E_3") ||
      !strcmp(val[i],"ECC_3")      || !strcmp(val[i],"ecc_3")) {
      output->e3 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->e3Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"pb_3") || !strcmp(val[i],"PB_3")) {
      output->Pb3 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->Pb3Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"om_3") || !strcmp(val[i],"OM_3")) {
      output->w03 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->w03Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "T0_3")){
      output->T03 = XLALTTMJDtoGPS(atof(val[i+1]));
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->T03Err = atof(val[i+3])*DAYSTOSECS; /* convert to seconds */
        j+=2;
      }
    }
    /* orbital frequency coefficients for BTX model (up to 12 FB coefficients), but
       only one orbit at the moment i.e. only a two body system */
    else if( val[i][0] == 'F' && val[i][1] == 'B' ){
      INT4 fbnum = 0;
      CHAR *loc;

      if (strlen(val[i])==2) fbnum = 0; /* only one coefficient */
      else{
        if( sscanf(val[i]+2,"%d",&fbnum) != 1 ){
          fprintf(stderr, "Error reading FB value from par file\n");
          exit(1);
        }
      }

      /* add to number of coefficients */
      if ( output->nfb < fbnum+1 ){
        output->fb = XLALRealloc(output->fb, (fbnum+1)*sizeof(REAL8));
        output->fbErr = XLALRealloc(output->fbErr, (fbnum+1)*sizeof(REAL8));
        output->nfb = fbnum+1;
      }

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[fbnum] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[fbnum] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[fbnum] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[fbnum] = atof(val[i+3]);
        }
        j+=2;
      }
    }
    /* read in pulsar gravitational wave parameters */
    else if( !strcmp(val[i],"h0") || !strcmp(val[i],"H0") ) {
      output->h0 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->h0Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"cosiota") || !strcmp(val[i],"COSIOTA") ||
      !strcmp(val[i],"ciota") || !strcmp(val[i],"CIOTA") ) {
      output->cosiota = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->cosiotaErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"psi") || !strcmp(val[i],"PSI") ) {
      output->psi = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->psiErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"phi0") || !strcmp(val[i],"PHI0") ) {
      output->phi0 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->phi0Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"aplus") || !strcmp(val[i],"APLUS") ) {
      output->Aplus = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->AplusErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"across") || !strcmp(val[i],"ACROSS") ) {
      output->Across = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->AcrossErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"i21") || !strcmp(val[i],"I21") ) {
      output->I21 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->I21Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"i31") || !strcmp(val[i],"I31") ) {
      output->I31 = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->I31Err = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"r") || !strcmp(val[i],"R") ) {
      output->r = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->rErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"lambda") || !strcmp(val[i],"LAMBDA") ) {
      output->lambda = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->lambdaErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i],"costheta") || !strcmp(val[i],"COSTHETA") ) {
      output->costheta = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->costhetaErr = atof(val[i+3]);
        j+=2;
      }
    }

    if(j==i){
      i++;
    }
    else{
      i+=(j-i);
    }

    if(i>=k)
      break;
  }

  /*fprintf(stderr, "Have I got to the end of LALReadPARFile.\n");*/
  fclose(fp);

  /* check linked parameters */
  if( output->sstr != NULL ){
    if( !strcmp(output->sstr, "KIN") || !strcmp(output->sstr, "kin") ){
      if ( output->kinset ) output->s = sin(output->kin);
      else{
        XLALPrintError("Error... KIN not set in .par file %s\n", pulsarAndPath);
        XLAL_ERROR_VOID( XLAL_EIO );
      }
    }
    else output->s = atof(output->sstr);
  }
}

/* function to print out to screen all the pulsar parameters and there associated errors */
void PrintPulsarParameters( BinaryPulsarParams params ){
  fprintf(stderr, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  fprintf(stderr, "PULSAR %s :\n", params.name);
  fprintf(stderr, "sky position:\tra %.7lf +/- %.3le rads, dec %.7lf +/- %.3le rads\n", params.ra,
params.raErr, params.dec, params.decErr);
  if(params.pmra != 0. || params.pmdec != 0.)
    fprintf(stderr, "proper motion:\tra %.4le +/- %.1le rads/s, dec %.4le +/- %.1le rads/s\n",
params.pmra, params.pmraErr, params.pmdec, params.pmdecErr);
  if(params.pepoch != 0. || params.posepoch != 0.)
    fprintf(stderr, "epochs:\tperiod %lf (GPS), position %lf (GPS)\n", params.pepoch,
params.posepoch);
  fprintf(stderr, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n");

  fprintf(stderr, "Frequency parameters\n");
  if(params.f0 != 0.)
    fprintf(stderr, "\tf0 = %.10lf +/- %.3le (Hz)\n", params.f0, params.f0Err);
  if(params.f1 != 0.)
    fprintf(stderr, "\tf1 = %.5le +/- %.3le (Hz/s)\n", params.f1, params.f1Err);
  if(params.f2 != 0.)
    fprintf(stderr, "\tf1 = %.5le +/- %.3le (Hz/s^2)\n", params.f2, params.f2Err);
  /* print binary parameters */
  if(params.model != NULL){
    fprintf(stderr, "\nBinary parameters:\tmodel %s\n", params.model);

    fprintf(stderr, "Keplarian parameters:-\n");
    if(params.Pb != 0.)
      fprintf(stderr, "\tperiod = %lf +/- %.3le (days)\n", params.Pb, params.PbErr);
    if(params.x != 0.)
      fprintf(stderr, "\tprojected semi-major axis = %lf +/- %.3le (light sec)\n", params.x,
params.xErr);
    if(params.e != 0.)
      fprintf(stderr, "\teccentricity = %lf +/- %.3le\n", params.e, params.eErr);
    if(params.w0 != 0.)
      fprintf(stderr, "\tlongitude of periastron = %lf +/- %.3lf (degs)\n", params.w0,
params.w0Err);
    if(params.T0 != 0.)
      fprintf(stderr, "\ttime of periastron = %lf +/- %.3lf (GPS)\n", params.T0, params.T0Err);
    if(params.Tasc != 0.)
      fprintf(stderr, "\ttime of ascending node = %lf +/- %.3lf (GPS)\n", params.Tasc,
params.TascErr);
    if(params.eps1 != 0.)
      fprintf(stderr, "\tfirst Laplace-Lagrange parameter (eps1) = %le +/- %.3le\n", params.eps1,
params.eps1Err);
    if(params.eps2 != 0.)
      fprintf(stderr, "\tsecond Laplace-Lagrange parameter (eps1) = %le +/- %.3le\n", params.eps2,
params.eps2Err);
    if(params.eps2 != 0.)
      fprintf(stderr, "\tsecond Laplace-Lagrange parameter (eps1) = %le +/- %.3le\n", params.eps2,
params.eps2Err);

    /*fprintf(stderr, "Post-Newtonian parameters:-\n");
    if(params.gamma != 0.)
      fprintf(stderr, "\tGravitational redshift parameter = %le +/- %.3le\n", params.gamma,
params.gammaErr);*/

  }
}

LALStringVector *XLALReadTEMPOCorFile( REAL8Array *cormat, CHAR *corfile )
/*void XLALReadTEMPOCorFile( REAL8Array *cormat, LALStringVector *params,
                           CHAR *corfile )*/{
  FILE *fp = NULL;
  CHAR *firstline = XLALStringDuplicate( "" );
  CHAR onechar[2];
  INT4 i = 0, numPars = 0, c = 1, sl = 0;
  LALStringVector *tmpparams = NULL; /* temporary parameter names */
  LALStringVector *params = NULL;
  UINT4Vector *dims = NULL;

  /* check the file exists */
  if( access(corfile, F_OK) != 0 ){
    XLALPrintError("Error... correlation matrix file does not exist!\n");
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* open file */
  if( (fp = fopen(corfile, "r")) == NULL ){
    XLALPrintError("Error... cannot open correlation matrix file!\n");
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* read in first line of the file */
  while( !strchr( fgets(onechar, 2, fp), '\n' ) )
    firstline = XLALStringAppend( firstline, onechar );

  sl = strlen(firstline);

  /* count the number of parameters */
  for ( i = 0; i < sl; i++ ){
    /* use isspace as delimiters could be unknown generic whitespace */
    if ( !isspace(firstline[i]) ){
      if ( c ){
        numPars++;
        c = 0;
      }
    }else
      c = 1;
  }

  /* parse the line and put into the params vector */
  rewind(fp); /* rewind to start of the file */
  for ( i = 0; i < numPars; i++ ){
    CHAR tmpStr[128];

    if( fscanf(fp, "%s", tmpStr) == EOF ){
      XLALPrintError("Error... Problem reading first line of correlation\
 matrix!\n");
      XLAL_ERROR_NULL(XLAL_EIO);
    }

    tmpparams = XLALAppendString2Vector( tmpparams, tmpStr );

    /* convert some parameter names to a more common convention */
    if ( !strcasecmp(tmpStr, "RAJ") ) /* convert RAJ to ra */
      params = XLALAppendString2Vector( params, "ra" );
    else if ( !strcasecmp(tmpStr, "DECJ") ) /* convert DECJ to dec */
      params = XLALAppendString2Vector( params, "dec" );
    else
      params = XLALAppendString2Vector( params, tmpStr );
  }

  dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = numPars;
  dims->data[1] = numPars;

  /* set the correlation matrix to the correct size */
  cormat = XLALResizeREAL8Array( cormat, dims );

  /* read through covariance values */
  for ( i = 0; i < numPars; i++ ){
    CHAR tmpStr[128];
    INT4 j = 0;

    if( fscanf(fp, "%s", tmpStr) == EOF ){
      XLALPrintError("Error... problem reading in correlation matrix!\n");
      XLAL_ERROR_NULL(XLAL_EIO);
    }

    if ( strcmp(tmpStr, tmpparams->data[i]) ){
      XLALPrintError("Error... problem reading in correlation matrix. \
Parameters not in consistent order!\n");
      XLAL_ERROR_NULL(XLAL_EIO);
    }

    for( j = 0; j < i+1; j++ ){
      REAL8 tmpval = 0.;

      if( fscanf(fp, "%lf", &tmpval) == EOF ){
        XLALPrintError("Error... problem reading in correlation matrix!\n");
        XLAL_ERROR_NULL(XLAL_EIO);
      }

      /* if off diagonal values are +/-1 set to +/- 0.99999 */
      if ( j != i && abs(tmpval) == 1. )
        tmpval *= 0.99999;

      cormat->data[i*numPars + j] = tmpval;

      /* set opposite elements */
      if( j != i )
        cormat->data[j*numPars + i] = tmpval;
    }
  }

  return params;
}


/* function to convert a string containing an angular coordinate in the format
 * degrees:minutues:seconds into radians */
REAL8 XLALdmsToRads( CHAR *dms ){
  REAL8 radians = 0., s = 0.;
  INT4 d = 0, m = 0, numitems = 0, negbutzero = 0;

  XLAL_CHECK_REAL8( dms != NULL, XLAL_EIO, "Angle string is NULL" );

  numitems = sscanf(dms, "%d:%d:%lf", &d, &m, &s);

  XLAL_CHECK_REAL8( numitems == 3, XLAL_EFUNC, "Angle string not in format 'degs:mins:secs'" );
  XLAL_CHECK_REAL8( m >= 0 && m < 60, XLAL_EFUNC, "Minutes is out of the 0 to 59 mins range" );
  XLAL_CHECK_REAL8( s >= 0. && s < 60., XLAL_EFUNC, "Seconds is out of the 0 to 60 secs range" );

  /* check if the string is negative in the case when the degrees value is zero */
  if( dms[0] == '-' && d == 0 ) negbutzero = 1;

  /* convert from dd:mm:ss to radians */
  radians = LAL_PI_180 * (REAL8)d;

  /* if dec is negative convert mins and secs to -ve numbers */
  if( d < 0 || negbutzero == 1 ){
    m = -m;
    s = -s;
  }

  radians += LAL_PI_180 * (REAL8)m / 60.;
  radians += LAL_PI_180 * s / 3600.;

  return radians;
}


/* function to convert a string containing an angular coordinate in the format
 * hours:minutues:seconds into radians */
REAL8 XLALhmsToRads( CHAR *hms ){
  REAL8 radians = 0., s = 0.;
  INT4 h = 0, m = 0, numitems = 0;

  REAL8 degsInHour = 360./24.;

  XLAL_CHECK_REAL8( hms != NULL, XLAL_EIO, "Angle string is NULL" );

  numitems = sscanf(hms, "%d:%d:%lf", &h, &m, &s);

  XLAL_CHECK_REAL8( numitems == 3, XLAL_EFUNC, "Angle string not in format 'hours:mins:secs'" );
  XLAL_CHECK_REAL8( h >= 0, XLAL_EFUNC, "Hours value must be positive" );
  XLAL_CHECK_REAL8( m >= 0 && m < 60, XLAL_EFUNC, "Minutes is out of the 0 to 59 mins range" );
  XLAL_CHECK_REAL8( s >= 0. && s < 60., XLAL_EFUNC, "Seconds is out of the 0 to 60 secs range" );

  radians = LAL_PI_180 * (REAL8)h * degsInHour;
  radians += LAL_PI_180 * ( (REAL8)m / 60.0 ) * degsInHour;
  radians += LAL_PI_180 * ( s / 3600. ) * degsInHour;

  return radians;
}


/* function converts dec or ra from format dd/hh:mm:ss.sss or format
   dd/hhmmss.ss to radians */
REAL8 LALDegsToRads(CHAR *degs, const CHAR *coord){
  REAL8 radians=0.;
  INT4 d, m;
  REAL8 s;
  CHAR dc[4]="", mc[3]="", *sc=NULL;
  CHAR *loc;
  INT4 n, negbutzero=0;

  /* if in format dd/hh:mm:ss.s do this*/
  /* locate first : */
  if((loc = strchr(degs, ':'))!=NULL){
    n = loc-degs;

    /* copy degrees part to dc */
    strncpy(dc, degs, n);
    d = atoi(dc);

    /* check if dec is negative but the degree part is zero */
    if((strchr(degs, '-') != NULL) && d == 0){
      negbutzero = 1;
    }

    /* copy minutes part to mc */
    strncpy(mc, loc+1, 2);
    m = atoi(mc);

    /* copy seconds part to sc */
    sc = XLALStringDuplicate(loc+4);
    s = atof(sc);
  }
  /* if in format hh/ddmmss.ss */
  else{
    /* find pos of decimal point . (ascii character 46) */
    loc = strchr(degs, '.');

    /* get seconds part */
    sc = XLALStringDuplicate(loc-2);
    s = atof(sc);

    /* get minutes part */
    strncpy(mc, loc-4, 2);
    m = atoi(mc);

    /* get hours or degs part part */
    /* check if first char is - (ascii character 45) */
    if(strchr(degs, '-') != NULL){
      /* first char is negative */
      strncpy(dc, loc-7, 3);
      d = atoi(dc);

      /* if dec is negative but the degrees part is zero set flag */
      negbutzero = 1;
    }
    else{
      strncpy(dc, loc-6, 2);
      d = atoi(dc);
    }
  }

  if(strstr(coord, "ra") || strstr(coord, "RA") || strstr(coord, "alpha")){
    /* convert from hh:mm:ss to radians */
    radians = LAL_PI_180*(REAL8)d*(360.0/24.0);
    radians += LAL_PI_180*((REAL8)m/60.0)*(360.0/24.0);
    radians += LAL_PI_180*(s/(60.0*60.0))*(360.0/24.0);
  }
  else if(strstr(coord, "dec") || strstr(coord, "DEC") || strstr(coord, "delta")){
    /* convert from dd:mm:ss to radians */
    radians = LAL_PI_180*(REAL8)d;

    /* if dec is negative convert mins and secs to -ve numbers */
    if(d<0 || negbutzero==1){
      m = -m;
      s = -s;
    }

    radians += LAL_PI_180*(REAL8)m/60.0;
    radians += LAL_PI_180*s/(60.0*60.0);
  }

  /* free mem */
  LALFree(sc);

  return radians;
}

/* functions for converting times given in Terrestrial time TT or TDB in MJD to
times in GPS - this is important for epochs given in .par files which are in
TDB. TT and GPS are different by a factor of 51.184 secs, this is just the
historical factor of 32.184 secs between TT and TAI (International Atomic Time)
and the other 19 seconds come from the leap seonds added between the TAI and
UTC up to the point of definition of GPS time at UTC 01/01/1980 (see
http://www.stjarnhimlen.se/comp/time.html for details) */

/* a very good paper describing the tranforms between different time systems
and why they are necessary can be found in Seidelmann and Fukushima, A&A 265
(1992) http://ukads.nottingham.ac.uk/abs/1992A%26A...265..833S */

/* This function converts a MJD format time corrected to Terrestrial Time (TT)
 * into an equivalent GPS time */
REAL8 XLALTTMJDtoGPS(REAL8 MJD){
  REAL8 GPS;

  /* Check not before the start of GPS time (MJD 44244) */
  XLAL_CHECK_REAL8 ( MJD >= GPS0MJD, XLAL_EDOM, "Input MJD time %.1f is not in\
 range, must be > %.1f.\n", MJD, GPS0MJD);

  /* there is the magical number factor of 32.184 + 19 leap seconds to the
   * start of GPS time */
  GPS = (MJD - GPS0MJD)*86400. - GPS_TDT;

  return GPS;
}


/* If you have an MJD arrival time on the Earth then this will convert it to
 * the equivalent GPS time in TDB (see Table 1 of Seidelmann and Fukushima,
 * Astronomy & Astrophysics, 265, 833-838 (1992).
 *
 * Note that LALBarycenter performs these TDBtoTT corrections (i.e. the
 * Einstein delay) when correcting a GPS time on the Earth to TDB. Also, for
 * TEMPO produced pulsar epochs given in MJD these are already in the TDB
 * system and an equivalent GPS time in the TDB can be calculated just using
 * XLALTTMJDtoGPS.
 */
REAL8 XLALTDBMJDtoGPS(REAL8 MJD){
  REAL8 GPS;
  REAL8 T, TDBtoTT;

  /* Check not before the start of GPS time */
  XLAL_CHECK_REAL8 ( MJD >= GPS0MJD, XLAL_EDOM, "Input MJD time %.1f is not in range, must be > %.1f.\n", MJD, GPS0MJD);

  /* use factors from Table 1 of Seidelmann and Fukushima, Astronomy &
   * Astrophysics, 265, 833-838 (1992) where TDB = TDT + P
   * and:
   * P = 0.0016568 sin(35999.37 degs x T + 357.5 degs) +
         0.0000224 sin(32964.5 degs x T + 246.0 degs) +
         0.0000138 sin(71998.7 degs x T + 355.0 degs) +
         0.0000048 sin(3034.9 degs x T + 25.0 degs) +
         0.0000047 sin(34777.3 degs x T + 230.0 degs)
   * and T is the elapsed time from J2000 (which has a Julian day date of
   * JD 2451545.0) in Julian centuries.*/
  T = MJD + (XLAL_MJD_REF - XLAL_EPOCH_J2000_0_JD);
  T /= 36525.; /* covert days to Julian centuries */

  /* time diff in seconds (the Einstein delay) */
  TDBtoTT = 0.0016568*sin((35999.37*T + 357.5) * LAL_PI_180) +
            0.0000224*sin((32964.5*T +  246.0) * LAL_PI_180) +
            0.0000138*sin((71998.7*T +  355.0) * LAL_PI_180) +
            0.0000048*sin((3034.9*T + 25.0) * LAL_PI_180) +
            0.0000047*sin((34777.3*T + 230.0) *LAL_PI_180);

  /* convert TDB to TT (TDB-TDBtoTT) and then convert TT to GPS */
  /* there is the magical number factor of 32.184 + 19 leap seconds to the
   * start of GPS time */
  GPS = (MJD - GPS0MJD)*86400. - GPS_TDT - TDBtoTT;

  return GPS;
}

/* If you have an MJD arrival time on the Earth then this will convert it to
 * the equivalent GPS time in TCB (see Table 1 of Seidelmann and Fukushima,
 * Astronomy & Astrophysics, 265, 833-838, 1992).
 *
 * Note that for default TEMPO2 produced pulsar epochs given in MJD these are
 * already in the TCB system and an equivalent GPS time in the TCB can be
 * calculated just using XLALTTMJDtoGPS. */
REAL8 XLALTCBMJDtoGPS(REAL8 MJD){
  REAL8 GPS;
  REAL8 Tdiff;
  REAL8 TCBtoTDB;

  /* Check not before the start of GPS time (MJD 44244) */
  XLAL_CHECK_REAL8 ( MJD >= GPS0MJD, XLAL_EDOM, "Input MJD time %.1f is not in\
 range, must be > %.1f.\n", MJD, GPS0MJD);

  /* from Seidelmann and Fukushima we have a linear drift term:
   * TCB - TDB = 1.550506e-8 x (JD - 2443144.5) x 86400
   */
  Tdiff = (MJD + XLAL_MJD_REF - 2443144.5)*86400.;
  TCBtoTDB = 1.550506e-8 * Tdiff;

  /* convert from TDB to GPS */
  GPS = XLALTDBMJDtoGPS(MJD);

  /* add extra factor as the MJD was really in TCB not TDB) */
  GPS -= TCBtoTDB;

  return GPS;
}
