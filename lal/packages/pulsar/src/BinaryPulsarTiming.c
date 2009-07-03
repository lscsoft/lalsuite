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
 * \ingroup pulsar
 * \brief Functions to calculate binary system time delays and read TEMPO pulsar parameter files
 *
 * $Id$
 *
 */

/* LAL functions to calculate the timing differences needed to
   take into account binary pulsar orbits
   Models are taken from Taylor and Weisberg (1989) and use the
   naming conventions therein and used by TEMPO */
   
/*   Also contains function to read TEMPO .par files to obtain parameters
  and errors on parameters (if available) */

/* <lalVerbatim file="BinaryPulsarTimingCV">
   Author: Pitkin, M. D.
   $Id$
   </lalVerbatim>
   
   <lalLaTeX>
   \subsection{Module \texttt{BinaryPulsarTiming.c}}
   \label{ss:BinaryPulsarTiming.c}
   
   Functions for calculating the timing delay to a signal from a pulsar in a
   binary system and reading pulsar parameters from TEMPO \cite{TEMPO} .par
   files.
   
   \subsubsection*{Prototypes}
   \vspace{0.1in}
   \input{BinaryPulsarTimingCP}
   \idx{LALBinaryPulsarDeltaT()}
   \idx{LALReadTEMPOParFile()}
   \idx{LALdegsToRads()}
   
   \subsubsection*{Description}
   
   The main function computes the time delay of a signal from a pulsar in a
   binary system due to doppler shifts and relativistic delays, 
   \begin{equation}
   \Delta{}t = t_{\rm Roemer} + t_{\rm Shapiro} + t_{\rm Einstein} + t_{\rm
   Abberation},
   \end{equation}
   where $t_{\rm Roemer}$ is the light travel time, $t_{\rm Shapiro}$ is the
   General relativistic time delay, $t_{\rm Einstein}$ is the special
   relativistic time delay, and $t_{\rm Abberation}$ is the delay caused by the
   pulsars' rotation. There are several models of the binary systems, described 
   in \cite{TaylorWeisberg:1989}, of which the four most common are so far
   implemented. The four models are the Blandford-Teukolsky model (BT)
   \cite{BlandfordTeukolsky:1976}, the low ellipticity model (ELL1)
   \cite{ChLangeetal:2001}, Damour-Deruelle model (DD) \cite{DamourDeruelle:1985}, 
   and the main sequence system model (MSS) \cite{Wex:1998}.
   These four models all use the five main binary parameters: the longitude of
   periastron $\omega_0$, the eccentricity of the orbit $e$, the orbital period
   $P$, the time of periastron/or the time of ascension of the first node 
   $T_0$/$T_{{\rm asc}}$, and the projected semi-major axis $a\sin{}i$. The are
   also many other model dependent parameters. These routines closely follow
   those used in the radio astronomy package TEMPO \cite{TEMPO}.
 
   Radio astronomers fit pulsar parameters using TEMPO which will output
   the parameters in a \verb+.par+ file. The values allowed in this file can be
   found in the TEMPO documentation. A function is included to extract these 
   parameters from the \verb+.par+ files and put them into a
   \verb+BinaryPulsarParams+ structure, it will set any unused parameters to
   zero or \texttt{NULL}. All parameters are in the units used by TEMPO with any
   conversion to SI units occuring within the binary timing routines. A function
   is also included which converts a string containing the right ascension or 
   declination in the format \texttt{ddd/hh:mm:ss.s} or \texttt{ddd/hhmmss.s} 
   (as is given in the \texttt{.par} file) into a \texttt{REAL8} value in 
   radians.
 
   \subsubsection*{Notes}
   
   \vfill{\footnotesize\input{BinaryPulsarTimingCV}}
   
   </lalLaTeX>
*/

/* Matt Pitkin 29/04/04 */

#include <lal/BinaryPulsarTiming.h>

#include <string.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/ComputeFstat.h>

#define DAYSTOSECS 86400.0 /* number of seconds in a day */

/******* DEFINE RCS ID STRING ************/
NRCSID( BINARYPULSARTIMINGC, "$Id$" );

/** Calculate the binary system time delay using the pulsar parameters in 
 *  \c params
 */
void
LALBinaryPulsarDeltaT( LALStatus            *status,
                       BinaryPulsarOutput   *output,
                       BinaryPulsarInput    *input,
                       BinaryPulsarParams   *params ){
  INITSTATUS(status, "LALBinaryPulsarDeltaT", BINARYPULSARTIMINGC);
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
         (!strcmp(params->model, "MSS")), status,
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
  const CHAR *fn = "XLALBinaryPulsarDeltaT()";

  REAL8 dt=0.; /* binary pulsar deltaT */
  REAL8 x, xdot;	/* x = asini/c */
  REAL8 w;  /* longitude of periastron */
  REAL8 e, edot;  /* eccentricity */
  REAL8 eps1, eps2;
  REAL8 eps1dot, eps2dot;
  REAL8 w0, wdot;
  REAL8 Pb, pbdot;
  REAL8 xpbdot;
  REAL8 T0, Tasc, tb=0.; /* time parameters */

  REAL8 s, r; /* Shapiro shape and range params */
  REAL8 gamma; /* time dilation and grav redshift */
  REAL8 dr, dth;

  REAL8 a0, b0;	/* abberation parameters */
  
  REAL8 M, m2;
  REAL8 c3 = (REAL8)LAL_C_SI*(REAL8)LAL_C_SI*(REAL8)LAL_C_SI;
  
  CHAR *model = params->model;

  /* Check input arguments */
  if( input == (BinaryPulsarInput *)NULL ){
    XLAL_ERROR_VOID( fn, BINARYPULSARTIMINGH_ENULLINPUT );
  }

  if( output == (BinaryPulsarOutput *)NULL ){
    XLAL_ERROR_VOID( fn, BINARYPULSARTIMINGH_ENULLOUTPUT );
  }

  if( params == (BinaryPulsarParams *)NULL ){
    XLAL_ERROR_VOID( fn, BINARYPULSARTIMINGH_ENULLPARAMS );
  }

  if((!strcmp(params->model, "BT")) &&
     (!strcmp(params->model, "BT1P")) &&
     (!strcmp(params->model, "BT2P")) &&
     (!strcmp(params->model, "BTX")) &&
     (!strcmp(params->model, "ELL1")) &&
     (!strcmp(params->model, "DD")) &&
     (!strcmp(params->model, "MSS"))){
    XLAL_ERROR_VOID( fn, BINARYPULSARTIMINGH_ENULLBINARYMODEL );
  }

  /* convert certain params to SI units */
  w0 = params->w0*LAL_PI_180; /* convert w to rads from degs */
  wdot = params->wdot*LAL_PI_180/(365.25*DAYSTOSECS); /* convert wdot to rads/s from degs/yr */

  Pb = params->Pb*DAYSTOSECS; /* covert period from days to secs */
  pbdot = params->Pbdot*1.0e-12;
  
  T0 = params->T0; /* these should be in TDB in seconds */
  Tasc = params->Tasc;	

  e = params->e;
  edot = params->edot*1.0e-12;
  eps1 = params->eps1;
  eps2 = params->eps2;
  eps1dot = params->eps1dot;
  eps2dot = params->eps2dot;

  x = params->x;
  xdot = params->xdot*1.0e-12;
  xpbdot = params->xpbdot*1.0e-12;

  gamma = params->gamma;
  s = params->s;
  dr = params->dr;
  dth = params->dth*1.0e-6;

  a0 = params->a0*1.0e-6; /* from microsecs to secs */
  b0 = params->b0*1.0e-6;

  M = params->M*LAL_MSUN_SI; /* from solar masses to kg */
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
    REAL8 du = 1.0; 
  
    INT4 nplanets=1; /* number of orbitting bodies in system */
    INT4 i=1, j=1;
    REAL8 fac=1.; /* factor in front of fb coefficients */  

    REAL8 su = 0., cu = 0.;
    REAL4 sw = 0., cw = 0.; /* phases from LUT */

    /* work out number of orbits i.e. have we got a BT1P or BT2P model */
    if(strstr(model, "BT1P") != NULL)
      nplanets = 2;
    if(strstr(model, "BT2P") != NULL)
      nplanets = 3;

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

        if( strstr(model, "BTX") != NULL ){
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
      du = 1.0;

      /* use numerical iteration to solve Kepler's eq for eccentric anomaly u */
      u = phase + e*sin(phase)*(1.0 + e*cos(phase));

      while(fabs(du) > 1.0e-12){
        du = (phase-(u-e*sin(u)))/(1.0-e*cos(u));
        u += du;
      }
      su = sin(u);
      cu = cos(u);

      /*fprintf(stderr, "Eccentric anomaly = %f, phase = %f.\n", u, phase);*/
      sin_cos_LUT(&sw, &cw, w);

      /* see eq 5 of Taylor and Weisberg (1989) */
      /**********************************************************/
      if( strstr(model, "BTX") != NULL ){
        /* dt += (x*sin(w)*(cos(u)-e) + (x*cos(w)*sqrt(1.0-e*e) +
          gamma)*sin(u))*(1.0 - params->fb[0]*(x*cos(w)*sqrt(1.0 -
          e*e)*cos(u) - x*sin(w)*sin(u))/(1.0 - e*cos(u))); */
        dt += (x*sw*(cu-e) + (x*cw*sqrt(1.0-e*e) +
          gamma)*su)*(1.0 - LAL_TWOPI*params->fb[0]*(x*cw*sqrt(1.0 -
          e*e)*cu - x*sw*su)/(1.0 - e*cu));
      }
      else{
        /* dt += (x*sin(w)*(cos(u)-e) + (x*cos(w)*sqrt(1.0-e*e) +
          gamma)*sin(u))*(1.0 - (LAL_TWOPI/Pb)*(x*cos(w)*sqrt(1.0 -
          e*e)*cos(u) - x*sin(w)*sin(u))/(1.0 - e*cos(u))); */
        dt += (x*sw*(cu-e) + (x*cw*sqrt(1.0-e*e) +
          gamma)*su)*(1.0 - (LAL_TWOPI/Pb)*(x*cw*sqrt(1.0 -
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
    q = alpha*(cos(u)-e) + (beta+gamma)*sin(u);
    r = -alpha*sin(u) + beta*cos(u);
    s = 1.0/(1.0-e*cos(u));
    dt = -(-q+(LAL_TWOPI/Pb)*q*r*s);*/
    /**********************************************************/
    /* There appears to be NO difference between either method */

    output->deltaT = -dt;
  }

  /* for ELL1 model (low eccentricity orbits so use eps1 and eps2) */
  /* see Appendix A, Ch. Lange etal, MNRAS (2001)                  */
  if(strstr(model, "ELL1") != NULL){
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

    REAL4 sp = 0., cp = 0., s2p = 0., c2p = 0.;

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
      REAL4 swint = 0., cwint = 0.;     

      ecc = sqrt(eps1*eps1 + eps2*eps2);
      ecc += edot*tt0;
      w_int = atan2(eps1, eps2);
      w_int = w_int + wdot*tt0;

      sin_cos_LUT(&swint, &cwint, w_int); 
      /* e1 = ecc*sin(w_int);
      e2 = ecc*cos(w_int); */
      e1 = ecc*swint;
      e2 = ecc*cwint;
    }

    sin_cos_LUT(&sp, &cp, phase);
    sin_cos_LUT(&s2p, &c2p, 2.*phase);

    /* this timing delay (Roemer + Einstein) should be most important in most cases */ 
    /* DRE = x*(sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase)));
    DREp = x*cos(phase);
    DREpp = -x*sin(phase); */
    DRE = x*(sp-0.5*(e1*c2p-e2*s2p));
    DREp = x*cp;
    DREpp = -x*sp;

    /* these params will normally be negligable */
    /* dlogbr = log(1.0-s*sin(phase)); */
    dlogbr = log(1.0-s*sp);
    DS = -2.0*r*dlogbr;
    /* DA = a0*sin(phase) + b0*cos(phase); */
    DA = a0*sp + b0*cp;

    Dbb = DRE*(1.0-nb*DREp+(nb*DREp)*(nb*DREp) + 0.5*nb*nb*DRE*DREpp) + DS + DA;

    output->deltaT = -Dbb;
    /********************************************************/
  }

  /* for DD model - code partly adapted from TEMPO bnrydd.f */
  /* also used for MSS model (Wex 1998) - main sequence star orbit - this only has two lines
different than DD model - TEMPO bnrymss.f */
  if(strstr(params->model, "DD") != NULL || strstr(params->model, "MSS") != NULL){
    REAL8 u, du=1.0;/* new eccentric anomaly */
    REAL8 Ae;       /* eccentricity parameter */
    REAL8 DRE;      /* Roemer delay + Einstein delay */
    REAL8	DREp, DREpp; /* see DD eqs 48 - 50 */
    REAL8 DS;       /* Shapiro delay */
    REAL8 DA;       /* aberation caused by pulsar rotation delay */
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
    
    REAL8 su = 0., cu = 0.;
    REAL4 sw = 0., cw = 0., swAe = 0., cwAe = 0.;

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

    /* use numerical iteration to solve Kepler's eq for eccentric anomaly u */
    u = phase + e*sin(phase)*(1.0 + e*cos(phase));

    while(fabs(du) > 1.0e-12){
      du = (phase-(u-e*sin(u)))/(1.0-e*cos(u));
      u += du;
    }
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
    if(strstr(params->model, "MSS") != NULL){
      x = x + xi*Ae; /* in bnrymss.f they also include a second time derivative of x (x2dot), but
this isn't defined for either of the two pulsars currently using this model */
    }
    else
      x = x + xdot*tt0;
    
    /* now compute time delays as in DD eqs 46 - 52 */

    /* calculate Einstein and Roemer delay */
    sin_cos_LUT(&sw, &cw, w);    
    /* sw = sin(w);
    cw = cos(w); */
    alpha = x*sw;
    beta = x*sqrt(1.0-eth*eth)*cw;
    bg = beta + gamma;
    DRE = alpha*(cu-er)+bg*su;
    DREp = -alpha*su + bg*cu;
    DREpp = -alpha*cu - bg*su;
    anhat = an/onemecu;

    /* calculate Shapiro and abberation delays DD eqs 26, 27 */
    sqr1me2 = sqrt(1.0-e*e);
    cume = cu-e;
    brace = onemecu-s*(sw*cume + sqr1me2*cw*su);
    dlogbr = log(brace);
    DS = -2.0*r*dlogbr;

    /* this abberation delay is prob fairly small */
    sin_cos_LUT(&swAe, &cwAe, (w+Ae));    
    /* DA = a0*(sin(w+Ae)+e*sw) + b0*(cos(w+Ae)+e*cw); */
    DA = a0*(swAe+e*sw) + b0*(cwAe+e*cw);

    /* timing difference */
    Dbb = DRE*(1.0 - anhat*DREp+anhat*anhat*DREp*DREp + 0.5*anhat*anhat*DRE*DREpp - 
          0.5*e*su*anhat*anhat*DRE*DREp/onemecu) + DS + DA;

    output->deltaT = -Dbb;
  }

  /* for DDGR model */  

  /* for Epstein-Haugan (EH) model - see Haugan, ApJ (1985) eqs 69 and 71 */
  
  /* check that the returned value is not a NaN */
  if( isnan(output->deltaT) ){
    XLAL_ERROR_VOID( fn, BINARYPULSARTIMINGH_ENAN );
  }
}


void
LALReadTEMPOParFile(  LALStatus *status,
                      BinaryPulsarParams *output,
                      CHAR      *pulsarAndPath )
{
  INITSTATUS(status, "LALReadTEMPOParFile", BINARYPULSARTIMINGC);
  ATTATCHSTATUSPTR(status);
  
  ASSERT(output != (BinaryPulsarParams *)NULL, status, 
  BINARYPULSARTIMINGH_ENULLOUTPUT, BINARYPULSARTIMINGH_MSGENULLOUTPUT);

  XLALReadTEMPOParFile( output, pulsarAndPath );

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
XLALReadTEMPOParFile( BinaryPulsarParams *output,
                      CHAR      *pulsarAndPath )
{
  const CHAR *fn = "XLALReadTEMPOParFile()";

  FILE *fp=NULL;
  CHAR val[500][40]; /* string array to hold all the read in values 
                        500 strings of max 40 characters is enough */
  INT4 i=0, j=1, k;

  if( output == (BinaryPulsarParams *)NULL ){
    XLAL_ERROR_VOID( fn, XLAL_EFAULT );
  }

  output->model = NULL; /* set binary model to null - incase not a binary */

  /* set all output params to zero*/
  output->e=0.0;      /* orbital eccentricity */
  output->Pb=0.0;     /* orbital period (days) */
  output->w0=0.0;     /* logitude of periastron (deg) */
  output->x=0.0;      /* projected semi-major axis/speed of light (light secs) */
  output->T0=0.0;     /* time of orbital perisastron as measured in TDB (MJD) */

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

  for(i=0;i<12;i++){
    output->fb[i] = 0.;
    output->fbErr[i] = 0.;
  }

  output->nfb=0;

  output->wdot=0.0;   /* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
  output->gamma=0.0;  /* gravitational redshift and time dilation parameter (s)*/
  output->Pbdot=0.0;  /* rate of change of Pb (dimensionless 10^-12) */
  output->xdot=0.0;   /* rate of change of x(=asini/c) - optional (10^-12)*/
  output->edot=0.0;   /* rate of change of e (10^-12)*/

  output->s=0.0;      /* Shapiro 'shape' parameter sin i */

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

  output->ra=0.0;
  output->dec=0.0;
  output->pmra=0.0;
  output->pmdec=0.0;

  output->px=0.;    /* parallax (mas) */
  output->dist=0.;  /* distance (kpc) */

  output->DM=0.;    /* dispersion measure */
  output->DM1=0.;   /* first derivative of dispersion measure */

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

  if((fp = fopen(pulsarAndPath, "r")) == NULL){
    XLAL_ERROR_VOID( fn, XLAL_EIO );
  }

  /* read all the pulsar data into the string array */
  while(!feof(fp)){
    /* make sure val[i] is clear first */
    sprintf(val[i], "%s", "");
    
    fscanf(fp, "%s", val[i]);
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
    else if(!strcmp(val[i],"ra") || !strcmp(val[i],"RA") || !strcmp(val[i],"RAJ")){
      /* this can be in form hh:mm:ss.ss or hhmmss.ss */
      output->ra = LALDegsToRads(val[i+1], "RA");
      j++;

      /* only try to get error if one exists */
      if(atoi(val[i+2])==1 && i+2<k){
        /* assuming at the moment that error is in arcsec */
        output->raErr = LAL_TWOPI*atof(val[i+3])/(24.0*60.0*60.0);
        j+=2;
      }
    }
    else if(!strcmp(val[i],"dec") || !strcmp(val[i],"DEC") || !strcmp(val[i],"DECJ")) {
      output->dec = LALDegsToRads(val[i+1], "DEC");
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
      output->pepoch = LALTTMJDtoGPS(atof(val[i+1])); /* convert all epochs to
        from MJD to GPS seconds in TDB */
      j++;

    }
    else if( !strcmp(val[i],"posepoch") || !strcmp(val[i],"POSEPOCH")){
      output->posepoch = LALTTMJDtoGPS(atof(val[i+1]));
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
    else if( !strcmp(val[i],"binary") || !strcmp(val[i],"BINARY")) {
      /*sprintf(output->model, "%s", val[j+1]);*/
      output->model = XLALStringDuplicate(val[i+1]);
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
      output->T0 = LALTTMJDtoGPS(atof(val[i+1]));
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->T0Err = atof(val[i+3])*DAYSTOSECS; /* convert to seconds */
        j+=2;
      }
    }
    else if( !strcmp(val[i], "Tasc") || !strcmp(val[i], "TASC")){
      output->Tasc = LALTTMJDtoGPS(atof(val[i+1]));
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
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eps1dotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "eps2dot") || !strcmp(val[i], "EPS2DOT")){
      output->eps2dot = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->eps2dotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "xpbdot") || !strcmp(val[i], "XPBDOT")){
      output->xpbdot = atof(val[i+1]);
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
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->PbdotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "xdot") || !strcmp(val[i], "XDOT")){
      output->xdot = atof(val[i+1]);
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->xdotErr = atof(val[i+3]);
        j+=2;
      }
    }
    else if( !strcmp(val[i], "edot") || !strcmp(val[i], "EDOT")){
      output->edot = atof(val[i+1]);
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
      output->s = atof(val[i+1]);
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

    /* parameters for distance */
    else if( !strcmp(val[i],"px") || !strcmp(val[i],"PX") ) { /* parallax */
      output->px = atof(val[i+1]); /* in milliarcsecs */
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
      output->T02 = LALTTMJDtoGPS(atof(val[i+1]));
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
      output->T03 = LALTTMJDtoGPS(atof(val[i+1]));
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        output->T03Err = atof(val[i+3])*DAYSTOSECS; /* convert to seconds */
        j+=2;
      }
    }

    /* orbital frequency coefficients for BTX model (up to 12 coefficients), but
       only one orbit at the moment i.e. only a two body system */
    else if( !strcmp(val[i], "fb0") || !strcmp(val[i], "FB0") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[0] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[0] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[0] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[0] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++; /* add to number of coefficients */
    }
    else if( !strcmp(val[i], "fb1") || !strcmp(val[i], "FB1") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[1] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[1] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[1] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[1] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb2") || !strcmp(val[i], "FB2") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[2] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[2] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[2] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[2] = atof(val[i+3]);
        }
        j+=2;
      }
    
      output->nfb++;
    }
    else if( !strcmp(val[i], "fb3") || !strcmp(val[i], "FB3") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[3] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[3] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[3] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[3] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb4") || !strcmp(val[i], "FB4") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[4] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[4] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[4] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[4] = atof(val[i+3]);
        }
        j+=2;
      }
      
      output->nfb++;
    }
    else if( !strcmp(val[i], "fb5") || !strcmp(val[i], "FB5") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[5] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[5] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[5] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[5] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb6") || !strcmp(val[i], "FB6") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[6] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[6] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[6] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[6] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb7") || !strcmp(val[i], "FB7") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[7] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[7] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[7] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[7] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb8") || !strcmp(val[i], "FB8") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[8] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[8] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[8] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[8] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb9") || !strcmp(val[i], "FB9") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[9] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[9] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[9] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[9] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb10") || !strcmp(val[i], "FB10") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[10] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[10] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[10] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[10] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
    }
    else if( !strcmp(val[i], "fb11") || !strcmp(val[i], "FB11") ){
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if((loc = strstr(val[i+1], "D"))!=NULL || (loc = strstr(val[i+1], "d"))!=NULL){
        output->fb[11] = atof(val[i+1])*pow(10, atof(loc+1));
      }
      else{
        output->fb[11] = atof(val[i+1]);
      }
      j++;

      if(atoi(val[i+2])==1 && i+2<k){
        /* check if exponent contains e/E or d/D or neither */
        if((loc = strstr(val[i+3], "D"))!=NULL || (loc = strstr(val[i+3], "d"))!=NULL){
          output->fbErr[11] = atof(val[i+3])*pow(10, atof(loc+1));
        }
        else{
          output->fbErr[11] = atof(val[i+3]);
        }
        j+=2;
      }

      output->nfb++;
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

/* Matt - you have tested these function using the radio pulsar data that Michael 
Kramer sent you and they are correct and are especially needed for the binary system 
epochs */

/* a very good paper describing the tranforms between different time systems
and why they are necessary can be found in Seidelmann and Fukushima, A&A 265
(1992) http://ukads.nottingham.ac.uk/abs/1992A%26A...265..833S */

REAL8 LALTTMJDtoGPS(REAL8 MJD){
  REAL8 GPS;

  /* Check not before the start of GPS time (MJD 44244) */
  if(MJD < 44244.){
    fprintf(stderr, "Input time is not in range.\n");
    exit(0);
  } 

  /* there is the magical number factor of 32.184 + 19 leap seconds to the
start of GPS time */
  GPS = (MJD-44244.)*86400. - 51.184;

  return GPS;
}

REAL8 LALTDBMJDtoGPS(REAL8 MJD){
  REAL8 GPS;
  REAL8 Tdiff, TDBtoTT;

  /* Check not before the start of GPS time (MJD 44244) */
  if(MJD < 44244.){
    fprintf(stderr, "Input time is not in range.\n");
    exit(0);
  } 
  
  /* use factors from Seidelmann and Fukushima (their factors in the sin terms
     are ~36525 times larger than what I use here as the time T they use is in
     Julian centuries rather than Julian days) */
  Tdiff = MJD + (2400000.5-2451545.0);

  /* time diff in seconds */
  TDBtoTT = 0.0016568*sin((357.5 + 0.98560028*Tdiff) * LAL_PI_180) +
            0.0000224*sin((246.0 + 0.90251882*Tdiff) * LAL_PI_180) +
            0.0000138*sin((355.0 + 1.97121697*Tdiff) * LAL_PI_180) +
            0.0000048*sin((25.0 + 0.08309103*Tdiff) * LAL_PI_180) + 
            0.0000047*sin((230.0 + 0.95215058*Tdiff) *LAL_PI_180);

  /* convert TDB to TT (TDB-TDBtoTT) and then convert TT to GPS */
  /* there is the magical number factor of 32.184 + 19 leap seconds to the
start of GPS time */
  GPS = (MJD-44244.)*86400. - 51.184 - TDBtoTT;

  return GPS;
}

/*TEMPO2 uses TCB rather than TDB so have another function for that conversion*/
REAL8 LALTCBMJDtoGPS(REAL8 MJD){
  REAL8 GPS;
  REAL8 Tdiff;
  REAL8 TCBtoTDB;

  /* Check not before the start of GPS time (MJD 44244) */
  if(MJD < 44244.){
    fprintf(stderr, "Input time is not in range.\n");
    exit(0);
  }

  /* from Seidelmann and Fukushima */
  Tdiff = (MJD + 2400000.5 - 2443144.5)*86400.;
  TCBtoTDB = 1.550506e-8 * Tdiff;

  /* convert from TDB to GPS */
  GPS = LALTDBMJDtoGPS(MJD);
  
  /* add extra factor as the MJD was really in TCB not TDB) */
  GPS -= TCBtoTDB;

  return GPS;
}
