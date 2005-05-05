/* LAL functions to calculate the timing differences needed to
	 take into account binary pulsar orbits
	 Models are taken from Taylor and Weisberg (1989) and use the
	 naming conventions therein and used by TEMPO									*/
	 
/*	 Also contains function to read TEMPO .par files to obtain parameters
	and errors on parameters (if available) */

/* <lalVerbatim file="BinaryPulsarTiming">
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
	 \idx{degsToRads()}
	 
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
	 in \cite{TaylorWeisberg:1989}, of which the three most common a so far
	 implemented. The three models are Blandford-Teukolsky (BT)
	 \cite{BlandfordTeukolsky:1976}, the low ellipticity model (ELL1)
	 \cite{ChLangeetal:2001} and Damour-Deruelle (DD) \cite{DamourDeruelle:1985}.
	 These three models all use the five main binary parameters: the longitude of
	 periastron $\omega_0$, the eccentricity of the orbit $e$, the orbital period
	 $P$, the time of periastron/or the time of ascension of the first node 
	 $T_0$/$T_{{\rm asc}}$, and the projected semi-major axis $a\sin{}i$. The are
	 also many other model dependent parameters. These routines closely follow
	 those used in the radio astronomy package TEMPO \cite{TEMPO}.
	 
	 Radio astronomer will fit pulsar parameters using TEMPO which will output the
	 parameters in a \verb+.par+ file. The values allowed in this file can be
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

#include "BinaryPulsarTiming.h"

#include <string.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>

#define DAYSTOSECS 86400.0 /* number of seconds in a day */

/******* DEFINE RCS ID STRING ************/
NRCSID( BINARYPULSARTIMINGC, "$Id$" );

void
LALBinaryPulsarDeltaT( LALStatus						*status,
											 BinaryPulsarOutput		*output,
											 BinaryPulsarInput		*input,
											 BinaryPulsarParams		*params )
{
	REAL8 dt;	/* binary pulsar deltaT */
	REAL8 x, xdot;	/* x = asini/c */
	REAL8 w;	/* longitude of periastron */
	REAL8 e, edot;	/* eccentricity */
	REAL8 eps1, eps2;
	REAL8 eps1dot, eps2dot;
	REAL8 w0, wdot;
	REAL8 Pb, pbdot;
	REAL8 xpbdot;
	REAL8 T0, Tasc, tb; /* time parameters */
	
	REAL8 s, r; /* Shapiro shape and range params */
	REAL8 gamma; /* time dilation and grav redshift */
	REAL8 dr, dth;
	
	REAL8 a0, b0;	/* abberation parameters */
	
	REAL8 M, m2;
	REAL8 c3 = (REAL8)LAL_C_SI*(REAL8)LAL_C_SI*(REAL8)LAL_C_SI;
	
	CHAR *model = params->model;
	
	INITSTATUS(status, "LALBinaryPulsarDeltaT", BINARYPULSARTIMINGC);
	ATTATCHSTATUSPTR(status);
	
	/* Check input arguments */
	ASSERT(input != (BinaryPulsarInput *)NULL, status, 
	BINARYPULSARTIMINGH_ENULLINPUT, BINARYPULSARTIMINGH_MSGENULLINPUT);
	
	ASSERT(output != (BinaryPulsarOutput *)NULL, status, 
	BINARYPULSARTIMINGH_ENULLOUTPUT, BINARYPULSARTIMINGH_MSGENULLOUTPUT);
	
	ASSERT(params != (BinaryPulsarParams *)NULL, status, 
	BINARYPULSARTIMINGH_ENULLPARAMS, BINARYPULSARTIMINGH_MSGENULLPARAMS);
	
	ASSERT((strstr(params->model, "BT") != NULL) || 
				 (strstr(params->model, "ELL1") != NULL) ||
				 (strstr(params->model, "DD") != NULL), status,
				 BINARYPULSARTIMINGH_ENULLBINARYMODEL, 
				 BINARYPULSARTIMINGH_MSGNULLBINARYMODEL);
	
	/* convert certain params to SI units */
	w0 = params->w0*LAL_PI_180; /* convert w to rads from degs */
	wdot = params->wdot*LAL_PI_180/(365.25*DAYSTOSECS); /* convert wdot to rads/s from degs/yr */
	
	Pb = params->Pb*DAYSTOSECS; /* covert period from days to secs */
	pbdot = params->Pbdot*1.0e-12;
	
	T0 = (params->T0 - 44244.0)*DAYSTOSECS; /* covert T0 from MJD to UTC */
	Tasc = (params->Tasc - 44244.0)*DAYSTOSECS; /* convert Tasc from MJD to UTC */
	
	if(strstr(input->tbflag, "MJD") != NULL){
		tb = (input->tb - 44244.0)*DAYSTOSECS;
	}
	else if(strstr(input->tbflag, "GPS") != NULL){
		tb = input->tb - input->leapSecs;
	}
	else{
		ASSERT((strstr(input->tbflag, "GPS") != NULL) || 
		(strstr(input->tbflag, "MJD") != NULL), status, 
		BINARYPULSARTIMINGH_ENULLTBFLAG, BINARYPULSARTIMINGH_MSGENULLTBFLAG);
	}
	 
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
	if(T0 == 0.0 && Tasc != 0.0){
		REAL8 fe, uasc, Dt; /* see TEMPO tasc2t0.f */
			
		fe = sqrt((1.0 - e)/(1.0 + e));
		uasc = 2.0*atan(fe*tan(w0/2.0));
		Dt = (Pb/LAL_TWOPI)*(uasc-e*sin(uasc));
			
		T0 = Tasc + Dt;
	}
			
	/* for BT model */
	if(strstr(model, "BT") != NULL){
		REAL8 tt0;
		REAL8 orbits; 
		INT4 norbits;
		REAL8 phase; /* same as mean anomaly */
		REAL8 u = 0.0; /* eccentric anomaly */
		REAL8 du = 1.0; 
		
		/* set some vars for bnrybt.f (TEMPO) method */
		/*REAL8 tt;
		REAL8 som;
		REAL8 com;
		REAL8 alpha, beta;*/
		/*REAL8 q, r, s;*/
		
		/*fprintf(stderr, "You are using the Blandford-Teukolsky (BT) binary model.\n");*/		
		
		tt0 = tb - T0;
				
		x = x + xdot*tt0;
		e = e + edot*tt0;
		w = w0 + wdot*tt0; /* calculate w */
		
		orbits = tt0/Pb - 0.5*(pbdot+xpbdot)*(tt0/Pb)*(tt0/Pb);
		norbits = (INT4)floor(orbits);
		
		if(orbits < 0)
			norbits = norbits - 1;
			
		phase = LAL_TWOPI*(orbits - norbits); /* called phase in TEMPO */
		/*phase = LAL_TWOPI*(orbits);*/
		du = 1.0;
		
		/* use numerical iteration to solve Kepler's eq for eccentric anomaly u */
		u = phase + e*sin(phase)*(1 + e*cos(phase));
		while(fabs(du) > 1.0e-12){
			du = (phase-(u-e*sin(u)))/(1.0-e*cos(u));
			u += du;
		}
		
		/*fprintf(stderr, "Eccentric anomaly = %f, phase = %f.\n", u, phase);*/
		
		/* see eq 5 of Taylor and Weisberg (1989) */
		/**********************************************************/
		dt = (x*sin(w)*(cos(u)-e) + (x*cos(w)*sqrt(1.0-e*e) +
		gamma)*sin(u))*(1.0 -
		(LAL_TWOPI/Pb)*(x*cos(w)*sqrt(1.0-e*e)*cos(u) -
		x*sin(w)*sin(u))/(1.0 - e*cos(u)));
		/**********************************************************/
		
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
	/* see Appendix A, Ch. Lange etal, MNRAS (2001)									 */
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
		norbits = (INT4)floor(orbits);				
		if(orbits < 0.0)
			norbits = norbits - 1.0;
			
		phase=LAL_TWOPI*(orbits - norbits);
		/*phase=LAL_TWOPI*(orbits);*/
		
		x = x + xdot*tt0;
						
		/* depending on whether we have eps derivs or w time derivs */
		/* calculate e1 and e2 accordingly													*/
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
		
		/* this timing delay (Roemer + Einstein) should be most important 
				in most cases																									*/	 
		DRE = x*(sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase)));
		DREp = x*cos(phase);
		DREpp = -x*sin(phase);
		
		/* these params will normally be negligable */
		dlogbr = log(1.0-s*sin(phase));
		DS = -2.0*r*dlogbr;
		DA = a0*sin(phase) + b0*cos(phase);
		
		Dbb = DRE*(1.0-nb*DREp+(nb*DREp)*(nb*DREp) + 0.5*nb*nb*DRE*DREpp) + 
					DS + DA;
		
		output->deltaT = -Dbb;
		/********************************************************/
	}
	
	/* for DD model - code partly adapted from TEMPO bnrydd.f */
	if(strstr(params->model, "DD") != NULL){
		REAL8 u, du=1.0;/* new eccentric anomaly */
		REAL8 Ae;				/* eccentricity parameter */
		REAL8 DRE; 			/* Roemer delay + Einstein delay */
		REAL8	DREp, DREpp; /* see DD eqs 48 - 50 */
		REAL8 DS;				/* Shapiro delay */
		REAL8 DA;				/* aberation caused by pulsar rotation delay */
		REAL8 tt0;
		/* various variable use during calculation */
		REAL8 er, eth, an, k;
		REAL8 orbits, phase;
		INT4 norbits;
		REAL8 w, su, cu, onemecu, cae, sae;
		REAL8 sw, cw, alpha, beta, bg;
		REAL8 anhat, sqr1me2, cume, brace, dlogbr;
		REAL8 Dbb;		/* Delta barbar in DD eq 52 */
		
		/* fprintf(stderr, "You are using the Damour-Deruelle (DD) binary model.\n");*/
		
		/* part of code adapted from TEMPO bnrydd.f */
		an = LAL_TWOPI/Pb;
		k = wdot/an;
		
		tt0 = tb - T0;
		x = x + xdot*tt0;
		e = e + edot*tt0;
		er = e*(1.0+dr);
		eth = e*(1.0+dth);
		
		orbits = (tt0/Pb) - 0.5*(pbdot+xpbdot)*(tt0/Pb)*(tt0/Pb);
		norbits = (INT4)floor(orbits);
		
		if(orbits < 0.0)
			norbits = norbits - 1;
			
		phase = LAL_TWOPI*(orbits - norbits);
		/*phase = LAL_TWOPI*(orbits);*/
		
		/* use numerical iteration to solve Kepler's eq for eccentric anomaly u */
		u = phase + e*sin(phase)*(1 + e*cos(phase));
		while(fabs(du) > 1.0e-12){
			du = (phase-(u-e*sin(u)))/(1.0-e*cos(u));
			u += du;
		}
		
		/* compute Ae as in TEMPO bnrydd.f */
		su = sin(u);
		cu = cos(u);
		onemecu = 1.0 - e*cu;
		cae = (cu - e)/onemecu;
		sae = sqrt(1.0 - e*e)*su/onemecu;
			
		Ae = atan2(sae,cae);
			
		if(Ae < 0.0)
			Ae = Ae + LAL_TWOPI;
				
		Ae = LAL_TWOPI*orbits + Ae - phase;
				
		w = w0 + k*Ae; /* add corrections to omega */
		
		/* now compute time delays as in DD eqs 46 - 52 */
		
		/* calculate Einstein and Roemer delay */
		sw = sin(w);
		cw = cos(w);
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
		DA = a0*(sin(w+Ae)+e*sw) + b0*(cos(w+Ae)+e*cw);
		
		/* timing difference */
		Dbb = DRE*(1.0 - anhat*DREp+anhat*anhat*DREp*DREp +
					0.5*anhat*anhat*DRE*DREpp - 
					0.5*e*su*anhat*anhat*DRE*DREp/onemecu) + DS + DA;
					
		output->deltaT = -Dbb;
	}
	
	/* for DDGR model */
	
	/* for Epstein-Haugan (EH) model - see Haugan, ApJ (1985) eqs 69 and 71 */
	
	/* other models, e.g. BT2P for 1 binary pulsar in current catalogue */
	
	DETATCHSTATUSPTR(status);
	RETURN(status);
		
}

void
LALReadTEMPOParFile(	LALStatus *status,
											BinaryPulsarParams *output,
											CHAR			*pulsarAndPath )
{
	FILE *fp;
	CHAR val[500][40]; /* string array to hold all the read in values 
												500 strings of max 40 characters is enough */	
	INT4 i=0, j=1, k;
	
	INITSTATUS(status, "LALReadTEMPOParFile", BINARYPULSARTIMINGC);
	ATTATCHSTATUSPTR(status);
	
	ASSERT(output != (BinaryPulsarParams *)NULL, status, 
	BINARYPULSARTIMINGH_ENULLOUTPUT, BINARYPULSARTIMINGH_MSGENULLOUTPUT);
	
	output->model = NULL; /* set binary model to null - incase not a binary */
	
	/* set all output params to zero*/
	output->e=0.0;			/* orbital eccentricity */
	output->Pb=0.0;			/* orbital period (days) */
	output->w0=0.0;			/* logitude of periastron (deg) */
	output->x=0.0;			/* projected semi-major axis/speed of light (light secs) */
	output->T0=0.0;			/* time of orbital perisastron as measured in TDB (MJD) */
	
	output->xpbdot=0.0;	/* (10^-12) */
	
	output->eps1=0.0;				/* e*sin(w) */
	output->eps2=0.0;				/* e*cos(w) */
	output->eps1dot=0.0;
	output->eps2dot=0.0;
	output->Tasc=0.0;		/* time of the ascending node (used rather than T0) */
	
	output->wdot=0.0;		/* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
	output->gamma=0.0;	/* gravitational redshift and time dilation parameter (s)*/
	output->Pbdot=0.0;	/* rate of change of Pb (dimensionless 10^-12) */
	output->xdot=0.0;		/* rate of change of x(=asini/c) - optional (10^-12)*/
	output->edot=0.0;		/* rate of change of e (10^-12)*/
	
	output->s=0.0;			/* Shapiro 'shape' parameter sin i */
	
	/*output.r=0.0;			/* Shapiro 'range' parameter */
	output->dr=0.0;			
	output->dth=0.0;		/* (10^-6) */
	output->a0=0.0;
	output->b0=0.0;	/* abberation delay parameters */
	
	output->M=0.0;			/* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
	output->m2=0.0;			/* companion mass */
	
	output->f0=0.0;
	output->f1=0.0;
	output->f2=0.0;
	output->f3=0.0;
	
	output->ra=0.0;
	output->dec=0.0;
	output->pmra=0.0;
	output->pmdec=0.0;
	
	output->raErr=0.0;
	output->decErr=0.0;
	output->pmraErr=0.0;
	output->pmdecErr=0.0;	
	
	output->posepoch=0.0;
	output->pepoch=0.0;
	
	output->posepochErr=0.0;
	output->pepochErr=0.0;
	
	/* set all errors on params to zero */
	output->xpbdotErr=0.0;	/* (10^-12) */
	
	output->eps1Err=0.0;				/* e*sin(w) */
	output->eps2Err=0.0;				/* e*cos(w) */
	output->eps1dotErr=0.0;
	output->eps2dotErr=0.0;
	output->TascErr=0.0;		/* time of the ascending node (used rather than T0) */
	
	output->wdotErr=0.0;		/* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
	output->gammaErr=0.0;	/* gravitational redshift and time dilation parameter (s)*/
	output->PbdotErr=0.0;	/* rate of change of Pb (dimensionless 10^-12) */
	output->xdotErr=0.0;		/* rate of change of x(=asini/c) - optional (10^-12)*/
	output->edotErr=0.0;		/* rate of change of e (10^-12)*/
	
	output->sErr=0.0;			/* Shapiro 'shape' parameter sin i */
	
	/*output->rErr=0.0;			/* Shapiro 'range' parameter */
	output->drErr=0.0;			
	output->dthErr=0.0;		/* (10^-6) */
	output->a0Err=0.0;
	output->b0Err=0.0;	/* abberation delay parameters */
	
	output->MErr=0.0;			/* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
	output->m2Err=0.0;			/* companion mass */
	
	output->f0Err=0.0;
	output->f1Err=0.0;
	output->f2Err=0.0;
	output->f3Err=0.0;
	
	output->eErr =0.0;
	output->w0Err=0.0;
	output->PbErr=0.0;
	output->xErr=0.0;
	output->T0Err=0.0;
	
	fp = fopen(pulsarAndPath, "r");
	
	ASSERT(fp!=NULL, status, BINARYPULSARTIMINGH_EPARFILEERROR,
				 BINARYPULSARTIMINGH_MSGEPARFILEERROR);
		
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
	while(1){
		j=i;
		if(!strcmp(val[i], "NAME") || !strcmp(val[i], "name")){
			output->name = val[i+1];
			
			j++;
		}
		else if(!strcmp(val[i],"ra") || !strcmp(val[i],"RA") || !strcmp(val[i],"RAJ")){
			/* this can be in form hh:mm:ss.ss or hhmmss.ss */
			output->ra = LALDegsToRads(val[i+1], "RA");
			j++;
			
			/* only try to get error if one exists */
			if(atoi(val[i+2])==1){
				/* assuming at the moment that error is in arcsec */
				output->raErr = LAL_TWOPI*atof(val[i+3])/(24.0*60.0*60.0);
				j+=2;
			}
		}
		else if(!strcmp(val[i],"dec") || !strcmp(val[i],"DEC") || !strcmp(val[i],"DECJ")) {
			output->dec = LALDegsToRads(val[i+1], "DEC");
			j++;
			
			if(atoi(val[i+2])==1){
				/* assuming at the moment that error is in arcsec */
				output->decErr = LAL_TWOPI*atof(val[i+3])/(360.0*60.0*60.0);
				j+=2;
			}
		}
		else if(!strcmp(val[i],"pmra") || !strcmp(val[i],"PMRA")) {
			/* convert pmra from mas/year to rads/sec */
			output->pmra = LAL_PI_180*atof(val[i+1])/(60.0*60.0*1000*LAL_YRSID_SI);
			j++;
			if(atoi(val[i+2])==1){
				output->pmraErr = LAL_PI_180*atof(val[i+3])/(60.0*60.0*1000*LAL_YRSID_SI);
				j+=2;
			}
		}
		else if(!strcmp(val[i],"pmdec") || !strcmp(val[i],"PMDEC")) {
			/* convert pmdec from mas/year to rads/sec */
			output->pmdec = LAL_PI_180*atof(val[j+1])/(60.0*60.0*1000*LAL_YRSID_SI);
			j++;
			
			if(atoi(val[i+2])==1){
				output->pmdecErr = LAL_PI_180*atof(val[i+3])/(60.0*60.0*1000*LAL_YRSID_SI);
				j+=2;
			}
		}
		else if(!strcmp(val[i],"pepoch") || !strcmp(val[i],"PEPOCH")) {
			output->pepoch = atof(val[i+1]);
			j++;
			
		}
		else if( !strcmp(val[i],"posepoch") || !strcmp(val[i],"POSEPOCH")){
			output->posepoch = atof(val[i+1]);
			j++;
			/* position epoch in MJD */
		}
		else if( !strcmp(val[i],"f0") || !strcmp(val[i],"F0")) {
			/* in .par files exponents sometimes shown as D/d rather than e/E
				 need way to check this as atof will not convert D (but will 
				 work for e/E (if a d/D is present atof will convert the number 
				 before the d/D but not the exponent */
			CHAR *loc;
				
			output->f0 = atof(val[i+1]);
			j++;
			
			if(atoi(val[i+2])==1){
				/* check if exponent contains e/E or d/D or neither */
				if((loc = strstr(val[i+3], "D"))!=NULL || 
						(loc = strstr(val[i+3], "d"))!=NULL){
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
			if((loc = strstr(val[i+1], "D"))!=NULL || 
					(loc = strstr(val[i+1], "d"))!=NULL){
				output->f1 = atof(val[i+1])*pow(10, atof(loc+1));
			}
			else{
				output->f1 = atof(val[i+1]);
			}
			j++;
			
			if(atoi(val[i+2])==1){
				/* check if exponent contains e/E or d/D or neither */
				if((loc = strstr(val[i+3], "D"))!=NULL || 
						(loc = strstr(val[i+3], "d"))!=NULL){
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
			if((loc = strstr(val[i+1], "D"))!=NULL || 
					(loc = strstr(val[i+1], "d"))!=NULL){
				output->f2 = atof(val[i+1])*pow(10, atof(loc+1));
			}
			else{
				output->f2 = atof(val[i+1]);
			}
			j++;
			
			if(atoi(val[i+2])==1){
				/* check if exponent contains e/E or d/D or neither */
				if((loc = strstr(val[i+3], "D"))!=NULL || 
						(loc = strstr(val[i+3], "d"))!=NULL){
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
			if((loc = strstr(val[i+1], "D"))!=NULL || 
					(loc = strstr(val[i+1], "d"))!=NULL){
				output->f3 = atof(val[i+1])*pow(10, atof(loc+1));
			}
			else{
				output->f3 = atof(val[i+1]);
			}
			j++;
				
			if(atoi(val[i+2])==1){
				/* check if exponent contains e/E or d/D or neither */
				if((loc = strstr(val[i+3], "D"))!=NULL || 
						(loc = strstr(val[i+3], "d"))!=NULL){
					output->f3Err = atof(val[i+3])*pow(10, atof(loc+1));
				}
				else{
					output->f3Err = atof(val[i+3]);
				}
				j+=2;
			}
    }
		else if( !strcmp(val[i],"binary") || !strcmp(val[i],"BINARY")) {
			/*sprintf(output->model, "%s", val[j+1]);*/
			output->model = val[i+1];
			
			j++;
    }
    else if( !strcmp(val[i],"a1") || !strcmp(val[i],"A1")) {
			while(atof(val[j+1])!=0){
				if(j==i)
					output->x = atof(val[j+1]);
				if(j==i+2)
					output->xErr = atof(val[j+1]);
				j++;
			}
		}
    else if( !strcmp(val[i],"e") || !strcmp(val[i],"E") || !strcmp(val[i],"ECC") ||
		!strcmp(val[i],"ecc")) {
			while(atof(val[j+1])!=0){
				if(j==i)
					output->e = atof(val[j+1]);
				if(j==i+2)
					output->eErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i],"pb") || !strcmp(val[i],"PB")) {
			while(atof(val[j+1])!=0){
				if(j==i)
					output->Pb = atof(val[j+1]);
				if(j==i+2)
					output->PbErr = atof(val[j+1]);
				j++;
			}
    }
    else if( !strcmp(val[i],"om") || !strcmp(val[i],"OM")) {
			while(atof(val[j+1])!=0){
				if(j==i)
					output->w0 = atof(val[j+1]);
				if(j==i+2)
					output->w0Err = atof(val[j+1]);
				j++;
			}
    }
		else if( !strcmp(val[i], "T0")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->T0 = atof(val[j+1]);
				if(j==i+2)
					output->T0Err = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "Tasc") || !strcmp(val[i], "TASC")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->Tasc = atof(val[j+1]);
				if(j==i+2)
					output->TascErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "eps1") || !strcmp(val[i], "EPS1")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->eps1 = atof(val[j+1]);
				if(j==i+2)
					output->eps1Err = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "eps2") || !strcmp(val[i], "EPS2")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->eps2 = atof(val[j+1]);
				if(j==i+2)
					output->eps2Err = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "eps1dot") || !strcmp(val[i], "EPS1DOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->eps1dot = atof(val[j+1]);
				if(j==i+2)
					output->eps1dotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "eps2dot") || !strcmp(val[i], "EPS2DOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->eps2dot = atof(val[j+1]);
				if(j==i+2)
					output->eps2dotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "xpbdot") || !strcmp(val[i], "XPBDOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->xpbdot = atof(val[j+1]);
				if(j==i+2)
					output->xpbdotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "omdot") || !strcmp(val[i], "OMDOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->wdot = atof(val[j+1]);
				if(j==i+2)
					output->wdotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "pbdot") || !strcmp(val[i], "PBDOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->Pbdot = atof(val[j+1]);
				if(j==i+2)
					output->PbdotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "xdot") || !strcmp(val[i], "XDOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->xdot = atof(val[j+1]);
				if(j==i+2)
					output->xdotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "edot") || !strcmp(val[i], "EDOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->edot = atof(val[j+1]);
				if(j==i+2)
					output->edotErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "gamma") || !strcmp(val[i], "GAMMA")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->gamma = atof(val[j+1]);
				if(j==i+2)
					output->gammaErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "sini") || !strcmp(val[i], "SINI")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->s = atof(val[j+1]);
				if(j==i+2)
					output->sErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "mtot") || !strcmp(val[i], "MTOT")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->M = atof(val[j+1]);
				if(j==i+2)
					output->MErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "m2") || !strcmp(val[i], "M2")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->m2 = atof(val[j+1]);
				if(j==i+2)
					output->m2Err = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "a0") || !strcmp(val[i], "A0")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->a0 = atof(val[j+1]);
				if(j==i+2)
					output->a0Err = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "b0") || !strcmp(val[i], "B0")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->b0 = atof(val[j+1]);
				if(j==i+2)
					output->b0Err = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "dr") || !strcmp(val[i], "DR")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->dr = atof(val[j+1]);
				if(j==i+2)
					output->drErr = atof(val[j+1]);
				j++;
			}
		}
		else if( !strcmp(val[i], "dtheta") || !strcmp(val[i], "DTHETA")){
			while(atof(val[j+1])!=0){
				if(j==i)
					output->dth = atof(val[j+1]);
				if(j==i+2)
					output->dthErr = atof(val[j+1]);
				j++;
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
		
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

/* function converts dec or ra from format dd/hh:mm:ss.sss or format 
   dd/hhmmss.ss to radians */
REAL8 LALDegsToRads(CHAR *degs, CHAR *coord){
	REAL8 radians;
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
		sc = strdup(loc+4);
		s = atof(sc);
	}
	/* if in format hh/ddmmss.ss */
	else{
		/* find pos of decimal point . (ascii character 46) */
		loc = strchr(degs, '.');
						
		/* get seconds part */
		sc = strdup(loc-2);
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
	
	return radians;
}

	
