/* LAL code to take into account binary pulsar motion */
/* 	for definitions of different models see Taylor and Weisberg (1989)
		and TEMPO software documentation 
		
		Also contains function to read TEMPO .par files to obtain parameters
		and errors on parameters (if available) */

/* Matt Pitkin 29/04/04 */

#ifndef _BINARYPULSARTIMING_H
#define _BINARYPULSARTIMING_H

#include <lal/LALStdlib.h>

/**** DEFINE STRUCTURES ****/

typedef struct
tagBinaryPulsarParams
{
	CHAR *model; 	/* TEMPO binary model e.g. BT, DD, ELL1 */
	REAL8 tb;			/* Time of arrival (TOA) at the SSB */
	CHAR *tbflag;	/* flag is "MJD" if tb in MJD (no leap secs needed)
									 flag is "GPS" if tb in GPS (need to subtract leap secs) */
	UINT4 leapSecs;/* number of leap seconds at GPS time tb */
	REAL8 t0;			/* epoch of data */
	
	REAL8 f0;
	REAL8 f1;
	REAL8 f2;
	REAL8 f3;
	
	REAL8 ra;
	REAL8 dec;
	REAL8 pmra;
	REAL8 pmdec;
	
	REAL8 posepoch;
	REAL8 pepoch;
	
	/* all parameters will be in the same units as used in TEMPO */
	
	/* Keplerian parameters */
	REAL8 e;			/* orbital eccentricity */
	REAL8 Pb;			/* orbital period (days) */
	REAL8 w0;			/* logitude of periastron (deg) */
	REAL8 x;			/* projected semi-major axis/speed of light (light secs) */
	REAL8 T0;			/* time of orbital perisastron as measured in TDB (MJD) */
	
	REAL8 xpbdot;	/* (10^-12) */
	
	/* for low eccentricity orbits (ELL1 model) use Laplace parameters */
	/* (eps1 = e*sin(w), eps2 = e*cos(w)) instead of e, w.						 */
	/* see Appendix A, Ch. Lange etal, MNRAS (2001)										 */
	REAL8 eps1;				/* e*sin(w) */
	REAL8 eps2;				/* e*cos(w) */
	REAL8 eps1dot;
	REAL8 eps2dot;
	REAL8 Tasc;		/* time of the ascending node (used rather than T0) */
	UINT4 nEll;			/* set to zero if have eps time derivs (default)
											 set to 1 if have wdot */
	
	/* Post-Keplarian parameters */
	/* for Blandford-Teukolsky (BT) model */
	REAL8 wdot;		/* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
	REAL8 gamma;	/* gravitational redshift and time dilation parameter (s)*/
	REAL8 Pbdot;	/* rate of change of Pb (dimensionless 10^-12) */
	REAL8 xdot;		/* rate of change of x(=asini/c) - optional (10^-12)*/
	REAL8 edot;		/* rate of change of e (10^-12)*/
	
	/* for Epstein-Haugan (EH) model */
	REAL8 s;			/* Shapiro 'shape' parameter sin i */
	
	/* for Damour-Deruelle (DD) model */
	/*REAL8 r;	 Shapiro 'range' parameter - defined internally as Gm2/c^3 */
	REAL8 dr;			
	REAL8 dth;		/* (10^-6) */
	REAL8 a0, b0;	/* abberation delay parameters */
	
	/* for DD (General Relativity) (DDGR) - assumes GR is correct model */
	/* 		we do not need wdot, gamma, Pbdot, s, r, xdot and edot				*/
	REAL8 M;			/* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
	REAL8 m2;			/* companion mass */
	
	/******** errors read in from a .par file **********/
	REAL8 f0Err;
	REAL8 f1Err;
	REAL8 f2Err;
	REAL8 f3Err;
	
	REAL8 pepochErr;
	REAL8 posepochErr;
	
	REAL8 raErr;
	REAL8 decErr;
	REAL8 pmraErr;
	REAL8 pmdecErr;
	
	REAL8 eErr;			
	REAL8 PbErr;			
	REAL8 w0Err;			
	REAL8 xErr;			
	REAL8 T0Err;			
	
	REAL8 xpbdotErr;	
	
	REAL8 eps1Err;				
	REAL8 eps2Err;				
	REAL8 eps1dotErr;
	REAL8 eps2dotErr;
	REAL8 TascErr;			

	REAL8 wdotErr;		
	REAL8 gammaErr;	
	REAL8 PbdotErr;	
	REAL8 xdotErr;		
	REAL8 edotErr;		
	
	REAL8 sErr;
	
	/*REAL8 rErr;	 Shapiro 'range' parameter - defined internally as Gm2/c^3 */
	REAL8 drErr;			
	REAL8 dthErr;		
	REAL8 a0Err, b0Err;	
	
	REAL8 MErr;			
	REAL8 m2Err;			
}BinaryPulsarParams;

typedef struct
tagBinaryPulsarOutput
{
	REAL8 deltaT;	/* deltaT to add to TDB in order to account for binary */
}BinaryPulsarOutput;

/**** DEFINE FUNCTIONS ****/
void
LALBinaryPulsarDeltaT( LALStatus						*status,
											 BinaryPulsarOutput		*output,
											 BinaryPulsarParams		*params );
											 
void
LALReadTEMPOParFile(	LALStatus							*status,
											BinaryPulsarParams		*output,
											CHAR									*pulsarAndPath );
											
/* define a function to convert RA and Dec on format dd:mm:ss.ss
	 or ddmmss.ss 
	 into the number of degrees as a float
	 degs is the string containing the dd/hh:mm:ss.sss
	 coords is either ra/RA or dec/DEC															*/
double
degsToRads(char *degs, char *coords);
											 
#endif
