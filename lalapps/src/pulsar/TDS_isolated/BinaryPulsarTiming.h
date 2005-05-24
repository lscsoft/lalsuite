/* LAL code to take into account binary pulsar motion */
/* 	for definitions of different models see Taylor and Weisberg (1989)
		and TEMPO software documentation 
		
		Also contains function to read TEMPO .par files to obtain parameters
		and errors on parameters (if available) */
/* <lalVerbatim file="BinaryPulsarTimingHV">
	 Author: Pitkin, M. D.
	 $Id$
*/
/* Matt Pitkin 29/04/04 */

/* <lalLaTeX>
	 \section{Header \texttt{BinaryPulsarTiming.h}}
	 label{ss:BinaryPulsarTiming.h}
	 
	 Calculates time delay to a signal from a pulsar in a binary system.
	 
	 \subsection*{Synopsis}
	 \begin{verbatim}
	 #include <lal/BinaryPulsarTiming.h>
	 \end{verbatim}
	 
	 \noindent This header covers routines for calculating the time delay to a
	 signal from a pulsar in a binary system. The are also routines for reading
	 pulsar data from TEMPO .par file formats.
	 </lalLaTeX> */

#ifndef _BINARYPULSARTIMING_H
#define _BINARYPULSARTIMING_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* <lalLaTeX>
	 \subsection*{Error conditions}
	 </lalLaTeX> */
	 
/* <lalErrTable> */
#define BINARYPULSARTIMINGH_ENULLINPUT 1
#define BINARYPULSARTIMINGH_ENULLOUTPUT 2
#define BINARYPULSARTIMINGH_ENULLPARAMS 3
#define BINARYPULSARTIMINGH_ENULLBINARYMODEL 4
#define BINARYPULSARTIMINGH_ENULLTBFLAG 5
#define BINARYPULSARTIMINGH_ENULLPULSARANDPATH 6
#define BINARYPULSARTIMINGH_EPARFILEERROR 7

#define BINARYPULSARTIMINGH_MSGENULLINPUT "Input was Null"
#define BINARYPULSARTIMINGH_MSGENULLOUTPUT "Output was Null"
#define BINARYPULSARTIMINGH_MSGENULLPARAMS "Params was Null"
#define BINARYPULSARTIMINGH_MSGNULLBINARYMODEL "Binary model is Null or not specified - you should not be in the binary timing routine" 
#define BINARYPULSARTIMINGH_MSGENULLTBFLAG "tbflag is Null"
#define BINARYPULSARTIMINGH_MSGENULLPULSARANDPATH "Path to pulsar.par file not specified"
#define BINARYPULSARTIMINGH_MSGEPARFILEERROR ".par file path is wrong or file doesn't exist"
/* </lalErrTable> */

/* <lalLaTeX>
	 \subsection*{Structures}
	 \idx[Type]{BinaryPulsarParams}
	 \idx[Type]{BinaryPulsarInput}
	 \idx[Type]{BinaryPulsarOutput}
	 
	 \begin{verbatim}
	 typedef struct tagBinaryPulsarParams BinaryPulsarParams;
	 \end{verbatim}
	 
	 This structure contains all the pulsars parameters. The structure does not
	 have to be used for a binary pulsar, but can just contain the parameters for
	 an isolated pulsar. All parameters are in the same units as given by TEMPO.
	 
	 \begin{verbatim}
	 typedef struct tagBinaryPulsarInput BinaryPulsarInput;
	 \end{verbatim}
	 
	 This structure contains the input time at which the binary correction needs
	 to be calculated.
	 
	 \begin{verbatim}
	 typedef struct tagBinaryPulsarOutput BinaryPulsarOutput;
	 \end{verbatim}
	 
	 This structure contains the binary time delay for the input time.
	 
	 </lalLaTeX> */

/**** DEFINE STRUCTURES ****/

typedef struct
tagBinaryPulsarParams
{
	CHAR *name;		/* pulsar name */
	
	CHAR *model; 	/* TEMPO binary model e.g. BT, DD, ELL1 */
	/*REAL8 tb;			 Time of arrival (TOA) at the SSB 
	CHAR *tbflag;	/* flag is "MJD" if tb in MJD (no leap secs needed)
									 flag is "GPS" if tb in GPS (need to subtract leap secs) 
	UINT4 leapSecs;/* number of leap seconds at GPS time tb 
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
tagBinaryPulsarInput
{
	REAL8 tb;			/* Time of arrival (TOA) at the SSB */
	CHAR *tbflag;	/* flag is "MJD" if tb in MJD (no leap secs needed)
									 flag is "GPS" if tb in GPS (need to subtract leap secs) */
	UINT4 leapSecs;/* number of leap seconds at GPS time tb */
	REAL8 t0;			/* epoch of data */
}BinaryPulsarInput;

typedef struct
tagBinaryPulsarOutput
{
	REAL8 deltaT;	/* deltaT to add to TDB in order to account for binary */
}BinaryPulsarOutput;

/* <lalLaTeX>
	 \newpage\input{BinaryPulsarTimingC}
	 </lalLaTeX> */

/**** DEFINE FUNCTIONS ****/
void
LALBinaryPulsarDeltaT( LALStatus						*status,
											 BinaryPulsarOutput		*output,
											 BinaryPulsarInput		*input,
											 BinaryPulsarParams		*params );
											 
void
LALReadTEMPOParFile(	LALStatus							*status,
											BinaryPulsarParams		*output,
											CHAR									*pulsarAndPath );
											
/* define a function to convert RA and Dec in format dd:mm:ss.ss
	 or ddmmss.ss into the number of degrees as a float
	 degs is the string containing the dd/hh:mm:ss.sss
	 coords is either ra/RA or dec/DEC															*/
REAL8
LALDegsToRads(CHAR *degs, CHAR *coords);

/* function for converting times given in Terrestrial time (TT) or TDB in MJD to
times in GPS - this is important for epochs given in .par files which are in
TDB. TT and GPS are different by a factor of 51.184 secs, this is just the
historical factor of 32.184 secs between TT and TAI (International Atomic Time)
and the other 19 seconds come from the leap seonds added between the TAI and
UTC up to the point of definition of GPS time at UTC 01/01/1980 (see
http://www.stjarnhimlen.se/comp/time.html for details) */
REAL8
LALTTtoGPS(REAL8 TT);

#ifdef __cplusplus
}
#endif
											 
#endif /* _BINARYPULSARTIMING_H */
