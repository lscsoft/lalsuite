/***************************************************************/
/* Version of HeterodynePulsar.h with added binary parameters  */
/* Matt Pitkin 05/05/04                                        */
/***************************************************************/   


/********************************* <lalVerbatim file="HeterodynePulsarHV">
Author: Dupuis, R. J.
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{HeterodynePulsar.h}}

Provides routines to heterodyne, average, and resample the data as required for time domain known pulsar search.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/HeterodynePulsar.h>
\end{verbatim}
The gravitational wave signal from a non-precessing pulsar at twice its frequency can be modeled as

\begin{equation}
h(t) = F_{+}(t;\psi)h_{0}(1 + \cos^{2}\iota)\cos 2 \Psi(t) + 2 F_{\times} h_{0}\cos \iota \sin 2 \Psi(t)
\end{equation}

where $F_{+,\times}$ are the amplitude responses of the detectors, $\psi$ is the polarization angle,
$\iota$ describes the inclination of the pulsar with respect to the line of sight, and $\Psi(t) =
\phi_{0} + \phi(t)$ describes the phase of the pulsar.

The phase $\Psi(t)$ of the pulsar is calculated as 
\begin{equation}
\Psi(t) = \phi_{0} +2\pi \left(f_{0}(T - T_{0}) + 
\frac{1}{2}\dot{f_{0}} (T - T_{0})^{2} + 
\frac{1}{6}\ddot{f_{0}}(T - T_{0})^{3}\right)
\end{equation}
where 
\begin{equation}
T = t + \delta t= t + \frac{\vec{r} \cdot \vec{n}}{c}  + \Delta_{E\odot}
\end{equation}
where {\emph T} is the time in a frame inertial with respect to the pulsar and $\phi_{0}$ is the phase of at time $T_{0}$.
 The time difference $\delta t$ due to the motion of the earth in the solar system is calculated using
 \texttt{LALBarycenter()}.


More documentation later.

\vfill{\footnotesize\input{HeterodynePulsarHV}}
\newpage\input{HeterodynePulsarC}
\newpage\input{HeterodynePulsarTestC}

********************************** </lalLaTeX> */

#ifndef _HETERODYNEPULSAR_H
#define _HETERODYNEPULSAR_H

#include <lal/LALStdlib.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h> 
#include <lal/AVFactories.h>

/* include the BinaryPulsarTiming.h header (only local in my directory at the mo*/
#include "BinaryPulsarTiming.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (HETERODYNEPULSARH, "$Id$");

/******************************** <lalErrTable file="HeterodynePulsarHE"> */
#define HETERODYNEPULSARH_ENULLINPUT 1
#define HETERODYNEPULSARH_ENULLOUTPUT 2
#define HETERODYNEPULSARH_ENULLPARAMS 3
#define HETERODYNEPULSARH_ERLENGTH 4

#define HETERODYNEPULSARH_MSGENULLINPUT "Input was Null"
#define HETERODYNEPULSARH_MSGENULLOUTPUT "Output was Null"
#define HETERODYNEPULSARH_MSGENULLPARAMS "Params was Null"
#define HETERODYNEPULSARH_MSGERLENGTH "Vector has wrong length"

/************************************ </lalErrTable> */

/****** DEFINE NEW STRUCTURES AND TYPES ************/
typedef struct
tagCalibrateFineInput
{
   COMPLEX16 	B;  /* uncalibrated band limited data */  
   COMPLEX16 var;
   REAL8 alpha; /* alpha value at time tgps */
   REAL8 beta; /* beta value at time tgps */
} CalibrateFineInput;

typedef struct
tagCalibrateFineOutput
{
    COMPLEX16 	B;  /* calibrated fine heterodyne output */  
    COMPLEX16   var;
    COMPLEX16 	R;  /* response of detector at frequency f */
} CalibrateFineOutput;

typedef struct
tagCalibrateFineParams
{
    REAL8 phaseH0; /* open loop gain phase at f*/
    REAL8 H0; /*  open loop gain amplitude at f*/
    REAL8 phaseC0; /* sensing function phase at f*/
    REAL8 C0;  /* sensing function amplitude at f*/
} CalibrateFineParams;


typedef struct
tagHeterodyneInput
{
  REAL8TimeSeries   V;             /* raw data */
  REAL8			f0;             /* frequency of the signal */
  REAL8			f1;             /* first time derivative of frequency */
  REAL8			f2;             /* second time derivative of frequency */
  REAL8 		fEpochGPS;      /* epoch of the frequency at SSB */
  SkyPosition           source;         /* location of pulsar in sky - equatorial coordinate system */
  REAL8 		pmRA;           /* proper motion RA radians / year */
  REAL8			pmDEC;		/* proper motion DEC radians / year*/
  REAL8 		posEpochGPS;    /* epoch of RA and DEC */
} HeterodyneInput;

typedef struct
tagHeterodyneOutput
{
  COMPLEX16 	B;		/* bin value */
  COMPLEX16 	var;  		/* variance -- not currently used */
} HeterodyneOutput;

typedef struct
tagHeterodyneParams
{  
  EphemerisData *edat;
  LALDetector detector;
  REAL8IIRFilter        *iirFilter1Re;   /* IIR filter to be applied to real part of complex heterodyned data */
  REAL8IIRFilter        *iirFilter1Im;   /* IIR filter to be applied to imaginary part of complex heterodyned data */
  REAL8IIRFilter        *iirFilter2Re;   /* IIR filter to be applied to real part of complex heterodyned data */
  REAL8IIRFilter        *iirFilter2Im;   /* IIR filter to be applied to imaginary part of complex heterodyned data */  
  REAL8IIRFilter        *iirFilter3Re;   /* IIR filter to be applied to real part of complex heterodyned data */
  REAL8IIRFilter        *iirFilter3Im;   /* IIR filter to be applied to imaginary part of complex heterodyned data */
	
	/* adding binary structure to the HeterodyneParams structure */
	BinaryPulsarParams binaryParams;
	BinaryPulsarInput binaryInput;
} HeterodyneParams;


/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/	     

void
LALCalibrateFineHeterodyne ( LALStatus                      *status,
		  	    CalibrateFineOutput		    *output,
		   	    CalibrateFineInput 	   	    *input,
		   	    CalibrateFineParams	            *params );

void
LALHeterodyneToPulsar       ( LALStatus                   *status,
		  	    HeterodyneOutput        *output,
		   	    HeterodyneInput         *input,
		   	    HeterodyneParams        *params );	
#ifdef  __cplusplus
}
#endif

#endif /* _HETERODYNEPULSAR_H */
