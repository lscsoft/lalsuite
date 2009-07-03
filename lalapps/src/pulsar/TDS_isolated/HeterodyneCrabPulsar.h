/*
*  Copyright (C) 2007 Matt Pitkin
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

/* <lalVerbatim file="HeterodyneCrabPulsarHV">
	 Author: Pitkin, M. D.
	 $Id$
	 </lalVerbatim> */
	 
/* <lalLaTeX>
	 \section{Header \texttt{HeterodyneCrabPulsar.h}}
	 label{ss:HeterodyneCrabPulsar.h}
	
	 Calculates and heterodynes the timing noise component of the Crab pulsar
	 phase.
	 
	 \subsection*{Synopsis}
	 \begin{verabtim}
	 #include <lal/HeterodyneCrabPulsar.h>
	 \end{verbatim}
	 
	 \noindent This header covers routines for reading in the Crab pulsar
	 ephemeris, calculating the spline interpolation of the phase between
	 ephemeris points and heterodyning data to remove the phase difference caused
	 by timing noise.
	 </lalLaTeX> */
	 
/* Matt Pitkin 10/03/04 - TNHeterodyne struct changed */

/* Matt Pitkin 26/03/04 - function arguments have been changed to conform to the
   LAL Spec */

#ifndef _HETERODYNECRABPULSAR_H
#define _HETERODYNECRABPULSAR_H

/* lal headers */
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h> 
#include <lal/DetectorSite.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/LALRCSID.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>

#ifdef __cplusplus
extern "C" {
#endif

/* DEFINE RCS ID STRING */
NRCSID (HETERODYNECRABPULSARH, "$Id$"); 

/* <lalLaTeX>
	 \subsection*{Error conditions}
	 </lalLaTeX> */

/* <lalErrTable> */
#define HETERODYNECRABPULSARH_ENULLINPUT 1
#define HETERODYNECRABPULSARH_ENULLOUTPUT 2
#define HETERODYNECRABPULSARH_EEPHEMERISFILENAME 3
#define HETERODYNECRABPULSARH_ENUMEPHEMERISDATA 4
#define HETERODYNECRABPULSARH_EINVALIDF0 5

#define HETERODYNECRABPULSARH_MSGENULLINPUT "Input was Null"
#define HETERODYNECRABPULSARH_MSGENULLOUTPUT "Output was Null"
#define HETERODYNECRABPULSARH_MSGEEPHEMERISFILENAME "No ephemeris filename given"
#define HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA "Number of ephemeris data points must be greater than zero"
#define HETERODYNECRABPULSARH_MSGEINVALIDF0 "F0 must be greater than 0"
/* </lalErrTable> */

/* <lalLaTex>
	 \subsection*{Structures}
	 \idx[Type]{GetCrabEphemerisInput}
	 \idx[Type]{CrabSpindownParamsInput}
	 \idx[Type]{CrabSpindownParamsOutput}
	 \idx[Type]{ParamsForHeterodyne}
	 \idx[Type]{TNHeterodyneInput}
	 \idx[Type]{TNHeterodyneOutput}
	 
	 \begin{verbatim}
	 typedef struct tagGetCrabEphemerisInput GetCrabEphemerisInput;
	 \end{verbatim}
	 
	 This structure contains the file name of the Crab ephemeris, nominally
	 \verb+crab_ephemeris.txt+.
	 
	 \begin{verbatim}
	 typedef struct tagCrabSpindownParamsInput CrabSpindownParamsInput;
	 \end{verbatim}
	 
	 This structure contains vectors of the values of time of arrival, $f_0$ (Hz)
	 and $\dot{f_0}$ ($10^{-10}$\,{\rm Hz}^2) extracted from the Crab ephemeris 
	 file.
	 
	 \begin{verbatim}
	 typedef struct tagCrabSpindownParamsOutput CrabSpindownParamsOutput;
	 \end{verbatim}
	 
	 This structure contains vectors of the values of time of arrival, $f_0$ (Hz)
	 and $\dot{f_0}$ (${\rm Hz}^2$) from the ephemeris file and values 
	 of $f$ time derivatives up to fourth order interpolated from the ephemeris 
	 values.
	 
	 \begin{verbatim}
	 typedef struct tagParamsForHeterodyne ParamsForHeterodyne;
	 \end{verbatim}
	 
	 This structure contains the single values of $f$ and its time derivative up 
	 to fourth order and their epoch.
	 
	 \begin{verbatim}
	 typedef struct tagTNHeterodyneInput TNHeterodyneInput;
	 \end{verbatim}
	 
	 This structure contains values of $f$, $\dot{f}$ and $\ddot{f}$ used in the
	 initial Crab heterodyne, the complex heterodyned data point, the epoch of the
	 intial heterodyne parameters $t_0$, and the \verb+LIGOTimeGPS+ time of the
	 data point.
	 
	 \begin{verbatim}
	 typedef struct tagTNHeterodyneOutput TNHeterodyneOutput;
	 \end{verbatim}

	 This structure contains the final complex heterodyned output and the phase
	 difference between the initial heterodyne and the timing noise heterodyne.
	 
	 </lalLaTeX> */

/* define new structures and type for crab ephemeris reading - to be put in header file later */
typedef struct
tagGetCrabEphemerisInput
{
	CHAR *filename; /* file name of the Crab ephemeris (crab_ephemeris.txt) */
} GetCrabEphemerisInput;

typedef struct
tagCrabSpindownParamsInput
{
  REAL8Vector *tArr;	/* vector of pulse arrival times from Crab ephemeris (GPS seconds) */
  REAL8Vector *f0;	/* vector of pulsar frequencies at each t_arr */
  REAL8Vector *f1;	/* vector of first frequecy derivs at each t_arr */
	UINT4 numOfData;  /* num of data in ephemeris file */
} CrabSpindownParamsInput;

typedef struct
tagCrabSpindownParamsOutput
{
  REAL8Vector *tArr;	/* vector of pulse arrival times from Crab ephemeris (GPS seconds) */
  REAL8Vector *f0;	/* vector of pulsar frequencies at each t_arr */
  REAL8Vector *f1;	/* vector of first frequecy derivs at each t_arr */
  REAL8Vector *f2;	/* vectors of second, third and fourth frequency */
  REAL8Vector *f3;	/* derivs of the Crab pulsar calculated at each  */
  REAL8Vector *f4;	/* t_arr, length will be one less than t_arr in input    */
	UINT4 numOfData;  /* num of data in ephemeris file */
} CrabSpindownParamsOutput;

/*typedef struct
tagParamsForHeterodyne
{
  REAL8 f0;
  REAL8 f1;
  REAL8 f2;
  REAL8 f3;
  REAL8 f4;
  REAL8 epoch;
  UINT4 startPos;
  EphemerisData *edat;
  LALDetector detector;
} ParamsForHeterodyne;*/

/* ParamsForHeterodyne struct changed */
typedef struct
tagParamsForHeterodyne
{
  REAL8 f0;
  REAL8 f1;
  REAL8 f2;
  REAL8 f3;
  REAL8 f4;
  REAL8 epoch;
} ParamsForHeterodyne;

/*typedef struct
tagTNHeterodyneInput
{
  REAL8			f0;
  REAL8			f1;
  REAL8			f2;
  COMPLEX16TimeSeries	Vh;
  REAL8			t0;
  REAL8			phase;
  SkyPosition           source;  
} TNHeterodyneInput;*/

/*TNHeterodyne structure that doesn't need to heterodyne at the SSB */
typedef struct
tagTNHeterodyneInput
{
  REAL8			f0;		/* the *initial* (i.e. that performed prior to the timing
	noise heterodyne) heterodyne frequency */
  REAL8			f1;		/* the initial heterdyne f1 */
  REAL8			f2;		/* the initial heterodyne f2 */
  REAL8                 f3;             /* the initial heterodyne f3 */
  COMPLEX16 Vh;		/* the complex output of the initial heterodyne */
	REAL8			t0;		/* the epoch of the initial heterodyne frequency */
	LIGOTimeGPS epoch;  /* the time of the data point Vh */
} TNHeterodyneInput; 

typedef struct
tagTNHeterodyneOutput
{
  COMPLEX16	Vh;
	REAL8 Dphase;
	REAL8 phi0;
	REAL8 phi1;
} TNHeterodyneOutput;

/* function definitions */

/* <lalLaTeX>
	 \newpage\input{HeterodyneCrabPulsarC}
	 </lalLaTeX> */

void 
LALGetCrabEphemeris	( LALStatus			*status,
			  CrabSpindownParamsInput	*output,
				GetCrabEphemerisInput *input );
			  
void
LALComputeFreqDerivatives	( LALStatus			*status,
				  CrabSpindownParamsOutput	*ouput,
					CrabSpindownParamsInput	*input );			  

void
LALSetSpindownParams	( LALStatus			*status,
			  ParamsForHeterodyne		*output,
				CrabSpindownParamsOutput	*input,
				LIGOTimeGPS			epoch );

void
LALTimingNoiseHeterodyne	( LALStatus		*status,
				  TNHeterodyneOutput	*output,
					TNHeterodyneInput	*input,
				  ParamsForHeterodyne	*params,
                                  EphemerisData *edat,
                                  BarycenterInput baryinput,
                                  EarthState earth );

#ifdef __cplusplus
}
#endif
				  
#endif /* _HETERODYNECRABPULSAR_H */
