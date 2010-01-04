/*
*  Copyright (C) 2007 Jolien Creighton
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

/********************************* <lalVerbatim file="FitToPulsarHV">
Author: Dupuis, R.J.
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{FitToPulsar.h}}
Provides routines for finding the best fit of the measured data to the strain expected from
non-precessing pulsar.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/FitToPulsar.h>
\end{verbatim}

The model to be fitted to the data after \texttt{LALFineHeterodyneToPulsar} has been applied is
\begin{equation}
y(t;\textrm{{\bf a}}) = F_{+}(t;\psi)h_{0} (1 + \cos^{2}\iota)e^{i2\phi_{0}} - 2iF_{\times}(t;\psi) h_{0} \cos\iota e^{i2\phi_{0}}
\end{equation}


The reduced set of data points is fitted to this model by minimizing
 $\chi^2$ over $h_{0}$, $\phi_{0}$, $\iota$, and $\psi$.

\begin{equation}
\chi^2(\textrm{{\bf a}}) = \sum_{k}\left|\frac{B_{k} - y(t;\textrm{{\bf a}})}{\sigma_{k}^{2}}\right|^2
\end{equation}

The minimization of $\chi^2$ is done in two steps \texttt{LALCoarseFitToPulsar()} and \texttt{LALFineFitToPulsar()}.


More documentation soon.
\subsection*{Error conditions}
\input{FitToPulsarHE}

\subsection*{Structures}

\subsubsection*{Structure \texttt{CoarseFitInput}}
\idx[Type]{CoarseFitInput}
\noindent This structure stores locked data to be fitted by model.
\begin{description}
\item[\texttt{COMPLEX16Vector *B}] heterodyned, averaged and resampled data
\item[\texttt{COMPLEX16Vector *var}] variance of the rFactor points that were averaged
\item[\texttt{LIGOTimeGPS *t}] time stamp for each data point (not necessarily with equal time steps)
\end{description}

\subsubsection*{Structure \texttt{CoarseFitOutput}}
\idx[Type]{CoarseFitOutput}
\noindent This structure stores the results from the coarse fit of parameters.
\begin{description}
\item[\texttt{REAL8 h0}] best fit h0
\item[\texttt{REAL8 eh0[3]}] standard error for h0, min standard error, max standard error
\item[\texttt{REAL8 cosIota}] best fit cosIota
\item[\texttt{REAL8 phase}] best fit phase
\item[\texttt{REAL8 psi}] best fit psi (polarization angle)
\item[\texttt{REAL8 chiSquare}] min value of chi square
\item[\texttt{REAL8Vector *mChiSquare}] matrix with chi square values
\end{description}

\subsubsection*{Structure \texttt{CoarseFitParams}}
\idx[Type]{CoarseFitParams}
\noindent This structure stores the parameters for the coarse fit.
\begin{description}
\item[\texttt{REAL8 meshH0[3]}] min h0, delta h0, number of steps
\item[\texttt{REAL8 meshCosIota[3]}] min cosIota, delta cosIota, number of steps
\item[\texttt{REAL8 meshPhase[3]}] min phase, delta phase, number of steps
\item[\texttt{REAL8 meshPsi[3]}] min psi, delta psi, number of steps
\item[\texttt{LALSource pulsarSrc}] describes sky position of pulsar
\item[\texttt{LALDetector detector}] detector
\end{description}

\vfill{\footnotesize\input{FitToPulsarHV}}
\newpage\input{FitToPulsarC}
\newpage\input{FitToPulsarTestC}

********************************** </lalLaTeX> */

#ifndef _FITTOPULSAR_H
#define _FITTOPULSAR_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/Units.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (FITTOPULSARH, "$Id$");

/******************************** <lalErrTable file="FitToPulsarHE"> */
#define FITTOPULSARH_ENULLINPUT 1
#define FITTOPULSARH_ENULLOUTPUT 2
#define FITTOPULSARH_ENULLPARAMS 3
#define FITTOPULSARH_EVECSIZE 4
#define FITTOPULSARH_EMESH 5
#define FITTOPULSARH_EVAR 6
#define FITTOPULSARH_EMAXCHI 7
#define FITTOPULSARH_EDIVZERO 8

#define FITTOPULSARH_MSGENULLINPUT "Input was Null"
#define FITTOPULSARH_MSGENULLOUTPUT "Output was Null"
#define FITTOPULSARH_MSGENULLPARAMS "Params was Null"
#define FITTOPULSARH_MSGEVECSIZE "Input vectors were not the same length"
#define FITTOPULSARH_MSGEMESH "Mesh paramaters supplied were invalid"
#define FITTOPULSARH_MSGEVAR "Variance vector in Input had invalid values"
#define FITTOPULSARH_MSGEMAXCHI "The minimum value of chiSquare was greater than INICHISQU"
#define FITTOPULSARH_MSGEDIVZERO "Attempted to divide by zero"

/************************************ </lalErrTable> */

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

/****** DEFINE NEW STRUCTURES AND TYPES ************/
typedef struct
tagCoarseFitInput
{
  COMPLEX16Vector *B;     /* heterodyned, averaged and resampled data */
  COMPLEX16Vector *var;   /* variance of the rFactor points that were averaged */
  LIGOTimeGPS *t;        /* time stamp for each data point (not necessarily with equal time steps)*/

} CoarseFitInput;

typedef struct
tagCoarseFitOutput
{
  REAL8 h0;              /* best fit h0 */
  REAL8	eh0[3];		 /* standard error for h0, min standard error, max standard error */
  REAL8 cosIota;         /* best fit cosIota */
  REAL8 phase;           /* best fit phase */
  REAL8 psi;             /* best fit psi (polarization angle) */
  REAL8 chiSquare;       /* value of min chi square */
  REAL8Vector *mChiSquare;  /* matrix with chi square values*/
} CoarseFitOutput;

typedef struct
tagCoarseFitParams
{
  REAL8 meshH0[3];	  /* min h0, delta h0, number of steps */
  REAL8 meshCosIota[3];   /* min cosIota, delta cosIota, number of steps */
  REAL8 meshPhase[3];     /* min phase, delta phase, number of steps */
  REAL8 meshPsi[3];       /* min psi, delta psi, number of steps */
  LALSource pulsarSrc;    /* describes sky position of pulsar */
  LALDetector detector;   /* detector */
} CoarseFitParams;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALCoarseFitToPulsar	(	LALStatus                       *status,
		   		CoarseFitOutput                 *output,
		   		CoarseFitInput 			*input,
		   		CoarseFitParams                 *params);

#ifdef  __cplusplus
}
#endif

#endif /* _FITTOPULSAR_H */
