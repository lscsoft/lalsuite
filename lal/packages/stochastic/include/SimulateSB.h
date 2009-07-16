/*
*  Copyright (C) 2007 Jolien Creighton, Tania Regimbau, Teviet Creighton, John Whelan
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

/*********************** <lalVerbatim file="SimulateSBHV">
Author: Sukanta Bose
$Id$
*********************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{SimulateSB.h}}
\label{inject:s:SimulateSB.h}

Provides prototype and error code information for the modules needed
to simulate a stochastic background signal (whitened, if desired) in a pair of
detectors, given the appropriate representations of the
detector transfer function in each detector.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SimulateSB.h>
\end{verbatim}

\noindent

\subsection*{Error conditions}
\input{SimulateSBHE}

\subsection*{Structures}

*********************************************************** </lalLaTeX> */

#ifndef _SIMULATESB_H
#define _SIMULATESB_H


#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/RealFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

  NRCSID( SIMULATESBH, "$Id$" );

/***************** <lalErrTable file="SimulateSBHE"> */

#define SIMULATESBH_ENULLP          1
#define SIMULATESBH_ENONPOSLEN      2
#define SIMULATESBH_ENONPOSDELTAF   3
#define SIMULATESBH_ENONPOSDELTAT   4
#define SIMULATESBH_ENEGFMIN        5
#define SIMULATESBH_EMMTIME         6
#define SIMULATESBH_EMMHETERO       7
#define SIMULATESBH_EMMFMIN         8
#define SIMULATESBH_EMMDELTAF       9
#define SIMULATESBH_EMMLEN         10
#define SIMULATESBH_EOORFREF       11
#define SIMULATESBH_ENONPOSOMEGA   12
#define SIMULATESBH_EALOC          13
#define SIMULATESBH_ENONZEROHETERO 14
#define SIMULATESBH_EWRONGUNITS    15
#define SIMULATESBH_ECOMPTIME      16
#define SIMULATESBH_ENOTYETHETERO 255

#define SIMULATESBH_MSGENULLP         "Null pointer"
#define SIMULATESBH_MSGENONPOSLEN     "Negative or zero length for data member of time series"
#define SIMULATESBH_MSGENONPOSDELTAF  "Negative or zero frequency spacing"
#define SIMULATESBH_MSGENONPOSDELTAT  "Negative or zero time spacing"
#define SIMULATESBH_MSGENEGFMIN       "Negative start frequency"
#define SIMULATESBH_MSGEMMTIME        "Mismatch in epochs"
#define SIMULATESBH_MSGEMMHETERO      "Mismatch in heterodyning frequencies"
#define SIMULATESBH_MSGEMMFMIN        "Mismatch in start frequencies"
#define SIMULATESBH_MSGEMMDELTAF      "Mismatch in frequency spacings"
#define SIMULATESBH_MSGEMMLEN         "Mismatch in sequence lengths"
#define SIMULATESBH_MSGEOORFREF       "Out of range reference frequency"
#define SIMULATESBH_MSGENONPOSOMEGA   "Negative stochastic background strength"
#define SIMULATESBH_MSGEALOC         "Memory allocation error"
#define SIMULATESBH_MSGENONZEROHETERO "Non-zero heterodyning frequency specified for real time series"
#define SIMULATESBH_MSGEWRONGUNITS    "Inconsistent input units"
#define SIMULATESBH_MSGECOMPTIME      "Time domain data complex instead of real"
#define SIMULATESBH_MSGENOTYETHETERO  "Non-zero heterodyning frequency not yet implemented"

/************************************ </lalErrTable> */

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             SimulateSB.c                  *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures associated with
  \texttt{SimulateSB.c}
  (Sec.~\ref{inject:ss:SimulateSB.c})}

\subsubsection*{\texttt{struct SSSimStochBGOutput}}
\idx[Type]{SSSimStochBGOutput}

\noindent Contains the output data produced by
\texttt{LALSSSimStochBGTimeSeries()}. It comprises of a pair of
(real) time-series simulated stochastic background signal in the outputs of
a given pair of detectors. The fields are:

\begin{description}
\item[\texttt{REAL4TimeSeries *SSimStochBG1}]
Simulated stochastic background signal in the output of
the first detector.
\item[\texttt{REAL4TimeSeries *SSimStochBG2}]
Simulated stochastic background signal in the output of
the second detector.
\end{description}

*********************************************************** </lalLaTeX> */

  typedef struct tagSSSimStochBGOutput {
    REAL4TimeSeries    *SSimStochBG1;
    REAL4TimeSeries    *SSimStochBG2;
  } SSSimStochBGOutput;

/*********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct SSSimStochBGInput}}
\idx[Type]{SSSimStochBGInput}

\noindent Contains the input data needed by
\texttt{LALSSSimStochBGTimeSeries()}
to calculate the whitened stochastic background signal in the output of
a detector.
The fields are:

\begin{description}
\item[\texttt{REAL4FrequencySeries *omegaGW}] The spectrum
$\Omega_{\scriptstyle{\rm GW}}(f)$ of the stochastic gravitational-wave
background.
\item[\texttt{COMPLEX8FrequencySeries *whiteningFilter1}]
The frequency-domain response function $\tilde{R}_1(f)$ for the first detector.
\item[\texttt{COMPLEX8FrequencySeries *whiteningFilter2}]
The frequency-domain response function $\tilde{R}_2(f)$ for the second detector.
\end{description}

*********************************************************** </lalLaTeX> */

  typedef struct tagSSSimStochBGInput {
    REAL4FrequencySeries     *omegaGW;
    COMPLEX8FrequencySeries  *whiteningFilter1;
    COMPLEX8FrequencySeries  *whiteningFilter2;
  } SSSimStochBGInput;

  typedef struct tagSSSimStochBGStrainInput {
    REAL4FrequencySeries     *omegaGW;
  } SSSimStochBGStrainInput;


/*********************************************************** <lalLaTeX>


\subsubsection*{\texttt{struct SSSimStochBGParams}}
\idx[Type]{SSSimStochBGParams}

\noindent Contains the parameters used by \texttt{LALSSSimStochBGTimeSeries()}
to compute the whitened stochastic background signal in the output of an
interferometric detector. The fields are:

\begin{description}
\item[\texttt{UINT4 length}]
The number of points in the output time series.

\item[\texttt{REAL8 deltaT}]
The temporal spacing of the output time series.

\item[\texttt{INT4 seed}]
The random number seed for the stochastic simulation.

\item[\texttt{LALDetector *detector1}]
The site location and orientation information of first detector involved in
the stochastic background search.

\item[\texttt{LALDetector *detector2}]
The site location and orientation information of second detector involved in
the stochastic background search.

\item[\texttt{LALUnit SSimStochBGTimeSeries1Unit}]
The unit field of the stochastic background, expressed as a Real4
time series, in detector 1.

\item[\texttt{LALUnit SSimStochBGTimeSeries2Unit}]
The unit field of the stochastic background, expressed as a Real4
time series, in detector 2.

\end{description}
*********************************************************** </lalLaTeX> */

  typedef struct tagSSSimStochBGParams {
    UINT4        length;   /* time length of output vector data samples */
    REAL8        deltaT;   /* time spacing */
    INT4         seed;     /* for random numbers x, y */
    LALDetector  detectorOne;
    LALDetector  detectorTwo;
    LALUnit      SSimStochBGTimeSeries1Unit;
    LALUnit      SSimStochBGTimeSeries2Unit;
  } SSSimStochBGParams;

  typedef struct tagSSSimStochBGStrainParams {
    UINT4        length1,length2;   /* time length of output vector data samples */
    REAL8        deltaT1, deltaT2;   /* time spacing */
    INT4         seed;     /* for random numbers x, y */
    LALDetector  detectorOne;
    LALDetector  detectorTwo;
    LALUnit      SSimStochBGTimeSeries1Unit;
    LALUnit      SSimStochBGTimeSeries2Unit;
  } SSSimStochBGStrainParams;


  void
  LALSSSimStochBGTimeSeries( LALStatus                  *status,
			     SSSimStochBGOutput           *output,
			     SSSimStochBGInput            *input,
			     SSSimStochBGParams           *params );

  void
  LALSSSimStochBGStrainTimeSeries( LALStatus                  *status,
			     SSSimStochBGOutput           *output,
			     SSSimStochBGStrainInput            *input,
			     SSSimStochBGStrainParams           *params );

#ifdef  __cplusplus
}
#endif

#endif /* _SIMULATESB_H */

/********************************************************** <lalLaTeX>

\vfill{\footnotesize\input{SimulateSBHV}}

\newpage\input{SimulateSBC}
%\newpage\input{SimulateSBTestC}
*********************************************************** </lalLaTeX> */
