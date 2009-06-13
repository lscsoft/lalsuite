/*
*  Copyright (C) 2007 Jolien Creighton, Kaice T. Reilly, Tania Regimbau, John Whelan
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

/*********************** <lalVerbatim file="StochasticCrossCorrelationHV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)
$Id$
*********************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{StochasticCrossCorrelation.h}}
\label{stochastic:s:StochasticCrossCorrelation.h}

Provides prototype and error code information for the modules needed
to calculate the standard optimally-filtered cross-correlation
statistic for stochastic background searches, given a pair of data
segments, along with appropriate representations of the detector
transfer function and the (uncalibrated) power spectral density of the
noise in each detector.  The relationship among these modules is
illustrated in Fig.~\ref{stochastic:f:CrossCorrFlowchart}.

\begin{figure}[htb!]
\begin{center}
\begin{picture}(410,250)(-32,-110)
\put(168,123){\vector(1,0){15}}
\put(168,123){\line(1,0){25}}
\put(195,120){$Y$}
\put(65,110)
{
  \framebox(100,30)
  {
  \texttt{CrossCorr}
  }
}
\put(80,85){\vector(0,1){15}}
\put(80,85){\line(0,1){25}}
\put(82,95){$\widetilde{\bar{h}}{}_{1,2}$}
\put(-37,65){$h_{1,2}$}
\put(-17,68){\vector(1,0){15}}
\put(-17,68){\line(1,0){25}}
\put(5,55)
{
  \framebox(90,30)
  {
    \texttt{ZeroPadAndFFT}
  }
}
\put(150,85){\vector(0,1){15}}
\put(150,85){\line(0,1){25}}
\put(152,95){$\widetilde{Q}$}
\put(100,55)
{
  \framebox(190,30)
  {
    \texttt{OptimalFilter}
  }
}
\put(260,30){\vector(0,1){15}}
\put(260,30){\line(0,1){25}}
\put(262,40){$\lambda$}
\put(343,13){\vector(1,0){15}}
\put(343,13){\line(1,0){25}}
\put(370,10){$\sigma^2/T$}
\put(190,0)
{
  \framebox(150,30)
  {
    \texttt{Normalization}
  }
}
\put(200,-30){\vector(0,1){25}}
\put(200,-30){\line(0,1){30}}
\put(202,-10){$\Omega_{\scriptstyle{\rm GW}}$}
\put(200,-10){\line(-1,0){35}}
\put(165,-10){\line(0,1){65}}
\put(165,-10){\vector(0,1){55}}
\put(167,40){$\Omega_{\scriptstyle{\rm GW}}$}
\put(170,-60)
{
  \framebox(70,30)
  {
    \texttt{OmegaGW}
  }
}
\put(280,-30){\vector(0,1){20}}
\put(280,-30){\line(0,1){30}}
\put(282,-15){$\gamma$}
\put(280,-20){\line(-1,0){130}}
\put(150,-20){\line(0,1){75}}
\put(150,-20){\vector(0,1){55}}
\put(152,30){$\gamma$}
\put(250,-60)
{
  \framebox(70,30)
  {
    \texttt{Overlap}
  }
}
\put(125,-80){\line(0,1){135}}
\put(125,-80){\vector(0,1){100}}
%\put(107,40){$({P^{\scriptstyle{\rm HC}}_{1,2}})^{-1}$}
\put(127,20){$\frac{1}{P^{\scriptstyle{\rm HC}}_{1,2}}$}
%\put(117,40){$\frac{\tilde{R}_{1,2}}{P^{\scriptstyle{\rm C}}_{1,2}}$}
\put(330,-80){\line(0,1){80}}
\put(330,-80){\vector(0,1){45}}
% \put(142,40){${P^{\scriptstyle{\rm C}}_{1,2}}^{-1}$}
\put(332,-35){$\frac{1}{P^{\scriptstyle{\rm C}}_{1,2}}$}
%\put(142,40){$\frac{|\tilde{R}_{1,2}|^2}{P^{\scriptstyle{\rm C}}_{1,2}}$}
\put(115,-110)
{
  \framebox(225,30)
  {
    \texttt{InverseNoise}
  }
}
\put(73,-108){$\tilde{R}_{1,2}$}
\put(93,-105){\vector(1,0){15}}
\put(93,-105){\line(1,0){25}}
\put(73,-90){$P_{1,2}$}
\put(93,-87){\vector(1,0){15}}
\put(93,-87){\line(1,0){25}}
\end{picture}
\end{center}
\caption{\label{stochastic:f:CrossCorrFlowchart} Relationship among
  the modules dependent on \texttt{StochasticCrossCorrelation.h},
  which are used to calculate the cross-correlation statistic $Y$
  and its theoretical variance per unit time $\sigma^2/T$ from
  (uncalibrated) stretches of data $h_1(t)$,
  $h_2(t)$, from two detectors, using metadata
  on the power spectral densities $P_1(f)$,
  $P_2(f)$ and transfer functions $\tilde{R}_1(f)$,
  $\tilde{R}_2(f)$ for each detector.
  \texttt{CrossCorr} represents the module
  \texttt{StochasticCrossCorrelation.c}
  (Sec.~\ref{stochastic:ss:StochasticCrossCorrelation.c})
  containing the functions
  \texttt{LALStochasticCrossCorrelationStatistic()},
  \texttt{LALStochasticHeterodynedCrossCorrelationStatistic()},
  and \texttt{LALStochasticCrossCorrelationSpectrum()},
  \texttt{ZeroPadAndFFT} represents the module
  \texttt{ZeroPadAndFFT.c} (Sec.~\ref{stochastic:ss:ZeroPadAndFFT.c})
  containing the functions
  \texttt{LALSZeroPadAndFFT()} and \texttt{LALCZeroPadAndFFT()};
  \texttt{OptimalFilter} represents the module
  \texttt{StochasticOptimalFilter.c}
  (Sec.~\ref{stochastic:ss:StochasticOptimalFilter.c})
  containing the function
  \texttt{LALStochasticOptimalFilter()};
  \texttt{Normalization} represents the module
  \texttt{StochasticOptimalFilterNormalization.c}
  (Sec.~\ref{stochastic:ss:StochasticOptimalFilterNormalization.c})
  containing the function
  \texttt{LALStochasticOptimalFilterNormalization()};
  \texttt{InverseNoise} represents the module
  \texttt{StochasticInverseNoise.c}
  (Sec.~\ref{stochastic:ss:StochasticInverseNoise.c})
  containing the function
  \texttt{LALStochasticInverseNoise()};
  \texttt{OmegaGW} represents the module
  \texttt{StochasticOmegaGW.c}
  (Sec.~\ref{stochastic:ss:StochasticOmegaGW.c})
  containing the function
  \texttt{LALStochasticOmegaGW()};
  \texttt{Overlap} represents the module
  \texttt{OverlapReductionFunction.c}
  (Sec.~\ref{stochastic:ss:OverlapReductionFunction.c})
  containing the function
  \texttt{OverlapReductionFunction()}.
 }
\end{figure}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/StochasticCrossCorrelation.h>
\end{verbatim}

\noindent

\subsection*{Error conditions}
\input{StochasticCrossCorrelationHE}

\subsection*{Structures}

*********************************************************** </lalLaTeX> */

#ifndef _STOCHASTICCROSSCORRELATION_H
#define _STOCHASTICCROSSCORRELATION_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( STOCHASTICCROSSCORRELATIONH,
        "$Id$" );

/****************** <lalErrTable file="StochasticCrossCorrelationHE"> */

#define STOCHASTICCROSSCORRELATIONH_ENULLPTR        1
#define STOCHASTICCROSSCORRELATIONH_ESAMEPTR        2
#define STOCHASTICCROSSCORRELATIONH_EZEROLEN        3
#define STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF   4
#define STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT   5
#define STOCHASTICCROSSCORRELATIONH_ENEGFMIN        6
#define STOCHASTICCROSSCORRELATIONH_EMMTIME         7
#define STOCHASTICCROSSCORRELATIONH_EMMHETERO       8
#define STOCHASTICCROSSCORRELATIONH_EMMFMIN         9
#define STOCHASTICCROSSCORRELATIONH_EMMDELTAF      10
#define STOCHASTICCROSSCORRELATIONH_EMMLEN         11
#define STOCHASTICCROSSCORRELATIONH_EOORFREF       12
#define STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA   13
#define STOCHASTICCROSSCORRELATIONH_ENONSYMDIJ     14
#define STOCHASTICCROSSCORRELATIONH_ENONZEROHETERO 15
#define STOCHASTICCROSSCORRELATIONH_EWRONGUNITS    16
#define STOCHASTICCROSSCORRELATIONH_ENONPOSWIN     17
#define STOCHASTICCROSSCORRELATIONH_EMEMORY        18
#define STOCHASTICCROSSCORRELATIONH_ENOTYETHETERO 255

#define STOCHASTICCROSSCORRELATIONH_MSGENULLPTR    "Null pointer"
#define STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR    "Input and Output pointers the same"
#define STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN    "Zero length for data member of series"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF "Negative or zero frequency spacing"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT "Negative or zero time spacing"
#define STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN "Negative start frequency"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMTIME     "Mismatch in epochs"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMHETERO   "Mismatch in heterodyning frequencies"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN     "Mismatch in start frequencies"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF   "Mismatch in frequency spacings"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMLEN      "Mismatch in sequence lengths"
#define STOCHASTICCROSSCORRELATIONH_MSGEOORFREF    "Out of range reference frequency"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA "Negative stochastic background strength"
#define STOCHASTICCROSSCORRELATIONH_MSGENONSYMDIJ   "Non-symmetric response tensor"
#define STOCHASTICCROSSCORRELATIONH_MSGENONZEROHETERO "Non-zero heterodyning frequency specified for real time series"
#define STOCHASTICCROSSCORRELATIONH_MSGEWRONGUNITS "Inconsistent input units"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSWIN  "Zero or negative total for window functions"
#define STOCHASTICCROSSCORRELATIONH_MSGEMEMORY     "Memory error"
#define STOCHASTICCROSSCORRELATIONH_MSGENOTYETHETERO   "Non-zero heterodyning frequency not yet implemented"

/************************************ </lalErrTable> */

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             StochasticCrossCorrelation.c                  *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures and prototypes associated with
  \texttt{StochasticCrossCorrelation.c}
  (Sec.~\ref{stochastic:ss:StochasticCrossCorrelation.c})}

\subsubsection*{Prototypes}

\idx{LALStochasticCrossCorrelationStatistic()}
\idx{LALStochasticHeterodynedCrossCorrelationStatistic()}
\idx{LALStochasticCrossCorrelationSpectrum()}
\input{StochasticCrossCorrelationHPCC}

\subsubsection*{\texttt{struct REAL4WithUnits}}
\idx[Type]{REAL4WithUnits}

\noindent
Represents a dimensionful number as a 4-byte float with an associated
units structure, which is the output of
\texttt{LALStochasticCrossCorrelationStatistic()}.  The fields are:

\begin{description}
\item[\texttt{REAL4 value}]
The numerical value.

\item[\texttt{LALUnit units}]
The units.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagREAL4WithUnits {
  REAL4     value;
  LALUnit   units;
} REAL4WithUnits;

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct COMPLEX8WithUnits}}
\idx[Type]{COMPLEX8WithUnits}

\noindent
Represents a dimensionful number as a single-precision (8-byte) complex
number with an associated
units structure, which is the output of
\texttt{LALStochasticHeterodynedCrossCorrelationStatistic()}.  The fields are:

\begin{description}
\item[\texttt{COMPLEX8 value}]
The numerical value.

\item[\texttt{LALUnit units}]
The units.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagCOMPLEX8WithUnits {
  COMPLEX8  value;
  LALUnit   units;
} COMPLEX8WithUnits;

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct StochasticCrossCorrelationInput}}
\idx[Type]{StochasticCrossCorrelationInput}

\noindent Contains the input data needed by
\texttt{LALStochasticCrossCorrelationStatistic()}
to calculate the value of the standard optimally-filtered
cross-correlation statistic.  The fields are:

\begin{description}
\item[\texttt{COMPLEX8FrequencySeries  *hBarTildeOne}]
Fourier transform of the first zero-padded data stream.

\item[\texttt{COMPLEX8FrequencySeries  *hBarTildeTwo}]
Fourier transform of the second zero-padded data stream.

\item[\texttt{COMPLEX8FrequencySeries  *optimalFilter}]
Optimal filter function in the frequency domain.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticCrossCorrelationInput {
  COMPLEX8FrequencySeries  *hBarTildeOne;
  COMPLEX8FrequencySeries  *hBarTildeTwo;
  COMPLEX8FrequencySeries  *optimalFilter;
} StochasticCrossCorrelationInput;

typedef struct tagStochasticCrossCorrelationStrainInput {
  COMPLEX8FrequencySeries  *hBarTildeOne;
  COMPLEX8FrequencySeries  *hBarTildeTwo;
  REAL4FrequencySeries  *optimalFilter;
} StochasticCrossCorrelationStrainInput;


typedef struct tagStochasticCrossCorrelationCalInput {
  COMPLEX8FrequencySeries  *hBarTildeOne;
  COMPLEX8FrequencySeries  *hBarTildeTwo;
  REAL4FrequencySeries     *optimalFilter;
  COMPLEX8FrequencySeries  *responseFunctionOne;
  COMPLEX8FrequencySeries  *responseFunctionTwo;
} StochasticCrossCorrelationCalInput;


/********** <lalVerbatim file="StochasticCrossCorrelationHPCC"> *********/

void
LALStochasticCrossCorrelationStatistic(
            LALStatus                              *status,
            REAL4WithUnits                         *output,
            const StochasticCrossCorrelationInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticHeterodynedCrossCorrelationStatistic(
            LALStatus                              *status,
            COMPLEX8WithUnits                      *output,
            const StochasticCrossCorrelationInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationSpectrum(
            LALStatus                              *status,
            COMPLEX8FrequencySeries                *output,
            const StochasticCrossCorrelationInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationStatisticStrain(
            LALStatus                              *status,
            REAL4WithUnits                         *output,
            const StochasticCrossCorrelationStrainInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticHeterodynedCrossCorrelationStatisticStrain(
            LALStatus                              *status,
            COMPLEX8WithUnits                      *output,
            const StochasticCrossCorrelationStrainInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationSpectrumStrain(
            LALStatus                              *status,
            COMPLEX8FrequencySeries                *output,
            const StochasticCrossCorrelationStrainInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationStatisticCal(
            LALStatus                                 *status,
            REAL4WithUnits                            *output,
            const StochasticCrossCorrelationCalInput  *input,
            BOOLEAN                                    epochsMatch);

void
LALStochasticHeterodynedCrossCorrelationStatisticCal(
            LALStatus                                 *status,
            COMPLEX8WithUnits                         *output,
            const StochasticCrossCorrelationCalInput     *input,
            BOOLEAN                                    epochsMatch);

void
LALStochasticCrossCorrelationSpectrumCal(
            LALStatus                                 *status,
            COMPLEX8FrequencySeries                   *output,
            const StochasticCrossCorrelationCalInput  *input,
            BOOLEAN                                   epochsMatch);

/********** </lalVerbatim> *********/

  /*************************************************************
   *                                                           *
   * Structures and prototypes associated with ZeroPadAndFFT.c *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct SZeroPadAndFFTParameters}}
\idx[Type]{SZeroPadAndFFTParameters}

\noindent
Contains the parameters of \texttt{LALSZeroPadAndFFT()}.
The fields are:

\begin{description}
\item[\texttt{RealFFTPlan *fftPlan}]
The FFT plan to be used by FFTW
\item[\texttt{REAL4Vector *window}]
The window which is to be applied to the data
\item[\texttt{UINT4 length}]
The length of the data after zero-padding
\end{description}

\subsubsection*{\texttt{struct CZeroPadAndFFTParameters}}
\idx[Type]{CZeroPadAndFFTParameters}

\noindent
Contains the parameters of \texttt{LALCZeroPadAndFFT()}.
The fields are:

\begin{description}
\item[\texttt{ComplexFFTPlan *fftPlan}]
The FFT plan to be used by FFTW
\item[\texttt{REAL4Vector *window}]
The window which is to be applied to the data
\item[\texttt{UINT4 length}]
The length of the data after zero-padding
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagSZeroPadAndFFTParameters {
  RealFFTPlan           *fftPlan;
  REAL4Vector           *window;
  UINT4                  length;
} SZeroPadAndFFTParameters;

typedef struct tagCZeroPadAndFFTParameters {
  ComplexFFTPlan        *fftPlan;
  REAL4Vector           *window;
  UINT4                  length;
} CZeroPadAndFFTParameters;

/********************************************************** <lalLaTeX>

\subsubsection*{Prototypes associated with
  \texttt{ZeroPadAndFFT.c}
  (Sec.~\ref{stochastic:ss:ZeroPadAndFFT.c})}

\idx{LALSZeroPadAndFFT()}
\idx{LALCZeroPadAndFFT()}
\input{StochasticCrossCorrelationHPZP}

********** </lalLaTeX> *********/

/********** <lalVerbatim file="StochasticCrossCorrelationHPZP"> *********/

void
LALSZeroPadAndFFT(LALStatus                *status,
                  COMPLEX8FrequencySeries  *output,
                  const REAL4TimeSeries    *input,
                  SZeroPadAndFFTParameters *parameters);

void
LALCZeroPadAndFFT(LALStatus                *status,
                  COMPLEX8FrequencySeries  *output,
                  const COMPLEX8TimeSeries *input,
                  CZeroPadAndFFTParameters *parameters);

/********** </lalVerbatim> *********/

  /*************************************************************
   *                                                           *
   *   Structures and prototypes associated with               *
   *               StochasticOptimalFilter.c                   *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures and protoypes associated with
  \texttt{StochasticOptimalFilter.c}
  (Sec.~\ref{stochastic:ss:StochasticOptimalFilter.c})}

\subsubsection*{Prototypes}

\idx{LALStochasticOptimalFilter()}
\input{StochasticCrossCorrelationHPOF}

\subsubsection*{\texttt{struct StochasticOptimalFilterInput}}
\idx[Type]{StochasticOptimalFilterInput}

\noindent
Contains the inputs of \texttt{LALStochasticOptimalFilter()}.
The fields are:

\begin{description}
\item[\texttt{REAL4FrequencySeries *overlapReductionFunction}]
The overlap reduction function $\gamma(f)$ describing the pair of detector
sites.
\item[\texttt{REAL4FrequencySeries *omegaGW}] The spectrum
$\Omega_{\scriptstyle{\rm GW}}(f)$ of the stochastic gravitational-wave
background.
\item[\texttt{COMPLEX8FrequencySeries *halfCalibratedInverseNoisePSD1}]
 The reciprocal
$1/P_1^{\scriptstyle{\rm HC}}(f)
=1/(\tilde{R_1}(f)P_1(f))
=\tilde{R_1}(f)^* / P_1(f)$ of the
half-calibrated noise power spectral density for the first detector.
\item[\texttt{COMPLEX8FrequencySeries *halfCalibratedInverseNoisePSD2}]
 The reciprocal
$1/P_2^{\scriptstyle{\rm HC}}(f)
=1/(\tilde{R_2}(f)P_2(f))
=\tilde{R_2}(f)^* / P_2(f)$ of the
half-calibrated noise power spectral density for the second detector.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticOptimalFilterInput {
  REAL4FrequencySeries     *overlapReductionFunction;
  REAL4FrequencySeries     *omegaGW;
  COMPLEX8FrequencySeries  *halfCalibratedInverseNoisePSD1;
  COMPLEX8FrequencySeries  *halfCalibratedInverseNoisePSD2;
} StochasticOptimalFilterInput;

typedef struct tagStochasticOptimalFilterCalInput {
  REAL4FrequencySeries     *overlapReductionFunction;
  REAL4FrequencySeries     *omegaGW;
  REAL4FrequencySeries     *calibratedInverseNoisePSD1;
  REAL4FrequencySeries     *calibratedInverseNoisePSD2;
} StochasticOptimalFilterCalInput;

/********** <lalVerbatim file="StochasticCrossCorrelationHPOF"> *********/

void
LALStochasticOptimalFilter(
            LALStatus                                *status,
            COMPLEX8FrequencySeries                  *optimalFilter,
            const StochasticOptimalFilterInput       *input,
            const REAL4WithUnits                     *lambda);

void
LALStochasticOptimalFilterCal(
            LALStatus                                *status,
            REAL4FrequencySeries                     *optimalFilter,
            const StochasticOptimalFilterCalInput    *input,
            const REAL4WithUnits                     *lambda);

/********** </lalVerbatim> *********/

  /*************************************************************
   *                                                           *
   *        Structures and prototypes associated with          *
   *         StochasticOptimalFilterNormalization.c            *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures and protoypes associated with
  \texttt{StochasticOptimalFilterNormalization.c}
  (Sec.~\ref{stochastic:ss:StochasticOptimalFilterNormalization.c})}

\subsubsection*{Prototypes}

\idx{LALStochasticOptimalFilterNormalization()}
\input{StochasticCrossCorrelationHPON}

\subsubsection*{\texttt{struct StochasticOptimalFilterNormalizationOutput}}
\idx[Type]{StochasticOptimalFilterNormalizationOutput}

\noindent
Contains the outputs of \texttt{LALStochasticOptimalFilterNormalization()}.
The fields are:

\begin{description}
\item[\texttt{REAL4WithUnits *normalization}]
The normalization parameter $\lambda$.
\item[\texttt{REAL4WithUnits *variance}]
The variance per unit time $\sigma^2/T$ of the cross-correlation statistic.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticOptimalFilterNormalizationOutput {
  REAL4WithUnits           *normalization;
  REAL4WithUnits           *variance;
} StochasticOptimalFilterNormalizationOutput;

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct StochasticOptimalFilterNormalizationInput}}
\idx[Type]{StochasticOptimalFilterNormalizationInput}

\noindent
Contains the inputs of \texttt{LALStochasticOptimalFilterNormalization()}.
The fields are:

\begin{description}
\item[\texttt{REAL4FrequencySeries *overlapReductionFunction}]
The overlap reduction function $\gamma(f)$ describing the pair of detector
sites.
\item[\texttt{REAL4FrequencySeries *omegaGW}] The spectrum
$\Omega_{\scriptstyle{\rm GW}}(f)$ of the stochastic gravitational-wave
background.
\item[\texttt{REAL4FrequencySeries *inverseNoisePSD1}]
 The reciprocal
$1/P_1(f)=|\tilde{R_1}(f)|^2/P_1(f)$ of the
ununcalibrated noise power spectral density for the first detector.
\item[\texttt{REAL4FrequencySeries *inverseNoisePSD2}]
 The reciprocal
$1/P_2(f)=|\tilde{R_2}(f)|^2/P_2(f)$ of the
ununcalibrated noise power spectral density for the second detector.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticOptimalFilterNormalizationInput {
  REAL4FrequencySeries     *overlapReductionFunction;
  REAL4FrequencySeries     *omegaGW;
  REAL4FrequencySeries     *inverseNoisePSD1;
  REAL4FrequencySeries     *inverseNoisePSD2;
} StochasticOptimalFilterNormalizationInput;

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct StochasticOptimalFilterNormalizationParameters}}
\idx[Type]{StochasticOptimalFilterNormalizationParameters}

\noindent
Contains the parameters of \texttt{LALStochasticOptimalFilterNormalization()}.
The fields are:

\begin{description}
\item[\texttt{REAL8 fRef}]
The reference frequency used in defining the normalization.
\item[\texttt{BOOLEAN heterodyned}]
Indicates whether the filter is to be used on heterodyned data or not.
\item[\texttt{REAL4Vector window1}]
The windowing function with which the first data stream was windowed
\item[\texttt{REAL4Vector window2}]
The windowing function with which the second data stream was windowed
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticOptimalFilterNormalizationParameters {
  REAL8               fRef;
  BOOLEAN             heterodyned;
  REAL4Vector        *window1;
  REAL4Vector        *window2;
} StochasticOptimalFilterNormalizationParameters;

/********** <lalVerbatim file="StochasticCrossCorrelationHPON"> *********/

void
LALStochasticOptimalFilterNormalization(
            LALStatus                                            *status,
            StochasticOptimalFilterNormalizationOutput           *output,
            const StochasticOptimalFilterNormalizationInput      *input,
            const StochasticOptimalFilterNormalizationParameters *parameters);

/********** </lalVerbatim> *********/

  /*************************************************************
   *                                                           *
   * Structures and prototypes associated with StochasticInverseNoise.c  *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures and prototypes associated with
  \texttt{StochasticInverseNoise.c}
  (Sec.~\ref{stochastic:ss:StochasticInverseNoise.c})}

\subsubsection*{Prototypes}

\idx{LALStochasticInverseNoise()}
\input{StochasticCrossCorrelationHPIN}

\subsubsection*{\texttt{struct StochasticInverseNoiseOutput}}
\idx[Type]{StochasticInverseNoiseOutput}

\noindent
Contains the outputs of \texttt{LALStochasticInverseNoise()}.
The fields are:

\begin{description}
\item[\texttt{REAL4FrequencySeries *calibratedInverseNoisePSD}]
The reciprocal
$1/P^{\scriptstyle{\rm C}}(f)=|\tilde{R}(f)|^2/P(f)$ of the
ununcalibrated noise power spectral density.

\item[\texttt{COMPLEX8FrequencySeries *halfCalibratedInverseNoisePSD}]
The reciprocal \\
$1/P^{\scriptstyle{\rm HC}}(f)=\tilde{R}(f)^* / P(f)$
of the half-calibrated noise power spectral density.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticInverseNoiseOutput {
  REAL4FrequencySeries     *calibratedInverseNoisePSD;
  COMPLEX8FrequencySeries  *halfCalibratedInverseNoisePSD;
} StochasticInverseNoiseOutput;

typedef struct tagStochasticInverseNoiseCalOutput {
  REAL4FrequencySeries     *calibratedInverseNoisePSD;
} StochasticInverseNoiseCalOutput;

/********************************************************** <lalLaTeX>
\subsubsection*{\texttt{struct StochasticInverseNoiseInput}}
\idx[Type]{StochasticInverseNoiseInput}

\noindent
Contains the inputs to \texttt{LALStochasticInverseNoise()}.
The fields are:

\begin{description}
\item[\texttt{REAL4FrequencySeries *unCalibratedNoisePSD}]
The power spectral density $P(f)$ of the noise
contribution to the detector output.

\item[\texttt{COMPLEX8FrequencySeries *responseFunction}]
The frequency-domain reponse function $\tilde{R}(f)$.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticInverseNoiseInput {
  REAL4FrequencySeries     *unCalibratedNoisePSD ;
  COMPLEX8FrequencySeries  *responseFunction;
} StochasticInverseNoiseInput;

/********** <lalVerbatim file="StochasticCrossCorrelationHPIN"> *********/

void
LALStochasticInverseNoise(
            LALStatus                             *status,
            StochasticInverseNoiseOutput          *output,
            const StochasticInverseNoiseInput     *input);

void
LALStochasticInverseNoiseCal(
            LALStatus                             *status,
            StochasticInverseNoiseCalOutput       *output,
            const StochasticInverseNoiseInput     *input);

/********** </lalVerbatim> *********/

  /*************************************************************
   *                                                           *
   *    Structures and prototypes associated with StochasticOmegaGW.c    *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures and protypes associated with
  \texttt{StochasticOmegaGW.c}
  (Sec.~\ref{stochastic:ss:StochasticOmegaGW.c})}

\subsubsection*{Prototypes}

\idx{LALStochasticOmegaGW()}
\input{StochasticCrossCorrelationHPOG}

\subsubsection*{\texttt{struct StochasticOmegaGWParameters}}
\idx[Type]{StochasticOmegaGWParameters}

\noindent
Contains the parameters used by \texttt{LALStochasticOmegaGW()} to define a
power law: $\Omega_{\scriptstyle{\rm GW}}(f)
= \Omega_{\scriptstyle{\rm R}} (f/f_{\scriptstyle{\rm R}})^\alpha$.
The fields are:

\begin{description}
\item[\texttt{REAL4 alpha}] The power-law exponent.

\item[\texttt{REAL8 fRef}] The reference frequency $f_{\scriptstyle{\rm
R}}$ used to define the normalization.

\item[\texttt{REAL4 omegaRef}] The amplitude
$\Omega_{\scriptstyle{\rm R}}
=\Omega_{\scriptstyle{\rm GW}}(f_{\scriptstyle{\rm R}})$
at reference frequency.

\item[\texttt{UINT4 length}]
The number of points in the output frequency series.

\item[\texttt{REAL8 f0}]
The start frequency of the output frequency series.

\item[\texttt{REAL8 deltaF}]
The frequency spacing of the output frequency series.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagStochasticOmegaGWParameters {
  REAL4     alpha;    /* exponent in power law: omegaGW(f) = f^alpha */
  UINT4     length;   /* length of vector containing omegaGW(f) values */
  REAL8     f0;       /* start frequency */
  REAL8     deltaF;   /* frequency spacing */
  REAL8     fRef;    /* reference normalization frequency */
  REAL4     omegaRef; /* refenence omega coefficent for normalization */
}
StochasticOmegaGWParameters;

/********** <lalVerbatim file="StochasticCrossCorrelationHPOG"> *********/

void
LALStochasticOmegaGW (
            LALStatus                          *status,
            REAL4FrequencySeries               *output,
            const StochasticOmegaGWParameters  *parameters);

/********** </lalVerbatim> *********/

  /*************************************************************
   *                                                           *
   *      Structures and prototypes associated with            *
   *            OverlapReductionFunction.c                     *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures and prototypes associated with
  \texttt{OverlapReductionFunction.c}
  (Sec.~\ref{stochastic:ss:OverlapReductionFunction.c})}

\subsubsection*{Prototypes}

\idx{LALOverlapReductionFunction()}
\input{StochasticCrossCorrelationHPOR}

\subsubsection*{\texttt{struct OverlapReductionFunctionParameters}}
\idx[Type]{OverlapReductionFunctionParameters}

\noindent Contains the parameters used by
\texttt{LALOverlapReductionFunction()} to determine the format of its
output for the overlap reduction function.  The fields are:

\begin{description}
\item[\texttt{UINT4 length}]
The number of points in the output frequency series.

\item[\texttt{REAL8 f0}]
The start frequency of the output frequency series.

\item[\texttt{REAL8 deltaF}]
The frequency spacing of the output frequency series.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagOverlapReductionFunctionParameters {
  UINT4     length;   /* length of vector containing overlap red function */
  REAL8     f0;       /* start frequency */
  REAL8     deltaF;   /* frequency spacing for overlap reduction function */
}
OverlapReductionFunctionParameters;

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct LALDetectorPair}}
\idx[Type]{LALDetectorPair}

\noindent Holds structures defining the location and orientation of a
pair of gravitational wave detectors.  This is the input to
\texttt{LALOverlapReductionFunction()}.  The fields are:

\begin{description}
\item[\texttt{LALDetector detectorOne}]
The first interferometer.

\item[\texttt{LALDetector detectorTwo}]
The second interferometer.
\end{description}

*********************************************************** </lalLaTeX> */

typedef struct tagLALDetectorPair {
  LALDetector    detectorOne;
  LALDetector    detectorTwo;
}
LALDetectorPair;

/********** <lalVerbatim file="StochasticCrossCorrelationHPOR"> *********/

void
LALOverlapReductionFunction(
                   LALStatus                                  *status,
                   REAL4FrequencySeries                       *output,
                   const LALDetectorPair                      *detectors,
                   const OverlapReductionFunctionParameters   *parameters);

/********** </lalVerbatim> *********/

#ifdef  __cplusplus
}
#endif /* C++ protection */

#endif /* _STOCHASTICCROSSCORRELATION_H */

/********************************************************** <lalLaTeX>

\vfill{\footnotesize\input{StochasticCrossCorrelationHV}}

\newpage\input{StochasticCrossCorrelationC}
\newpage\input{StochasticCrossCorrelationStatisticTestC}
\newpage\input{StochasticHeterodynedCrossCorrelationStatisticTestC}
\newpage\input{StochasticCrossCorrelationSpectrumTestC}
\newpage\input{ZeroPadAndFFTC}
\newpage\input{SZeroPadAndFFTTestC}
\newpage\input{CZeroPadAndFFTTestC}
\newpage\input{StochasticOptimalFilterC}
\newpage\input{StochasticOptimalFilterTestC}
\newpage\input{StochasticOptimalFilterNormalizationC}
\newpage\input{StochasticOptimalFilterNormalizationTestC}
\newpage\input{StochasticInverseNoiseC}
\newpage\input{StochasticInverseNoiseTestC}
\newpage\input{StochasticOmegaGWC}
\newpage\input{StochasticOmegaGWTestC}
\newpage\input{OverlapReductionFunctionC}
\newpage\input{OverlapReductionFunctionTestC}

*********************************************************** </lalLaTeX> */
