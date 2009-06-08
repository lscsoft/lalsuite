/*
*  Copyright (C) 2007 Chad Hanna, David Churches, Duncan Brown, Jolien Creighton, Benjamin Owen, B.S. Sathyaprakash, Anand Sengupta, Craig Robinson , Thomas Cokelaer, Evan Ochsner
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

/* <lalVerbatim file="LALInspiralBankHV">

Author: Churches, D.K. and Sathyaprakash, B.S., Cokelaer, T.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALInspiralBank.h}}
\label{s:LALInspiralBank.h}

Header file for the template placement codes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALInspiralBank.h>
\end{verbatim}

\noindent This header file covers routines that are used in template placement.

</lalLaTeX> */

#ifndef _LALINSPIRALBANK_H
#define _LALINSPIRALBANK_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/LALNoiseModels.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID(LALINSPIRALBANKH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error codes}

</lalLaTeX>  */

/* <lalErrTable> */
#define LALINSPIRALBANKH_ENULL      1
#define LALINSPIRALBANKH_EMEM       2
#define LALINSPIRALBANKH_ECHOICE    3
#define LALINSPIRALBANKH_EDIV0      4
#define LALINSPIRALBANKH_ESIZE      8
#define LALINSPIRALBANKH_EFRANGE    16
#define LALINSPIRALBANKH_EORDER     32
#define LALINSPIRALBANKH_EGRIDSPACING 64
#define LALINSPIRALBANKH_EHEXAINIT 128
#define LALINSPIRALBANKH_EFCUT      5
#define LALINSPIRALBANKH_EFHIGH     6
#define LALINSPIRALBANKH_ENUMFCUT   7

#define LALINSPIRALBANKH_MSGENULL   "Null pointer"
#define LALINSPIRALBANKH_MSGEMEM    "Memory allocation failure"
#define LALINSPIRALBANKH_MSGECHOICE "Invalid choice for an input parameter"
#define LALINSPIRALBANKH_MSGEDIV0   "Division by zero"
#define LALINSPIRALBANKH_MSGESIZE   "Invalid input range"
#define LALINSPIRALBANKH_MSGEFRANGE "Limits outside range of frequency series"
#define LALINSPIRALBANKH_MSGEORDER  "Inappropriate PN order"
#define LALINSPIRALBANKH_MSGEGRIDSPACING "Inappropriate grid spacing parameter [SquareNotOriented or Hexagonal]"
#define LALINSPIRALBANKH_MSGEHEXAINIT "Empty bank. abnormal behaviour in HexaBank generation."
#define LALINSPIRALBANKH_MSGEFCUT   "Inappropriate cutoff frequency [SchwarzISCO, BKLISCO, LightRing, ERD, FRD or LRD]"
#define LALINSPIRALBANKH_MSGEFHIGH  "Final frequency is less than the low frequency cutoff."
#define LALINSPIRALBANKH_MSGENUMFCUT "Number of fcut must be greater or equal to 1"

/* </lalErrTable> */

/* <lalLaTeX>

\subsection*{Enums}
\begin{enumerate}

\item \texttt{CoordinateSpace:}
\input{LALCoordinateSpaceH}
Choose templates either in the $(\tau_0,\tau_2)$ or $(\tau_0,\tau_3)$
space.  This is one of the members of the InspiralCoarseBankIn structure.

This enum allows users to choose template bank either in the $(\tau_0, \tau_2)$
space of chirptimes (the choice made by \texttt{Tau0Tau2}) or in the
$(\tau_0, \tau_3)$ space of chirptimes (the choice made by \texttt{Tau0Tau3}).
This was implemented in releases before May 25, 2002. On May 25 we migrated to a
new, slightly faster, computation of the metric in which, at present, only the
choice \texttt{Tau0Tau3} can be made. Since October 2003 a new choice {\tt Psi0Psi3}
was added to handle BCV templates. In November 2007 two new choices were addded:
{\tt PTFIntrinctic} is a PTF metric in only the intrinsic parameters (a $4
\times 4$ matrix), and {\tt PTFFull} is the PTF metric in the full parameter
space (intrinsic and extrinsic parameters).

\item\texttt{InspiralBankMassRange:}
\input{LALInspiralBankMassRangeH}
An enum that appears in the \texttt{InspiralCoarseBankIn} structure
which fixes the way templates are chosen: The choice
\texttt{MinComponentMassMaxTotalMass} means the minimum of the
component masses will be given by \texttt{mMin} and maximum total
mass is given by \texttt{MMax} of the \texttt{InspiralBankCoarseIn} structure.
The choice \texttt{MinMaxComponentMass} means the minimum of the
components masses will be again fixed by \texttt{mMin} and the
maximum of the component masses is fixed by \texttt{mMax} of the
\texttt{InspiralCoarseIn} structure below.

\item\texttt{FreqCut:}
\input{LALFreqCutH}
An enum that lists all the formulas that can be used to specify an upper
frequency cutoff. From lowest to highest, the choices are:
\begin{itemize}
\item \texttt{FreqCut\_SchwarzISCO},
the innermost stable circular orbit (ISCO) for a test particle orbiting a
Schwarzschild black hole.
\item \texttt{FreqCut\_BKLISCO},
a mass ratio dependent ISCO derived from
estimates of the final spin of a merged black found in a paper by Buonanno,
Kidder, and Lehner (arXiv:0709.3839).
\item \texttt{FreqCut\_LightRing},
the unstable circular orbit
for photons orbiting a Schwarzschild black hole.
\item \texttt{FreqCut\_FRD},
the "Fundamental
RingDown" frequency which is calculated from the Berti, Cardoso and Will
(arXiv:gr-qc/0512160) value for the $\omega_{220}$ QNM frequency using mass
ratio dependent fits to the final BH mass and spin from Buonanno et al
(arXiv:0706.3732).
\item \texttt{FreqCut\_ERD},
an effective ringdown
frequency studied in Pan et al (arXiv:0704.1964) that was found to give good
fit between stationary-phase templates and  numerical relativity waveforms.
\item \texttt{FreqCut\_LRD},
the "Lorentzian RingDown" frequency = 1.2*FRD which captures part of the
Lorentzian tail from the decay of the QNMs.
\end{itemize}

\item\texttt{GridSpacing:}
\input{LALGridSpacingH}
This enum is set by the user to specify the type of placement requested. It can be
\texttt{Square, Hexagonal, SquareNotOriented, HexagonalNotOriented, S2BCV}. The two first
align the ellipse along the eigen-vectors whereas the two next do not. The last is a
square placement which was being used during S2 and is therefore obsolete and should
not be used (feel free to remove it). Historically, we used the \texttt{SquareNotOriented}
placement until S4. Then, in S5, we switched to the \texttt{Hexagonal} placement,
which should be used for future searches.

\item\texttt{Position:}
\input{LALPositionH}
This enum can take the following values \texttt{In, Out, Below, Edge, Above} and is used
 \emph{only} by the Hexagonal placement.  It simply specifies
the place of a point with respect to the parameter space. Edge, means that the ellipse
covers two boundaries(upper and lower).

\item\texttt{InsidePolygon:}
\input{LALInsidePolygonH}
This enum is set to true or false, it is just a boolean variable for the
 purpose of BCV placement but can be used in an other context..

\item\texttt{Generation}
\input{LALGenerationH}
This enum is either \texttt{fertile,sterile}, and is a boolean expression used \emph{only}
 by the Hexagonal placement.



\end{enumerate}

\subsection*{Structures}
\begin {enumerate}
\item \texttt{InspiralMetric}
Structure to store metric at various points the signal manifold.
\input{LALInspiralMetricH}
We store the diagonalized metric together with the angle theta
between the $\tau_0$-axis and the semi-major axis of the ambiguity ellipse.
The members of this structure are:
\begin{itemize}
\item \texttt{G00}: 00-component of the metric in $(\tau_0,\tau_{2(3)})$ coordinates.
\item \texttt{G11}: 11-component of the metric in $(\tau_0,\tau_{2(3)})$ coordinates.
\item \texttt{G01}: 01-component of the metric in $(\tau_0,\tau_{2(3)})$ coordinates.
\item \texttt{g00}: 00-component of the diagonalised metric.
\item \texttt{g11}: 11-component of the diagonalised metric.
\item \texttt{theta}:  Angle from tau0 to semi-major axis of the ellipse.
\item \texttt{Gamma[6]}: 3d metric co-efficients in $(t_C, \tau_0,\tau_{2(3)})$ coordinates.
\item \texttt{space}:  The enum describing the coordinate space in which
the metric is computed.
\end{itemize}

\item \texttt{InspiralCoarseBankIn:}
Input for choosing a template bank. This is the structure that must
	be filled by a routine calling the code \texttt{InspiralCreateCoarseBank} or \texttt{InspiralCreateBCVBank}.
Unless BCV template bank is needed (that is, \texttt{InspiralCreateBCVBank})  then one can ignore the
parameters \texttt{psi0Min, psi0Max, psi3Min, psi3Max, alpha, numFcutTemplates.}

\input{LALInspiralCoarseBankH}

\begin{itemize}
\item \texttt{massRange}:   enum that determines whether templates should be
	chosen using fixed ranges for component masses or
	to use minimum component mass and maximum totalmass.
\item \texttt{space}: enum that decides whether to use $(\tau_0,\tau_2)$
        or $(\tau_0,\tau_3)$ in constructing the template bank
\item \texttt{alpha}: 	the BCV amplitude correction parameter
\item \texttt{psi0Min}: minimum value of the parameter $\psi_0$
\item \texttt{psi0Max}: maximum value of the parameter $\psi_0$
\item \texttt{psi3Min}: minimum value of the parameter $\psi_3$
\item \texttt{psi3Max}: maximum value of the parameter $\psi_3$
\item \texttt{mMin}: 	minimum mass of components to search for
\item \texttt{mMax}: 	maximum mass of components to search for
\item \texttt{MMax}:   	alternatively, maximum total mass of binary to search for
\item \texttt{mmCoarse}:Coarse grid minimal match
\item \texttt{mmFine}:  Fine grid minimal match
\item \texttt{fLower}:  Lower frequency cutoff
\item \texttt{fUpper}:  Upper frequency cutoff
\item \texttt{tSampling}:  Sampling rate
\item \texttt{etamin}: 	minimum value of eta in our search
\item \texttt{shf}: 	Frequency series containing the PSD
\item \texttt{iflso}: 	(currently not implemented) flso will be used as an
\item \texttt{numFcutTemplates}: number of templates in the {\tt fcut} direction

The next two members are used in setting up the InspiralTemplate
parameter structure but not in creating the template bank.

\item \texttt{order}: Post-Newtonian order of the waveform
\item \texttt{approximant}: Approximant of the waveform
\item \texttt{numFreqCut}: Number of different upper frequency cutoffs (spaced evenly between minFreqCut and maxFreqCut) to use when creating a template bank.
\item \texttt{maxFreqCut}: largest upper frequency cutoff to use
\item \texttt{minFreqCut}: smallest upper frequency cutoff to use
\end{itemize}

\item \texttt{InspiralFineBankIn}
Structure needed by the function \texttt{LALInspiralCreateFineBank}.
	which computes a finer mesh around a given lattice point
	using the value of the fine-mesh minimal match, coarse-mesh
	minimal match and the metric at the current lattice point.
\input{LALInspiralFineBankInH}
\begin{itemize}
\item {templateList:} A list containing all the fine-mesh templates
\item {coarseIn:} input structure that contains useful necessary parameters
to construct a fine-mesh.
\end{itemize}

\item \texttt{InspiralTemplateList}
A grid of inspiral templates (i.e., a template list).

\input{LALInspiralTemplateListH}
Structure returned by the coarse and fine bank generation routines.
Currently we generate an array of type \texttt{InspiralTemplateList}
which contains the coordinate markers (the parameter structure
\texttt{InspiralTemplate} defined in the \texttt{inspiral} package)
and the metric at each of those points. There is a desire to make this
a truly linked list at some time in the future. The member of this
structure are:
\begin{itemize}
\item \texttt{ID}: An unique integer ID of the template
\item \texttt{params}: Value of the parameters at the lattice point
\item \texttt{metric}:  metric at the lattice point
\item \texttt{*next}:  pointer to next lattice point; but this is currently
not filled by the bank code.
\end{itemize}

\item \texttt{InspiralBankParams:}
This is a structure needed in the inner workings
of the \texttt{LALInspiralCreateCoarseBank} code.
\input{LALInspiralBankParamsH}
\begin{itemize}
\item \texttt{nparams}: Number of parameters (currently fixed at 2, so this
		is as of now unused)
\item \texttt{x0}: the first coordinate, chosen to be always $\tau_0$
\item \texttt{x1}: the second coordinate, chosen to be either $\tau_2$ or $\tau_3$
\item \texttt{dx0}: increment in the x0-direction
\item \texttt{dx1}: increment in the x1-direction
\item \texttt{x0Min}: minimum value of the first coordinate as
defined by the search region
\item \texttt{x0Max}: maximum value of the first coordinate as
defined by the search region
\item \texttt{x1Min}: minimum value of the second coordinate as
defined by the search region
\item \texttt{x1Max}: maximum value of the second coordinate as
defined by the search region
\item \texttt{*metric}: pointer to the metric at the current location.
\end{itemize}


\item \texttt{InspiralMomentsIn}
Inputs to the function that computes the moments of the PSD.
	The moment is defined as:
	$$I(p) \equiv \int_{x_{\rm min}}^{x_{\rm max}}
\frac{x^{-p}}{S_h(x)} dx,$$
	where $x=f/f_0$ is a scaled frequency, $f_0$
	being a fiducial frequency, taken in these routines
	as the user supplied lower cutoff of the detector
	response.
\input{LALInspiralMomentsInH}
\begin{itemize}
\item \texttt{xmin}: lower limit of the integral $x_{\rm min}$
\item \texttt{xmax}: upper limit of the integral $x_{\rm max}$
\item \texttt{ndx}: index $p$ (without the negative sign) in the moment integral as above
\item \texttt{norm}: norm to be used in computing the moment, the returned value is
the above integral divided by the norm.
\item \texttt{*shf}: the frequency series containing the noise psd.
\end{itemize}


\item \texttt{InspiralMomentsEtc}
Parameter structure that holds the moments of the PSD and other useful
	constants required in the computation of the metric.
\input{LALInspiralMomentsEtcH}
\begin{itemize}
\item {a01, a21, \ldots:} Coefficients in the expansion of the phase
	of the Fourier transform of an inspiral waveform computed
	in the stationary phase approximation. See documentation under
	the function \texttt{LALInspiralComputeMetric} later in this
	Section for a description of these coefficients.
\item\texttt{j[18]:} The required moments are all computed once and
stored in this array. The required moments are from J(1) to J(17)
(except J(2), J(3) and J(16) that are not required at 2PN order,
 however, they are computed since future extensions, planned in the
 near future, will require them). However, in C we need an array size
18 to use an array that has an index 18. To ease the notation we have
therefore defined an over sized (by one element) array.
\end{itemize}

\item {\texttt{RectangleIn} and \texttt{RectangleOut}:}
Input and output structures to function LALRectangleVertices.
\input{LALRectangleInH}
\input{LALRectangleOutH}



\item\texttt{HexaGridParam}
This is a structure needed in the inner workings of the \texttt{LALInspiralHexagonalBank} code.
\input{LALHexaGridParamH}
It contains some part of CoarseBankIn and some other standard parameters.  It provides the
parameter space boundaries with the minimum and maximum values of mass parameters, the
minimal match, the space, massRange and gridSpacing parameter.

\item\texttt{CellEvolution}
This is a structure needed in the inner workings of the \texttt{LALInspiralHexagonalBank} code.
\input{LALCellEvolutionH}
This structure checks the status of the placement. \texttt{fertile} tells if the
placement is still evolving or not. \texttt{nTemplateMax} is the number of maximum templates allowed,
 which can be resized. And \texttt{nTemplate} is the number of template set. nTemplate can not
 be higher than nTemplateMax.

\item\texttt{CellList}
This is a structure needed in the inner workings of the \texttt{LALInspiralHexagonalBank} code.
\input{LALCellListH}
Similarly to the square placement, which uses InspiralList, we used a
linked list for the hexagonal placement. A different structure has been
implemented so as to simplify the complexity of the algorithm. It also set
an id to each cell which has been created. This id is unique to each
cell/template.

\item\texttt{InspiralCell}
This is a structure needed in the inner workings of the \texttt{LALInspiralHexagonalBank} code.
\input{LALInspiralCellH}
Each cell is defined by this structure, which contains the position of
each cell in the tau0/tau3 parameter space, the metric at that point, and various
information such as the status of the cell. Is it still fertile ? what is its position
with respect to the parameter space and so on. child is a 6-length array with a link
to the 6 templates (hexagonal) around the current template that we are dealing with.

\end{enumerate}
</lalLaTeX> */



/* <lalVerbatim file="LALComputeMomentsH"> */
typedef enum
{
  disable,
  enable
}
ComputeMoments;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{ComputeMoments}
</lalLaTeX>  */


/* <lalVerbatim file="LALCoordinateSpaceH"> */
typedef enum
{
  Tau0Tau2,
  Tau0Tau3,
  Psi0Psi3,
  PTFIntrinsic,
  PTFFull
}
CoordinateSpace;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{CoordinateSpace}
</lalLaTeX>  */


/* <lalVerbatim file="LALGridSpacingH"> */
typedef enum
{
  SquareNotOriented,
  Square,
  HexagonalNotOriented,
  Hexagonal,
  HybridHexagonal,
  S2BCV
}
GridSpacing;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{GridSpacing}
</lalLaTeX>  */

/* <lalVerbatim file="LALPositionH"> */
typedef enum
{
  In,
  Above,
  Below,
  Out,
  Edge
}
Position;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{Position}
</lalLaTeX>  */

/* <lalVerbatim file="LALInsidePolygonH"> */
typedef enum
{
  False,
  True
}
InsidePolygon;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InsidePolygon}
</lalLaTeX>  */


/* <lalVerbatim file="LALGenerationH"> */
typedef enum
{
  Sterile,
  Fertile
}
Generation;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{Generation}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralBankMassRangeH"> */
typedef enum
{
  MinComponentMassMaxTotalMass,
  MinMaxComponentMass,
  MinMaxComponentTotalMass
}
InspiralBankMassRange;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralBankMassRange}
</lalLaTeX>  */

/* <lalVerbatim file="LALFreqCutH"> */
typedef enum
{
  FreqCut_SchwarzISCO,
  FreqCut_BKLISCO,
  FreqCut_LightRing,
  FreqCut_ERD,
  FreqCut_FRD,
  FreqCut_LRD
}
FreqCut;
/* </lalVerbatim> */
/* <lalLaTeX>
\idx[Type]{FreqCut}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralMetricH"> */
typedef struct
tagInspiralMetric
{
  REAL8            G00;
  REAL8            G11;
  REAL8            G01;

  REAL8            g00;
  REAL8            g11;
  REAL8            theta;

  /* Gamma[] is a vector that stores the upper triangular part of the metric in
   * the space of parameters. For time domain searches, Gamma[0,...,5] stores
   * the following information :
   *    Gamma[0] -> (tc,tc) metric component
   *    Gamma[1] -> (tc,t0) metric component
   *    Gamma[2] -> (tc,t3) metric component
   *    Gamma[3] -> (t0,t0) metric component
   *    Gamma[4] -> (t0,t3) metric component
   *    Gamma[5] -> (t3,t3) metric component
   * For spinBCV searches, (in 4 dimensions) Gamma[0,...,9] would be required.
   */
  REAL4            Gamma[10];

  CoordinateSpace  space;
}
InspiralMetric;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralMetric}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralTemplateListH"> */
typedef struct
tagInspiralTemplateList
{
  INT4              ID;
  InspiralTemplate  params;
  InspiralMetric    metric;
  UINT4             nLayer;
  struct tagInspiralTemplateList *next;
}
InspiralTemplateList;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralTemplateList}
</lalLaTeX>  */


/* <lalVerbatim file="LALHexaGridParamH"> */
typedef struct
tagHexaGridParam
{
  REAL4 x0Min;
  REAL4 x1Min;
  REAL4 x0Max;
  REAL4 x1Max;
  REAL4 mm;
  REAL4 mMin;
  REAL4 mMax;
  REAL4 etaMin;
  REAL4 MMin;
  REAL4 MMax;
  REAL4 fLower;
  GridSpacing 		gridSpacing;
  InspiralBankMassRange massRange;
  CoordinateSpace       space;
}
HexaGridParam;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{HexaGridParam}
</lalLaTeX>  */

/* <lalVerbatim file="LALCellEvolutionH"> */
typedef struct
tagCellEvolution
{
  INT4 nTemplateMax;
  INT4 nTemplate;
  INT4 fertile;
}
CellEvolution;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{CellEvolution}
</lalLaTeX>  */

/* <lalVerbatim file="LALCellListH"> */
typedef struct
tagCellList
{
  INT4 id;
  struct tagCellList *next;
}
CellList;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{CellList}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralCellH"> */
typedef struct
tagInspiralCell
{
  INT4  ID;
  INT4  in;
  INT4  child[6];
  REAL4 t0;
  REAL4 t3;
  REAL4 dx0;
  REAL4 dx1;
  Generation status;
  Position position;
  Position RectPosition[5];
  InspiralMetric metric;
}
InspiralCell;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralCell}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankParamsH"> */
typedef struct
tagInspiralBankParams
{
  INT4           nparams;
  REAL8          minimalMatch;
  REAL8          x0;
  REAL8          x1;
  REAL8          dx0;
  REAL8          dx1;
  REAL8          x0Min;
  REAL8          x0Max;
  REAL8          x1Min;
  REAL8          x1Max;
  InspiralMetric *metric;
}
InspiralBankParams;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralBankParams}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralCoarseBankH"> */
typedef struct
tagInspiralCoarseBankIn
{
  InspiralBankMassRange         massRange;
  CoordinateSpace               space;

  REAL8                         mMin;
  REAL8                         mMax;
  REAL8                         MMax;
  REAL8                         MMin;
  REAL8                         alpha;
  REAL8                         psi0Min;
  REAL8                         psi0Max;
  REAL8                         psi3Min;
  REAL8                         psi3Max;
  REAL8                         mmCoarse;
  REAL8                         mmFine;
  REAL8                         fLower;
  REAL8                         fUpper;
  REAL8                         tSampling;
  REAL8                         etamin;
  REAL8				betaMin;
  REAL8				betaMax;
  REAL8                         chiMin;
  REAL8                         chiMax;
  REAL8                         kappaMin;
  REAL8                         kappaMax;
  INT4                          nPointsChi;
  INT4                          nPointsKappa;
  REAL8FrequencySeries          shf;
  /* Maximum size of the power spectral density array for use in
   * the computation of the metric in SBBH; typical values that
   * assures that the code runs quickly are 1024-8192.
   */
  UINT4				ShMaxSz;

  /* See for random number generation in RandomBank algorithm */
  UINT4				iseed;
  /* nTIni is an estimate for the number of templates that might
   * be required; this is used in the random bank generation
   * routine with a seed number of templates = nTIni*sqrt(nTIni)
   */
  UINT4				nTIni;
  /* iflso is an integer that tells whether to compute the moments
   * using an upper limit defined by flso; this is not used anywhere
   * at the moment
   */
  INT4                          iflso;
  /* spinBank=0:use Owen+Hanna bank*/
  /* spinBank=1:use extended bank by AEI/Cardiff/Osaka */
  /* spinBank=2:use random bank algorithm */
  INT4                          spinBank;
  /* Number of templates required in the fCut (upper cutoff)
   * dimension and the value of upper and lower cutoffs
   */

  UINT4                         numFcutTemplates;
  REAL4				HighGM;
  REAL4				LowGM;

  /* Type of gridspacing required:
  1=SquareNotOriented,
  2=Square,
  3=HexagonalNotOriented,
  4=Hexagonal,
  5=HybridHexagonal,
  6=S2BCV
  */
  GridSpacing                   gridSpacing;

  /* post-Newtonian order( phase), approximant, and amplitude PN order */
  LALPNOrder                    order;
  Approximant                   approximant;
  LALPNOrder                    ampOrder;

  /* parameters for different/multiple freq. cutoffs */
  INT4                          numFreqCut;
  FreqCut                       maxFreqCut;
  FreqCut                       minFreqCut;

  InsidePolygon                 insidePolygon;
  ComputeMoments                computeMoments;
  /* ComputeMoments tells whether to re-compute the moments
   * using an upper limit defined by flso; This is done after
   * the template bank is gnerated
   */
}
InspiralCoarseBankIn;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralCoarseBankIn}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralMomentsInH"> */
typedef struct
{
  REAL8                xmin;
  REAL8                xmax;
  REAL8                ndx;
  REAL8	             norm;
  REAL8FrequencySeries *shf;
}
InspiralMomentsIn;
/* </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralMomentsIn}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralFineBankInH"> */
typedef struct
tagInspiralFineBankIn
{
  InspiralTemplateList templateList;
  InspiralCoarseBankIn coarseIn;
}
InspiralFineBankIn;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralFineBankIn}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralMomentsEtcH"> */
typedef struct
tagInspiralMomentsEtc
{
  REAL8 a01, a21, a22, a31, a41, a42, a43;
  REAL8 j[18];
}
InspiralMomentsEtc;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralMomentsEtc}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralMomentsEtcBCVH"> */
typedef struct
tagInspiralMomentsEtcBCV
{
  REAL8 n0, n15;
  REAL8 j[9];
  REAL8 i[23];
  REAL8 alpha;
  REAL8 fcut;

  REAL8 M1[2][2];
  REAL8 M2[2][2];
  REAL8 M3[2][2];
}
InspiralMomentsEtcBCV;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\idx[Type]{InspiralMomentsEtcBCV}
</lalLaTeX>  */

/*  <lalVerbatim file="LALRectangleInH"> */
typedef struct
tagRectangleIn
{
  REAL8 x0, y0, dx, dy, theta;
}
RectangleIn;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\index{\texttt{RectangleIn}}
</lalLaTeX>  */

/*  <lalVerbatim file="LALRectangleOutH"> */
typedef struct
tagRectangleOut
{
  REAL8 x1;
  REAL8 y1;
  REAL8 x2;
  REAL8 y2;
  REAL8 x3;
  REAL8 y3;
  REAL8 x4;
  REAL8 y4;
  REAL8 x5;
  REAL8 y5;
}
RectangleOut;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\index{\texttt{RectangleOut}}
</lalLaTeX>  */

/*  <lalVerbatim file="LALHexagonOutH"> */
typedef struct
tagHexagonOut
{
  REAL8 x1;
  REAL8 y1;
  REAL8 x2;
  REAL8 y2;
  REAL8 x3;
  REAL8 y3;
  REAL8 x4;
  REAL8 y4;
  REAL8 x5;
  REAL8 y5;
  REAL8 x6;
  REAL8 y6;
  REAL8 x7;
  REAL8 y7;
}
HexagonOut;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\index{\texttt{HexagonOut}}
</lalLaTeX>  */


typedef struct{
REAL4 ct;
REAL4 b;
}
PRIN;


/*  <lalLaTeX>
\vfill{\footnotesize\input{LALInspiralBankHV}}
</lalLaTeX>  */





/* Function prototypes */
/* <lalLaTeX>
\newpage\input{LALInspiralCreateCoarseBankC}
</lalLaTeX>  */

void
LALInspiralCreateCoarseBank (
    LALStatus              *status,
    InspiralTemplateList   **list,
    INT4                   *nlist,
    InspiralCoarseBankIn   bankIn
    );

void
LALInspiralCreatePNCoarseBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

/* <lalLaTeX>
\newpage\input{LALInspiralBCVBankC}
</lalLaTeX>  */
void
LALInspiralCreateBCVBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralCreateFlatBankS3S4 (
    LALStatus            *status,
    REAL4VectorSequence  *list,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn coarseIn
    );

void
LALExcludeTemplate(
    LALStatus            *status,
    INT4 *valid,
    InspiralBankParams   *bankParams,
    REAL4 x,
    REAL4 y
    );

void
LALInspiralBCVBankFcutS3S4 (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                *NList,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralBCVFcutBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                *NList,
    InspiralCoarseBankIn coarseIn
    );

void
PSItoMasses (
    LALStatus            *status,
    InspiralTemplate 	                *params,
    UINT4 				*valid,
    REAL4 				highGM
);


void
LALEmpiricalPSItoMassesConversion(
    LALStatus            *status,
    InspiralTemplate    *params,
    UINT4               *valid,
    REAL4               lightring
);

void
LALPSItoMasses(
    LALStatus            *status,
    InspiralTemplate    *params,
    UINT4               *valid,
    REAL4               thisfreq
);

void
LALInspiralBCVRegularFcutBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                *NList,
    InspiralCoarseBankIn coarseIn);

void
LALNudgeTemplatesToConstantTotalMassLine(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 nlist,
    InspiralCoarseBankIn coarseIn
    );

/* <lalLaTeX>
   \newpage\input{LALInspiralBankUtilsC}
   </lalLaTeX>  */
REAL4
XLALInspiralTau3FromTau0AndEqualMassLine(
    REAL4               tau0,
    REAL4               fL
    );

REAL4
XLALInspiralTau3FromNonEqualMass(
  REAL4               	m1,
  REAL4 		m2,
  REAL4			fL
);

REAL4
XLALInspiralTau0FromMEta(
  REAL4              	M,
  REAL4 		eta,
  REAL4			fL
);


REAL8
XLALInspiralBissectionLine (
  REAL8 x,
  REAL8 fL,
  REAL8 mMin,
  REAL8 mMax);

REAL8
XLALInspiralMFromTau0AndNonEqualMass(
  REAL8 tau0,
  REAL8 extremMass,
  REAL8 fL);

/* <lalLaTeX>
   \newpage\input{InspiralSpinBankC}
   </lalLaTeX>  */

void
LALInspiralSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable    **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    );

void
LALInspiralBCVSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable    **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    );

#if 0
void
LALInspiralSpinBankBoundary(
    LALStatus            *status,
    NDTemplateBankInput  *input,
    NDTemplateBankOutput *output,
    INT2                 *flag
    );

void
LALInspiralSpinBankMetric(
    LALStatus           *status,
    NDTemplateBankInput *input,
    REAL4Array          *metric
    );
#endif

void
LALInspiralBankGeneration(
    LALStatus            *status,
    InspiralCoarseBankIn *in,
    SnglInspiralTable    **out,
    INT4                 *count
    );

void
LALInspiralCreateFlatBank (
    LALStatus            *status,
    REAL4VectorSequence  *list,
    InspiralBankParams   *bankParams
    );

/* <lalLaTeX>
   \newpage\input{LALInspiralCreateFineBankC}
   </lalLaTeX>  */

void
LALInspiralCreateFineBank (
    LALStatus              *status,
    InspiralTemplateList   **outlist,
    INT4                   *nlist,
    InspiralFineBankIn     fineIn
    );

/* <lalLaTeX>
   \newpage\input{LALInspiralComputeMetricC}
   </lalLaTeX>  */
void
LALInspiralComputeMetric (
    LALStatus           *status,
    InspiralMetric      *metric,
    InspiralTemplate    *params,
    InspiralMomentsEtc  *moments

    );
void
LALInspiralComputeMetricBCV
(
 LALStatus             *status,
 InspiralMetric        *metric,
 REAL8FrequencySeries  *psd,
 InspiralTemplate      *params
);

/* <lalLaTeX>
\newpage\input{LALInspiralLongestTemplateInBankC}
</lalLaTeX>  */
void
LALInspiralLongestTemplateInBank (
    LALStatus            *status,
    UINT4                *templateLength,
    InspiralCoarseBankIn *coarseIn
    );

/* <lalLaTeX>
\newpage\input{LALInspiralMomentsC}
</lalLaTeX>  */
void
LALGetInspiralMoments (
    LALStatus            *status,
    InspiralMomentsEtc   *moments,
    REAL8FrequencySeries *psd,
    InspiralTemplate     *params
    );

void
LALGetInspiralMomentsBCV (
    LALStatus               *status,
    InspiralMomentsEtcBCV   *moments,
    REAL8FrequencySeries    *psd,
    InspiralTemplate        *params
    );

void
LALInspiralMoments (
    LALStatus         *status,
    REAL8             *moment,
    InspiralMomentsIn pars
    );

/* <lalLaTeX>
\newpage\input{LALInspiralMomentsIntegrandC}
</lalLaTeX>  */
void
LALInspiralMomentsIntegrand
(
   LALStatus *status,
   REAL8  *integrand,
   REAL8  f,
   void   *pars
   );

/*  <lalLaTeX>
\newpage\input{LALInspiralSetSearchLimitsC}
</lalLaTeX>  */
void
LALInspiralSetSearchLimits (
    LALStatus            *status,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn coarseIn
    );

/* <lalLaTeX>
\newpage\input{LALInspiralNextTemplateC}
</lalLaTeX>  */
void
LALInspiralNextTemplate (
    LALStatus          *status,
    InspiralBankParams *bankPars,
    InspiralMetric      metric
    );

/*  <lalLaTeX>
\newpage\input{LALInspiralComputeParamsC}
</lalLaTeX>  */
void
LALInspiralComputeParams (
    LALStatus            *status,
    InspiralTemplate     *pars,
    InspiralBankParams   bankParams,
    InspiralCoarseBankIn coarseIn
    );

/* <lalLaTeX>
\newpage\input{LALInspiralValidParamsC}
</lalLaTeX>  */
void
LALInspiralValidParams (
    LALStatus            *status,
    INT4                 *valid,
    InspiralBankParams   bankParams,
    InspiralCoarseBankIn coarseIn
    );

/* <lalLaTeX>
\newpage\input{LALInspiralValidTemplateC}
</lalLaTeX>  */
void
LALInspiralValidTemplate(
   LALStatus            *status,
   INT4                 *valid,
   InspiralBankParams   bankParams,
   InspiralCoarseBankIn coarseIn
   );

/* <lalLaTeX>
\newpage\input{LALInspiralUpdateParamsC}
</lalLaTeX>  */
void
LALInspiralUpdateParams (
    LALStatus          *status,
    InspiralBankParams *bankParams,
    InspiralMetric     metric,
    REAL8              minimalMatch
    );

/* <lalLaTeX>
\newpage\input{LALMatrixTransformC}
</lalLaTeX>  */
void
LALMatrixTransform (
    LALStatus *status,
    INT4      Dim,
    REAL8     **trans,
    REAL8     **buff1,
    REAL8     **mm3
    );

/* <lalLaTeX>
\newpage\input{LALDeterminant3C}
</lalLaTeX>  */
void
LALDeterminant3 (
    LALStatus *status,
    REAL8  *determinant,
    REAL8  **matrix
    );

/* <lalLaTeX>
\newpage\input{LALInverse3C}
</lalLaTeX>  */
void
LALInverse3(
        LALStatus *status,
        REAL8     **inverse,
        REAL8     **matrix
);

/* <lalLaTeX>
\newpage\input{LALInspiralSetParamsC}
</lalLaTeX>  */
void
LALInspiralSetParams (
    LALStatus            *status,
    InspiralTemplate     *tempPars,
    InspiralCoarseBankIn coarseIn
    );

/* <lalLaTeX>
\newpage\input{LALRectangleVerticesC}
</lalLaTeX>  */

void
LALRectangleVertices(
   LALStatus *status,
   RectangleOut *out,
   RectangleIn *in
   );

void
LALHexagonVertices(
   LALStatus *status,
   HexagonOut *out,
   RectangleIn *in
   );


/* <lalLaTeX>
\newpage\input{LALInsidePolygonC}
</lalLaTeX>  */
void
LALInsidePolygon(
   LALStatus            *status,
   REAL4                *inputx,
   REAL4               *inputy,
   INT4                 n,
   REAL4                x,
   REAL4                y,
   INT4                 *valid
   );

/* <lalLaTeX>
   \newpage\input{LALInspiralHybridHexagonalBankC}
   </lalLaTeX>  */
void
LALInspiralCreatePNCoarseBankHexa(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralCreatePNCoarseBankHybridHexa(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void
LALInitHexagonalBank(
    LALStatus         *status,
    InspiralCell      **cell,
    INT4              id,
    InspiralMomentsEtc        *moments,
    InspiralTemplate          *paramsIn,
    HexaGridParam             *gridParam,
    CellEvolution             *cellEvolution,
    CellList                  **cellList
    );


void
LALPopulateCell(
    LALStatus               *status,
    InspiralMomentsEtc      *moments,
    InspiralCell            **cell,
    INT4                    l,
    InspiralTemplate        *paramsIn,
    HexaGridParam           *gridParam,
    CellEvolution           *cellEvolution,
    CellList **cellList
    );

void
LALFindPosition(
    LALStatus               *status,
    REAL4                   dx0,
    REAL4                   dx1,
    Position                *position,
    InspiralTemplate        *paramsIn,
    HexaGridParam           *gridParam
    );


void
LALSPAValidPosition(
    LALStatus           *status,
    InspiralCell        **cell,
    INT4                id1,
    InspiralMomentsEtc  *moments,
    CellEvolution *cellEvolution,
    CellList **list
    );

void
GetPositionRectangle(
    LALStatus *status,
    InspiralCell **cell,
    INT4 id,
    InspiralTemplate *params,
    HexaGridParam           *gridParam,
    CellEvolution *cellEvolution,
    CellList **cellList,
    INT4 *valid
    );



void
LALListAppend(
    CellList ** headRef,
    INT4 id
    );

UINT4
GetIdFromIndex(
    CellList *head,
    INT4 index
    );

void
LALListDelete(
    CellList ** headRef,
    INT4 id
    );

UINT4
LALListLength(
    CellList *head
    );

void
LALSPAF(
	LALStatus 	*status,
	REAL4 		*result,
	REAL4 		x,
	void 		*t3
);

/* <lalLaTeX>
\newpage\input{CoarseTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{CoarseTest2C}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{ChirpSpaceC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{GetOrientationEllipseC}
</lalLaTeX> */

/*<lalLaTeX>
\newpage\input{InspiralSpinBankTestC}
</lalLaTeX>*/

/* <lalLaTeX>
\newpage\input{SpaceCoveringC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{LALInspiralComputePTFMetricC}
</lalLaTeX> */
INT4 XLALInspiralComputePTFIntrinsicMetric (
    InspiralMetric             *metric,
	REAL8Vector				   *fullmetric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    );

INT4 XLALInspiralComputePTFFullMetric (
    InspiralMetric             *metric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    );

INT4 XLALInspiralComputePTFWaveform (
    REAL8Vector				   *ptfwave,
    InspiralTemplate           *params
    );

INT4 XLALInspiralComputePTFWDeriv (
    COMPLEX16Vector			   *Wderiv,
	REAL8FrequencySeries       *psd,
    InspiralTemplate           *params,
	INT4					   paramid,
	REAL8					   initdelta,
	REAL8					   tolerance
    );

INT4 XLALInspiralComputePTFQDeriv (
    REAL8VectorSequence		   *Qderiv,
    InspiralTemplate           *params
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRALBANK_H */
