/* <lalVerbatim file="LALInspiralBankHV">

Author: Churches, D.K. and Sathyaprakash, B.S.
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
#define LALINSPIRALBANKH_MSGENULL   "Null pointer"
#define LALINSPIRALBANKH_MSGEMEM    "Memory allocation failure"
#define LALINSPIRALBANKH_MSGECHOICE "Invalid choice for an input parameter"
#define LALINSPIRALBANKH_MSGEDIV0   "Division by zero"
#define LALINSPIRALBANKH_MSGESIZE   "Invalid input range"
#define LALINSPIRALBANKH_MSGEFRANGE "Limits outside range of frequency series"
#define LALINSPIRALBANKH_MSGEORDER  "Inappropriate PN order"
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
was added to handle BCV templates.

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
\item \texttt{alpha}: the BCV amplitude correction parameter
\item \texttt{psi0Min}: minimum value of the parameter $\psi_0$
\item \texttt{psi0Max}: maximum value of the parameter $\psi_0$
\item \texttt{psi3Min}: minimum value of the parameter $\psi_3$
\item \texttt{psi3Max}: maximum value of the parameter $\psi_3$
\item \texttt{mMin}: minimum mass of components to search for 
\item \texttt{mMax}: maximum mass of components to search for
\item \texttt{MMax}:   alternatively, maximum total mass of binary to search for
\item \texttt{mmCoarse}:  Coarse grid minimal match 
\item \texttt{mmFine}:  Fine grid minimal match 
\item \texttt{fLower}:  Lower frequency cutoff
\item \texttt{fUpper}:  Upper frequency cutoff 
\item \texttt{tSampling}:  Sampling rate
\item \texttt{etamin}: minimum value of eta in our search 
\item \texttt{shf}: Frequency series containing the PSD 
\item \texttt{iflso}: (currently not implemented) flso will be used as an 
\item \texttt{numFcutTemplates}: number of templates in the {\tt fcut} direction

The next two members are used in setting up the InspiralTemplate
parameter structure but not in creating the template bank. 

\item \texttt{order}: Post-Newtonian order of the waveform 
\item \texttt{approximant}: Approximant of the waveform 
\end{itemize}

\item \texttt{InspiralFineBankIn}
Structre needed by the function \texttt{LALInspiralCreateFineBank}.
	which computes a finer mesh around a given lattice point
	using the value of the fine-mesh minimal match, coarse-mesh
	minimal match and the metric at the current lattice point.
\input{LALInspiralFineBankInH}
\begin{itemize}
\item {templateList:} A list contianing all the fine-mesh templates
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
\input{LALInspiralParamsH}
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
therefore defined an oversized (by one element) array.
\end{itemize}

\item {\texttt{RectangleIn} and \texttt{RectangleOut}:}
Input and ouput structures to function LALRectangleVertices.
\input{LALRectangleInH}
\input{LALRectangleOutH}

\end{enumerate}
</lalLaTeX>  */


/* <lalVerbatim file="LALCoordinateSpaceH"> */
typedef enum
{
  Tau0Tau2, 
  Tau0Tau3, 
  Psi0Psi3
}
CoordinateSpace;
/*  </lalVerbatim>  */
/*  <lalLaTeX> 
\idx[Type]{CoordinateSpace} 
</lalLaTeX>  */


/* <lalVerbatim file="LALGridTypeH"> */
typedef enum
{
  Square, 
  OrientedSquare,
  Hexagonal,
  OrientedHexagonal
}
GridType;
/*  </lalVerbatim>  */
/*  <lalLaTeX> 
\idx[Type]{GridType} 
</lalLaTeX>  */






/* <lalVerbatim file="LALInspiralBankMassRangeH"> */
typedef enum
{
  MinComponentMassMaxTotalMass,
  MinMaxComponentMass
} 
InspiralBankMassRange;
/*  </lalVerbatim>  */
/*  <lalLaTeX> 
\idx[Type]{InspiralBankMassRange} 
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

/*  <lalVerbatim file="LALInspiralParamsH"> */
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

  REAL8FrequencySeries          shf;

  INT4                          iflso;          
  UINT4                         numFcutTemplates;
  REAL4				HighGM;
  REAL4				LowGM;

  GridType                      gridType;
  Order                         order;        
  Approximant                   approximant;  
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
  REAL8 x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;
}
RectangleOut;
/*  </lalVerbatim>  */
/*  <lalLaTeX>
\index{\texttt{RectangleOut}}
</lalLaTeX>  */

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

void 
LALInspiralCreateBCVBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void 
LALInspiralBCVBankFcutS3 (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                *NList, 
    InspiralCoarseBankIn coarseIn);

void
LALInspiralBCVFcutBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                *NList, 
    InspiralCoarseBankIn coarseIn);

void 
LALInspiralBCVRegularFcutBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                *NList, 
    InspiralCoarseBankIn coarseIn);


/* <lalLaTeX>
   \newpage\input{InspiralSpinBankC}
   </lalLaTeX>  */

void
LALInspiralSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable   **tiles,
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
LALInspiralValidTemplate
(
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
LALInverse3
(
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
LALRectangleVertices
(
   LALStatus *status, 
   RectangleOut *out,
   RectangleIn *in
);


void 
LALEmpiricalPSI2MassesConversion(
    InspiralTemplate    *params,
    UINT4               *valid,
    REAL4               lightring
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

/*<lalLaTeX>
\newpage\input{InspiralSpinBankTestC}
</lalLaTeX>*/

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRALBANK_H */
