/*************************** <lalVerbatim file="SimulateSBCV">
Author: Sukanta Bose (Adapted from a non-LAL code written by Bruce Allen)
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{SimulateSB.c}}
\label{stochastic:ss:SimulateSB.c}

Simulates whitened time-domain signal in a pair 
of interferometric detectors that arises purely from  an isotropic and
unpolarized stochastic background of gravitational radiation with the
desired power spectrum, $\Omega_{\scriptstyle{\rm GW}}(f)$.

\subsubsection*{Prototypes}
\input{SimulateSBCP}
\index{\texttt{LALSimulateSB()}}

\subsubsection*{Description}

The frequency domain strains $\widetilde{h}_1(f_i)$ 
and $\widetilde{h}_2(f_j)$ caused by
the stochastic background in two detectors are random variables that have
zero mean and that obey  \cite{simulatesb:Allen:1999}:
  \begin{equation}
    \langle\widetilde{h}_1^*(f_i)\widetilde{h}_1(f_j)\rangle
    = \frac{3H_0^2T}{20\pi^2}\delta_{ij}f_i^{-3}
    \Omega_{\scriptstyle{\rm GW}}(|f|)
  \end{equation}
  and
   \begin{equation}
    \langle\widetilde{h}_2^*(f_i)\widetilde{h}_2(f_j)\rangle
    = \frac{3H_0^2T}{20\pi^2}\delta_{ij}f_i^{-3}
    \Omega_{\scriptstyle{\rm GW}}(|f|)
  \end{equation}
  and
  \begin{equation}
    \langle\widetilde{h}_1^*(f_i)\widetilde{h}_2(f_j)\rangle
    = \frac{3H_0^2T}{20\pi^2}\delta_{ij}f_i^{-3}\gamma(f_i)
    \Omega_{\scriptstyle{\rm GW}}(|f|) \ ,
  \end{equation}
where $\langle\rangle$ denotes ensemble average, $T$ is the time of
observation, and $\gamma$ is the overlap reduction function
\cite{simulatesb:Flanagan:1993}. Above, $\widetilde{h}_1(f_i)$ and 
$\widetilde{h}_2(f_j)$ are the Fourier components of the gravitational strains
$h_1(t)$ and $h_2(t)$ at the two detectors. 

The Fourier components that 
obey the above relations are
  \begin{equation}
    \widetilde{h}_1(f_i) = \sqrt{3H_0^2T \over 40\pi^2}f_i^{-3/2}
    \Omega^{1/2}_{\scriptstyle{\rm GW}}(|f|) (x_{1i} + i y_{1i})
    \,
  \end{equation}
  and
  \begin{equation}
    \widetilde{h}_2(f_i) = \widetilde{h}_1(f_i)\gamma(f_i) +
    \sqrt{3H_0^2T \over 40\pi^2}f_i^{-3/2}
    \Omega^{1/2}_{\scriptstyle{\rm GW}}(|f|) 
    \sqrt{1-\gamma^2(f_i)} (x_{2i} + i y_{2i})
    \,
  \end{equation}
where $x_{1i}$, $y_{1i}$, $x_{2i}$, and $y_{2i}$ are statistically 
independent real Gaussian random variables, each of zero mean and unit 
variance. 

The routine assumes as inputs the data sample length, temporal spacing,
stochastic background characteristics, detector locations, the appropriate 
representations of the detector response function in each detector, etc.
The (frequency domain) response functions, $\widetilde{R}_1(f_i)$ and
$\widetilde{R}_2(f_i)$  are used to whiten the strains
$\widetilde{h}_1(f_i)$ and 
$\widetilde{h}_2(f_i)$, respectively, to obtain the whitened 
Fourier components:
  \begin{equation}
    \widetilde{o}_1(f_i) = \widetilde{R}_1(f_i)\widetilde{h}_1(f_i)
    \,
  \end{equation}
  and
  \begin{equation}
    \widetilde{o}_2(f_i) = \widetilde{R}_2(f_i)\widetilde{h}_2(f_i)
    \ .
  \end{equation}
To obtain the whitened (real) 
outputs $o_1(t_i)$ and $o_2(t_i)$ in the time domain, the inverse 
Fourier transforms of the above frequency series are taken.

\subsubsection*{Algorithm}

The routine \texttt{LALSimulateSB()} first inputs the frequency series 
describing the overlap reduction function consistent with the specified
pair of detectors. It uses this information, the specified power
spectrum for the stochastic background, and a random number generator (of
zero mean, unit variance Gaussian distributions) to generate 
$\widetilde{h}_1(f_i)$ and $\widetilde{h}_2(f_i)$. The 
response functions of the two detectors are then used to whiten the two 
strains in the Fourier domain. Their inverse transform is taken to obtain
at each detector the whitened simulated signal in the time domain.

\subsubsection*{Uses}

\begin{verbatim}
LALOverlapReductionFunction()
LALStochasticOmegaGW()
LALReverseRealFFT()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item This routine does not yet support non-zero heterodyning frequencies.
\end{itemize}

\vfill{\footnotesize\input{SimulateSBCV}}

******************************************************* </lalLaTeX> */ 

/**************************** <lalLaTeX file="SimulateSBCB">
\bibitem{simulate:Allen:1999}
  B.~Allen and J.~D.~Romano, ``Detecting a stochastic background of
  gravitational radiation: Signal processing strategies and
  sensitivities''
  Phys.\ Rev.\ D {\bf 59}, 102001 (1999);
  \href{http://www.arXiv.org/abs/gr-qc/9710117}{gr-qc/9710117}
  \bibitem{simulate:Flanagan:1993}
  E.~Flanagan, ``The sensitivity of the laser interferometer gravitational 
  wave observatory (LIGO) to a stochastic background, and its dependence on
  the detector orientations''
  Phys.\ Rev.\ D {\bf 48}, 2389 (1993);
  \href{http://www.arXiv.org/abs/astro-ph/9305029}{astro-ph/9305029}
******************************************************* </lalLaTeX> */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StochasticCrossCorrelation.h> 
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>
#include <lal/PrintVector.h>
#include <lal/Random.h>
#include <lal/SimulateSB.h> 

NRCSID (SIMULATESBC, "$Id$");

/* <lalVerbatim file="SimulateSBCP"> */
void
LALCreateSimulateSBOutput (
			   LALStatus                     *status,
			   SimulateSBOutput             **output,
			   SimulateSBInitParams          *params
			   )
     /* </lalVerbatim> */
{
  SimulateSBOutput          *outputPtr;

  INITSTATUS( status, "LALCreateSimulateSBOutput", SIMULATESBC );
  ATTATCHSTATUSPTR( status );

  /* check that the arguments are reasonable */
  /* parameter structure */
  ASSERT(params != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  ASSERT(params->length > 0, status, 
         SIMULATESBH_ENONPOSLEN,
         SIMULATESBH_MSGENONPOSLEN);

  /* output structure */
  ASSERT(output != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* create the output structure */
  outputPtr = *output = (SimulateSBOutput *)
    LALCalloc( 1, sizeof(SimulateSBOutput) );
  if ( ! outputPtr ) 
    { 
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
  /* create memory for the whitenedSimulatedSB1 structure */
  outputPtr->whitenedSimulatedSB1 = (REAL4TimeSeries *)
    LALCalloc( 1, sizeof(REAL4TimeSeries) );
  if ( ! outputPtr->whitenedSimulatedSB1 ) 
    { 
      LALFree( *output );
      *output = NULL;
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }

  /* create memory for the whitenedSimulatedSB1 data */
  LALSCreateVector(status->statusPtr, &(outputPtr->whitenedSimulatedSB1->data),
		   params->length);
  BEGINFAIL( status )
    {
      LALFree( outputPtr->whitenedSimulatedSB1 ); 
      outputPtr->whitenedSimulatedSB1 = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );

  /* create memory for the whitenedSimulatedSB2 structure */
  outputPtr->whitenedSimulatedSB2 = (REAL4TimeSeries *)
    LALCalloc( 1, sizeof(REAL4TimeSeries) );
  if ( ! outputPtr->whitenedSimulatedSB2 ) 
    { 
      LALFree( outputPtr->whitenedSimulatedSB1->data );
      LALFree( outputPtr->whitenedSimulatedSB1 );
      outputPtr->whitenedSimulatedSB1 = NULL;
      LALFree( *output );
      *output = NULL;
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
    /* create memory for the whitenedSimulatedSB2 data */
  LALSCreateVector(status->statusPtr, &(outputPtr->whitenedSimulatedSB2->data),
		   params->length);
  BEGINFAIL( status )
    {
      LALFree( outputPtr->whitenedSimulatedSB2 );
      LALFree( outputPtr->whitenedSimulatedSB1->data );
      LALFree( outputPtr->whitenedSimulatedSB1 );
      outputPtr->whitenedSimulatedSB1 = NULL;
      outputPtr->whitenedSimulatedSB2 = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="SimulateSBCP"> */
void
LALDestroySimulateSBOutput (
			    LALStatus                      *status,
			    SimulateSBOutput              **output
			    )
/* </lalVerbatim> */
{
  SimulateSBOutput          *outputPtr;

  INITSTATUS( status, "LALDestroySimulateSBOutput", SIMULATESBC );
  ATTATCHSTATUSPTR( status );
  
  /* check that the arguments are reasonable */
  
  /* output structure */
  ASSERT(output != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* destroy the output data storage */
  outputPtr = *output;
  LALSDestroyVector(status->statusPtr, 
		    &(outputPtr->whitenedSimulatedSB1->data) );
  CHECKSTATUSPTR( status );

  LALSDestroyVector(status->statusPtr, 
		    &(outputPtr->whitenedSimulatedSB2->data) );
  CHECKSTATUSPTR( status );

  /* destroy the output whitenedSimulatedSB structures */
  LALFree( outputPtr->whitenedSimulatedSB1 );
  LALFree( outputPtr->whitenedSimulatedSB2 );

  /* destroy the output structure */
  LALFree( outputPtr );
  *output = NULL;
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="SimulateSBCP"> */
void
LALSimulateSBInit (
		   LALStatus                     *status,
		   SimulateSBParams             **output,
		   SimulateSBInitParams          *params
		   )
     /* </lalVerbatim> */
{
  SimulateSBParams          *outputPtr;
  
  INITSTATUS( status, "LALSimulateSBInit", SIMULATESBC );
  ATTATCHSTATUSPTR( status );

  /* check that the arguments are reasonable */

  /* output structure */
  ASSERT(output != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
 
  /* parameter structure */
  ASSERT(params != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* ensure that number of points in time segment is positive */
  ASSERT(params->length > 0 , status, 
         SIMULATESBH_ENONPOSLEN,
         SIMULATESBH_MSGENONPOSLEN);

  /* create the output structure */
  outputPtr = *output = (SimulateSBParams *)
    LALCalloc( 1, sizeof(SimulateSBParams) );
  if ( ! outputPtr ) 
    { 
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="SimulateSBCP"> */
void
LALSimulateSBFinalize (
		       LALStatus                     *status,
		       SimulateSBParams             **output
		       )
     /* </lalVerbatim> */
{
  SimulateSBParams          *outputPtr;
  
  INITSTATUS( status, "LALSimulateSBFinalize", SIMULATESBC );
  ATTATCHSTATUSPTR( status );
  
  /* check that the arguments are reasonable */
  
  /* output structure */
  ASSERT(output != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* local pointer to output structure */
  outputPtr = *output;
  
  /* destroy the output structure */
  LALFree( outputPtr );
  *output = NULL;
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="SimulateSBCP"> */
void
LALCreateSimulateSBInput (
			  LALStatus                     *status,
			  SimulateSBInput              **output,
			  SimulateSBInitParams          *params
			  )
     /* </lalVerbatim> */
{
  SimulateSBInput          *outputPtr;
  
  INITSTATUS( status, "LALCreateSimulateSBInput", SIMULATESBC );
  ATTATCHSTATUSPTR( status );
  
  /* check that the arguments are reasonable */
  
  /* output structure */
  ASSERT(output != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* create the output structure */
  outputPtr = *output = (SimulateSBInput *)
    LALCalloc( 1, sizeof(SimulateSBInput) );
  if ( ! outputPtr ) 
    { 
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
  /* create memory for the overlapReductionFunction structure */ 
  outputPtr->overlapReductionFunction = (REAL4FrequencySeries *)
    LALCalloc( 1, sizeof(REAL4FrequencySeries) );
  if ( ! outputPtr->overlapReductionFunction ) 
    { 
      LALFree( *output );
      *output = NULL;
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
  /* create memory for the omegaGW structure */
  outputPtr->omegaGW = (REAL4FrequencySeries *)
    LALCalloc( 1, sizeof(REAL4FrequencySeries) );
  if ( ! outputPtr->omegaGW ) 
    { 
      LALFree( *output );
      *output = NULL;
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
   
  /* create memory for the overlapReductionFunction data */
  LALSCreateVector(status->statusPtr, 
		   &(outputPtr->overlapReductionFunction->data),
		   ((params->length)/2 + 1));
  BEGINFAIL( status )
    {
      LALFree( outputPtr->overlapReductionFunction );
      outputPtr->overlapReductionFunction = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );
  
  /* create memory for the omegaGW data */
  LALSCreateVector(status->statusPtr, &(outputPtr->omegaGW->data),
		   ((params->length)/2 + 1) );
  BEGINFAIL( status )
    {
      LALFree( outputPtr->omegaGW );
      outputPtr->omegaGW = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );
  
  /* create memory for the whiteningFilter1 structure */
  outputPtr->whiteningFilter1 = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! outputPtr->whiteningFilter1 ) 
    { 
      LALFree( *output );
      *output = NULL;
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
  /* create memory for the whiteningFilter2 structure */
  outputPtr->whiteningFilter2 = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! outputPtr->whiteningFilter2 ) 
    { 
      LALFree( *output );
      *output = NULL;
      ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC );
    }
  
  
  /* create memory for the whiteningFilter1 data */
  LALCCreateVector(status->statusPtr, &(outputPtr->whiteningFilter1->data),
		   ((params->length)/2 + 1));
  BEGINFAIL( status )
    {
      LALFree( outputPtr->whiteningFilter1 );
      outputPtr->whiteningFilter1 = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );
  
  /* create memory for the whiteningFilter2 data */
  LALCCreateVector(status->statusPtr, &(outputPtr->whiteningFilter2->data),
		   ((params->length)/2 + 1) );
  BEGINFAIL( status )
    {
      LALFree( outputPtr->whiteningFilter2 );
      outputPtr->whiteningFilter2 = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="SimulateSBCP"> */
void
LALDestroySimulateSBInput (
			   LALStatus                      *status,
			   SimulateSBInput               **output
			   )
     /* </lalVerbatim> */
{
  SimulateSBInput          *outputPtr;
  
  INITSTATUS( status, "LALDestroySimulateSBInput", SIMULATESBC );
  ATTATCHSTATUSPTR( status );
  
  /* check that the arguments are reasonable */
  
  /* output structure */
  ASSERT(output != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* destroy the output data storage */
  outputPtr = *output;
  LALSDestroyVector(status->statusPtr, 
		    &(outputPtr->overlapReductionFunction->data) );
  CHECKSTATUSPTR( status );

  LALSDestroyVector(status->statusPtr, 
		    &(outputPtr->omegaGW->data) );
  CHECKSTATUSPTR( status );

  LALCDestroyVector(status->statusPtr, 
		    &(outputPtr->whiteningFilter1->data) );
  CHECKSTATUSPTR( status );

  LALCDestroyVector(status->statusPtr, 
		    &(outputPtr->whiteningFilter2->data) );
  CHECKSTATUSPTR( status );

  /* destroy the output whiteningFilter structures */
  LALFree( outputPtr->overlapReductionFunction );
  LALFree( outputPtr->omegaGW );
  LALFree( outputPtr->whiteningFilter1 );
  LALFree( outputPtr->whiteningFilter2 );
  
  /* destroy the output structure */
  LALFree( outputPtr );
  *output = NULL;
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="SimulateSBCP"> */
void
LALSimulateSB( LALStatus                  *status,
	       SimulateSBOutput          **output,
	       SimulateSBInput            *input,
	       SimulateSBParams           *params 
	       )
     /* </lalVerbatim> */
{
  /* output */
  SimulateSBOutput *outputPtr=NULL;
  REAL4Vector      *outData[2]={NULL,NULL};
  
  /* parameters */
  UINT4             length;   /* (time) length of output vector data samples */
  UINT4             freqlen;
  REAL8             deltaT;   /* time spacing */
  REAL8             f0;       /* start frequency */
  REAL8             fRef;    /* reference normalization frequency */
  
  /* counters */ 
  UINT4             i;
  UINT4             j;
  
  /* other variables used */
  REAL8             deltaF;
  RandomParams     *randParams=NULL;
  INT4              seed=123;         
  
  /* vector for storing random numbers */ 
  REAL4Vector      *gaussdevs=NULL;
  
  /* IFO output counts in freq domain : */
  COMPLEX8Vector  *ccounts[2]={NULL,NULL};
  COMPLEX8Vector  *ccountsTmp[2]={NULL,NULL};
  
  /* Plan for reverse FFTs */ 
  RealFFTPlan      *invPlan=NULL;  
  
  /* initialize status pointer */
  INITSTATUS(status, "LALSimulateSB", SIMULATESBC);
  ATTATCHSTATUSPTR(status);
  
  /*
   *
   *ERROR CHECKING
   *
   */
  
  /***** check input/output structures exist *****/

  /* output structure */
  ASSERT(output, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* input structure */
  ASSERT(input != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* overlap member of input */
  ASSERT(input->overlapReductionFunction != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* omega member of input */
  ASSERT(input->omegaGW != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* detectorpair member of input */
  ASSERT(input->detectors != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);
  
  /* First detector's complex response (whitening filter) part of input */
  ASSERT(input->whiteningFilter1 != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* Second detector's complex response (whitening filter) part of input */
  ASSERT(input->whiteningFilter2 != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* data member of overlap */
  ASSERT(input->overlapReductionFunction->data != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* data member of omega */
  ASSERT(input->omegaGW->data != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* data member of first detector's response (whitening filter) */
  ASSERT(input->whiteningFilter1->data != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* data member of second detector's response (whitening filter) */
  ASSERT(input->whiteningFilter2->data != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* data-data member of first detector's response (whitening filter) */
  ASSERT(input->whiteningFilter1->data->data != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /* data-data member of second detector's response (whitening filter) */
  ASSERT(input->whiteningFilter2->data->data != NULL, status, 
         SIMULATESBH_ENULLP,
         SIMULATESBH_MSGENULLP);

  /************* check parameter structures ***********/

  /*No. of discrete time samples (length) is non-zero in each detector output*/
  ASSERT(params->length > 0, status, 
         SIMULATESBH_ENONPOSLEN,
         SIMULATESBH_MSGENONPOSLEN);

  /* time-interval between successive samples is non-zero in each output */
  ASSERT(params->deltaT > 0, status, 
         SIMULATESBH_ENONPOSLEN,
         SIMULATESBH_MSGENONPOSLEN);

  /************* done with null pointers *****************/
  
  
  /**** check for legality ****/
  
  /* start frequency must not be negative */
  f0 = params->f0;
  if (f0 < 0)
    {
      ABORT( status,
	     SIMULATESBH_ENEGFMIN,
	     SIMULATESBH_MSGENEGFMIN );
    }
  
  
  /** check for mismatches **/
  /* frequency length = length/2 +1  */
  length = params->length;
  freqlen = length/2 +1;
  if (input->overlapReductionFunction->data->length != (length/2 +1)) 
    {
      ABORT(status,
	    SIMULATESBH_EMMLEN,
	    SIMULATESBH_MSGEMMLEN);
    }
  if (input->omegaGW->data->length != (length/2 +1)) 
    {
      ABORT(status,
	    SIMULATESBH_EMMLEN,
	    SIMULATESBH_MSGEMMLEN);
    }
  if (input->whiteningFilter1->data->length != (length/2 +1)) 
    {
      ABORT(status,
	    SIMULATESBH_EMMLEN,
	    SIMULATESBH_MSGEMMLEN);
    }
  if (input->whiteningFilter2->data->length != (length/2 +1)) 
    {
      ABORT(status,
	    SIMULATESBH_EMMLEN,
	    SIMULATESBH_MSGEMMLEN);
    }
  if (input->whiteningFilter1->f0 != f0) 
    {
      ABORT(status,
	    SIMULATESBH_EMMFMIN,
	    SIMULATESBH_MSGEMMFMIN);
    }
  if (input->whiteningFilter2->f0 != f0) 
    {
      ABORT(status,
	    SIMULATESBH_EMMFMIN,
	    SIMULATESBH_MSGEMMFMIN);
    }

  /* frequency spacing */
  deltaT = params->deltaT;
  deltaF = 1/(deltaT*length);
  if (input->whiteningFilter1->deltaF != deltaF) 
    {
      ABORT(status,
	    SIMULATESBH_EMMDELTAF,
	    SIMULATESBH_MSGEMMDELTAF);
    }
  if (input->whiteningFilter2->deltaF != deltaF) 
    {
      ABORT(status,
	    SIMULATESBH_EMMDELTAF,
	    SIMULATESBH_MSGEMMDELTAF);
    }
  
  
  /** check for reference frequency lower and upper limits **/ 
  fRef = params->fRef;
  if ( fRef < deltaF )
    {
      ABORT(status,
	    SIMULATESBH_EOORFREF,
	    SIMULATESBH_MSGEOORFREF);
    }
  if ( fRef > ((length-1)*deltaF) )
    {
      ABORT(status,
	    SIMULATESBH_EOORFREF,
	    SIMULATESBH_MSGEOORFREF);
    }

  /*
   *
   *EVERYHTING OKAY HERE 
   *
   */
  

  /******** create fft plans and workspace vectors *****/

  /* fft plan */
  LALCreateReverseRealFFTPlan(status->statusPtr,&invPlan,length,0);
  CHECKSTATUSPTR( status); 
  
  LALSCreateVector( status->statusPtr, 
		    &gaussdevs, 4*length ); 
  CHECKSTATUSPTR( status); 
  
  /* create parameters for generating random numbers from seed */
  LALCreateRandomParams( status->statusPtr, 
			 &randParams, seed ); 
  CHECKSTATUSPTR( status);
  
  /* create random numbers from parameters */
  LALNormalDeviates( status->statusPtr, 
		     gaussdevs, randParams ); 
  CHECKSTATUSPTR( status);
  
  for (i=0;i<2;i++)
    {
    LALCCreateVector(status->statusPtr, &ccounts[i],freqlen);
    LALCCreateVector(status->statusPtr, &ccountsTmp[i],freqlen);
  }
  CHECKSTATUSPTR( status); 
  
  if (f0 == 0)
    {
      REAL4    gamma;
      REAL4    omega;
      REAL8    freq;
      REAL8    factor;
      REAL8    factor2;
      COMPLEX8 wFilter1;
      COMPLEX8 wFilter2;
       
      /* loop over frequencies; will do DC and Nyquist below */
      j=0;
      for (i = 1; i < freqlen; ++i)
	{
	  gamma = input->overlapReductionFunction->data->data[i];
	  omega = input->omegaGW->data->data[i];
	  freq  = i*deltaF;	  
	  factor = deltaF * sqrt(3.0L * length * deltaT * omega / 
				 (40.0L *freq*freq*freq)
				 )* LAL_H0FAC_SI / LAL_PI;
	  factor2 = sqrt(1-gamma*gamma)*factor;
	  
	  wFilter1 = input->whiteningFilter1->data->data[i];
	  wFilter2 = input->whiteningFilter2->data->data[i];
	  
	  ccountsTmp[0]->data[i].re=factor*gaussdevs->data[j++];
	  ccountsTmp[0]->data[i].im=factor*gaussdevs->data[j++];
	  ccountsTmp[1]->data[i].re=ccountsTmp[0]->data[i].re*gamma+
	    factor2*gaussdevs->data[j++];
	  ccountsTmp[1]->data[i].im=ccountsTmp[0]->data[i].im*gamma+
	    factor2*gaussdevs->data[j++];
	  
	  ccounts[0]->data[i].re = wFilter1.re * ccountsTmp[0]->data[i].re -
	    wFilter1.im * ccountsTmp[0]->data[i].im;
	  ccounts[0]->data[i].im = wFilter1.re * ccountsTmp[0]->data[i].im +
	    wFilter1.im * ccountsTmp[0]->data[i].re;
	  ccounts[1]->data[i].re = wFilter2.re * ccountsTmp[1]->data[i].re -
	    wFilter2.im * ccountsTmp[1]->data[i].im;
	  ccounts[1]->data[i].im = wFilter2.re * ccountsTmp[1]->data[i].im +
	    wFilter2.im * ccountsTmp[1]->data[i].re;
	}
      
      /* Set DC, Nyquist (imaginary) components to zero */
      for (i=0;i<2;++i)
	{
	  ccounts[i]->data[0].re=0.0;
	  ccounts[i]->data[0].im=0.0;
	  ccountsTmp[i]->data[length/2].im=0.0;
	}
      
      /* Compute the whitened Nyquist (real) component */
      gamma = input->overlapReductionFunction->data->data[length/2];
      omega = input->omegaGW->data->data[length/2];
      freq = deltaF*length/2;
      
      wFilter1 = input->whiteningFilter1->data->data[length/2];
      wFilter2 = input->whiteningFilter2->data->data[length/2];
      
      /* check that whitening filter is real in time domain */
      if (wFilter1.im != 0) 
	{
	  ABORT(status,
		SIMULATESBH_ECOMPTIME,
		SIMULATESBH_MSGECOMPTIME);
	};
      if (wFilter2.im != 0) 
	{
	  ABORT(status,
		SIMULATESBH_ECOMPTIME,
		SIMULATESBH_MSGECOMPTIME);
	};
      
      factor = deltaF * sqrt(3.0L * length * deltaT * omega / 
			     (40.0L *freq*freq*freq)
			     )* LAL_H0FAC_SI / LAL_PI;
      factor2 = sqrt(1-gamma*gamma)*factor;
      
      ccountsTmp[0]->data[length/2].re=factor*gaussdevs->data[j++];
      
      ccountsTmp[1]->data[length/2].re=
	(ccountsTmp[0]->data[length/2].re*gamma + 
	 factor2*gaussdevs->data[j++]);
      
      ccounts[0]->data[length/2].re = 
	(wFilter1.re * ccountsTmp[0]->data[length/2].re - 
	 wFilter1.im * ccountsTmp[0]->data[length/2].im);
      ccounts[0]->data[length/2].im = 0;
      
      
      ccounts[1]->data[length/2].re = 
	(wFilter2.re * ccountsTmp[1]->data[length/2].re - 
	 wFilter2.im * ccountsTmp[1]->data[length/2].im);
      ccounts[1]->data[length/2].im = 0;
      
      for (i=0;i<2;i++)
	{
	  LALSCreateVector(status->statusPtr, 
			   &outData[i],length);
	}
      CHECKSTATUSPTR( status );
     
      /* InvFFT from freq to time domain & get output (no detector noise) */
      for (i=0;i<2;i++)
	{
	  LALReverseRealFFT(status->statusPtr,outData[i],
			    ccounts[i],invPlan);
	}
      CHECKSTATUSPTR( status );
      
      /*
       * 
       * assign parameters to output 
       *
       */
      
      outputPtr = *output = (SimulateSBOutput *)
	LALCalloc(1, sizeof(SimulateSBOutput));
      if( ! outputPtr )
	{ 
	  ABORT( status, SIMULATESBH_EALOC, SIMULATESBH_MSGEALOC);
	}
      
      outputPtr->whitenedSimulatedSB1->f0                   = f0;
      outputPtr->whitenedSimulatedSB1->deltaT               = deltaT;
      outputPtr->whitenedSimulatedSB1->epoch.gpsSeconds     = 0;
      outputPtr->whitenedSimulatedSB1->epoch.gpsNanoSeconds = 0;
      outputPtr->whitenedSimulatedSB1->data                 = outData[0];
      outputPtr->whitenedSimulatedSB2->f0                   = 0; 
      outputPtr->whitenedSimulatedSB2->deltaT               = deltaT;
      outputPtr->whitenedSimulatedSB2->epoch.gpsSeconds     = 0;
      outputPtr->whitenedSimulatedSB2->epoch.gpsNanoSeconds = 0;
      outputPtr->whitenedSimulatedSB2->data                 = outData[1];
      
      
      strncpy( outputPtr->whitenedSimulatedSB1->name, 
	       "Simulated stochastic background One", LALNameLength );
      strncpy( outputPtr->whitenedSimulatedSB2->name, 
	       "Simulated stochastic background Two", LALNameLength );
      
    } /* if (f0 == 0) */
  else
    {
      
      /*****This abort should be replaced with the correct *****/
      /***** non-zero heterodyne frequency procedure******/
      ABORT(status, SIMULATESBH_ENOTYETHETERO,
	    SIMULATESBH_MSGENOTYETHETERO); 
    }
  
  
  /* normal exit */
  DETATCHSTATUSPTR(status);
  RETURN(status);
  
}
