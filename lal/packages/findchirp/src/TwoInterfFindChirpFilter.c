/*************************** <lalVerbatim file="TwoInterfFindChirpFilterCV">
Author: Sukanta Bose
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{TwoInterfFindChirpFilter.c}}
\label{twointerffindchirp:ss:TwoInterfFindChirpFilter.c}

Computes the search statistic for an inspiral waveform in the outputs of a 
pair of interferometers.

\subsubsection*{Prototypes}
%\input{TwoInterfFindChirpFilterCP}
\idx{LALTwoInterfFindChirpFilter()}

\subsubsection*{Description}

This package computes the coherent statistic 
\cite{twointerffindchirp:BPD:2000} as a function of time, given a 
post-Newtonian template and coincident data from two interferometers.
\vfill{\footnotesize\input{TwoInterfFindChirpFilterCV}}

For a network of $M$ detectors, the data consist of $M$ data trains, 
$\{x^I(t)|\>I=1,2,...,M$ and $t \in [0,T]\}$. The network matched-template can 
be obtained naturally by the maximum-likelihood method, where the 
decision whether
the signal is present or not is made by evaluating the likelihood ratio (LR)
for the network \cite{twointerffindchirp:Hels:1968}. 
For independent detector noises, it is straightforward to show that the 
network LR, denoted by $\lambda$, is just a product of the individual 
detector LRs. (For non-independent detector noises, see Ref. 
\cite{twointerffindchirp:Finn:2001}.)
In addition, for Gaussian noise, the logarithmic likelihood 
ratio (LLR) for the network is just the sum of the LLRs of the individual 
detectors \cite{twointerffindchirp:BDP:1999},
\begin{equation}
\ln\lambda = \sum_{I=1}^M \ln\lambda_{(I)} \ \ ,
\end{equation}
where 
\begin{equation}
\ln\lambda_{(I)} = \langle s^I, x^I \rangle_{(I)} -
{1 \over 2} \langle s^I, s^I \rangle_{(I)} \,.
\end{equation}
Above, we have denoted the cross-correlation between two vectors by the
inner-product:
\begin{equation}
\label{twointerffindchirp:e:innerProduct}
\langle a ,\>b \rangle_{(I)} = 2 \Re \int_{0}^{\infty} \! df\>
{\tilde{a}^* (f) \tilde{b}(f) \over s_{h(I)}(f)} \ \ , 
\end{equation}
where $\tilde{a}(f)$ and $\tilde{b}(f)$ are the Fourier
transforms of $a(t)$ and $b(t)$, respectively, and $s_{h(I)}(f)$ is the
two-sided noise power spectral density of the $I$-th detector.

The network LLR takes a compact form in terms of the network inner-product,
\begin{equation}
\langle {\sf s},{\sf x}\rangle_{NW} = \sum_{I=1}^M \langle s^{I}(t), x^{I}
(t)\rangle_{(I)} \ \ ,
\end{equation}
where 
\begin{equation}
\label{twointerffindchirp:e:netSignal}
{\sf s} (t) =\left(s^{1}(t), s^{2}(t),\cdots , s^{M} (t) \right)
\end{equation}
is the network template-vector, which comprises of individual 
detector-templates as its components, and
\begin{equation}
\label{twointerffindchirp:e:netData}
{\sf x} (t) =\left(x^{1} (t), x^{2}(t),\cdots , x^{M} (t) \right)
\end{equation}
is the network data-vector. 

It can be shown by using the Schwarz inequality 
that the network template, ${\sf s}$, defined 
above yields the maximum signal-to-noise (SNR) amongst all linear templates 
and, hence, is the matched template. As shown in Ref. 
\cite{twointerffindchirp:BDP:1999}, in terms of 
the above definitions, the network LLR takes the following simple form:
\begin{equation}\label{twointerffindchirp:e:netMLR}
\ln\lambda = \langle {\sf s},{\sf x}\rangle_{NW} -
{1 \over 2} \langle {\sf s}, {\sf s}\rangle_{NW} \ \ ,
\end{equation}
which is a function of the source parameters that determine ${\sf s}$. Given 
${\sf s}$, different selections of source-parameter values and, therefore,
different values of ${\sf s}$ result in varying magnitudes of the LLR. The
selection that gives the maximum value stands the best chance for beating the
pre-set threshold on the LLR. Since scanning the complete source-parameter 
manifold for the maximum of LLR is computationally very expensive
\cite{twointerffindchirp:PDB:2001}, it is wiser
to perform its maximization analytically over as many parameters as possible.
Such a maximization can be done over four parameters, namely, the 
luminosity distance to the source, the waveform phase at the time of 
final coalescence, the inclination angle of the binary's orbit, and the
polarization-ellipse angle. This package computes for a pair of 
interferometers the value of the reduced statistic 
that is derived from the above statistic by such a maximization
\cite{twointerffindchirp:BPD:2000}.

\subsubsection*{Uses}

\begin{verbatim}
LALFindChirpChisqVeto()
LALCOMPLEX8VectorFFT()
\end{verbatim}

\vfill{\footnotesize\input{TwoInterfFindChirpFilterCV}}

******************************************************* </lalLaTeX> */ 

/**************************** <lalLaTeX file="TwoInterfFindChirpFilterCB">
\bibitem{twointerffindchirp:BPD:2000}
 S.~Bose, A.~Pai, and S.~V.~Dhurandhar, \/ Int.\ J.\ Mod.\ Phys. \ D\
 {\bf 9}, 325 (2000).   
 \href{http://www.arXiv.org/abs/gr-qc/0002010}{gr-qc/0002010}
\bibitem{twointerffindchirp:Hels:1968} 
 C.~W.~Helstrom, \textit{Statistical Theory of Signal Detection}, 
 Pergamon Press, London (1968).
\bibitem{twointerffindchirp:Finn:2001} 
 L.~S.~Finn, Phys.\ Rev.\ D\ {\bf 63}, 102001 (2001).
 \href{http://www.arXiv.org/abs/gr-qc/0010033}{gr-qc/0010033}
\bibitem{twointerffindchirp:BDP:1999}
 S.~Bose, S.~V.~Dhurandhar, and A.~Pai, \/ Pramana\ {\bf 53}, 1125 (1999). 
 \href{http://www.arXiv.org/abs/gr-qc/9906064}{gr-qc/9906064}
\bibitem{twointerffindchirp:PDB:2001}
 A.~Pai, S.~V.~Dhurandhar, and S.~Bose, \/ Phys.\ Rev.\ D{\bf 64}, 042004
 (2001).
 \href{http://www.arXiv.org/abs/gr-qc/0009078}{gr-qc/0009078}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h> 
#include <lal/TwoInterfFindChirp.h> 

double rint(double x);

NRCSID (TWOINTERFFINDCHIRPFILTERC, "$Id$");

static REAL4 cartesianInnerProduct(REAL4 x[3], REAL4 y[3])
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
void
LALCreateTwoInterfFindChirpFilterInput (
					LALStatus                      *status,
					TwoInterfFindChirpFilterInput **output,
					TwoInterfFindChirpInitParams   *params
					)
{
  TwoInterfFindChirpFilterInput         *outputPtr;
  
  INITSTATUS( status, "LALCreateTwoInterfFindChirpFilterInput", 
	      TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * ensure that the arguments are reasonable
   *
   */
  

  /* check that the output handle exists, but points to a null pointer */
  ASSERT( output, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );  
  ASSERT( !*output, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );

  /* check that the parameter structure exists */
  ASSERT( params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* check that the parameter sub-structures exists */
  ASSERT( params->initParams1, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params->initParams2, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* check that the number of points is positive */
  ASSERT( params->initParams1->numPoints > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );

  ASSERT( params->initParams2->numPoints > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );

  /** check for mismatches **/

  if ( params->initParams1->numPoints !=  params->initParams2->numPoints ) 
    {
      ABORT(status,
	    TWOINTERFFINDCHIRPH_ENUMF,
	    TWOINTERFFINDCHIRPH_MSGENUMF);
    }
  

  /*
   *
   * create the twointerffindchirp filter input structure
   *
   */


  /* create the output structure */
  outputPtr = *output = (TwoInterfFindChirpFilterInput *)
    LALCalloc( 1, sizeof(TwoInterfFindChirpFilterInput) );
  if ( ! outputPtr )
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  /* create memory for the chirp template structure */
  outputPtr->fcTmplt = (FindChirpTemplate *)
    LALCalloc( 1, sizeof(FindChirpTemplate) );
  if ( !outputPtr->fcTmplt )
    {
      LALFree( *output );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  /* create memory for the chirp template data */
  LALCCreateVector (status->statusPtr, &(outputPtr->fcTmplt->data), 
		    params->numPoints);
  BEGINFAIL( status )
    {
      LALFree( outputPtr->fcTmplt );
      outputPtr->fcTmplt = NULL;
      LALFree( *output );
      *output = NULL;
    }
  ENDFAIL( status );
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALDestroyTwoInterfFindChirpInput (
				   LALStatus                          *status,
				   TwoInterfFindChirpFilterInput     **output
				   )
{
  TwoInterfFindChirpFilterInput         *outputPtr;

  INITSTATUS( status, "LALDestroyTwoInterfFindChirpFilterInput", TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  
  /*
   *
   * ensure that the arguments are reasonable
   *
   */


  /* check that handle is non-null and points to a non-null pointer */
  ASSERT( output, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *output, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  
  /*
   *
   * destroy the twointerffindchirp input structure
   *
   */
  

  /* local pointer to output */
  outputPtr = *output;
  
  /* destroy the chirp template data storage */
  LALCDestroyVector( status->statusPtr, &(outputPtr->fcTmplt->data) );
  CHECKSTATUSPTR( status );
  
  /* destroy the chirp template structure */
  LALFree( outputPtr->fcTmplt );
  
  /* destroy the filter input structure */
  LALFree( outputPtr );
  *output = NULL;
  

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALTwoInterfFindChirpFilterInit (
				 LALStatus                           *status,
				 TwoInterfFindChirpFilterParams     **output,
				 TwoInterfFindChirpInitParams        *params
				 )
{
  TwoInterfFindChirpFilterParams        *outputPtr;
  
  INITSTATUS( status, "LALTwoInterfFindChirpFilterInit", TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  
  
  /* check that the output handle exists, but points to a null pointer */
  ASSERT( output, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*output, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );
  
  /* check that the parameter structure exists */
  ASSERT( params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the parameter sub-structure exists */
  ASSERT( params->initParams1, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  ASSERT( params->initParams2, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0,  status, 
	  TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );
  
  ASSERT( params->initParams1->numPoints > 0,  status, 
	  TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );
  
  ASSERT( params->initParams2->numPoints > 0,  status, 
	  TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );
  
  /** check for mismatches **/
  
  if ( params->initParams1->numPoints !=  params->initParams2->numPoints ) 
    {
      ABORT(status,
	    TWOINTERFFINDCHIRPH_ENUMF,
	    TWOINTERFFINDCHIRPH_MSGENUMF);
    }
  
  
  if ( params->numPoints !=  params->initParams1->numPoints ) 
    {
      ABORT(status,
	    TWOINTERFFINDCHIRPH_ENUMF,
	    TWOINTERFFINDCHIRPH_MSGENUMF);
    }
  
  
  /*
   *
   * allocate memory for the TwoInterfFindChirpFilterParams
   *
   */
  
  
  /* create the output structure */
  outputPtr = *output = (TwoInterfFindChirpFilterParams *)
    LALCalloc( 1, sizeof(TwoInterfFindChirpFilterParams) );
  if ( ! outputPtr )
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  /* create memory for twointerffindchirp filter parameters */
  outputPtr->filterParams1 = (FindChirpFilterParams *)
    LALCalloc( 1, sizeof(FindChirpFilterParams) );
  if ( !outputPtr->filterParams1 )
    {
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  outputPtr->filterParams2 = (FindChirpFilterParams *)
    LALCalloc( 1, sizeof(FindChirpFilterParams) );
  if ( !outputPtr->filterParams2 )
    {
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }

  /* create memory for the chisq parameters */
  outputPtr->filterParams1->chisqParams = (FindChirpChisqParams *)
    LALCalloc( 1, sizeof(FindChirpChisqParams) );
  if ( !outputPtr->filterParams1->chisqParams )
    {
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }

  outputPtr->filterParams2->chisqParams = (FindChirpChisqParams *)
    LALCalloc( 1, sizeof(FindChirpChisqParams) );
  if ( !outputPtr->filterParams2->chisqParams )
    {
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  /* create memory for the chisq input */
  outputPtr->filterParams1->chisqInput = (FindChirpChisqInput *)
    LALCalloc( 1, sizeof(FindChirpChisqInput) );
  if ( !outputPtr->filterParams1->chisqInput )
  {
    LALFree( outputPtr->filterParams1->chisqParams );
    LALFree( outputPtr->filterParams1 );
    LALFree( outputPtr );
    *output = NULL;
    ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
  }

  outputPtr->filterParams2->chisqInput = (FindChirpChisqInput *)
    LALCalloc( 1, sizeof(FindChirpChisqInput) );
  if ( !outputPtr->filterParams2->chisqInput )
    {
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  
  /*
   *
   * create fft plans and workspace vectors
   *
   */
  
  
  /* create plans for optimal filter */
  LALCreateReverseComplexFFTPlan( status->statusPtr, 
				  &(outputPtr->filterParams1->invPlan), 
				  params->initParams1->numPoints, 0); 
  
  BEGINFAIL( status )
    {
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;
    }
  ENDFAIL( status );
  
  LALCreateReverseComplexFFTPlan( status->statusPtr, 
				  &(outputPtr->filterParams2->invPlan), 
				  params->initParams2->numPoints, 0); 
  
  BEGINFAIL( status )
    {
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;
    }
  ENDFAIL( status );
  
  /* create workspace vectors for optimal filter: time domain */
  LALCCreateVector( status->statusPtr, &(outputPtr->filterParams1->qVec), 
		    params->initParams1->numPoints );
  BEGINFAIL( status )
    {
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams2->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;
    }
  ENDFAIL( status );
  
  LALCCreateVector( status->statusPtr, &(outputPtr->filterParams2->qVec), 
		    params->initParams2->numPoints );
  BEGINFAIL( status )
    {
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams1->qVec) ), 
	   status );
      
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams2->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;
    }
  ENDFAIL( status );

  /* create workspace vectors for optimal filter: freq domain */
  LALCCreateVector( status->statusPtr, &(outputPtr->filterParams1->qtildeVec), 
		    params->initParams1->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, 
			    &(outputPtr->filterParams2->qVec) ), 
	 status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams1->qVec) ), 
	   status );
      
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams2->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;    

  }
  ENDFAIL( status );

  LALCCreateVector( status->statusPtr, &(outputPtr->filterParams2->qtildeVec), 
		    params->initParams2->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, 
			    &(outputPtr->filterParams1->qtildeVec) ), 
	 status );
    TRY( LALCDestroyVector( status->statusPtr, 
			    &(outputPtr->filterParams2->qVec) ), 
	 status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams1->qVec) ), 
	   status );
      
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams2->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;    
  }
  ENDFAIL( status );
  
  /* create workspace vector for chisq filter */
  LALCreateVector (status->statusPtr, &(outputPtr->filterParams1->chisqVec), 
      params->initParams1->numPoints);
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, 
			    &(outputPtr->filterParams2->qtildeVec) ), 
	 status );
    TRY( LALCDestroyVector( status->statusPtr, 
			    &(outputPtr->filterParams1->qtildeVec) ), 
	 status );
    TRY( LALCDestroyVector( status->statusPtr, 
			    &(outputPtr->filterParams2->qVec) ), 
	 status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams1->qVec) ), 
	   status );
      
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams2->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;    
  }
  ENDFAIL( status );

  LALCreateVector (status->statusPtr, &(outputPtr->filterParams2->chisqVec), 
		   params->initParams2->numPoints);
  BEGINFAIL( status )
    {
      TRY( LALDestroyVector( status->statusPtr, 
			     &(outputPtr->filterParams1->chisqVec) ), 
	   status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams2->qtildeVec) ), 
	   status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams1->qtildeVec) ), 
	   status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams2->qVec) ), 
	   status );
      TRY( LALCDestroyVector( status->statusPtr, 
			      &(outputPtr->filterParams1->qVec) ), 
	   status );
      
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams1->invPlan) ), 
	   status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
				     &(outputPtr->filterParams2->invPlan) ), 
	   status );
      
      LALFree( outputPtr->filterParams1->chisqInput );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqInput );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr->filterParams1->chisqParams );
      LALFree( outputPtr->filterParams1 );
      LALFree( outputPtr->filterParams2->chisqParams );
      LALFree( outputPtr->filterParams2 );
      LALFree( outputPtr );
      *output = NULL;    
    }
  ENDFAIL( status );
  
  
  /*
   *
   * create vectors to store snrsq, if required
   *
   */
  
  
  if ( params->initParams1->createRhosqVec )
    {
      LALCreateVector (status->statusPtr, 
		       &(outputPtr->filterParams1->rhosqVec), 
		       params->initParams1->numPoints);
      BEGINFAIL( status )
	{
	  TRY( LALDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->chisqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->chisqVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->qVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->qVec) ), 
	       status );
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(outputPtr->filterParams1->invPlan)),
	       status );
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(outputPtr->filterParams2->invPlan)),
	       status );
	  
	  LALFree( outputPtr->filterParams1->chisqInput );
	  LALFree( outputPtr->filterParams1 );
	  LALFree( outputPtr->filterParams2->chisqInput );
	  LALFree( outputPtr->filterParams2 );
	  LALFree( outputPtr->filterParams1->chisqParams );
	  LALFree( outputPtr->filterParams1 );
	  LALFree( outputPtr->filterParams2->chisqParams );
	  LALFree( outputPtr->filterParams2 );
	  LALFree( outputPtr );
	  *output = NULL;    
	}
      ENDFAIL( status );
    }    
  
  if ( params->initParams2->createRhosqVec )
    {
      LALCreateVector (status->statusPtr, 
		       &(outputPtr->filterParams2->rhosqVec), 
		       params->initParams2->numPoints);
      BEGINFAIL( status )
	{
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams1->rhosqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams2->chisqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams2->chisqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams1->chisqVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->qVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->qVec) ), 
	       status );
	  
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(outputPtr->filterParams1->invPlan)),
	       status );
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(outputPtr->filterParams2->invPlan)),
	       status );
	  
	  LALFree( outputPtr->filterParams1->chisqInput );
	  LALFree( outputPtr->filterParams1 );
	  LALFree( outputPtr->filterParams2->chisqInput );
	  LALFree( outputPtr->filterParams2 );
	  LALFree( outputPtr->filterParams1->chisqParams );
	  LALFree( outputPtr->filterParams1 );
	  LALFree( outputPtr->filterParams2->chisqParams );
	  LALFree( outputPtr->filterParams2 );
	  LALFree( outputPtr );
	  *output = NULL;    
	}
      ENDFAIL( status );
    }    
  
  if ( params->createTwoInterfRhosqVec )
    {
      LALCreateVector (status->statusPtr, 
		       &(outputPtr->twoInterfRhosqVec), 
		       params->numPoints);
      BEGINFAIL( status )
	{
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams2->rhosqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams1->rhosqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams2->chisqVec) ), 
	       status );
	  TRY( LALDestroyVector( status->statusPtr, 
				 &(outputPtr->filterParams1->chisqVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams2->qVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, 
				  &(outputPtr->filterParams1->qVec) ), 
	       status );
	  
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(outputPtr->filterParams1->invPlan)),
	       status );
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(outputPtr->filterParams2->invPlan)),
	       status );
	  
	  LALFree( outputPtr->filterParams1->chisqInput );
	  LALFree( outputPtr->filterParams1 );
	  LALFree( outputPtr->filterParams2->chisqInput );
	  LALFree( outputPtr->filterParams2 );
	  LALFree( outputPtr->filterParams1->chisqParams );
	  LALFree( outputPtr->filterParams1 );
	  LALFree( outputPtr->filterParams2->chisqParams );
	  LALFree( outputPtr->filterParams2 );
	  LALFree( outputPtr );
	  *output = NULL;    
	}
      ENDFAIL( status );
    }    
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALTwoInterfFindChirpFilterFinalize (
				     LALStatus                        *status,
				     TwoInterfFindChirpFilterParams  **output
				     )
{
  TwoInterfFindChirpFilterParams        *outputPtr;
  
  INITSTATUS( status, "LALTwoInterfFindChirpFilterFinalize", 
	      TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  
  /*
   *
   * check that the arguments are reasonable
   *
   */

  
  /* check that handle is non-null and points to a non-null pointer */
  ASSERT( output, status, TWOINTERFFINDCHIRPH_ENULL, 
	  TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *output, status, TWOINTERFFINDCHIRPH_ENULL, 
	  TWOINTERFFINDCHIRPH_MSGENULL );
  

  /*
   *
   * destroy the filter parameter structure
   *
   */
  
  
  /* local pointer to output structure */
  outputPtr = *output;
  
  
  /*
   *
   * destroy fft plans and workspace vectors
   *
   */
  
  
  /* destroy plans for optimal filter */
  LALDestroyComplexFFTPlan( status->statusPtr, 
			    &(outputPtr->filterParams1->invPlan) );
  CHECKSTATUSPTR( status );
  
  LALDestroyComplexFFTPlan( status->statusPtr, 
			    &(outputPtr->filterParams2->invPlan) );
  CHECKSTATUSPTR( status );
  
  /* destroy workspace vectors for optimal filter: freq domain */
  LALCDestroyVector( status->statusPtr, 
		     &(outputPtr->filterParams1->qVec) );
  CHECKSTATUSPTR( status );
  
  LALCDestroyVector( status->statusPtr, &(outputPtr->filterParams2->qVec) );
  CHECKSTATUSPTR( status );
  
  /* destroy workspace vectors for optimal filter: freq domain */
  LALCDestroyVector( status->statusPtr, 
		     &(outputPtr->filterParams1->qtildeVec) );
  CHECKSTATUSPTR( status );
  
  LALCDestroyVector( status->statusPtr, 
		     &(outputPtr->filterParams2->qtildeVec) );
  CHECKSTATUSPTR( status );
  
  /* destroy workspace vectors for chisq filter */
  LALDestroyVector( status->statusPtr, &(outputPtr->filterParams1->chisqVec) );
  CHECKSTATUSPTR( status );
  
  LALDestroyVector( status->statusPtr, &(outputPtr->filterParams2->chisqVec) );
  CHECKSTATUSPTR( status );
  
  /*
   *
   * free the chisq structures
   *
   */
  
  
  /* parameter structure */
  LALFree( outputPtr->filterParams1->chisqParams );
  LALFree( outputPtr->filterParams1 );
  LALFree( outputPtr->filterParams2->chisqParams );
  LALFree( outputPtr->filterParams2 );
  
  /* input structure */  
  LALFree( outputPtr->filterParams1->chisqInput );
  LALFree( outputPtr->filterParams1 );
  LALFree( outputPtr->filterParams2->chisqInput );
  LALFree( outputPtr->filterParams2 );
  
  
  /*
   *
   * destroy vectors to store snrsq, if they exist
   *
   */
  
  
  if ( outputPtr->twoInterfRhosqVec )
    {
      LALDestroyVector( status->statusPtr, &(outputPtr->twoInterfRhosqVec) );
      CHECKSTATUSPTR( status );
    }    
  
  if ( outputPtr->filterParams1->rhosqVec )
    {
      LALDestroyVector( status->statusPtr, 
			&(outputPtr->filterParams1->rhosqVec) );
      CHECKSTATUSPTR( status );
    }    
  
  if ( outputPtr->filterParams2->rhosqVec )
    {
      LALDestroyVector( status->statusPtr, 
			&(outputPtr->filterParams2->rhosqVec) );
      CHECKSTATUSPTR( status );
    }    
  
  
  /*
   *
   * free memory for the TwoInterfFindChirpFilterParams
   *
   */
  
  
  LALFree( outputPtr );
  *output = NULL;
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALTwoInterfFindChirpFilterSegment (
				    LALStatus                       *status,
				    TwoInterfInspiralEvent         **eventList,
				    TwoInterfFindChirpFilterInput   *input,
				    TwoInterfFindChirpFilterParams  *params
				    )
{
  UINT4                          j, k, m;
  UINT4                          deltaEventIndex;
  UINT4                          twoInterfEventId = 0; 
  REAL4                          deltaT; 
  REAL4                          maxLag; /*lightTravel time between dets 1,2*/
  REAL4                          s[3];
  REAL4                          distance;
  REAL4                          norm1;
  REAL4                          norm2;
  REAL4                          modqsq1;
  REAL4                          modqsq2;
  REAL4                          rhosq = 0;
  REAL4                          twoInterfRhosqThresh = 0;
  BOOLEAN                        haveChisq    = 0;
  UINT4                          numPoints; 
  COMPLEX8                      *q1tilde      = NULL;
  COMPLEX8                      *q1           = NULL;
  COMPLEX8                      *inputData1   = NULL;
  COMPLEX8                      *q2tilde      = NULL;
  COMPLEX8                      *q2           = NULL;
  COMPLEX8                      *inputData2   = NULL;
  COMPLEX8                      *tmpltSignal  = NULL;
  TwoInterfInspiralEvent        *thisEvent    = NULL; 
  LALDetectorPair               *detectors = NULL;
  
  
  INITSTATUS( status, "LALTwoInterFindChirpFilter", 
	      TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  
  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, 
	  TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );
  
  /* check that the parameter structure exists */
  ASSERT( params, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the parameter sub-structures exist */
  ASSERT( params->filterParams1, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the filter parameters are reasonable */
  ASSERT( params->maxLag > 0, status,
	  TWOINTERFFINDCHIRPH_EMLZO, TWOINTERFFINDCHIRPH_MSGEMLZO );
  ASSERT( params->twoInterfRhosqThresh > 0, status,
	  TWOINTERFFINDCHIRPH_ERHOT, TWOINTERFFINDCHIRPH_MSGERHOT );
  ASSERT( params->filterParams1->deltaT > 0, status,
	  TWOINTERFFINDCHIRPH_EDTZO, TWOINTERFFINDCHIRPH_MSGEDTZO );
  ASSERT( params->filterParams2->deltaT > 0, status,
	  TWOINTERFFINDCHIRPH_EDTZO, TWOINTERFFINDCHIRPH_MSGEDTZO );
  ASSERT( params->filterParams1->chisqThresh > 0, status,
	  TWOINTERFFINDCHIRPH_ECHIT, TWOINTERFFINDCHIRPH_MSGECHIT ); 
  ASSERT( params->filterParams2->chisqThresh > 0, status,
	  TWOINTERFFINDCHIRPH_ECHIT, TWOINTERFFINDCHIRPH_MSGECHIT ); 
  
  

  /* check that the fft plan exists */
  ASSERT( params->filterParams1->invPlan, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2->invPlan, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the workspace vectors exist for detector 1 */
  ASSERT( params->filterParams1->qVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams1->qVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams1->qtildeVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams1->qtildeVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the workspace vectors exist for detector 2 */
  ASSERT( params->filterParams2->qVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2->qVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2->qtildeVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2->qtildeVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the chisq parameter and input structures exist */
  ASSERT( params->filterParams1->chisqParams, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams1->chisqInput, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2->chisqParams, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->filterParams2->chisqInput, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->twoInterfRhosqVec ) 
    {
      ASSERT( params->twoInterfRhosqVec->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
  if ( params->filterParams1->rhosqVec ) 
    {
      ASSERT( params->filterParams1->rhosqVec->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
  if ( params->filterParams2->rhosqVec ) 
    {
      ASSERT( params->filterParams2->rhosqVec->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
 

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->filterParams1->chisqVec ) 
    {
      ASSERT( params->filterParams1->chisqVec->data, status,
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
  if ( params->filterParams2->chisqVec ) 
    {
      ASSERT( params->filterParams2->chisqVec->data, status,
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
  
  /* check that the input structure exists */
  ASSERT( input, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the input structure contains some input */
  ASSERT( input->tmplt, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->fcTmplt, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->detectors, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->segment1, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->segment2, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL ); /* introduced this line */
  
  /*
   *
   * point local pointers to input and output pointers
   *
   */

  
  /* workspace vectors for detector 1 */
  q1 = params->filterParams1->qVec->data;
  q1tilde = params->filterParams1->qtildeVec->data;
  
  /* workspace vectors for detector 2 */
  q2 = params->filterParams2->qVec->data;
  q2tilde = params->filterParams2->qtildeVec->data;
  
  /* template and data */
  inputData1 = input->segment1->data->data->data; 
  inputData2 = input->segment2->data->data->data;
  tmpltSignal = input->fcTmplt->data->data;
  deltaT = params->filterParams1->deltaT; 
 
  /* calculate separation vector between sites */
  for ( j=0; j<3; j++) {
    s[j] = (REAL4) (detectors->detectorOne.location[j] -
		    detectors->detectorTwo.location[j]);
  }
  /* calculate distance between sites (in meters) */
  distance = sqrt( cartesianInnerProduct(s,s) );
  
  /* light-travel time between sites (in seconds) */
  maxLag = distance / LAL_C_SI;
  
  /* number of points in a segment from either detector 1 or detector 2 */
  numPoints = params->filterParams1->qVec->length;
  

  /*
   *
   * calculate the minimum distance between distinct events
   *
   */
  

  /* This should be the length of a chirp, as in the FindChirp package, */
  /* As in FindChirpFilter.c, it is hard wired to sixty seconds at the moment*/
  deltaEventIndex = (UINT4) rint( (60.0 / deltaT) + 1.0 );
  
  
  
  /*
   *
   * compute q1tilde and q1
   *
   */
  
  
  memset( q1tilde, 0, numPoints * sizeof(COMPLEX8) );
  
  /* q1tilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
    {
      REAL4 r1 = inputData1[k].re;
      REAL4 s1 = inputData1[k].im;
      REAL4 x = tmpltSignal[k].re;
      REAL4 y = 0 - tmpltSignal[k].im;       /* note complex conjugate */
      
      q1tilde[k].re = r1*x - s1*y;
      q1tilde[k].im = r1*y + s1*x;
    }
  
  /* q1tilde negative frequency only: not DC or nyquist */
  if ( params->filterParams1->computeNegFreq )
    {
      for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
	{
	  REAL4 r1 = inputData1[k].re;
	  REAL4 s1 = inputData1[k].im;
	  REAL4 x = tmpltSignal[k].re;
	  REAL4 y = 0 - tmpltSignal[k].im;     /* note complex conjugate */
	  
	  q1tilde[k].re = r1*x - s1*y;
	  q1tilde[k].im = r1*y + s1*x;
	}
    }
  
  /* inverse fft to get q1 */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->filterParams1->qVec, 
			params->filterParams1->qtildeVec, 
			params->filterParams1->invPlan );
  CHECKSTATUSPTR( status );


  /*
   *
   * compute q2tilde and q2
   *
   */
  

  memset( q2tilde, 0, numPoints * sizeof(COMPLEX8) );
  
  /* q2tilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
    {
      REAL4 r2 = inputData2[k].re;
      REAL4 s2 = inputData2[k].im;
      REAL4 x = tmpltSignal[k].re;
      REAL4 y = 0 - tmpltSignal[k].im;       /* note complex conjugate */
      
      q2tilde[k].re = r2*x - s2*y;
      q2tilde[k].im = r2*y + s2*x;
    }
  
  /* q2tilde negative frequency only: not DC or nyquist */
  if ( params->filterParams1->computeNegFreq )
    {
      for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
	{
	  REAL4 r2 = inputData2[k].re;
	  REAL4 s2 = inputData2[k].im;
	  REAL4 x = tmpltSignal[k].re;
	  REAL4 y = 0 - tmpltSignal[k].im;     /* note complex conjugate */
	  
	  q2tilde[k].re = r2*x - s2*y;
	  q2tilde[k].im = r2*y + s2*x;
	}
    }
  
  /* inverse fft to get q2 */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->filterParams2->qVec, 
			params->filterParams2->qtildeVec, 
			params->filterParams2->invPlan );
  CHECKSTATUSPTR( status );

  
  /* 
   *
   * calculate signal to noise squared 
   *
   */
  
  
  /* if full snrsq vector for the network is required, set it to zero */
  if ( params->twoInterfRhosqVec )
    memset( params->twoInterfRhosqVec->data, 0, numPoints * sizeof( REAL4 ) );
  
  /* if full snrsq vector for detector 1 is required, set it to zero */
  if ( params->filterParams1->rhosqVec )
    memset( params->filterParams1->rhosqVec->data, 0, 
	    numPoints * sizeof( REAL4 ) );
  
  /* normalization */
  norm1 = 4.0 * (deltaT / (REAL4)numPoints) / input->segment1->segNorm;
  norm2 = 4.0 * (deltaT / (REAL4)numPoints) / input->segment2->segNorm;
  
  m = ( maxLag / deltaT ); 
  
  for ( j = 0; j < numPoints; ++j ) 
    {
      
      for ( k = ( j < m ) ? 0 : ( j - m + 1) ; k <= j + m - 1 ; ++k )
	{
	  modqsq1 = ( q1[j].re * q1[j].re + q1[j].im * q1[j].im );
	  modqsq2 = ( q2[k].re * q2[k].re + q2[k].im * q2[k].im );
	  
	  rhosq = norm1 * modqsq1 + norm2 * modqsq2; 
	  
	  /* if full snrsq vector is required, store the snrsq */
	  if ( params->twoInterfRhosqVec ) 
	    params->twoInterfRhosqVec->data[j] = rhosq;
	  
	  /* if snrsq exceeds threshold at any point */
	  if ( rhosq > twoInterfRhosqThresh )
	    {
	      /* compute chisq vector for detector 1 if it does not exist */
	      if ( ! haveChisq )
		{
		  memset( params->filterParams1->chisqVec->data, 0, 
			  params->filterParams1->chisqVec->length 
			  * sizeof(REAL4) );

		  /* pointers to chisq input for detector 1 */
		  params->filterParams1->chisqInput->qtildeVec 
		    = params->filterParams1->qtildeVec;
		  params->filterParams1->chisqInput->qVec      
		    = params->filterParams1->qVec;
		  
		  /* pointer to the chisq bin vector in segment of detector1 */
		  params->filterParams1->chisqParams->chisqBinVec 
		    = input->segment1->chisqBinVec;
		  params->filterParams1->chisqParams->chisqNorm   
		    = sqrt( norm1 );
		  
		  /* compute the chisq threshold for detector 1 */
		  LALFindChirpChisqVeto( status->statusPtr, 
					 params->filterParams1->chisqVec, 
					 params->filterParams1->chisqInput, 
					 params->filterParams1->chisqParams );
		  CHECKSTATUSPTR (status); 
		  
		  /* compute chisqVector for detector 2 if it does not exist */
		  
		  memset( params->filterParams2->chisqVec->data, 0, 
			  params->filterParams2->chisqVec->length 
			  * sizeof(REAL4) );
		  
		  /* pointers to chisq input for detector 2 */
		  params->filterParams2->chisqInput->qtildeVec 
		    = params->filterParams2->qtildeVec;
		  params->filterParams2->chisqInput->qVec      
		    = params->filterParams2->qVec;
		  
		  /* pointer to the chisq bin vector in segment of detector2 */
		  params->filterParams2->chisqParams->chisqBinVec 
		    = input->segment2->chisqBinVec;
		  params->filterParams2->chisqParams->chisqNorm   
		    = sqrt( norm2 );
		  
		  /* compute the chisq threshold for detector 2 */
		  LALFindChirpChisqVeto( status->statusPtr, 
					 params->filterParams2->chisqVec, 
					 params->filterParams2->chisqInput, 
					 params->filterParams2->chisqParams );
		  CHECKSTATUSPTR (status); 
		}
	      
	      /* if chisq for either detector drops below threshold,
		 then start processing events */
	      if ( ( params->filterParams1->chisqVec->data[j] 
		     < params->filterParams1->chisqThresh ) 
		   || ( params->filterParams2->chisqVec->data[j] 
			< params->filterParams2->chisqThresh ) )
		{
		  if ( ! *eventList )
		    {
		      /* if this is the first event, start the list */
		      thisEvent = *eventList = (TwoInterfInspiralEvent *) 
			LALCalloc( 1, sizeof(TwoInterfInspiralEvent) );
		      if ( ! thisEvent )
			{
			  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, 
				 TWOINTERFFINDCHIRPH_MSGEALOC );
			}
		      thisEvent->twoInterfId = twoInterfEventId++;
		      
		      /* stick minimal data into the event */
		      thisEvent->eventIn1->timeIndex = j;
		      thisEvent->eventIn2->timeIndex = k;
		      thisEvent->snrsq = rhosq;
		    }
		  else if ((j < thisEvent->eventIn1->timeIndex+deltaEventIndex 
			    || 
			    k <thisEvent->eventIn2->timeIndex+deltaEventIndex) 
			   && rhosq > thisEvent->snrsq )
		    {
		      /* if this is the same event, update the maximum */
		      thisEvent->eventIn1->timeIndex = j;
		      thisEvent->eventIn2->timeIndex = k;
		      thisEvent->snrsq = rhosq;
		    }
		  else if ( j > thisEvent->eventIn1->timeIndex+deltaEventIndex
			    || 
			    k >thisEvent->eventIn2->timeIndex+ deltaEventIndex)
		    {
		      /* clean up this event */
		      TwoInterfInspiralEvent *lastEvent;
		      INT8                    timeNS1;
		      INT8                    timeNS2;
		      
		      /* set the event LIGO GPS time of the event 
			 for detector 1 */
		      timeNS1 = 1000000000L * 
			(INT8) (input->segment1->data->epoch.gpsSeconds);
		      timeNS1 += 
			(INT8) (input->segment1->data->epoch.gpsNanoSeconds);
		      timeNS1 += 
			(INT8) (1e9 * (thisEvent->eventIn1->timeIndex)*deltaT);
		      thisEvent->eventIn1->time.gpsSeconds = 
			(INT4) (timeNS1/1000000000L);
		      thisEvent->eventIn1->time.gpsNanoSeconds = 
			(INT4) (timeNS1%1000000000L);
 
		      /* set the event LIGO GPS time of the event 
			 for detector 2 */
		      timeNS2 = 1000000000L * 
			(INT8) (input->segment2->data->epoch.gpsSeconds);
		      timeNS2 += 
			(INT8) (input->segment2->data->epoch.gpsNanoSeconds);
		      timeNS2 += 
			(INT8) (1e9 * (thisEvent->eventIn2->timeIndex)*deltaT);
		      thisEvent->eventIn2->time.gpsSeconds = 
			(INT4) (timeNS2/1000000000L);
		      thisEvent->eventIn2->time.gpsNanoSeconds = 
			(INT4) (timeNS2%1000000000L);
		      
		      /* copy the template into the event */
		      memcpy( &(thisEvent->tmplt), input->tmplt, 
			      sizeof(InspiralTemplate) );
		      thisEvent->tmplt.next = NULL;
		      thisEvent->tmplt.fine = NULL;
		      
		      /* set snrsq, chisq, sigma and effDist for this event */
		      thisEvent->eventIn1->chisq =  
			params->filterParams1->chisqVec->data[thisEvent->eventIn1->timeIndex];
		      thisEvent->eventIn2->chisq = 
			params->filterParams2->chisqVec->data[thisEvent->eventIn2->timeIndex];          
		      thisEvent->eventIn1->sigma   = norm1;
		      thisEvent->eventIn2->sigma   = norm2;
		      thisEvent->snrsq  = rhosq;
		      
		      /* allocate memory for the newEvent */
		      lastEvent = thisEvent;
		      
		      lastEvent->next = thisEvent =(TwoInterfInspiralEvent *) 
			LALCalloc( 1, sizeof(TwoInterfInspiralEvent) );
		      if ( ! lastEvent->next )
			{
			  ABORT( status, FINDCHIRPH_EALOC, 
				 FINDCHIRPH_MSGEALOC );
			}
		      thisEvent->twoInterfId = twoInterfEventId++;
		      
		      /* stick minimal data into the event */
		      thisEvent->eventIn1->timeIndex = j;
		      thisEvent->eventIn2->timeIndex = k;
		      thisEvent->snrsq = rhosq; 
		    }
		}
	    }
	}
    }    
  /* 
   *
   * clean up the last event if there is one
   *
   */
  
  
  if ( thisEvent )
    {
      INT8           timeNS1;
      INT8           timeNS2;
      
      /* set the event LIGO GPS time of the event in detector 1 */
      timeNS1 = 1000000000L * 
	(INT8) (input->segment1->data->epoch.gpsSeconds);
      timeNS1 += (INT8) (input->segment1->data->epoch.gpsNanoSeconds);
      timeNS1 += (INT8) (1e9 * (thisEvent->eventIn1->timeIndex) * deltaT);
      thisEvent->eventIn1->time.gpsSeconds = (INT4) (timeNS1/1000000000L);
      thisEvent->eventIn1->time.gpsNanoSeconds 
	= (INT4) (timeNS1%1000000000L);
      
      /* set the event LIGO GPS time of the event in detector 2 */
      timeNS2 = 1000000000L * 
	(INT8) (input->segment2->data->epoch.gpsSeconds);
      timeNS2 += (INT8) (input->segment2->data->epoch.gpsNanoSeconds);
      timeNS2 += (INT8) (1e9 * (thisEvent->eventIn2->timeIndex) * deltaT);
      thisEvent->eventIn2->time.gpsSeconds = (INT4) (timeNS2/1000000000L);
      thisEvent->eventIn2->time.gpsNanoSeconds 
	= (INT4) (timeNS2%1000000000L);
      
      /* copy the template into the event */
      memcpy( &(thisEvent->tmplt), input->tmplt, sizeof(InspiralTemplate));
      
      /* set snrsq, chisq, sigma and effDist for this event */
      thisEvent->eventIn1->chisq =  
	params->filterParams1->chisqVec->data[thisEvent->eventIn1->timeIndex];
      thisEvent->eventIn2->chisq = 
	params->filterParams2->chisqVec->data[thisEvent->eventIn2->timeIndex];          
      thisEvent->eventIn1->sigma   = norm1;
      thisEvent->eventIn2->sigma   = norm2;
      thisEvent->snrsq  = rhosq;
    }    
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
  
}
