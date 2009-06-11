/*
*  Copyright (C) 2007 Duncan Brown, Gareth Jones
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBCVSpinFilter.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVSpinFilterCV">
Author: Brown, D. A., Spinning BCV-Modifications: Jones, G.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinFilter.c}}
\label{ss:FindChirpBCVSpinFilter.c}

Provides functions to filter data for spinning BCV templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpBCVSpinFilterCP}
\idx{LALFindChirpBCVSpinFilter()}

The function \texttt{LALFindChirpBCVSpinFilterSegment()} filters data for
spinning BCV templates as described by the algorithm below.

\subsubsection*{Algorithm}

Using the quantities calculated in {\tt LALFindChirpBCVSpinTemplate()} we
construct a template and filter our data producing a clustered
signal-to-noise ratio time series.
We filter our data in 256 second data segments.
We first calculate the following functions in the frequency domain:
\begin{eqnarray}
{\tt qtilde}         & = & \frac {\mathcal{\widehat{A}}_1(f)
e^{i \psi_{NM} (f)} s^* (f)} {S_h (f)} \nonumber\\
{\tt qtildeBCVSpin1} & = & \frac {\mathcal{\widehat{A}}_2(f)
e^{i \psi_{NM} (f)} s^* (f)} {S_h (f)} \nonumber\\
{\tt qtildeBCVSpin2} & = & \frac {\mathcal{\widehat{A}}_3(f)
e^{i \psi_{NM} (f)} s^* (f)} {S_h (f)}
\end{eqnarray}
where $\mathcal{\widehat{A}}_1(f)$, $\mathcal{\widehat{A}}_2(f)$
and $\mathcal{\widehat{A}}_3(f)$ are the orthonormal amplitude functions
and $\psi_{NM} (f)$ is the non-modulational phase of our template. These
quantitites were calculated in {\tt LALFindChirpBCVSpinTemplate()}. $s^*$
is the complex conjugate of our (detector) data in the frequency domain
and $S_h (f)$ is our estimate of the power spectral density of the detector
data estimated over a 2048 second ``blocks''.
Using inverse FFTs we construct the complex time domain quantities
{\tt q}, {\tt qBCVSpin1} and {\tt qBCVSpin2}.
We then calculate signal-to-noise ratio as
\begin{eqnarray}
\rho(t)^2 & = & {\tt q.re}^2
           + {\tt q.im}^2
           + {\tt qBCVSpin1.re}^2
           + {\tt qBCVSpin1.im}^2 + \nonumber\\
       &   & {\tt qBCVSpin2.re}^2
           + {\tt qBCVSpin2.im}^2.
\end{eqnarray}
We then look for values of $\rho(t)$ above our threshold - note that the
$\beta = 0$ threshold is currently hardcoded. We do not calculate
signal-to-noise ratio for the 64 second stretch at the beginning and end
of each data segment to avoid edge-effects. These times are picked up by
overlapping our 256 second data segments.
For times for which signal-to-noise ratio is calculated we have the option
of clustering our output using the {\tt --cluster-method window} option in
{\tt lalapps\_inspiral} with an appropriate choice of cluster length.
For events that pass the signal-to-noise ratio threshold and survive
clustering we store the template parameters $\psi_0$, $\psi_3$, $\beta$
and $f_{final}$ as well as 6 $\alpha$ values which encode the relative
contribution of the {\tt q}, {\tt qBCVSpin1} and {\tt qBCVSpin2} functions
to the overall signal-to-noise ratio. These are simply calculated as
\begin{eqnarray}
\alpha_1 & = & {\tt q.re} / \rho \nonumber \\
\alpha_2 & = & {\tt qBCVSpin1.re} / \rho \nonumber \\
\alpha_3 & = & {\tt qBCVSpin2.re} / \rho \nonumber \\
\alpha_4 & = & {\tt q.im} / \rho \nonumber \\
\alpha_5 & = & {\tt qBCVSpin1.im} / \rho \nonumber \\
\alpha_6 & = & {\tt qBCVSpin2.im} / \rho.
\end{eqnarray}
These obey $\sum_{i=1}^6 \alpha_i = 1$ and might prove useful in future
signal based vetoe studies.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpBCVSpinFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCVSpin.h>

#define rint(x) (floor((x)+0.5))

NRCSID (FINDCHIRPBCVSPINFILTERC, "$Id$");

/* <lalVerbatim file="FindChirpBCVSpinFilterCP"> */
void
LALFindChirpBCVSpinFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    FindChirpDataParams        *fcDataParams
  )
/* </lalVerbatim> */
{
  UINT4                 j, k;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT;
  REAL4                 deltaF;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  UINT4                 eventStartIdx    = 0;
  REAL4                 chirpTime        = 0;
  COMPLEX8             *qtilde           = NULL;
  COMPLEX8             *qtildeBCVSpin1   = NULL;
  COMPLEX8             *qtildeBCVSpin2   = NULL;
  COMPLEX8             *q                = NULL;
  COMPLEX8             *qBCVSpin1        = NULL;
  COMPLEX8             *qBCVSpin2        = NULL;
  COMPLEX8             *tmpltSignal      = NULL;
  SnglInspiralTable    *thisEvent        = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  COMPLEX8             *wtilde;
  COMPLEX8             *inputData1;
  REAL4 		rho = 0.0;
  REAL4                 invRho = 0.0;
  REAL4 		m;
  REAL4                 normFac;
  REAL4                 normFacSq;
  REAL4   		alpha1hat;
  REAL4                 alpha2hat;
  REAL4                 alpha3hat;
  REAL4                 alpha4hat;
  REAL4                 alpha5hat;
  REAL4                 alpha6hat;
  REAL8                *A1Vec            = NULL;
  REAL8                *A2Vec            = NULL;
  REAL8                *A3Vec            = NULL;
  REAL4                 normData;
  REAL4                 invRootNormData;

  REAL8                 I;
  REAL8                 J;
  REAL8                 K;

  REAL8			rootI;
  REAL8                 denominator;
  REAL8                 rootDenominator;
  REAL8                 numerator1;
  REAL8                 denominator1;

/* REAL4                 A1A1             = 0.0;
 REAL4                 A2A2             = 0.0;
 REAL4                 A3A3             = 0.0;
 REAL4                 A1A2             = 0.0;
 REAL4                 A1A3             = 0.0;
 REAL4                 A2A3             = 0.0;*/



  INITSTATUS( status, "LALFindChirpBCVSpinFilter", FINDCHIRPBCVSPINFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   * may need to remove asserts regarding chisq
   *
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh > 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
 /* ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL,
	  FINDCHIRPH_MSGENULL);*/

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if we are doing a chisq, check we can store the data */
  if ( input->segment->chisqBinVec->length )
  {
    ASSERT( params->chisqVec, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the template and the segment are both BCVSpin */
  ASSERT( input->fcTmplt->tmplt.approximant == BCVSpin, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCVSpin, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );

  /* template and data */
  inputData1    = input->segment->data->data->data; /* data */
  tmpltSignal   = input->fcTmplt->data->data;       /* expPsi */
  wtilde        = fcDataParams->wtildeVec->data;    /* inverse psd */

  numPoints     = params->qVec->length;
  normFac       = 4./numPoints;
  normFacSq     = pow(normFac, 2);
  deltaT        = params->deltaT;
  deltaF        =  1.0/((REAL4)numPoints * deltaT);

  /* amplitude vectors calculated in LALFindChirpBCVSpinTemplate() */
  A1Vec = input->fcTmplt->A1BCVSpin->data;
  A2Vec = input->fcTmplt->A2BCVSpin->data;
  A3Vec = input->fcTmplt->A3BCVSpin->data;

  I           = input->fcTmplt->momentI;
  J           = input->fcTmplt->momentJ;
  K           = input->fcTmplt->momentK;

  rootI           = input->fcTmplt->rootMomentI;
  denominator     = input->fcTmplt->numFactor;
  rootDenominator = input->fcTmplt->numFactor1;
  numerator1      = input->fcTmplt->numFactor2;
  denominator1    = input->fcTmplt->numFactor3;

  /* workspace vectors */
  q              = params->qVec->data;
  qBCVSpin1      = params->qVecBCVSpin1->data;
  qBCVSpin2      = params->qVecBCVSpin2->data;

  qtilde         = params->qtildeVec->data;
  qtildeBCVSpin1 = params->qtildeVecBCVSpin1->data;
  qtildeBCVSpin2 = params->qtildeVecBCVSpin2->data;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;

  /* finding cross product of data with itself,
     to be used for normalisation later  */

 /* checking orthonormalisation of A vectors */


  /*  {
           fprintf (stdout, "Checking orthonormalisation of amplitude vectors \n");

           for (k=1; k < numPoints/2; ++k)
           {

	            A1A1 += A1Vec[k] * A1Vec[k] * wtilde[k].re;
	            A2A2 += A2Vec[k] * A2Vec[k] * wtilde[k].re;
	            A3A3 += A3Vec[k] * A3Vec[k] * wtilde[k].re;
	            A1A2 += A1Vec[k] * A2Vec[k] * wtilde[k].re;
	            A1A3 += A1Vec[k] * A3Vec[k] * wtilde[k].re;
	            A2A3 += A2Vec[k] * A3Vec[k] * wtilde[k].re;
	   }

	   A1A1 *= 4 * deltaF;
	   A2A2 *= 4 * deltaF;
	   A3A3 *= 4 * deltaF;
	   A1A2 *= 4 * deltaF;
	   A1A3 *= 4 * deltaF;
	   A2A3 *= 4 * deltaF;

	  fprintf (stdout, "A1hat cross A1hat %e\n", A1A1);
	  fprintf (stdout, "A2hat cross A2hat %e\n", A2A2);
	  fprintf (stdout, "A3hat cross A3hat %e\n", A3A3);
	  fprintf (stdout, "A1hat cross A2hat %e\n", A1A2);
	  fprintf (stdout, "A1hat cross A3hat %e\n", A1A3);
	  fprintf (stdout, "A2hat cross A3hat %e\n\n", A2A3);
  } */



  	normData = 0.;

  	for (k = 0; k < (numPoints/2)+1; ++k )
  	{
  		normData += ((inputData1[k].re * inputData1[k].re)
 		 	+ (inputData1[k].im * inputData1[k].im))
             	 	* wtilde[k].re;
  	}

  	normData *= deltaT * normFac;

  	/*fprintf (stdout, "Cross product of input data with itself = %e\n",
			normData);*/

  	invRootNormData = pow(normData,-0.5);


        /*UNCOMMENT LOOP BELOW TO NORMALISE INPUT  */
    	/*fprintf (stdout,
		"Normalising input data (should not generally be used) \n");

  	for (k = 0; k < (numPoints/2)+1; ++k )
  	{
  		inputData1[k].re *= invRootNormData;
		inputData1[k].im *= invRootNormData;
  	}*/
        /*UNCOMMENT LOOP ABOVE TO NORMALISE INPUT  */

  memset( qtilde,         0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin1, 0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin2, 0, numPoints * sizeof(COMPLEX8) );

  /*
   *
   * Compute qtilde, qtildeBCVSpin1, qtildeBCVSpin2
   *
   */


  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
  	REAL4 r        =   inputData1[k].re;
    	REAL4 s        =   inputData1[k].im;

    	REAL4 x =  tmpltSignal[k].re;
    	REAL4 y =  0. - tmpltSignal[k].im;

    	qtilde[k].re        = r * x - s * y ;
    	qtilde[k].im        = s * x + r * y ;

/*      	qtilde[k].re *= wtilde[k].re;
      	qtilde[k].im *= wtilde[k].re; */

    	qtildeBCVSpin1[k] =  qtilde[k];
    	qtildeBCVSpin2[k] =  qtilde[k];

    	/* real parts */

     	qtilde[k].re         *= A1Vec[k];
     	qtildeBCVSpin1[k].re *= A2Vec[k];
     	qtildeBCVSpin2[k].re *= A3Vec[k];

    	/* imaginary parts */

     	qtilde[k].im         *= A1Vec[k];
     	qtildeBCVSpin1[k].im *= A2Vec[k];
     	qtildeBCVSpin2[k].im *= A3Vec[k];
  }


  /*
   *
   * inverse fft to get q, qBCVSpin1 and qBCVSpin2
   *
   */

  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec,
  	params->qtildeVec, params->invPlan );
  CHECKSTATUSPTR( status );

  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCVSpin1,
	params->qtildeVecBCVSpin1, params->invPlan );
  CHECKSTATUSPTR( status );

  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCVSpin2,
	params->qtildeVecBCVSpin2, params->invPlan );
  CHECKSTATUSPTR( status );

  /*
   *
   * calculate signal to noise squared
   *
   */

  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
  	memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec )
  {
   	/* fprintf (stdout, "Inside rhosq loop \n "); */

    	memcpy( params->rhosqVec->name, input->segment->data->name,
        	LALNameLength * sizeof(CHAR) );
    	memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch),
        	sizeof(LIGOTimeGPS) );
    	params->rhosqVec->deltaT = input->segment->deltaT;

	for ( j = 0; j < numPoints; ++j)
	{

		REAL4 	rhoSq =
			( ( q[j].re * q[j].re + q[j].im * q[j].im ) +
   	                ( qBCVSpin1[j].re * qBCVSpin1[j].re
			+ qBCVSpin1[j].im * qBCVSpin1[j].im ) +
   	                ( qBCVSpin2[j].re * qBCVSpin2[j].re
			+ qBCVSpin2[j].im * qBCVSpin2[j].im ) )
                        * normFacSq;

		params->rhosqVec->data->data[j] = rhoSq;
	}

  }

  /*
   * Looking for event in output
   * Writing out to SnglInspiralTable
   *
   */


  /* will need to set up ignoreIndex and deltaEventIndex */
  /* for time being... */
  /* temporarily set chirpTime equal to 0.5 seconds */
 /* chirpTime = 0.5;
  deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );*/

  /* Calculate deltaEventIndex : the only acceptable clustering */
  /* is "window" method, for BCV                                */
  if ( params->clusterMethod == FindChirpClustering_window )
  {
	 deltaEventIndex=(UINT4) rint((params->clusterWindow/params->deltaT)+1.0);
  }
  else if ( params->clusterMethod == FindChirpClustering_tmplt )
  {
        ABORT( status, FINDCHIRPBCVSPINH_ECLUW, FINDCHIRPBCVSPINH_MSGECLUW );
  }

  /* ignore corrupted data at start and end */
  ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

  /* fprintf (stdout, "ignoreIndex = %d\n", ignoreIndex);*/

/*  if ( ignoreIndex > numPoints / 4 )
  {
  	ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
  }*/
  /* XXX reset ignoreIndex to one quarter of a segment XXX */
  ignoreIndex = numPoints / 4;

  /* REMOVE THIS */
  /*  ignoreIndex = 0; */
  /* REMOVE THIS */


 /* fprintf (stdout, "ignoreIndex = %d\n", ignoreIndex);*/


  rhosqThresh = params->rhosqThresh;
  modqsqThresh = rhosqThresh;

/*  fprintf (stdout, "beta = %e\n", input->fcTmplt->tmplt.beta);*/

  if (input->fcTmplt->tmplt.beta == 0.0)
  {
          rhosqThresh = 125;
	  modqsqThresh = 125;
  	/*  fprintf (stdout, "beta = 0 so changing rhosq thresh = %e\n ", rhosqThresh);*/
  }

/*  fprintf (stdout, "modqsqThresh = %e\n", modqsqThresh);
  fprintf (stdout, "rhosqThresh = %e\n",  rhosqThresh); */


  /* look for an event in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
         REAL4 rhoSq = (
		  ( q[j].re * q[j].re ) +
		  ( q[j].im * q[j].im ) +
                  ( qBCVSpin1[j].re * qBCVSpin1[j].re ) +
		  ( qBCVSpin1[j].im * qBCVSpin1[j].im ) +
                  ( qBCVSpin2[j].re * qBCVSpin2[j].re ) +
		  ( qBCVSpin2[j].im * qBCVSpin2[j].im ) )
                  * normFacSq;

	rho    = pow(rhoSq, 0.5);

	invRho = 1/rho;



	/* if snrsq exceeds threshold at any point */
   	if ( rhoSq > modqsqThresh )
    	{

      /*  fprintf (stdout, "\n");
        fprintf (stdout, "normFac = %e\n", normFac);
        fprintf (stdout, "normFacSq = %e\n", normFacSq);
	fprintf (stdout, "rhoSq = %e\n", rhoSq);
        fprintf (stdout, "q[j].re = %e\n", q[j].re );
	fprintf (stdout, "q[j].im = %e\n", q[j].im );
	fprintf (stdout, "qBCVSpin1[j].re = %e\n", qBCVSpin1[j].re );
	fprintf (stdout, "qBCVSpin1[j].im = %e\n", qBCVSpin1[j].im );
	fprintf (stdout, "qBCVSpin2[j].re = %e\n", qBCVSpin2[j].re );
	fprintf (stdout, "qBCVSpin2[j].im = %e\n", qBCVSpin2[j].im );
        fprintf (stdout, "rho = %e\n", rho);
        fprintf (stdout, "invRho = %e\n", invRho);*/

		if (! *eventList )    /* if *eventlist is empty */
                {
 			/* store the start of the crossing */
          		eventStartIdx = j;

          		/* if this is the first event, start the list */
          		thisEvent = *eventList = (SnglInspiralTable *)
            		LALCalloc( 1, sizeof(SnglInspiralTable) );
          		if ( ! thisEvent )
          		{
            			ABORT( status, FINDCHIRPH_EALOC,
				FINDCHIRPH_MSGEALOC );
          		}

          		/* record the data that we need
			for the clustering algorithm */
          		thisEvent->end_time.gpsSeconds = j;
          		thisEvent->snr = rho;

			alpha1hat = q[j].re * invRho * normFac;
 	 	 	alpha4hat = q[j].im * invRho * normFac;
			alpha2hat = qBCVSpin1[j].re * invRho * normFac;
			alpha5hat = qBCVSpin1[j].im * invRho * normFac;
			alpha3hat = qBCVSpin2[j].re * invRho * normFac;
			alpha6hat = qBCVSpin2[j].im * invRho * normFac;
		/*													                           		     fprintf (stdout, "alpha1hat = %e\n", alpha1hat);
		        fprintf (stdout, "alpha2hat = %e\n", alpha2hat);
			fprintf (stdout, "alpha3hat = %e\n", alpha3hat);
			fprintf (stdout, "alpha4hat = %e\n", alpha4hat);
			fprintf (stdout, "alpha5hat = %e\n", alpha5hat);
			fprintf (stdout, "alpha6hat = %e\n", alpha6hat);*/
                			                   									                                     thisEvent->alpha1 = alpha1hat;
		        thisEvent->alpha2 = alpha2hat;
			thisEvent->alpha3 = alpha3hat;
			thisEvent->alpha4 = alpha4hat;
			thisEvent->alpha5 = alpha5hat;
			thisEvent->alpha6 = alpha6hat;

		}

                /* check to see if snr>threshold
		within interval defined by
                deltaEventIndex */

		else if ( ! params->clusterMethod == FindChirpClustering_none &&
	        j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
		rho > thisEvent->snr )
        	{
          		/* this is the same event so update maximum */
          		thisEvent->end_time.gpsSeconds = j;
          		thisEvent->snr = rho;

                         alpha1hat = q[j].re * invRho * normFac;
			 alpha4hat = q[j].im * invRho * normFac;
			 alpha2hat = qBCVSpin1[j].re * invRho * normFac;
			 alpha5hat = qBCVSpin1[j].im * invRho * normFac;
			 alpha3hat = qBCVSpin2[j].re * invRho * normFac;
			 alpha6hat = qBCVSpin2[j].im * invRho * normFac;
			                                                                                                                                          /*    fprintf (stdout, "alpha1hat = %e\n", alpha1hat);
			 fprintf (stdout, "alpha2hat = %e\n", alpha2hat);
			 fprintf (stdout, "alpha3hat = %e\n", alpha3hat);
			 fprintf (stdout, "alpha4hat = %e\n", alpha4hat);
			 fprintf (stdout, "alpha5hat = %e\n", alpha5hat);
			 fprintf (stdout, "alpha6hat = %e\n", alpha6hat);*/
			                                                                                                                                              thisEvent->alpha1 = alpha1hat;
			 thisEvent->alpha2 = alpha2hat;
			 thisEvent->alpha3 = alpha3hat;
			 thisEvent->alpha4 = alpha4hat;
			 thisEvent->alpha5 = alpha5hat;
			 thisEvent->alpha6 = alpha6hat;

                }

                else if (j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
		         params->clusterMethod == FindChirpClustering_none )
        	{
          		/* clean up this event */
          		SnglInspiralTable *lastEvent;
          		INT8 timeNS;
          		INT4 timeIndex = thisEvent->end_time.gpsSeconds;

          		/* set the event LIGO GPS time of the event */
          		timeNS = 1000000000L *
            		(INT8) (input->segment->data->epoch.gpsSeconds);
          		timeNS +=
			(INT8) (input->segment->data->epoch.gpsNanoSeconds);
          		timeNS += (INT8) (1e9 * timeIndex * deltaT);
          		thisEvent->end_time.gpsSeconds
				= (INT4) (timeNS/1000000000L);
          		thisEvent->end_time.gpsNanoSeconds
				= (INT4) (timeNS%1000000000L);
          		LALGPStoGMST1( status->statusPtr,
				&(thisEvent->end_time_gmst),
              			&(thisEvent->end_time), &gmstUnits );
          		CHECKSTATUSPTR( status );

          		/* set the impulse time for the event */
          		thisEvent->template_duration = (REAL8) chirpTime;

			/* record the ifo and channel name for the event */
          		strncpy( thisEvent->ifo, input->segment->data->name,
              			2 * sizeof(CHAR) );
          		strncpy( thisEvent->channel,
				input->segment->data->name + 3,
              			(LALNameLength - 3) * sizeof(CHAR) );
         		thisEvent->impulse_time = thisEvent->end_time;

                	/* record the beta value */
               		/* eventually beta will be provided FROM
				the template bank */
                	thisEvent->beta   = input->fcTmplt->tmplt.beta;
          		thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;
          		thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;
          		/* chirp mass in units of M_sun */
          		thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
            		pow( 3.0 / 128.0 / input->fcTmplt->tmplt.psi0 , 3.0/5.0 );
          		m =  fabs(thisEvent->psi3) /
	                (16.0 * LAL_MTSUN_SI * LAL_PI
			* LAL_PI * thisEvent->psi0) ;
          	        thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 *
              		pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );
          		thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal ;

          		/* set the type of the template used in the analysis */
          		LALSnprintf( thisEvent->search,
				LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              			"FindChirpBCVSpin" );

                	/* commented out all chisq stuff */

          		/* compute the time since the snr crossing */
          		thisEvent->event_duration =
            		(REAL8) timeIndex - (REAL8) eventStartIdx;
          		thisEvent->event_duration *= (REAL8) deltaT;

          		/* store the start of the crossing */
          		eventStartIdx = j;

			/* allocate memory for the newEvent */
          		lastEvent = thisEvent;

          		lastEvent->next = thisEvent = (SnglInspiralTable *)
           		LALCalloc( 1, sizeof(SnglInspiralTable) );
          		if ( ! lastEvent->next )
          		{
            		ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
         		}

          		/* stick minimal data into the event */
          		thisEvent->end_time.gpsSeconds = j;
          		thisEvent->snr = rho;

  			alpha1hat = q[j].re * invRho * normFac;
                        alpha4hat = q[j].im * invRho * normFac;
			alpha2hat = qBCVSpin1[j].re * invRho * normFac;
			alpha5hat = qBCVSpin1[j].im * invRho * normFac;
			alpha3hat = qBCVSpin2[j].re * invRho * normFac;
			alpha6hat = qBCVSpin2[j].im * invRho * normFac;
															                  /*                         fprintf (stdout, "alpha1hat = %e\n", alpha1hat);
			fprintf (stdout, "alpha2hat = %e\n", alpha2hat);
			fprintf (stdout, "alpha3hat = %e\n", alpha3hat);
			fprintf (stdout, "alpha4hat = %e\n", alpha4hat);
			fprintf (stdout, "alpha5hat = %e\n", alpha5hat);
			fprintf (stdout, "alpha6hat = %e\n", alpha6hat);*/
															                          		     thisEvent->alpha1 = alpha1hat;
			thisEvent->alpha2 = alpha2hat;
			thisEvent->alpha3 = alpha3hat;
			thisEvent->alpha4 = alpha4hat;
			thisEvent->alpha5 = alpha5hat;
			thisEvent->alpha6 = alpha6hat;

  	        }
        }
  }

 /*
  *
  * clean up the last event if there is one
  *
  */

 if (thisEvent)
 {
 	/* clean up this event */
	INT8 timeNS;
	INT4 timeIndex = thisEvent->end_time.gpsSeconds;

	/* set the event LIGO GPS time of the event */
	timeNS = 1000000000L *
	    (INT8) (input->segment->data->epoch.gpsSeconds);
	timeNS +=
	   (INT8) (input->segment->data->epoch.gpsNanoSeconds);
	timeNS += (INT8) (1e9 * timeIndex * deltaT);
	thisEvent->end_time.gpsSeconds
	        = (INT4) (timeNS/1000000000L);
	thisEvent->end_time.gpsNanoSeconds
	        = (INT4) (timeNS%1000000000L);
	LALGPStoGMST1( status->statusPtr,
	          &(thisEvent->end_time_gmst),
	          &(thisEvent->end_time), &gmstUnits );
	          CHECKSTATUSPTR( status );

	/* set the impulse time for the event */
	thisEvent->template_duration = (REAL8) chirpTime;

	/* record the ifo and channel name for the event */
	strncpy( thisEvent->ifo, input->segment->data->name,
	         2 * sizeof(CHAR) );
	         strncpy( thisEvent->channel,
	         input->segment->data->name + 3,
	         (LALNameLength - 3) * sizeof(CHAR) );
	         thisEvent->impulse_time = thisEvent->end_time;
																																	                          thisEvent->beta   = input->fcTmplt->tmplt.beta;
	thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;
	thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;
																                     /* chirp mass in units of M_sun */
	thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
	               pow( 3.0 / 128.0 / input->fcTmplt->tmplt.psi0 , 3.0/5.0 );
	               m =  fabs(thisEvent->psi3) /
	               (16.0 * LAL_MTSUN_SI * LAL_PI
	               * LAL_PI * thisEvent->psi0) ;

	thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 *
		       pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );
        thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal ;

       /* set the type of the template used in the analysis */
       LALSnprintf( thisEvent->search,
                    LIGOMETA_SEARCH_MAX * sizeof(CHAR),
                   "FindChirpBCVSpin" );

        /* compute the time since the snr crossing */
        thisEvent->event_duration =
                   (REAL8) timeIndex - (REAL8) eventStartIdx;
                   thisEvent->event_duration *= (REAL8) deltaT;

    }




  DETATCHSTATUSPTR( status );
  RETURN( status );
  }
