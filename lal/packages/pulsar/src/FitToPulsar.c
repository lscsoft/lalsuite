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

/************************************ <lalVerbatim file="FitToPulsarCV">
Author: Dupuis, R. J.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{FitToPulsar.c}}
Calculates the best fit parameters for a GW signal originating from a
non-precessing pulsar.

\subsubsection*{Prototypes}
\idx{LALCoarseFitToPulsar()}
\idx{LALFineFitToPulsar}
\input{FitToPulsarCP}
\begin{verbatim}
void
LALFineFitToPulsar	(            LALStatus            *status,
                                FineFitOutput        *output,
                                FineFitInput         *input,
                                FineFitParams        *params )
\end{verbatim}

\subsubsection*{Description}

This routine calculates the best fit of parameters by minimizing $\chi^2$ by going through
fixed grid for $\iota, \psi, \phi_{0}$ and  $h_{0}$.  The best fit parameters
returned by \texttt{LALCoarseFitToPulsar()} are then used as initial parameters for
\texttt{LALFineFitToPulsar()}.

The function \texttt{LALFineFitToPulsar()} refines the fit using the Levenberg-Marquardt method for nonlinear fitting. This is
done by calculating the Hessian and the gradient of $\chi^2$ ...

\subsubsection*{Algorithm}
To be completed.

\subsubsection*{Uses}
\begin{verbatim}
LALSCreateVector()
LALSDestroyVector()
LALComputeDetAMResponse()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FitToPulsarCV}}

******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */
#include <math.h>

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/FitToPulsar.h>

/******* DEFINE RCS ID STRING ************/
NRCSID( FITTOPULSARC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/
#define INICHISQU 1.e50
#define INIMAXEH 0
#define INIMINEH 1e10

/******* DECLARE LOCAL (static) FUNCTIONS ************/
/* (definitions can go here or at the end of the file) */

/******* DEFINE GLOBAL FUNCTIONS ************/

/* <lalVerbatim file="FitToPulsarCP"> */
void
LALCoarseFitToPulsar	( 	LALStatus            *status,
		    		CoarseFitOutput      *output,
		    		CoarseFitInput       *input,
		    		CoarseFitParams      *params )

/* </lalVerbatim> */
{
  /******* DECLARE VARIABLES ************/

	UINT4 			n,i;		/* integer indices */
	REAL8 			psi;
	REAL8 			h0/*,h*/;
	REAL8 			cosIota;
	REAL8			phase;
	REAL8   		chiSquare;
	LALDetector 		detector;
	LALSource 		pulsar;
	LALDetAndSource 	detAndSource;
 	LALDetAMResponse 	computedResponse;
	REAL4Vector 		*Fp, *Fc;	/* pointers to amplitude responses of the detector */
	REAL8Vector		*var;
	REAL8 			cos2phase,sin2phase;
 	COMPLEX16		Xp, Xc;
	REAL8 			Y, cosIota2;
	COMPLEX16		A,B;
  	REAL8			sumAA, sumBB, sumAB;
	REAL8			h0BestFit=0,phaseBestFit=0, psiBestFit=0, cosIotaBestFit=0;
  	REAL8			minChiSquare;
	COMPLEX16		eh0;
	REAL8 oldMinEh0, oldMaxEh0, weh0;
	UINT4			iH0, iCosIota, iPhase, iPsi, arg;

  INITSTATUS(status, "LALCoarseFitToPulsar", FITTOPULSARC);
  ATTATCHSTATUSPTR(status);

  /******* CHECK VALIDITY OF ARGUMENTS  ************/
  ASSERT(output != (CoarseFitOutput *)NULL, status,
         FITTOPULSARH_ENULLOUTPUT, FITTOPULSARH_MSGENULLOUTPUT);

  ASSERT(params != (CoarseFitParams *)NULL, status,
         FITTOPULSARH_ENULLPARAMS, FITTOPULSARH_MSGENULLPARAMS);

  ASSERT(input != (CoarseFitInput *)NULL, status,
         FITTOPULSARH_ENULLINPUT, FITTOPULSARH_MSGENULLINPUT);

  ASSERT(input->B->length == input->var->length, status,
  FITTOPULSARH_EVECSIZE , FITTOPULSARH_MSGEVECSIZE );

  /******* EXTRACT INPUTS AND PARAMETERS ************/

  n = input->B->length;

  for (i = 0; i < n; i++)
     ASSERT(input->var->data[i].re > 0.0 && input->var->data[i].im > 0.0, status,
         FITTOPULSARH_EVAR, FITTOPULSARH_MSGEVAR);

  detector = params->detector;
  detAndSource.pDetector = &detector;

  pulsar = params->pulsarSrc;

  /******* ALLOCATE MEMORY TO VECTORS **************/

  Fp = NULL;
  LALSCreateVector(status->statusPtr, &Fp, n);
  Fp->length = n;

  Fc = NULL;
  LALSCreateVector(status->statusPtr, &Fc, n);
  Fc->length = n;

  var = NULL;
  LALDCreateVector(status->statusPtr, &var, n);
  var->length = n;
  /******* INITIALIZE VARIABLES *******************/

  minChiSquare = INICHISQU;
  oldMaxEh0 = INIMAXEH;
  oldMinEh0 = INIMINEH;


  for (i=0;i<n;i++)
  {
    var->data[i] = input->var->data[i].re + input->var->data[i].im;
  }
  /******* DO ANALYSIS ************/

  for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
  {
    psi = params->meshPsi[0] + iPsi*params->meshPsi[1];
    pulsar.orientation = psi;
    detAndSource.pSource = &pulsar;

   /* create vectors containing amplitude response of detector for specified times */
   for (i=0;i<n;i++)
   {
     LALGPSandAcc gpsAndAcc;
     gpsAndAcc.gps = input->t[i];
     gpsAndAcc.accuracy = LALLEAPSEC_STRICT; /* FIXME: check */
     LALComputeDetAMResponse(status->statusPtr, &computedResponse,&detAndSource, &gpsAndAcc);
     Fp->data[i] = (REAL8)computedResponse.plus;
     Fc->data[i] = (REAL8)computedResponse.cross;
   }

   for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
   {
     phase = params->meshPhase[0] + iPhase*params->meshPhase[1];
     cos2phase = cos(2.0*phase);
     sin2phase = sin(2.0*phase);
     for (iCosIota = 0; iCosIota < params->meshCosIota[2]; iCosIota++)
     {
       cosIota = params->meshCosIota[0] + iCosIota*params->meshCosIota[1];
       cosIota2 = 1.0 + cosIota*cosIota;
       Xp.re = cosIota2 * cos2phase;
       Xp.im = cosIota2 * sin2phase;

       Y = 2.0*cosIota;

       Xc.re = Y*sin2phase;
       Xc.im = -Y*cos2phase;

       sumAB = 0.0;
       sumAA = 0.0;
       sumBB = 0.0;
       eh0.re=0.0;
       eh0.im=0.0;

       for (i = 0; i < n; i++)
       {
	 B.re = input->B->data[i].re;
	 B.im = input->B->data[i].im;

	 A.re = Fp->data[i]*Xp.re + Fc->data[i]*Xc.re;
	 A.im = Fp->data[i]*Xp.im + Fc->data[i]*Xc.im;

	 sumBB += (B.re*B.re + B.im*B.im) / var->data[i];
	 sumAA += (A.re*A.re + A.im*A.im) / var->data[i];
	 sumAB += (B.re*A.re + B.im*A.im) / var->data[i];

         /**** calculate error on h0 **********/
	 eh0.re += (Fp->data[i]*cosIota2*cos2phase
                + 2.0*Fc->data[i]*cosIota*sin2phase)
	        * (Fp->data[i]*cosIota2*cos2phase
                + 2.0*Fc->data[i]*cosIota*sin2phase) / input->var->data[i].re;

         eh0.im += (Fp->data[i]*cosIota2*sin2phase
                - 2.0*Fc->data[i]*cosIota*cos2phase)
	        * (Fp->data[i]*cosIota2*sin2phase
                - 2.0*Fc->data[i]*cosIota*cos2phase) / input->var->data[i].im;
        }

	for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
        {
	  h0 = params->meshH0[0] + iH0*params->meshH0[1];
	  chiSquare = sumBB - 2.0*h0*sumAB + h0*h0*sumAA;

	  if (chiSquare<minChiSquare)
	  {
	    minChiSquare = chiSquare;
	    h0BestFit = h0;
	    cosIotaBestFit = cosIota;
	    psiBestFit = psi;
	    phaseBestFit = phase;
	    output->eh0[0] = 1.0 /sqrt(sqrt(eh0.re*eh0.re + eh0.im*eh0.im));
	  }

	  weh0 = 1.0 /sqrt(sqrt(eh0.re*eh0.re + eh0.im*eh0.im));

	  if (weh0>oldMaxEh0)
	  {
	    output->eh0[2] = weh0;
	    oldMaxEh0 = weh0;
	  }

	  if (weh0<oldMinEh0)
	  {
	    output->eh0[1] = weh0;
	    oldMinEh0 = weh0;
          }
	  arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi));
	  output->mChiSquare->data[arg] = chiSquare;
        }
      }
    }
  }

  /******* CONSTRUCT OUTPUT ************/

  output->h0 = h0BestFit;
  output->cosIota = cosIotaBestFit;
  output->phase = phaseBestFit;
  output->psi = psiBestFit;
  output->chiSquare = minChiSquare;


  ASSERT(minChiSquare < INICHISQU,  status,
  FITTOPULSARH_EMAXCHI , FITTOPULSARH_MSGEMAXCHI);

  LALSDestroyVector(status->statusPtr, &Fp);
  LALSDestroyVector(status->statusPtr, &Fc);
  LALDDestroyVector(status->statusPtr, &var);
  DETATCHSTATUSPTR(status);
  RETURN(status);


}
