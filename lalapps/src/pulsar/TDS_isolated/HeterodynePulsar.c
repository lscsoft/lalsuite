/***********************************************************************/
/* Version of HeterodynePulsar.c code with binary timing added to      */
/* the heterodyne - Matt Pitkin 05/05/04                               */
/***********************************************************************/

/************************************ <lalVerbatim file="HeterodynePulsarCV">
Author: Dupuis, R. J.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{HeterodynePulsar.c}}

The routines in this module heterodyne, average and resample a time series for known pulsars.

\subsubsection*{Prototypes}
\input{HeterodynePulsarCP}
\idx{LALHeterodyneToPulsar()}
\idx{LALCalibrateFineHeterodyne()}
\idx{LALCoarseHeterodyne()}
\idx{LALFineHeterodyneToPulsar()}

\subsubsection*{Description}

\noindent The function \texttt{LALHeterodyneToPulsar()} takes as input 60 seconds of raw data
from the frame files. This data is then heterodyned with current phase phi(t) of the gravitational wave 
signal (other than an offset \phi_0).  Before averaging (binning) these 60*16384 points to obtain our 
point estimate Bk we apply three Butterworth IIR filters to prevent aliasing noise in from other frequencies.
Three low order filters are used instead of one higher order filter to prevent running into numerical errors (see 
the LAL tdfilter package for details).

\noindent The function \texttt{LALCalibrateFineHeterodyne()} simply calibrates the Bk's using the ASQ sensing
function, the open loop gain, and the alpha/beta parameters at the base frequecy of the gravitational wave signal.

\noindent The function \texttt{LALCoarseHeterodyne()} ... (not currently in use for S2/S3 analysis)

\noindent The function \texttt{LALFineHeterodyneToPulsar()} ... (not currently in use for S2/S3 analysis)

\subsubsection*{Algorithm}

To be completed.

\subsubsection*{Uses}
\begin{verbatim}
LALDIIRFilter()
LALZCreateVector()
LALZDestroyVector()
LALBarycenterEarth()
LALBarycenter()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{HeterodynePulsarCV}}

******************************************************* </lalLaTeX> */ 

/******* INCLUDE ANY LAL HEADERS ************/
//#include <lal/HeterodynePulsar.h>
#include "HeterodynePulsar.h"  /* use local version until up-to-date version is commited to lal */

#include <lal/LALConstants.h>
/******* DEFINE RCS ID STRING ************/
NRCSID( HETERODYNEPULSARC, "$Id$" );

/* <lalVerbatim file="HeterodynePulsarCP"> */
void
LALHeterodyneToPulsar     ( LALStatus                      *status,
		  	    HeterodyneOutput	   *output,
		   	    HeterodyneInput 	   *input,
		   	    HeterodyneParams	   *params )
/* </lalVerbatim> */
{
  UINT4		n; 
  UINT4		i;
  REAL8Vector   *phase;
  REAL8		f0,f1,f2;
  REAL8 	t,T0;
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime  emit;    
  REAL8		deltaT;
  COMPLEX16Vector   *Xh, *Yh;
  REAL8 tdt,ph;
  COMPLEX16 sum;
	
	BinaryPulsarOutput binOut; 
    
  INITSTATUS( status, "LALHeterodyneToPulsar", HETERODYNEPULSARC );
  ATTATCHSTATUSPTR (status); 

 /******* CHECK VALIDITY OF ARGUMENTS  ************/
 
  ASSERT(input != (HeterodyneInput *)NULL, status,
	HETERODYNEPULSARH_ENULLINPUT, HETERODYNEPULSARH_MSGENULLINPUT);
	
  ASSERT(output != (HeterodyneOutput *)NULL, status,
	HETERODYNEPULSARH_ENULLOUTPUT, HETERODYNEPULSARH_MSGENULLOUTPUT);

  ASSERT(params != (HeterodyneParams *)NULL, status,
	HETERODYNEPULSARH_ENULLPARAMS, HETERODYNEPULSARH_MSGENULLPARAMS); 

  ASSERT(( 983040 == input->V.data->length), status, 
	HETERODYNEPULSARH_ERLENGTH, HETERODYNEPULSARH_MSGERLENGTH);
	 
	 
  n = input->V.data->length; 
  
 /******* EXTRACT INPUTS AND PARAMETERS ************/
  
  deltaT = (REAL8)input->V.deltaT;
  
  /* set up spin parameters */
  f0 = input->f0;
  f1 = input->f1;
  f2 = input->f2;
 
 /* set up location of detector */
  baryinput.site.location[0] = params->detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = params->detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = params->detector.location[2]/LAL_C_SI;
 
 /* set up epoch, RA, DEC, and distance variables for LALBarycenter*/
  baryinput.tgps = input->V.epoch;
  baryinput.alpha = input->source.longitude;
  baryinput.delta = input->source.latitude;
  baryinput.dInv = 0.0;  /* no parallax */

/******* ALLOCATE MEMORY *************/
 
  phase = NULL;
  LALDCreateVector(status->statusPtr, &phase, n);
 
  Xh = NULL;
  LALZCreateVector(status->statusPtr, &Xh, n);
  
  Yh = NULL;
  LALZCreateVector(status->statusPtr, &Yh, n);
 /******* DO ANALYSIS ************/
 
 /* calculate instantaneous phase */
    for (i=0;i<n;i++)
    {    
      t = (REAL8)input->V.epoch.gpsSeconds + (REAL8)input->V.epoch.gpsNanoSeconds*1e-9 + (REAL8)i*deltaT;  
      baryinput.tgps.gpsSeconds = (UINT8)floor(t);
      baryinput.tgps.gpsNanoSeconds = (UINT8)floor((fmod(t,1.0)*1.e9));	 
      LALBarycenterEarth(status->statusPtr, &earth, &baryinput.tgps, params->edat); 
      LALBarycenter(status->statusPtr, &emit, &baryinput, &earth);   
			
			T0 = input->fEpochGPS;      
      
			/* if binary pulsar add extra time delay */
			if(!strcmp(params->binaryParams.model, "BT") || 
				 !strcmp(params->binaryParams.model, "ELL1") ||
				 !strcmp(params->binaryParams.model, "DD")){
				/* input SSB time into binary timing function */
			  params->binaryInput.tb = t + (REAL8)emit.deltaT;
			
				/* calculate binary time delay */
				LALBinaryPulsarDeltaT(status->statusPtr, &binOut, &params->binaryInput, &params->binaryParams);
			
				/* add binary time delay */
				tdt = t - T0 + (REAL8)emit.deltaT +  binOut.deltaT;
			}
			else{
			  tdt = t - T0 + (REAL8)emit.deltaT;
			}
			
      ph = f0*tdt + 0.5*f1*tdt*tdt + f2*tdt*tdt*tdt/6.0;
      /* store phase of time stamp i */
      phase->data[i] =  2.0*LAL_PI*fmod(ph, 1.0);
    }
  

  /* do heterodyning, iir filtering, and resampling */
 for (i = 0; i < n; i++)
 { 
   Xh->data[i].re = (REAL8)input->V.data->data[i]* cos(-phase->data[i]);	
	 Yh->data[i].re = LALDIIRFilter( Xh->data[i].re, params->iirFilter1Re ); 		    	    
	 Xh->data[i].re = LALDIIRFilter( Yh->data[i].re, params->iirFilter2Re ); 
   Yh->data[i].re = LALDIIRFilter( Xh->data[i].re, params->iirFilter3Re ); 

	 
   Xh->data[i].im = (REAL8)input->V.data->data[i] * sin(-phase->data[i]);
   Yh->data[i].im = LALDIIRFilter( Xh->data[i].im, params->iirFilter1Im ); 		    	    
   Xh->data[i].im = LALDIIRFilter( Yh->data[i].im, params->iirFilter2Im ); 
   Yh->data[i].im = LALDIIRFilter( Xh->data[i].im, params->iirFilter3Im );  
 }

 sum.re = 0.0;
 sum.im = 0.0;

/* sum (bin) heterodyned/filtered data */
 for (i=0;i<n;i++)
 {
   sum.re += Yh->data[i].re;
   sum.im += Yh->data[i].im;
 }
 
 /* divide by n to get average */   
 output->B.re = sum.re / (REAL8)n;
 output->B.im = sum.im / (REAL8)n;
 
 /****** CLEAN UP *******/
 LALDDestroyVector(status->statusPtr, &phase); 
 LALZDestroyVector(status->statusPtr, &Xh); 
 LALZDestroyVector(status->statusPtr, &Yh); 
 DETATCHSTATUSPTR(status); 
 RETURN(status);
}

/*****************************************************************************/


/* <lalVerbatim file="HeterodynePulsarCP"> */
void
LALCalibrateFineHeterodyne ( LALStatus                      *status,
		  	    CalibrateFineOutput		    *output,
		   	    CalibrateFineInput 	   	    *input,
		   	    CalibrateFineParams	            *params )
/* </lalVerbatim> */
{
  COMPLEX16 R;
    
  /* simply calculate imaginary response of IFO at this specific frequency bin 
  and adjust the Bk's correspondingly */ 
    
  R.re = (cos(params->phaseC0) + 
  input->alpha*input->beta*params->H0*cos(params->phaseH0 - params->phaseC0))/(input->alpha*params->C0);
  
  R.im = (-sin(params->phaseC0) + 
  input->alpha*input->beta*params->H0*sin(params->phaseH0 - params->phaseC0))/(input->alpha*params->C0);

  output->R.re = R.re;
  output->R.im = R.im;
  
  output->B.re = input->B.re*R.re - input->B.im*R.im;
  output->B.im = input->B.re*R.im + input->B.im*R.re;
  
  DETATCHSTATUSPTR(status); 
  RETURN(status);
}
/*****************************************************************************/
