/*----------------------------------------------------------------------- 
 * 
 * File Name: GenerateInspiral.c
 *
 * Author: Thomas Cokelaer 
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="GenerateInspiralCV">
Author: Thomas Cokelaer
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{GenerateInspiral.c}}
\label{ss:GenerateInspiral.c}

\noindent Provides an interface between \texttt{inspiral} package as well as the 
\texttt{inject} package to be used in \texttt{FindChirpSimulation} in \texttt{findchirp}
 packages. Basically, it is used  for injecting chirps into data.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateInspiralCP}
\idx{(GetOrderAndApproximant)}
\idx{ComputeSpin()}

\begin{description}
\item[\texttt{LALGenerateInspiral()}] create an inspiral  waveform and fill a 
\texttt{CoherentGW} structure. It uses the \texttt{thisEvent} structure and/or
\textt{PPNParamsStruct} to fill input structure needed by inspiral package and/or 
inject package.

\item[\texttt{LALGetApproximantAndOrder()}] convert the string  in the \texttt{CoherentGW}
structure into a integer in order to get the order and approximant. 


\item[\texttt{LALComputeSpin()}] might be used later for computing some spin parameters. 

\end{description}

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Notes}
\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpSimulationCV}}
</lalLaTeX> 
#endif


#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>



NRCSID( GENERATEINSPIRALC, "$Id$" );

#define GENERATEINSPIRALC_ENORM 0
#define GENERATEINSPIRALC_ESUB  1
#define GENERATEINSPIRALC_EARG  2
#define GENERATEINSPIRALC_EVAL  3
#define GENERATEINSPIRALC_EFILE 4
#define GENERATEINSPIRALC_EMEM  5

#define GENERATEINSPIRALC_EAPPROXIMANT  6
#define GENERATEINSPIRALC_EORDER  7

#define GENERATEINSPIRALC_MSGENORM "Normal exit"
#define GENERATEINSPIRALC_MSGESUB  "Subroutine failed"
#define GENERATEINSPIRALC_MSGEARG  "Error parsing arguments"
#define GENERATEINSPIRALC_MSGEVAL  "Input argument out of valid range"
#define GENERATEINSPIRALC_MSGEFILE "Could not open file"
#define GENERATEINSPIRALC_MSGEMEM  "Out of memory"

#define GENERATEINSPIRALC_MSGEAPPROXIMANT  "Error parsing order"
#define GENERATEINSPIRALC_MSGEORDER        "Error parsing approximant"



#define GENERATEINSPIRAL_APPROXIMANT PPN
#define GENERATEINSPIRAL_M1          1.4
#define GENERATEINSPIRAL_M2         1.4
#define GENERATEINSPIRAL_DISTANCE   10  /*mpc*/
#define GENERATEINSPIRAL_TC          0.
#define GENERATEINSPIRAL_INCLINATION 0.
#define GENERATEINSPIRAL_DEC         0.
#define GENERATEINSPIRAL_RA          0.
#define GENERATEINSPIRAL_PSI         0.
#define GENERATEINSPIRAL_PHIC        0.
#define GENERATEINSPIRAL_FUPPER   ((GENERATEINSPIRAL_SAMPLING-1)/2.)   /*must be < sampling/2*/
#define GENERATEINSPIRAL_ZETA2       0.
#define GENERATEINSPIRAL_OMEGAS      0.
#define GENERATEINSPIRAL_THETA       0.
#define GENERATEINSPIRAL_T0          0.
#define GENERATEINSPIRAL_T1          1.
#define GENERATEINSPIRAL_T2          2.
#define GENERATEINSPIRAL_SAMPLING 2048.
#define GENERATEINSPIRAL_A1          1000.
#define GENERATEINSPIRAL_A2          1000.
#define GENERATEINSPIRAL_PHI0        0.
#define GENERATEINSPIRAL_F0         100.
#define GENERATEINSPIRAL_ARG         0
#define GENERATEINSPIRAL_UDOT        0.5
#define GENERATEINSPIRAL_RP          0.5
#define GENERATEINSPIRAL_E           0.   
#define GENERATEINSPIRAL_ALPHA       0.
#define GENERATEINSPIRAL_FBCV        GENERATEINSPIRAL_FUPPER
#define GENERATEINSPIRAL_DELTAT      0  /*default use deltaT = 1/sampling*/
#define GENERATEINSPIRAL_LENGTH   1024  
#define GENERATEINSPIRAL_LENGTHIN  1e9

#define GENERATEINSPIRAL_SOURCETHETA 1
#define GENERATEINSPIRAL_SOURCEPHI   2
#define GENERATEINSPIRAL_SPIN1       0.3
#define GENERATEINSPIRAL_SPIN2       0.7
#define GENERATEINSPIRAL_THETA1      0.86
#define GENERATEINSPIRAL_THETA2      0.86

#define GENERATEINSPIRAL_PHI1      1.53
#define GENERATEINSPIRAL_PHI2      1.53


void
LALGenerateInspiral(
		    LALStatus             *status,
		    CoherentGW            *waveform,
		    SimInspiralTable      *thisEvent,
		    PPNParamStruc         *ppnParamsInputOutput
		    )
{
  REAL8  zeta2		= GENERATEINSPIRAL_ZETA2;				/* EOB 3PN parameter	*/
  REAL8  omegas   	= GENERATEINSPIRAL_OMEGAS;				/* EOB 3PN parameter 	*/
  REAL8  theta    	= GENERATEINSPIRAL_THETA;
  /* We must retrieve order and approximant from thisEvent structure*/
  UINT4  order       	= -1;
  UINT4	 approximant 	= -1;
  UINT4  n;
  /* We need to parse some structure to the inspiral generation functions.
   * For the time being there are two strcutures related to either the 
   * inspiral package or the inject package.
   * */
  PPNParamStruc         *ppnParamsInput = NULL;
  InspiralTemplate      inspiralParams;
  

  
  INITSTATUS(status, "LALGenerateInspiral",GENERATEINSPIRALC);
  ATTATCHSTATUSPTR(status);


  /* Get the order and approximant from the thisEvent table. In this table 
   * order and approximant are in a string format. I need them as an integer
   * to be used in the switch and to used the standard enumerate structure
   * of the inspiral package
   * */
  LALGetApproximantAndOrder(status->statusPtr, thisEvent->waveform, &order, &approximant);  
  CHECKSTATUSPTR(status);

  /* Now that we have the order and the approximant, let's switch 
   * to the correct approximant and fill the needed structure
   * */
  switch (approximant)
    {
    case PPN:
      ppnParamsInput = ppnParamsInputOutput;
      /*
      ppnParamsInput.mTot = thisEvent->mass1 + thisEvent->mass2;
      ppnParamsInput.eta  = thisEvent->eta;
      ppnParamsInput.d    = thisEvent->distance * 1.0e6 * LAL_PC_SI;
      ppnParamsInput.inc  = thisEvent->inclination; 
      ppnParamsInput.phi  = thisEvent->coa_phase;   
      ppnParamsInput.fStartIn = GENERATEINSPIRAL_FLOWER;
      ppnParamsInput.fStopIn  = -1.0 / 
	(6.0 * sqrt(6.0) * LAL_PI * ppnParamsInput.mTot * LAL_MTSUN_SI); 


      ppnParamsInput.position.longitude   = thisEvent->longitude;
      ppnParamsInput.position.latitude    = thisEvent->latitude;  
      ppnParamsInput.position.system      = COORDINATESYSTEM_EQUATORIAL;
      ppnParamsInput.psi                  = thisEvent->polarization; 
      ppnParamsInput.epoch.gpsSeconds     = 0;
      ppnParamsInput.epoch.gpsNanoSeconds = 0;
      ppnParamsInput.ppn = NULL;
      ppnParamsInput.lengthIn = 0;
      ppnParamsInput.deltaT = ppnParamsInputOutput->deltaT;
      */
      break;

    case EOB:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
      inspiralParams.ieta        =  1;
      /* Masses */
      inspiralParams.mass1       =  thisEvent->mass1;
      inspiralParams.mass2       =  thisEvent->mass2;
      /* padding */
      inspiralParams.startTime   =  0.0;
      inspiralParams.startPhase  =  thisEvent->coa_phase;
      inspiralParams.nStartPad   =  0;
      inspiralParams.nEndPad     =  0;
      /* distance and amplitude factor*/
      inspiralParams.signalAmplitude = 1.; 
      inspiralParams.distance    =  thisEvent->distance * LAL_PC_SI * 1e6; 
      /* lower and upper frequency*/
      inspiralParams.fLower      =  ppnParamsInputOutput->fStartIn;
      inspiralParams.fCutoff     = 1./ (ppnParamsInputOutput->deltaT)/2.-1;;
      /* Order and approximant */
      inspiralParams.order       =  order;
      inspiralParams.approximant = approximant;
      /* orbit plane? */
      inspiralParams.inclination =  thisEvent->inclination ;
      inspiralParams.massChoice=m1Andm2;
      /* Sampling */
      inspiralParams.tSampling = 1./ (ppnParamsInputOutput->deltaT);
      /* For EOB */
      inspiralParams.OmegaS      = omegas;
      inspiralParams.Theta       = theta;
      inspiralParams.Zeta2       = zeta2;
      break;
    case SpinTaylorT3:
      inspiralParams.ieta        = 1;
      inspiralParams.mass1       = thisEvent->mass1;
      inspiralParams.mass2       = thisEvent->mass2;
      inspiralParams.startTime   = 0.0;
      inspiralParams.startPhase  = thisEvent->coa_phase;
      inspiralParams.nStartPad   = 0;
      inspiralParams.nEndPad     = 0;
      inspiralParams.signalAmplitude = 1.; 
      inspiralParams.distance    = thisEvent->distance * LAL_PC_SI * 1e6; 
      inspiralParams.fLower      = ppnParamsInputOutput->fStartIn;
      inspiralParams.fCutoff     = 1./ (ppnParamsInputOutput->deltaT)/2.-1;;
      inspiralParams.order       = order;
      inspiralParams.inclination = thisEvent->inclination ;
/*      inspiralParams.massChoice  = m1Andm2;*/
      inspiralParams.tSampling   = 1./ (ppnParamsInputOutput->deltaT);
      inspiralParams.OmegaS      = omegas;
      inspiralParams.Theta       = theta;
      inspiralParams.Zeta2       = zeta2;
      inspiralParams.approximant = approximant;
      inspiralParams.sourceTheta = GENERATEINSPIRAL_SOURCETHETA;
      inspiralParams.sourcePhi   = GENERATEINSPIRAL_SOURCEPHI;
      inspiralParams.spin1[0]    = GENERATEINSPIRAL_SPIN1; 
      inspiralParams.spin2[0]    = GENERATEINSPIRAL_SPIN2;
      inspiralParams.spin1[1]    = GENERATEINSPIRAL_THETA1;
      inspiralParams.spin2[1]    = GENERATEINSPIRAL_THETA2;
      inspiralParams.spin1[2]    = GENERATEINSPIRAL_PHI1;
      inspiralParams.spin2[2]    = GENERATEINSPIRAL_PHI2;

      /* Extra computation for spin parameters: */
      ComputeSpin(&inspiralParams);      

      break;
    }


  /* IF we're going to use inspiral package, we need to use this initialization 
   * functions 
   * */
  switch  (approximant){
  case EOB:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case PadeT1:
  case SpinTaylorT3:
    LALInspiralWaveLength(status->statusPtr, &n, inspiralParams);
    CHECKSTATUSPTR(status);      
    LALInspiralParameterCalc(status->statusPtr, &inspiralParams);
    CHECKSTATUSPTR(status);                  
    break;
  }
  
  /* Here is the main switch to get the waveform */
  switch(approximant)
    {
    case TaylorT1: 
      LALInspiralWave1ForInjection(status->statusPtr, waveform, &inspiralParams);
      CHECKSTATUSPTR(status);      
      break;
    case TaylorT2: 
      LALInspiralWave2ForInjection(status->statusPtr, waveform, &inspiralParams);
      CHECKSTATUSPTR(status);      
      break;
    case TaylorT3: 
      LALInspiralWave3ForInjection(status->statusPtr, waveform, &inspiralParams);
      CHECKSTATUSPTR(status);      
      break;
    case TaylorF1:
    case TaylorF2:      
    case PadeT1:
    case PadeF1:
    case BCV:
    case BCVSpin:
    case SpinTaylorT3:
      /*    LALInspiralSpinModulatedWaveForInjection(status->statusPtr, 
					       waveform,
					       &inspiralParams);
					       CHECKSTATUSPTR(status);      */
      break;
    case EOB:       
      LALEOBWaveformForInjection(status->statusPtr, waveform, &inspiralParams);
      CHECKSTATUSPTR(status);      
      break;
    case PPN:
      LALGeneratePPNInspiral(status->statusPtr, waveform, ppnParamsInput);	    	
      CHECKSTATUSPTR(status);   
      break;
      /*      case SpinOrbitCW:	   
	      LALGenerateSpinOrbitCW(status->statusPtr, waveform, &params->socw);
		CHECKSTATUSPTR(status);      
		break;
		case TaylorCW:
		LALGenerateTaylorCW(status->statusPtr, waveform, &params->taylorcw);
		CHECKSTATUSPTR(status);      
	break;*/
    default: 
      fprintf(stderr,"nothing to do (bad approximant?)");
      break;
    }

  /* Finally let;s fill the output ppn structure used in FindChirpSimulation.c */
  switch(approximant)
    {
    case TaylorT1: 
    case TaylorT2: 
    case TaylorT3: 
    case TaylorF1:
    case TaylorF2:      
    case PadeT1:
    case PadeF1:
    case BCV:
    case BCVSpin:
    case SpinTaylorT3:
    case EOB:       
      ppnParamsInputOutput->dfdt   = inspiralParams.tSampling;
      ppnParamsInputOutput->tc     = inspiralParams.tC;
      ppnParamsInputOutput->length = inspiralParams.nStartPad;
      break;
    case PPN:
      /* ppnParamsInputOutput->dfdt   = ppnParamsInput->dfdt;
      ppnParamsInputOutput->tc     = ppnParamsInput->tc;
      ppnParamsInputOutput->length = ppnParamsInput->length;*/
      break;
    default: 
      fprintf(stderr,"nothing to do (bad approximant?)");
      break;
    }
  
  
  DETATCHSTATUSPTR( status );
  RETURN(status);
}




/* <lalVerbatim file="GetApproximantAndOrder"> */
void
LALGetApproximantAndOrder(
			  LALStatus *status,
			  CHAR *waveform,
			  UINT4 *order,
			  UINT4 *approximant)
{
  /* Function to transform the string of the PN order as well
     as the approximant into an integer value. */
  CHAR  *orderptr 	= NULL;
  CHAR  *approxptr 	= NULL; 
  
  INITSTATUS( status, "LALGetApproximantAndOrder", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( waveform, status, 1,"Null pointer");
  
  if (  (orderptr = strstr(waveform, "newtonian") ) ){
    *order = newtonian;}
  else if (  (orderptr = strstr(waveform, "oneHalfPN") ) ){
    *order = oneHalfPN;}
  else if (  (orderptr = strstr(waveform,  "onePN") ) ){
    *order = onePN;}
  else if (  (orderptr =  strstr(waveform,  "onePointFivePN") ) ){
    *order = onePointFivePN;}
  else if (  (orderptr = strstr(waveform, "twoPN") ) ){
    *order = twoPN;}
  else if (  (orderptr = strstr(waveform,  "twoPointFive") ) ){
    *order = twoPointFivePN;}
  else if (  (orderptr = strstr(waveform, "threePN") ) ){
    *order = threePN;}
  else if (  (orderptr = strstr(waveform, "threePointFivePN") ) ){
    *order = threePointFivePN;}
  else {
      ABORT(status, GENERATEINSPIRALC_EORDER, GENERATEINSPIRALC_MSGEORDER);
  }


  if (  (approxptr = strstr(waveform, "TaylorT1" ) ) ){
    *approximant = TaylorT1;}
  else if (  (approxptr = strstr(waveform, "TaylorT2" ) ) ){
    *approximant = TaylorT2;}
  else if (  (approxptr = strstr(waveform, "SpinTaylorT3" ) ) ){
    *approximant = SpinTaylorT3;}
  else if (  (approxptr = strstr(waveform, "TaylorT3" ) ) ){
    *approximant = TaylorT3;}
  else if (  (approxptr = strstr(waveform, "TaylorF1" ) ) ){
    *approximant = TaylorF1;}
  else if (  (approxptr = strstr(waveform, "TaylorF2" ) ) ){
    *approximant = TaylorF2;}
  else if (  (approxptr = strstr(waveform, "PadeT1" ) ) ){
    *approximant = PadeT1;}
  else if (  (approxptr = strstr(waveform, "PadeF1" ) ) ){
    *approximant = PadeF1;}
  else if (  (approxptr = strstr(waveform, "EOB" ) ) ){
    *approximant = EOB;}
  else if (  (approxptr = strstr(waveform, "BCV" ) ) ){
    *approximant = BCV;}
  else if (  (approxptr = strstr(waveform, "GeneratePPN" ) ) ){
    *approximant = PPN;}
  else if (  (approxptr = strstr(waveform, "TaylorCW" ) ) ){
    *approximant = TaylorCW;}
  else if (  (approxptr = strstr(waveform, "SpinOrbitCW" ) ) ){
    *approximant = SpinOrbitCW;}
  else {
      ABORT(status, GENERATEINSPIRALC_EAPPROXIMANT, GENERATEINSPIRALC_MSGEAPPROXIMANT);
  }


  DETATCHSTATUSPTR( status );
  RETURN( status );
}





void ComputeSpin(InspiralTemplate *params)
{


  REAL8 mass1Sq;
  REAL8 mass2Sq;
  REAL8 spin1Frac;
  REAL8 spin2Frac;
  REAL8 spin1Phi;
  REAL8 spin2Phi;
  REAL8 spin1Theta;
  REAL8 spin2Theta;

  mass1Sq = pow(params->mass1*LAL_MTSUN_SI,2.L); 
  mass2Sq = pow(params->mass2*LAL_MTSUN_SI,2.L); 
  spin1Frac = params->spin1[0];  
  spin2Frac = params->spin2[0];

  params->sourceTheta = LAL_PI/6.L;
  params->sourcePhi   = LAL_PI/6.L;

  spin1Theta = params->spin1[1];
  spin2Theta = params->spin2[1];
  spin1Phi   = params->spin1[2];
  spin2Phi   = params->spin2[2]; 
  
  params->spin1[0] =  mass1Sq * spin1Frac * sin(spin1Theta) * cos(spin1Phi);
  params->spin1[1] =  mass1Sq * spin1Frac * sin(spin1Theta) * sin(spin1Phi);
  params->spin1[2] =  mass1Sq * spin1Frac * cos(spin1Theta);	 	 
  params->spin2[0] =  mass2Sq * spin2Frac * sin(spin2Theta) * cos(spin2Phi);
  params->spin2[1] =  mass2Sq * spin2Frac * sin(spin2Theta) * sin(spin2Phi);
  params->spin2[2] =  mass2Sq * spin2Frac * cos(spin2Theta);  

}
