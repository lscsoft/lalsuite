/* <lalVerbatim file="GenerateInspiralCV">
 * Author: Thomas Cokelaer
 * $Id$
 * </lalVerbatim>  */

/*
 *   <lalLaTeX>
 *     	\subsection{Module \texttt{GenerateInspiral.c}}
 *     	\label{ss:GenerateInspiral.c}
        \noindent Generates a CoherentGW inspiral waveform for injection.


 *
 *     	\subsubsection*{Prototypes}
 *     	\vspace{0.1in}
 *     	\input{LALGenerateInspiralCP}         
        \input{LALGenerateInspiralGetApproxFromStringCP}
	\input{LALGenerateInspiralGetOrderFromStringCP}
	\input{LALGenerateInspiralGetModelFromStringCP}
	\input{LALGenerateInspiralPopulatePPNCP}
	\input{LALGenerateInspiralPopulateInspiralCP}
	\input{LALGenerateInspiralPopulateInspiralSpinCP}


 *     	\idx{LALGenerateInspiral}
        \idx{LALGenerateInspiralGetApproxFromString}
	\idx{LALGenerateInspiralGetOrderFromString}
	\idx{LALGenerateInspiralGetModelFromString}
	\idx{LALGenerateInspiralPopulatePPN}
	\idx{LALGenerateInspiralPopulateInspiral}
	\idx{LALGenerateInspiralPopulateInspiralSpin}
 *

 *	\begin{description}
 * 	\item[\texttt{LALGenerateInspiral()}] create an inspiral binary waveform generated 
 either by the \texttt{inspiral} package (EOB, PadeT1, TaylorT1, TaylorT2, TaylorT3, 
 SpinTaylor) or the \texttt{inject} package (GeneratePPN). It is used in the module 
 \texttt{FindChirpSimulation} in \texttt{findchirp} package. 
 *
 *      There are three  parsed arguments
 \begin{itemize}
 \item a \texttt{CoherentGW}  structure which stores amplitude, frequency and phase of the 
 waveform (output)
 \item a \texttt{thisEvent}  structure which provides some waveform parameters (input)
 \item a \texttt{PPNParamStruc} which gives some input parameters needed by the PPN waveform 
 generation. That arguments is also used as an output by all the different approximant (output/input).
 \end{itemize}

 The input must be composed of a valid thisEvent structure as well as the variable deltaT of the 
 PPNparamsStruct. All others variables of the PPNParamStruc are populated within that function.


 
 *  	\item[\texttt{LALGenerateInspiralGetModelFromString()}] convert a string  
 *  	provided by the \texttt{CoherentGW} structure in order to retrieve the order 
 *  	and approximant of the waveform to generate. 

 *  	\item[\texttt{LALGenerateInspiralPopulatePPN()}] Populate the PPNParamsStruc
 with the input argument \texttt{thisEvent}. That structure is used by both inspiral waveforms
 inject waveforms

 *  	\item[\texttt{LALGenerateInspiralPopulateInspiral()}]  Populate the InspiralTemplate strucuture
if the model chosen belongs to the inspiral package.

 *  	\item[\texttt{LALGenerateInspiralPopulateInspiralSpin()}] Populate some extra InspiralTemplate 
parameters related to spin if the model is SpinTaylor. That functions might be merge with 
\texttt{LALGenerateInspiralPopulatePPN()}
\end{description}
 *  

 *  
 *  	\subsubsection*{Algorithm}
 *  	\noindent None.
 *  
 *  	\subsubsection*{Notes}
 *  	Inject only time-domain waveforms for the time being such as PPN, 
 *  	TaylorT1, TaylorT2, TaylorT3, PadeT1 and EOB , Spintaylor..
 *  	\subsubsection*{Uses}
 *  	\begin{verbatim}
 *  	
 LALInspiralWaveLength()            
 LALInspiralParameterCalc()            
 LALInspiralWave1ForInjection()
 LALInspiralWave2ForInjection()
 LALInspiralWave3ForInjection()
 LALEOBWaveformForInjection()
 LALSTPNWaveformForInjection()
 LALCalloc()
 LALFree()
      
 *  	\end{verbatim}
 *  
 *  	\vfill{\footnotesize\input{GenerateInspiralCV}}
 *  	</lalLaTeX> 
 * */


#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>



NRCSID( GENERATEINSPIRALC, "$Id$" );



/* <lalVerbatim file="LALGenerateInspiralCP"> */
void LALGenerateInspiral(LALStatus             *status,
			 CoherentGW            *waveform,
			 SimInspiralTable      *thisEvent,
			 PPNParamStruc         *ppnParams)
/* </lalVerbatim> */
{

  UINT4  order       	= -1;  /* order of the model        */
  UINT4	 approximant 	= -1;  /* and its approximant value */
  UINT4  n;                    /* length variable*/

  InspiralTemplate      inspiralParams; /* if we want to use the inspiral package, we need
					   to fill that structure*/
  
  INITSTATUS(status, "LALGenerateInspiral",GENERATEINSPIRALC);
  ATTATCHSTATUSPTR(status);

  /* NO ASSERT FOR THE TIME BEING. WE ASSUME THAT ASSERTS ARE PRESENT IN THE
   * FINDCHIRP FILE OR INSPIRAL PACKAGE. */
  
  
  /* Get the order and approximant from the thisEvent table. In this table 
   * order and approximant are in a string format.  we need them as integers
   * to be used in the switch 
   * */
  LALGenerateInspiralGetModelFromString(status->statusPtr,					
					thisEvent->waveform,
					&order,
					&approximant);  


  /* --- switch for the model approximant  --- */
  if  (approximant == PPN )
    {
      LALGenerateInspiralPopulatePPN(status->statusPtr, ppnParams, thisEvent);   
      CHECKSTATUSPTR(status);    
      LALGeneratePPNInspiral( status->statusPtr, waveform, ppnParams );  
      CHECKSTATUSPTR(status);    
    }
  else
    {
      
      inspiralParams.approximant = approximant;
      inspiralParams.order       = order;

      LALGenerateInspiralPopulatePPN(status->statusPtr, ppnParams, thisEvent);
      CHECKSTATUSPTR(status);      

      LALGenerateInspiralPopulateInspiral(status->statusPtr, &inspiralParams, thisEvent, ppnParams);
      CHECKSTATUSPTR(status);      
      
      LALInspiralWaveLength(status->statusPtr, &n, inspiralParams);
      CHECKSTATUSPTR(status);      
      
      LALInspiralParameterCalc(status->statusPtr, &inspiralParams);
      CHECKSTATUSPTR(status);         
      
      switch (approximant){
      case TaylorT1: 
      case PadeT1:
	LALInspiralWave1ForInjection(status->statusPtr, 
				     waveform,
				     &inspiralParams,
				     ppnParams);
	CHECKSTATUSPTR(status);      
	break;
      case TaylorT2: 
	LALInspiralWave2ForInjection(status->statusPtr,
				     waveform, 
				     &inspiralParams,
				     ppnParams);
	CHECKSTATUSPTR(status);      
	break;
      case TaylorT3: 
	LALInspiralWave3ForInjection(status->statusPtr,
				     waveform, 
				     &inspiralParams,
				     ppnParams);
	CHECKSTATUSPTR(status);      
	break;    
      case SpinTaylor:
	LALGenerateInspiralPopulateInspiralSpin(status ->statusPtr, 
						&inspiralParams);		
	CHECKSTATUSPTR(status);      
	LALSTPNWaveformForInjection(status->statusPtr,
				    waveform,
				    &inspiralParams,
				    ppnParams);
	CHECKSTATUSPTR(status);      
	break;
      case EOB:       
	LALEOBWaveformForInjection(status->statusPtr, 
				   waveform, 
				   &inspiralParams,
				   ppnParams);
	CHECKSTATUSPTR(status);      
	break;
      }
    }


  /*  if (approximant==SpinTaylorT3)  ComputeSpin(&inspiralParams);*/



  
  DETATCHSTATUSPTR( status );
  RETURN(status);
}





/* <lalVerbatim file="LALGenerateInspiralGetModelFromStringCP"> */
void LALGenerateInspiralGetModelFromString(LALStatus *status,
					   CHAR *message,
					   UINT4 *order,
					   UINT4 *model) 
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALGenerateInspiralGetOrderFromString", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( message, status, GENERATEINSPIRALH_ENULL, GENERATEINSPIRALH_MSGENULL);
  
  LALGenerateInspiralGetOrderFromString( status->statusPtr,
					 message,
					 order);
  CHECKSTATUSPTR(status);      
  
  LALGenerateInspiralGetApproxFromString( status->statusPtr,
					  message,
					  model);
  CHECKSTATUSPTR(status);      
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  
}






/* <lalVerbatim file="LALGenerateInspiralGetOrderFromStringCP"> */
void LALGenerateInspiralGetOrderFromString(LALStatus *status,
					   CHAR *message,
					   UINT4 *order)
/* </lalVerbatim> */
{
  
  /* Function to transform the string of the PN order as well
     as the approximant into an integer value. */
  CHAR  *ptr 	= NULL;
  
  INITSTATUS( status, "LALGenerateInspiralGetOrderFromString", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( message, status, GENERATEINSPIRALH_ENULL, GENERATEINSPIRALH_MSGENULL);
  
  if (  (ptr = strstr(message, 		"newtonian") ) )
    {*order = newtonian;}
  else if (  (ptr = strstr(message, 	"oneHalfPN") ) )
    {*order = oneHalfPN;}
  else if (  (ptr = strstr(message,  	"onePN") ) )
    {*order = onePN;}
  else if (  (ptr =  strstr(message,  	"onePointFivePN") ) )
    {*order = onePointFivePN;}
  else if (  (ptr = strstr(message, 	"twoPN") ) )
    {*order = twoPN;}
  else if (  (ptr = strstr(message,  	"twoPointFivePN") ) )
    {*order = twoPointFivePN;}
  else if (  (ptr = strstr(message, 	"threePN") ) )
    {*order = threePN;}
  else if (  (ptr = strstr(message, 	"threePointFivePN") ) )
    {*order = threePointFivePN;}
  else {
    ABORT(status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);
  }
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="LALGenerateInspiralGetApproxFromStringCP"> */
void LALGenerateInspiralGetApproxFromString(LALStatus *status,
					    CHAR *message,
					    UINT4 *approximant)  
     /* </lalVerbatim> */
{
  
  /* Function to search for the approximant into a string */
  CHAR  *ptr 	= NULL;
  
  INITSTATUS( status, " LALGenerateInspiralGetApproxFromString", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( message, status, GENERATEINSPIRALH_ENULL, GENERATEINSPIRALH_MSGENULL);
 

  if (  (ptr = strstr(message, 		"TaylorT1" ) ) )
    {*approximant = TaylorT1;}
  else if (  (ptr = strstr(message, 	"TaylorT2" ) ) )
    {*approximant = TaylorT2;}
  else if (  (ptr = strstr(message, 	"SpinTaylorT3" ) ) )
    {*approximant = SpinTaylorT3;}
  else if (  (ptr = strstr(message, 	"TaylorT3" ) ) )
    {*approximant = TaylorT3;}
  else if (  (ptr = strstr(message, 	"EOB" ) ) )
    {*approximant = EOB;}
  else if (  (ptr = strstr(message, 	"SpinTaylor" ) ) )
    {*approximant = SpinTaylor;}
  else if (  (ptr = strstr(message, 	"PadeT1" ) ) )
    {*approximant = PadeT1;}
 else if (  (ptr = strstr(message, 	"GeneratePPN" ) ) )
    {*approximant = PPN;}

  /* those are not implemented for time injection */
  /*  else if (  (ptr = strstr(message, 	"TaylorF1" ) ) )
      {*approximant = TaylorF1;}
      else if (  (ptr = strstr(message, 	"TaylorF2" ) ) )
      {*approximant = TaylorF2;}
      else if (  (ptr = strstr(message, 	"PadeT1" ) ) )
      {*approximant = PadeT1;}
      else if (  (ptr = strstr(message, 	"PadeF1" ) ) )
      {*approximant = PadeF1;}
      else if (  (ptr = strstr(message, 	"BCV" ) ) )
      {*approximant = BCV;}
      else if (  (ptr = strstr(message, 	"TaylorCW" ) ) )
      {*approximant = TaylorCW;}
      else if (  (ptr = strstr(message, 	"SpinOrbitCW" ) ) )
      {*approximant = SpinOrbitCW;}
  */
  else {
    ABORT(status, GENERATEINSPIRALH_EAPPROX, GENERATEINSPIRALH_MSGEAPPROX);
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


/* <lalVerbatim file="LALGenerateInspiralPopulatePPNCP"> */
void  LALGenerateInspiralPopulatePPN(LALStatus             *status,
				     PPNParamStruc         *ppnParams,
				     SimInspiralTable      *thisEvent)
/* </lalVerbatim> */
{

  INITSTATUS( status, " LALGenerateInspiralGetApproxFromString", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );

  /* input fields */
  ppnParams->mTot = thisEvent->mass1 + thisEvent->mass2;
  ppnParams->eta  = thisEvent->eta;
  ppnParams->d    = thisEvent->distance * 1.0e6 * LAL_PC_SI;
  ppnParams->inc  = thisEvent->inclination;
  ppnParams->phi  = thisEvent->coa_phase;
  
  /* frequency cutoffs */
  ppnParams->fStartIn = 40.0;
  ppnParams->fStopIn  = -1.0 / 
    (6.0 * sqrt(6.0) * LAL_PI * ppnParams->mTot * LAL_MTSUN_SI);
  
  /* passed fields */
  ppnParams->position.longitude   = thisEvent->longitude;
  ppnParams->position.latitude    = thisEvent->latitude;
  ppnParams->position.system      = COORDINATESYSTEM_EQUATORIAL;
  ppnParams->psi                  = thisEvent->polarization;
  ppnParams->epoch.gpsSeconds     = 0;
  ppnParams->epoch.gpsNanoSeconds = 0;
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  
}


/* <lalVerbatim file="LALGenerateInspiralPopulateInspiralCP"> */
void LALGenerateInspiralPopulateInspiral(LALStatus             *status,
					 InspiralTemplate      *inspiralParams,
					 SimInspiralTable      *thisEvent,
					 PPNParamStruc        *ppnParams)
     
     /* </lalVerbatim> */
{

  INITSTATUS( status, " LALGenerateInspiralGetApproxFromString", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );

  /* --- Let's fill the inspiral structure now --- */
  inspiralParams->mass1	     =  thisEvent->mass1;  /* masses 1 */
  inspiralParams->mass2	     =  thisEvent->mass2;  /* masses 2 */
  inspiralParams->fLower	     =  ppnParams->fStartIn; /* lower cutoff frequency */
  inspiralParams->fCutoff	 = 1./ (ppnParams->deltaT)/2.-1; /* -1 to be in agreement with
								    the inspiral assert.*/
  inspiralParams->tSampling	 = 1./ (ppnParams->deltaT); /* sampling*/
  inspiralParams->signalAmplitude = 1.; 
  inspiralParams->distance	 =  thisEvent->distance * LAL_PC_SI * 1e6; 
  inspiralParams->startTime	 =  0.0;
  inspiralParams->startPhase	 =  thisEvent->coa_phase;
    
  inspiralParams->OmegaS	 = GENERATEINSPIRAL_OMEGAS;  /* EOB at 3PN contribution */
  inspiralParams->Theta	 = GENERATEINSPIRAL_THETA;   /* EOB at 3PN contribution */
  inspiralParams->Zeta2	 = GENERATEINSPIRAL_ZETA2;   /* EOB at 3PN contribution */

  inspiralParams->alpha	 = -1.;      /* bcv useless for the time being */
  inspiralParams->psi0	 = -1.;      /* bcv useless for the time being */
  inspiralParams->psi3	 = -1.;      /* bcv useless for the time being */
  inspiralParams->alpha1	 = -1.;      /* bcv useless for the time being */
  inspiralParams->alpha2	 = -1.;      /* bcv useless for the time being */
  inspiralParams->beta	 = -1.;      /* bcv useless for the time being */

  inspiralParams->inclination =  thisEvent->inclination ; /*inclination of the binary*/

  inspiralParams->ieta	   =  1;
  inspiralParams->nStartPad =  0;
  inspiralParams->nEndPad   =  0;
  
  inspiralParams->massChoice  = m1Andm2;      
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  
}



/* <lalVerbatim file="LALGenerateInspiralPopulateInspiralSpinCP"> */
void LALGenerateInspiralPopulateInspiralSpin(LALStatus             *status,
					     InspiralTemplate      *inspiralParams)
     /* </lalVerbatim> */
{
  INITSTATUS( status, "LALGenerateInspiralPopulateInspiralSpin", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );

  /* not needed ???*/
  inspiralParams->sourceTheta = GENERATEINSPIRAL_SOURCETHETA;
  inspiralParams->sourcePhi   = GENERATEINSPIRAL_SOURCEPHI;
  inspiralParams->spin1[0]    = GENERATEINSPIRAL_SPIN1; 
  inspiralParams->spin2[0]    = GENERATEINSPIRAL_SPIN2;
  inspiralParams->spin1[1]    = GENERATEINSPIRAL_THETA1;
  inspiralParams->spin2[1]    = GENERATEINSPIRAL_THETA2;
  inspiralParams->spin1[2]    = GENERATEINSPIRAL_PHI1;
  inspiralParams->spin2[2]    = GENERATEINSPIRAL_PHI2;
  inspiralParams->orbitTheta0 = 1.5;
  inspiralParams->orbitPhi0   = 0.0;

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}
