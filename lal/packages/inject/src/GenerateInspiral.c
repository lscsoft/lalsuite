/* <lalVerbatim file="GenerateInspiralCV">
 * Author: Thomas Cokelaer
 * $Id$
 * </lalVerbatim>  */

/*
 *   <lalLaTeX>
 *     	\subsection{Module \texttt{GenerateInspiral.c}}
 *     	\label{ss:GenerateInspiral.c}
 *     	\noindent Provides an interface between \texttt{inspiral} package as well as the 
 *     	\texttt{inject} package with \texttt{FindChirpSimulation} in \texttt{findchirp}
 *     	packages. Basically, it is used  for injecting chirps into data through the 
 *     	FindChirpSimulation function of the \texttt{findchirp} package.
 *
 *     	\subsubsection*{Prototypes}
 *     	\vspace{0.1in}
 *     	\input{GenerateInspiralCP}         
 *     	\input{LALGetApproximantAndOrderFromStringCP}
 *     	\idx{LALGetApproximantAndOrderFromString}
 *     	\idx{ComputeSpin}
 *     	\idx{GenerateInspiral}
 *
 *	\begin{description}
 * 	\item[\texttt{LALGenerateInspiral()}] create an inspiral  waveform and fill a 
 *  	\texttt{CoherentGW} structure. It uses the \texttt{thisEvent} structure and/or
 *  	\texttt{PPNParamsStruct} to fill input structure needed by inspiral package and/or 
 *  	inject package.
 *  
 *  	\item[\texttt{LALGetApproximantAndOrderFromString()}] convert a string  
 *  	provided by the \texttt{CoherentGW} structure into  integers in order 
 *  	to get the order and approximant of the waveform to generate. 
 *  
 *  	\item[\texttt{LALComputeSpin()}] might be used later for computing some spin parameters. 
 *  	\end{description}
 *  
 *  	\subsubsection*{Algorithm}
 *  	\noindent None.
 *  
 *  	\subsubsection*{Notes}
 *  	Inject only time-domain waveforms for the time being such as PPN, 
 *  	TaylorT1, TaylorT2, TaylorT3, PadeT1 and EOB.
 *  	\subsubsection*{Uses}
 *  	\begin{verbatim}
 *  	LALCalloc()
 *  	LALFree()
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


/* Those default values might be removed or used later.
 let's keep same for the time being */
#define GENERATEINSPIRAL_TC          0.
#define GENERATEINSPIRAL_INCLINATION 0.
#define GENERATEINSPIRAL_DEC         0.
#define GENERATEINSPIRAL_RA          0.
#define GENERATEINSPIRAL_PSI         0.
#define GENERATEINSPIRAL_PHIC        0.

/* idem for CW waveform*/
#define GENERATEINSPIRAL_F0         100.
#define GENERATEINSPIRAL_ARG         0
#define GENERATEINSPIRAL_UDOT        0.5
#define GENERATEINSPIRAL_RP          0.5
#define GENERATEINSPIRAL_E           0.   
#define GENERATEINSPIRAL_ALPHA       0.

/* For the spinning case. might be changed later or include 
   in the injection itself ? */
#define GENERATEINSPIRAL_SOURCETHETA 1.
#define GENERATEINSPIRAL_SOURCEPHI   2.
#define GENERATEINSPIRAL_SPIN1       0.3
#define GENERATEINSPIRAL_SPIN2       0.7
#define GENERATEINSPIRAL_THETA1      0.86
#define GENERATEINSPIRAL_THETA2      0.86

#define GENERATEINSPIRAL_PHI1        1.53
#define GENERATEINSPIRAL_PHI2        1.53


/* <lalVerbatim file="GenerateInspiralCP"> */
void LALGenerateInspiral(LALStatus             *status,
			 CoherentGW            *waveform,
			 SimInspiralTable      *thisEvent,
			 PPNParamStruc         *ppnParams)
/* </lalVerbatim> */
{

  /* We must retrieve order and approximant from thisEvent structure */
  UINT4  order       	= -1;
  UINT4	 approximant 	= -1;

  /* length variable*/
  UINT4  n;

  /* We need to parse some structure to the inspiral generation functions.
   * For the time being there are two strcutures related to either the 
   * inspiral package or the inject package.
   * */
  InspiralTemplate      inspiralParams;
  

  /* NO ASSERT FOR THE TIME BEING. WE ASSUME THAT ASSERTS ARE PRESENT IN THE
   * FINDCHIRP FILE OR INSPIRAL PACKAGE. */
  INITSTATUS(status, "LALGenerateInspiral",GENERATEINSPIRALC);
  ATTATCHSTATUSPTR(status);


  /* Get the order and approximant from the thisEvent table. In this table 
   * order and approximant are in a string format. I need them as an integer
   * to be used in the switch and to used the standard enumerate structure
   * of the inspiral package.
   * We probably should add that function in a standard place to be used by 
   * everybody...
   * */
  LALGetApproximantAndOrderFromString(status->statusPtr,
				      thisEvent->waveform,
				      &order,
				      &approximant);  
  CHECKSTATUSPTR(status);


  /* Now that we have the order and the approximant, let's switch 
   * to the approximant and fill the needed structure.
   * */

  inspiralParams.massChoice = m1Andm2; /* common to all time-domain inspiral models */

  /* --- Let's fill inspiralTemplate structure given the ppn strcuture of 
   * FindChirpSimulation --- */
  switch (approximant)
    {
      /* if PPN then nothing special to do for the time being since the findchirpSimulation
       has already filled the structure. Take care of it: findchirpSimulation might erase the 
       ppn structure later. Then we'll have to add the ppn strucutre here.*/
    case PPN:
      break;
    case EOB:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case PadeT1:
      /* --- Let's fill the inspiral structure now --- */
      /* approximant and order already filled */
      inspiralParams.order       = order;
      inspiralParams.approximant = approximant;
      /* masses */
      inspiralParams.mass1	 =  thisEvent->mass1;
      inspiralParams.mass2	 =  thisEvent->mass2;
      /* frequency cutoff */
      inspiralParams.fLower	 =  ppnParams->fStartIn;
      inspiralParams.fCutoff	 = 1./ (ppnParams->deltaT)/2.-1; /* -1 to be in agreement with
								    the inspiral assert.*/
      /* sampling */
      inspiralParams.tSampling	 = 1./ (ppnParams->deltaT);
      /* distance and amplitude factor*/
      inspiralParams.signalAmplitude = 1.; 
      inspiralParams.distance	 =  thisEvent->distance * LAL_PC_SI * 1e6; 
      /* initial phi and time */
      inspiralParams.startTime	 =  0.0;
      inspiralParams.startPhase	 =  thisEvent->coa_phase;
      /* EOB at 3PN contribution */
      inspiralParams.OmegaS	 = GENERATEINSPIRAL_OMEGAS;
      inspiralParams.Theta	 = GENERATEINSPIRAL_THETA;
      inspiralParams.Zeta2	 = GENERATEINSPIRAL_ZETA2;      
      /* bcv useless for the time being */
      inspiralParams.alpha	 = -1.;
      inspiralParams.psi0	 = -1.;   
      inspiralParams.psi3	 = -1.;                
      inspiralParams.alpha1	 = -1.;    
      inspiralParams.alpha2	 = -1.;    
      inspiralParams.beta	 = -1.;      
      /* orientation */
      inspiralParams.inclination =  thisEvent->inclination ;
      /* others */
      inspiralParams.ieta	 =  1;
      inspiralParams.nStartPad	 =  0;
      inspiralParams.nEndPad	 =  0;
      /* orbit plane ? */
      inspiralParams.inclination = thisEvent->inclination ;
      inspiralParams.massChoice  = m1Andm2;      
      break;

    case SpinTaylorT3:
    case SpinTaylor:
      /* --- Let's fill the inspiral structure now --- */
      inspiralParams.order       = order;
      inspiralParams.approximant = approximant;  
      /* masses */
      inspiralParams.mass1	 =  thisEvent->mass1;
      inspiralParams.mass2	 =  thisEvent->mass2;
      /* frequency cutoff */
      inspiralParams.fLower	 =  ppnParams->fStartIn;
      inspiralParams.fCutoff	 = 1./ (ppnParams->deltaT)/2.-1; /* -1 to be in agreement with
								    the inspiral assert.*/
      /* sampling */
      inspiralParams.tSampling	 = 1./ (ppnParams->deltaT);
      /* distance and amplitude factor*/
      inspiralParams.signalAmplitude = 1.; 
      inspiralParams.distance	 =  thisEvent->distance * LAL_PC_SI * 1e6; 
      /* initial phi and time */
      inspiralParams.startTime	 =  0.0;
      inspiralParams.startPhase	 =  thisEvent->coa_phase;
      /* EOB at 3PN contribution */
      inspiralParams.OmegaS	 = GENERATEINSPIRAL_OMEGAS;
      inspiralParams.Theta	 = GENERATEINSPIRAL_THETA;
      inspiralParams.Zeta2	 = GENERATEINSPIRAL_ZETA2;      
      /* bcv useless for the time being */
      inspiralParams.alpha	 = -1.;
      inspiralParams.psi0	 = -1.;   
      inspiralParams.psi3	 = -1.;                
      inspiralParams.alpha1	 = -1.;    
      inspiralParams.alpha2	 = -1.;    
      inspiralParams.beta	 = -1.;      
      /* orientation */
      inspiralParams.inclination =  thisEvent->inclination ;
      /* others */
      inspiralParams.ieta	 =  1;
      inspiralParams.nStartPad	 =  0;
      inspiralParams.nEndPad	 =  0;
      /* orbit plane ? */
      inspiralParams.inclination =  thisEvent->inclination ;
      inspiralParams.massChoice  = m1Andm2;      
      /* not needed ???*/
      inspiralParams.sourceTheta = GENERATEINSPIRAL_SOURCETHETA;
      inspiralParams.sourcePhi   = GENERATEINSPIRAL_SOURCEPHI;
      inspiralParams.spin1[0]    = GENERATEINSPIRAL_SPIN1; 
      inspiralParams.spin2[0]    = GENERATEINSPIRAL_SPIN2;
      inspiralParams.spin1[1]    = GENERATEINSPIRAL_THETA1;
      inspiralParams.spin2[1]    = GENERATEINSPIRAL_THETA2;
      inspiralParams.spin1[2]    = GENERATEINSPIRAL_PHI1;
      inspiralParams.spin2[2]    = GENERATEINSPIRAL_PHI2;
      inspiralParams.orbitTheta0 = 1.5;
      inspiralParams.orbitPhi0 	 = 0.0;

      break;

      /* must not reach that point */
    case PadeF1:
    case TaylorF1:
    case TaylorF2:
    case BCV:
    case BCVSpin:
      break;
    }

  if (approximant==SpinTaylorT3)
    ComputeSpin(&inspiralParams);
  /*otherwise it is computed inside the STPN spim taylor*/


  /* IF we're going to use inspiral package, we need to use those
   * additional initialization functions 
   * */
  switch  (approximant){
  case EOB:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case PadeT1:
  case SpinTaylorT3:
  case SpinTaylor:
    LALInspiralWaveLength(status->statusPtr, &n, inspiralParams);
    CHECKSTATUSPTR(status);      
    LALInspiralParameterCalc(status->statusPtr, &inspiralParams);
    CHECKSTATUSPTR(status);                  
    break;
    /* must not reach that point  */
  case PadeF1:
  case TaylorF1:
  case TaylorF2:
  case BCV:
  case BCVSpin:
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
    case PadeT1:
      LALInspiralWave1ForInjection(status->statusPtr, waveform, &inspiralParams);
      CHECKSTATUSPTR(status);    
      break;
    case SpinTaylorT3:
      /* not implemted yet
	 LALInspiralSpinModulatedWaveForInjection(status->statusPtr, 
	 waveform,
	 &inspiralParams);
	 CHECKSTATUSPTR(status);      */
      break;
    case SpinTaylor:
      LALSTPNWaveformForInjection(status->statusPtr,  waveform, &inspiralParams);
      CHECKSTATUSPTR(status);      
      break;
    case EOB:       
      LALEOBWaveformForInjection(status->statusPtr, waveform, &inspiralParams);
      CHECKSTATUSPTR(status);      
      break;
    case PPN:
      LALGeneratePPNInspiral(status->statusPtr, waveform, ppnParams);	    	
      CHECKSTATUSPTR(status);   
      break;
    case PadeF1:
    case TaylorF1:
    case TaylorF2:
    case BCV:
    case BCVSpin:
    default: /* must not reach that point*/
      break;
    }


  /* Finally let;s fill the output ppn structure used in FindChirpSimulation.c 
   1 - dfdt is computed by using the output of the GWcohehrent structure
   2 - tc is the same as tC given by inspiral package
   3 - However, length has no equivalent in the inspiral strucutre, therefore I 
   used nStartPad which seems to be useless for inspiral package. */
  switch(approximant)
    {
    case TaylorT1: 
    case TaylorT2: 
    case TaylorT3: 
    case PadeT1:
    case SpinTaylorT3:
    case SpinTaylor:
    case EOB:       

      ppnParams->dfdt   = 
	((REAL4)(waveform->f->data->data[inspiralParams.nStartPad-1] 
		 - waveform->f->data->data[inspiralParams.nStartPad-2]))
	/ inspiralParams.tSampling;    

      ppnParams->tc     = inspiralParams.tC;
      ppnParams->length = inspiralParams.nStartPad;
      inspiralParams.nStartPad     = 0;      
      break;
    case PPN:
      break;
    /* must not reach that point */
    case PadeF1:
    case TaylorF1:
    case TaylorF2:
    case BCV:
    case BCVSpin:
    default: 
      break;
    }
  
  
  DETATCHSTATUSPTR( status );
  RETURN(status);
}


                      
/* <lalVerbatim file="LALGetApproximantAndOrderFromStringCP"> */
void LALGetApproximantAndOrderFromString(LALStatus *status,
					 CHAR *waveform,
					 UINT4 *order,
					 UINT4 *approximant)    
{  /* </lalVerbatim> */

  /* Function to transform the string of the PN order as well
     as the approximant into an integer value. */
  CHAR  *orderptr 	= NULL;
  CHAR  *approxptr 	= NULL; 
  
  INITSTATUS( status, "LALGetApproximantAndOrder", GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( waveform, status, GENERATEINSPIRALH_ENULL, GENERATEINSPIRALH_MSGENULL);
  
  if (  (orderptr = strstr(waveform, 		"newtonian") ) )
  {*order = newtonian;}
  else if (  (orderptr = strstr(waveform, 	"oneHalfPN") ) )
  {*order = oneHalfPN;}
  else if (  (orderptr = strstr(waveform,  	"onePN") ) )
  {*order = onePN;}
  else if (  (orderptr =  strstr(waveform,  	"onePointFivePN") ) )
  {*order = onePointFivePN;}
  else if (  (orderptr = strstr(waveform, 	"twoPN") ) )
  {*order = twoPN;}
  else if (  (orderptr = strstr(waveform,  	"twoPointFivePN") ) )
  {*order = twoPointFivePN;}
  else if (  (orderptr = strstr(waveform, 	"threePN") ) )
  {*order = threePN;}
  else if (  (orderptr = strstr(waveform, 	"threePointFivePN") ) )
  {*order = threePointFivePN;}
  else {
    ABORT(status, GENERATEINSPIRALH_EORDER, GENERATEINSPIRALH_MSGEORDER);
  }
  
  
  if (  (approxptr = strstr(waveform, 		"TaylorT1" ) ) )
  {*approximant = TaylorT1;}
  else if (  (approxptr = strstr(waveform, 	"TaylorT2" ) ) )
  {*approximant = TaylorT2;}
  else if (  (approxptr = strstr(waveform, 	"SpinTaylorT3" ) ) )
  {*approximant = SpinTaylorT3;}
  else if (  (approxptr = strstr(waveform, 	"TaylorT3" ) ) )
  {*approximant = TaylorT3;}
  else if (  (approxptr = strstr(waveform, 	"TaylorF1" ) ) )
  {*approximant = TaylorF1;}
  else if (  (approxptr = strstr(waveform, 	"TaylorF2" ) ) )
  {*approximant = TaylorF2;}
  else if (  (approxptr = strstr(waveform, 	"PadeT1" ) ) )
  {*approximant = PadeT1;}
  else if (  (approxptr = strstr(waveform, 	"PadeF1" ) ) )
  {*approximant = PadeF1;}
  else if (  (approxptr = strstr(waveform, 	"EOB" ) ) )
  {*approximant = EOB;}
  else if (  (approxptr = strstr(waveform, 	"BCV" ) ) )
  {*approximant = BCV;}
  else if (  (approxptr = strstr(waveform, 	"SpinTaylor" ) ) )
  {*approximant = SpinTaylor;}
  else if (  (approxptr = strstr(waveform, 	"GeneratePPN" ) ) )
  {*approximant = PPN;}
  else if (  (approxptr = strstr(waveform, 	"TaylorCW" ) ) )
  {*approximant = TaylorCW;}
  else if (  (approxptr = strstr(waveform, 	"SpinOrbitCW" ) ) )
  {*approximant = SpinOrbitCW;}
  else {
      ABORT(status, GENERATEINSPIRALH_EAPPROXIMANT, GENERATEINSPIRALH_MSGEAPPROXIMANT);
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
