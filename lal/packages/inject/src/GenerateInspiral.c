#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>



NRCSID( INJECTINSPIRALC, "$Id$" );

#define INJECTINSPIRALC_ENORM 0
#define INJECTINSPIRALC_ESUB  1
#define INJECTINSPIRALC_EARG  2
#define INJECTINSPIRALC_EVAL  3
#define INJECTINSPIRALC_EFILE 4
#define INJECTINSPIRALC_EMEM  5
#define INJECTINSPIRALC_EAPPROXIMANT  6
#define INJECTINSPIRALC_EORDER  7

#define INJECTINSPIRALC_MSGENORM "Normal exit"
#define INJECTINSPIRALC_MSGESUB  "Subroutine failed"
#define INJECTINSPIRALC_MSGEARG  "Error parsing arguments"
#define INJECTINSPIRALC_MSGEVAL  "Input argument out of valid range"
#define INJECTINSPIRALC_MSGEFILE "Could not open file"
#define INJECTINSPIRALC_MSGEMEM  "Out of memory"
#define INJECTINSPIRALC_MSGEAPPROXIMANT  "Error parsing order"
#define INJECTINSPIRALC_MSGEORDER  "Error parsing approximant"



#define INJECTINSPIRAL_APPROXIMANT PPN
#define INJECTINSPIRAL_M1          1.4
#define INJECTINSPIRAL_M2         1.4
#define INJECTINSPIRAL_DISTANCE   10  /*mpc*/
#define INJECTINSPIRAL_TC          0.
#define INJECTINSPIRAL_INCLINATION 0.
#define INJECTINSPIRAL_DEC         0.
#define INJECTINSPIRAL_RA          0.
#define INJECTINSPIRAL_PSI         0.
#define INJECTINSPIRAL_PHIC        0.
#define INJECTINSPIRAL_FLOWER     40.
#define INJECTINSPIRAL_FUPPER   ((INJECTINSPIRAL_SAMPLING-1)/2.)   /*must be < sampling/2*/
#define INJECTINSPIRAL_ZETA2       0.
#define INJECTINSPIRAL_OMEGAS      0.
#define INJECTINSPIRAL_THETA       0.
#define INJECTINSPIRAL_T0          0.
#define INJECTINSPIRAL_T1          1.
#define INJECTINSPIRAL_T2          2.
#define INJECTINSPIRAL_SAMPLING 2048.
#define INJECTINSPIRAL_A1          1000.
#define INJECTINSPIRAL_A2          1000.
#define INJECTINSPIRAL_PHI0        0.
#define INJECTINSPIRAL_F0         100.
#define INJECTINSPIRAL_ARG         0
#define INJECTINSPIRAL_UDOT        0.5
#define INJECTINSPIRAL_RP          0.5
#define INJECTINSPIRAL_E           0.   
#define INJECTINSPIRAL_ALPHA       0.
#define INJECTINSPIRAL_FBCV        INJECTINSPIRAL_FUPPER
#define INJECTINSPIRAL_DELTAT      0  /*default use deltaT = 1/sampling*/
#define INJECTINSPIRAL_LENGTH   1024  
#define INJECTINSPIRAL_LENGTHIN  1e9

#define INJECTINSPIRAL_SOURCETHETA 1
#define INJECTINSPIRAL_SOURCEPHI   2
#define INJECTINSPIRAL_SPIN1       0.3
#define INJECTINSPIRAL_SPIN2       0.7
#define INJECTINSPIRAL_THETA1      0.86
#define INJECTINSPIRAL_THETA2      0.86

#define INJECTINSPIRAL_PHI1      1.53
#define INJECTINSPIRAL_PHI2      1.53






void
LALInjectInspiral(
		    LALStatus             *status,
		    CoherentGW            *waveform,
		    SimInspiralTable      *thisEvent,
		    REAL4 *buffer1,
		    REAL8 *buffer2,
		    UINT4 *buffer3
		    )
{

  /*
  REAL8 
    tc       = INJECTINSPIRAL_TC,
    m1       = INJECTINSPIRAL_M1,
    m2       = INJECTINSPIRAL_M2,
    d        = INJECTINSPIRAL_DISTANCE,
    inc      = INJECTINSPIRAL_INCLINATION,
    ra       = INJECTINSPIRAL_RA,
    dec      = INJECTINSPIRAL_DEC,
    psi      = INJECTINSPIRAL_PSI,
    phic     = INJECTINSPIRAL_PHIC,
    fi       = INJECTINSPIRAL_FLOWER,
    ff       = INJECTINSPIRAL_FUPPER,*/
REAL8
    zeta2    = INJECTINSPIRAL_ZETA2,
    omegas   = INJECTINSPIRAL_OMEGAS,
  theta    = INJECTINSPIRAL_THETA;

/*
    t0       = INJECTINSPIRAL_T0,
    t1       = INJECTINSPIRAL_T1,
    t2       = INJECTINSPIRAL_T2,
    a1       = INJECTINSPIRAL_A1,
    a2       = INJECTINSPIRAL_A2,
    phi0     = INJECTINSPIRAL_PHI0,
    f0       = INJECTINSPIRAL_F0,
    arg      = INJECTINSPIRAL_ARG,
    udot     = INJECTINSPIRAL_UDOT,
    rp       = INJECTINSPIRAL_RP,
    e        = INJECTINSPIRAL_E,
    alpha    = INJECTINSPIRAL_ALPHA,
    fbcv     = INJECTINSPIRAL_FBCV,
    sampling = INJECTINSPIRAL_SAMPLING,
    deltaT   = INJECTINSPIRAL_DELTAT, 
    length   = INJECTINSPIRAL_LENGTH;
  
*/


  UINT4 
    order  = -1,
    approximant = -1;
    /*    start=0,
	  tail=0;*/
  
  /*
  UINT4 i,
    nspin=0;
    


  REAL8
    fdata[128];
  */


  PPNParamStruc         ppnParams;
  InspiralTemplate      inspiralParams;
  
  /* System-derived constants. */
  

  
  
  UINT4 n;

  /*  LALAssert(approximant);*/

  INITSTATUS(status, "LALGenerateInspiral",INJECTINSPIRALC);
  ATTATCHSTATUSPTR(status);

  LALGetApproximantAndOrder(status->statusPtr, 
			    thisEvent->waveform, 
			    &order ,
			    &approximant);
  
  CHECKSTATUSPTR(status);
  
  if (order == -1)
    {
      ABORT(status, INJECTINSPIRALC_EORDER, INJECTINSPIRALC_MSGEORDER);
    }
  if (approximant == -1)
    {
      ABORT(status, INJECTINSPIRALC_EAPPROXIMANT, INJECTINSPIRALC_MSGEAPPROXIMANT);
    }
  

  /*  dt = 1./params->inspiral.tSampling;*/
  


  /* fill structures ppn ou inspiral */ 
  switch (approximant)
    {
    case PPN:
      ppnParams.mTot = thisEvent->mass1 + thisEvent->mass2;
      ppnParams.eta  = thisEvent->eta;
      ppnParams.d    = thisEvent->distance * 1.0e6 * LAL_PC_SI;
      ppnParams.inc  = thisEvent->inclination; 
      ppnParams.phi  = thisEvent->coa_phase;   
      ppnParams.fStartIn = INJECTINSPIRAL_FLOWER;
      ppnParams.fStopIn  = -1.0 / 
	(6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot * LAL_MTSUN_SI); 


      ppnParams.position.longitude   = thisEvent->longitude;
      ppnParams.position.latitude    = thisEvent->latitude;  
      ppnParams.position.system      = COORDINATESYSTEM_EQUATORIAL;
      ppnParams.psi                  = thisEvent->polarization; 
      ppnParams.epoch.gpsSeconds     = 0;
      ppnParams.epoch.gpsNanoSeconds = 0;
      ppnParams.ppn = NULL;
      ppnParams.lengthIn = 0;
      ppnParams.deltaT = *buffer1;

      break;

    case EOB:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
      inspiralParams.ieta        = 1;
      inspiralParams.mass1       = thisEvent->mass1;
      inspiralParams.mass2       = thisEvent->mass2;
      inspiralParams.startTime   = 0.0;
      inspiralParams.startPhase  = 0.0;
      inspiralParams.nStartPad   = 0;
      inspiralParams.nEndPad     = 0;
      inspiralParams.signalAmplitude = 1.; 
      inspiralParams.distance    = thisEvent->distance * LAL_PC_SI * 1e6; 
      inspiralParams.fLower      = INJECTINSPIRAL_FLOWER;
      inspiralParams.fCutoff     = INJECTINSPIRAL_FUPPER;
      inspiralParams.order       = order;
      inspiralParams.inclination = thisEvent->inclination ;
      inspiralParams.massChoice=m1Andm2;
      inspiralParams.tSampling = 1./ (*buffer1);
      inspiralParams.OmegaS = omegas;
      inspiralParams.Theta  = theta;
      inspiralParams.Zeta2  = zeta2;
      inspiralParams.approximant  = approximant;
      break;
    case SpinTaylorT3:
      inspiralParams.ieta        = 1;
      inspiralParams.mass1       = thisEvent->mass1;
      inspiralParams.mass2       = thisEvent->mass2;
      inspiralParams.startTime   = 0.0;
      inspiralParams.startPhase  = 0.0;
      inspiralParams.nStartPad   = 0;
      inspiralParams.nEndPad     = 0;
      inspiralParams.signalAmplitude = 1.; 
      inspiralParams.distance    = thisEvent->distance * LAL_PC_SI * 1e6; 
      inspiralParams.fLower      = INJECTINSPIRAL_FLOWER;
      inspiralParams.fCutoff     = INJECTINSPIRAL_FUPPER;
      inspiralParams.order       = order;
      inspiralParams.inclination = thisEvent->inclination ;
      inspiralParams.massChoice  = m1Andm2;
      inspiralParams.tSampling   = 1./ (*buffer1);
      inspiralParams.OmegaS      = omegas;
      inspiralParams.Theta       = theta;
      inspiralParams.Zeta2       = zeta2;
      inspiralParams.approximant = approximant;
      inspiralParams.sourceTheta = INJECTINSPIRAL_SOURCETHETA;
      inspiralParams.sourcePhi   = INJECTINSPIRAL_SOURCEPHI;
      inspiralParams.spin1[0]    = INJECTINSPIRAL_SPIN1; 
      inspiralParams.spin2[0]    = INJECTINSPIRAL_SPIN2;
      inspiralParams.spin1[1]    = INJECTINSPIRAL_THETA1;
      inspiralParams.spin2[1]    = INJECTINSPIRAL_THETA2;
      inspiralParams.spin1[2]    = INJECTINSPIRAL_PHI1;
      inspiralParams.spin2[2]    = INJECTINSPIRAL_PHI2;




      /* Extra computation for spin parameters: */
      ComputeSpin(&inspiralParams);      

      break;
    }


  /* si inspiral on affecte memoire */
  switch  (approximant){
  case EOB:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
    LALInspiralWaveLength(status->statusPtr, &n, inspiralParams);
    CHECKSTATUSPTR(status);      
    LALInspiralParameterCalc(status->statusPtr, &inspiralParams);
    CHECKSTATUSPTR(status);                  
    break;
  }

  /*veritable computation c ici*/
  switch(approximant)
    {
    case TaylorT1: 
      LALInspiralWave1ForInjection(status->statusPtr, 
				   waveform,
				   &inspiralParams);
      *buffer1 = inspiralParams.tSampling;
      *buffer2 = inspiralParams.tC;
      *buffer3 = inspiralParams.nStartPad;
      CHECKSTATUSPTR(status);      
      break;
    case TaylorT2: 
      LALInspiralWave2ForInjection(status->statusPtr, 
				   waveform,
				     &inspiralParams);
      *buffer1 = inspiralParams.tSampling;
      *buffer2 = inspiralParams.tC;
      *buffer3 = inspiralParams.nStartPad;
      CHECKSTATUSPTR(status);      
      break;
    case TaylorT3: 
      LALInspiralWave3ForInjection(status->statusPtr, 
				   waveform,
				     &inspiralParams);
      *buffer1 = inspiralParams.tSampling;
      *buffer2 = inspiralParams.tC;
      *buffer3 = inspiralParams.nStartPad;
      CHECKSTATUSPTR(status);      
	break;
    case TaylorF1:
    case TaylorF2:      
    case PadeT1:
    case PadeF1:
    case BCV:
    case BCVSpin:
    case SpinTaylorT3:
      LALInspiralSpinModulatedWaveForInjection(status->statusPtr, 
					       waveform,
					       &inspiralParams);
      *buffer1 = inspiralParams.tSampling;
      *buffer2 = inspiralParams.tC;
      *buffer3 = inspiralParams.nStartPad;
      CHECKSTATUSPTR(status);      
      break;
    case EOB:       
      LALEOBWaveformForInjection(status->statusPtr, 
					       waveform,
					       &inspiralParams);
      *buffer1 = inspiralParams.tSampling;
      *buffer2 = inspiralParams.tC;
      *buffer3 = inspiralParams.nStartPad;
      CHECKSTATUSPTR(status);      
      break;
    case PPN:
      LALGeneratePPNInspiral(status->statusPtr, waveform, &ppnParams);	    	
      CHECKSTATUSPTR(status);   
      *buffer1 = ppnParams.dfdt;
      *buffer2 = ppnParams.tc;
      *buffer3 = ppnParams.length;
	
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
  char *orderptr = NULL;
  char *approxptr = NULL; 
  
  INITSTATUS( status, "LALGetApproximantAndOrder", INJECTINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( waveform, status, 1,"Null pointer");
  
  if (  (orderptr = strstr(waveform, "newtonian") ) )
    *order = newtonian;
  else if (  (orderptr = strstr(waveform, "oneHalfPN") ) )
    *order = oneHalfPN;
  else if (  (orderptr = strstr(waveform,  "onePN") ) )
    *order = onePN;
  else if (  (orderptr =  strstr(waveform,  "onePointFivePN") ) )
    *order = onePointFivePN;
  else if (  (orderptr = strstr(waveform, "twoPN") ) )
    *order = twoPN;
  else if (  (orderptr = strstr(waveform,  "twoPointFive") ) )
    *order = twoPointFivePN;
  else if (  (orderptr = strstr(waveform, "threePN") ) )
    *order = threePN;
  else if (  (orderptr = strstr(waveform, "threePointFivePN") ) )
    *order = threePointFivePN;


  if (  (approxptr = strstr(waveform, "TaylorT1" ) ) )
    *approximant = TaylorT1;
  else if (  (approxptr = strstr(waveform, "TaylorT2" ) ) )
    *approximant = TaylorT2;
  else if (  (approxptr = strstr(waveform, "SpinTaylorT3" ) ) )
    *approximant = SpinTaylorT3;
  else if (  (approxptr = strstr(waveform, "TaylorT3" ) ) )
    *approximant = TaylorT3;
  else if (  (approxptr = strstr(waveform, "TaylorF1" ) ) )
    *approximant = TaylorF1;
  else if (  (approxptr = strstr(waveform, "TaylorF2" ) ) )
    *approximant = TaylorF2;
  else if (  (approxptr = strstr(waveform, "PadeT1" ) ) )
    *approximant = PadeT1;
  else if (  (approxptr = strstr(waveform, "PadeF1" ) ) )
    *approximant = PadeF1;
  else if (  (approxptr = strstr(waveform, "EOB" ) ) )
    *approximant = EOB;
  else if (  (approxptr = strstr(waveform, "BCV" ) ) )
    *approximant = BCV;
  else if (  (approxptr = strstr(waveform, "GeneratePPN" ) ) )
    *approximant = PPN;
  else if (  (approxptr = strstr(waveform, "TaylorCW" ) ) )
    *approximant = TaylorCW;
  else if (  (approxptr = strstr(waveform, "SpinOrbitCW" ) ) )
    *approximant = SpinOrbitCW;
  

  DETATCHSTATUSPTR( status );
  RETURN( status );
}





void ComputeSpin(InspiralTemplate *params)
{


  double mass1Sq, mass2Sq, spin1Frac, spin2Frac, spin1Phi, spin2Phi, spin1Theta, spin2Theta;

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
