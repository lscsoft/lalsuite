#include <lal/LALInspiral.h>

NRCSID (LALBCVWAVEFORMC, "$Id$");

void LALBCVWaveform(LALStatus *status, REAL4Vector *signal, InspiralTemplate *params)
{

  REAL8 f, df;
  REAL8 shift, phi, psi, amp0, amp;
  REAL8 Sevenby6, Fiveby3, Twoby3, alpha, totalMass;
  INT4 n, i;

  INITSTATUS(status, "LALBCVWaveform", LALBCVWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  n = signal->length;
  Twoby3 = 2./3.;
  Sevenby6 = 7.0/6.0;
  Fiveby3 = 5./3.;

  df = params->tSampling/(REAL8)n;
  totalMass = params->totalMass * LAL_MTSUN_SI;
  alpha = params->alpha / pow(params->fendBCV, Twoby3);

  shift = LAL_TWOPI * (params->tC - (float)n /params->tSampling 
		  - params->startTime - params->nStartPad/params->tSampling);

  phi = - params->startPhase + LAL_PI/4.;
  amp0 = params->signalAmplitude * pow(5./(384.*params->eta), 0.5) * 
	   totalMass * pow(LAL_PI * totalMass,-Sevenby6) * 
	   params->tSampling * (2. / signal->length); 

  /*  Computing BCV waveform */

  signal->data[0] = 0.0;
  signal->data[n/2] = 0.0;

/*    printf("test\n"); */


  for(i=1;i<n/2;i++)
  {
	  f = i*df;
    
	  if (f < params->fLower || f > params->fendBCV)
	  {
	  /* 
	   * All frequency components below params->fLower and above fn are set to zero  
	   */
              signal->data[i] = 0.;
              signal->data[n-i] = 0.;
            
	  }
       	  else
	  {
          /* What shall we put for sign phi? for uspa it must be "-" */
              psi =  (shift*f + phi + params->psi0*pow(f,-Fiveby3) + params->psi3*pow(f,-Twoby3));
	      amp = amp0 * (1. - alpha * pow(f,Twoby3)) * pow(f,-Sevenby6);
	      
              signal->data[i] = (REAL4) (amp * cos(psi));
              signal->data[n-i] = (REAL4) (-amp * sin(psi));
          }
  }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
