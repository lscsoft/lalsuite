/**** <lalVerbatim file="LALBCVWaveformCV">
 * Author: B.S. Sathyaprakash
 * $Id$  
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{LALBCVWaveform.c}}
 *
 * This module contains a single function {\it LALBCVWaveform.}
 * \subsubsection*{Prototypesc}
 * \input{LALBCVWaveformCP}
 * \idx{LALLALBCVWaveform()}
 * 
 * \subsubsection*{Description}
 * 
 * This module can be used to generate {\it detection template family}
 * of Buonanno, Chen and Vallisneri (Phys. Rev. D Phys.Rev. D67 (2003) 024016). 
 * The only function needed is {\it LALBCVWaveform.}
 * 
 * \subsubsection*{Algorithm}
 * 
 * %% A description of the method used to perform the calculation.
 * 
 * \subsubsection*{Uses}
 *  None.
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * \end{verbatim}
 * 
 * \subsubsection*{Notes}
 * 
 * %% Any relevant notes.
 * 
 * \vfill{\footnotesize\input{LALBCVWaveformCV}}
 * 
 **** </lalLaTeX> */ 

#include <lal/LALInspiral.h>

/** DEFINE RCS ID STRING **/
NRCSID (LALBCVWAVEFORMC, "$Id$");


/*  <lalVerbatim file="LALBCVWaveformCP"> */
void LALBCVWaveform(
		LALStatus *status, 
		REAL4Vector *signal, 
		InspiralTemplate *params)
 { /* </lalVerbatim>  */

  REAL8 f, df;
  REAL8 shift, phi, psi, amp0, amp;
  REAL8 Sevenby6, Fiveby3, Twoby3, alpha, totalMass;
  INT4 n, i;

  INITSTATUS(status, "LALBCVWaveform", LALBCVWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

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
