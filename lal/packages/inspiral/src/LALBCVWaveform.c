/**** <lalVerbatim file="LALBCVWaveformCV">
 * Author: B.S. Sathyaprakash
 * $Id$  
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{LALBCVWaveform.c}}
 *
 * This module contains a single function {\tt LALBCVWaveform.}
 * \subsubsection*{Prototypesc}
 * \input{LALBCVWaveformCP}
 * \idx{LALLALBCVWaveform()} 
 * \begin{itemize} 
 * \item {\tt signal:} Output containing the \emph {Fourier transform} of the inspiral waveform.  
 * \item {\tt params:} Input containing binary chirp parameters; 
 * it is necessary and sufficent to specify the following parameters
 * of the {\tt params} structure: 
 * {\tt psi0, psi3, alpha, fendBCV(fFinal), nStartPad, fLower, tSampling}.  
 * All other parameters in {\tt params} are ignored.  \end{itemize} 
 * \input{LALBCVSpinWaveformCP}
 * \idx{LALLALBCVSpinWaveform()} 
 * \begin{itemize} 
 * \item {\tt signal:} Output containing the \emph {Fourier transform} of the inspiral waveform.  
 * \item {\tt params:} Input containing binary chirp parameters; 
 * it is necessary and sufficent to specify the following parameters
 * of the {\tt params} structure: 
 * {\tt psi0, psi3, alpha1, alpha2, beta, fendBCV(fFinal), nStartPad, fLower, tSampling}.  
 * All other parameters in {\tt params} are ignored.  \end{itemize} * 
 * 
 * \subsubsection*{Description}
 * This module can be used to generate {\it detection template 
 * family} of Buonanno, Chen and Vallisneri \cite{BCV03,BCV03b}.
 * There are two modules: {\tt LALBCVWaveform.} and {\tt LALBCVSpinWaveform.}
 * The former can be used to generate non-spinning waveforms and the DTF
 * it implements is given in Sec.~\ref{sec:BCV} and 
 * Eq.(\ref{eq:BCV:NonSpinning}) and the latter
 * to generate spinning waveforms (Eq.~\ref{eq:BCV:Spinning}).
 *
 * \subsubsection*{Algorithm}
 * A straightforward implementation of the formula. Note that the routine returns
 * {\bf Fourier transform} of the signal as opposed to most other modules in this 
 * package which return time-domain signals. Also, the amplitude is quite arbitrary.
 * 
 * \subsubsection*{Uses}
 * \begin{verbatim} 
 * ASSERT 
 * ATTATCHSTATUSPTR
 * DETATCHSTATUSPTR 
 * INITSTATUS 
 * RETURN
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
void 
LALBCVWaveform(
   LALStatus        *status, 
   REAL4Vector      *signal, 
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

  REAL8 f, df;
  REAL8 shift, phi, psi, amp0, amp;
  REAL8 Sevenby6, Fiveby3, Twoby3, alpha;
  INT4 n, i;

  INITSTATUS(status, "LALBCVWaveform", LALBCVWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  ASSERT (signal,  status,       LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params,  status,       LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->fLower > 0, status,     LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->tSampling > 0, status,  LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  n = signal->length;
  Twoby3 = 2.L/3.L;
  Sevenby6 = 7.L/6.L;
  Fiveby3 = 5.L/3.L;

  df = params->tSampling/(REAL8)n;
  params->fFinal = params->fCutoff; 
  alpha = params->alpha / pow(params->fCutoff, Twoby3);


  /* to do : check that a_f in [0,1]??*/

  /*
  totalMass = params->totalMass * LAL_MTSUN_SI;
  shift = LAL_TWOPI * (params->tC - (float)n /params->tSampling 
		  - params->startTime - params->nStartPad/params->tSampling);
  */
  /*
  shift = LAL_TWOPI * ((double) params->nStartPad/params->tSampling + params->startTime);
  phi = params->startPhase + params->psi0 * pow(params->fLower, -Fiveby3) 
	  + params->psi3 * pow(params->fLower, -Twoby3);
   */
  shift = -LAL_TWOPI * ((double)params->nStartPad/params->tSampling);
  phi = - params->startPhase;

  /*
  amp0 = params->signalAmplitude * pow(5./(384.*params->eta), 0.5) * 
	   totalMass * pow(LAL_PI * totalMass,-Sevenby6) * 
	   params->tSampling * (2. / signal->length); 

   */

  amp0 = params->signalAmplitude;
  /*  Computing BCV waveform */

  signal->data[0] = 0.0;
  signal->data[n/2] = 0.0;

  for(i=1;i<n/2;i++)
  {
	  f = i*df;
    
	  if (f < params->fLower || f > params->fFinal)
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




/*  <lalVerbatim file="LALBCVSpinWaveformCP"> */

void 
LALBCVSpinWaveform(
   LALStatus        *status, 
   REAL4Vector      *signal, 
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

  REAL8 f, df;
  REAL8 shift, phi, psi, amp0, amp, beta, modphase;
  REAL8 Sevenby6, Fiveby3, Twoby3, alpha1, alpha2;
  INT4 n, i;

  INITSTATUS(status, "LALBCVSpinWaveform", LALBCVWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  n = signal->length;
  Twoby3 = 2.L/3.L;
  Sevenby6 = 7.L/6.L;
  Fiveby3 = 5.L/3.L;




  df = params->tSampling/(REAL8)n;

 
  alpha1 = params->alpha1;
  alpha2 = params->alpha2;

  /*
  totalMass = params->totalMass * LAL_MTSUN_SI;
  shift = LAL_TWOPI * (params->tC - (float)n /params->tSampling 
		  - params->startTime - params->nStartPad/params->tSampling);
  */
  shift = -LAL_TWOPI * (params->nStartPad/params->tSampling);
      

  phi = - params->startPhase;    
  beta = params->beta; 


  /*
  amp0 = params->signalAmplitude * pow(5./(384.*params->eta), 0.5) * 
	   totalMass * pow(LAL_PI * totalMass,-Sevenby6) * 
	   params->tSampling * (2. / signal->length); 

   */
  amp0 = 1.0L;
  /*  Computing BCV waveform */

  signal->data[0] = 0.0;
  signal->data[n/2] = 0.0;



  for(i=1;i<n/2;i++)
  {
	  f = i*df;
	  if (f < params->fLower || f > params->fCutoff)
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
           
	    psi = (shift*f + phi + params->psi0*pow(f,-Fiveby3) + params->psi3*pow(f,-Twoby3)); 
	    
	    modphase = beta * pow(f,-Twoby3);  
	    
	    amp = amp0 * (1. + (alpha1 * cos(modphase)) + (alpha2 * sin(modphase))) * pow(f,-Sevenby6);

	      	      
              signal->data[i] = (REAL4) (amp * cos(psi));
              signal->data[n-i] = (REAL4) (-amp * sin(psi));
          }
  }

  params->fFinal = params->fCutoff;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
