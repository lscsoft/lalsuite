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
 * {\tt psi0, psi3, alpha, fendBCV, nStartPad, fLower, tSampling}.  
 * All other parameters in {\tt params} are ignored.  \end{itemize} * 
 * \subsubsection*{Description}
 * 
 * This module can be used to generate {\it detection template family}
 * of Buonanno, Chen and Vallisneri \cite{BCV03}.
 * The only function needed is {\tt LALBCVWaveform.}
 *
 * The Fourier transform of a chirp waveform in the restricted post-Newtonian
 * approximation in * the stationary phase approximation is given, for 
 * positive frequencies $f,$ by
 * \begin{equation}
 * \tilde h(f) = h_0 f^{-7/6} \exp \left [ \sum_k \psi_k f^{(k-5)/3} \right ],
 * \end{equation}
 * where $h_0$ is a constant for a given system and $psi_k$ are parameters that
 * depend on the two masses of the binary. Since the time-domain waveform 
 * is terminated at when the instantaneous GW frequency reaches a certain value
 * $F_{\rm cut}$ (which is either the last stable
 * orbit or the light-ring defined by the model) and since the contribution to
 * a Fourier component comes mainly from times when the GW instantaneous frequency
 * reaches that value, it is customery to terminate the Fourier transform at the
 * same frequency, namely $f_{\rm cut} = F_{\rm cut}.$ In otherwords, the Fourier
 * transform is taken to be
 * \begin{equation}
 * \tilde h(f) = h_0 f^{-7/6} \theta(f-f_{\rm cut}) \exp \left [ \sum_k \psi_k f^{(k-5)/3} \right ],
 * \end{equation}
 * where $\theta(x<0)=0$ and $\theta(x\ge 0) =1.$ 
 *
 * There are different post-Newtonian models such as the standard post-Newtonian,
 * P-approximants, effective one-body, and their overlaps with one another is not
 * as good as we would like them to be. The main reason for this is that matched 
 * filtering is sensitive to the phasing of the waves. It is not clear which model
 * best describes the true GW signal from a compact binary inspiral although some,
 * like the EOB, are more robust in their predictions than others. Thus, Buonanno,
 * Chen and Vallisneri proposed a {\it phenomenological} model as a detection template
 * family (DTF) based on the above expression for the Fourier transform. Indeed, they
 * proposed to use a DTF that depends on four parameters $(\psi_0, \psi_3, f_{\rm cut}, \alpha)$.
 * The proposed waveform has structure similar to the one above: 
 * \begin{equation}
 * \tilde h(f) = h_0 f^{-7/6} \left (1 - \alpha f^{2/3} \right) \theta(f-f_{\rm cut}) 
 * \exp \left [ \psi_0 f^{-5/3} + \psi_3 f^{-2/3} \right ],
 * \end{equation}
 * where the motivation to include an amplitude correction term $\alpha$ is based on the
 * fact that the first post-Newtonian correction to the amplitude would induce a term like this.
 * Note carefully that the phasing does not include the full post-Newtonian expansion 
 * but only the Newtonian and 1.5 post-Newtonian terms. It turns out this four-parameter
 * family of waveforms has good overlap with the two-parameter family of templates
 * corresponding to different post-Newtonian models and their improvments.
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
  Twoby3 = 2.L/3.L;
  Sevenby6 = 7.L/6.L;
  Fiveby3 = 5.L/3.L;

  df = params->tSampling/(REAL8)n;
  alpha = params->alpha / pow(params->fendBCV, Twoby3);

  /*
  totalMass = params->totalMass * LAL_MTSUN_SI;
  shift = LAL_TWOPI * (params->tC - (float)n /params->tSampling 
		  - params->startTime - params->nStartPad/params->tSampling);
  */
  shift = -LAL_TWOPI * (params->nStartPad/params->tSampling);
  phi = - params->startPhase + LAL_PI/4.;
  /*
  amp0 = params->signalAmplitude * pow(5./(384.*params->eta), 0.5) * 
	   totalMass * pow(LAL_PI * totalMass,-Sevenby6) * 
	   params->tSampling * (2. / signal->length); 

   */
  amp0 = 1.0L;
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
