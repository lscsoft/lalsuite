/*  <lalVerbatim file="LALTappRpnTdomFreqTemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomFreqTemplates.c}}

Module which generates two inspiral waveforms which differ in phase by $\pi/2$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTappRpnTdomFreqTemplatesCP}
\index{\verb&LALTappRpnTdomFreqTemplates()&}

\subsubsection*{Description}

Module which generates two inspiral waveforms which differ in phase by $\pi/2$.
This module is identical to \texttt{LALTappRpnTdomFreq} except that is generates two waveforms instead of one.
This is so that when we need to maximise over the phase at time of arrival by seperately performing
correlations of the signal with a zero and $\pi/2$ phase waveform, we do not need to call the waveform
generation routine twice.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc}\\
\texttt{LALDBisectionFindRoot}\\
\texttt{LALTappRpnTdomFreqPhase}

\subsubsection*{Notes}

See the ducumentation for the function \texttt{LALTappRpnTdomFreq} for further details.

\vfill{\footnotesize\input{LALTappRpnTdomFreqTemplatesCV}}

</lalLaTeX>  */




#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>



NRCSID (LALTAPPRPNTDOMFREQC, "$Id$");

/*
void (*LALTappRpnTdomFreqPhase) (LALStatus *status,
                                 REAL8 *phase, 
                                 InspiralPhasesInput *params);

void (*LALTappRpnTdomFreqTofF) (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params);
*/

/*  <lalVerbatim file="LALTappRpnTdomFreqTemplatesCP"> */
void LALTappRpnTdomFreqTemplates(LALStatus *status, 
                                 REAL4Vector *signal1, 
                                 REAL4Vector *signal2, 
                                 InspiralTemplate *params)
{ /* </lalVerbatim>  */

  void (*LALTappRpnTdomFreqPhase) (LALStatus *, REAL8 *, InspiralPhasesInput *);
  void (*LALTappRpnTdomFreqTofF) (LALStatus *, REAL8 *, REAL8, void *);

  REAL8 dt, fs, fu, fsPi, fHigh, phase0;
  REAL8 phase, v, totalMass, fLso, freq;
  INT4 i, startShift, endShift, count;
  InspiralPhasesInput phaseIn;
  DFindRootIn rootIn;
  InspiralToffInput toffIn;
  void *funcParams;


  INITSTATUS (status, "LALTappRpnTdomFreq", LALTAPPRPNTDOMFREQC);
  ATTATCHSTATUSPTR(status);


  ASSERT(signal1,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(signal1->data,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(signal2,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(signal2->data,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(params,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);


  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  dt = 1.0/(params->tSampling);
  fs = params->fLower;
  fu = params->fCutoff;
  fsPi = fs*LAL_PI;
  startShift = params->nStartPad;
  endShift = params->nEndPad;
  phase0 = params->startPhase;



  switch (params->order) {
     case newtonian:
     case oneHalfPN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase0PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF0PN;
          rootIn.function = LALTappRpnTdomFreqTofF0PN;
          break;
     case onePN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase2PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF2PN;
          rootIn.function = LALTappRpnTdomFreqTofF2PN;
          break;
     case onePointFivePN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase3PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF3PN;
          rootIn.function = LALTappRpnTdomFreqTofF3PN;
          break;
     case twoPN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase4PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF4PN;
          rootIn.function = LALTappRpnTdomFreqTofF4PN;
          break;
     default:
          fprintf(stderr, "LALTappRpnTdomFreq: No order selected ... exiting\n");
          exit(0);
     }

/* Calculate the three unknowns from (m1,m2,M,eta,mu) from the two which
   are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);


  ASSERT(params->totalMass > 0.4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->totalMass < 100, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta <=0.25, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->mu >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);



  toffIn.t0 = params->t0;
  toffIn.t2 = params->t2;
  toffIn.t3 = params->t3;
  toffIn.t4 = params->t4;
  toffIn.tc = params->tC;

  phaseIn.p0 = 3.2 * fsPi * params->t0;
  phaseIn.p2 = 4.0 * fsPi * params->t2;
  phaseIn.p3 = 5.0 * fsPi * params->t3;
  phaseIn.p4 = 8.0 * fsPi * params->t4;
  phaseIn.pc = phaseIn.p0 + phaseIn.p2 - phaseIn.p3 + phaseIn.p4;


  totalMass = params->totalMass*LAL_MTSUN_SI;

  fLso = 1.0/(LAL_PI*totalMass*pow(6.0, 1.5));

  if (fu)
  fHigh = (fu < fLso) ? fu/fs : fLso/fs;
  else
  fHigh = fLso/fs;


  ASSERT(fHigh*fs < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh*fs > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  toffIn.t = 0.0;

/* Initialize the members of the structure which will be used as the input to
   the function which will calcute the frequency of the wave  */


  rootIn.xmax = 1.1*fu/fs;
  rootIn.xacc = 1.0e-8;
  rootIn.xmin = 0.999999;


  i=0;
  while (i<startShift) 
	signal1->data[i] = signal2->data[i++] = 0.0;



/* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
  of type void * rather than InspiralToffInput  */

  funcParams = (void *) &toffIn;

  LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
  CHECKSTATUSPTR(status);

  
  phaseIn.f=freq;
  count=1;

  while (freq < fHigh) {
    v = pow(freq*fs*LAL_PI*totalMass, oneby3);
    LALTappRpnTdomFreqPhase(status->statusPtr, &phase, &phaseIn);
    CHECKSTATUSPTR(status);
    signal1->data[i] = (REAL4) params->signalAmplitude * v*v * cos(phase+phase0);
    signal2->data[i++] = (REAL4) params->signalAmplitude * v*v * cos(phase+phase0+LAL_PI_2);
    toffIn.t=count*dt;
    ++count;
    funcParams = (void *) &toffIn;
    LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
    CHECKSTATUSPTR(status);
    phaseIn.f=freq;
  }

  while (i < (int)signal1->length) 
	signal1->data[i]= signal2->data[i++]=0.0;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}
