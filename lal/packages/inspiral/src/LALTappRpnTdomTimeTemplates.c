/*  <lalVerbatim file="LALTappRpnTdomTimeTemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomTimeTemplates.c}}

Module to generate two inspiral waveforms simultaneously, which differ in phase by $\pi/2$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTappRpnTdomTimeTemplatesCP}
\index{\verb&LALTappRpnTdomTimeTemplates()&}

\subsubsection*{Description}

Module to generate two inspiral waveforms simultaneously, which differ in 
phase by $\pi/2$.  This is so that when we need to maximise over the phase 
at time of arrival by seperately correlating a signal with a zero and $\pi/2$ 
phase waveform, we do not need to call the waveform generation function
twice.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc} \\
\texttt{LALTappRpnTdomTimePhase} \\
\texttt{LALTappRpnTdomTimeFrequency}. \\

\subsubsection*{Notes}

See the documentation for the module \texttt{LALTappRpnTdomTime} for further details.

\vfill{\footnotesize\input{LALTappRpnTdomTimeTemplatesCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALTAPPRPNTDOMTIMEC, "$Id$");


/*  <lalVerbatim file="LALTappRpnTdomTimeTemplatesCP"> */
void LALTappRpnTdomTimeTemplates (LALStatus *status,
                                  REAL4Vector *signal1, 
                                  REAL4Vector *signal2, 
                                  InspiralTemplate *params)
{ /* </lalVerbatim>  */
  void (*LALTappRpnTdomTimeFrequency) (LALStatus *,
      InspiralwaveFrequencyOutput *, InspiralwaveFrequencyInput *);

  void (*LALTappRpnTdomTimePhase) (LALStatus *, InspiralwavePhaseOutput *,
      InspiralwavePhaseInput *);

  INT4 i, startShift, endShift, count;
  REAL8 dt, fu, eta, tc, totalMass, t, c1, phi0, phi;
  REAL8 v, fLso, fHigh, amp, tmax;


  InspiralwavePhaseInput input1;       /*Input structure for LALTappRpnTdomTimePhase()*/
  InspiralwavePhaseOutput output1;     /*Output structure for LALTappRpnTdomTimePhase()*/
  InspiralwaveFrequencyInput input2;   /*Input structure for LALTappRpnTdomTimeFrequency()*/
  InspiralwaveFrequencyOutput output2; /*Output structure for LALTappRpnTdomTimeFrequency()*/ 


  INITSTATUS (status, "LALTappRpnTdomTime", LALTAPPRPNTDOMTIMEC);
  ATTATCHSTATUSPTR(status);

  ASSERT(signal1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signal1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signal2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signal2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  dt = 1.0/(params->tSampling);          /* The sampling rate  */
  fu = params->fCutoff;            /* The upper frequency cutoff  */
  phi = params->startPhase;        /* The initial phase  */
  startShift = params->nStartPad;  /* The number of zeros at the start of the wave  */
  endShift = params->nEndPad;      /* The number of zeros at the end of the wave  */

/* Calculate the three unknowns in (m1,m2,M,eta,mu) from the two which are given. */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);


  ASSERT(params->totalMass > 0.4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->totalMass < 100, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->mu >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  ASSERT(params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->order >= 4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  switch (params->order) {
     case newtonian:
     case oneHalfPN:
          LALTappRpnTdomTimePhase = &LALTappRpnTdomTimePhase0PN;
          LALTappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency0PN;
          break;
     case onePN:
          LALTappRpnTdomTimePhase = &LALTappRpnTdomTimePhase2PN;
          LALTappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency2PN;
          break;
     case onePointFivePN:
          LALTappRpnTdomTimePhase = &LALTappRpnTdomTimePhase3PN;
          LALTappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency3PN;
          break;
     case twoPN:
          LALTappRpnTdomTimePhase = &LALTappRpnTdomTimePhase4PN;
          LALTappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency4PN;
          break;
     }
  
 
  tc=params->tC;              /* Instant of coalescence of the compact objects */
  eta = params->eta;                              /* Symmetric mass ratio  */
  totalMass = (params->totalMass)*LAL_MTSUN_SI;   /* The mass of the system in seconds */



/* Calculate the frequency of the last stable orbit flso. If flso is less 
   than the user inputted upper frequency cutoff fu, then the waveforn is 
   truncated at f=flso.  If fu is less than flso, then we truncate the 
   waveform at fu. */

  fLso = 1.0/(LAL_PI*totalMass*pow(6.0, 1.5));

  if (fu)
  fHigh = (fu < fLso) ? fu : fLso;
  else
  fHigh = fLso;

/* Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum */

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Initialize members of the structures which get fed into LALTappRpnTdomTimePhase()
   and LALTappRpnTdomTimeFrequency(). */

  input1.etaby2 = 0.5 * eta;
  input1.a2 = 3715./8064. + 55.*eta/96.;
  input1.a3 = 0.75*LAL_PI;
  input1.a4 = 9275495./14450688. + 284875.*eta/258048. + 1855.*pow(eta,2.0)/2048.;

  input2.eightPiM = 8.*LAL_PI*totalMass;
  input2.a2 = 743./2688.+(11.*eta)/32.;
  input2.a3 = 0.3*LAL_PI;
  input2.a4 =  1855099./14450688. +  56975.*eta/258048. + 371.*pow(eta,2.0)/2048;

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);
  i=0; 
  while (i<startShift) 
	 signal1->data[i] = signal2->data[i++] = 0.0;

  t=0.0;
  input1.td = input2.td = c1*(tc-t);
  LALTappRpnTdomTimePhase(status->statusPtr, &output1, &input1);
  CHECKSTATUSPTR(status);
  LALTappRpnTdomTimeFrequency(status->statusPtr, &output2, &input2);
  CHECKSTATUSPTR(status);
  phi0=-output1.phase+phi;

  count = 0;
  tmax = tc - dt;
  while (output2.frequency<fHigh && t<tmax) 
  {
    v = pow(output2.frequency*LAL_PI*totalMass, oneby3);
    amp = v*v; 
    signal1->data[i] = (REAL4) params->signalAmplitude * amp * cos(output1.phase+phi0);
    signal2->data[i++] = (REAL4) params->signalAmplitude * amp * cos(output1.phase+phi0+LAL_PI_2);
    ++count;
    t=count*dt;
    input1.td = input2.td = c1*(tc-t);
    LALTappRpnTdomTimeFrequency(status->statusPtr, &output2, &input2);
    CHECKSTATUSPTR(status); 
    LALTappRpnTdomTimePhase(status->statusPtr, &output1, &input1);
    CHECKSTATUSPTR(status);
  }

  while (i < (int)signal1->length) 
	signal1->data[i]= signal2->data[i++]=0.0;


  DETATCHSTATUSPTR(status);
  RETURN(status);
}
