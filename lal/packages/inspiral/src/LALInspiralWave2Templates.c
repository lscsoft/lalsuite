/*  <lalVerbatim file="LALInspiralWave2TemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave2Templates.c}}

This module is exactly the same as \texttt{LALInspiralWave2}
except that it generates two waveforms that differ by a
phase difference of $\pi/2.$

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave2TemplatesCP}
\index{\verb&LALInspiralWave2Templates()&}

\subsubsection*{Description}

\texttt{LALInspiralParameterCalc}\\
\texttt{LALDBisectionFindRoot}\\
\texttt{LALInspiralPhasing2}\\

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave2TemplatesCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

NRCSID (LALINSPIRALWAVE2TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave2TemplatesCP"> */

void LALInspiralWave2Templates(
  LALStatus *status, 
  REAL4Vector *output1, 
  REAL4Vector *output2, 
  InspiralTemplate *params)

{ /* </lalVerbatim>  */

  REAL8 amp, eta, dt, fs, fu, fHigh, phase0, phase1, tC;
  REAL8 phase, v, totalMass, fLso, freq, fOld;
  INT4 i, startShift, count;
  DFindRootIn rootIn;
  InspiralToffInput toffIn;
  void *funcParams;
  expnCoeffs ak;
  expnFunc func;

  INITSTATUS (status, "LALInspiralWave2Templates", LALINSPIRALWAVE2TEMPLATESC);
  ATTATCHSTATUSPTR(status);

  ASSERT(output1,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(output2,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(output1->data,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(output2->data,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(params,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  LALInspiralSetup(status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);

  dt = 1.0/(params->tSampling);   /* sampling interval */
  fs = params->fLower;            /* lower frequency cutoff */
  fu = params->fCutoff;           /* upper frequency cutoff */
  startShift = params->nStartPad; /* number of bins to pad at the beginning */
  phase0 = params->startPhase;    /* initial phasea */
  phase1 = phase0 + LAL_PI_2;

  rootIn.function = func.timing2; /* function to solve for v, given t: 
                                     timing2(v;tC,t)=0.  */

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc(status->statusPtr, params);
  CHECKSTATUSPTR(status);

  ASSERT(params->totalMass > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta <=0.25, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  eta = params->eta;
  totalMass = params->totalMass*LAL_MTSUN_SI; /* solar mass in seconds */

  toffIn.tN = ak.tvaN;
  toffIn.t2 = ak.tva2;
  toffIn.t3 = ak.tva3;
  toffIn.t4 = ak.tva4;
  toffIn.t5 = ak.tva5;
  toffIn.t6 = ak.tva6;
  toffIn.t7 = ak.tva7;
  toffIn.tl6 = ak.tvl6;
  toffIn.piM = ak.totalmass * LAL_PI;

/* Determine the total chirp-time tC: the total chirp time is 
   timing2(v0;tC,t) with t=tc=0*/

  toffIn.t = 0.;
  toffIn.tc = 0.;
  funcParams = (void *) &toffIn;
  func.timing2(status->statusPtr, &tC, fs, funcParams);
  CHECKSTATUSPTR(status);
/* Reset chirp time in toffIn structure */
  toffIn.tc = -tC;

/* Determine the initial phase: it is phasing2(v0) with ak.phiC=0 */
  v = pow(fs * LAL_PI * totalMass, oneby3);
  ak.phiC = 0.0;
  func.phasing2(status->statusPtr, &phase, v, &ak);
  CHECKSTATUSPTR(status);
  ak.phiC = -phase;

/* 
   If flso is less than the user inputted upper frequency cutoff fu, 
*/

  fLso = ak.fn;
  if (fu) 
     fHigh = (fu < fLso) ? fu : fLso; 
  else 
     fHigh = fLso;

/* Is the sampling rate large enough? */

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  rootIn.xmax = 1.1*fu;
  rootIn.xacc = 1.0e-8;
  rootIn.xmin = 0.999999*fs;

  i=0;
  while (i<startShift) output1->data[i] = output2->data[i++] = 0.0;

/* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
  of type void * rather than InspiralToffInput  */

  funcParams = (void *) &toffIn;

  toffIn.t = 0.0;
  freq = fs;
  count=1;
  do
  {
    fOld = freq;
    v = pow(freq*toffIn.piM, oneby3);
    func.phasing2(status->statusPtr, &phase, v, &ak); /* phase at given v */
    CHECKSTATUSPTR(status);
    amp = params->signalAmplitude*v*v;
    output1->data[i]=(REAL4)(amp*cos(phase+phase0));
    output2->data[i]=(REAL4)(amp*cos(phase+phase1));
    i++;
    toffIn.t=count*dt;
    ++count;
/* 
   Determine the frequency at the current time by solving timing2(v;tC,t)=0 
*/
    LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
    CHECKSTATUSPTR(status);
  } while (freq < fHigh && freq > fOld && toffIn.t < -tC);

  while (i<output1->length) output1->data[i]= output2->data[i++]=0.0;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}
