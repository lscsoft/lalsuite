/*  <lalVerbatim file="LALInspiralWave2CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave2.c} and \texttt{LALInspiralWave2Templates.c}}

These modules generate a time-domain chirp waveform of type {\tt TaylorT2}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave2CP}
\index{\verb&LALInspiralWave2()&}
\begin{itemize}
\item {\tt output:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}
\vspace{0.1in}
\input{LALInspiralWave2TemplatesCP}
\index{\verb&LALInspiralWave2Templates()&}
\begin{itemize}
\item {\tt output1:} Output containing the 0-phase inspiral waveform.
\item {\tt output2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

\texttt{LALInspiralWave2} generates {\tt TaylorT2} approximant wherein 
the phase of the waveform is given as an implicit function of time
as in Equation (\ref{eq:InspiralWavePhase2}). A template is required
to be sampled at equal intervals of time. Thus, first of the equations
in Equation (\ref{eq:InspiralWavePhase2}) is solved for $v$ at equally
spaced values of the time steps
$t_k$ and the resulting value of $v_k$ is used in the second equation to
obtain the phase $\phi_k$. 

\texttt{LALInspiralWave2Templates} is exactly the same as \texttt{LALInspiralWave2}
except that it generates two waveforms that differ in phase by $\pi/2.$

\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc}\\
\texttt{LALDBisectionFindRoot}\\
\texttt{LALInspiralPhasing2}\\

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave2CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

NRCSID (LALINSPIRALWAVE2C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave2CP"> */

void 
LALInspiralWave2(
   LALStatus        *status, 
   REAL4Vector      *output, 
   InspiralTemplate *params
   )
{ /* </lalVerbatim>  */

  REAL8 eta, dt, fs, fu, fHigh, phase0, tC;
  REAL8 phase, v, totalMass, fLso, freq, fOld;
  INT4 i, startShift;
  UINT4 count;
  DFindRootIn rootIn;
  InspiralToffInput toffIn;
  void *funcParams;
  expnCoeffs ak;
  expnFunc func;

  INITSTATUS (status, "LALInspiralWave2", LALINSPIRALWAVE2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(output,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(output->data,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(params,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT((INT4)params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((REAL8)params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((REAL8)params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((INT4)params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((INT4)params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/*  params->nStartPad = 0;*/
  
  LALInspiralSetup(status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);

  dt = 1.0/(params->tSampling);   /* sampling interval */
  fs = params->fLower;            /* lower frequency cutoff */
  fu = params->fCutoff;           /* upper frequency cutoff */
  startShift = params->nStartPad; /* number of bins to pad at the beginning */
  phase0 = params->startPhase;    /* initial phasea */

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
  rootIn.xmin = 0.99999*fs;
  
  i=0;
  while (i<startShift) 
    {
      output->data[i] = 0.0;
      i++;
  }
  
  /* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
     of type void * rather than InspiralToffInput  */
  
  funcParams = (void *) &toffIn;
  fprintf(stderr,"fu=%f fs=%f fHigh=%f\n",fu, fs  ,fHigh); 
  toffIn.t = 0.0;
  freq = fs;
  count=1;
  do
    {
    /*
    Check we're not writing past the end of the vector
    */
    if (i >= output->length)
    {
	ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
    }

    fOld = freq;
    v = pow(freq*toffIn.piM, oneby3);
    func.phasing2(status->statusPtr, &phase, v, &ak); /* phase at given v */
    CHECKSTATUSPTR(status);
    output->data[i++]=(REAL4)(params->signalAmplitude* v*v *cos(phase+phase0));
    toffIn.t=count*dt;
    ++count;
    /* 
       Determine the frequency at the current time by solving timing2(v;tC,t)=0 
    */
    LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
    if (status->statusPtr->statusCode!=0)
    {
	    freq = fHigh+1;
	    status->statusPtr->statusCode = 0;
	    fprintf(stderr,"put a warning here here\n");
	    CHECKSTATUSPTR(status);
    }
    
    } while (freq < fHigh && freq > fOld && toffIn.t < -tC);
  params->fFinal = fOld;
  params->tC = toffIn.t;
  
  while (i<(INT4)output->length) 
    {
      output->data[i]=0.0;
      i++;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

NRCSID (LALINSPIRALWAVE2TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave2TemplatesCP"> */

void 
LALInspiralWave2Templates(
			  LALStatus        *status, 
			  REAL4Vector      *output1, 
			  REAL4Vector      *output2, 
			  InspiralTemplate *params
			  )

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
  ASSERT((INT4)params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((REAL8)params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((REAL8)params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((INT4)params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((INT4)params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  
  LALInspiralSetup(status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);
 
  params->nStartPad = 0;          /* that value must be zero for template generation */ 
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
  while (i<startShift) 
    {
      output1->data[i] = output2->data[i] = 0.0;
      i++;
  }
  
  /* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
     of type void * rather than InspiralToffInput  */
  
  funcParams = (void *) &toffIn;
  
  toffIn.t = 0.0;
  freq = fs;
  count=1;
  do
    {
    /*
    Check we're not writing past the end of the vector
    */
    if (i >= output1->length)
    {
        ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
    }

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
  params->fFinal = fOld;

  while (i<(INT4)output1->length) 
    {
      output1->data[i]= output2->data[i]=0.0;
      i++;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}




NRCSID (LALINSPIRALWAVE2FORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave2ForInjectionCP"> */

void 
LALInspiralWave2ForInjection(
   LALStatus        *status, 
   CoherentGW *waveform,   
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams			     
   )

{ /* </lalVerbatim>  */

  REAL8  eta, dt, fs, fu, fHigh, phase0, phase1, tC, omega;
  REAL8 phase, v, totalMass, fLso, freq, fOld;
  INT4  startShift;
  UINT4 count,i;

  REAL4Vector *a   = NULL;/* pointers to generated amplitude  data */
  REAL4Vector *ff  = NULL ;/* pointers to generated  frequency data */
  REAL8Vector *phi = NULL;/* pointer to generated phase data */

  CreateVectorSequenceIn in;
  DFindRootIn rootIn;
  InspiralToffInput toffIn;
  void *funcParams;
  expnCoeffs ak;
  expnFunc func;
 
  REAL8 unitHz;
  REAL8 f2a;
  REAL8 mu; 
  REAL8 mTot;
  REAL8 cosI;/* cosine of system inclination */
  REAL8 etab;
  REAL8 fFac; /* SI normalization for f and t */
  REAL8 f2aFac;/* factor multiplying f in amplitude function */
  REAL8 apFac, acFac;/* extra factor in plus and cross amplitudes */
  REAL8 phiC;/* phase at coalescence */

  CHAR message[256];
  
  InspiralInit paramsInit;  

  /** -- -- */
  INITSTATUS (status, "LALInspiralWave2ForInjection", LALINSPIRALWAVE2FORINJECTIONC);
  ATTATCHSTATUSPTR(status);

  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);  
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
   
  
  params->nStartPad=0;
  
  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);   
  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
      
    }
  func = paramsInit.func;
  ak = paramsInit.ak;
  
  fHigh = params->fCutoff;
  dt = 1./params->tSampling;

  /** -- some aliases -- */
  mTot   =  params->mass1 + params->mass2;
  etab   =  params->mass1 * params->mass2;
  etab  /= mTot;
  etab  /= mTot;
  unitHz = (mTot) *LAL_MTSUN_SI*(REAL8)LAL_PI;
  cosI   = cos( params->inclination );
  mu     = etab * mTot;  
  fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
  f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;   
  apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
  apFac *= 1.0 + cosI*cosI;
  acFac *= 2.0 * cosI;

  /* Now we can allocate memory and vector for coherentGW structure*/     
  LALSCreateVector(status->statusPtr, &ff, paramsInit.nbins);
  CHECKSTATUSPTR(status);   
  LALSCreateVector(status->statusPtr, &a, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);   
  LALDCreateVector(status->statusPtr, &phi, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  
  /* By default the waveform is empty */
  for (count = 0; count < paramsInit.nbins; count++) 
    {
      ff->data[count]           = 0.;
      a->data[2*count+1]        = 0.;
      phi->data[count]         = 0.;
      a->data[2*count]          = 0.;
    }
  count = 0;


  
  fs = params->fLower;            /* lower frequency cutoff */
  fu = params->fCutoff;           /* upper frequency cutoff */
  startShift = params->nStartPad; /* number of bins to pad at the beginning */
  phase0 = params->startPhase;    /* initial phasea */
  phase1 = phase0 + LAL_PI_2;
  
  rootIn.function = func.timing2; /* function to solve for v, given t: 
                                     timing2(v;tC,t)=0.  */
  
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
  ak.phiC = LAL_PI/2;
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
  
  i = 0;
  
  /* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
     of type void * rather than InspiralToffInput  */

  funcParams = (void *) &toffIn;
  
  toffIn.t = 0.0;
  freq     = fs;
  count    = 0;
  do
    {
    /* If attempting to write beyond the end of vectors,
       free up memory and abort.
    */
    if (count >= ff->length)
    {
	LALSDestroyVector(status->statusPtr, &ff);
	CHECKSTATUSPTR(status);
	LALSDestroyVector(status->statusPtr, &a);
	CHECKSTATUSPTR(status);
	LALDDestroyVector(status->statusPtr, &phi);
	CHECKSTATUSPTR(status);
	ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
    }

    fOld = freq;
    v = pow(freq*toffIn.piM, oneby3);
    func.phasing2(status->statusPtr, &phase, v, &ak); /* phase at given v */
    CHECKSTATUSPTR(status);

    omega = v*v*v;

    ff->data[count]= (REAL4)(omega/unitHz);
    f2a = pow (f2aFac * omega, 2./3.);
    a->data[2*count]          = (REAL4)(4.*apFac * f2a);
    a->data[2*count+1]        = (REAL4)(4.*acFac * f2a);
    phi->data[count]          = (REAL8)(phase);
    
    ++count; /*has to be before toffIn */
    toffIn.t = count * dt;
    
    /* 
       Determine the frequency at the current time by solving timing2(v;tC,t)=0 
    */
    LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
    CHECKSTATUSPTR(status);
    } while (freq < fHigh && freq > fOld && toffIn.t < -tC);
  params->fFinal = fOld;

  
  if ( fabs(phi->data[count-1]/2.)/LAL_PI < 2. ){
        sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.", 
	       (double)(fabs(phi->data[count-1]/2.)/LAL_PI) );
    LALWarning(status, message);


  }
  else
    {
      /*wrap the phase vector*/
      phiC =  phi->data[count-1] ;
      for (i=0; i<count;i++)
	{
	  phi->data[i] =  phi->data[i] -phiC + ppnParams->phi;
	}
      /* Allocate the waveform structures. */
      if ( ( waveform->a = (REAL4TimeVectorSeries *)
	     LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );
      if ( ( waveform->f = (REAL4TimeSeries *)
	     LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->f, 0, sizeof(REAL4TimeSeries) );
      if ( ( waveform->phi = (REAL8TimeSeries *)
	     LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	LALFree( waveform->f ); waveform->f = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );
      
      
      in.length = (UINT4)count;
      in.vectorLength = 2;
      LALSCreateVectorSequence( status->statusPtr,
				&( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);      
      LALSCreateVector( status->statusPtr,
			&( waveform->f->data ), count);
      CHECKSTATUSPTR(status);      
      LALDCreateVector( status->statusPtr,
			&( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);        
      
      memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));
      
      
      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
	= 1./params->tSampling;
      
      waveform->a->sampleUnits   = lalStrainUnit;
      waveform->f->sampleUnits   = lalHertzUnit;
      waveform->phi->sampleUnits = lalDimensionlessUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      LALSnprintf( waveform->a->name, LALNameLength,   "T2 inspiral amplitude" );
      LALSnprintf( waveform->f->name, LALNameLength,   "T2 inspiral frequency" );
      LALSnprintf( waveform->phi->name, LALNameLength, "T2 inspiral phase" );
      
      
      /* --- fill some output ---*/
      ppnParams->tc     = (double)(count-1) / params->tSampling ;
      ppnParams->length = count;
      ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1] 
				   - waveform->f->data->data[count-2]))
	* ppnParams->deltaT;
      ppnParams->fStop  = params->fFinal;
      ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
      ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
      
      ppnParams->fStart   = ppnParams->fStartIn;
    }/*end of coherentGW storage */

  
  /* --- free memory --- */
  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phi);
  CHECKSTATUSPTR(status);
  

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
