/*  <lalVerbatim file="LALInspiralWave3CV">

Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave3.c} and \texttt{LALInspiralWave3Templates.c}}
These modules generate a time-domain chirp waveform of type {\tt TaylorT3}.


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave3CP}
\index{\verb&LALInspiralWave3()&}
\begin{itemize}
\item {\tt output:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}
\vspace{0.1in}
\input{LALInspiralWave3TemplatesCP}
\index{\verb&LALInspiralWave3Templates()&}
\begin{itemize}
\item {\tt output1:} Output containing the 0-phase inspiral waveform.
\item {\tt output2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}


\subsubsection*{Description}
{\tt LALInspiralWave3} generates {\tt TaylorT3} approximant which
corresponds to the case wherein 
the phase of the waveform is given as an explicit function of time
as in Equation (\ref{eq:InspiralWavePhase3}). 

{\tt LALInspiralWave3Templates} simultaneously generates 
two inspiral waveforms and the two differ in 
phase by $\pi/2$.  
 
\subsubsection*{Algorithm}

\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc} \\
\texttt{LALInspiralChooseModel} \\
\texttt{LALInspiralSetup} \\
\texttt{LALInspiralPhasing3} (via expnFunc)\\
\texttt{LALInspiralFrequency3}. (via expnFunc)\\


\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave3CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>


typedef struct
{
	void (*func)(LALStatus *status, REAL8 *f, REAL8 tC, expnCoeffs *ak);
	expnCoeffs ak;
} 
ChirptimeFromFreqIn;

static void LALInspiralFrequency3Wrapper(LALStatus *status, REAL8 *f, REAL8 tC, void *pars);

NRCSID (LALINSPIRALWAVE3C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3CP"> */

void 
LALInspiralWave3 (
   LALStatus        *status,
   REAL4Vector      *output, 
   InspiralTemplate *params
   )

{ /* </lalVerbatim>  */

  INT4 i, startShift, count;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase;
  DFindRootIn rootIn;
  expnFunc func;
  expnCoeffs ak;
  ChirptimeFromFreqIn timeIn;
  void *pars;


  INITSTATUS (status, "LALInspiralWave3", LALINSPIRALWAVE3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  params->nStartPad = 0; 

  LALInspiralSetup (status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  ASSERT(params->totalMass >0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  eta = params->eta;   /* Symmetric mass ratio  */
  ieta = params->ieta;
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */
  /* constant that appears in the definition of the time parameter Theta */
  c1 = eta/(5.*totalMass);
 
  /*
   * In Jan 2003 we realized that the tC determined as a sum of chirp times is
   * not quite the tC that should enter the definition of Theta in the expression
   * for the frequency as a function of time (see DIS3, 2000). This is because
   * chirp times are obtained by inverting t(f). Rather tC should be obtained by
   * solving the equation f0 - f(tC) = 0. This is what is implemented below.
   */

  timeIn.func = func.frequency3;
  timeIn.ak = ak;
  rootIn.function = &LALInspiralFrequency3Wrapper;
  rootIn.xmin = c1*params->tC/2.;
  rootIn.xmax = c1*params->tC*2.;
  rootIn.xacc = 1.e-6;
  pars = (void*) &timeIn;
  /* tc is the instant of coalescence */
  LALDBisectionFindRoot (status->statusPtr, &tc, &rootIn, pars);
  CHECKSTATUSPTR(status);

  tc /= c1;

  tc += params->startTime;       /* Add user given startTime to instant of 
				     coalescence of the compact objects */

/* 
   If flso is less than the user inputted upper frequency cutoff fu, 
   then the waveforn is truncated at f=flso.  
*/

  if (fu) 
     fHigh = (fu < ak.flso) ? fu : ak.flso; 
  else 
     fHigh = ak.flso;

/* 
   Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum 
*/

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Here's the part which calculates the waveform */

  i=0; while (i<startShift) 
  {
      output->data[i] = 0.0;
      i++;
  }

  t=0.0;
  td = c1*(tc-t);
  func.phasing3(status->statusPtr, &phase, td, &ak);
  CHECKSTATUSPTR(status);
  func.frequency3(status->statusPtr, &f, td, &ak);
  CHECKSTATUSPTR(status);
  phi0=-phase+phi;

  /*
  fprintf(stderr, "Starting frequency=%e\n", f);
   */

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */

  while (f<fHigh && t<tmax && f>fOld) 
  {
    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = v*v; 
    output->data[i++] = (REAL4) (params->signalAmplitude * amp * cos(phase+phi0));
    ++count;
    t=count*dt;
    td = c1*(tc-t);
    func.phasing3(status->statusPtr, &phase, td, &ak);
    CHECKSTATUSPTR(status);
  
    func.frequency3(status->statusPtr, &f, td, &ak);
    CHECKSTATUSPTR(status); 
  }
  params->fFinal = fOld;
  params->tC = t;
  
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/
  while (i < (int)output->length) 
  {
      output->data[i]=0.0;
      i++;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

static void LALInspiralFrequency3Wrapper(LALStatus *status, REAL8 *f, REAL8 tC, void *pars)
{

  ChirptimeFromFreqIn *in;
  REAL8 freq;
  INITSTATUS (status, "LALInspiralWave3", LALINSPIRALWAVE3C);
  ATTATCHSTATUSPTR(status);

  in = (ChirptimeFromFreqIn *) pars;
  in->func(status->statusPtr, &freq, tC, &(in->ak));
  *f = freq - in->ak.f0;

  /*
  fprintf(stderr, "Here freq=%e f=%e tc=%e f0=%e\n", freq, *f, tC, in->ak.f0);
   */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

NRCSID (LALINSPIRALWAVE3TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3TemplatesCP"> */

void 
LALInspiralWave3Templates (
   LALStatus        *status,
   REAL4Vector      *output1, 
   REAL4Vector      *output2, 
   InspiralTemplate *params
   )

{ /* </lalVerbatim>  */

  INT4 i, startShift, count;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi1, phi;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase;
  expnFunc func;
  expnCoeffs ak;


  INITSTATUS (status, "LALInspiralWave3Templates", LALINSPIRALWAVE3TEMPLATESC);
  ATTATCHSTATUSPTR(status);

  ASSERT(output1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  LALInspiralSetup (status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  ASSERT(params->totalMass >0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  tc=params->tC;       /* Instant of coalescence of the compact objects */
  eta = params->eta;   /* Symmetric mass ratio  */
  ieta = params->ieta;
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */


/* 
   If flso is less than the user inputted upper frequency cutoff fu, 
   then the waveforn is truncated at f=flso.  
*/

  if (fu) 
     fHigh = (fu < ak.flso) ? fu : ak.flso; 
  else 
     fHigh = ak.flso;

/* 
   Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum 
*/

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);
  i=0; while (i<startShift) 
  {
      output1->data[i] = output2->data[i] = 0.0;
      i++;
  }

  t=0.0;
  td = c1*(tc-t);
  func.phasing3(status->statusPtr, &phase, td, &ak);
  CHECKSTATUSPTR(status);
  func.frequency3(status->statusPtr, &f, td, &ak);
  CHECKSTATUSPTR(status);
  phi0=-phase+phi;
  phi1=phi0+LAL_PI_2;

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */

  while (f<fHigh && t<tmax && f>fOld) 
  {
    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = params->signalAmplitude * v*v; 
    output1->data[i] = (REAL4) amp * cos(phase+phi0);
    output2->data[i] = (REAL4) amp * cos(phase+phi1);
    ++i;
    ++count;
    t=count*dt;
    td = c1*(tc-t);
    func.phasing3(status->statusPtr, &phase, td, &ak);
    CHECKSTATUSPTR(status);
    func.frequency3(status->statusPtr, &f, td, &ak);
    CHECKSTATUSPTR(status); 
  }
  params->fFinal = fOld;
  
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/
  while (i < (int)output1->length) 
  {
      output1->data[i]=output2->data[i]=0.0;
      i++;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}





NRCSID (LALINSPIRALWAVE3FORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3ForInjectionCP"> */

void 
LALInspiralWave3ForInjection (
   LALStatus        *status,
   CoherentGW       *waveform,   
   InspiralTemplate *params
   )

{ /* </lalVerbatim>  */

  INT4 i, startShift, count, length;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi1, phi,omega;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase;
  expnFunc func;
  expnCoeffs ak;
  REAL4Vector *a=NULL;
  REAL4Vector *ff=NULL ;
  REAL8Vector *phiv=NULL;
  CreateVectorSequenceIn in;
  REAL8 unitHz,f2a, mu, mTot, cosI, etab, fFac,  f2aFac, apFac, acFac, phiC;
    
  mTot   =  params->mass1 + params->mass2;
  etab   =  params->mass1 * params->mass2;
  etab  /= mTot;
  etab  /= mTot;
  unitHz = (mTot) *LAL_MTSUN_SI*(REAL8)LAL_PI;
  cosI   = cos( params->inclination );
  mu     = etab * mTot;  
  fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
  dt     = -1. * etab / ( params->tSampling * 5.0*LAL_MTSUN_SI*mTot );      
  f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;   
  apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
  apFac *= 1.0 + cosI*cosI;
  acFac *= 2.0 * cosI;

  


  INITSTATUS (status, "LALInspiralWave3ForInjection", LALINSPIRALWAVE3FORINJECTIONC);
  ATTATCHSTATUSPTR(status);
  

/* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL,
	  LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);  
  /* Make sure waveform fields don't exist. */
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL,
	  LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL,
	  LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL,
	  LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status, LALINSPIRALH_ENULL,
	  LALINSPIRALH_MSGENULL );
   
  

  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  params->nStartPad = 0;

  /*Compute some parameters*/
   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralWaveLength(status->statusPtr, &length, *params);
   CHECKSTATUSPTR(status);

  
   /*Now we can allocate memory and vector for coherentGW structure*/     
   LALSCreateVector(status->statusPtr, &ff, length);
   CHECKSTATUSPTR(status);   
   LALSCreateVector(status->statusPtr, &a, 2*length);
   CHECKSTATUSPTR(status);   
   LALDCreateVector(status->statusPtr, &phiv, length);
   CHECKSTATUSPTR(status);
  

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  ASSERT(params->totalMass >0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  tc=params->tC;       /* Instant of coalescence of the compact objects */
  eta = params->eta;   /* Symmetric mass ratio  */
  ieta = params->ieta;
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */


/* 
   If flso is less than the user inputted upper frequency cutoff fu, 
   then the waveforn is truncated at f=flso.  
*/

  if (fu) 
     fHigh = (fu < ak.flso) ? fu : ak.flso; 
  else 
     fHigh = ak.flso;

/* 
   Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum 
*/

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);
  i=0;


  t=0.0;
  td = c1*(tc-t);
  func.phasing3(status->statusPtr, &phase, td, &ak);
  CHECKSTATUSPTR(status);
  func.frequency3(status->statusPtr, &f, td, &ak);
  CHECKSTATUSPTR(status);
  phi0=-phase+phi;
  phi1=phi0+LAL_PI_2;

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */
  while (f<fHigh && t<tmax && f>fOld) 
  {
    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = params->signalAmplitude * v*v;

    omega = v*v*v;

    ff->data[count]= (REAL4)(omega/unitHz);
    f2a = pow (f2aFac * omega, 2./3.);
    a->data[2*count]          = (REAL4)(4.*apFac * f2a);
    a->data[2*count+1]        = (REAL4)(4.*acFac * f2a);
    phiv->data[count]          = (REAL8)(phase);


    ++i;
    ++count;
    t=count*dt;
    td = c1*(tc-t);
    func.phasing3(status->statusPtr, &phase, td, &ak);
    CHECKSTATUSPTR(status);
    func.frequency3(status->statusPtr, &f, td, &ak);
    CHECKSTATUSPTR(status); 
  }
  params->fFinal = fOld;
  
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/
/*wrap the phase vector*/
   phiC =  phiv->data[count-1] ;
   for (i=0; i<count;i++)
     {
       phiv->data[i] =  phiv->data[i] -phiC;
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
  memcpy(waveform->phi->data->data ,phiv->data, count*(sizeof(REAL8)));


 
  dt = -1. * etab / ( params->tSampling * 5.0*LAL_MTSUN_SI*mTot );      

  waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
    = 1./params->tSampling;

  waveform->a->sampleUnits = lalStrainUnit;
  waveform->f->sampleUnits = lalHertzUnit;
  waveform->phi->sampleUnits = lalDimensionlessUnit;
 LALSnprintf( waveform->a->name, LALNameLength, "T3 inspiral amplitudes" );
  LALSnprintf( waveform->f->name, LALNameLength, "T3 inspiral frequency" );
  LALSnprintf( waveform->phi->name, LALNameLength, "T3  inspiral phase" );


  params->tC = count / params->tSampling ;
  params->nStartPad = count;

  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phiv);
  CHECKSTATUSPTR(status);



  DETATCHSTATUSPTR(status);
  RETURN(status);
}
