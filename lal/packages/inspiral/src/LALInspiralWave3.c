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

static void
LALInspiralWave3Engine(
                LALStatus        *status,
                REAL4Vector      *output1,
                REAL4Vector      *output2,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                );

NRCSID (LALINSPIRALWAVE3C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3CP"> */

void 
LALInspiralWave3 (
   LALStatus        *status,
   REAL4Vector      *output, 
   InspiralTemplate *params
   )

{ /* </lalVerbatim>  */

  
  UINT4 count;
  InspiralInit paramsInit;

  INITSTATUS (status, "LALInspiralWave3", LALINSPIRALWAVE3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  
  LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &(paramsInit.func),
					 &(paramsInit.ak), params);
  CHECKSTATUSPTR(status);


  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  memset( output->data, 0, output->length * sizeof(REAL4) );

  /* Call the engine function */
  LALInspiralWave3Engine(status->statusPtr, output, NULL, NULL, 
			NULL, NULL, &count, params, &paramsInit);

  CHECKSTATUSPTR(status);

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

  UINT4 count;

  InspiralInit paramsInit;

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

  LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &(paramsInit.func), 
					&(paramsInit.ak), params);
  CHECKSTATUSPTR(status);

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  /* Initialise the waveforms to zero */
  memset(output1->data, 0, output1->length * sizeof(REAL4));
  memset(output2->data, 0, output2->length * sizeof(REAL4));

  /* Call the engine function */
  LALInspiralWave3Engine(status->statusPtr, output1, output2, NULL, 
			    NULL, NULL, &count, params, &paramsInit);
  CHECKSTATUSPTR(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}





NRCSID (LALINSPIRALWAVE3FORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3ForInjectionCP"> */

void 
LALInspiralWave3ForInjection (
			      LALStatus        *status,
			      CoherentGW       *waveform,   
			      InspiralTemplate *params,
			      PPNParamStruc  *ppnParams)
     
     /* </lalVerbatim>  */
{
  
  UINT4 count, i;
  expnFunc func;
  expnCoeffs ak;
  REAL4Vector *a=NULL;
  REAL4Vector *ff=NULL ;
  REAL8Vector *phiv=NULL;
  CreateVectorSequenceIn in;

  REAL8 phiC;/* phase at coalescence */

  
  InspiralInit paramsInit;  
  
 /** -- -- */
  INITSTATUS (status, "LALInspiralWave3ForInjection", LALINSPIRALWAVE3FORINJECTIONC);
  ATTATCHSTATUSPTR(status);
  
  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);  
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  
  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);   
  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
      
    }
  
  /* Now we can allocate memory and vector for coherentGW structure*/     
  LALSCreateVector(status->statusPtr, &ff, paramsInit.nbins);
  CHECKSTATUSPTR(status);   
  LALSCreateVector(status->statusPtr, &a, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);   
  LALDCreateVector(status->statusPtr, &phiv, paramsInit.nbins);
  CHECKSTATUSPTR(status);

 /* By default the waveform is empty */
 
  memset(ff->data, 0, ff->length * sizeof(REAL4));
  memset(a->data,  0, a->length * sizeof(REAL4));
  memset(phiv->data, 0, phiv->length * sizeof(REAL8));

  /* Call the engine function */
  LALInspiralWave3Engine(status->statusPtr, NULL, NULL, a, ff, phiv, &count, params, &paramsInit);
  BEGINFAIL( status ) {
     LALSDestroyVector(status->statusPtr, &ff);
     CHECKSTATUSPTR(status);
     LALSDestroyVector(status->statusPtr, &a);
     CHECKSTATUSPTR(status);
     LALDDestroyVector(status->statusPtr, &phiv);
     CHECKSTATUSPTR(status);
  }
  ENDFAIL(status);

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phiv->length; i++)
  {
    if (phiv->data[i] != 0.0) break;
    /* If the waveform returned is empty, return now */
    if (i == phiv->length - 1)
    {
      LALSDestroyVector(status->statusPtr, &ff);
      CHECKSTATUSPTR(status);
      LALSDestroyVector(status->statusPtr, &a);
      CHECKSTATUSPTR(status);
      LALDDestroyVector(status->statusPtr, &phiv);
      CHECKSTATUSPTR(status);

      DETATCHSTATUSPTR(status);
      RETURN(status);
    }
  }

  /*  if ( (phase/2./LAL_PI) < 2. ){
    sprintf(message, "The waveform has only %lf cycles; we don't keep waveform with less than 2 cycles.", 
	       (double)phase/2./(double)LAL_PI );
    LALWarning(status, message);


  }
  else*/
    {
      
      /*wrap the phase vector*/
      phiC =  phiv->data[count-1] ;
      for (i=0; i<count;i++)
	{
	  phiv->data[i] =  phiv->data[i] -phiC + ppnParams->phi;
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
      
      
      
      
      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
	= 1./params->tSampling;
      
      waveform->a->sampleUnits = lalStrainUnit;
      waveform->f->sampleUnits = lalHertzUnit;
      waveform->phi->sampleUnits = lalDimensionlessUnit;  
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      LALSnprintf( waveform->a->name, LALNameLength, "T3 inspiral amplitudes" );
      LALSnprintf( waveform->f->name, LALNameLength, "T3 inspiral frequency" );
      LALSnprintf( waveform->phi->name, LALNameLength, "T3  inspiral phase" );
      
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
    }    /*end of coherentGW storage */  


  /* --- free memory --- */
  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phiv);
  CHECKSTATUSPTR(status);



  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* Engine function used to generate the waveforms */
static void
LALInspiralWave3Engine(
                LALStatus        *status,
                REAL4Vector      *output1,
                REAL4Vector      *output2,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phiv,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )

{
  INT4 i, startShift, count;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi1, phi;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase, omega;
  DFindRootIn rootIn;
  expnFunc func;
  expnCoeffs ak;
  ChirptimeFromFreqIn timeIn;
  void *pars;

  REAL8 temp, tempMax=0, tempMin = 0;
  
  /* Only used in injection case */
  REAL8 unitHz;
  REAL8 f2a;
  REAL8 mu;
  REAL8 mTot;
  REAL8 cosI;/* cosine of system inclination */
  REAL8 etab;
  REAL8 fFac; /* SI normalization for f and t */
  REAL8 f2aFac;/* factor multiplying f in amplitude function */
  REAL8 apFac, acFac;/* extra factor in plus and cross amplitudes */


  INITSTATUS (status, "LALInspiralWave3Engine", LALINSPIRALWAVE3TEMPLATESC);
  ATTATCHSTATUSPTR(status);

  ak   = paramsInit->ak;
  func = paramsInit->func;

  if (output2 || a)
      params->nStartPad = 0; /* must be zero for templates and injections */

  /* Only in injection case.. */
  if (a) {
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
  }

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */


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

  i = startShift;

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

  rootIn.xmax = c1*params->tC*3 + 5.; /* we add 5 so that if tC is small then xmax
                                         is always greater than a given value (here 5)*/
                                                                                                                             
  /* for x in [rootIn.xmin, rootIn.xmax], we search the value which gives the max frequency.
   and keep the corresponding rootIn.xmin. */
                                                                                                                             
                                                                                                                             
                                                                                                                             
  for (tc = c1*params->tC/1000.; tc < rootIn.xmax; tc+=c1*params->tC/1000.){
    LALInspiralFrequency3Wrapper(status->statusPtr, &temp,  tc , pars);
    if (temp > tempMax) {
      rootIn.xmin = tc;
      tempMax = temp;
    }
    if (temp < tempMin) {
      tempMin = temp;
    }
  }
                                                                                                                             
  /* if we have found a value positive then everything should be fine in the
     BissectionFindRoot function */
  if (tempMax > 0  &&  tempMin < 0){
    LALDBisectionFindRoot (status->statusPtr, &tc, &rootIn, pars);
    CHECKSTATUSPTR(status);
  }
  else if (a)
  {
    /* Otherwise we return an empty waveform for injection */
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }
  else
    /* Or abort if not injection */
    ABORT( status, LALINSPIRALH_EROOTINIT, LALINSPIRALH_MSGEROOTINIT );
                                                                                                                             
  tc /= c1;
                                                                                                                             
  tc += params->startTime;       /* Add user given startTime to instant of
                                     coalescence of the compact objects */

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

  while (f < fHigh && t < tmax && f > fOld) 
  {
    /* Check we don't write past the end of the vector */
    if ((output1 && (i >= output1->length)) || (ff && (count >= ff->length)))
    {
        ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
    }

    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = params->signalAmplitude * v*v;
    if (output1)
    { 
      output1->data[i] = (REAL4) amp * cos(phase+phi0);
      if (output2)
        output2->data[i] = (REAL4) amp * cos(phase+phi1);
    }
    else if (a) {
      omega = v*v*v;
      
      ff->data[count] = (REAL4)(omega/unitHz);
      f2a = pow(f2aFac * omega, 2. / 3.);

      a->data[2*count]          = (REAL4)(4.*apFac * f2a);
      a->data[2*count+1]        = (REAL4)(4.*acFac * f2a);
      phiv->data[count]          = (REAL8)(phase);
    }
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
  if (output1 && !output2) params->tC = t;  
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/
  *countback = count;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
