/*  <lalVerbatim file="LALInspiralWave1CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave1.c} and \texttt{LALInspiralWave1Templates.c}}

The code \texttt{LALInspiralWave1} generates an time-domain inspiral waveform corresponding to the 
\texttt{approximant} \texttt{TaylorT1} and \texttt{PadeT1} as outlined in the
documentation for the function \texttt{LALInspiralWave}. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave1CP}
\index{\verb&LALInspiralWave1()&}
\begin{itemize}
\item {\tt signal:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALInspiralWave1TemplatesCP}
\index{\verb&LALInspiralWave1Templates()&}
\begin{itemize}
\item {\tt signal1:} Output containing the 0-phase inspiral waveform.
\item {\tt signal2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

\texttt{LALInspiralWave1} is called if the user has specified the 
\texttt{enum} \texttt{approximant} to be
either \texttt{TaylorT1} or \texttt{PadeT1}.
{\tt LALInspiralWave1Templates} is exactly the same as \texttt{LALInspiralWave1,} except that
it generates two templates one for which the starting phase is 
\texttt{params.startPhase} and the other for which the phase is
\texttt{params.startPhase + $\pi/2$}.


\subsubsection*{Algorithm}
This code uses a fourth-order Runge-Kutta algorithm to solve the ODEs 
in Equation (\ref{eq:ode2}).

\subsubsection*{Uses}

\texttt{LALInspiralSetup}\\
\texttt{LALInspiralChooseModel}\\
\texttt{LALInspiralVelocity}\\
\texttt{LALInspiralPhasing1}\\
\texttt{LALInspiralDerivatives}\\
\texttt{LALRungeKutta4}.
 

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave1CV}}

</lalLaTeX>  */

/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

NRCSID (LALINSPIRALWAVE1C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1CP"> */
void 
LALInspiralWave1(
   LALStatus        *status,
   REAL4Vector      *signal,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 n=2, count;
   REAL8 m, dt, t, v, p, h, f, fu, fHigh, piM;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;
   
   INITSTATUS(status, "LALInspiralWave1", LALINSPIRALWAVE1C);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);
   
   values.length = dvalues.length = valuesNew.length =
   yt.length = dym.length = dyt.length = n;
   dummy.length = n * 6;
   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * n * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[n];
   valuesNew.data = &dummy.data[2*n];
   yt.data = &dummy.data[3*n];
   dym.data = &dummy.data[4*n];
   dyt.data = &dummy.data[5*n];
   
   params->nStartPad = 0;
   params->nEndPad 	= 0.0;

   m = ak.totalmass;
   dt = 1./params->tSampling;
   fu = params->fCutoff;
   if (fu) 
      fHigh = (fu < ak.flso) ? fu : ak.flso; 
   else 
      fHigh = ak.flso;

/* 
    Check that the highest frequency is less than half 
    the sampling frequency - the Nyquist theorem 
*/

   ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


   ASSERT(ak.totalmass > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   t = 0.0;
   in1.t = t;
   in1.t0=ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;

   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams = (void *) &in3;

   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

   f = (v*v*v)/(LAL_PI*m);

   LALInspiralPhasing1(status->statusPtr, &p, v, &in2);
   CHECKSTATUSPTR(status);


   *(values.data) = v; 
   *(values.data+1) = p;

   in4.function = LALInspiralDerivatives;
   in4.x = t;
   in4.y = &values;
   in4.h = dt;
   in4.n = n;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   count = 0;
  /* while (count < params->nStartPad) 
   {
       *(signal->data + count) = 0.;
       count++;
   }*/

   t = 0.0;
   piM = LAL_PI * m;
   do {
      /* Check we're not writing beyond the end of the vector */
      if (count >= signal->length) {
          LALFree(dummy.data);
          ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }

      h = params->signalAmplitude * v*v * cos(p);
      LALInspiralDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = &dvalues;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);
      *(signal->data+count) = (REAL4) h;
      t = (++count-params->nStartPad) * dt;
      f = v*v*v/piM;
   } while (t < ak.tn &&  f<fHigh);
   params->fFinal = f;
   params->tC = t;       

   while (count < (int)signal->length) 
   {
       *(signal->data + count) = 0.;
       count++;
   }

   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}



/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVE1TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1TemplatesCP"> */
void 
LALInspiralWave1Templates(
   LALStatus        *status,
   REAL4Vector      *signal1,
   REAL4Vector      *signal2,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 n=2, count;
   REAL8 amp, m, dt, t, v, p, h1, h2, f, fu, fHigh, piM;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;

   INITSTATUS(status, "LALInspiralWave1Templates", LALINSPIRALWAVE1TEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   values.length = dvalues.length = valuesNew.length =
   yt.length = dym.length = dyt.length = n;
   dummy.length = n * 6;
   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * n * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[n];
   valuesNew.data = &dummy.data[2*n];
   yt.data = &dummy.data[3*n];
   dym.data = &dummy.data[4*n];
   dyt.data = &dummy.data[5*n];

   m = ak.totalmass;
   dt = 1./params->tSampling;

   ASSERT(ak.totalmass > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   t = 0.0;
   in1.t = t;
   in1.t0=ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;

   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams = (void *) &in3;

   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

   piM = LAL_PI * m;
   f = (v*v*v)/piM;

   fu = params->fCutoff;
   if (fu) 
      fHigh = (fu < ak.flso) ? fu : ak.flso; 
   else 
      fHigh = ak.flso;
   f = (v*v*v)/(LAL_PI*m);

   LALInspiralPhasing1(status->statusPtr, &p, v, &in2);
   CHECKSTATUSPTR(status);


   *(values.data) = v; 
   *(values.data+1) = p;

   in4.function = LALInspiralDerivatives;
   in4.x = t;
   in4.y = &values;
   in4.h = dt;
   in4.n = n;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   count = 0;
   params->nStartPad = 0;
   /*while (count < params->nStartPad) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }*/

   t = 0.0;
   do {
      /* Free up memory and abort if writing beyond the end of vector*/
      if (count >= signal1->length) {
          LALFree(dummy.data);
          ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }
      amp = params->signalAmplitude * v*v;
      h1 = amp * cos(p);
      h2 = amp * cos(p+LAL_PI_2);
      *(signal1->data + count) = (REAL4) h1;
      *(signal2->data + count) = (REAL4) h2;
      LALInspiralDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = &dvalues;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);
      t = (++count-params->nStartPad) * dt;
      f = v*v*v/piM;
   } while (t < ak.tn && f<fHigh);
   
   params->vFinal = p;
   params->fFinal = f;
   params->tC = t;
   
   while (count < (int)signal1->length) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}

/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms for injection packages T.Cokelaer sept 2003
*/

NRCSID (LALINSPIRALWAVE1FORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1ForInjectionCP"> */
void 
LALInspiralWave1ForInjection(
			     LALStatus        *status,
			     CoherentGW       *waveform,
			     InspiralTemplate *params,
			     PPNParamStruc  *ppnParams			     
			     )
{ /* </lalVerbatim>  */
  
  INT4 n = 2  ;
  UINT4 count, i; 
  REAL8  dt, t, v, p, f, fu, fHigh, piM, omega;
  REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
  
  REAL4Vector *a	= NULL;/* pointers to generated amplitude  data */
  REAL4Vector *ff	= NULL ;/* pointers to generated  frequency data */
  REAL8Vector *phi	= NULL;/* pointer to generated phase data */
  
  CreateVectorSequenceIn in;
  TofVIn in1;
  InspiralPhaseIn in2;
  InspiralDerivativesIn in3;
  rk4In in4;
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

  INITSTATUS(status, "LALInspiralWave1ForInjection", LALINSPIRALWAVE1TEMPLATESC);
  ATTATCHSTATUSPTR(status);
  
  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);  
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  
  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);   
  
  if (paramsInit.nbins == 0){
      DETATCHSTATUSPTR(status);
      RETURN (status);
  }
  func 	= paramsInit.func;
  ak   	= paramsInit.ak;
  dt 	= 1./params->tSampling;


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
  acFac *= 2.0*cosI;
  params->nStartPad = 0;

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
      phi->data[count]          = 0.;
      a->data[2*count]          = 0.;
    }
  count = 0;

  /* some others parameters for inspiral computation*/
  values.length = dvalues.length = valuesNew.length =
    yt.length = dym.length = dyt.length = n;
  dummy.length = n * 6;
  if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * n * 6))) {
    ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
  }
  
  values.data    = &dummy.data[0];
  dvalues.data   = &dummy.data[n];
  valuesNew.data = &dummy.data[2*n];
  yt.data        = &dummy.data[3*n];
  dym.data       = &dummy.data[4*n];
  dyt.data       = &dummy.data[5*n];
  

  dt = 1./params->tSampling;

  
 /* --- actual computation start here --- */
  t 		= 0.0;
  in1.t 	= t;
  in1.t0	= ak.t0;
  in1.v0 	= ak.v0;
  in1.vlso 	= ak.vlso;
  in1.totalmass = ak.totalmass;
  in1.dEnergy 	= func.dEnergy;
  in1.flux 	= func.flux;
  in1.coeffs 	= &ak;
  LALInspiralVelocity(status->statusPtr, &v, &in1);
  CHECKSTATUSPTR(status);
  
  in2.v0 	= ak.v0;
  in2.phi0 	= params->startPhase;
  in2.dEnergy 	= func.dEnergy;
  in2.flux 	= func.flux;
  in2.coeffs 	= &ak;
  LALInspiralPhasing1(status->statusPtr, &p, v, &in2);
  CHECKSTATUSPTR(status);
  
  in3.totalmass = ak.totalmass;
  in3.dEnergy 	= func.dEnergy;
  in3.flux 	= func.flux;
  in3.coeffs 	= &ak;
  funcParams 	= (void *) &in3;
  
  piM = LAL_PI * ak.totalmass;
  f = (v*v*v)/piM;
  
  /* a meetre ailleurs. Au debut avec les asserts?*/
  fu = params->fCutoff;
  if (fu) 
    fHigh = (fu < ak.flso) ? fu : ak.flso; 
  else 
    fHigh = ak.flso;
  f = (v*v*v)/(LAL_PI*ak.totalmass);
  
  *(values.data) = v; 
  *(values.data+1) = p;
  
  in4.function 	= LALInspiralDerivatives;
  in4.x 	= t;
  in4.y 	= &values;
  in4.h 	= dt;
  in4.n 	= n;
  in4.yt 	= &yt;
  in4.dym 	= &dym;
  in4.dyt 	= &dyt;
  
  count = 0;
  t    	= 0.0;
  
  /* --- Main loop --- */   
  do {
    /* If trying to write beyond the end of a vector
       free up memory and abort. */
    if (count >= ff->length) {
        LALSDestroyVector(status->statusPtr, &ff);
        CHECKSTATUSPTR(status);
        LALSDestroyVector(status->statusPtr, &a);
        CHECKSTATUSPTR(status);
        LALDDestroyVector(status->statusPtr, &phi);
        CHECKSTATUSPTR(status);

        LALFree(dummy.data);
        ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
    }  
  
    omega = v*v*v;
    
    ff->data[count]           = (REAL4)(omega/unitHz);
    f2a                       = pow (f2aFac * omega, 2./3.);
    a->data[2*count]          = (REAL4)(4.*apFac * f2a);
    a->data[2*count+1]        = (REAL4)(4.*acFac * f2a);
    phi->data[count]          = (REAL8)(p);
    
    LALInspiralDerivatives(&values, &dvalues, funcParams);
    CHECKSTATUSPTR(status);
    
    in4.dydx = &dvalues;
    in4.x=t;
    
    LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
    CHECKSTATUSPTR(status);
    
    *(values.data) = v = *(valuesNew.data);
    *(values.data+1) = p = *(valuesNew.data+1);
    
    t = (++count-params->nStartPad) * dt;
    f = v*v*v/piM;
  } while (t < ak.tn && f<fHigh);
  
  params->fFinal = ff->data[count-1];
  sprintf(message, "cycles = %f", p/2./3.14159);
  LALInfo(status, message);

  if ( (p/2./LAL_PI) < 2 ){
    sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.", 
	       p/2./(double)LAL_PI );
    LALWarning(status, message);
  }
  else {

      /*wrap the phase vector*/
      phiC =  phi->data[count-1] ;
      for (i = 0; i < count; i++)
	{
	  phi->data[i] =  phi->data[i] - phiC + ppnParams->phi;
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
      
      
      
      in.length = (UINT4)(count);
      in.vectorLength = 2;
      
      LALSCreateVectorSequence( status->statusPtr, &( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);      
      
      LALSCreateVector( status->statusPtr, &( waveform->f->data ), count);
      CHECKSTATUSPTR(status);      
      
      LALDCreateVector( status->statusPtr, &( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);        
      
      memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));
      
      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
	= ppnParams->deltaT;
      
      waveform->a->sampleUnits    = lalStrainUnit;
      waveform->f->sampleUnits    = lalHertzUnit;
      waveform->phi->sampleUnits  = lalDimensionlessUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      LALSnprintf( waveform->a->name, LALNameLength,   "T1 inspiral amplitude" );
      LALSnprintf( waveform->f->name, LALNameLength,   "T1 inspiral frequency" );
      LALSnprintf( waveform->phi->name, LALNameLength, "T1 inspiral phase" );
      
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
    } /*end of coherentGW storage */


  /* --- free memory --- */
  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phi);
  CHECKSTATUSPTR(status);

  LALFree(dummy.data);
   
  DETATCHSTATUSPTR(status);
  RETURN (status);
}
