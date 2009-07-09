/**
*  Copyright (C) 2007
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*  <lalVerbatim file="LALInspiralEccentricityCV">
Author: Devanka Pathak and Thomas Cokelaer
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralEccentricity.c} and \texttt{LALInspiralEccentricityTemplates.c}}

The code \texttt{LALInspiralEccentricity} generates a time-domain inspiral waveform corresponding to the
\texttt{approximant} \texttt{Eccentricity} as outlined PRD 60 for the Newtonian case.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralEccentricityCP}
\index{\verb&LALInspiralEccentricity()&}
\begin{itemize}
\item {\tt signalvec:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters and eccentricity.
\end{itemize}

\input{LALInspiralEccentricityTemplatesCP}
\index{\verb&LALInspiralEccentricityTemplates()&}
\begin{itemize}
\item {\tt signalvec1:} Output containing the 0-phase inspiral waveform.
\item {\tt signalvec2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

\texttt{LALInspiralEccentricity} is called if the user has specified the
\texttt{enum} \texttt{approximant} to be
either \texttt{TaylorT1} or \texttt{PadeT1}.
{\tt LALInspiralEccentricityTemplates} is exactly the same as \texttt{LALInspiralEccentricity,} except that
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

\vfill{\footnotesize\input{LALInspiralEccentricityCV}}

</lalLaTeX>  */

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

/* macro to "use" unused function parameters */
#define UNUSED(expr) do { (void)(expr); } while (0)

/* structure to provide M and eta. */
typedef struct
tagecc_CBC_ODE_Struct
{
  double totalMass; /* M=m1+m2 */
  double eta;       /* eta=m1 m2/M^2 */
}
ecc_CBC_ODE_Input;


void
LALInspiralEccentricityDerivatives (
   REAL8Vector *values,
   REAL8Vector *dvalues,
   void        *params
   );

static void
LALInspiralEccentricityEngine(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   REAL4Vector      *a,
   REAL4Vector      *ff,
   REAL8Vector      *phi,
   INT4             *countback,
   InspiralTemplate *params
   );


NRCSID (LALINSPIRALECCENTRICITYC, "$Id$");

/*  <lalVerbatim file="LALInspiralEccentricityCP"> */
void
LALInspiralEccentricity(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 count;

   INITSTATUS(status, "LALInspiralEccentricity", LALINSPIRALECCENTRICITYC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   /* Initially the waveform is empty*/
   memset(signalvec->data, 0, signalvec->length*sizeof(REAL4));

   /*Call the engine function*/
   LALInspiralEccentricityEngine(status->statusPtr, signalvec, NULL, NULL, NULL, NULL, &count, params);
   CHECKSTATUSPTR(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}



/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/

NRCSID (LALINSPIRALECCENTRICITYTEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralEccentricityTemplatesCP"> */
void
LALInspiralEccentricityTemplates(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 count;

   INITSTATUS(status, "LALInspiralEccentricityTemplates", LALINSPIRALECCENTRICITYTEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   /* Initially the waveforms are empty */
   memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
   memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

   /* Call the engine function */
   LALInspiralEccentricityEngine(status->statusPtr, signalvec1, signalvec2, NULL, NULL, NULL, &count, params);
   CHECKSTATUSPTR(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms for injection packages T.Cokelaer sept 2003
*/

NRCSID (LALINSPIRALECCENTRICITYFORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralEccentricityForInjectionCP"> */
void
LALInspiralEccentricityForInjection(
			     LALStatus        *status,
			     CoherentGW       *waveform,
			     InspiralTemplate *params,
			     PPNParamStruc  *ppnParams
			     )
{ /* </lalVerbatim>  */

  INT4        count, i;
  REAL8       p, phiC;

  REAL4Vector a;           /* pointers to generated amplitude  data */
  REAL4Vector ff;          /* pointers to generated  frequency data */
  REAL8Vector phi;         /* generated phase data */

  CreateVectorSequenceIn in;

  CHAR message[256];

  InspiralInit paramsInit;

  INITSTATUS(status, "LALInspiralEccentricityForInjection", LALINSPIRALECCENTRICITYTEMPLATESC);
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

  /* Now we can allocate memory and vector for coherentGW structure*/

  ff.length  = paramsInit.nbins;
  a.length   = 2* paramsInit.nbins;
  phi.length = paramsInit.nbins;

  ff.data = (REAL4 *) LALCalloc(paramsInit.nbins, sizeof(REAL4));
  a.data  = (REAL4 *) LALCalloc(2 * paramsInit.nbins, sizeof(REAL4));
  phi.data= (REAL8 *) LALCalloc(paramsInit.nbins, sizeof(REAL8));

  /* Check momory allocation is okay */
  if (!(ff.data) || !(a.data) || !(phi.data))
  {
    if (ff.data)  LALFree(ff.data);
    if (a.data)   LALFree(a.data);
    if (phi.data) LALFree(phi.data);

    ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
  }

  count = 0;

  /* Call the engine function */
  LALInspiralEccentricityEngine(status->statusPtr, NULL, NULL, &a, &ff, &phi, &count, params);
  BEGINFAIL( status )
  {
    LALFree(ff.data);
    LALFree(a.data);
    LALFree(phi.data);
  }
  ENDFAIL( status );

  p = phi.data[count-1];

  params->fFinal = ff.data[count-1];
  sprintf(message, "cycles = %f", p/(double)LAL_TWOPI);
  LALInfo(status, message);

  if ( (INT4)(p/LAL_TWOPI) < 2 ){
    sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.",
	       p/(double)LAL_TWOPI );
    XLALPrintError(message);
    LALWarning(status, message);
  }

      /*wrap the phase vector*/
      phiC =  phi.data[count-1] ;
      for (i = 0; i < count; i++)
	{
	  phi.data[i] =  phi.data[i] - phiC + ppnParams->phi;
	}

      /* Allocate the waveform structures. */
      if ( ( waveform->a = (REAL4TimeVectorSeries *)
	     LALCalloc(1, sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      if ( ( waveform->f = (REAL4TimeSeries *)
	     LALCalloc(1, sizeof(REAL4TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      if ( ( waveform->phi = (REAL8TimeSeries *)
	     LALCalloc(1, sizeof(REAL8TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	LALFree( waveform->f ); waveform->f = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }


      in.length = (UINT4)(count);
      in.vectorLength = 2;

      LALSCreateVectorSequence( status->statusPtr, &( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);

      LALSCreateVector( status->statusPtr, &( waveform->f->data ), count);
      CHECKSTATUSPTR(status);

      LALDCreateVector( status->statusPtr, &( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);

      memcpy(waveform->f->data->data , ff.data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a.data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi.data, count*(sizeof(REAL8)));

      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
	= ppnParams->deltaT;

      waveform->a->sampleUnits    = lalStrainUnit;
      waveform->f->sampleUnits    = lalHertzUnit;
      waveform->phi->sampleUnits  = lalDimensionlessUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      snprintf( waveform->a->name, LALNameLength,   "T1 inspiral amplitude" );
      snprintf( waveform->f->name, LALNameLength,   "T1 inspiral frequency" );
      snprintf( waveform->phi->name, LALNameLength, "T1 inspiral phase" );

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

  /* --- free memory --- */
  LALFree(ff.data);
  LALFree(a.data);
  LALFree(phi.data);

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/*
 *  Engine function for use by other LALInspiralEccentricity* functions
 *  Craig Robinson April 2005
 */

NRCSID (LALINSPIRALECCENTRICITYENGINEC, "$Id$");

void
LALInspiralEccentricityEngine(
		LALStatus        *status,
		REAL4Vector      *signalvec1,
		REAL4Vector      *signalvec2,
		REAL4Vector      *a,
		REAL4Vector      *ff,
		REAL8Vector      *phi,
		INT4             *countback,
		InspiralTemplate *params)
{
   INT4 number_of_diff_equations = 3;
   INT4 count = 0;
   ecc_CBC_ODE_Input in3;
   REAL8 phase;
   REAL8 orbital_element_p,orbital_element_e_squared;
   REAL8 orbital_element_e;
   REAL8 twoPhim2Beta = 0;
   REAL8 threePhim2Beta = 0;
   REAL8 phim2Beta = 0;
   REAL8 twoBeta;
   REAL8 rbyM=1e6, rbyMFlso=6.;
   REAL8 sin2Beta,cos2Beta,iota,onepCosSqI, SinSqI, cosI, e0, f_min, beta, p0;



   REAL8 amp, m, dt, t,  h1, h2, f,  fHigh, piM, fu;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   INT4 done=0;
   rk4In in4;
   rk4GSLIntegrator *integrator;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;

#if 0
   REAL8 mTot = 0;
   REAL8 unitHz = 0;
   REAL8 f2a = 0;
   REAL8 mu = 0;
   REAL8 etab = 0;
   REAL8 fFac = 0; /* SI normalization for f and t */
   REAL8 f2aFac = 0;/* factor multiplying f in amplitude function */
   REAL8 apFac = 0, acFac = 0;/* extra factor in plus and cross amplitudes */
#endif

   /* ff and phi are unused in this function */
   UNUSED(ff);
   UNUSED(phi);

   INITSTATUS(status, "LALInspiralEccentricityEngine", LALINSPIRALECCENTRICITYENGINEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   m = ak.totalmass = params->mass1+params->mass2;

   values.length = dvalues.length = valuesNew.length =
   yt.length = dym.length = dyt.length = number_of_diff_equations;
   dummy.length = number_of_diff_equations * 6;
   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * number_of_diff_equations * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[number_of_diff_equations];
   valuesNew.data = &dummy.data[2*number_of_diff_equations];
   yt.data = &dummy.data[3*number_of_diff_equations];
   dym.data = &dummy.data[4*number_of_diff_equations];
   dyt.data = &dummy.data[5*number_of_diff_equations];

/*   m = ak.totalmass;*/
   dt = 1./params->tSampling;

   if (a)
   {
     /*to be implemented. */


/*      mTot   = params->mass1 + params->mass2;
      etab   = params->mass1 * params->mass2;
      etab  /= mTot;
      etab  /= mTot;
      unitHz = mTot *LAL_MTSUN_SI*(REAL8)LAL_PI;
      cosI   = cos( params->inclination );
      mu     = etab * mTot;
      fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
      f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
      apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
      apFac *= 1.0 + cosI*cosI;
      acFac *= 2.0*cosI;
      params->nStartPad = 0;
      */
   }

   ASSERT(ak.totalmass > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   t = 0.0;



   in3.totalMass = (params->mass1 + params->mass2) * LAL_MTSUN_SI ;
   in3.eta = (params->mass1 * params->mass2) /(params->mass1 + params->mass2) / (params->mass1 + params->mass2);
   funcParams = (void *) &in3;


   piM = LAL_PI * m * LAL_MTSUN_SI;
/*   f = (v*v*v)/piM;

   fu = params->fCutoff;
   if (fu)
      fHigh = (fu < ak.flso) ? fu : ak.flso;
   else
      fHigh = ak.flso;
   f = (v*v*v)/(LAL_PI*m);

   ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
*/


   /* e0 is set at f_min */
   e0 = params->eccentricity;

   /* the second harmonic will start at fLower*2/3 */
   f_min = params->fLower;
   iota = params->inclination; /*overwritten later */

   beta = 0.;
   twoBeta = 2.* beta;
   cos2Beta = cos(twoBeta);
   sin2Beta = sin(twoBeta);
   iota = LAL_PI/4.;
   onepCosSqI = 1. + cos(iota) * cos(iota);
   SinSqI = sin(iota) * sin(iota);
   cosI = cos(iota);

   p0 = (1. - e0*e0)/pow(2. * LAL_PI * m * LAL_MTSUN_SI* f_min/3. , 2./3.);

   *(values.data) = orbital_element_p = p0;
   *(values.data+1) = phase = params->startPhase;
   *(values.data+2) = orbital_element_e = e0;




   in4.function = LALInspiralEccentricityDerivatives;
   in4.x = t;
   in4.y = &values;
   in4.h = dt;
   in4.n = number_of_diff_equations;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   xlalErrno = 0;
   /* Initialize GSL integrator */
   if (!(integrator = XLALRungeKutta4Init(number_of_diff_equations, &in4)))
   {
     INT4 errNum = XLALClearErrno();
     LALFree(dummy.data);

     if (errNum == XLAL_ENOMEM)
       ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
     else
       ABORTXLAL( status );
   }

   count = 0;
   if (signalvec2) {
   params->nStartPad = 0;} /* for template generation, that value must be zero*/

   else if (signalvec1) {
     count = params->nStartPad;
   }

   t = 0.0;


   fu = params->fCutoff;
   if (fu)
      fHigh = (fu < ak.flso) ? fu : ak.flso;
   else
      fHigh = ak.flso;

   f = 1./(pow(orbital_element_p, 3./2.))/piM;


   /*fprintf(stderr, "fFinal = %f %f %f %f\n", fu,fHigh,f,ak.flso);*/
 done = 0;
   do {
      /* Free up memory and abort if writing beyond the end of vector*/
      /*if ((signalvec1 && (UINT4)count >= signalvec1->length) || (ff && (UINT4)count >= ff->length))*/
      if ((signalvec1 && (UINT4)count >= signalvec1->length))
      {
          XLALRungeKutta4Free( integrator );
          LALFree(dummy.data);
          ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }

      /* Non-injection case */
      if (signalvec1)
      {
        twoPhim2Beta = 2.* phase - twoBeta;
        phim2Beta = phase - twoBeta;
        threePhim2Beta = 3.* phase - twoBeta;
        orbital_element_e_squared = orbital_element_e * orbital_element_e;
        amp = params->signalAmplitude / orbital_element_p;

/*        fprintf(stderr, "%e %e %e %e %e\n", twoBeta, twoPhim2Beta, phim2Beta, threePhim2Beta, orbital_element_e_squared);*/


        h1 = amp * ( ( 2. * cos(twoPhim2Beta) + 2.5 * orbital_element_e * cos(phim2Beta)
          + 0.5 * orbital_element_e * cos(threePhim2Beta) + orbital_element_e_squared * cos2Beta) * onepCosSqI +
          + ( orbital_element_e * cos(orbital_element_p) + orbital_element_e_squared) * SinSqI);
        if ((f >= params->fLower) && (done == 0))
        {
        /*fprintf(stderr, "freq=%e p=%e, e=%e, phase = %e\n", f,orbital_element_p, orbital_element_e, phase);fflush(stderr);
*/
        params->alpha1 =  orbital_element_e;
        done = 1;
         }
         /*if (f>=params->fLower)*/
        {
          *(signalvec1->data + count) = (REAL4) h1;
         }

	 if (signalvec2)
	 {
            h2 = amp * ( ( 4. * sin(twoPhim2Beta) + 5 * orbital_element_e * sin(phim2Beta)
          + orbital_element_e * sin(threePhim2Beta) - 2. * orbital_element_e_squared * sin2Beta) * cosI);
/*           if (f>=params->fLower)*/
           {
              *(signalvec2->data + count) = (REAL4) h2;
            }
	 }
      }

      /* Injection case */
      else if (a)
      {
        /*to be done*/
        /*
        omega = v*v*v;

          ff->data[count]       = (REAL4)(omega/unitHz);
          f2a                   = pow (f2aFac * omega, 2./3.);
          a->data[2*count]      = (REAL4)(4.*apFac * f2a);
          a->data[2*count+1]    = (REAL4)(4.*acFac * f2a);
          phi->data[count]      = (REAL8)(p);
          */
      }

      LALInspiralEccentricityDerivatives(&values, &dvalues,funcParams);
      CHECKSTATUSPTR(status);

      in4.dydx = &dvalues;
      in4.x=t;

      LALRungeKutta4(status->statusPtr, &valuesNew, integrator, funcParams);
      CHECKSTATUSPTR(status);

      *(values.data) = orbital_element_p = *(valuesNew.data);
      *(values.data+1) = phase = *(valuesNew.data+1);
      *(values.data+2) = orbital_element_e = *(valuesNew.data+2);

      t = (++count-params->nStartPad) * dt;
      /* v^3 is equal to p^(3/2)*/
      f = 1./(pow(orbital_element_p, 3./2.))/piM;
/*      fprintf(stderr, "p=%e, e=%e, phase = %e\n", orbital_element_p, orbital_element_e, phase);
      fflush(stderr);*/
      rbyM = orbital_element_p/(1.+orbital_element_e * cos(phase));
      /*fprintf(stderr, "rbyM=%e rbyMFlso=%e t=%e, ak.tn=%e f=%e e=%e\n",rbyM, rbyMFlso, t, ak.tn,f,orbital_element_e);
      fflush(stderr);*/

   } while ( (t < ak.tn) && (rbyM>rbyMFlso) && (f<fHigh));


   /*fprintf(stderr, "t=%e ak.tn=%e rbyM=%e rbyMFlso=%e f=%f fHigh=%f", t, ak.tn, rbyM, rbyMFlso, f, fHigh);*/

   /*also need to add for the fupper nyquist frequency**/
   params->vFinal = orbital_element_p;
   params->fFinal = f;
   params->tC = t;
   *countback = count;

   XLALRungeKutta4Free( integrator );
   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}




void
LALInspiralEccentricityDerivatives (
   REAL8Vector *values,
   REAL8Vector *dvalues,
   void        *params
   )
 { /* </lalVerbatim> */

  ecc_CBC_ODE_Input *par;
  double M, eta, mu, c1, e0, e2, p0, p2, p3, p4, phi;
  par = (ecc_CBC_ODE_Input *) params;

  /* affectation */
  M = par->totalMass; /* in MTSUN*/
  eta = par->eta;  /*dimensionless*/
  mu = eta * M; /*therefore in MTSUN*/
  phi = values->data[1];
  e0 = values->data[2];
  e2 = values->data[2] * values->data[2];
  p0 = values->data[0];
  p2 = values->data[0] * values->data[0];
  p3 = p2 * values->data[0];
  p4 = p3 * values->data[0];

  c1 = pow(1. - e2, 1.5);

  /* y[0] is p
     y[1] is phase
   * y[2] is e
   */

/*  fprintf(stderr, "before=%e %e %e\n", p0, e0, phi);*/
  /* Eq 6 */
  dvalues->data[1] = 1./(M * sqrt (p3)) * ( 1. + e0 * cos(phi) ) * ( 1. + e0 * cos(phi) );

  /* Eq 8 */
  dvalues->data[0] = (-64./5.) * mu / M / M * c1 / p3 * (1 + 7./8. * e2);

  /* Eq 9*/
  dvalues->data[2] = (-304./15.) * mu / M / M / p4 * c1 * e0 * (1. + 121./304. * e2);

/*  fprintf(stderr, "new values p=%e, e=%e, phase = %e\n",
      dvalues->data[0],
      dvalues->data[2],
      dvalues->data[1]);*/

  return;
}



