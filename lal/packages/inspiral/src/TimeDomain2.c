/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include "Inspiral.h"
#include "LALStdlib.h"

NRCSID (TIMEDOMAIN2C, "$Id$");

void LALTimeDomain2(LALStatus *status,
		 REAL8Vector *signal,
		 InspiralTemplate *params)
 {

   INT4 n=2, count;
   double m, dt, t, v, p, h, f;
   REAL8Vector *vandp;
   REAL8Vector *dvanddp;
   REAL8Vector *vandpnew;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   REAL8Vector *yt;
   REAL8Vector *dym;
   REAL8Vector *dyt;
   int ndx;
   double length;
   expnCoeffs ak;
   expnFunc func;


   INITSTATUS(status, "timedomain2", TIMEDOMAIN2C);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal,  status, TIMEDOMAIN2_ENULL, TIMEDOMAIN2_MSGENULL);
   ASSERT (params,  status, TIMEDOMAIN2_ENULL, TIMEDOMAIN2_MSGENULL);


   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);

   LALChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);



      length = (ak.tn) * params->tSampling + params->nStartPad + params->nEndPad;

      ndx = ceil(log10(length)/log10(2.));
      signal->length = pow(2, ndx);
      signal->data = (REAL8 *)LALMalloc(sizeof(REAL8) * signal->length);

/*      fprintf(stderr,"chirp time=%e, length=%e, power=%d\n", ak.tn,length,ndx);
*/


       vandp = (REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
       dvanddp = (REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
       vandpnew = (REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
       vandp->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
       dvanddp->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
       vandpnew->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);



   m = ak.totalmass;
   dt = 1./params->tSampling;

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

   f = pow(v,3.)/(LAL_PI*m);

   LALInspiralPhase(status->statusPtr, &p, v, &in2);
   CHECKSTATUSPTR(status);

/*   fprintf (stderr, "Initial velocity %e, frequency %e and phase %e\n", v,f,p);
*/

   *(vandp->data) = v; 
   *(vandp->data+1) = p;



	dym = (REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
	dyt = (REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
	yt = (REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
	dym->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
	dyt->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
	yt->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);

   in4.function = LALInspiralDerivatives;
   in4.x = t;
   in4.y = vandp;
   in4.h = dt;
   in4.n = n;
   in4.yt = yt;
   in4.dym = dym;
   in4.dyt = dyt;

 

   count = 0;
   while (count < params->nStartPad) *(signal->data + count++) = 0.;


   t = 0.0;
   do {
      h = params->signalAmplitude * v*v * cos(p);
      LALInspiralDerivatives(status->statusPtr, vandp, dvanddp, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = dvanddp;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, vandpnew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(vandp->data) = v = *(vandpnew->data);
      *(vandp->data+1) = p = *(vandpnew->data+1);
      *(signal->data+count) = h;
      t = (++count-params->nStartPad) * dt;
   } while (t < ak.tn);

   while (count < signal->length) *(signal->data + count++) = 0.;


	LALFree(yt->data);
	LALFree(yt);
	LALFree(dyt->data);
	LALFree(dyt);
	LALFree(dym->data);
	LALFree(dym);


   DETATCHSTATUSPTR(status);
   RETURN (status);

}
