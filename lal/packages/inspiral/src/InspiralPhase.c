#include <math.h>
#include "LALStdlib.h"
#include "Inspiral.h"
#include "Integrate.h"

NRCSID (INSPIRALPHASEC, "$Id$"); 


void LALInspiralPhase (LALStatus *status,
	            REAL8 *phiofv,
	            REAL8 v,
	            void *params)
{

   void *funcParams;
   DIntegrateIn intinp;
   PhiofVIntegrandIn in2;
   InspiralPhaseIn *in1;
   REAL8 sign;
   REAL8 answer;

  INITSTATUS (status, "LALInspiralPhase", INSPIRALPHASEC);
  ATTATCHSTATUSPTR (status);

  ASSERT (phiofv, status, INSPIRALPHASE_ENULL, INSPIRALPHASE_MSGENULL);
  ASSERT (params, status, INSPIRALPHASE_ENULL, INSPIRALPHASE_MSGENULL);

   sign = 1.0;
  

   in1 = (InspiralPhaseIn *) params;


   intinp.function = LALPhiofVIntegrand;
   intinp.xmin = in1->v0;
   intinp.xmax = v;
   intinp.type = ClosedInterval;

   in2.dEnergy = in1->dEnergy;
   in2.flux = in1->flux;
   in2.coeffs = in1->coeffs;

   funcParams = (void *) &in2;



   if (v==in1->v0)
     {
     *phiofv = in1->phi0;
     DETATCHSTATUSPTR (status);
     RETURN (status);
     }

   if(in1->v0 > v)
       {
	intinp.xmin = v;
	intinp.xmax = in1->v0;
	sign = -1.0;
	}


   LALDRombergIntegrate (status->statusPtr, &answer, &intinp, funcParams);
   CHECKSTATUSPTR (status);

   printf("inside InspiralPhase, in1->phi0,p= %e %e\n",in1->phi0,answer);

   *phiofv = in1->phi0 - 2.0*sign*answer;




   DETATCHSTATUSPTR (status);
   RETURN (status);


}

