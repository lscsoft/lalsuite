#include <math.h>
#include "LALStdlib.h"
#include "Inspiral.h"
#include "Integrate.h"

NRCSID (TOFVC, "$Id$");


void TofV (Status *status,
	   REAL8 *tofv,
	   REAL8 v,
	   void *params)
{

   void *funcParams;
   DIntegrateIn intinp;
   TofVIntegrandIn in2;
   TofVIn *in1;
   REAL8 answer;
   REAL8 sign;


   INITSTATUS (status, "TofV", TOFVC);
   ATTATCHSTATUSPTR(status);

   ASSERT (tofv, status, TOFV_ENULL, TOFV_MSGENULL);
   ASSERT (params, status, TOFV_ENULL, TOFV_MSGENULL);

   sign = 1.0;


   in1 = (TofVIn *) params;

   intinp.function = TofVIntegrand;
   intinp.xmin = in1->v0;
   intinp.xmax = v;
   intinp.type = ClosedInterval;


   in2.dEnergy = in1->dEnergy;
   in2.flux = in1->flux;
   in2.coeffs = in1->coeffs;

   funcParams = (void *) &in2;

   if (v==in1->v0)
     {
     *tofv = in1->t - in1->t0;
     DETATCHSTATUSPTR(status);
     RETURN (status);
     }

   if(in1->v0 > v)
       {
	intinp.xmin = v;
	intinp.xmax = in1->v0;
	sign = -1.0;
	}

	
   DRombergIntegrate (status->statusPtr, &answer, &intinp, funcParams);
   CHECKSTATUSPTR(status);
   *tofv = in1->t - in1->t0 + in1->totalmass*answer*sign;

   DETATCHSTATUSPTR(status);
   RETURN (status);


}

