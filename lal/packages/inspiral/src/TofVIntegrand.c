#include <math.h>
#include "LALStdlib.h"
#include "Inspiral.h"

NRCSID (TOFVINTEGRANDC, "$Id$"); 


void TofVIntegrand (Status *status,
		   REAL8 *integrand,
		   REAL8 v,
		   void *params)
{

   TofVIntegrandIn *ak;

   INITSTATUS (status, "TofVIntegrand", TOFVINTEGRANDC);

   ASSERT (integrand, status, TOFVINTEGRAND_ENULL, TOFVINTEGRAND_MSGENULL);
   ASSERT (params, status, TOFVINTEGRAND_ENULL, TOFVINTEGRAND_MSGENULL);



   ak = (TofVIntegrandIn *) params;

   *integrand = ak->dEnergy(v, ak->coeffs)/ak->flux(v, ak->coeffs);

   RETURN (status);


}
