#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>

NRCSID (TOFVINTEGRANDC, "$Id$"); 


void LALTofVIntegrand (LALStatus *status,
		   REAL8 *integrand,
		   REAL8 v,
		   void *params)
{

   TofVIntegrandIn *ak;

   INITSTATUS (status, "LALTofVIntegrand", TOFVINTEGRANDC);

   ASSERT (integrand, status, TOFVINTEGRAND_ENULL, TOFVINTEGRAND_MSGENULL);
   ASSERT (params, status, TOFVINTEGRAND_ENULL, TOFVINTEGRAND_MSGENULL);



   ak = (TofVIntegrandIn *) params;

   *integrand = ak->dEnergy(v, ak->coeffs)/ak->flux(v, ak->coeffs);

   RETURN (status);


}
