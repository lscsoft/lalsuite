#include <math.h>
#include "LALStdlib.h"
#include "Inspiral.h"

NRCSID (PHIOFVINTEGRANDC, "$Id$"); 

void PhiofVIntegrand (Status *status,
		      REAL8 *integrand,
		      REAL8 v,
		      void *params)
{
  PhiofVIntegrandIn *in;

  INITSTATUS (status, "PhiofVIntegrand", PHIOFVINTEGRANDC);

  ASSERT (integrand, status, PHIOFVINTEGRAND_ENULL, PHIOFVINTEGRAND_MSGENULL);
  ASSERT (params, status, PHIOFVINTEGRAND_ENULL, PHIOFVINTEGRAND_MSGENULL);

  in = (PhiofVIntegrandIn *) params;

  *integrand = pow (v, 3.) * in->dEnergy(v,in->coeffs)/in->flux(v,in->coeffs);

  RETURN(status);


}

