#include "Inspiral.h"
#include "LALStdlib.h"

NRCSID (INSPIRALDERIVATIVESC, "$Id$");

void LALInspiralDerivatives (LALStatus *status,
			  REAL8Vector *values,
			  REAL8Vector *dvalues,
			  void *params)
 {

  InspiralDerivativesIn *ak;

  INITSTATUS(status, "LALInspiralDerivatives", INSPIRALDERIVATIVESC);

  ASSERT (values, status, INSPIRALDERIVATIVES_ENULL, INSPIRALDERIVATIVES_MSGENULL);
  ASSERT (values->data, status, INSPIRALDERIVATIVES_ENULL, INSPIRALDERIVATIVES_MSGENULL);
  ASSERT (dvalues, status, INSPIRALDERIVATIVES_ENULL, INSPIRALDERIVATIVES_MSGENULL);
  ASSERT (dvalues->data, status, INSPIRALDERIVATIVES_ENULL, INSPIRALDERIVATIVES_MSGENULL);
  ASSERT (params, status, INSPIRALDERIVATIVES_ENULL, INSPIRALDERIVATIVES_MSGENULL);

  ak = (InspiralDerivativesIn *) params;


   *(dvalues->data) = -ak->flux(*(values->data), ak->coeffs)/(ak->totalmass*ak->dEnergy(*(values->data), ak->coeffs));
   *(dvalues->data+1) = 2.*pow(*(values->data),3.)/ak->totalmass;

   RETURN (status);



}
