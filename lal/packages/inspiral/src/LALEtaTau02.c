#include <lal/LALInspiral.h>

void LALEtaTau02(LALStatus *status, REAL8 *x, REAL8 eta, void *p) {
   EtaTau02In *q;
/* 
   In the notation of June 4, pp 1 
   Given tau0 and tau2 determine eta by solving EtaTau02 = 0.
   p[0] = t2
   p[1] = A2 * pow(t0/A0, 3/5)
   p[2] = B2
*/
   status = NULL;
   q = (EtaTau02In *) p;
   *x = -q->t2 + q->A2/pow(eta,0.4) * (1. + q->B2*eta);
}
