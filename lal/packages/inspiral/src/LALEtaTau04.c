#include <lal/LALInspiral.h>

void LALEtaTau04(LALStatus *status, REAL8 *x, REAL8 eta, void *p) {
   EtaTau04In *q;
/* 
   In the notation of June 4, pp 1 
   Given tau0 and tau4 determine eta by solving EtaTau04 = 0.
   p[0] = t4
   p[1] = A4 * pow(t0/A0, 1/5)
   p[2] = B4
   p[3] = C4
*/
   status = NULL;
   q = (EtaTau04In *) p;
   *x = -q->t4 + q->A4/pow(eta,0.8) * (1. + q->B4*eta + q->C4*eta*eta);
}
