#include <stdio.h>
#include "lal/LALRCSID.h"

NRCSID (MAIN, "$Id$");

#include <lal/TFCThresholds.h>
#define CHKST if(status.statusCode != 0) return -1

int lalDebugLevel = LALMSGLVL3;

int main(void) {

  static LALStatus status;
  
  UINT4 i, nFreq;
  REAL4 *rho;
  RiceThresholdParams params;

  nFreq = 16;

  rho = (REAL4 *)LALMalloc(nFreq * sizeof(REAL4));
  if(!rho) return 1;

  params.nFreq = nFreq;
  params.bpp = 0.1;
  params.eGoal = 1e-3;

  params.P0 = (REAL8 *)LALMalloc(nFreq * sizeof(REAL8));
  if(!params.P0) return 1;

  params.Q = (REAL8 *)LALMalloc(nFreq * sizeof(REAL8));
  if(!params.Q) return 1;

  
  /* set dummy values */
  for(i=0;i<nFreq;i++) {
    params.P0[i] = (REAL4)(i+1);
    params.Q[i] = 0.25 * (REAL4)(1+i*i);
  }

  LALTFCRiceThreshold(&status, rho, &params);

  if(status.statusCode != 0) return 1;

  for(i=0;i<nFreq;i++) printf("%g\t%g\t->\t%g\n",params.P0[i], params.Q[i],rho[i]);
  
  return 0;
}
