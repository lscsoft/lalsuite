#include <stdio.h>
#include "lal/LALRCSID.h"

NRCSID (MAIN, "$Id$");

#include <lal/TFCThresholds.h>
#define CHKST if(status.statusCode != 0) return -1

int lalDebugLevel = 0;

int main(int argc, char* argv[]) {

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

  params.meanRe = (REAL4 *)LALMalloc(nFreq * sizeof(REAL4));
  if(!params.meanRe) return 1;

  params.meanIm = (REAL4 *)LALMalloc(nFreq * sizeof(REAL4));
  if(!params.meanIm) return 1;

  params.varRe = (REAL4 *)LALMalloc(nFreq * sizeof(REAL4));
  if(!params.varRe) return 1;

  params.varIm = (REAL4 *)LALMalloc(nFreq * sizeof(REAL4));
  if(!params.varIm) return 1;

  
  /* set dummy values */
  for(i=0;i<nFreq;i++) {
    params.meanRe[i] = (REAL4)(i+1);
    params.meanIm[i] = (REAL4)(nFreq-i);
    params.varRe[i] = 0.25 * (REAL4)(1+i*i);
    params.varIm[i] = 0.25 * (REAL4)(nFreq*nFreq+i*i);
  }

  LALTFCRiceThreshold(&status, rho, &params);

  if(status.statusCode != 0) return 1;

  for(i=0;i<nFreq;i++) printf("%g\t%g\t%g\t%g\t->\t%g\n",params.meanRe[i], params.meanIm[i],params.varRe[i],params.varIm[i],rho[i]);
  
  return 0;
}
