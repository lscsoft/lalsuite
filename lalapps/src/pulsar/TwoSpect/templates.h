

#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__

#include <lal/RealFFT.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "candidates.h"
#include "TwoSpect.h"

typedef struct
{
   REAL4 far;
   REAL4 distMean;
   REAL4 distSigma;
} farStruct;


farStruct * new_farStruct(void);
void free_farStruct(farStruct *farstruct);
void estimateFAR(farStruct *out, REAL4Vector *weights, topbinsStruct *topbinsstruct, REAL4 thresh, REAL4Vector *ffplanenoise);

void makeTemplateGaussians(ffdataStruct *out, candidate *in);
void makeTemplate(ffdataStruct *out, candidate *in, REAL4FFTPlan *plan);




#endif

