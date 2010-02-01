

#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "TwoSpect.h"

typedef struct
{
   REAL4 fsig; /* 0 value means candidate not valid */
   REAL4 period;
   REAL4 moddepth;
   REAL4 Tobs;
   REAL4 Tcoh;
   REAL4 fmin;
   REAL4 fspan;
   REAL4 stat;
   REAL4 snr;
} candidate;


candidate * new_candidate(void);
void free_candidate(candidate *cand);
void loadCandidateData(candidate *out, REAL4 fsig, REAL4 period, REAL4 moddepth, REAL4 Tobs, REAL4 Tcoh, REAL4 fmin, REAL4 fspan, REAL4 stat, REAL4 snr);
void clusterCandidates(candidate *out[], candidate *in[], ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, INT4 numofcandidates, INT4 option);





#endif


