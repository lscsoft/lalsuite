#ifndef COINCIDENCESH
#define COINCIDENCESH

/*
#include <lal/EPSearch.h>
*/

#include <BurstProcessTypeDef.h>

typedef struct tagCoincidence2Params {

  double dt;
  double df;
  double dbw;
  double ddur;
  double damp;
  double dconf;
  double dsnr;

} Coincidence2Params;

int Get2Coincidences(
		     SnglBurstTableC *output,
		     SnglBurstTableC *input1,
		     SnglBurstTableC *input2,
		     Coincidence2Params *cparams
		     );
		     

#endif
