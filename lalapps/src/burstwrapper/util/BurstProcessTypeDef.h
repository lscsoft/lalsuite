#ifndef BurstProcessTypeDefh
#define BurstProcessTypeDefh

#include <lal/LALDatatypes.h>

typedef struct tagSnglBurstTableC {

  struct tagSnglBurstTableC *next;

  LIGOTimeGPS start_time;
  REAL4         duration;
  REAL4         central_freq;
  REAL4         bandwidth;
  REAL4         amplitude;
  REAL4         snr;
  REAL4         confidence;

} SnglBurstTableC;



#endif
