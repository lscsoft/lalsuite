#include <lal/LALDatatypes.h>

typedef enum { Time, Freq } domain;

#define IS_TIME( domain_ ) ( domain_ == Time )
#define IS_FREQ( domain_ ) ( domain_ == Freq )

struct series
{
  const char *name;
  LIGOTimeGPS tbeg;
  LIGOTimeGPS tend;
  domain      dom;
  int         type;
  float       step;
  const char *unit;
  size_t      size;
  float      *data;
};

double epoch_diff( const LIGOTimeGPS *t2, const LIGOTimeGPS *t1 );
void epoch_add( LIGOTimeGPS *t1, LIGOTimeGPS *t0, double dt );
int write_ilwd( const char *fname, const struct series *ser );
struct FrameH *fr_add_proc_data( struct FrameH *frame, const struct series *ser );
