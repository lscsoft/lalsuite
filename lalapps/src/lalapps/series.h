#ifndef SERIES_H_
#define SERIES_H_

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( SERIESH, "$Id$" );

typedef enum { Time, Freq, Trans } domain;

#define IS_TIME( domain_ ) ( domain_ == Time )
#define IS_FREQ( domain_ ) ( ( domain_ == Freq ) || ( domain_ == Trans ) )
#define IS_TRANS( domain_ ) ( domain_ == Trans )

struct series
{
  const char *name;
  LIGOTimeGPS tbeg;
  LIGOTimeGPS tend;
  domain      dom;
  int         type;
  double      step;
  const char *unit;
  size_t      size;
  float      *data;
};

double epoch_diff( const LIGOTimeGPS *t2, const LIGOTimeGPS *t1 );
void epoch_add( LIGOTimeGPS *t1, LIGOTimeGPS *t0, double dt );
int write_ilwd( const char *fname, const struct series *ser );
struct FrameH *fr_add_proc_data( struct FrameH *frame, const struct series *ser );

#ifdef  __cplusplus
}
#endif

#endif
