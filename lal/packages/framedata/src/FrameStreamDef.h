#ifndef _FRAMESTREAMDEF_H
#define _FRAMESTREAMDEF_H
#include <FrameL.h>
#include <lal/LALDatatypes.h>

/* Useful macros */
#define SECNAN_TO_I8TIME( sec, nan ) \
  ((INT8)1000000000*(INT8)(sec)+(INT8)(nan))
/* Dangerous!!!: */
#define EPOCH_TO_I8TIME( epoch ) \
  SECNAN_TO_I8TIME( (epoch).gpsSeconds, (epoch).gpsNanoSeconds )
#define SET_EPOCH( pepoch, i8time ) \
  do { INT8 t=(i8time); LIGOTimeGPS *pe=(pepoch); \
    pe->gpsSeconds=t/(INT8)1000000000; pe->gpsNanoSeconds=t%(INT8)1000000000; \
  } while( 0 )

typedef struct
tagFrFileInfo
{
  INT4  ind;
  CHAR *url;
  INT4  t0;
  INT4  dt;
}
FrFileInfo;

/* Definition of FrStream */
struct
tagFrStream
{
  FrFileInfo     *filelist;
  UINT4           numfiles;
  UINT4           filenum;
  struct FrFile  *frfile;
  struct FrameH  *frame;
  LIGOTimeGPS     epoch;
  INT4            end;
  INT4            err;
  INT4            gap;
};
#endif /* _FRAMESTREAMDEF_H */
