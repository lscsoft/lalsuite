#include <stdio.h>
#include <lal/LALDatatypes.h>

typedef struct tagLALDataFileNameFields
{
  CHAR  site[FILENAME_MAX];
  CHAR  description[FILENAME_MAX];
  INT4  tstart;
  INT4  duration;
  CHAR  extension[FILENAME_MAX];
}
LALDataFileNameFields;

typedef struct tagLALCalRefFileNameDescriptionFields
{
  CHAR ifo[FILENAME_MAX];
  CHAR channelPostfix[FILENAME_MAX];
  CHAR run[FILENAME_MAX];
  INT4 version;
}
LALCalRefFileNameDescriptionFields;

typedef struct tagLALCalFacFileNameDescriptionFields
{
  CHAR ifo[FILENAME_MAX];
  CHAR run[FILENAME_MAX];
  INT4 version;
  INT4 deltaT;
}
LALCalFacFileNameDescriptionFields;

int XLALDataFileNameParse( LALDataFileNameFields *fields, const char *fname );

int XLALCalRefFileNameDescriptionParse( LALCalRefFileNameDescriptionFields *fields, const char *description );
int XLALCalFacFileNameDescriptionParse( LALCalFacFileNameDescriptionFields *fields, const char *description );

int XLALASCIIFileCountRows( const char *fname );
REAL8VectorSequence * XLALASCIIFileReadColumns( INT4 ncol, const char *fname );

REAL4 XLALASCIIFileReadCalFacHeader( const char *fname );
REAL4 XLALASCIIFileReadCalRefHeader( const char *fname );
int XLALASCIIFileReadCalFac( REAL4TimeSeries **alpha, REAL4TimeSeries **gamma, const char *fname );
int XLALASCIIFileReadCalRef( COMPLEX8FrequencySeries **series, REAL8 *duration, const char *fname );
