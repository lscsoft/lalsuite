#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <glob.h>
#include <getopt.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/StackMetric.h>
#include <lal/MatrixUtils.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>

#define MAXFILES 50000

typedef struct Cornertag {
  REAL8 x;
  REAL8 y;
} Corner;

typedef struct Linetag {
  REAL8 minangle;
  REAL8 maxangle;
  REAL8 gradient;
  REAL8 intercept;
} Line;

typedef struct Resulttag {
  REAL8 freq;
  REAL8 RA;
  REAL8 dec;
  REAL8 sma;
  REAL8 period;
  LIGOTimeGPS tp;
  REAL8 ecc;
  REAL8 argp;
  INT4 ncluster;
  REAL8 meantwoF;
  REAL8 stdtwoF;
  REAL8 twoF;
} Result;

typedef struct Resultstag {
  INT4 Nresults;
  Result *result;
} Results;

typedef struct Significancetag {
  REAL8 log10oneminusp_sig;
  REAL8 log10oneminuss_sig;
  REAL8 log10oneminusco_sig;
} Significance;

typedef struct CoResultstag {
  INT4 Nresults;
  Result *primary_result;
  Result *secondary_result;
  Significance *significance;
} CoResults;

typedef struct FreqMeshtag {
  BinaryMeshFileHeader p_header;
  BinaryMeshFileHeader s_header;
  REAL8 f_min;
  REAL8 f_max;
  REAL8 f_band;
} FreqMesh;

typedef struct FreqMeshestag {
  INT4 Nheaders;
  FreqMesh *freqmesh;
} FreqMeshes;
