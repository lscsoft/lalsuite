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
#include <lal/FlatMesh.h>
#include <lal/FindRoot.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/StackMetric.h>
#include <lal/MatrixUtils.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>


#define DIM 2         /* Current hardcoded number of dimensions in parameter space  */
#define ECC 0.0       /* circular orbit currently harcoded in */
#define ARGP 0.0

typedef struct CLargstag {
  REAL8 f_max;
  REAL8 tspan;
  LIGOTimeGPS tstart;
  REAL8 RA;
  REAL8 dec;
  REAL8 period;
  REAL8 sma_0;
  LIGOTimeGPS tperi_0;
  REAL8 sma_MIN;
  REAL8 sma_MAX;
  LIGOTimeGPS tperi_MIN;
  LIGOTimeGPS tperi_MAX;
  REAL8 mismatch;
  CHAR ephemdir[256];
  CHAR yr[256];
  CHAR ifo[256];
  CHAR meshfile[256];
} CLargs;
  
typedef struct BinaryMeshFileHeadertag {
  REAL8 f_max;
  REAL8 tspan;
  LIGOTimeGPS tstart;
  UINT4 Nfilters;
  REAL8 mismatch;
  REAL8 sma_0;
  REAL8 sma_MIN;
  REAL8 sma_MAX;
  LIGOTimeGPS tperi_0;
  LIGOTimeGPS tperi_MIN;
  LIGOTimeGPS tperi_MAX;
  REAL8 ecc_MIN;
  REAL8 ecc_MAX;
  REAL8 argp_MIN;
  REAL8 argp_MAX;
  REAL8 period_MIN;
  REAL8 period_MAX;
  REAL8 metric_XX;
  REAL8 metric_XY;
  REAL8 metric_YY;
  CHAR version[256];
  CHAR det[256];
  REAL8 RA;
  REAL8 dec;
} BinaryMeshFileHeader;

typedef struct XYparameterspacetag {
  REAL8 X_MIN;
  REAL8 X_MAX;
  REAL8 X_0;
  REAL8 Y_MIN;
  REAL8 Y_MAX;
  REAL8 Y_0;
} XYparameterspace;
  
typedef struct RTparameterspacetag {
  REAL8 sma_MIN;
  REAL8 sma_MAX;
  REAL8 sma_0;
  LIGOTimeGPS tperi_MIN;
  LIGOTimeGPS tperi_MAX;
  LIGOTimeGPS tperi_0;
} RTparameterspace;

typedef struct Metrictag {
  REAL8 X_0;
  REAL8 Y_0;
  REAL8 determinant;
  REAL8Array *element;
  REAL8Array *eigenvec;
  REAL8Vector *eigenval;
} Metric;

typedef struct RTMeshtag {
  UINT4 length;
  REAL8Vector *sma;
  LIGOTimeGPS *tperi;
} RTMesh;

typedef struct RTPLocationtag {
  REAL8 sma;
  REAL8 period;
  LIGOTimeGPS tperi;
  LIGOTimeGPS tstartSSB;
  REAL8 ecc;
  REAL8 argp;
} RTPLocation;

typedef struct XYLocationtag {
  REAL8 X;
  REAL8 Y;
  REAL8 period;
  LIGOTimeGPS tstartSSB;
  REAL8 ecc;
  REAL8 argp;
} XYLocation;

int ConvertXYtoRTperi(XYLocation *, RTPLocation *);
int WriteMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader);
int ReadMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader);
int ConvertRTperitoXY(RTPLocation *, XYLocation *, REAL8 *);



