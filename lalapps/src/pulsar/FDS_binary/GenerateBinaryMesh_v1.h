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
#include <sys/stat.h>
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
  REAL8 fmax;
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

int FreeMem();
int GenerateMesh(REAL4VectorSequence **,XYparameterspace *,Metric *);
int SetupPspaceParams(RTparameterspace *,XYparameterspace *);
int GenMetricComp(REAL8 *,REAL8 *,Metric *);
int CheckRTBoundary(REAL8 *,LIGOTimeGPS *,RTparameterspace *);
int ConvertMesh(REAL4VectorSequence **,RTMesh *,RTparameterspace *);
int OutputRTMesh(RTMesh *,Metric *, RTparameterspace *);
int ReadCommandLine(int argc,char *argv[]);
int ConvertXYtoRTperi(REAL8 *, REAL8 *, REAL8 *, LIGOTimeGPS *);
int ConvertTperitoPhase();
int SetupBaryInput();
int CheckInput();
int GetSSBTime(LIGOTimeGPS *, LIGOTimeGPS *);
int ConvertRTperitoXY(REAL8 *,LIGOTimeGPS *,REAL8 *,REAL8 *,REAL8 *);
static void OrbPhaseFunc(LALStatus *, REAL8 *, REAL8, void *);
