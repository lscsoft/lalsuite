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

typedef struct SingleSplitULDatatag {
  REAL8 *conf;
  REAL8 est_h0;
} SingleSplitULData;

typedef struct SingleULDatatag {
  REAL8 f_min;
  REAL8 f_max;
  REAL8 p_loudestsig;
  REAL8 s_loudestsig;
  REAL8 co_loudestsig;
  INT4 N_per_split;
  SingleSplitULData *SSULdata;
  REAL8 *final_conf;
  REAL8 *final_conf_err;
  REAL8 *h0;
  REAL8 final_h0;
  REAL8 final_h0_err;
} SingleULData;

typedef struct ULDatatag {
  INT4 N_UL;
  INT4 N_h0;
  INT4 N_split;
  SingleULData *UL;
} ULData;

typedef struct Loudesttag {
  REAL8 co_sig;
  REAL8 p_sig;
  REAL8 s_sig;
} Loudest;
