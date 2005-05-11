#ifndef _GENRANDINPUT_H
#define _GENRANDINPUT_H

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
#include <getopt.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>

typedef struct REAL8range {
  REAL8 min;
  REAL8 max;
  REAL8 value;
} REAL8range;

typedef struct INT4range {
  INT4 min;
  INT4 max;
  INT4 value;
} INT4range;

typedef struct LIGOTimeGPSrange {
  LIGOTimeGPS min;
  LIGOTimeGPS max;
  LIGOTimeGPS value;
} LIGOTimeGPSrange;

typedef struct RandomParameterstag {
  INT4 start;
  char detector[256];
  char stampsfile[256];
  char noisedir[256];
  char sftbase[256];
  REAL8range *ra;
  REAL8range *dec;
  REAL8range *sma;
  LIGOTimeGPSrange *tperi;
  REAL8range *period;
  REAL8range *ecc;
  REAL8range *argp;
  REAL8range *freq;
  REAL8range *phi;
  REAL8range *psi;
  REAL8range *cosiota;
  REAL8range *h0;
  INT4range *det;
} RandomParameters;


#endif
