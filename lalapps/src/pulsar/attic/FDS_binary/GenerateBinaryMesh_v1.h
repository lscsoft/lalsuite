/*
*  Copyright (C) 2007 Chris Messenger, Reinhard Prix
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _GENBINMESH_H
#define _GENBINMESH_H

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
#include <lal/Random.h>
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
#include <glob.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>

#define DIM 2         /* Current hardcoded number of dimensions in parameter space  */
#define ECC 0.0       /* circular orbit currently harcoded in */
#define ARGP 0.0
#define MAXSFTS 50000

typedef struct CLargstag {
  CHAR sourcefile[256];
  CHAR source[256];
  CHAR datadir[256];
  LIGOTimeGPS tstart;
  REAL8 tspan;
  REAL8 mismatch;
  REAL8 band;
  CHAR ephemdir[256];
  CHAR yr[256];
  CHAR ifo[256];
  CHAR meshdir[256];
  BOOLEAN datadirflag;
  BOOLEAN mismatchedflag;
  BOOLEAN exactflag;
} CLargs;

typedef struct GlobVartag {
  INT4 nband;
  REAL8 *f_max;
  INT4 band;
  REAL8 tspan;
  LIGOTimeGPS tstart;
  LIGOTimeGPS tstartSSB;
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
  CHAR ifo[256]; 
  CHAR ephemdir[256];
  CHAR yr[256];
  CHAR meshdir[256];
  CHAR sourcefile[256];
  CHAR source[256];
  BOOLEAN mismatchedflag;
  BOOLEAN exactflag;
} GlobVar;
  
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

typedef struct BinaryTemplatetag {             /* BINARY-MOD - structure to store a single binary signal template */
  REAL8       ProjSMaxis;
  REAL8       Period;
  LIGOTimeGPS TperiSSB;
  REAL8       Eccentricity;
  REAL8       ArgPeri;
} BinaryTemplate;

typedef struct BinaryTemplateBanktag {
  BinaryMeshFileHeader BMFheader;
  BinaryTemplate *BTB;       
} BinaryTemplateBank;

/*struct sftheadertag {
  REAL8 endian;
  INT4 gps_sec;
  INT4 gps_nsec;
  REAL8 tbase;
  INT4 firstfreqindex;
  INT4 nsamples;
  } sftheader;*/

int ConvertXYtoRTperi(XYLocation *, RTPLocation *);
int WriteMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader);
int ReadMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader);
int ConvertRTperitoXY(RTPLocation *, XYLocation *, REAL8 *);
int PeriapseShift(LIGOTimeGPS, LIGOTimeGPS *,LIGOTimeGPS, REAL8,INT4 *);
int PeriapseShiftBack(LIGOTimeGPS, LIGOTimeGPS,LIGOTimeGPS,LIGOTimeGPS *, REAL8,INT4);

#endif /* end the if over GENBINMESH_H */

