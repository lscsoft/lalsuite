/*
*  Copyright (C) 2007 Chris Messenger
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
#include "ComputeFStatisticBinary_v2.h" 
#include "ReadSourceFile_v1.h"

typedef struct datasettag {
  INT4 tsft;
  REAL8 tobs;
  INT4 nbins;
  INT4 nsamples;
  INT4 sftno;
  LIGOTimeGPS *stamps;
  LIGOTimeGPS start;
  LIGOTimeGPS end;
} dataset;

typedef struct sensitivityparamstag {
  INT4 overlap;
  REAL8 tspan;
  INT4 nchunks;
  LIGOTimeGPS fullstart;
  LIGOTimeGPS fullend;
  dataset *dataparams;
} sensitivityparams;

typedef struct sensresultstag {
  REAL8 A;
  REAL8 B;
  REAL8 ShAV[2];
  REAL8 ShAVweight;
  REAL8 Q;
} sensresults;


