/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>

#include <lal/Matrix.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/Velocity.h>
#include <lal/Units.h>
#include <lal/VectorOps.h>
#include <lal/SkyCoordinates.h>

#include <lal/PrintFTSeries.h>
#include <lal/StreamOutput.h>

#include "config.h"

/* globals */
extern int lalDebugLevel;

int main(int argc, char **argv)
{
  static LALStatus status;
  REAL4Array      *skygrid = NULL;
  UINT4            num_ra, num_dec;
  UINT4Vector     *grid_dims = NULL;
  UINT4            i, j;
  
  num_ra  = 4;
  num_dec = 2;
  
  LALU4CreateVector(&status, &grid_dims, 2);
  grid_dims->data[0] = num_ra;
  grid_dims->data[1] = num_dec;
  
  printf("%d\n", grid_dims->length);
  printf("%d\n", grid_dims->data[0]);
  printf("%d\n", grid_dims->data[1]);

  LALSCreateArray(&status, &skygrid, grid_dims);
  
  printf("%d\n", skygrid->dimLength->length);
  printf("%d\n", skygrid->dimLength->data[0]);
  printf("%d\n", skygrid->dimLength->data[1]);
  
  LALSDestroyArray(&status, &skygrid);
  LALU4DestroyVector(&status, &grid_dims);
  
  return 0;
}
