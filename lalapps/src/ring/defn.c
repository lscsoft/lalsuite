#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <glob.h>
#include <lalapps.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>

#include "defn.h"

REAL4 ampl( REAL4 M, REAL4 Q, REAL4 a, REAL4 eps, REAL4 D ){
  REAL4 h0;
  h0  = sqrt( 5.0 / 2.0 * eps )
    * ( LAL_G_SI * M * LAL_MSUN_SI / pow( LAL_C_SI,2.0) /
        D / LAL_PC_SI /1000000.0 ) * pow( Q, -0.5 )
    * pow( 1.0 + 7.0 / 24.0 / pow( Q, 2.0), -0.5 )
    * pow(  1.0 - 0.63 * pow( 1.0 - a,0.3 ), -0.5);

  return h0;
}


