/*----------------------------------------------------------------------- 
 * 
 * File Name: olapredfcn.h
 *
 * Author: John T. Whelan
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/DetectorSite.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>

#include "lalapps.h"

#define OLAPREDFCNH_TRUE     1
#define OLAPREDFCNH_FALSE    0
#define OLAPREDFCNH_OOR    1e+5

/***************************** <lalErrTable file="olapredfcnHE"> */
#define OLAPREDFCNH_ENOM 0
#define OLAPREDFCNH_EARG 1
#define OLAPREDFCNH_ESUB 2
#define OLAPREDFCNH_MSGENOM "Nominal exit"
#define OLAPREDFCNH_MSGEARG "Invalid command-line arguments"
#define OLAPREDFCNH_MSGECHK "LAL Subroutine Returned Error"
/***************************** </lalErrTable> */

enum
{
  LALBarIndexAURIGA,
  LALBarIndexNAUTILUS,
  LALBarIndexEXPLORER,
  LALBarIndexALLEGROOLD,
  LALBarIndexNIOBE,
  LALBarIndexALLEGRONEW,
  LALNumCachedBars
};

extern const LALFrDetector lalCachedBars[LALNumCachedBars];

void
olapredfcn_usage (
    const char *program, 
    int         exitcode
    );

void
olapredfcn_parse_options (
    int         argc, 
    char       *argv[]
    );

void
olapredfcn_print_options ( 
    void 
    );
