/*----------------------------------------------------------------------- 
 * 
 * File Name: findchirp.h
 *
 * Author: Brown, D. A.
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
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpEngine.h>

#include "lalapps.h"

void
findchirp_usage (
    const char *program, 
    int         exitcode
    );

void
findchirp_parse_options (
    int         argc, 
    char       *argv[]
    );

void
findchirp_print_options ( 
    void 
    );
