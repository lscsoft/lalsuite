/*----------------------------------------------------------------------- 
 * 
 * File Name: tmpltbank.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef TMPLTBANK_H_
#define TMPLTBANK_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <errno.h>
#include <lalapps.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>

#define PROGRAM_NAME "tmpltbank"

#define USAGE \
"go fish"

int snprintf(char *str, size_t size, const  char  *format, ...);
int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#endif
