/*----------------------------------------------------------------------- 
 * 
 * File Name: inspiral.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lalapps.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>

#define PROGRAM_NAME "lalapps_inspiral"

#define USAGE \
"lalapps_inspiral is a stand alone code for performing matched filtering\n" \
"of LIGO data for graviational wave signals and Monte Carlo analysis.\n\n" \
"Usage: lalapps_inspiral [OPTIONS]\n" \
"Program options:\n" \
"   --verbose                   verbose operation\n" \
"   --debug                     trun on debugging information\n\n" \
"Input data options:\n" \
"   --gps-start-time            GPS start of data to be filtered\n" \
"   --gps-stop-time             GPS stop time of data to be filtered\n"

