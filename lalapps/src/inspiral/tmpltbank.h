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
"lalapps_tmpltbank is a stand alone code for generating inspiral template\n" \
"banks for LIGO data with the LAL bank package.\n" \
"\n" \
"The code generates a calibrated power spectrum at the specified time for\n" \
"the requested channel and uses this to compute the template bank.\n" \
"\n" \
"See the LAL bank package documentation for detailed information on the\n" \
"algorithms used to generate the banks.\n" \
"\n" \
"Usage: lalapps_tmpltbank [OPTIONS]\n" \
"\n" \
"If the option takes an argument, the agument type is shown after the option.\n" \
"Some of the options are required for operation. If a required option is\n" \
"not specified, the program will exit with an error message reporting the\n" \
"missing option. Sanity checking is performed on option arguments.\n"\
"\n" \
"\n" \
"PROGRAM OPTIONS.\n" \
"\n" \
"The following options control non-sceintific program behaviour, such as\n" \
"enabling debugging information:\n"\
"\n" \
"   --help                      display this message\n" \
"   --debug-level lstring       set the LAL debug level to the specified\n" \
"                                 value. Useful values are: NDEBUG, ERROR,\n" \
"                                 WARNING, INFO, TRACE, MEMINFO and MEMDBG\n" \
"   --verbose                   verbose operation\n" \
"\n" \
"\n" \
"POWER SPECTRUM GENERATION OPTIONS.\n" \
"\n" \
"The following options contol reading of the raw data and generation of \n" \
"the power spectum. The raw data is read from IGWD frame files.\n" \
"\n" \
"   --gps-start-time                GPS start (seconds) of data in PSD\n" \
"   --gps-stop-time                 GPS stop time (seconds) of data in PSD\n" \
"   --channel-name lstring          name of the channel (e.g. H1:LSC-AS_Q)\n" \
"   --segment-length int_4u         length of each individual segment\n" \
"   --number-of-segments int_4u     numer of segments to break the data into\n" \
"   --sample-rate int_4u            resample the raw data to the specified rate\n" \
"   --enable-high-pass real_4       high pass the data above the specified\n"\
"                                     frequency before computing PSD\n" \
"   --disable-high-pass             do not high pass the data\n" \
"   --spectrum-type lstring         method of psd computation (mean or median)\n" \
"   --frame-cache lstring           name of the cache file containing input data\n" \
"   --calibration-cache lstring     location of a frame cache file containing\n"\
"                                     calibration information.\n" \
"\n" \
"\n" \
"TEMPLATE BANK OPTIONS.\n" \
"\n" \
"The following options are using to generate the bank. A brief overview of the\n" \
"paramerers is given here. The LAL bank package contains should be consulted\n" \
"for the full documentation on the parameters\n" \
"\n" \
"   --minimum-mass real_4           minimum (component) mass in the bank\n" \
"   --maximum-mass real_4           maximum (component) mass in the bank\n" \
"   --minimum-match real_4          minimal match of the template bank\n" \
"   --low-frequency-cutoff real_4   high frequency cutoff\n"\
"   --high-frequency-cutoff real_4  upper frequency cutoff\n"\
"   --order lstring                 post-Newtonian order of the waveform\n" \
"                                    (newtonian|oneHalfPN|onePN|onePointFivePN|\n" \
"                                    twoPN|twoPointFive|threePN|threePointFivePN)\n" \
"   --approximant lstring           approximant of the waveform\n" \
"                                    (TaylorT1|TaylorT2|TaylorT3|\n" \
"                                    TaylorF1|TaylorF2|PadeT1|PadeT2|\n" \
"                                    EOB|BCV|SpinTaylorT3)\n" \
"   --space lstring                 space in which to lay down the bank\n" \
"                                     (Tau0Tau2|Tau0Tau3)\n" \
"\n" \
"\n" \
"OUTPUT OPTIONS.\n" \
"\n" \
"The following options control writing any result data to disk. The bank is\n" \
"written as LIGO lightweigh XML file (LIGOLw), data as IGWD frame files.\n" \
"\n" \
"   --comment lstring           add a comment to the process params table\n" \
"   --write-raw-data            write the raw data read in as a frame\n" \
"   --write-response            write the response function used as a frame\n" \
"   --write-spectrum            write the power spectrum used as a frame\n" \
"\n" \
"\n" \
"EXAMPLE USEAGE.\n" \
"\n" \
"Generates a typical template bank as used in the S1 analysis:\n" \
"\n" \
"lalapps_tmpltbank --gps-start-time 715482185 \n" \
"                  --gps-end-time 715482441 \n" \
"                  --segment-length 262144 --number-of-segments 7 \n" \
"                  --sample-rate 4096 \n" \
"                  --enable-high-pass 30.0 --low-frequency-cutoff 40.0 \n" \
"                  --spectrum-type mean \n" \
"                  --channel-name \"L1:LSC-AS_Q\" \n" \
"                  --frame-cache s1_bank_data.catalog \n" \
"                  --calibration-cache s1_bank_cal.catalog \n" \
"                  --minimum-mass 1.0 --maximum-mass 3.0 --minimum-match 0.97 \n" \
"                  --high-frequency-cutoff 1024.0 \n" \
"                  --order twoPN --approximant TaylorT1 --space Tau0Tau3 \n" \
"                  --write-raw-data --write-spectrum \n" \
"                  --write-response --write-strain-spectrum\n" \
"\n"


int snprintf(char *str, size_t size, const  char  *format, ...);
int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#endif
