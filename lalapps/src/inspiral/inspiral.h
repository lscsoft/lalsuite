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

#ifndef INSPIRAL_H_
#define INSPIRAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <FrameL.h>

#include <lalapps.h>
#include <series.h>
#include <processtable.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/FrameCalibration.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpEngine.h>

#define PROGRAM_NAME "inspiral"

#define USAGE \
"lalapps_inspiral is a stand alone code for performing matched filtering\n" \
"of LIGO data for graviational wave signals and Monte Carlo analysis.\n" \
"\n" \
"Usage: lalapps_inspiral [OPTIONS]\n" \
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
"INPUT DATA OPTIONS.\n" \
"\n" \
"The following options contol reading of the raw data that the code will filter.\n" \
"The inspiral code reads in a block of data from the desired channel, starting" \
"at the specified GPS start time and ending at the last sample before the GPS\n" \
"end time (i.e. start <= data < end). It breaks this into a specified number of \n" \
"segments of a given length and overlap for filtering. Sanity checking is \n" \
"performed to make sure that the specified block of data can be appropriately\n" \
"segmented.\n" \
"\n" \
"The raw data is read from IGWD frame files and the template bank is read from a\n" \
"LIGO lightweight XML file containing the template parameters.\n" \
"\n" \
"   --gps-start-time             GPS start of data to be filtered\n" \
"   --gps-stop-time              GPS stop time of data to be filtered\n" \
"   --channel-name lstring       name of the channel to filter (e.g. H1:LSC-AS_Q)\n" \
"   --segment-length int_4u      length of each individual segment\n" \
"   --number-of-segments int_4u  numer of segments to break the data into\n" \
"   --segment-overlap int_4u     overlap between consecutive segments\n" \
"   --start-template int_4u      start parsing templates at this row number\n" \
"   --stop-template int_4u       stop parsing templates at this row number\n" \
"\n" \
"\n" \
"DATA CONDITIONING OPTIONS.\n" \
"\n" \
"The following options control any pre-conditioning applied to the raw data\n"\
"before it is passed to the filtering code as well as parameters to generate\n"\
"the one-sided power specral density of the data. It is possible to resample\n"\
"and/or apply a time-domain high-pass filter to the data before the power\n" \
"spectrum is computed\n" \
"\n" \
"   --sample-rate int_4u            resample the raw data to the specified rate\n" \
"   --enable-high-pass real_4       high pass the data above the specified\n"\
"                                      frequency\n" \
"   --disable-high-pass             do not high pass the data\n" \
"   --low-frequency-cutoff real_4   specifiy a low frequency cutoff for the\n"\
"                                      filter. All frequencies below this value\n"\
"                                      are ignored. This takes effect after\n"\
"                                      the high-pass filter.\n" \
"   --spectrum-type lstring         method of psd computation (mean or median)\n" \
"   --inverse-spec-length int_4u    truncate the inverse power spectal density in\n" \
"                                   time domain. Specify desired length of the\n"\
"                                   inverse specrum in seconds.\n" \
"   --dynamic-range-exponent real_4 to keep the data in range a dynamic range\n"\
"                                      scale factor is applied to the data as\n" \
"                                      described in the findchirp documentation.\n"\
"   --calibration string            location of a frame cache file containing\n"\
"                                      calibration information.\n" \
"\n" \
"\n" \
"FILTERING OPTIONS.\n" \
"\n" \
"The following options are using by the filtering routines to control the\n" \
"parameters of the matched filter and chi squared veto. The thresholds for\n" \
"signal to noise and chi squared are specified as a comma separated list\n" \
"(e.g. 9.0,7.0,6.0) with the first value being used for the coarse search,\n" \
"the second value for the first hierarchical search and so on. If the\n" \
"numer of chi squared bins is set to 0 then the chi squared veto is disabled.\n" \
"\n" \
"   --chisq-bins int_4u         number chi squared veto frequency bins\n" \
"   --rhosq-thresholds real_4   list of signal to noise thresholds\n" \
"   --chisq-thresholds real_4   list of chi squared thresholds\n" \
"   --enable-event-cluster      enable the event clustering algorithm\n" \
"   --disable-event-cluster     disable the event clustering algorithm\n" \
"\n" \
"\n" \
"OUTPUT OPTIONS.\n" \
"\n" \
"The following options control writing any result data to disk. Results are\n" \
"written as LIGO lightweigh XML file (LIGOLw) or as IGWD frame files (frames)\n" \
"\n" \
"   --enable-output             write inspiral triggers found as LIGOLw\n" \
"   --disable-output            do not write triggers as LIGOLw \n" \
"   --comment lstring           add a comment to the process params table\n" \
"   --write-raw-data            write the raw data read in as a frame\n" \
"   --write-filter-data         write the pre-conditioned data as a frame\n" \
"   --write-response            write the response function used as a frame\n" \
"   --write-spectrum            write the power spectrum used as a frame\n" \
"   --write-rhosq               write the last vector of snr^2 data as a frame\n" \
"   --write-chisq               write the last vector of chi^2 data as a frame\n" \
"\n" \
"\n" \
"EXAMPLE USEAGE.\n" \
"\n" \
"Requires 256 seconds of LHO frame data in the current directory.\n" \
"\n" \
"lalapps_inspiral --enable-output --enable-high-pass 50.0\n" \
"                 --enable-event-cluster --gps-start-time 600000000\n" \
"                 --gps-end-time 600000256 --segment-length 262144\n" \
"                 --number-of-segments 7 --segment-overlap 131072\n" \
"                 --sample-rate 4096 --chisq-bins 8\n" \
"                 --low-frequency-cutoff 100.0 --spectrum-type mean\n" \
"                 --inverse-spec-length 0 --dynamic-range-exponent 69.0\n" \
"                 --channel-name H1:LSC-AS_Q --rhosq-threshold 1.0\n" \
"                 --chisq-threshold 0.001 --verbose --debug-level NDEBUG\n" \
"\n"

#define BANK_FILE "tmpltbank.xml"

int snprintf(char *str, size_t size, const  char  *format, ...);
char *strsep(char **stringp, const char *delim);
int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#endif
