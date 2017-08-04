#ifndef _SPLINTER_PULSAR_H
#define _SPLINTER_PULSAR_H

#define _GNU_SOURCE   /* for alphasort() and scandir() */

#include <sys/time.h>

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>

/* LAL headers */
#include <lal/BinaryPulsarTiming.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/FrequencySeries.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/SFTfileIO.h>
#include <lal/LALCache.h>
/* lalapps header */
#include <lalapps.h>

/* Usage format string. */
#define USAGE \
"Usage: lalapps_SplInter [options]\n\n"\
" --help (-h)           display this message\n"\
" \nRequired Input Variables:\n\n"\
" --ifo (-i)            name of ifo e.g. L1, H1, H2, G1\n"\
" --start-freq (-S)     Start frequency of the SFTs\n"\
" --end-freq (-E)       End frequency of the SFTs\n"\
" --ephem-dir (-e)      directory containing the ephemeris files\n"\
" --output-dir (-o)     directory for output data files\n"\
" --seg-file (-l)       name of file containing science segment list\n"\
"\n Plus one of the following:\n\n"\
" --param-file (-P)     name of file containing initial pulsar\n\
                        parameters when using just one file as input\n\
                        (.par file)\n"\
" --param-dir (-d)      name of directory containing pulsar parameter\n\
                        files when using multiple input files (.par\n\
                        files)\n"\
"\n and one of the following:\n\n"\
" --sft-cache (-F)      location of list of SFTs for use, can use\n\
                        list:listsfts.txt or /location/*.sft\n"\
" --sft-loc (-L)        directory of files containing a list of SFTs\n\
                        for each segment, located in\n\
                        /sft-loc/Segment_Start-End.sftcache, where\n\
                        Start and End are the start and end of the\n\
                        segments in the seg-file.\n"\
" --sft-lalcache (-C)   file containing a list of SFTs in the LALCache\n\
                        format.\n"\
"\nThe following are not required but are set to defaults if not\n\
specified:\n\n"\
" --starttime (-s)      start time of the analysis (default 0 i.e. from\n\
                        the start of the available data/segment list)\n"\
" --finishtime (-f)     finish time of the analysis (default Infinity\n\
                        i.e. to the end of the available data/segment\n\
                        list)\n"\
" --psr-name (-N)       set the name of the pulsar for the output file.\n\
                        This option will overwrite any name set in\n\
                        the parfile. The name is only set if the\n\
                        param-file option is used. (default: not set)\n"\
" --stddev-thresh (-T)  Set the number of standard deviations to use as\n\
                        the threshold for removing outliers. Set to\n\
                        zero for no outlier removal to be performed.\n\
                        (default zero)\n"\
" --freq-factor (-m)    factor which multiplies the pulsar spin\n\
                        frequency (default 2.0 i.e. a signal from a\n\
                        triaxial pulsar)\n"\
" --bandwidth (-b)      width of frequency band around central\n\
                        frequency to use when performing the\n\
                        interpolation. (default 0.3Hz)\n"\
" --min-seg-length (-M) Minimum length of segments (default 1800s)\n"\
" --max-seg-length (-Z) Maximum length of segments (default INFINITY).\n\
                        Use this to reduce memory requirement if\n\
                        trying to read in SFTs for long segement.\n"\
"\nThe following flags are used in testing only, but are included here\n\
for completeness (defaults for all are not to set the flag) \n"\
" --geocentreFlag (-g)  Flag to set the position of the ifo to be at\n\
                        the geocentre.\n"\
" --baryFlag (-B)       Flag to set the position of the ifo to be at\n\
                        the solar system barycentre.\n"\
" --output-timing (-t)  Flags whether to print timing information to\n\
                        stderr\n"\
"\n"

#define XLAL_FRESNEL_EPS 6.0e-8
#define XLAL_FRESNEL_MAXIT 100
#define XLAL_FRESNEL_FPMIN 1.0e-30
#define XLAL_FRESNEL_XMIN 1.5

#define FILENAME_MAXLEN 1024

typedef struct tagInputParams{
  CHAR ifo[3];

  CHAR pulsarDir[FILENAME_MAXLEN];
  CHAR paramfile[FILENAME_MAXLEN];

  CHAR PSRname [12];
  UINT4 nameset;

  UINT4 parfileflag;
  UINT4 pardirflag;

  REAL8 stddevthresh;
  REAL8 minSegLength;
  REAL8 maxSegLength;

  CHAR ephemdir[FILENAME_MAXLEN];

  CHAR filePattern[FILENAME_MAXLEN];

  CHAR outputdir[FILENAME_MAXLEN];
  CHAR segfile[FILENAME_MAXLEN];

  REAL8 freqfactor;
  REAL8 bandwidth;

  UINT4 geocentre;
  UINT4 baryFlag;
  UINT4 Timing;
  UINT4 cacheDir;
  UINT4 cacheFile;
  UINT4 lalcacheFile;

  REAL8 startF;
  REAL8 endF;

  REAL8 startTime;
  REAL8 endTime;

}InputParams;

typedef struct tagSplInterParams{
  LALDetector detector;

  CHAR timefile[FILENAME_MAXLEN];
  CHAR earthfile[FILENAME_MAXLEN];
  CHAR sunfile[FILENAME_MAXLEN];
}SplInterParams;

void get_input_args(InputParams *inputParam, int argc, char *argv[]);

INT4 remove_outliers_using_running_median_data(REAL8Vector *redata, REAL8Vector *imdata,  REAL8Vector *rermdata,
  REAL8Vector *imrmdata, REAL8Vector *times, REAL8 stddevthresh);

REAL8Vector *subtract_running_median( REAL8Vector *data, REAL8Vector *timeStamp , UINT4 npoints);

INT4 XLALFresnel(REAL8 *C, REAL8 *S, REAL8 x);

#endif /* _SPLINTER_PULSAR_H */
