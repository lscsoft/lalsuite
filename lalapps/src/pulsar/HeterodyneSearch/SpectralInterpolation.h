#define _GNU_SOURCE   /* for alphasort() and scandir() */

#include <sys/time.h>
#include <lalapps.h>
#include <lal/LALConstants.h>

#include <lal/BinaryPulsarTiming.h>
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
#include <dirent.h>
#include <lal/BinaryPulsarTiming.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

/* LAL headers */
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/SFTfileIO.h>
/* lalapps header */
#include <lalapps.h>

/* frame headers */
#include <FrIO.h>
#include <lal/LALFrameL.h>

/* Usage format string. */
 #define USAGE \
"Usage: lalapps_SplInter [options]\n\n"\
" --help (-h)              display this message\n"\
" \nRequired Input Variables:\n\n"\
" --ifo (-i)               name of ifo e.g. L1, H1, H2, G1\n"\
" --start-freq (-S)        Start frequency of the SFTs\n"\
" --end-freq (-E)          End frequency of the SFTs\n"\
" --ephem-dir (-e)         directory containing the ephemeris files\n"\
" --output-dir (-o)        directory for output data files\n"\
" --seg-file (-l)          name of file containing science segment list\n"\
"\n Plus one of the following:\n\n"\
" --param-file (-P)        name of file containing initial pulsar parameters\n\
                             when using just one file as input (.par file)\n"\
" --param-dir (-d)         name of directory containing pulsar parameter files\n\
                             when using multiple input files (.par files)\n"\
"\n and one of the following:\n\n"\
" --sft-cache (-F)         location of list of SFTs for use, can use \n\
                             list:listsfts.txt or /location/*.sft\n"\
" --sft-loc (-L)           directory of files containing a list of SFTs for \n\
                             each segment, located in \n\
                             /sft-loc/Segment_Start-End.sftcache, where Start and \n\
                             End are the start and end of the segments in the seg-file. \n"\
" \nThe following are not required but are set to defaults if not specified:\n\n"\
" --starttime (-s)         start time of the analysis \n\
                             (default 0 - i.e. whenever parfile(s) are applicable)\n"\
" --finishtime (-f)        finish time of the analysis \n\
                             (default Infinity - i.e. whenever parfile(s) are applicable)\n"\
" --psr-name (-N)          set the name of the pulsar for the output file.\n\
                             This option will overwrite any name set in the parfile.\n\
                             The name is only set if the param-file option is used.\n\
                             (default: not set) \n"\
" --stddev-thresh (-T)     Set the number of standard deviations to use as the\n\
                             threshold for removing outliers. Set to zero for no\n\
                             outlier removal to be performed. (default zero)\n"\
" --freq-factor (-m)       factor which multiplies the pulsar spin frequency\n\
                             (default 2.0 i.e. a signal from a triaxial pulsar)\n"\
" --bandwidth (-b)         width of frequency band around central frequency to\n\
                             use when performing the interpolation.\n\
                             (default 0.3Hz)\n"\
" --min-seg-length (-M)    Minimum length of segments (default 5400s)\n"\
" \nThe following are 'flags', and require no associated value:\n\n"\
" --noise-output (-n)      flags whether to print noise output into the B_k output file\n\
                             (default not to output)\n"\
" --output-timing (-t)     flags whether to print timing information to stderr\n\
                             (default not to output)\n"\
" \nThe following are used in testing only, but are included here for completeness:\n\n"\
" --geocentreFlag (-g)     Flag to set the position of the ifo to be \n\
                             at the geocentre.\n"\
" --baryFlag (-B)          Flag to set the position of the ifo to be \n\
                             at the solar system barycentre.\n"\
"\n"

#define XLAL_FRESNEL_EPS 6.0e-8
#define XLAL_FRESNEL_MAXIT 100
#define XLAL_FRESNEL_FPMIN 1.0e-30
#define XLAL_FRESNEL_XMIN 1.5


typedef struct tagInputParams{
  CHAR ifo[3];

  CHAR pulsarDir[256];
  CHAR paramfile[256];

  CHAR PSRname [12];
  UINT4 nameset;

  UINT4 parfileflag;
  UINT4 pardirflag;

  REAL8 stddevthresh;
  REAL8 minSegLength;

  CHAR ephemdir[256];

  CHAR filePattern[256];

  CHAR outputdir[256];
  CHAR segfile[256];

  REAL8 freqfactor;
  REAL8 bandwidth;

  UINT4 geocentre;
  UINT4 baryFlag;
  UINT4	Timing;
  UINT4	cacheDir;
  UINT4	cacheFile;

  REAL8 startF;
  REAL8 endF;

  REAL8 startTime;
  REAL8 endTime;

}InputParams;

typedef struct tagSplInterParams{
  LALDetector detector;

  CHAR timefile[256];
  CHAR earthfile[256];
  CHAR sunfile[256];
}SplInterParams;

void get_input_args(InputParams *inputParam, int argc, char *argv[]);

INT4 get_segment_list(INT4Vector *start, INT4Vector *stop, REAL8 minDur, REAL8 startSegs, REAL8 endSegs, CHAR *seglistfile);

INT4 remove_outliers_using_running_median_data(REAL8Vector *redata, REAL8Vector *imdata,  REAL8Vector *rermdata,
  REAL8Vector *imrmdata, REAL8Vector *times, REAL8 stddevthresh);

REAL8Vector *subtract_running_median( REAL8Vector *data, REAL8Vector *timeStamp , UINT4 npoints);

INT4 XLALFresnel(REAL8 *C, REAL8 *S, REAL8 x);
