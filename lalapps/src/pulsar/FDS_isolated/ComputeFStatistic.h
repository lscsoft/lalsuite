#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <getopt.h>
#include <time.h>
#include <errno.h>
#include <lal/AVFactories.h>
#include <lal/ComputeSky.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDemod.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define BUFFERSIZE 1024                                                                   

/* Maximum fractional doppler shift */
#define DOPPLERMAX 1.e-4

/* (Half the ) number of terms to keep in the Dirichlet kernel sum */
#define NTERMS 32

#define MAXFILES 40000         /* Maximum # of files in a directory  */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/* The command line arguments are not the actual search parameters that will be used, 
just the ones that the user requested.  The actual parameters that
will be used are computed as exact integers, which correspond as closely as
possible to the users request, but are not the same.  Don't use these variables! */

struct CommandLineArgsTag {
  INT4 Dterms ;
  INT4 IFO;
  INT4 noise;
  BOOLEAN EstimSigParam;
  REAL8 Freq;
  REAL8 dFreq;
  REAL8 FreqBand;
  REAL8 Alpha;
  REAL8 dAlpha;
  REAL8 AlphaBand;
  REAL8 Delta;
  REAL8 dDelta;
  REAL8 DeltaBand;
  REAL8 Spin;
  REAL8 dSpin;
  REAL8 SpinBand;
  REAL8 Fthreshold;
  CHAR *DataDir;
  CHAR *EphemDir;
  CHAR *EphemYear;
  CHAR *BaseName;
} CommandLineArgs;

typedef struct GlobalVariablesTag {
  char EphemEarth[256];
  char EphemSun[256];
  REAL8 Freq;
  REAL8 dFreq;
  REAL8 FreqBand;
  INT4 FreqImax;  // number of computed F values: F[0]....F[FreqImax-1]
  REAL8 Alpha;
  REAL8 dAlpha;
  REAL8 AlphaBand;
  //  INT4 AlphaImax;  // FIXME: these two are now in ScanSky.c
  //  INT4 DeltaImax;
  REAL8 Delta;
  REAL8 dDelta;
  REAL8 DeltaBand;

  REAL8 Spin;
  REAL8 dSpin;
  REAL8 SpinBand;
  INT4 SpinImax;
  INT4 ifmax;
  INT4 ifmin;
  INT4 Dterms;
  REAL8 df;
  REAL8 tsft;
  INT4 SFTno;
  INT4 noise;
  INT4 nsamples;
  INT4 Ti;        /* GPS seconds of first SFT */
  INT4 Tf;        /* GPS seconds of last SFT */
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];
  INT4 IFO;
  BOOLEAN EstimSigParam;
  REAL8 Fthreshold;

  INT2 useMetric;	/* use metric grid or "manual" stepping : 0 = manual, 1 = PtoleMetric, 2 = CoherentMetric */
  REAL8 metricMismatch;	/* maximum allowed mismatch for metric grid */
  BOOLEAN flipTiling;	/* use non-standard internal grid order? ORDER_DELTA_ALPHA */
  LALDetector Detector;              /* Our detector*/
} GlobalVariables;
  
struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
  
/* Function Prototypes */

INT4 ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
INT4 SetGlobalVariables(struct CommandLineArgsTag CLA);
INT4 ReadSFTData(void);
INT4 CreateDemodParams(void);
INT4 NormaliseSFTData(void);
INT4 CreateDetector(LALDetector *Detector);
INT4 AllocateMem(void);
INT4 Freemem(void);

INT4 EstimateSignalParameters(INT4 * maxIndex);
INT4 EstimatePSDLines(void);
INT4 EstimateFLines(void);
INT4 writeFLines(INT4 *maxIndex);
INT4 NormaliseSFTDataRngMdn(void);
INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN);
INT4 EstimateFloor(REAL8Vector *Sp, INT2 windowSize, REAL8Vector *SpFloor);
int compare(const void *ip, const void *jp);
INT4 MedianBias(INT2 * BlockSize, REAL8 * medianbias);
INT4 writeFaFb(INT4 *maxIndex);










