#include <unistd.h>
#include <sys/types.h>
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
#include <lal/ComputeSkyBinary.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDemod.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define BUFFERSIZE 1024                                                                   

/* Maximum fractional doppler shift */
#define DOPPLERMAX 4.0e-4  /* this value is double the maximum expected value for sco X1 */

/* (Half the ) number of terms to keep in the Dirichlet kernel sum */
#define NTERMS 32

#define MAXFILES 100000         /* Maximum # of files in a directory  */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/* The command line arguments are not the actual search parameters that will be used, 
just the ones that the user requested.  The actual parameters that
will be used are computed as exact integers, which correspond as closely as
possible to the users request, but are not the same.  Don't use these variables! */

struct CommandLineArgsTag {
  int Dterms ;
  int IFO;
  int noise;
  int EstimSigParam;
  double Freq;
  double dFreq;
  double FreqBand;
  double skyAlpha;
  double skyDelta;
  int spinDwnOrder;
  double Fthreshold;
  const char *DataDir;
  const char *EphemDir;
  const char *EphemYear;
  const char *filters;
} CommandLineArgs;

typedef struct GlobalVariablesTag {
  char EphemEarth[256];
  char EphemSun[256];
  double Freq;
  double dFreq;
  double FreqBand;
  int FreqImax;
  double skyAlpha;
  double skyDelta;
  int spinDwnOrder;
  double SpinDownOne;
  int ifmax;
  int ifmin;
  int Dterms;
  double df;
  double tsft;
  int SFTno;
  int noise;
  int nsamples;
  int Ti;        /* GPS seconds of first SFT */
  int Tf;        /* GPS seconds of last SFT */
  char filelist[MAXFILES][MAXFILENAMELENGTH];
  const char *filters;
  double SemiMajorAxis;
  double OrbitalPeriod;
  LIGOTimeGPS TperiapseSSB;
  double OrbitalEccentricity;
  double ArgPeriapse;
  int IFO;
  int EstimSigParam;
  double Fthreshold;
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

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int SetGlobalVariables(struct CommandLineArgsTag CLA);
int ReadSFTData(void);
int CreateDemodParamsBinary(struct CommandLineArgsTag CLA);
int NormaliseSFTData(void);
int CreateDetector(LALDetector *Detector);
int AllocateMem(void);
int Freemem(void);
int EstimateSignalParameters(INT4 * maxIndex);
int EstimatePSDLines(void);
int EstimateFLines(void);
int writeFLines(INT4 *maxIndex);
int NormaliseSFTDataRngMdn(LALStatus *status);
int PrintTopValues(double TwoFthr, int ReturnMaxN);
int EstimateFloor(REAL8Vector *Sp, INT2 windowSize, REAL8Vector *SpFloor);
int compare(const void *ip, const void *jp);
int MedianBias(INT2 * BlockSize, double * medianbias);
int writeFaFb(INT4 *maxIndex);










