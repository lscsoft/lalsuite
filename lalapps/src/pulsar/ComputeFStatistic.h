#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
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
  int    Dterms ;
  int IFO;
  int noise;
  double startingdemodfreq;
  double demodband;
  double demodfreqres;
  REAL8 skyalpha;
  REAL8 skydelta;
  char *directory;
  char *efiles;
  int spinDwnOrder;
  REAL8 spinParams[4];    /*  For now maximum of 4 spindown params  */
} CommandLineArgs;

typedef struct GlobalVariablesTag {
  double startingdemodfreq;
  double demodband;
  double demodfreqres;
  int imax;
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
int CreateDemodParams(struct CommandLineArgsTag CLA);
int NormaliseSFTData();
int CreateDetector(LALDetector *Detector);
int Freemem();

















