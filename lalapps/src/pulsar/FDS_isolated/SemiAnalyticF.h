#include <unistd.h>
#include <sys/types.h>
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

#define BUFFERSIZE 1024                                                                   

/* Maximum fractional doppler shift */
#define DOPPLERMAX 1.e-4

/* (Half the ) number of terms to keep in the Dirichlet kernel sum */
#define NTERMS 32

#define MAXFILES 40000         /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/* The command line arguments are not the actual search parameters that will be used, 
just the ones that the user requested.  The actual parameters that
will be used are computed as exact integers, which correspond as closely as
possible to the users request, but are not the same.  Don't use these variables! */

struct CommandLineArgsTag {
  REAL8 skyalpha;
  REAL8 skydelta;
  REAL8 tsft;
  int nTsft;
  int detector;
  char *timestamps;
  int gpsStart;
  char *efiles;
  REAL8 phi;
  REAL8 psi;
  REAL8 h0;
  REAL8 cosiota;
  REAL8 sqrtSh;
} CommandLineArgs;

/* Function Prototypes */

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadTimeStamps(struct CommandLineArgsTag CLA);
int MakeTimeStamps(struct CommandLineArgsTag CLA);
int ComputeF(struct CommandLineArgsTag CLA);
int CreateDetector(LALDetector *Detector);
int Freemem(void);

















