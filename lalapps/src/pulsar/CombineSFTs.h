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
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <errno.h>

#define MAXFILES 60000         // Maximum # of files in a directory 
#define MAXFILENAMELENGTH 256   // Maximum # of characters of a SFT filename

/* The command line arguments are not the actual search parameters that will be used, 
just the ones that the user requested.  The actual parameters that
will be used are computed as exact integers, which correspond as closely as
possible to the users request, but are not the same.  Don't use these variables! */

struct CommandLineArgsTag {
   int Dterms;
  char *inputdirectory;
  char *outputdirectory;
  int number;
} CommandLineArgs;
  
struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
  
typedef struct FFTTag 
{
  COMPLEX8FrequencySeries *fft;   
} FFT;

/* Function Prototypes */

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int CreateFileList(struct CommandLineArgsTag CLA);
int ReadSFTs(struct CommandLineArgsTag CLA, int i);
int CSFTs(struct CommandLineArgsTag CLA);
int AllocateMem(struct CommandLineArgsTag CLA);














