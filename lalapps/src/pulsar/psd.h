#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALDemod.h>

#define MAXFILES 60000         /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

struct CommandLineArgsTag 
{
  char *directory;
  char *outputfile;
} CommandLineArgs;

struct headertag 
{
  REAL8 endian;
  INT4  gps_sec;
  INT4  gps_nsec;
  REAL8 tbase;
  INT4  firstfreqindex;
  INT4  nsamples;
} header;

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadSFTDirectory(struct CommandLineArgsTag CLA);
int ComputePSD(struct CommandLineArgsTag CLA);
int Freemem();


