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

#define MAXFILES 70000          /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define MAXLINESRS   4000000     /* Maximum # of lines in a Response or Sensing file */
#define MAXLINESF    100000     /* Maximum # of lines in a Factors file */

struct CommandLineArgsTag 
{
  char *directory;
  char *caldirectory;
  char *run;
  char *IFO; 
  char *outputdirectory;
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

typedef struct ResponseFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Response;

typedef struct SensingFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Sensing;

struct FactorsTag
{
  INT4  t[MAXLINESF];      /* in GPS seconds */
  REAL4 alpha[MAXLINESF];
  REAL4 alpha_beta[MAXLINESF];
} Factors;

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadSFTDirectory(struct CommandLineArgsTag CLA);
int ReadCalibrationFiles(struct CommandLineArgsTag CLA);
int ComputeInitialRSFunctions(struct CommandLineArgsTag CLA);
int CalibrateSfts(struct CommandLineArgsTag CLA);
int Freemem();


COMPLEX8 tmpa, tmpb, tmpc; 
REAL4 tmpx, tmpy;

#define cmul( a, b ) \
( tmpa = (a), tmpb = (b), \
  tmpc.re = tmpa.re * tmpb.re - tmpa.im * tmpb.im, \
  tmpc.im = tmpa.re * tmpb.im + tmpa.im * tmpb.re, \
  tmpc )

#define cdiv( a, b ) \
( tmpa = (a), tmpb = (b), \
  fabs( tmpb.re ) >= fabs( tmpb.im ) ? \
    ( tmpx = tmpb.im / tmpb.re, \
      tmpy = tmpb.re + tmpx * tmpb.im, \
      tmpc.re = ( tmpa.re + tmpx * tmpa.im ) / tmpy, \
      tmpc.im = ( tmpa.im - tmpx * tmpa.re ) / tmpy, \
      tmpc ) : \
    ( tmpx = tmpb.re / tmpb.im, \
      tmpy = tmpb.im + tmpx * tmpb.re, \
      tmpc.re = ( tmpa.re * tmpx + tmpa.im ) / tmpy, \
      tmpc.im = ( tmpa.im * tmpx - tmpa.re ) / tmpy, \
      tmpc ) )

