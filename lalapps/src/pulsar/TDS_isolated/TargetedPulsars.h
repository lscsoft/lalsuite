/*
$Id$
*/
// standard C headers
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

// LAL headers
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>

// new LAL headers
#include "HeterodynePulsar.h" /* use this local version until up-to-date version commited to lal*/

// Frame headers
#include <FrIO.h>
#include <FrameL.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

INT4 lalDebugLevel = LALWARNING | LALINFO;

#define SRATE 16384  // sample rate
#define MAXSTART 5000 // maximum number of start times from input file
#define MAXLINES 150000 //maximum number of lines in alpha/beta file (num of min)
#define MAXFINE 2000 // maximum number of fine heterodyned data from each node

#define TBASE 60

// filter info
#define fc1 0.50 /* cutoff frequency of first IIR */
#define fc2 0.50  /* cutoff frequency of second IIR */
#define fc3 0.50  /* cutoff frequency of third IIR */

#define SQRT3_2 0.8660254037844386467637231707529361L


/* this should be changed before S3 final results if necessary*/
#define ALPHAMAX 2.3  /* alpha calibration cut */
#define ALPHAMIN 0.25


#define S2_GPS_SECS 729273613.0
#define S2_MJD 52684.66666666665

#define S3_GPS_SECS 751680013.0
#define S3_MJD 52944.0

#define S4_GPS_SECS 793130413.0
#define S4_MJD 53423.75

#define MAXPSRS 150

typedef struct
tagFrFileInfo
{
  INT4  ind;
  CHAR *url;
  INT4  t0;
  INT4  dt;
}
FrFileInfo;

/* Definition of FrStream */
struct
tagFrStream
{
  FrFileInfo     *filelist;
  UINT4           numfiles;
  UINT4           filenum;
  struct FrFile  *frfile;
  struct FrameH  *frame;
  LIGOTimeGPS     epoch;
  INT4            end;
  INT4            err;
  INT4            gap;
};

FILE* tryopen(char *name, char *mode){
  int count=0;
  FILE *fp;

  while (!(fp=fopen(name, mode))){
    fprintf(stderr,"Unable to open file %s in mode %s\n", name, mode);
    fflush(stderr);
    if (count++<10){
      continue; /* slight change as sleep() might effect Condor checkpointing */
			/* sleep(10); */
		}
    else
      exit(3);
  }
  return fp;
}
