#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <lal/AVFactories.h>
#include <lal/ComputeSky.h> 
#include <lal/ComputeSkyBinary.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDemod.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>
#ifndef NOGLOB
#include <glob.h>
#endif
#define BUFFERSIZE 1024                                                                   

/* Maximum fractional doppler shift */  /* BINARY-MOD - Need to sort this */
#define DOPPLERMAX 1.e-4

/* (Half the ) number of terms to keep in the Dirichlet kernel sum */
#define NTERMS 32

#define MAXFILES 40000         /* Maximum # of files in a directory  */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/* The command line arguments are not the actual search parameters that will be used, 
just the ones that the user requested.  The actual parameters that
will be used are computed as exact integers, which correspond as closely as
possible to the users request, but are not the same.  Don't use these variables! */



/* 
 *this structure holds all configuration-settings for the code, including the 
 * user-input variables, but also derived ones 
 */
typedef struct {
  CHAR EphemEarth[MAXFILENAMELENGTH];
  CHAR EphemSun[MAXFILENAMELENGTH];
  INT4 FreqImax;  /* number of computed F values: F[0]....F[FreqImax-1] */

  INT4 SpinImax;
  INT4 ifmax;
  INT4 ifmin;

  REAL8 dFreq;
  REAL8 tsft;
  INT4 SFTno;

  INT4 nsamples;
  INT4 Ti;        /* GPS seconds of first SFT */
  INT4 Tf;        /* GPS seconds of last SFT */
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];

  LALDetector Detector;         /* Our detector*/
  EphemerisData *edat;		/* ephemeris data */

  CHAR *skyRegion;		/* sky-region (polygon) to search over */
} ConfigVariables;
  
struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
  
typedef struct BinaryTemplatetag {             /* BINARY-MOD - structure to store a single binary signal template */
    REAL8       ProjSMaxis;
    REAL8       Period;
    LIGOTimeGPS TperiSSB;
    REAL8       Eccentricity;
    REAL8       ArgPeri;
} BinaryTemplate;

typedef struct BinaryTemplateBanktag {
    REAL8          BTBfmax;
    REAL8          BTBTspan;
    LIGOTimeGPS    BTBTobsStart;
    INT4           BTBNFilters;
    REAL8          BTBMismatch;
    REAL8          BTBProjSMaxisMIN;
    REAL8          BTBProjSMaxisMAX;
    LIGOTimeGPS    BTBTperiSSBMIN;
    LIGOTimeGPS    BTBTperiSSBMAX;
    REAL8          BTBEccMIN;
    REAL8          BTBEccMAX;
    REAL8          BTBArgPeriMIN;
    REAL8          BTBArgPeriMAX;
    REAL8          BTBPeriodMIN;
    REAL8          BTBPeriodMAX;
    CHAR           BTBversion[256];
    BinaryTemplate *BTB;       
} BinaryTemplateBank;

