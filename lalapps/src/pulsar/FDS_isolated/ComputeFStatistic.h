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
#include <lal/LALDatatypes.h>
#include <lal/LALDemod.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>
#ifndef NOGLOB
#include <glob.h>
#endif
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



/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  CHAR EphemEarth[MAXFILENAMELENGTH];	/**< filename of earth-ephemeris data */
  CHAR EphemSun[MAXFILENAMELENGTH];	/**< filename of sun-ephemeris data */
  INT4 FreqImax;  		/**< number of frequency-bins to compute F-stat for */
  INT4 SpinImax;		/**< number of spindown-bins to compute F for */
  INT4 ifmax;			/**< highest frequency-bin needed in calculation */
  INT4 ifmin;			/**< lowest frequency-bin needed */
  REAL8 dFreq;			/**< search frequency spacing  */
  REAL8 tsft;			/**< length of an SFT in seconds */
  INT4 SFTno;			/**< number of SFTs in input */
  INT4 nsamples;		/**< number of frequency-bins in an SFT */
  INT4 Ti;	        	/**< GPS seconds of start of observation */
  INT4 Tf;      	  	/**< GPS-time (seconds) of end of observation */
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH]; /**< array of filenames to load SFTs from */
  LALDetector Detector;         /**< Our detector*/
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  CHAR *skyRegion;		/**< sky-region to search (polygon defined by list of points) */
} ConfigVariables;
  
struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
  
