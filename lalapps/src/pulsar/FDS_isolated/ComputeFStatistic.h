#ifndef _COMPUTEFSTATISTIC_H
#define _COMPUTEFSTATISTIC_H


#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>

#include "DopplerScan.h"

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#define BUFFERSIZE 1024                                                                   

/* Maximum fractional doppler shift */
#define DOPPLERMAX 1.046e-4

/* (Half the ) number of terms to keep in the Dirichlet kernel sum */
#define NTERMS 32

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
  INT4 ifmax;			/**< highest frequency-bin needed in calculation */
  INT4 ifmin;			/**< lowest frequency-bin needed */
  REAL8 dFreq;			/**< search frequency spacing  */
  REAL8 df1dot;			/**< spindown resolution (f1 = df/dt!!) */
  REAL8 tsft;			/**< length of an SFT in seconds */
  INT4 SFTno;			/**< number of SFTs in input */
  INT4 nsamples;		/**< number of frequency-bins in an SFT */
  INT4 Ti;	        	/**< GPS seconds of start of observation */
  INT4 Tf;      	  	/**< GPS-time (seconds) of end of observation */
  CHAR **filelist; 		/**< array of filenames to load SFTs from */
  LALDetector Detector;         /**< detector of data to be searched */
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  DopplerRegion searchRegion;	/**< parameter-space region to search over */
  INT4 FreqImax;		/**< number of frequency-bins to run LALDemod for */
  INT4 SpinImax;		/**< number of spindown-values */
} ConfigVariables;
  
struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
  

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
