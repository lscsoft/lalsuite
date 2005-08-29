#ifndef _COMPUTEFSTATISTIC_H
#define _COMPUTEFSTATISTIC_H

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDemod.h>

#include "DopplerScan.h"

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#define BUFFERSIZE 1024                                                                   

/*----------------------------------------------------------------------*/
/* Error-codes */

#define COMPUTEFSTATC_ENULL             1
#define COMPUTEFSTATC_ESYS              2
#define COMPUTEFSTATC_EINPUT            3
#define COMPUTEFSTATC_EMEM              4
#define COMPUTEFSTATC_ECHECKPOINT       5
#define COMPUTEFSTATC_ECLUSTER          6
#define COMPUTEFSTATC_EXLAL		50

#define COMPUTEFSTATC_MSGENULL          "Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGESYS           "System call failed (probably file IO)"
#define COMPUTEFSTATC_MSGEINPUT         "Invalid input"
#define COMPUTEFSTATC_MSGEMEM           "Out of memory. Bad."
#define COMPUTEFSTATC_MSGECHECKPOINT    "Illegal checkpoint"
#define COMPUTEFSTATC_MSGECLUSTER       "Unspecified error in cluster-related routine"
#define COMPUTEFSTATC_MSGEXLAL		"XLALFunction-call failed"

/*----------------------------------------------------------------------*/
/* Exit values */
#define COMPUTEFSTAT_EXIT_OK              0  /* normal exit */
#define COMPUTEFSTAT_EXIT_USAGE           7  /* user requested help */
#define COMPUTEFSTAT_EXIT_READSFTFAIL     8  /* ReadSFT failed */
#define COMPUTEFSTAT_EXIT_OPENFMAX        9  /* error opening Fmax file */
#define COMPUTEFSTAT_EXIT_OPENFSTAT      10  /* error opening FStats file */
#define COMPUTEFSTAT_EXIT_OPENFSTAT2     11  /* error opening FStats file for append (chkpt) */
#define COMPUTEFSTAT_EXIT_WRITEFSTAT     12  /* error writing FStats file */
#define COMPUTEFSTAT_EXIT_WRITEFAFB      13  /* writeFaFb failed */
#define COMPUTEFSTAT_EXIT_ESTSIGPAR      14  /* EstimateSignalParameters failed */
#define COMPUTEFSTAT_EXIT_NOMEM          15  /* out of memory */
#define COMPUTEFSTAT_EXIT_CANTZIP        16  /* unable to zip Fstats file */
#define COMPUTEFSTAT_EXIT_CANTUNZIP      17  /* unable to zip Fstats file */
#define COMPUTEFSTAT_EXIT_CANTRENAME     18  /* unable to zip Fstats file */
#define COMPUTEFSTAT_EXIT_NOPOLKADEL     19  /* no // found in command line */
#define COMPUTEFSTAT_EXIT_USER           20  /* user asked for exit */
#define COMPUTEFSTAT_EXIT_DEMOD          21  /* error in LAL-Demod */
#define COMPUTEFSTAT_EXIT_SIGNAL         22  /* exited on signal */
#define COMPUTEFSTAT_EXIT_BOINCRESOLVE   23  /* boinc_resolve_filename failed */
#define COMPUTEFSTAT_EXIT_DLOPEN         24  /* problems with dynamic lib */
#define COMPUTEFSTAT_EXIT_WORKER         25  /* can't start worker-thread */

#define COMPUTEFSTAT_EXIT_LALCALLERROR  100  /* added to LAL status for BOINC exit value */

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
  REAL8 DeltaFreqRef;		/**< difference f(tRef) - f(tStart) */
} ConfigVariables;


/** Type to hold the fields that will be output in unclustered output file  */
typedef struct {
  REAL8 Freq;			/**< Frequency at maximum (?) of the cluster */
  REAL8 f1dot;			/**< spindown value f1dot = df/dt */
  REAL8 Alpha; 			/**< Skyposition: longitude in equatorial coords, radians */
  REAL8 Delta;			/**< skyposition: latitude */
  REAL8 Fstat;			/**< value of 2F */
} FstatOutputEntry;

  
/* LALDemod functions now put into CFSLALDemod.c */
extern void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params);

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
