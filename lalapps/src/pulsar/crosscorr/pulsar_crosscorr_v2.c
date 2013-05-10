/*
 *  Copyright (C) 2013 Badri Krishnan, Shane Larson, John Whelan
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \author B.Krishnan, S.Larson, J.T.Whelan
 * \date 2013
 * \file
 * \ingroup pulsarApps
 * \brief Perform CW cross-correlation search - version 2
 *
 * 
 *


/ lalapps includes */
#include <lalapps.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>
#include <lal/DopplerScan.h>
#include <lal/ExtrapolatePulsarSpins.h>


extern int lalDebugLevel;

/* user input variables */
typedef struct{
  BOOLEAN help; /**< if the user wants a help message */
  INT4    startTime;   /**< desired start GPS time of search */ 
  INT4    endTime;     /**< desired end GPS time */
  REAL8   fStart;      /**< start frequency */
  REAL8   fBand;       /**< frequency band to search over */
  REAL8   fdotStart;   /**< starting value for first spindown */
  REAL8   fdotBand;    /**< range of first spindown to search over */
  REAL8   refTime;     /**< reference time for pulsar phase definition */
  CHAR    *sftLocation;/**< location of SFT data */
  CHAR    *ephemYear;  /**< range of years for ephemeris file */
} UserInput_t;

/* struct to store useful stuff */
typedef struct{
  EphemerisData *edat; /**< ephemeris data */
} ConfigVariables;


#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
#define EPHEM_YEARS "00-19-DE405"

#define TRUE (1==1)
#define FALSE (1==0)
#define MAXFILENAMELENGTH 512

/* empty user input struct for initialization */
UserInput_t empty_UserInput;

/* local function prototypes */
int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitializeConfigVars (ConfigVariables *config, const UserInput_t *uvar);

int main(int argc, char *argv[]){
  /* LALStatus pointer */
  static LALStatus status;  

  UserInput_t uvar = empty_UserInput;

  /* sft related variables */ 
  MultiSFTVector *inputSFTs = NULL;
  LIGOTimeGPS firstTimeStamp, lastTimeStamp;
  REAL8 tObs;
  SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;
  REAL8 fMin, fMax; /* min and max freuencies read from SFTs */
  



  /* read lal-debug leve with short-option -v */
   if ( XLALGetDebugLevel ( argc, argv, 'v') != XLAL_SUCCESS )
     XLAL_ERROR ( XLAL_EFUNC );

   /* initialize and register user variables */
   if ( XLALInitUserVars( &uvar ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALInitUserVars() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }

  /* read user input from the command line or config file */
  if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALUserVarReadAllInput() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }

  if (uvar.help)	/* if help was requested, then exit */
    return 0;

  
  /* set sft catalog constraints */
  constraints.detector = NULL;
  constraints.timestamps = NULL;
  XLALGPSSet( constraints.startTime, uvar.startTime, 0);
  XLALGPSSet( constraints.endTime, uvar.endTime,0); 
  if ( (constraints.startTime == NULL)&& (constraints.startTime == NULL) ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALGPSSet() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }

  

  /* FIXME: need to correct fMin and fMax for Doppler shift, rndmedian bins and spindown range */
  /* this is essentially just a place holder for now */
  fMin = uvar.fStart;
  fMin = uvar.fStart + uvar.fBand;

  /* get catalog of SFTs and load them */
  if ((catalog = XLALSFTdataFind (uvar.sftLocation, &constraints)) == NULL){ 
    LogPrintf ( LOG_CRITICAL, "%s: XLALSFTdataFind() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }

  if ((inputSFTs = XLALLoadMultiSFTs ( catalog, fMin, fMax)) == NULL){ 
    LogPrintf ( LOG_CRITICAL, "%s: XLALLoadSFTs() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }


  


  /* /\* get SFT parameters so that we can initialise search frequency resolutions *\/ */
  /* /\* calculate deltaF_SFT *\/ */
  /* deltaF_SFT = catalog->data[0].header.deltaF;  /\* frequency resolution *\/ */
  /* timeBase= 1.0/deltaF_SFT; /\* sft baseline *\/ */

  /* /\* catalog is ordered in time so we can get start, end time and tObs *\/ */
  /* firstTimeStamp = catalog->data[0].header.epoch; */
  /* lastTimeStamp = catalog->data[catalog->length - 1].header.epoch; */
  /* tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase; */

  /* /\*set pulsar reference time *\/ */
  /* if (LALUserVarWasSet ( &uvar_refTime )) { */
  /*   XLALGPSSetREAL8(&refTime, uvar_refTime); */
  /* }  */
  /* else {	/\*if refTime is not set, set it to midpoint of sfts*\/ */
  /*   XLALGPSSetREAL8(&refTime, (0.5*tObs) + XLALGPSGetREAL8(&firstTimeStamp));  */
  /* } */

  /* /\* set frequency resolution defaults if not set by user *\/ */
  /* if (!(LALUserVarWasSet (&uvar_fResolution))) { */
  /*   uvar_fResolution = 1/tObs; */
  /* } */

  /* /\*  set up ephemeris  *\/ */
  /* if(uvar_ephemDir) { */
  /*   snprintf(EphemEarth, MAXFILENAMELENGTH, "%s/earth%s.dat", */
  /* 		uvar_ephemDir, uvar_ephemYear); */
  /*   snprintf(EphemSun, MAXFILENAMELENGTH, "%s/sun%s.dat", */
  /* 		uvar_ephemDir, uvar_ephemYear); */
  /* } else { */
  /*   snprintf(EphemEarth, MAXFILENAMELENGTH, "earth%s.dat", uvar_ephemYear); */
  /*   snprintf(EphemSun, MAXFILENAMELENGTH, "sun%s.dat", uvar_ephemYear); */
  /* } */

  /* EphemEarth[MAXFILENAMELENGTH-1] = 0; */
  /* EphemSun[MAXFILENAMELENGTH-1] = 0; */

  /* edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData)); */
  /* (*edat).ephiles.earthEphemeris = EphemEarth; */
  /* (*edat).ephiles.sunEphemeris = EphemSun; */

  /* LAL_CALL( LALInitBarycenter( &status, edat), &status); */


  /* { */
  /*   /\* block for calculating frequency range to read from SFTs *\/ */
  /*   /\* user specifies freq and fdot range at reftime */
  /*      we translate this range of fdots to start and endtime and find */
  /*      the largest frequency band required to cover the  */
  /*      frequency evolution  *\/ */
  /*   PulsarSpinRange spinRange_startTime; /\**< freq and fdot range at start-time of observation *\/ */
  /*   PulsarSpinRange spinRange_endTime;   /\**< freq and fdot range at end-time of observation *\/ */
  /*   PulsarSpinRange spinRange_refTime;   /\**< freq and fdot range at the reference time *\/ */

  /*   REAL8 startTime_freqLo, startTime_freqHi, endTime_freqLo, endTime_freqHi, freqLo, freqHi; */

  /*   REAL8Vector *fdotsMin=NULL; */
  /*   REAL8Vector *fdotsMax=NULL; */

  /*   UINT4 k; */

  /*   fdotsMin = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector)); */
  /*   fdotsMin->length = N_SPINDOWN_DERIVS; */
  /*   fdotsMin->data = (REAL8 *)LALCalloc(fdotsMin->length, sizeof(REAL8)); */

  /*   fdotsMax = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector)); */
  /*   fdotsMax->length = N_SPINDOWN_DERIVS; */
  /*   fdotsMax->data = (REAL8 *)LALCalloc(fdotsMax->length, sizeof(REAL8)); */

  /*   INIT_MEM(spinRange_startTime); */
  /*   INIT_MEM(spinRange_endTime); */
  /*   INIT_MEM(spinRange_refTime); */

  /*   spinRange_refTime.refTime = refTime; */
  /*   spinRange_refTime.fkdot[0] = uvar_f0; */
  /*   spinRange_refTime.fkdotBand[0] = uvar_fBand; */
  /* } */

  /* LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status); */


  XLALDestroySFTCatalog (catalog);


  /* de-allocate memory for user input variables */
  XLALDestroyUserVars();

  /* check memory leaks if we forgot to de-allocate anything */
  LALCheckMemoryLeaks();

  return 0;

} /* main */


/* initialize and register user variables */
int XLALInitUserVars (UserInput_t *uvar)
{

  /* initialize with some defaults */
  uvar->help = FALSE;
  uvar->startTime = 814838413;	/* 1 Nov 2005, ~ start of S5 */
  uvar->endTime = uvar->startTime + (INT4) round ( LAL_YRSID_SI ) ;	/* 1 year of data */
  uvar->fStart = 100.0; 
  uvar->fBand = 0.1;
  uvar->fdotStart = 0.0;
  uvar->fdotBand = 0.0;

  /* default for reftime is in the middle */
  uvar->refTime = 0.5*(uvar->startTime + uvar->endTime);

  uvar->ephemYear = XLALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

  uvar->sftLocation = XLALCalloc(1, MAXFILENAMELENGTH+1);

  /* register  user-variables */
  XLALregBOOLUserStruct ( help, 	   'h',  UVAR_HELP, "Print this message");  
  
  XLALregINTUserStruct   ( startTime,       0,  UVAR_OPTIONAL, "Desired start time of analysis in GPS seconds");
  XLALregINTUserStruct   ( endTime,         0,  UVAR_OPTIONAL, "Desired end time of analysis in GPS seconds");
  XLALregREALUserStruct  ( fStart,          0,  UVAR_OPTIONAL, "Start frequency in Hz");
  XLALregREALUserStruct  ( fBand,           0,  UVAR_OPTIONAL, "Frequency band to search over in HZ ");
  XLALregREALUserStruct  ( fdotStart,       0,  UVAR_OPTIONAL, "Start value of spindown in Hz/s");
  XLALregREALUserStruct  ( fdotBand,        0,  UVAR_OPTIONAL, "Band for spindown values in Hz/s");
  XLALregSTRINGUserStruct( ephemYear,       0,  UVAR_OPTIONAL, "String Ephemeris year range");
  XLALregSTRINGUserStruct( sftLocation,     0,  UVAR_REQUIRED, "Filename pattern for locating SFT data");

  if ( xlalErrno ) {
    XLALPrintError ("%s: user variable initialization failed with errno = %d.\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  return XLAL_SUCCESS;
}



/* initialize and register user variables */
int XLALInitializeConfigVars (ConfigVariables *config, const UserInput_t *uvar)
{
  
  return XLAL_SUCCESS;

} /* XLALInitializeConfigVars() */
