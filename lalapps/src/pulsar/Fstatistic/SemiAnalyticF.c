/*
*  Copyright (C) 2007 Chris Messenger, Iraj Gholami,  Holger Pletsch, Reinhard Prix, Xavier Siemens
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
 * \file
 * \ingroup lalapps_pulsar_Fstatistic
 * \author Chris Messenger, Iraj Gholami,  Holger Pletsch, Reinhard Prix, Xavier Siemens
 * \brief Semi-Analytic calculation of the F-statistic
 */

/************************************************************************
 ** WARNING: THIS CODE OUTPUTS 'F' INSTEAD OF '2F' WHICH WOULD BE OUR USUAL
 **          CONVENTION FOR THE "F-STATISTIC"
 ************************************************************************
 **
 ** An attempt to change this to 2F to make it consistent with the other
 ** pulgroup codes was rejected by pulgroup. This warning has been put in place
 ** instead: http://blip.phys.uwm.edu/twiki/bin/save/CW/ProposedCodePatches
 **
 ************************************************************************/


/*********************************************************************************/
/*                 Semi-Analytic calculation of the F-statistic                  */
/*                                                                               */
/*			           X. Siemens                                    */
/*                                                                               */
/*                             UWM - February 2003                               */
/*********************************************************************************/

#include <errno.h>

#include <lal/LALString.h>
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/LALStdio.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALComputeAM.h>

#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>

#include <lalapps.h>

/*---------- error-codes ---------- */
#define SEMIANALYTIC_ENORM 		0
#define SEMIANALYTIC_ESUB  		1
#define SEMIANALYTIC_EINPUT  		2
#define SEMIANALYTIC_EBAD  		3
#define SEMIANALYTIC_EFILE 		4
#define SEMIANALYTIC_ENOARG 		5
#define SEMIANALYTIC_EMEM 		6
#define SEMIANALYTIC_EREADFILE 		8

#define SEMIANALYTIC_MSGENORM 		"Normal exit"
#define SEMIANALYTIC_MSGESUB  		"Subroutine failed"
#define SEMIANALYTIC_MSGEINPUT 		"Invalid input"
#define SEMIANALYTIC_MSGEBAD  		"Bad argument values"
#define SEMIANALYTIC_MSGEFILE 		"File IO error"
#define SEMIANALYTIC_MSGENOARG 		"Missing argument"
#define SEMIANALYTIC_MSGEMEM 		"Out of memory..."
#define SEMIANALYTIC_MSGEREADFILE 	"Error reading in file"

/*---------- defines ---------- */
#define TRUE (1==1)
#define FALSE (1==0)

#define SQ(x) ((x)*(x))

/*---------- local types ---------- */
struct CommandLineArgsTag {
  REAL8 Alpha;
  REAL8 Delta;
  REAL8 Tsft;
  INT4  nTsft;
  CHAR  *IFO;
  CHAR  *timestamps;
  INT4  gpsStart;
  REAL8 phi0;
  REAL8 psi;
  REAL8 sqrtSh;
  REAL8 duration;
  CHAR  *ephemEarth;
  CHAR  *ephemSun;
  REAL8 aPlus;
  REAL8 aCross;
  REAL8 h0;
  REAL8 cosi;
  /* ----- deprecated ----- */
  REAL8 cosiota;
  CHAR *detector;
} CommandLineArgs;

/*---------- global variables ---------- */
LIGOTimeGPSVector *timestamps = NULL;
AMCoeffs amc;

extern int vrbflg;

/*---------- local prototypes ---------- */
void InitUserVars (LALStatus *status, struct CommandLineArgsTag *CLA);
void ReadUserInput (LALStatus *, struct CommandLineArgsTag *CLA, int argc,char *argv[]);
void Freemem( LALStatus *);
void Initialize (LALStatus *status, struct CommandLineArgsTag *CLA);
void ComputeF(LALStatus *, struct CommandLineArgsTag CLA);

void CheckUserInput (LALStatus *,  struct CommandLineArgsTag *CLA );

/*---------- function definitions ---------- */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */

  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (InitUserVars(&status, &CommandLineArgs), &status);	  

  /* read cmdline & cfgfile  */	
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN (XLALUserVarReadAllInput(&should_exit, argc, argv) == XLAL_SUCCESS, XLAL_EFUNC);
  if (should_exit)
    return EXIT_FAILURE;
  LAL_CALL ( CheckUserInput (&status, &CommandLineArgs), &status);

  LAL_CALL ( Initialize (&status, &CommandLineArgs), &status);

  /*---------- central function: compute F-statistic ---------- */
  LAL_CALL ( ComputeF(&status, CommandLineArgs), &status); 

  
  /* Free remaining memory */
  LAL_CALL ( LALSDestroyVector(&status, &(amc.a) ), &status);
  LAL_CALL ( LALSDestroyVector(&status, &(amc.b) ), &status);
  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return 0;

} /* main() */


void 
ComputeF( LALStatus *status, struct CommandLineArgsTag CLA)
{

  REAL8 A,B,C,A1,A2,A3,A4, To,Sh,F;
  REAL8 aPlus, aCross;
  REAL8 twopsi, twophi;

  INITSTATUS(status);
  ATTATCHSTATUSPTR ( status);
  
  A = amc.A;
  B = amc.B;
  C = amc.C;

  twophi = 2.0 * CLA.phi0;
  twopsi = 2.0 * CLA.psi;

  aPlus = CLA.aPlus;
  aCross = CLA.aCross;
  
  A1 = aPlus * cos(twopsi) * cos(twophi) - aCross * sin(twopsi) * sin(twophi);
  A2 = aPlus * sin(twopsi) * cos(twophi) + aCross * cos(twopsi) * sin(twophi);
  A3 =-aPlus * cos(twopsi) * sin(twophi) - aCross * sin(twopsi) * cos(twophi);
  A4 =-aPlus * sin(twopsi) * sin(twophi) + aCross * cos(twopsi) * cos(twophi);
  
  To = CLA.nTsft * CLA.Tsft;
  
  Sh=pow(CLA.sqrtSh,2);

  F = A * ( SQ(A1) + SQ(A3) ) + 2.0 * C * (A1 * A2 + A3 * A4 ) + B * ( SQ(A2) + SQ(A4) );
  F *= To / (4.0 * Sh);

  /* Note: the expectation-value of 2F is 4 + lambda ==> add 2 to Fstat*/
  F += 2.0;
  
  fprintf( stdout, "%g\n", F );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeF() */


/**
 * register all our "user-variables"
 */
void
InitUserVars (LALStatus *status, struct CommandLineArgsTag *CLA)
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Initialize default values */
  CLA->Tsft=1800;
  CLA->nTsft=0;            
  CLA->timestamps=NULL;
  CLA->gpsStart=-1;
  CLA->sqrtSh=1.0;
  
  /** Default year-span of ephemeris-files to be used */
  CLA->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  CLA->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");
  
  /* ---------- register all our user-variable ---------- */
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->Alpha),      "Alpha",          REAL8,  'a', OPTIONAL,  "Sky position Alpha (equatorial coordinates) in radians") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->Alpha),      "longitude",      REAL8,  0,   DEVELOPER, "[DEPRECATED] Use --Alpha instead!") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->Delta),      "Delta",          REAL8,  'd', OPTIONAL,  "Sky position Delta (equatorial coordinates) in radians") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->Delta),      "latitude",       REAL8,  0,   DEVELOPER, "[DEPRECATED] Use --Delta instead!") == XLAL_SUCCESS, XLAL_EFUNC);
 
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->phi0),       "phi0",           REAL8,  'Q', OPTIONAL,  "Phi_0: Initial phase in radians") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->psi),        "psi",            REAL8,  'Y', OPTIONAL,  "Polarisation in radians") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->cosi),       "cosi",           REAL8,  'i', OPTIONAL,  "Cos(iota)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->cosiota),    "cosiota",        REAL8,  0,   DEVELOPER, "[DEPRECATED] Use --cosi instead") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->h0),         "h0",             REAL8,  's', OPTIONAL,  "Strain amplitude h_0") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->sqrtSh),     "sqrtSh",         REAL8,  'N', OPTIONAL,  "Noise floor: one-sided sqrt(Sh) in 1/sqrt(Hz)") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->timestamps), "timestampsFile", STRING, 'T', OPTIONAL,  "Name of timestamps file") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->gpsStart),   "startTime",      INT4,   'S', OPTIONAL,  "GPS start time of continuous observation") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->Tsft),       "Tsft",           REAL8,  't', OPTIONAL,  "Length of an SFT in seconds") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->nTsft),      "nTsft",          INT4,   'n', OPTIONAL,  "Number of SFTs") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->IFO),        "IFO",            STRING, 'D', OPTIONAL,  "Detector: H1, H2, L1, G1, ... ") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->detector),   "detector",       STRING, 0,   DEVELOPER, "[DEPRECATED] Use --IFO instead!") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->ephemEarth), "ephemEarth",     STRING, 0,   OPTIONAL,  "Earth ephemeris file to use") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_LAL( status, XLALRegisterNamedUvar(&(CLA->ephemSun),   "ephemSun",       STRING, 0,   OPTIONAL,  "Sun ephemeris file to use") == XLAL_SUCCESS, XLAL_EFUNC);

  /* ----- added for mfd_v4 compatibility ---------- */
  XLAL_CHECK_LAL ( status, XLALRegisterNamedUvar(&(CLA->duration),  "duration",       REAL8, 0,    OPTIONAL,  "Duration of requested signal in seconds") == XLAL_SUCCESS, XLAL_EFUNC);
  
  XLAL_CHECK_LAL ( status, XLALRegisterNamedUvar(&(CLA->aPlus),     "aPlus",          REAL8, 0,    OPTIONAL,  "Plus polarization amplitude aPlus") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL ( status, XLALRegisterNamedUvar(&(CLA->aCross),    "aCross",         REAL8, 0,    OPTIONAL,  "Cross polarization amplitude aCross") == XLAL_SUCCESS, XLAL_EFUNC);


  DETATCHSTATUSPTR (status);
  RETURN(status);
} /* InitUserVars() */
  
/**
 * Handle user-input and check its validity.
 * Load ephemeris and calculate AM-coefficients (stored globally)
 */
void
Initialize (LALStatus *status, struct CommandLineArgsTag *CLA)
{
  EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  EarthState earth;
  AMCoeffsParams *amParams;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  LALDetector *Detector;              /* Our detector*/
  INT4 k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  if ( XLALUserVarWasSet ( &(CLA->nTsft) ) )
    CLA->duration = 1.0 * CLA->nTsft * CLA->Tsft;

  /* read or generate SFT timestamps */
  if ( XLALUserVarWasSet(&(CLA->timestamps)) ) 
    { 
      XLAL_CHECK_LAL ( status, ( timestamps = XLALReadTimestampsFile ( CLA->timestamps ) ) != NULL, XLAL_EFUNC );
      if ( (CLA->nTsft > 0) && ( (UINT4)CLA->nTsft < timestamps->length ) )	/* truncate if required */
	timestamps->length = CLA->nTsft;
      
      CLA->nTsft = timestamps->length;
    } /* if have_timestamps */
  else 
    {
      LIGOTimeGPS tStart;
      tStart.gpsSeconds = CLA->gpsStart;
      tStart.gpsNanoSeconds = 0;

      XLAL_CHECK_LAL ( status, ( timestamps = XLALMakeTimestamps( tStart, CLA->duration, CLA->Tsft, 0 ) ) != NULL, XLAL_EFUNC );
      CLA->nTsft = timestamps->length;

    } /* no timestamps */

  /*---------- initialize detector ---------- */
  {
    BOOLEAN have_IFO       = XLALUserVarWasSet ( &CLA->IFO );
    BOOLEAN have_detector  = XLALUserVarWasSet ( &CLA->detector );
    CHAR *IFO;

    if ( !have_IFO  && !have_detector ) {
      fprintf (stderr, "\nNeed to specify the detector (--IFO) !\n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }
    if ( have_IFO )
      IFO = CLA->IFO;
    else
      IFO = CLA->detector;

    if ( ( Detector = XLALGetSiteInfo ( IFO ) ) == NULL ) {
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }
  }

  /* ---------- load ephemeris-files ---------- */
  {
    edat = XLALInitBarycenter( CLA->ephemEarth, CLA->ephemSun );
    if ( !edat ) {
      XLALPrintError("XLALInitBarycenter failed: could not load Earth ephemeris '%s' and Sun ephemeris '%s'\n", CLA->ephemEarth, CLA->ephemSun);
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }
  } /* ephemeris-reading */


  /* ---------- calculate AM-coefficients ---------- */

  /* prepare call to barycentering routing */
  baryinput.site.location[0] = Detector->location[0]/LAL_C_SI;
  baryinput.site.location[1] = Detector->location[1]/LAL_C_SI;
  baryinput.site.location[2] = Detector->location[2]/LAL_C_SI;
  baryinput.alpha = CLA->Alpha;
  baryinput.delta = CLA->Delta;
  baryinput.dInv = 0.e0;

  /* amParams structure to compute a(t) and b(t) */

  /* Allocate space for amParams stucture */
  /* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  /* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = edat;
  amParams->das->pDetector = Detector; 
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->das->pSource->equatorialCoords.longitude = CLA->Alpha;
  amParams->das->pSource->equatorialCoords.latitude = CLA->Delta;
  amParams->das->pSource->orientation = 0.0;

  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  
  /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  TRY ( LALSCreateVector(status->statusPtr, &(amc.a), (UINT4)  CLA->nTsft), status);
  TRY ( LALSCreateVector(status->statusPtr, &(amc.b), (UINT4)  CLA->nTsft), status);
  
  /* Mid point of each SFT */
  midTS = (LIGOTimeGPS *)LALCalloc(CLA->nTsft,sizeof(LIGOTimeGPS));
  for(k=0; k < CLA->nTsft; k++)
    {
      /* FIXME:  loss of precision; consider
      midTS[k] = timestamps->data[k];
      XLALGPSAdd(&midTS[k], 0.5*CLA->Tsft);
      */
      REAL8 teemp=0.0;
      teemp = XLALGPSGetREAL8(&(timestamps->data[k]));
      teemp += 0.5*CLA->Tsft;
      XLALGPSSetREAL8(&(midTS[k]), teemp);
    }
  
  TRY ( LALComputeAM(status->statusPtr, &amc, midTS, amParams), status);

  /* Free memory */
  XLALDestroyTimestampVector ( timestamps);

  LALFree(midTS);
  LALFree(Detector);
  XLALDestroyEphemerisData(edat);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* ParseUserInput() */


/**
 * Check validity of user-input
 */
void
CheckUserInput (LALStatus *status,  struct CommandLineArgsTag *CLA )
{

  /* set a few abbreviations */
  BOOLEAN have_timestamps= XLALUserVarWasSet (&(CLA->timestamps));
  BOOLEAN have_gpsStart = XLALUserVarWasSet  (&(CLA->gpsStart));
  BOOLEAN have_duration  = XLALUserVarWasSet (&(CLA->duration));
  BOOLEAN have_nTsft     = XLALUserVarWasSet (&(CLA->nTsft));
  
  INITSTATUS(status);
  
  if( have_timestamps && (have_gpsStart||have_duration) )
    {
      fprintf(stderr,"\nBoth start time/duration and timestamps file specified - just need one !!\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }   
  
  if( !have_timestamps && !have_gpsStart )
    {
      fprintf(stderr,"\nNeed to specify gpsStart time or a timestamps file !!\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }
  
  if ( have_duration && have_nTsft )
    {
      fprintf (stderr, "\nSpecify only one of {duration, nTsft}!\n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }
  if ( !have_duration && !have_nTsft && !have_timestamps )
    {
      fprintf (stderr, "\nDuration has not been specified! Use one of {duration, nTsft, timestamps}!\n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }

  /* now one can either specify {h0, cosiota} OR {aPlus, aCross} */
  {
    BOOLEAN have_h0 = XLALUserVarWasSet (&(CLA->h0));
    BOOLEAN have_cosi = XLALUserVarWasSet (&(CLA->cosi));
    BOOLEAN have_aPlus = XLALUserVarWasSet (&(CLA->aPlus));
    BOOLEAN have_aCross = XLALUserVarWasSet (&(CLA->aCross));
    
    if ( (have_h0 || have_cosi) && (have_aPlus || have_aCross) ) 
      {
	fprintf (stderr, "\nSpecify EITHER {h0/cosiota} OR {aPlus/aCross}\n\n");
	ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
      }

    if ( (have_h0 && !have_cosi) || ( !have_h0 && have_cosi ) )
      {
	fprintf (stderr, "\nYou need to specify both --h0 and --cosi\n\n");
	ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
      }
    if ( (have_aPlus && !have_aCross) || ( !have_aPlus && have_aCross ) )
      {
	fprintf (stderr, "\nYou need to specify both --aPlus and --aCross\n\n");
	ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
      }
    if ( (CLA->cosi < -1.0) || (CLA->cosi > 1.0) )
      {
        fprintf (stderr, "\nIncorrect value for cos(iota)\n\n");
        ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
      }

    /* hack, hack... */
    if ( have_h0 )	/* internally only use aPlus, aCross */
      {
	CLA->aPlus = 0.5 * CLA->h0 * ( 1.0 + SQ ( CLA->cosi ) );
	CLA->aCross = CLA->h0 * CLA->cosi;
      }
  }

  RETURN (status);

} /* CheckUserInput() */
