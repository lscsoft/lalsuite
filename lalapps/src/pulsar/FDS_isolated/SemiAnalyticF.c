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

#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/LALStdio.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALComputeAM.h>

#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>

#include <lalapps.h>


RCSID("$Id$");

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
  CHAR  *efiles;
  REAL8 phi0;
  REAL8 psi;
  REAL8 sqrtSh;
  REAL8 duration;
  CHAR  *ephemYear;
  REAL8 aPlus;
  REAL8 aCross;
  REAL8 h0;
  REAL8 cosi;
  BOOLEAN help;
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

  lalDebugLevel = 0;	/* default value */
  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /*----------  read user-input and set up shop ----------*/
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);

  /* register all user-variables */
  LAL_CALL (InitUserVars(&status, &CommandLineArgs), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);
  if (CommandLineArgs.help)	/* help was called, do nothing here.. */
    return (0);
  LAL_CALL ( CheckUserInput (&status, &CommandLineArgs), &status);

  LAL_CALL ( Initialize (&status, &CommandLineArgs), &status);

  /*---------- central function: compute F-statistic ---------- */
  LAL_CALL ( ComputeF(&status, CommandLineArgs), &status); 

  
  /* Free remaining memory */
  LAL_CALL ( LALSDestroyVector(&status, &(amc.a) ), &status);
  LAL_CALL ( LALSDestroyVector(&status, &(amc.b) ), &status);
  LAL_CALL ( LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  return 0;

} /* main() */


void 
ComputeF( LALStatus *status, struct CommandLineArgsTag CLA)
{

  REAL8 A,B,C,D, A1,A2,A3,A4, To,Sh,F;
  REAL8 aPlus, aCross;
  REAL8 twopsi, twophi;

  INITSTATUS (status, "ComputeF", rcsid );
  ATTATCHSTATUSPTR ( status);
  
  A = amc.A;
  B = amc.B;
  C = amc.C;
  D = amc.D; 

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

  INITSTATUS( status, "InitUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* Initialize default values */
  CLA->Tsft=1800;
  CLA->nTsft=0;            
  CLA->timestamps=NULL;
  CLA->gpsStart=-1;
  CLA->efiles=NULL;
  CLA->sqrtSh=1.0;
  
  /** Default year-span of ephemeris-files to be used */
#define EPHEM_YEARS  "00-04"
  CLA->ephemYear = LALCalloc(1, strlen(EPHEM_YEARS)+1);
  strcpy (CLA->ephemYear, EPHEM_YEARS);
  
#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  CLA->efiles = LALCalloc(1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (CLA->efiles, DEFAULT_EPHEMDIR);
  
  CLA->help = FALSE;
  
  /* ---------- register all our user-variable ---------- */
  TRY (LALRegisterBOOLUserVar(status->statusPtr, "help", 'h', UVAR_HELP, "Print this message",
			      &(CLA->help)), status); 

  TRY( LALRegisterREALUserVar(status->statusPtr, "Alpha", 'a', UVAR_OPTIONAL, 
			      "Sky position Alpha (equatorial coordinates) in radians", 
			      &(CLA->Alpha)), status);
  TRY( LALRegisterREALUserVar(status->statusPtr, "longitude",  0, UVAR_DEVELOPER, 
			      "[DEPRECATED] Use --Alpha instead!",  
			      &(CLA->Alpha)), status);

  TRY( LALRegisterREALUserVar(status->statusPtr, "Delta", 'd', UVAR_OPTIONAL, 
			      "Sky position Delta (equatorial coordinates) in radians", 
			      &(CLA->Delta)), status);
  TRY( LALRegisterREALUserVar(status->statusPtr, "latitude", 0, UVAR_DEVELOPER, 
			      "[DEPRECATED] Use --Delta instead!", 
			      &(CLA->Delta)), status);
 
  TRY( LALRegisterREALUserVar(status->statusPtr, "phi0",   'Q', UVAR_OPTIONAL, 
			     "Phi_0: Initial phase in radians", &(CLA->phi0)), status);

  TRY( LALRegisterREALUserVar(status->statusPtr, "psi",    'Y', UVAR_OPTIONAL, 
			     "Polarisation in radians", &(CLA->psi)), status);

  TRY( LALRegisterREALUserVar(status->statusPtr, "cosi",   'i', UVAR_OPTIONAL, 
			      "Cos(iota)", &(CLA->cosi)), status);
  TRY( LALRegisterREALUserVar(status->statusPtr, "cosiota", 0, UVAR_DEVELOPER, 
			      "[DEPRECATED] Use --cosi instead", &(CLA->cosiota)), status);

  TRY( LALRegisterREALUserVar(status->statusPtr, "h0", 's', UVAR_OPTIONAL, 
			      "Strain amplitude h_0", &(CLA->h0)), status);
  TRY( LALRegisterREALUserVar(status->statusPtr, "sqrtSh", 'N', UVAR_OPTIONAL, 
			      "Noise floor: one-sided sqrt(Sh) in 1/sqrt(Hz)", &(CLA->sqrtSh)), status);
  
  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "timestampsFile", 'T', UVAR_OPTIONAL, 
				"Name of timestamps file", &(CLA->timestamps)), status);
  
  TRY( LALRegisterINTUserVar(status->statusPtr, "startTime", 'S', UVAR_OPTIONAL, 
			     "GPS start time of continuous observation", &(CLA->gpsStart)), status);
  
  TRY( LALRegisterREALUserVar(status->statusPtr, "Tsft", 't', UVAR_OPTIONAL, 
			      "Length of an SFT in seconds", &(CLA->Tsft)), status);
  
  TRY( LALRegisterINTUserVar(status->statusPtr, "nTsft", 'n', UVAR_OPTIONAL, 
			     "Number of SFTs", &(CLA->nTsft)), status);
  
  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "ephemDir", 'E', UVAR_OPTIONAL, 
				"Directory where Ephemeris files are located", 
				&(CLA->efiles)), status);

  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "IFO", 'D', UVAR_OPTIONAL, 
				"Detector: H1, H2, L1, G1, ... ",
				&(CLA->IFO)), status);
  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "detector",  0, UVAR_DEVELOPER, 
				"[DEPRECATED] Use --IFO instead!",
				&(CLA->detector)), status);
  
  /* ----- added for mfd_v4 compatibility ---------- */
  TRY ( LALRegisterREALUserVar(status->statusPtr, "duration", 0, UVAR_OPTIONAL,
			       "Duration of requested signal in seconds", 
			       &(CLA->duration)), status); 
  
  TRY ( LALRegisterSTRINGUserVar(status->statusPtr, "ephemYear", 0, UVAR_OPTIONAL,
				 "Year (or range of years) of ephemeris files to be used",
				 &(CLA->ephemYear)), status);
  
  TRY ( LALRegisterREALUserVar(status->statusPtr, "aPlus", 0, UVAR_OPTIONAL, 
			       "Plus polarization amplitude aPlus", 
			       &(CLA->aPlus)), status);
  TRY ( LALRegisterREALUserVar(status->statusPtr, "aCross", 0, UVAR_OPTIONAL, 
			       "Cross polarization amplitude aCross", 
			       &(CLA->aCross)), status);


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
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  LALDetector *Detector;              /* Our detector*/
  INT4 k;

  INITSTATUS (status, "Initialize", rcsid);
  ATTATCHSTATUSPTR (status);

  if ( LALUserVarWasSet ( &(CLA->nTsft) ) )
    CLA->duration = 1.0 * CLA->nTsft * CLA->Tsft;

  /* read or generate SFT timestamps */
  if ( LALUserVarWasSet(&(CLA->timestamps)) ) 
    { 
      TRY ( LALReadTimestampsFile (status->statusPtr, &timestamps, CLA->timestamps ), status );
      if ( (CLA->nTsft > 0) && ( (UINT4)CLA->nTsft < timestamps->length ) )	/* truncate if required */
	timestamps->length = CLA->nTsft;
      
      CLA->nTsft = timestamps->length;
    } /* if have_timestamps */
  else 
    {
      LIGOTimeGPS tStart;
      tStart.gpsSeconds = CLA->gpsStart;
      tStart.gpsNanoSeconds = 0;

      TRY ( LALMakeTimestamps(status->statusPtr, &timestamps, tStart, CLA->duration, CLA->Tsft ), status );
      CLA->nTsft = timestamps->length;

    } /* no timestamps */

  /*---------- initialize detector ---------- */
  {
    BOOLEAN have_IFO       = LALUserVarWasSet ( &CLA->IFO );
    BOOLEAN have_detector  = LALUserVarWasSet ( &CLA->detector );
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
#define MAXFILENAME 256
  {
    CHAR filenameE[MAXFILENAME], filenameS[MAXFILENAME];
    INT4 leap;
    
    /* don't use absolute path if none was given, this
     * allows LAL to find the ephemeris in LAL_DATA_PATH */
    if ( LALUserVarWasSet (&(CLA->efiles)) ) 
      {
	LALSnprintf (filenameE, MAXFILENAME, "%s/earth%s.dat", CLA->efiles, CLA->ephemYear );
	LALSnprintf (filenameS, MAXFILENAME, "%s/sun%s.dat", CLA->efiles, CLA->ephemYear );
      }
    else
      {
	LALSnprintf (filenameE, MAXFILENAME, "earth%s.dat", CLA->ephemYear );
	LALSnprintf (filenameS, MAXFILENAME, "sun%s.dat", CLA->ephemYear );
      }
    filenameE[MAXFILENAME-1] = 0;
    filenameS[MAXFILENAME-1] = 0;

    edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
    (*edat).ephiles.earthEphemeris = filenameE;     
    (*edat).ephiles.sunEphemeris = filenameS;         

    TRY ( LALLeapSecs(status->statusPtr, &leap, &(timestamps->data[0]), &formatAndAcc), status);
    (*edat).leap=leap; 

    /* Reads in ephemeris files */
    TRY( LALInitBarycenter (status->statusPtr, edat), status );

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
  amParams->leapAcc = formatAndAcc.accuracy;
  
  /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  TRY ( LALSCreateVector(status->statusPtr, &(amc.a), (UINT4)  CLA->nTsft), status);
  TRY ( LALSCreateVector(status->statusPtr, &(amc.b), (UINT4)  CLA->nTsft), status);
  
  /* Mid point of each SFT */
  midTS = (LIGOTimeGPS *)LALCalloc(CLA->nTsft,sizeof(LIGOTimeGPS));
  for(k=0; k < CLA->nTsft; k++)
    {
      REAL8 teemp=0.0;
      teemp = XLALGPSGetREAL8(&(timestamps->data[k]));
      teemp += 0.5*CLA->Tsft;
      XLALGPSSetREAL8(&(midTS[k]), teemp);
    }
  
  TRY ( LALComputeAM(status->statusPtr, &amc, midTS, amParams), status);

  /* Free memory */
  TRY ( LALDestroyTimestampVector (status->statusPtr, &timestamps), status );

  LALFree(midTS);
  LALFree(Detector);
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

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
  BOOLEAN have_timestamps= LALUserVarWasSet (&(CLA->timestamps));
  BOOLEAN have_gpsStart = LALUserVarWasSet  (&(CLA->gpsStart));
  BOOLEAN have_duration  = LALUserVarWasSet (&(CLA->duration));
  BOOLEAN have_nTsft     = LALUserVarWasSet (&(CLA->nTsft));
  
  INITSTATUS (status, "CheckUserInput", rcsid);
  
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
    BOOLEAN have_h0 = LALUserVarWasSet (&(CLA->h0));
    BOOLEAN have_cosi = LALUserVarWasSet (&(CLA->cosi));
    BOOLEAN have_aPlus = LALUserVarWasSet (&(CLA->aPlus));
    BOOLEAN have_aCross = LALUserVarWasSet (&(CLA->aCross));
    
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
