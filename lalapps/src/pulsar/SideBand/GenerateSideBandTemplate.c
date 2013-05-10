/*
 * Copyright (C) 2006 C. Messenger
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

/*********************************************************************************/
/** \author C. Messenger
 * \file 
 * \ingroup pulsarApps
 * \brief
 * Generate a frequency domain template for a binary signal as recieved at the SSB.
 *                                                                          
 *********************************************************************************/

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/UserInput.h>

#include <lal/LALInitBarycenter.h>
#include <lal/LALComputeAM.h>
#include <lalapps.h>
#include "SideBand.h"

/*---------- DEFINES ----------*/
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Error-codes -----*/
#define GENERATESIDEBANDTEMPLATEC_ENULL 		1
#define GENERATESIDEBANDTEMPLATEC_ESYS     		2
#define GENERATESIDEBANDTEMPLATEC_EINPUT   		3
#define GENERATESIDEBANDTEMPLATEC_EXLAL		        4

#define GENERATESIDEBANDTEMPLATEC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define GENERATESIDEBANDTEMPLATEC_MSGESYS		"System call failed (probably file IO)"
#define GENERATESIDEBANDTEMPLATEC_MSGEINPUT   		"Invalid input"
#define GENERATESIDEBANDTEMPLATEC_MSGEXLAL		"XLALFunction-call failed"


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
REAL8 uvar_f0;
REAL8 uvar_freqband;
REAL8 uvar_minfreq;
REAL8 uvar_semimajoraxis;
REAL8 uvar_orbitalperiod;
REAL8 uvar_eccentricity;
REAL8 uvar_argperiapse;
LIGOTimeGPS uvar_timeofperiapse;
REAL8 uvar_h0;
REAL8 uvar_cosi;
REAL8 uvar_psi;
REAL8 uvar_alpha;
REAL8 uvar_delta;
CHAR* uvar_ephemdir;
CHAR* uvar_ephemyear;
CHAR* uvar_ifo;
LIGOTimeGPS uvar_tstart;
LIGOTimeGPS uvar_reftime;
REAL8 uvar_tobs;
REAL8 uvar_sqrtSh;
REAL8 uvar_freqresfactor;
REAL8 uvar_phi0;
CHAR* uvar_output;
BOOLEAN uvar_help;
BOOLEAN uvar_local;
CHAR *uvar_timestampsfile;
CHAR *uvar_windowfile;
REAL8 uvar_highresfactor;
INT4 uvar_windowrange;
INT4 uvar_tsft;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *);
void checkUserInputConsistency (LALStatus *);

/*----------------------------------------------------------------------*/
/* Function definitions start here */
/*----------------------------------------------------------------------*/

/** 
 * MAIN function of GenerateBinaryFDTemplate code.
 * Compute the waveform of binary signal in the frequency domain as recieved at the SSB 
 */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */

  FILE *fp = NULL;                      /* pointer to output file */
  SideBandTemplate *Template = NULL;         /* structure for complex frequency domain results */
  SideBandTemplateParams *TParams = NULL;
  BinarySourceParams *BSParams = NULL;
  UINT4 Nbins;                          /* number of frequency bins */
  INT4 Mmin;
  INT4 Mmax;
  INT4 i;
  EphemerisData *edat;			    /**< ephemeris data (from LALInitBarycenter()) */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  LALDetector Detector;              /* Our detector*/
  EarthState earth;
  ABCcoParams abcparams;
  ABCcoefficients *ABCco = NULL;
  REAL8 df;
  EmissionTime emit;

  lalDebugLevel = 8;  
  vrbflg = 1;	/* verbose error-messages */
  
  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LAL_CALL (initUserVars(&status), &status); 	

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc,argv), &status);	

  if (uvar_help)	/* if help was requested, we're done here */
    exit (0);
  
  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(&status), &status);

  /* allocate memory to the template parameter structures */
  TParams = (SideBandTemplateParams *)LALMalloc(sizeof(SideBandTemplateParams));

  /*********************************************************************************************/
  /* read in time stamps from timestamps file */
  /*********************************************************************************************/
  if (uvar_timestampsfile!=NULL) ReadTimeStamps(&status,uvar_timestampsfile,uvar_tsft,&TParams);
  else {
    printf("ERROR : timestamps file required !!!, exiting\n");
    exit(0);
  }
  TParams->dfwindow = 1.0/(TParams->T*uvar_highresfactor);
  df = 1.0/(TParams->T*uvar_freqresfactor);
  TParams->windowrange = uvar_windowrange;

  for (i=0;i<(INT4)TParams->gapvectorstart->length;i++) {
    printf("%d -> %d\n",TParams->gapvectorstart->data[i],TParams->gapvectorend->data[i]);
  }
  
  if (lalDebugLevel) printf ("\nFinished Reading timestamps file.\n");
  
  /*********************************************************************************************/
  /* Load ephemeris data */
  /*********************************************************************************************/
  
  /* load ephemeris-data */
  edat = LALCalloc(1,sizeof(EphemerisData));
  InitEphemeris(&status,edat,uvar_ephemdir,uvar_ephemyear);

  if (lalDebugLevel) printf ("Finished initialising ephemeris.\n");

  /* select the detector */
  if(!strcmp(uvar_ifo,"G1")) {
    Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
    printf("DET = G1\n");
  }
  else if(!strcmp(uvar_ifo,"L1")) {
    Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
    printf("DET = L1\n");
  }
  else if(!strcmp(uvar_ifo,"H1")) {
    Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
    printf("DET = H1\n");
  }
  else if(!strcmp(uvar_ifo,"H2")) {
    Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
    printf("DET = H2\n");
  }
  else {
    XLALPrintError ("\nUnknown detector (use either G1,L1,H1,H2)!\n\n");
  }
  
  /* Detector location */
  baryinput.site.location[0] = Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = Detector.location[2]/LAL_C_SI;
  baryinput.alpha = uvar_alpha;
  baryinput.delta = uvar_delta;
  baryinput.dInv = 0.e0;
  baryinput.tgps.gpsSeconds = TParams->tstart.gpsSeconds;
  baryinput.tgps.gpsNanoSeconds = TParams->tstart.gpsNanoSeconds;
  
  /* compute SSB start time */
  LALBarycenterEarth(&status,&earth,&(TParams->tstart),edat);
  LALBarycenter(&status,&emit,&baryinput,&earth);
  printf("SSB tstart = %d %d\n",emit.te.gpsSeconds,emit.te.gpsNanoSeconds);
  printf("SSB deltaT = %6.12f\n",emit.deltaT);

  /* define a few variables */
  Nbins = (UINT4)ceil(TParams->T*uvar_freqband*uvar_freqresfactor);          /* define the number of frequency bins in frequency band */
  Mmin = (-1.0)*ceil(LAL_TWOPI*uvar_f0*uvar_semimajoraxis);             /* compute the expected number of upper sidebands */
  Mmax = ceil(LAL_TWOPI*uvar_f0*uvar_semimajoraxis);                    /* compute the expected number of lower sidebands */
 
  /* make sure that we only generate harmonics in the band */
  if (abs(Mmin)>(floor((uvar_f0-uvar_minfreq)*uvar_orbitalperiod))) Mmin = (-1.0)*floor((uvar_f0-uvar_minfreq)*uvar_orbitalperiod);
  if (Mmax>(floor((uvar_minfreq+uvar_freqband-uvar_f0)*uvar_orbitalperiod))) Mmax = floor((uvar_minfreq+uvar_freqband-uvar_f0)*uvar_orbitalperiod);

  if (lalDebugLevel) {
    printf("Changed : Mmin = %d\n",Mmin);
    printf("Changed : Mmax = %d\n",Mmax);
  }

  /* allocate memory to the template parameter structures */
  BSParams = (BinarySourceParams *)LALMalloc(sizeof(BinarySourceParams));

  if (lalDebugLevel) printf ("\nAllocated memory to the template parameter structures.\n");

  /* fill in structures to compute ABC coefficients */
  abcparams.tstart.gpsSeconds = TParams->tstart.gpsSeconds;
  abcparams.tstart.gpsNanoSeconds = TParams->tstart.gpsNanoSeconds;
  abcparams.alpha = uvar_alpha;
  abcparams.delta = uvar_delta;

  /* allocate memory for abccoefficients */
  ABCco = (ABCcoefficients *)LALMalloc(sizeof(ABCcoefficients));
  
  /* compute ABCcoefficients */
  ComputeABCcoefficients(&status,&abcparams,&Detector,&ABCco);

  printf("a co = %f %f %f %f %f %f\n",ABCco->aco[0],ABCco->aco[1],ABCco->aco[2],ABCco->apco[0],ABCco->apco[1],ABCco->apco[2]);
  printf("b co = %f %f %f %f %f %f\n",ABCco->bco[0],ABCco->bco[1],ABCco->bco[2],ABCco->bpco[0],ABCco->bpco[1],ABCco->bpco[2]);

  /* fill in the input parameters for BinaryFDTemplate */  
  if (uvar_sqrtSh==0.0) TParams->sqrtSh = 1.0;
  else TParams->sqrtSh = uvar_sqrtSh;
  TParams->Mmin = Mmin;
  TParams->Mmax = Mmax;
  TParams->tstartSSB.gpsSeconds = emit.te.gpsSeconds;
  TParams->tstartSSB.gpsNanoSeconds = emit.te.gpsNanoSeconds;
  TParams->reftime.gpsSeconds = uvar_reftime.gpsSeconds;
  TParams->reftime.gpsNanoSeconds = uvar_reftime.gpsNanoSeconds;
  printf("reftime = %d %d\n",TParams->reftime.gpsSeconds,TParams->reftime.gpsNanoSeconds);    
  TParams->local = uvar_local;
  TParams->ABCco = ABCco;
  TParams->pmax = 4;

  /* fill in frequency samples */
  /* in this case we want a highly resolved frequency */
  TParams->freqsamples = XLALCreateREAL8Vector(Nbins);
  for (i=0;i<(INT4)Nbins;i++) TParams->freqsamples->data[i] = uvar_minfreq + i*df;

  BSParams->OrbitalSemiMajorAxis = uvar_semimajoraxis;
  BSParams->OrbitalPeriod = uvar_orbitalperiod;
  BSParams->OrbitalEccentricity = uvar_eccentricity;
  BSParams->ArgumentofPeriapse = uvar_argperiapse;
  BSParams->TimeofSSBPeriapsePassage.gpsSeconds = uvar_timeofperiapse.gpsSeconds;
  BSParams->TimeofSSBPeriapsePassage.gpsNanoSeconds = uvar_timeofperiapse.gpsNanoSeconds;
  BSParams->alpha = uvar_alpha;
  BSParams->delta = uvar_delta;
  BSParams->f0 = uvar_f0;
  BSParams->phi0 = uvar_phi0;
  BSParams->psi = uvar_psi;
  BSParams->cosi = uvar_cosi;
  BSParams->h0 = uvar_h0;

  if (lalDebugLevel) printf ("\nFilled in the template parameter structures.\n");

  /* allocate memory for the window */
  TParams->wa = NULL;
  TParams->wb = NULL;
  /* TParams->wa = XLALCreateCOMPLEX16Vector(2*TParams->windowrange);
     TParams->wb = XLALCreateCOMPLEX16Vector(2*TParams->windowrange); */
 
  /* compute the window function */
  ComputeSideBandWindow(&status,ABCco,uvar_windowfile,&TParams);

  Template = (SideBandTemplate *)LALMalloc(sizeof(SideBandTemplate));
  Template->fourier = NULL;
  Template->minfreq = TParams->freqsamples->data[0];
  Template->length = TParams->freqsamples->length;
  Template->epoch.gpsSeconds = TParams->tstart.gpsSeconds;
  Template->epoch.gpsNanoSeconds = TParams->tstart.gpsNanoSeconds;
  Template->fourier = XLALCreateCOMPLEX16Vector(TParams->freqsamples->length);

  if (lalDebugLevel) printf ("\nFinished allocating memory to the template results structures.\n");
  
  /* call the template generating function */
  GenerateSideBandTemplate(&status,BSParams,TParams,&Template);
  
  if (lalDebugLevel) printf ("\nFinished SideBandTemplate.\n");

  /* check if output file is OK */
  if ((fp = fopen(uvar_output,"w"))==NULL) {
    XLALPrintError("\nError opening file '%s' for writing..\n\n",uvar_output);
    return (GENERATESIDEBANDTEMPLATEC_ESYS);
  }

  /* output the results to file */
  for (i=0;i<(INT4)Nbins;i++) {
    fprintf(fp,"%16.12f %6.9f %6.9f\n",
	    TParams->freqsamples->data[i],
	    creal(Template->fourier->data[i]),
	    cimag(Template->fourier->data[i]));
  }
  fclose(fp);
    
  if (lalDebugLevel) printf ("\nFinished outputting to file.\n");

  /* Free memory */ 
  XLALDestroyREAL8Vector(TParams->freqsamples);
  XLALDestroyCOMPLEX16Vector(Template->fourier);
  LALZDestroyVector(&status,&(TParams->wa));
  LALZDestroyVector(&status,&(TParams->wb));
  LALFree(BSParams);
  LALFree(Template);
  LALFree(TParams);
  LALI4DestroyVector(&status,&(TParams->timestamps));
  LALI4DestroyVector(&status,&(TParams->gapvectorstart));
  LALI4DestroyVector(&status,&(TParams->gapvectorend));
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LALFree(ABCco);

  /* Free config-Variables and userInput stuff */
  LALDestroyUserVars(&status);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
    
  return 0;
  
} /* main() */


/** 
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar_semimajoraxis = 0.0;
  uvar_orbitalperiod = 0.0;
  uvar_timeofperiapse.gpsSeconds = 0;
  uvar_timeofperiapse.gpsNanoSeconds = 0;
  uvar_eccentricity = 0.0;
  uvar_argperiapse = 0.0;
  uvar_h0 = 0.0;
  uvar_cosi = 0.0;
  uvar_psi = 0.0;
  uvar_phi0 = 0.0;
  uvar_f0 = 0.0;
  uvar_freqband = 0.0;
  uvar_minfreq = 0.0;
  uvar_freqresfactor = 1.0;
  uvar_tobs = 0.0;
  uvar_tstart.gpsSeconds = 0;
  uvar_tstart.gpsNanoSeconds = 0;
  uvar_reftime.gpsSeconds = 0;
  uvar_reftime.gpsNanoSeconds = 0;
  uvar_alpha = 0.0;
  uvar_delta = 0.0;
  uvar_ephemdir = NULL;
  uvar_ephemyear = NULL;
  uvar_ifo = NULL;
  uvar_help = FALSE;
  uvar_output = NULL;
  uvar_windowfile = NULL;
  uvar_local = 0;
  uvar_sqrtSh = 0.0;
  uvar_highresfactor = 10.0;
  uvar_windowrange = 1000;
  uvar_tsft = 1800;

  /* register all our user-variables */
  LALregBOOLUserVar(status, 	help, 		                      'h', UVAR_HELP,     "Print this message"); 
  LALregREALUserVar(status, 	semimajoraxis, 		              'a', UVAR_REQUIRED, "Orbital semi-major axis normalised by c in Seconds");
  LALregREALUserVar(status, 	orbitalperiod, 	                      'b', UVAR_REQUIRED, "Orbital period in Seconds");
  LALregINTUserVar(status,      timeofperiapse.gpsSeconds,            'c', UVAR_REQUIRED, "Time of periapse passage (Seconds part)");
  LALregINTUserVar(status, 	timeofperiapse.gpsNanoSeconds, 	      'd', UVAR_REQUIRED, "Time of periapse passage (NanoSeconds part)");
  LALregREALUserVar(status, 	eccentricity, 	                      'e', UVAR_REQUIRED, "Orbital eccentricity");
  LALregREALUserVar(status, 	argperiapse, 	                      'f', UVAR_REQUIRED, "Orbital argument of periapse in radians");
  LALregREALUserVar(status, 	h0,   	                              'i', UVAR_REQUIRED, "Gravitational wave amplitude");
  LALregREALUserVar(status, 	cosi,   	                      'j', UVAR_REQUIRED, "The cosine of the source inclination angle");
  LALregREALUserVar(status, 	psi,    	                      'k', UVAR_REQUIRED, "The gravitational wave polarisation angle in radians");
  LALregREALUserVar(status, 	f0, 	                              'l', UVAR_REQUIRED, "Intrinsic gravitational wave frequency in Hz");
  LALregREALUserVar(status, 	minfreq, 	                      'm', UVAR_REQUIRED, "Start frequency of band in Hz");
  LALregREALUserVar(status, 	freqband, 		              'n', UVAR_REQUIRED, "Width of frequency band in Hz");
  LALregREALUserVar(status, 	tobs, 	                              'o', UVAR_REQUIRED, "Total observation span in Seconds");
  LALregREALUserVar(status, 	sqrtSh, 	                      'p', UVAR_REQUIRED, "Sqrt of Noise spectral density in Hz^{-1/2} (0 = no noise)");
  LALregINTUserVar(status,      tstart.gpsSeconds,                    'q', UVAR_REQUIRED, "Time of observation start (Seconds part)");
  LALregINTUserVar(status, 	tstart.gpsNanoSeconds,       	      'r', UVAR_REQUIRED, "Time of observation start (NanoSeconds part)");
  LALregINTUserVar(status,      reftime.gpsSeconds,                   's', UVAR_REQUIRED, "Reference time at which pulsar parameters are defined (Seconds part)");
  LALregINTUserVar(status, 	reftime.gpsNanoSeconds,       	      't', UVAR_REQUIRED, "Reference time at which pulsar parameters are defined (NanoSeconds part)");
  LALregREALUserVar(status, 	phi0, 	                              'u', UVAR_REQUIRED, "Gravitational wave phase at start of observation in radians");
  LALregREALUserVar(status, 	alpha, 	                              'v', UVAR_REQUIRED, "Source right ascension in radians");
  LALregREALUserVar(status, 	delta, 	                              'w', UVAR_REQUIRED, "Source declination in radians");
  LALregSTRINGUserVar(status, 	ephemdir, 	                      'x', UVAR_REQUIRED, "Location of ephemeris files");
  LALregSTRINGUserVar(status, 	ephemyear, 	                      'y', UVAR_REQUIRED, "Ephemeris file year");
  LALregSTRINGUserVar(status, 	ifo, 	                              'z', UVAR_REQUIRED, "Interferometer name (G1,L1,H1,H2)");
  LALregREALUserVar(status, 	freqresfactor, 		              'A', UVAR_OPTIONAL, "Factor by which to over-resolve in frequency");
  LALregSTRINGUserVar(status, 	output, 	                      'B', UVAR_REQUIRED, "Name of output file");
  LALregBOOLUserVar(status, 	local,                                'E', UVAR_OPTIONAL, "Set flag if you wish to use only contibution from local sideband");
  LALregSTRINGUserVar(status, 	timestampsfile,                       'F', UVAR_REQUIRED, "Name of input time stamps file");
  LALregSTRINGUserVar(status, 	windowfile,                           'K', UVAR_REQUIRED, "Name of output window file");
  LALregREALUserVar(status, 	highresfactor, 	                      'G', UVAR_REQUIRED, "Over resolution factor for window function calculation (default = 10)");
  LALregINTUserVar(status,      windowrange,                          'I', UVAR_REQUIRED, "Range in highly resolved bins used to compute the window (default 1000)");
  LALregINTUserVar(status,      tsft,                                 'J', UVAR_OPTIONAL, "Length of SFT in seconds (default 1800)");

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */


/*----------------------------------------------------------------------*/
/** Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *status)
{

  INITSTATUS(status);  

  /* don't allow negative orbital sem-major axis */
  if (uvar_semimajoraxis < 0.0) {
    XLALPrintError ("\nOrbital semi-major axis must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }      
  /* don't allow zero or negative orbital periods */
  if (uvar_orbitalperiod <= 0.0) {
    XLALPrintError ("\nOrbital period must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative periapse passage times */
  if (uvar_timeofperiapse.gpsSeconds < 0) {
    XLALPrintError ("\nTime of periapse passage must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative periapse passage nanosecond parts */
  if (uvar_timeofperiapse.gpsNanoSeconds < 0) {
    XLALPrintError ("\nTime of periapse passage nanosecond part must be >= 0!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative observation start times */
  if (uvar_tstart.gpsSeconds < 0) {
    XLALPrintError ("\nTime of observation start must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative observation start time nanosecond parts */
  if (uvar_tstart.gpsNanoSeconds < 0) {
    XLALPrintError ("\nTime of observation start nanosecond part must be >= 0!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative frequency */
  if (uvar_f0 <= 0) {
    XLALPrintError ("\nFrequency must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative frequency band */
  if (uvar_freqband <= 0) {
    XLALPrintError ("\nFrequency band must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative minimum frequency */
  if (uvar_minfreq <= 0) {
    XLALPrintError ("\nMinimum frequency must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow phi0 out of range */
  if ((uvar_phi0 < 0) || (uvar_phi0 >= LAL_TWOPI)) {
    XLALPrintError ("\nPhi0 must be in range [0,2*pi)!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative frequency resolution factor */
  if (uvar_freqresfactor <= 0) {
    XLALPrintError ("\nFrequency resolution factor must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow negative observation time */
  if (uvar_tobs <= 0) {
    XLALPrintError ("\nObservation time must be positive!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow unspecified output file */
  if (uvar_output == NULL) {
    XLALPrintError ("\nOutput file must be specified (option 'output')\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow right ascension out of range */
  if ((uvar_alpha < 0)||(uvar_alpha >= LAL_TWOPI)) {
    XLALPrintError ("\nRight ascension must be in range [0,2PI)!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  /* don't allow declination out of range */
  if ((uvar_delta < ((-1.0)*LAL_PI/2.0))||(uvar_delta > (LAL_PI/2.0))) {
    XLALPrintError ("\nDeclination must be in range [-PI/2,PI/2]!\n\n");
    ABORT (status, GENERATESIDEBANDTEMPLATEC_EINPUT, GENERATESIDEBANDTEMPLATEC_MSGEINPUT);
  }
  
  
  RETURN (status);
} /* checkUserInputConsistency() */





