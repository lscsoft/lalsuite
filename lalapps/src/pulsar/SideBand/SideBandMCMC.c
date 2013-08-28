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
/**
 * \author C. Messenger
 * \file
 * \ingroup pulsarApps
 * \brief
 * Generates posterior pdfs for a subset of the unknown orbital and nuisance
 * parameters given a set of candidate regions in frequency of demodulated Fourier
 * transform results.
 *
 */

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#include <time.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALComputeAM.h>
#include <lal/Random.h>
#include <lalapps.h>
#include "SideBand.h"

/*---------- DEFINES ----------*/
#define MAXUINT4 2147483647

/*----- Error-codes -----*/
#define SIDEBANDMCMCC_ENULL 		1
#define SIDEBANDMCMCC_ESYS     		2
#define SIDEBANDMCMCC_EINPUT   		3
#define SIDEBANDMCMCC_EXLAL	        4

#define SIDEBANDMCMCC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define SIDEBANDMCMCC_MSGESYS		"System call failed (probably file IO)"
#define SIDEBANDMCMCC_MSGEINPUT       	"Invalid input"
#define SIDEBANDMCMCC_MSGEXLAL		"XLALFunction-call failed"


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
CHAR* uvar_ephemdir;
CHAR* uvar_ephemyear;
CHAR* uvar_ifo;
REAL8 uvar_alpha;
REAL8 uvar_delta;
LIGOTimeGPS uvar_reftime;
CHAR* uvar_output;
CHAR* uvar_fourierfile;
CHAR* uvar_rangefile;
CHAR* uvar_windowfile;
BOOLEAN uvar_help;
INT4 uvar_seed;
INT4 uvar_Nchain;
INT4 uvar_Nburn;
REAL8 uvar_sqrtSh; 
CHAR* uvar_stampsfile; 
REAL8 uvar_highresfactor;
INT4 uvar_tsft;
BOOLEAN uvar_noneccentric;
BOOLEAN uvar_ampmod;

REAL8 jumptemp;

static SideBandMCMCVector empty_SideBandMCMCVector;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *);
void checkUserInputConsistency (LALStatus *);
void TestJump(LALStatus *,SideBandMCMCVector *,SideBandMCMCVector *,SideBandMCMCVector *,RandomParams *,INT4,INT4,INT4 *);
void MakeJump(LALStatus *,SideBandMCMCVector,SideBandMCMCVector *,SideBandMCMCJumpProbs,SideBandMCMCRanges,RandomParams *);
void InitialiseLambda(LALStatus *,SideBandMCMCVector *,SideBandMCMCRanges,RandomParams *);
void PriorRanges(LALStatus *,SideBandMCMCVector ,SideBandMCMCRanges,INT4 *);

/*----------------------------------------------------------------------*/
/* Function definitions start here */
/*----------------------------------------------------------------------*/

/**
 * MAIN function of SideBandMCMC code
 * Compute the posterior pdfs of the orbital and nuisance parameters of a binary signal
 * in Fstat form
 */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */
  FILE *fp = NULL;                      /* pointer to output file */
  SideBandTemplate *Template = NULL;          /* structure for complex frequency domain template */
  SideBandTemplateParams *TParams = NULL;       /* structure for template parameters */
  INT4 i;                           /* general loop indexes */
  EphemerisData *edat;			/* ephemeris data (from LALInitBarycenter()) */
  BarycenterInput baryinput;            /* Stores detector location and other barycentering data */
  LALDetector Detector;                 /* Our detector */
  EarthState earth;
  ABCcoParams abcparams;
  ABCcoefficients *ABCco = NULL;
  SideBandMCMCVector **MCMCchain = NULL;       /* a vector of MCMC parameter vectors where the Markov chain is stored */
  SideBandMCMCVector lambda;                   /* stores the current MCMC parameter vector */
  SideBandMCMCVector newlambda;                /* stores the prospective MCMC parameter vector */
  SideBandMCMCVector currentlambda = empty_SideBandMCMCVector; /* stores the last successful MCMC jump parameter vector */
  SideBandMCMCRanges ranges;                   /* used to store the prior ranges on all MCMC parameters */
  SideBandMCMCJumpProbs jumpsizes;                 /* used to store the jump sizes for all MCMC parameters */
  RandomParams *randparams = NULL;
  REAL8 safety;
  EmissionTime emit;
  INT4 pmax = 4; 
  INT4 out;
  REAL8 temp2;
  INT4 AC = 0;
  INT4 ACtemp = 0;
  INT4 flag;
  ReadSideBandDataParams RFDparams;
  SideBandDataSet *fulldata = NULL;
  SideBandDataSet *reddata = NULL;
  REAL8 minfreq,maxfreq;
  SelectSideBandFrequencyParams SFparams;
  EstimateSideBandNoiseParams ENparams;

  vrbflg = 0;	/* verbose error-messages */

  /**************************************************************************************/
  /* do some standard things to start */
  
  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (initUserVars(&status), &status); 	

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc,argv), &status);	

  if (uvar_help)	/* if help was requested, we're done here */
    exit (0);
  
  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(&status), &status);

  if (lalDebugLevel) printf ("\nFinished checking input.\n");

  /* allocate memory to the template parameter structures */
  TParams = (SideBandTemplateParams *)LALMalloc(sizeof(SideBandTemplateParams));
  
  /**************************************************************************************/
  /* read in time stamps from timestamps file */
  /**************************************************************************************/
  
  if (uvar_stampsfile!=NULL) ReadTimeStamps(&status,uvar_stampsfile,uvar_tsft,&TParams);
  else {
    printf("ERROR : timestamps file required !!!, exiting\n");
    exit(0);
  }
  TParams->dfwindow = 1.0/(TParams->T*uvar_highresfactor);
  if (uvar_ampmod) TParams->windowrange = 4.0/(LAL_DAYSID_SI*TParams->dfwindow);
  else TParams->windowrange = (INT4)uvar_highresfactor;
  
  for (i=0;i<(INT4)TParams->gapvectorstart->length;i++) {
    printf("%d -> %d\n",TParams->gapvectorstart->data[i],TParams->gapvectorend->data[i]);
  } 

  if (lalDebugLevel) printf ("\nFinished Reading timestamps file.\n");

  /**************************************************************************************/
  /* here we read in the MCMC parameter ranges from file */
  /**************************************************************************************/
  
  ReadSideBandPriors(&status,uvar_rangefile,&ranges,&jumpsizes);

  printf("f0min = %6.12f f0max = %6.12f\n",ranges.f0min,ranges.f0max);

   if (lalDebugLevel) printf ("\nFinished Reading Priors and Proposals file.\n");

  /**************************************************************************************/
  /* compute values of m and M given prior ranges */
  /**************************************************************************************/
  
  /* compute maximum range of M given amax */
  safety = 1.2;
  TParams->mmin = (INT4)-ceil(safety*LAL_TWOPI*ranges.f0max*(ranges.amax+ranges.amin)/2.0);
  TParams->mmax = (INT4)ceil(safety*LAL_TWOPI*ranges.f0max*(ranges.amax+ranges.amin)/2.0);
  TParams->Mmax = TParams->mmax - TParams->mmin + 1;
  minfreq = ranges.f0min + TParams->mmin/ranges.periodmin - 2.0/LAL_DAYSID_SI;
  maxfreq = ranges.f0max + TParams->mmax/ranges.periodmin + 2.0/LAL_DAYSID_SI;
  
  printf("mmin = %d mmax = %d Mmax = %d\n",TParams->mmin,TParams->mmax,TParams->Mmax);
  printf("minfreq = %6.12f maxfreq = %6.12f\n",minfreq,maxfreq);

  if (lalDebugLevel) printf("*** Finished computing number of spikes\n");

  /**************************************************************************************/
  /* here we read in the relevant Fourier data */
  /**************************************************************************************/

  RFDparams.minf = minfreq;
  RFDparams.maxf = maxfreq;
  RFDparams.Tobs = TParams->Tobs;
  RFDparams.nsft = TParams->nsft;
  sprintf(RFDparams.file,"%s",uvar_fourierfile);
  fulldata = (SideBandDataSet *)LALMalloc(sizeof(SideBandDataSet));
  ReadSideBandData(&status,&RFDparams,&fulldata);
    
  if (lalDebugLevel) printf("*** Finished reading Fourier data\n");

  /************************************************************************/
  /* now we select the data corresponding to each FM and AM sideband      */
  /************************************************************************/

  SFparams.ranges = ranges;
  SFparams.mmin = TParams->mmin;
  SFparams.mmax = TParams->mmax;
  SFparams.am = uvar_ampmod;
  SFparams.df = RFDparams.df;
  reddata = (SideBandDataSet *)LALMalloc(sizeof(SideBandDataSet));

  SelectSideBandFrequencies(&status,&fulldata,&SFparams,&reddata);

  if (lalDebugLevel) printf("*** Finished reducing Fourier data\n");

  /************************************************************************/
  /* now we estimate noise from data */
  /************************************************************************/

  ENparams.minfreq = minfreq;
  ENparams.maxfreq = maxfreq;
  ENparams.minf = SFparams.minf;
  ENparams.maxf = SFparams.maxf;
  ENparams.Nthresh = 10;
  ENparams.safedf = 2.0/TParams->T;

  EstimateSideBandNoise(&status,&fulldata,&ENparams); 

  if (lalDebugLevel) printf("*** Finished estimating noise\n");

  /************************************************************************/
  /* here we set up the ephemeris */
  /************************************************************************/

  /* load ephemeris-data */
  edat = LALCalloc(1,sizeof(EphemerisData));
  InitEphemeris(&status,edat,uvar_ephemdir,uvar_ephemyear);

  if (lalDebugLevel) printf ("Finished initialising ephemeris.\n");

  /************************************************************************/
  /* here we set up the detector */
  /************************************************************************/

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

  /************************************************************************************/
  /* here we compute the amplitude modulation coeeficients */
  /************************************************************************************/

  /* fill in structures to compute ABC coefficients */
  abcparams.tstart.gpsSeconds = TParams->tstart.gpsSeconds;
  abcparams.tstart.gpsNanoSeconds = TParams->tstart.gpsNanoSeconds;
  abcparams.alpha = uvar_alpha;
  abcparams.delta = uvar_delta;

  /* allocate memory for abccoefficients */
  ABCco = (ABCcoefficients *)LALMalloc(sizeof(ABCcoefficients));
  
  /* compute ABCcoefficients */
  ComputeABCcoefficients(&status,&abcparams,&Detector,&ABCco);

  if (lalDebugLevel) printf ("\nFinished computing amplitude modulation coefficients.\n");

  /************************************************************************************/
  /* here we compute the window functions */
  /************************************************************************************/

  /* allocate memory for the window */
  TParams->wa = NULL;
  TParams->wb = NULL;
  
  /* compute the window function */
  ComputeSideBandWindow(&status,ABCco,uvar_windowfile,&TParams);

  if (lalDebugLevel) printf ("\nFinished computing the window functions.\n");

  /************************************************************************************/
  /* here we set up the template parameter structure and allocate some memory */
  /************************************************************************************/

  /* fill in the input parameters for BinaryFDTemplate */
  /* these remain fixed for the duration of the MCMC (at present) */
  /* if (uvar_sqrtSh==0.0) TParams->sqrtSh = 1.0;
     else TParams->sqrtSh = uvar_sqrtSh; */
  TParams->sqrtSh = ENparams.sqrtSh;
  TParams->tstartSSB.gpsSeconds = emit.te.gpsSeconds;
  TParams->tstartSSB.gpsNanoSeconds = emit.te.gpsNanoSeconds;
  if (uvar_reftime.gpsSeconds+1e-9*uvar_reftime.gpsNanoSeconds==0.0) {
    uvar_reftime.gpsSeconds = TParams->tstartSSB.gpsSeconds;
    uvar_reftime.gpsNanoSeconds = TParams->tstartSSB.gpsNanoSeconds;
  } 
  TParams->reftime.gpsSeconds = uvar_reftime.gpsSeconds;
  TParams->reftime.gpsNanoSeconds = uvar_reftime.gpsNanoSeconds;
  printf("reftime = %d %d\n",TParams->reftime.gpsSeconds,TParams->reftime.gpsNanoSeconds);
  TParams->local = 1;
  if (uvar_noneccentric) pmax = 0;
  TParams->pmax = pmax;
  TParams->ABCco = ABCco;
  TParams->freqsamples = reddata->freq;
  
  if (lalDebugLevel) printf ("\nFinished setting up template parameter structure.\n");

  /************************************************************************************/
  /* allocate memory for template results */
  /************************************************************************************/

  Template = (SideBandTemplate *)LALMalloc(sizeof(SideBandTemplate));
  Template->fourier = NULL;
  Template->minfreq = TParams->freqsamples->data[0];
  Template->length = TParams->freqsamples->length;
  Template->epoch.gpsSeconds = TParams->tstart.gpsSeconds;
  Template->epoch.gpsNanoSeconds = TParams->tstart.gpsNanoSeconds;
  Template->fourier = XLALCreateCOMPLEX16Vector(TParams->freqsamples->length);
 
  /* allocate memory to the chain */
  MCMCchain = (SideBandMCMCVector **)LALCalloc(uvar_Nchain,sizeof(SideBandMCMCVector *));
  for (i=0;i<uvar_Nchain;i++) MCMCchain[i] = (SideBandMCMCVector *)LALCalloc(1,sizeof(SideBandMCMCVector));
 
  if (lalDebugLevel) printf ("\nFinished allocating memory for MCMC results.\n");

  /************************************************************************************/
  /* generate random number parameters */
  /************************************************************************************/
  
  randparams = (RandomParams *)LALMalloc(sizeof(RandomParams));
  randparams = NULL;
  LALCreateRandomParams(&status,&randparams,uvar_seed);

  if (lalDebugLevel) printf ("\nFinished setting up Template parameters.\n");

  /************************************************************************************/
  /* begin the MCMC part of the code */
  /************************************************************************************/

  /* define initial random parameter values */
  InitialiseLambda(&status,&lambda,ranges,randparams);
  
  printf("initialised lambda as [%f %f %f %d %d %f %f %e %f %f %f]\n",
	 lambda.f0,lambda.period,lambda.a,lambda.tp.gpsSeconds,
	 lambda.tp.gpsNanoSeconds,lambda.argp,lambda.e,lambda.h0,lambda.cosi,lambda.psi,lambda.phi0); 

  /* compute likelihood */
  ComputeSideBandLikelihood(&status,&lambda,reddata,&Template,TParams);

  /* newlambda.f0 = 200.0; */
  /* printf("likelihood = %e\n",lambda.logL); */
 
  temp2 = uvar_Nchain/100.0;

  /* record very first point */
  *(MCMCchain[0]) = lambda; 

  /* start a loop over number of links in the chain */
  for (i=1;i<(INT4)uvar_Nchain;i++) {
    
    /* jump to a new location */
    /* MakeJumpF(&status,lambda,&newlambda,jumpsizes,ranges,randparams); */
    MakeJump(&status,lambda,&newlambda,jumpsizes,ranges,randparams);
  
    /* printf("jumped lambda to [%f %f %d %d %f %f %f %f %f %f]\n",
	   newlambda.f0,newlambda.a,newlambda.tp.gpsSeconds,newlambda.tp.gpsNanoSeconds,
	   newlambda.argp,newlambda.e,newlambda.h0,newlambda.cosi,newlambda.psi,newlambda.phi0); */
    
    /* if we have jumped within our prior ranges */
    PriorRanges(&status,newlambda,ranges,&out);
    
    /* compute likelihood */
    ComputeSideBandLikelihood(&status,&newlambda,reddata,&Template,TParams);      
    
    /* Test the new location */
    flag = 0;
    TestJump(&status,&lambda,&newlambda,MCMCchain[i],randparams,i,out,&flag);
      
    /* if it's a good jump then move to new location and record acceptance increment */
    if (flag==1) {
      AC++; 
      ACtemp++;
      currentlambda = lambda;
    }
    
    /* printf("decided lambda is [%f %d %d %f %f %f %f %f %f]\n",
       MCMCchain[i]->a,MCMCchain[i]->tp.gpsSeconds,MCMCchain[i]->tp.gpsNanoSeconds,MCMCchain[i]->argp,MCMCchain[i]->e,MCMCchain[i]->h0,MCMCchain[i]->cosi,MCMCchain[i]->psi,MCMCchain[i]->phi0);  */
    
    if (fmod(i,temp2)==0) {
      printf("*** %6.3f %% complete *** [acceptance ratio = %6.3f(%6.3f)] logL = %f\n",
	     100*(REAL8)i/(REAL8)uvar_Nchain,(REAL8)AC/(REAL8)i,(REAL8)ACtemp/(REAL8)temp2,MCMCchain[i]->logL);
      printf("currently at [%16.12f %16.12f %f %d %d %f %f %e %f %f %f]\n",
	     currentlambda.f0,currentlambda.period,currentlambda.a,currentlambda.tp.gpsSeconds,currentlambda.tp.gpsNanoSeconds,
	     currentlambda.argp,currentlambda.e,currentlambda.h0,currentlambda.cosi,currentlambda.psi,currentlambda.phi0);
      ACtemp = 0;
    }

  }
  

  /************************************************************************************/
  /* output results to file */
  
  /* check if output file is OK */
  if ((fp = fopen(uvar_output,"w"))==NULL) {
    XLALPrintError("\nError opening file '%s' for writing..\n\n",uvar_output);
    return (SIDEBANDMCMCC_ESYS);
  }
  
  for (i=0;i<(INT4)uvar_Nchain;i++) {
    fprintf(fp,"%16.12f %16.12f %16.12f %d %d %6.12f %6.12f %6.12e %6.12f %6.12f %6.12f %6.12f\n",
	    MCMCchain[i]->f0,
	    MCMCchain[i]->period,
	    MCMCchain[i]->a,
	    MCMCchain[i]->tp.gpsSeconds,
	    MCMCchain[i]->tp.gpsNanoSeconds,
	    MCMCchain[i]->argp,
	    MCMCchain[i]->e,
	    MCMCchain[i]->h0,
	    MCMCchain[i]->cosi,
	    MCMCchain[i]->psi,
	    MCMCchain[i]->phi0,
	    MCMCchain[i]->logL);
  }
  fclose(fp);

  
  if (lalDebugLevel) printf ("\nFinished outputting to file.\n");
  
  /************************************************************************************/
  /* Free memory */ 
  
  XLALDestroyREAL8Vector(TParams->freqsamples);
  XLALDestroyREAL8Vector(reddata->freq);
  XLALDestroyCOMPLEX16Vector(reddata->fourier);
 
  LALFree(Template);
  LALFree(TParams);
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


/********************************************************************************************/
/* Register all our "user-variables" that can be specified from cmd-line and/or config-file.
* Here we set defaults for some user-variables and register them with the UserInput module.
*/
void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /* set a few defaults */
  uvar_alpha = 0.0;
  uvar_delta = 0.0;
  uvar_ephemdir = NULL;
  uvar_ephemyear = NULL;
  uvar_ifo = NULL;
  uvar_reftime.gpsSeconds = 0;
  uvar_reftime.gpsNanoSeconds = 0;
  uvar_output = NULL;
  uvar_fourierfile = NULL;
  uvar_seed = 0;
  uvar_Nchain = 0;
  uvar_sqrtSh = 1.0; 
  uvar_rangefile = NULL;
  uvar_Nburn = 0;
  uvar_stampsfile = NULL;
  uvar_windowfile = NULL;
  uvar_highresfactor = 100;
  uvar_tsft = 1800;
  uvar_noneccentric = 0;
  uvar_ampmod = 0;

  /* register all our user-variables */
  LALregBOOLUserVar(status, 	help, 		                      'h', UVAR_HELP,     "Print this message"); 
  LALregINTUserVar(status,      reftime.gpsSeconds,                   'A', UVAR_REQUIRED, "Reference time at which pulsar parameters are defined (Seconds part)");
  LALregINTUserVar(status, 	reftime.gpsNanoSeconds,       	      'B', UVAR_REQUIRED, "Reference time at which pulsar parameters are defined (NanoSeconds part)");
  LALregREALUserVar(status, 	alpha, 	                              'i', UVAR_REQUIRED, "Source right ascension in radians");
  LALregREALUserVar(status, 	delta, 	                              'j', UVAR_REQUIRED, "Source declination in radians");
  LALregSTRINGUserVar(status, 	ephemdir, 	                      'k', UVAR_REQUIRED, "Location of ephemeris files");
  LALregSTRINGUserVar(status, 	ephemyear, 	                      'l', UVAR_REQUIRED, "Ephemeris file year");
  LALregSTRINGUserVar(status, 	ifo, 	                              'm', UVAR_REQUIRED, "Interferometer name (G1,L1,H1,H2)");
  LALregSTRINGUserVar(status, 	output, 	                      'n', UVAR_REQUIRED, "Name of output file");
  LALregSTRINGUserVar(status, 	fourierfile, 	                      'o', UVAR_REQUIRED, "Name of input fstat file");
  LALregSTRINGUserVar(status, 	rangefile, 	                      'y', UVAR_REQUIRED, "Name of parameter range file");
  LALregSTRINGUserVar(status, 	windowfile, 	                      'J', UVAR_REQUIRED, "Name of output window file");
  LALregINTUserVar(status, 	seed,                                 'v', UVAR_OPTIONAL, "The seed for random number generation (0 uses clock)");
  LALregINTUserVar(status, 	Nchain,                               'w', UVAR_REQUIRED, "The number of links in the chain");
  LALregREALUserVar(status, 	sqrtSh, 	                      'x', UVAR_REQUIRED, "The sqrt of the noise spectral density (in Hz^{-1/2})"); 
  LALregINTUserVar(status, 	Nburn,                                'C', UVAR_REQUIRED, "The number of steps in the burn in stage");
  LALregSTRINGUserVar(status, 	stampsfile, 	                      'D', UVAR_OPTIONAL, "Name of timestamps file");
  LALregREALUserVar(status, 	highresfactor, 	                      'E', UVAR_OPTIONAL, "The over-resolution of the window function (default 100)");
  LALregINTUserVar(status, 	tsft,                                 'F', UVAR_REQUIRED, "The length of the SFTs used (default 1800)");
  LALregBOOLUserVar(status, 	noneccentric,                         'G', UVAR_OPTIONAL, "Set if using a non-eccentric orbit (e=0) (default 0)");
  LALregBOOLUserVar(status, 	ampmod,                               'I', UVAR_OPTIONAL, "Set if you wish to use amplitude modulation sidebands");
  
  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */


/*----------------------------------------------------------------------*/
/**
 * Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *status)
{

  INITSTATUS(status);  
  /* don't allow unspecified output file */
  if (uvar_output == NULL) {
    XLALPrintError ("\nOutput file must be specified (option 'output')\n\n");
    ABORT (status, SIDEBANDMCMCC_EINPUT, SIDEBANDMCMCC_MSGEINPUT);
  }
  /* don't allow right ascension out of range */
  if ((uvar_alpha < 0)||(uvar_alpha >= LAL_TWOPI)) {
    XLALPrintError ("\nRight ascension must be in range [0,2PI)!\n\n");
    ABORT (status, SIDEBANDMCMCC_EINPUT, SIDEBANDMCMCC_MSGEINPUT);
  }
  /* don't allow declination out of range */
  if ((uvar_delta < ((-1.0)*LAL_PI/2.0))||(uvar_delta > (LAL_PI/2.0))) {
    XLALPrintError ("\nDeclination must be in range [-PI/2,PI/2]!\n\n");
    ABORT (status, SIDEBANDMCMCC_EINPUT, SIDEBANDMCMCC_MSGEINPUT);
  }
    
 
  RETURN (status);
} /* checkUserInputConsistency() */


/***********************************************************************************/
/*----------------------------------------------------------------------*/
/* Simply initialises the MCMC parameter vector with some randomly chosen parameter values
 */
void 
InitialiseLambda(LALStatus *status,SideBandMCMCVector *lambda,SideBandMCMCRanges ranges,RandomParams *randparams)
{
  
  REAL4Vector *vector = NULL;          /* stores vector of random parameters */
  INT4 i;                              /* general loop index */
  REAL8 temp;                          /* temporary time used for setting tp */
  int compareGPS;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (lambda,status,SIDEBANDMCMCC_ENULL,SIDEBANDMCMCC_MSGENULL );

  /* generate random numbers */
  vector = XLALCreateREAL4Vector(NPARAMS+1);
  for (i=0;i<(INT4)vector->length;i++) {
    LALUniformDeviate(status->statusPtr,vector->data+i,randparams);
  }

  /* initialise the parameters */
  lambda->f0 = ranges.f0min + (REAL8)vector->data[0]*(ranges.f0max - ranges.f0min);
  lambda->a = ranges.amin + (REAL8)vector->data[1]*(ranges.amax - ranges.amin);
  compareGPS = XLALGPSCmp(&(ranges.tpmin),&(ranges.tpmax));
  if (compareGPS!=0) {
    temp = ranges.tpmin.gpsSeconds+1e-9*ranges.tpmin.gpsNanoSeconds + 
      (REAL8)vector->data[2]*(ranges.tpmax.gpsSeconds+1e-9*ranges.tpmax.gpsNanoSeconds 
			      - ranges.tpmin.gpsSeconds-1e-9*ranges.tpmin.gpsNanoSeconds);
    XLALGPSSetREAL8(&(lambda->tp),temp);
  }
  else {
    lambda->tp.gpsSeconds = ranges.tpmin.gpsSeconds;
    lambda->tp.gpsNanoSeconds = ranges.tpmin.gpsNanoSeconds;
  }
  lambda->argp = ranges.argpmin + (REAL8)vector->data[3]*(ranges.argpmax - ranges.argpmin);
  lambda->e = ranges.emin + (REAL8)vector->data[4]*(ranges.emax - ranges.emin);
  lambda->h0 = ranges.h0min + (REAL8)vector->data[5]*(ranges.h0max - ranges.h0min);
  lambda->period = ranges.periodmin + (REAL8)vector->data[6]*(ranges.periodmax - ranges.periodmin);
  lambda->cosi = ranges.cosimin + (REAL8)vector->data[7]*(ranges.cosimax - ranges.cosimin);
  lambda->psi = ranges.psimin + (REAL8)vector->data[8]*(ranges.psimax - ranges.psimin);
  lambda->phi0 = ranges.phi0min + (REAL8)vector->data[9]*(ranges.phi0max - ranges.phi0min);
  
  /* free memory */
  XLALDestroyREAL4Vector(vector);
  
  DETATCHSTATUSPTR (status);
  RETURN(status);
  
}
/***********************************************************************************/
/*----------------------------------------------------------------------*/
/* Makes an N-Dimensional random jump in parameter space
 */
void
MakeJump(LALStatus *status,SideBandMCMCVector lambda,SideBandMCMCVector *newlambda,SideBandMCMCJumpProbs jumpsizes,SideBandMCMCRanges ranges,RandomParams *randparams)
{
  
  REAL4Vector *rndm = NULL;           /* vector for storing random numbers */
  REAL8 temp;                           /* used to store a temporary random time jump */
  REAL8 temptp;
  INT4 tpflag = 1;
  REAL4Vector *randomscale = NULL;
  REAL8 x,newx;
  REAL8 y,newy;
  INT4 i;
  REAL8 f0jump,ajump,ejump,xjump,yjump,h0jump,psijump,cosijump,phi0jump,periodjump;
  int compare1,compare2,compare3;
  double delta1;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* ASSERT (newlambda,status,SIDEBANDMCMCC_ENULL,SIDEBANDMCMCC_MSGENULL ); */

  /* generate random numbers for actual jumps */
  rndm = XLALCreateREAL4Vector(NPARAMS+1);
  LALNormalDeviates(status->statusPtr,rndm,randparams);
  
  /* generate random numbers to determine jump scale */
  randomscale = XLALCreateREAL4Vector(NPARAMS+1);
  for (i=0;i<NPARAMS;i++) {
    LALUniformDeviate(status->statusPtr,&(randomscale->data[i]),randparams);
    /* printf("randomscale = %f\n",randomscale->data[i]);*/
  }

  /* select jump sizes */
  if (randomscale->data[0]<jumpsizes.f0.prob[0]) f0jump = jumpsizes.f0.jump[0];
  else if (randomscale->data[0]<jumpsizes.f0.prob[0]+jumpsizes.f0.prob[1]) {
    f0jump = jumpsizes.f0.jump[1];
    /* random->data[0] = random->data[0]/fabs(random->data[0]); */  /* try an exact jump of medium size */
  }
  else f0jump = jumpsizes.f0.jump[2];
  /* printf("fijump = %9.12f\n",fijump); */

  /* if (randomscale->data[1]<jumpsizes.a.prob[0]) ajump = jumpsizes.a.jump[0];
  else if (randomscale->data[1]<jumpsizes.a.prob[0]+jumpsizes.a.prob[1]) ajump = jumpsizes.a.jump[1];
  else ajump = jumpsizes.a.jump[2]; */
  ajump = jumpsizes.a.jump[0];
  /* printf("ajump = %9.12f\n",ajump); */

  if (randomscale->data[2]<jumpsizes.x.prob[0]) xjump = jumpsizes.x.jump[0];
  else if (randomscale->data[2]<jumpsizes.x.prob[0]+jumpsizes.x.prob[1]) xjump = jumpsizes.x.jump[1];
  else xjump = LAL_PI/(2.0*lambda.a*lambda.f0);  /* approximate ridge seperation in tp-argp plane */
  /* printf("xjump = %9.12f\n",xjump); */

  /* if (randomscale->data[3]<jumpsizes.y.prob[0]) yjump = jumpsizes.y.jump[0];
  else if (randomscale->data[3]<jumpsizes.y.prob[0]+jumpsizes.y.prob[1]) yjump = jumpsizes.y.jump[1];
  else yjump = jumpsizes.y.jump[2]; */
  yjump = jumpsizes.y.jump[0];
  /* printf("yjump = %9.12f\n",yjump); */

  /* if (randomscale->data[4]<jumpsizes.e.prob[0]) ejump = jumpsizes.e.jump[0];
  else if (randomscale->data[4]<jumpsizes.e.prob[0]+jumpsizes.e.prob[1]) ejump = jumpsizes.e.jump[1];
  else ejump = jumpsizes.e.jump[2]; */
  ejump = jumpsizes.e.jump[0];
  /* printf("ejump = %9.12f\n",ejump); */

  /* if (randomscale->data[5]<jumpsizes.h0.prob[0]) h0jump = jumpsizes.h0.jump[0];
  else if (randomscale->data[5]<jumpsizes.h0.prob[0]+jumpsizes.h0.prob[1]) h0jump = jumpsizes.h0.jump[1];
  else h0jump = jumpsizes.h0.jump[2]; */
  h0jump = jumpsizes.h0.jump[0];

  periodjump = jumpsizes.period.jump[0];

  /* if (randomscale->data[6]<jumpsizes.psi.prob[0]) psijump = jumpsizes.psi.jump[0];
  else if (randomscale->data[6]<jumpsizes.psi.prob[0]+jumpsizes.psi.prob[1]) psijump = jumpsizes.psi.jump[1];
  else psijump = jumpsizes.psi.jump[2]; */
  psijump = jumpsizes.psi.jump[0];

  /* if (randomscale->data[7]<jumpsizes.cosi.prob[0]) cosijump = jumpsizes.cosi.jump[0];
  else if (randomscale->data[7]<jumpsizes.cosi.prob[0]+jumpsizes.cosi.prob[1]) cosijump = jumpsizes.cosi.jump[1];
  else cosijump = jumpsizes.cosi.jump[2]; */
  cosijump = jumpsizes.cosi.jump[0];

  /* if (randomscale->data[8]<jumpsizes.phi0.prob[0]) phi0jump = jumpsizes.phi0.jump[0];
  else if (randomscale->data[8]<jumpsizes.phi0.prob[0]+jumpsizes.phi0.prob[1]) phi0jump = jumpsizes.phi0.jump[1];
  else phi0jump = jumpsizes.phi0.jump[2]; */
  phi0jump = jumpsizes.phi0.jump[0];

  /* decide whether to take a large scale jump in x */
  if (!uvar_noneccentric) {

    /********************************************************************/
    /* reparameterise tp and argp as x and y */
    /* tp = 0.5*(x + y) and argp = 0.5*(x - y) */
    /* x = tp + argp and y = tp - argp */
    x = (LAL_TWOPI/lambda.period)*(lambda.tp.gpsSeconds+1e-9*lambda.tp.gpsNanoSeconds) - lambda.argp;
    y = (LAL_TWOPI/lambda.period)*(lambda.tp.gpsSeconds+1e-9*lambda.tp.gpsNanoSeconds) + lambda.argp;
    /********************************************************************/
    /* printf("x = %6.12f y = %6.12f\n",x,y); */

    /* make the x,y jumps */
    /* printf("xjump = %6.12f yjump = %6.12f\n",xjump,yjump); */
    newx = x + (REAL8)rndm->data[1]*xjump;
    newy = y + (REAL8)rndm->data[2]*yjump; 
    /* printf("newx = %6.12f newy = %6.12f\n",newx,newy); */

    /* printf("old tp = %6.12f old argp = %6.12f\n",lambda.tp.gpsSeconds+1e-9*lambda.tp.gpsNanoSeconds,lambda.argp); */
    /* convert back to tp and argp */
    temptp = (0.5*lambda.period/LAL_TWOPI)*(newx + newy);
    XLALGPSSetREAL8(&(newlambda->tp),temptp);
    newlambda->argp = 0.5*(newy - newx); 
    /* printf("tp = %6.12f argp = %6.12f\n",temptp,newlambda->argp);*/
    
  }
  else {

    /* make the tp,argp jumps */
    temp = (REAL8)rndm->data[1]*xjump;
    newlambda->tp = lambda.tp;
    XLALGPSAdd(&(newlambda->tp), temp);
    newlambda->argp = lambda.argp + (REAL8)rndm->data[2]*yjump;
  }  


  /* printf("jumping fi %f\n",(REAL8)rndm->data[0]*fijump); */
  newlambda->f0 = lambda.f0 + (REAL8)rndm->data[0]*f0jump;
  newlambda->a = lambda.a + (REAL8)rndm->data[3]*ajump;
  newlambda->e = lambda.e + (REAL8)rndm->data[4]*ejump;
  newlambda->h0 = lambda.h0 + (REAL8)rndm->data[5]*h0jump;
  newlambda->period = lambda.period + (REAL8)rndm->data[6]*periodjump;
  newlambda->cosi = lambda.cosi + (REAL8)rndm->data[7]*cosijump;
  newlambda->psi = lambda.psi + (REAL8)rndm->data[8]*psijump;
  newlambda->phi0 = lambda.phi0 + (REAL8)rndm->data[9]*phi0jump;

  /* make sure we are within prior boundaries */

  /* wrap time of periapse around range */
  /* make comparisons between boundary values */
  compare3 = XLALGPSCmp(&(ranges.tpmin),&(ranges.tpmax));

  /* if boundary values not equal to each other */
  if (compare3!=0) {
    
    while (tpflag) {

      /* make comparisons between new tp and boundary values */
      compare1 = XLALGPSCmp(&(newlambda->tp),&(ranges.tpmin));
      compare2 = XLALGPSCmp(&(newlambda->tp),&(ranges.tpmax));

      /* if new tp is less than min boundary */
      if (compare1==-1) {
	/* printf("newlambda->tp = %d %d tpmin = %d %d tpmax = %d %d\n",newlambda->tp.gpsSeconds,newlambda->tp.gpsNanoSeconds,ranges.tpmin.gpsSeconds,ranges.tpmin.gpsNanoSeconds,ranges.tpmax.gpsSeconds,ranges.tpmax.gpsNanoSeconds);
	   printf(" new tp < min\n"); */
	/* find difference between new tp and min boundary */
        delta1 = XLALGPSDiff(&(ranges.tpmin),&(newlambda->tp));
	/* printf("interval = %d %d\n",delta1.seconds,delta1.nanoSeconds); */
	/* minus interval from max boundary (max boundary must be one orbital period from min boundary) */
        newlambda->tp = ranges.tpmax;
        XLALGPSAdd(&(newlambda->tp), -delta1);
	/* printf("changed newlambda->tp to %d %d\n",newlambda->tp.gpsSeconds,newlambda->tp.gpsNanoSeconds); */
      }
      else if (compare2==1) {
	/* printf("newlambda->tp = %d %d tpmin = %d %d tpmax = %d %d\n",newlambda->tp.gpsSeconds,newlambda->tp.gpsNanoSeconds,ranges.tpmin.gpsSeconds,ranges.tpmin.gpsNanoSeconds,ranges.tpmax.gpsSeconds,ranges.tpmax.gpsNanoSeconds);
	   printf(" new tp > max\n"); */
	/* find difference between new tp and max boundary */
        delta1 = XLALGPSDiff(&(newlambda->tp),&(ranges.tpmax));
	/* printf("interval = %d %d\n",delta1.seconds,delta1.nanoSeconds); */
	/* add interval to the min boundary (max boundary must be one orbital period from min boundary) */
        newlambda->tp = ranges.tpmin;
	XLALGPSAdd(&(newlambda->tp), delta1);
	/* printf("changed newlambda->tp to %d %d\n",newlambda->tp.gpsSeconds,newlambda->tp.gpsNanoSeconds); */
      }
      else tpflag = 0;
      
    }
   
  }
  else {
    newlambda->tp.gpsSeconds = lambda.tp.gpsSeconds;
    newlambda->tp.gpsNanoSeconds = lambda.tp.gpsNanoSeconds;
  }

  /* wrap argument of periapse around range */
  if (ranges.argpmax!=ranges.argpmin) {  
    if (newlambda->argp<ranges.argpmin) {
      temp = fmod((ranges.argpmin-newlambda->argp),(ranges.argpmax-ranges.argpmin));
      newlambda->argp = ranges.argpmax - temp;
    }
    else if (newlambda->argp>ranges.argpmax) {
      temp = fmod((newlambda->argp-ranges.argpmax),(ranges.argpmax-ranges.argpmin));
      newlambda->argp = ranges.argpmin + temp;
    }
  }
  else newlambda->argp = lambda.argp;

  /* bounce cos(i) off the boundaries */
  if (ranges.cosimax!=ranges.cosimin) {  
    if (newlambda->cosi<ranges.cosimin) {
      temp = fmod((ranges.cosimin-newlambda->cosi),(ranges.cosimax-ranges.cosimin));
      newlambda->cosi = ranges.cosimin +  temp;
    }
    else if (newlambda->cosi>ranges.cosimax) {
      temp = fmod((newlambda->cosi-ranges.cosimax),(ranges.cosimax-ranges.cosimin));
      newlambda->cosi = ranges.cosimax - temp;
    }
  }
  else newlambda->cosi = lambda.cosi; 

  /* wrap phi0 around range */
  if (newlambda->phi0<ranges.phi0min) {
    temp = fmod((ranges.phi0min-newlambda->phi0),(ranges.phi0max-ranges.phi0min));
    newlambda->phi0 = ranges.phi0max - temp;
  }
  else if (newlambda->phi0>ranges.phi0max) {
    temp = fmod((newlambda->phi0-ranges.phi0max),(ranges.phi0max-ranges.phi0min));
    newlambda->phi0 = ranges.phi0min + temp;
  }

  /* wrap psi around range and flip phi0 */
  if (newlambda->psi<ranges.psimin) {
    temp = fmod((ranges.psimin-newlambda->psi),(ranges.psimax-ranges.psimin));
    newlambda->psi = ranges.psimax - temp;
    newlambda->phi0 = fmod(newlambda->phi0+LAL_PI,LAL_TWOPI);
  }
  else if (newlambda->psi>ranges.psimax) {
    temp = fmod((newlambda->psi-ranges.psimax),(ranges.psimax-ranges.psimin));
    newlambda->psi = ranges.psimin + temp;
    newlambda->phi0 = fmod(newlambda->phi0+LAL_PI,LAL_TWOPI);
  }

  /* we are going to try bouncing eccentricity off zero */
  if (ranges.emax!=ranges.emin) {  
    if (newlambda->e<ranges.emin) {
      temp = fmod((ranges.emin-newlambda->e),(ranges.emax-ranges.emin));
      newlambda->e = ranges.emin +  temp;
    }
  }
  else newlambda->e = lambda.e; 
 
  /* free memory */
  XLALDestroyREAL4Vector(rndm);

  DETATCHSTATUSPTR (status);
  RETURN(status);
  
  
}


/***********************************************************************************/
/*----------------------------------------------------------------------*/
/* Tests the new jump likelihood and decides whether to jump there or not 
 */
void
TestJump(LALStatus *status,SideBandMCMCVector *lambda,SideBandMCMCVector *newlambda,SideBandMCMCVector *chain,RandomParams *randparams,INT4 n,INT4 out,INT4 *flag)
{
  
  REAL4 rndm;                        /* stores a single random number */
  REAL8 R;                             /* the ratio of likelihoods */
  REAL8 T = 1e6;                       /* starting temperature */
  REAL8 logrand;

  INITSTATUS(status); 
  ATTATCHSTATUSPTR (status); 
  
  ASSERT (lambda,status,SIDEBANDMCMCC_ENULL,SIDEBANDMCMCC_MSGENULL );
  ASSERT (newlambda,status,SIDEBANDMCMCC_ENULL,SIDEBANDMCMCC_MSGENULL ); 
  /* ASSERT (chain,status,SIDEBANDMCMCC_ENULL,SIDEBANDMCMCC_MSGENULL ); */
  
  /* if we have not jumped out of prior range */
  if (out==0) {

    /* generate random numbers */
    LALUniformDeviate(status->statusPtr,&rndm,randparams);
    
    /* printf("comparing (old) %f - (new) %f\n",lambda->logL,newlambda->logL); */
    
    if (n<=uvar_Nburn) R = log(T)*(1.0 - ((REAL8)n/(REAL8)uvar_Nburn)) + newlambda->logL - lambda->logL;
    else R = newlambda->logL - lambda->logL;
    /* R = newlambda->L/lambda->L;  */
    logrand = log((REAL8)rndm);

    /* printf("R = %f logrand = %f\n",R,logrand); 
       printf("rndm = %f\n",(REAL8)rndm); */
    

    /* if new likelihood is greater than old likelihood then defineitely jump there */
    if (R>logrand) {
      *lambda = *newlambda;
      lambda->accepted = 1;
      *flag = 1;
      /* printf("jumping\n"); */
    }
    else {
      lambda->accepted = 0;
      *flag = 0;
      /* printf("not jumping\n"); */
    }

    /* if (newlambda->forcejump) printf("FORCING JUMP!!!\n"); */
    
    
    /* testing FORCING JUMP inside range */
    /* if (newlambda->logL > -1e32) {
     *lambda=*newlambda;
     *flag = 1; 
     } */

  } 

  /* add to the chain */
  *chain = *lambda;
  /* chain->a = lambda->a;
  chain->tp.gpsSeconds = lambda->tp.gpsSeconds;
  chain->tp.gpsNanoSeconds = lambda->tp.gpsNanoSeconds;
  chain->h0 = lambda->h0;
  chain->cosi = lambda->cosi;
  chain->psi = lambda->psi;
  chain->phi0 = lambda->phi0; */

  /* printf("added to the chain\n"); */
  
 
  DETATCHSTATUSPTR (status); 
  RETURN(status);

}

/***********************************************************************************/
/*----------------------------------------------------------------------*/
/* determines whether we are inside prior ranges */
void
PriorRanges(LALStatus *status,SideBandMCMCVector lambda,SideBandMCMCRanges ranges, INT4 *out)
{

  int compare1, compare2;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  compare1 = XLALGPSCmp(&(lambda.tp),&(ranges.tpmin));
  compare2 = XLALGPSCmp(&(lambda.tp),&(ranges.tpmax));


  *out = 0;
  if ((lambda.f0<ranges.f0min)||(lambda.f0>ranges.f0max)) {
    /* printf("jumped out in findex\n"); */ 
    *out = 1;
  }
  else if ((lambda.h0<ranges.h0min)||(lambda.h0>ranges.h0max)) {
    /* printf("jumped out in h0\n"); */
    *out = 1;   
  }
  else if ((compare1==-1)||(compare2==1)) {
    /* printf("jumped out in tp\n"); */
    *out = 1;  
  }
  else if ((lambda.period<ranges.periodmin)||(lambda.period>ranges.periodmax)) {
    /* printf("jumped out in argp\n");  */
    *out = 1;   
  }
  else if ((lambda.argp<ranges.argpmin)||(lambda.argp>ranges.argpmax)) {
    /* printf("jumped out in argp\n"); */
    *out = 1;   
  }
  else if ((lambda.e<ranges.emin)||(lambda.e>ranges.emax)) {
    /* printf("jumped out in e\n"); */
    *out = 1;   
  }
  else if ((lambda.psi<ranges.psimin)||(lambda.psi>ranges.psimax)) {
    /* printf("jumped out in psi\n"); */
    *out = 1;   
  }
  else if ((lambda.cosi<ranges.cosimin)||(lambda.cosi>ranges.cosimax)) {
    /* printf("jumped out in cosi\n"); */
    *out = 1;   
  }
  else if ((lambda.phi0<ranges.phi0min)||(lambda.phi0>ranges.phi0max)) {
    /* printf("jumped out in phi0\n"); */
    *out = 1;   
  }
  else if ((lambda.a<ranges.amin)||(lambda.a>ranges.amax)) {
    /* printf("jumped out in a\n");  */
    *out = 1;  
  }
  

  DETATCHSTATUSPTR (status);
  RETURN(status);

 
}


