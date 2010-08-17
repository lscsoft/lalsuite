/*
 * Copyright (C) 2009 Letizia Sammut
 * Copyright (C) 2009 Chris Messenger
 * Copyright (C) 2007 Reinhard Prix
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

/** \author L.Sammut, C. Messenger
 * \file
 * \brief
 * Calculates the C-statistic for a given parameter-space of GW signals from binary sources with known sky position.
 * 
 * Uses outputFStat file of lalapps_ComputeFStatistic_v2.c as input.
 *
 */

/***********************************************************************************************/
#include "config.h"

/* System includes */
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <strings.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* LAL-includes */
#include <lal/LALConfig.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdio.h>
#include <lal/LALError.h>
#include <lal/LALStdlib.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>

#include <lal/AVFactories.h>
#include <lal/GSLSupport.h>

#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>

#include <lalapps.h>


RCSID( "$Id$");

/***********************************************************************************************/
/* defines */

#define MAXFILENAMELENGTH 256   	/* Maximum # of characters of a SFT filename */
#define MAXLINELENGTH 1024   		/* Maximum # of characters of in a line */
#define TRUE (1==1)
#define FALSE (1==0)

/***********************************************************************************************/
/* internal structures */

/* Structure for parameters to be input to Cstat */
typedef struct {
  REAL8 Alpha;				/**< sky position alpha in radians read in from data file 	*/
  REAL8 Delta;                          /**< sky position delta in radians read in from data file 	*/
  REAL8 *Fstat;				/**< values of Fstatistic read in from data file 		*/
  REAL8 fmin;				/**< min frequency read in from data file 			*/
  REAL8 fmax;				/**< max frequency read in from data file 			*/
  REAL8 df;				/**< size of frequency bin read in from data file 		*/
   
  REAL8 f0;				/**< guess frequency of search 					*/
  INT4 dfout;				/**< number of bins to exclude from the start of the Fstat file */
  INT4 dm;			        /**< number of bins excluded from Cstat output 			*/

  INT4 N;				/**< number of frequency bins to search in ComputeFStatistic 	*/
} ConfigVariables;


/* Structure for output of computeCStat */
typedef struct {
  INT4 sideband;			/**< number of sidebands 	*/
  REAL8 *Cstat; 		        /**< cstat array 		*/
} CstatOut;


/* User-variables: can be set from config-file or command-line */
typedef struct {
 
  REAL8 Freq;				/**< start user frequency band for output 	*/
  REAL8 FreqBand;			/**< user Frequency-band for output 		*/
  
  /* orbital parameters */
  REAL8 orbitPeriod;			/**< binary-system orbital period in s 		*/
  REAL8 orbitasini;			/**< amplitude of radial motion 		*/
 
  CHAR *inputFstat;			/**< filename of Fstat input data file to use 	*/

  BOOLEAN help;				/**< output help-string 			*/
  CHAR *outputCstat;			/**< filename to output Cstatistic in 		*/

  BOOLEAN version;			/**< output version information 		*/
} UserInput_t;

/***********************************************************************************************/
/* Global variables */
extern int vrbflg;			/**< defined in lalapps.c */

/* empty initializers */
static const ConfigVariables empty_ConfigVariables;
static UserInput_t empty_UserInput;
static CstatOut empty_CstatOut;

/* local prototypes */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *, UserInput_t *uvar);
int checkUserInputConsistency (const UserInput_t *uvar );
void OutputVersion ( void );
int getFStat(CHAR *filename, ConfigVariables *cfg);
int computeCStat(ConfigVariables *cfg, UserInput_t *uvar, CstatOut *cst);
void Freemem(LALStatus *,  ConfigVariables *cfg, CstatOut *cst);

const char *va(const char *format, ...);	/* little var-arg string helper function */

/***********************************************************************************************/
/* Function definitions */

/*--------------------------------------------------------------- */
/** getFstat function 
 * reads in Fstat file input and assigns parameter values 
 */
/*----------------------------------------------------------------*/
int getFStat(CHAR *filename, ConfigVariables *cfg)
{
  static const char *fn = __func__;             /* store function name for log output */

  CHAR line[MAXLINELENGTH];
  REAL8 dummy, *a, *d, *frequency;
  FILE *data = NULL;
  INT4 c=0, l=0;
   
  /* Open data file - check it is good */
  if ((data = fopen(filename, "r")) == NULL)	{
    LogPrintf (LOG_CRITICAL, "%s: Error opening file '%s' for reading.. Error %d\n",fn,filename,xlalErrno);
    return XLAL_EIO;
  }

  /* need to know number of data points - counting the number of data lines (ignoring comments) */
  while(fgets(line, MAXLINELENGTH, data) != NULL)	{
    if (strncmp(&line[0], "%",1) != 0 && isdigit(line[0]) != 0) {
      l++;
    }
  }
 
  fclose(data); /* close the file prior to exiting the routine */
   
  /* Check data has more than one Fstat entry */
  if (l==1)	{
    LogPrintf (LOG_CRITICAL, "%s: Must be more than one Fstat entry in data file %s. Error %d\n",fn,filename,xlalErrno);
    return XLAL_EDOM;
  }
  
  /* Allocate some memory according to number of frequency bins */  
  if ((frequency=(REAL8*)LALMalloc(l*sizeof(REAL8)))== NULL )	{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory.Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }
  if ((a=(REAL8*)LALMalloc(l*sizeof(REAL8)))== NULL )		{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }
  if ((d=(REAL8*)LALMalloc(l*sizeof(REAL8)))== NULL )		{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM; 
  }
  if ((cfg->Fstat=(REAL8*)LALMalloc(l*sizeof(REAL8)))== NULL )	{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }
 
  /* Open data file - check it is good */
  if ((data = fopen(filename, "r")) == NULL)	{
    LogPrintf (LOG_CRITICAL, "%s: Error opening file '%s' for reading.. Error %d\n",fn,filename,xlalErrno);
    return XLAL_EIO;
  }

  /* Get data, assign configuration variables */  
  while (fgets(line, MAXLINELENGTH, data) != NULL)	{
    if (strncmp(line,"%%",1) != 0 && isdigit(line[0]) != 0 ) 	{
      sscanf (line, "%lf %lf %lf %lf %lf %lf %lf", &frequency[c], &a[c], &d[c], &dummy, &dummy ,&dummy, &(cfg->Fstat[c]) );
      c++;
    }  
  }	

  /* Assign config variables alpha, delta, fmin, fmax and df */
  cfg-> Alpha 	= a[0];
  cfg-> Delta	= d[0];
  cfg-> fmin	= frequency[0];
  cfg-> fmax	= frequency[c-1];
  cfg-> df	= frequency[1]-frequency[0];
  
  /* check for constant frequency seperation and constant sky position */
  for (c=1; c<l; c++) {
    if (((frequency[c]-frequency[c-1]) - cfg->df) > 1e-8) 	{
      LogPrintf (LOG_CRITICAL, "%s: Frequency separation df in Fstat file '%s' must be constant.\n",fn,filename,xlalErrno);
      return XLAL_EINVAL;
    }
    if ((a[c]-a[0]) != 0) 	{
      LogPrintf (LOG_CRITICAL, "%s: Sky position (alpha) must be constant. Error %d\n",fn,xlalErrno);
      return XLAL_EINVAL;
    }
    if ((d[c]-d[0]) != 0) 	{
      LogPrintf (LOG_CRITICAL, "%s: Sky position (delta) must be constant. Error %d\n",fn,xlalErrno);
      return XLAL_EINVAL;
    }  
  }
 
  /* Free getFStat function specific memory */
  LALFree(a);
  LALFree(d);
  LALFree(frequency);
       
  return XLAL_SUCCESS;

} /* getFStat*/



/*--------------------------------------------------------------- */
/** computCStat function 
 * receives input of Fstat and outputs Cstat 
 */
/*----------------------------------------------------------------*/
int computeCStat(ConfigVariables *cfg, UserInput_t *uvar, CstatOut *cst)
{ 
  static const char *fn = __func__;             /* store function name for log output */ 

  REAL8 ufreq, fstart, fend;
  INT4 ufband, ufreq_bin, mm, i, ind=0;
  		
  /* Allocate memory for FFT plans and vectors */
  REAL8FFTPlan *pfwd = NULL;			
  REAL8FFTPlan *prev = NULL;
  REAL8Vector *fstat = NULL;
  REAL8Vector *comb = NULL;
  COMPLEX16Vector *fout = NULL;
  COMPLEX16Vector *cout = NULL;
  COMPLEX16Vector *out = NULL;
  REAL8Vector *cstat = NULL;

  /* create comb search parameters */
  cfg->f0 	= uvar->Freq + 0.5*uvar->FreqBand ; 			/* allocate guess frequency as centre of user input search frequency band */ 
  ufreq_bin	= floor(0.5 + (uvar->Freq-cfg->fmin)/cfg->df);		/* nearest bin to user starting frequency */
  ufreq		= cfg->fmin + (cfg->df*ufreq_bin);			/* frequency of Fstat bin nearest to user frequency */
  ufband	= floor(0.5 + (uvar->FreqBand/cfg->df));		/* user frequency band in bins */
  mm		= floor(0.5 + LAL_TWOPI*cfg->f0*uvar->orbitasini);	/* whole number of sidebands on either side of central spike */
  cfg->dm	= floor(0.5 + mm/(uvar->orbitPeriod*cfg->df));		/* exclusion region: number of bins to exclude after Cstat calculation */
  cfg->N	= ufband + 2*cfg->dm;					/* number of bins for cstat memory allocation */
  fstart	= ufreq - cfg->df*cfg->dm;				/* start frequency of search, half a comb width before user specified Freq */
  fend		= ufreq + cfg->df*cfg->N;				/* end search frequency band half a comb width after user specified Freq */
  cfg->dfout	= (fstart - cfg->fmin)/cfg->df;				/* number of bins to exclude from calculation of cstat */
   
  /* check search band plus exclusion region is within data frequency band */
  if ( (fstart < cfg->fmin) || (cfg->fmax < fend) )	{
      LogPrintf (LOG_CRITICAL, "%s: User input frequency range and/or exclusion range outside data limits. Error %d\n",fn,xlalErrno);
      return XLAL_EDOM;
  } 
         
  /* Create FFT plans and vectors */
  pfwd=XLALCreateREAL8FFTPlan(cfg->N, 1, 0);
  prev=XLALCreateREAL8FFTPlan(cfg->N, 0, 0);
  fstat=XLALCreateREAL8Vector(cfg->N );
  comb=XLALCreateREAL8Vector(cfg->N );
  fout=XLALCreateCOMPLEX16Vector(cfg->N/2 +1 );
  cout=XLALCreateCOMPLEX16Vector(cfg->N/2 +1 );
  out=XLALCreateCOMPLEX16Vector(cfg->N/2 +1 );
  cstat=XLALCreateREAL8Vector(cfg->N );  

  /* Assign fstat arrary values of Fstatistic from input file */
  for (i=0;i<(cfg->N);i++) {
    fstat->data[i] =cfg->Fstat[i+cfg->dfout];
  }
  
  /* Fourier transform fstat array */
  XLALREAL8ForwardFFT(fout, fstat, pfwd );
  
  /* Create comb template */
  /*  - zero all values first */
  for (i=0;i<cfg->N;i++) { 
    comb->data[i]=0;						
  }
  
  /*  - assign unity at spacings of 1/P for 1 zero spike + mm positive frequency spikes */
  for (i=0;i<=mm;i++)	{
    ind=floor(0.5+ i/(uvar->orbitPeriod*cfg->df));
    comb->data[ind]=1;						
  }
  
  /*  - assign unity at spacings of 1/P for mm negative frequency spikes */
  for (i=1;i<=mm;i++)	{
    ind=(cfg->N) -floor(0.5 +i/(uvar->orbitPeriod*cfg->df));
    comb->data[ind]=1;						
  }
  
  /*  - Fourier transform template for convolution with fstat Fourier transform */
  XLALREAL8ForwardFFT(cout, comb, pfwd );			


  /* Perform convolution of fstat with template by multiplication in Fourier time domain */
  for (i=0;i<(cfg->N/2 +1); i++)	{
    out->data[i].re = (fout->data[i].re * cout->data[i].re) - (fout->data[i].im * cout->data[i].im); /* real part of out */
    out->data[i].im = (fout->data[i].re * cout->data[i].im) + (fout->data[i].im * cout->data[i].re); /* imaginary part of out */
   }

  /* Inverse FFT back to frequency domain to retrieve Cstat */
  XLALREAL8ReverseFFT(cstat, out, prev );
  
  /* Fill out CstatOut struct */
  /*  - Allocate some memory */
  if ((cst->Cstat=(REAL8*)LALMalloc(cfg->N*sizeof(REAL8)))== NULL )	{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM; 
  }
  
  /*  - Allocate total number of sidebands */
  cst->sideband = 2*mm+1;
  							
  /*  - Assign Cstat arrary values of Cstatistic */
  for (i=0; i<cfg->N; i++) {
    cst->Cstat[i] = cstat->data[i];						
  }
   
  /* Destroy FFT plans and free computeCStat function specific memory */
  XLALDestroyREAL8FFTPlan( pfwd );
  XLALDestroyREAL8FFTPlan( prev );
  XLALDestroyREAL8Vector(fstat );
  XLALDestroyREAL8Vector(comb );
  XLALDestroyCOMPLEX16Vector( fout);
  XLALDestroyCOMPLEX16Vector( cout);
  XLALDestroyCOMPLEX16Vector( out);
  XLALDestroyREAL8Vector(cstat );

  return XLAL_SUCCESS;

} /* computeCStat */



/*----------------------------------------------------------------------*/ 
/**
 * MAIN function of sb_search code.
 * Calculate the C-statistic over a given portion of the parameter-space
 * and write output into a file(default: 'Cstats').
 */
/*----------------------------------------------------------------------*/ 
int main(int argc,char *argv[])
{ 
  static const char *fn = __func__;             	/* store function name for log output */
  LALStatus status = blank_status;			/* initialize status */
  UserInput_t uvar = empty_UserInput;			/* global container for user variables */
  ConfigVariables cfg = empty_ConfigVariables;		/* global container for various derived configuration settings */
  CstatOut cst = empty_CstatOut;			/* global container for cstat struct */
  FILE *Cstat_out = NULL;
  INT4 i;
  CHAR *cmdline = NULL; 
  

  lalDebugLevel = 0;
  vrbflg = 1;						/* verbose error-messages */


  /* register all user-variable */
  if (XLALGetDebugLevel(argc, argv, 'v')) {
    LogPrintf(LOG_CRITICAL,"%s : XLALGetDebugLevel() failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }

  LAL_CALL (initUserVars(&status, &uvar), &status); 	

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);	

  if (uvar.help)	/* if help was requested, we're done here */
    exit (0);

  if ( uvar.version )	{
    /*OutputVersion();*/
    exit (0);
  }
   
  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(&uvar), &status);
   
  /* call the function that reads the needed information from the data (Fstat) file */
  LAL_CALL ( getFStat(uvar.inputFstat, &cfg), &status);
  
  /* call the function to compute Cstatistic */
  LAL_CALL ( computeCStat( &cfg, &uvar, &cst ), &status);
  
  /*some programming checks, can removed from real code*/
  printf("\n number of sidebands = %d \n", cst.sideband);
  printf("\n c[0] = %f \n", cst.Cstat[0]/cfg.N);
  printf("\n c[1] = %f \n", cst.Cstat[1]/cfg.N);
  printf ("\n orbit period = %le \t", uvar.orbitPeriod);
  printf ("\n df = %le \t", cfg.df);
  printf ("\n number of lines (cfg) = %d \t", cfg.N);
  printf ("\n fmin = %f \t", cfg.fmin);
  printf ("\n f0 = %f \t", cfg.f0);
  printf ("\n df = %le \t", cfg.df);    
  printf ("\n alpha = %lf \t", cfg.Alpha);     
  printf ("\n delta = %f \t", cfg.Delta);   
  printf("\n fstat: = \t %lf  \n", cfg.Fstat[0]);
  
  printf ("\n\n cfg f0 = %lf \t", cfg.f0);     
  printf("\n cfg dfout = \t %d  \n", (cfg.dfout));
  printf("\n cfg dm = \t %d  \n", (cfg.dm));
  printf("\n fstart=f[dfout] = \t %lf  \n", cfg.fmin + (cfg.dfout)*cfg.df);
  printf("\n f[dm]=Freq = \t %lf  \n", cfg.fmin+((cfg.dfout)+(cfg.dm))*cfg.df);
  printf("\n f[N-dm]=Freq+FreqBand = \t %lf  \n", cfg.fmin+((cfg.N-(2*(cfg.dm)))+(cfg.dfout)+(cfg.dm))*cfg.df );
  
  /* open output file */
  //LogSetLevel ( lalDebugLevel );
  if (LALUserVarWasSet(&uvar.outputCstat)) 	{
    if ((Cstat_out = fopen(uvar.outputCstat, "wb")) == NULL) 	{
      LogPrintf (LOG_CRITICAL, "%s: Error opening file '%s' for reading.. Error %d\n",fn,uvar.outputCstat,xlalErrno);
      return XLAL_EIO;
    }
    
    /* get full commandline describing search */
    LAL_CALL(LALUserVarGetLog(&status, &cmdline, UVAR_LOGFMT_CMDLINE), &status);
    fprintf(Cstat_out, "%%%% %s\n%%%% %s\n", rcsid, cmdline);
    LALFree(cmdline);
    
    /* output user input to file */
    fprintf(Cstat_out,"%%input: f0 = %f \t asin = %f \t fmin = %f \t",cfg.f0,uvar.orbitasini,uvar.Freq);
    fprintf(Cstat_out,"\n%% fmin = \t %f \n%% df = \t %f  \n%% alpha = \t %lf \n%% delta = %f  \n%% M = %d \t",cfg.fmin, cfg.df, cfg.Alpha, cfg.Delta, cst.sideband);
    fprintf(Cstat_out,"\n%% \n%% i \t frequency \t\t Fstat \t\t Cstat \n ");
    
    /* output fstat inputs and cstat output to file */
    for (i=0; i< cfg.N-(2*cfg.dm) ; i++) 	{
      fprintf(Cstat_out,"%d\t%6.12f\t%lf\t%lf\n",i,cfg.fmin+(i+cfg.dfout+cfg.dm)*cfg.df,cfg.Fstat[i+cfg.dfout+cfg.dm],(cst.Cstat[i+cfg.dm])/cfg.N);
    } 
    
  }
  fclose(Cstat_out);
 
  /* Free memory */
  LAL_CALL ( Freemem(&status, &cfg, &cst), &status); 
   
  /* did we forget anything ? */
  LALCheckMemoryLeaks();
    
  return 0;
  
} /* main() */



/*----------------------------------------------------------------------*/ 
/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
/*----------------------------------------------------------------------*/ 
void
initUserVars (LALStatus *status, UserInput_t *uvar)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->FreqBand = 0.0;
  uvar->orbitasini = 0.0;			/* define default orbital semi-major axis */
  uvar->help = FALSE;
  uvar->version = FALSE;
  uvar->outputCstat = NULL;
 
  /* register all user-variables */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 

  LALregREALUserStruct(status, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz"); 
  LALregREALUserStruct(status, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
   
  LALregREALUserStruct(status, 	orbitPeriod, 	'P',  UVAR_OPTIONAL, "Orbital period in seconds");
  LALregREALUserStruct(status, 	orbitasini, 	'A',  UVAR_OPTIONAL, "Orbital projected semi-major axis (normalised by the speed of light) in seconds [Default: 0.0]");
  
  LALregSTRINGUserStruct(status,inputFstat, 	'D', UVAR_REQUIRED, "Filename specifying input Fstat file"); 
  LALregSTRINGUserStruct(status,outputCstat,	'C', UVAR_REQUIRED, "Output-file for C-statistic");
   
  LALregBOOLUserStruct(status, version,	'V', UVAR_SPECIAL,  "Output version information");

 
  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */


/*----------------------------------------------------------------------*/
/* Free all globally allocated memory. */
void
Freemem(LALStatus *status,  ConfigVariables *cfg, CstatOut *cst) 
{
  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  LALFree (cfg->Fstat);
  LALFree (cst->Cstat);
  
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* Freemem() */



/*----------------------------------------------------------------------*/
/**
 * Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
/*----------------------------------------------------------------------*/
int checkUserInputConsistency ( const UserInput_t *uvar )
{
  const CHAR *fn = __func__;      		/* store function name for log output */ 

  /* binary parameter checks */
  if ( XLALUserVarWasSet(&uvar->orbitPeriod) && (uvar->orbitPeriod <= 0) )	{
    LogPrintf (LOG_CRITICAL, "%s: Negative or zero value of orbital period not allowed! Error %d\n",fn,xlalErrno);
    return XLAL_EDOM;
  }
  if ( XLALUserVarWasSet(&uvar->orbitasini) && (uvar->orbitasini < 0) )	{
    LogPrintf (LOG_CRITICAL, "%s: Negative or zero value of projected orbital semi-major axis not allowed! Error %d\n",fn,xlalErrno);
    return XLAL_EDOM;
  }
  
  return XLAL_SUCCESS; 
} /* checkUserInputConsistency() */



