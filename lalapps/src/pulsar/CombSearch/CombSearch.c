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
 * \ingroup pulsarApps
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
#define LAL_USE_OLD_COMPLEX_STRUCTS
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


/***********************************************************************************************/
/* defines */

#define MAXLINELENGTH 1024   		/* Maximum # of characters of in a line */
#define TRUE (1==1)
#define FALSE (1==0)

/***********************************************************************************************/
/* internal structures */

/** A structure for frequency series 
 */
typedef struct {
  REAL8Vector* 	fvect;			/**< frequency series of data, REAL8Vector has components data and length */
  REAL8		fmin;			/**< start frequency of data in Fseries */
  REAL8		df;			/**< frequency separation in Fseries */
  INT4		dof;			/**< degrees of freedom of detection statistic represented by frequency series */
  CHAR*		comment;		/**< comment field */ 
} VectorStruct;

/** A structure for template parameters 
 */
typedef struct {
  REAL8 freq;				/**< user defined start frequency */
  REAL8	fband;				/**< user defined frequency band */
  REAL8 orbitPeriod;			/**< orbital period of source */
  REAL8 orbitasini;			/**< light travel time across semi-major axis of source */
  REAL8	f0;				/**< search frequency, centre of user band, used to determine width of template */
  INT4	unitspikes;			/**< number of unit spikes in template */
} ParamStruct;

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            	/**< trigger to output help string */
  BOOLEAN version;			/**< output version information */
  BOOLEAN tophat;			/**< tophat template flag */
  CHAR *inputFstat;			/**< filename of Fstat input data file to use */
  CHAR *outputCstat;			/**< filename to output Cstatistic */
  REAL8 Freq;				/**< user defined start frequency */
  REAL8 FreqBand;			/**< user defined frequency band */
  /* orbital parameters */
  REAL8 orbitPeriod;			/**< binary-system orbital period in s */
  REAL8 orbitasini;			/**< amplitude of radial motion */
} UserInput_t;

/** A structure for information on the exclusion region needed to avoid edge effects in calculation of Cstat 
 */
typedef struct { 
  INT4 mm;		            	/**< number of sidebands on either side of central spike */
  INT4 dm;				/**< number of bins spanned by mm */
  INT4 fbins;				/**< number of bins in Fstat */
  INT4 cbins;				/**< number of bins in Cstat */
  REAL8 df;				/**< size of bin */
}ExRegion;

/***********************************************************************************************/
/* Global variables */
extern int vrbflg;				/* defined in lalapps.c */

const char *va(const char *format, ...);	/* little var-arg string helper function */

/* local prototypes */
int main(int argc,char *argv[]);
int initUserVars(int argc, char *argv[], UserInput_t *uvar);
int checkUserInputConsistency (const UserInput_t *uvar);
int ReadInput(UserInput_t *uvar, ParamStruct *userParams, VectorStruct *Fstat, ExRegion *exr);
int createComb(ParamStruct *userParams, VectorStruct *comb, ExRegion *exr);
int createTopHat(ParamStruct *userParams, VectorStruct *tophat, ExRegion *exr);
int ComputeCstat(VectorStruct *template, VectorStruct *Fstat, VectorStruct *Cstat, ExRegion *exr);
int OutputCstats(UserInput_t *uvar, ParamStruct *userParams, VectorStruct *Fstat, VectorStruct *Cstat, ExRegion *exr);

/* empty initializers */
UserInput_t empty_UserInput;
VectorStruct empty_VectorStruct;
ParamStruct empty_ParamStruct;
ExRegion empty_ExRegion;

/***********************************************************************************************/
/* Function definitions */


/*----------------------------------------------------------------------*/ 
/**
 * MAIN function of sb_search code.
 * Calculate the C-statistic over a given portion of the parameter-space
 * and write output into a file(default: 'Cstats').
 */
/*----------------------------------------------------------------------*/ 
int main(int argc,char *argv[])
{ 
  static const char *fn = __func__;             /* store function name for log output */
  UserInput_t uvar = empty_UserInput;		/* global container for user variables */
  ParamStruct userParams = empty_ParamStruct;	/* initialise paramStruct for search parameters, set by user */
  VectorStruct Fstat = empty_VectorStruct;	/* initialise Fstat structure */ 
  VectorStruct Cstat = empty_VectorStruct; 	/* initialise Cstat structure */ 
  VectorStruct template = empty_VectorStruct;	/* initialise template structure */ 
  ExRegion exr = empty_ExRegion;		/* initialise the exclusion region structure */
	

  lalDebugLevel = 0;				/* lalDebugLevel default */
  vrbflg = 1;					/* verbose error-messages */


  /* setup LAL debug level */
  if (XLALGetDebugLevel(argc, argv, 'v')) {
    LogPrintf(LOG_CRITICAL,"%s : XLALGetDebugLevel() failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }
  LogSetLevel(lalDebugLevel);
  LogPrintf(LOG_DEBUG,"Debug level set to %d \n", lalDebugLevel);
  
  /* register and read all user-variables */
  if (initUserVars(argc, argv, &uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : initUserVars failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }
  LogPrintf(LOG_DEBUG,"registered all user variables \n"); 
  
  /* do some sanity checks on the user-input before we proceed */
  if (checkUserInputConsistency(&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : checkUserInputConsistency failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }
  LogPrintf(LOG_DEBUG,"checked user input consistency \n");
  
  /* call function to get parameters and data from user input and data file */
  if (ReadInput(&uvar, &userParams, &Fstat, &exr)) {
    LogPrintf(LOG_CRITICAL,"%s : ReadInput failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }
  LogPrintf(LOG_DEBUG,"userParams and Fstat retrieved from input \n");

  /* check if tophat flag was raised */
  if (uvar.tophat) {
    /* call function to create tophat template */
    if (createTopHat(&userParams, &template, &exr)) {
      LogPrintf(LOG_CRITICAL,"%s : createTemplate failed with error = %d\n",fn,xlalErrno);
      return XLAL_EFAULT;
    }
  }
  else { 
    /* call function to create comb template */
    if (createComb(&userParams, &template, &exr)) {
      LogPrintf(LOG_CRITICAL,"%s : createComb failed with error = %d\n",fn,xlalErrno);
      return XLAL_EFAULT;
    }
  } 
  LogPrintf(LOG_DEBUG,"template created \n"); 

  /* call function to compute Cstat */
  if (ComputeCstat(&template, &Fstat, &Cstat, &exr)) {
    LogPrintf(LOG_CRITICAL,"%s : createComb failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }  
  LogPrintf(LOG_DEBUG,"ComputeCstat done \n"); 


  /* call function to output Cstat to file */
  if (OutputCstats(&uvar, &userParams, &Fstat, &Cstat, &exr)) {
    LogPrintf(LOG_CRITICAL,"%s : OutputCstats failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }  
  LogPrintf(LOG_DEBUG,"OutputCstats done \n");   
  
  /********************************/
  /* Free memory */
  /********************************/

  /* free data */
  XLALFree(Fstat.fvect->data);
  XLALFree(Cstat.fvect->data);
  XLALFree(template.fvect->data);
  XLALFree(Fstat.comment);
  XLALFree(Cstat.comment);
  
  /* free vectors */
  XLALFree(Fstat.fvect);
  XLALFree(Cstat.fvect);
  XLALFree(template.fvect);

  /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();
     
  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",fn);

  /* we're done */
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Done. \n",fn);  
  return 0;
  
} /* main() */

/*----------------------------------------------------------------------*/ 
/**
 * initUserVars function
 * Register "user-variables" specified from cmd-line and/or config-file.
 * Set defaults for some user-variables and register them with the UserInput module.
 */
/*----------------------------------------------------------------------*/ 
int initUserVars(int argc,char *argv[],UserInput_t *uvar)
{
  const CHAR *fn = __func__;      /* store function name for log output */
  
  CHAR *version_string;

  /* set a few defaults */
  uvar->help = FALSE;
  uvar->version = FALSE;

 
  /* register all user-variables */
  XLALregBOOLUserStruct(help, 		'h', UVAR_HELP,     "Print this message"); 
  XLALregREALUserStruct(Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz"); 
  XLALregREALUserStruct(FreqBand, 	'b', UVAR_REQUIRED, "Search frequency band in Hz"); 
  XLALregREALUserStruct(orbitPeriod, 	'P', UVAR_REQUIRED, "Orbital period in seconds");
  XLALregREALUserStruct(orbitasini, 	'A', UVAR_REQUIRED, "Light travel time of orbital projected semi-major axis, in seconds");
  XLALregSTRINGUserStruct(inputFstat, 	'D', UVAR_REQUIRED, "Filename specifying input Fstat file"); 
  XLALregSTRINGUserStruct(outputCstat,	'C', UVAR_REQUIRED, "Output-file for C-statistic");
  XLALregBOOLUserStruct(tophat,		't', UVAR_OPTIONAL, "Perform search with tophat template");
  XLALregBOOLUserStruct(version,	'V', UVAR_SPECIAL,  "Output version information");
  
  /* do ALL cmdline and cfgfile handling */
  if (XLALUserVarReadAllInput(argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput failed with error = %d\n",fn,xlalErrno);
    return XLAL_EFAULT;
  }
 
  /* if help was requested, we're done here */
  if (uvar->help) exit(0);

  if ((version_string = XLALGetVersionString(0)) == NULL) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }
  
  if (uvar->version) {
    printf("%s\n",version_string);
    exit(0);
  }
  XLALFree(version_string);
 
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn); 
  return XLAL_SUCCESS;

} /* initUserVars() */


/*----------------------------------------------------------------------*/
/**
 * Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
/*----------------------------------------------------------------------*/
int checkUserInputConsistency (const UserInput_t *uvar)
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
 
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn); 
  return XLAL_SUCCESS; 
  
} /* checkUserInputConsistency() */


/*--------------------------------------------------------------- */
/** ReadInput function 
 * reads user input and Fstat file and assigns parameter values and Fstat array 
 */
/*----------------------------------------------------------------*/
int ReadInput(UserInput_t *uvar, ParamStruct *userParams, VectorStruct *Fstat, ExRegion *exr)
{
  const CHAR *fn = __func__;      			/* store function name for log output */

  REAL8 *frequency = NULL, *fstats = NULL;   		/* temp storage for frequency and fstat vectors */
  REAL8 dummy, a=0, d=0, alpha=0, delta=0, fstat=0; 	/* temp variables to store data from file */ 
  INT4 c=0, l=0, i=0;					/* declare and initialise counters */
  UINT4 totlen=0;					/* recursive length counter for Fstat header */

  FILE *data = NULL;					/* data file */
  CHAR *filename = uvar->inputFstat;			/* Fstat input filename */
  CHAR line[MAXLINELENGTH];   				/* string to store Fstat header information */
  Fstat->comment = XLALCalloc(1,sizeof(CHAR));          /* comment field for Fstat header information - allocate initial memory */
  


  /* Open data file - check it is good */
  if ((data = fopen(filename, "r")) == NULL)	{
    LogPrintf (LOG_CRITICAL, "%s: Error opening file '%s' for reading.. Error %d\n",fn,filename,xlalErrno);
    return XLAL_EIO;
  }
  LogPrintf(LOG_DEBUG,"opened inputFstat file %s for reading \n",filename);

  /* go through file, store header information, count number of data lines (ignores comments), and check validity of data */
  while(fgets(line, MAXLINELENGTH, data) != NULL)	{ 
    /* get Fstat header information */
    if (strncmp(&line[0], "%%",2) == 0 ) {  
      UINT4 len = strlen(line);
      CHAR* tempstring = XLALCalloc(len+1,sizeof(CHAR));
      totlen += len;
      Fstat->comment = XLALRealloc(Fstat->comment,(totlen+1)*sizeof(CHAR));
      strncpy (tempstring, line, len);  
      strcat(Fstat->comment, tempstring);
      XLALFree(tempstring);
    }
    /* count data lines, check fstat validity and store last sky position entry */
    if (strncmp(&line[0], "%",1) != 0 && isdigit(line[0]) != 0) {
      sscanf (line, "%lf %lf %lf %lf %lf %lf %lf", &dummy, &alpha, &delta, &dummy, &dummy ,&dummy, &fstat);      
      l++;
      /*check if any fstats are negative, NaN or Inf*/
      if (fstat < 0 || isinf(fstat) || isnan(fstat)) {
        LogPrintf (LOG_CRITICAL, "%s: Fstat file %s contains negative, Inf or NaN values %d\n",fn,filename,xlalErrno);
        return XLAL_EDOM;
      }
    } 
  }
  /* close the file prior to exiting the routine */
  fclose(data); 

   /* Check data has more than one Fstat entry */
  if (l==1)	{
    LogPrintf (LOG_CRITICAL, "%s: Must be more than one Fstat entry in data file %s. Error %d\n",fn,filename,xlalErrno);
    return XLAL_EDOM;
  }
  LogPrintf(LOG_DEBUG,"checked for multiple Fstats, lines l= %d \n", l); 
  
  /* Allocate some memory according to number of frequency bins in input Fstat file */  
  if ((frequency=XLALCalloc(l,sizeof(REAL8))) == NULL )	{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory.Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }
  if ((fstats=XLALCalloc(l,sizeof(REAL8)))== NULL )		{
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }
  
  /* Open data file - check it is good */
  if ((data = fopen(filename, "r")) == NULL)	{
    LogPrintf (LOG_CRITICAL, "%s: Error opening file '%s' for reading.. Error %d\n",fn,filename,xlalErrno);
    return XLAL_EIO;
  } 
  
  /* Get data and assign some parameter variables */  
  while ( (fgets(line, MAXLINELENGTH, data)) != NULL)	{
    if ( ((strncmp(line,"%%",1)) != 0) && ((isdigit(line[0])) != 0) ) 	{
      sscanf (line, "%lf %lf %lf %lf %lf %lf %lf", &frequency[c], &a, &d, &dummy, &dummy ,&dummy, &fstats[c] );
      c++;      
      if ((alpha != a) || (delta != d)) {
        LogPrintf (LOG_CRITICAL, "%s: Sky position (alpha and delta) must be constant. Error %d\n",fn,xlalErrno);
        return XLAL_EINVAL;
      }
    }
  }
  fclose(data); /* close the file prior to exiting the routine */ 
  LogPrintf(LOG_DEBUG,"Acquired Fstats and search parameter information. File closed. \n"); 

   
  /* Assign all local variables */ 
  REAL8 f0, a0, df, Porb, freq, fband, fstart, fend;
  INT4 ufreq_bin, ufreqmax_bin, ufband, mm, dm, fbins, fstart_bin, fend_bin;

  /* Assign variables for easier reading */ 
  df 	= (frequency[c-1]-frequency[0])/(c-1); 			/* determine frequency resolution from Fstat file */
  Porb 	= userParams->orbitPeriod	= uvar->orbitPeriod;	/* orbital period of source */
  a0 	= userParams->orbitasini 	= uvar->orbitasini;	/* semi-major axis */
  freq 	= userParams->freq 		= uvar->Freq;		/* user input frequency */
  fband = userParams->fband 		= uvar->FreqBand;	/* user input frequency search band */


  /* get parameters - fill out userParams*/
  f0		= freq + 0.5*fband ; 				/* allocate guess frequency as centre of user input search frequency band */ 
  ufreq_bin	= floor(0.5 +(freq-frequency[0])/df);		/* nearest bin to user starting frequency */
  ufreqmax_bin	= floor(0.5 +(freq+fband-frequency[0])/df);  	/* nearest bin to user end frequency */
  ufband	= ufreqmax_bin - ufreq_bin +1 ; 		/* number of bins in user frequency band */
  mm		= floor(0.5 + LAL_TWOPI*f0*a0);			/* whole number of sidebands on either side of central spike */
  dm		= floor(0.5 + mm/(Porb*df));			/* exclusion region: number of Fstat bins to include on either side of user fband for Cstat calculation */
  fbins		= ufband + 2*dm;				/* length of Fstat vector requried for search*/ 
  fstart_bin    = ufreq_bin - dm ;				/* frequency bin of first Fstat required for search, half a comb width before user specified Freq */
  fstart	= frequency[0] + fstart_bin*df;			/* start frequency of search, half a comb width before user specified Freq */				
  fend_bin	= fstart_bin + ufband + 2*dm;			/* frequency bin corresponding to end of search, half a comb width after user specified Freq */	
  fend		= frequency[0] + fend_bin*df;			/* end frequency of search, half a comb width after user specified Freq */
  
    
  /* check search band plus exclusion region is within data frequency band */
  if ( (fstart < frequency[0]) || (frequency[c-1] < fend) )	{
      LogPrintf (LOG_CRITICAL, "%s: User input frequency range and/or exclusion range outside data limits. Error %d\n",fn,xlalErrno);
      return XLAL_EDOM;
  }
  LogPrintf(LOG_DEBUG,"Calculated exlusion region - within data limits. Good. \n"); 
  
  /* fill in exRegion struct exr */
  exr->mm			= mm;				/* sidebands either side of central spike - defines width of exlusion region */
  exr->dm			= dm;				/* number of bins in exclusion region */
  exr->fbins			= fbins;			/* number of bins required in Fstat */
  exr->cbins			= ufband;			/* number of bins to output Cstat */
  exr->df			= df;				/* frequency resolution */

  /* fill out rest of userParam struct */
  userParams->f0 		= f0;				/* guess frequency of search, taken as centre of user input search frequency band */
  if (uvar->tophat) {
    userParams->unitspikes	= floor(0.5+(2*mm+1)/(Porb*df)); 		/* number of unit amplitude spikes if tophat template flag is raised */
  }
  else userParams->unitspikes 	= 2*mm+1;			/* number of unit amplitude spikes in comb template (default) */


  /* Allocate some memory according to number of Fstat frequency bins required to complete search */  
  if ( (Fstat->fvect = XLALCreateREAL8Vector(fbins)) == NULL ){
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory to Fstat->fvect. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }

  /* get Fstats required for search - only need Fstats defined by user defined search band plus the exlusion region */
  for (i = 0; i < fbins; i++) {
    Fstat->fvect->data[i] = fstats[fstart_bin+i];
  }
  /* fill out the rest of the Fstat struct*/
  Fstat->fmin = fstart; 
  Fstat->df = df;
  Fstat->dof = 4;
 
  /* Free ReadInput function specific memory */
  XLALFree(frequency);
  XLALFree(fstats);    
  
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn);   
  return XLAL_SUCCESS;

} /* ReadInput */


/*--------------------------------------------------------------- */
/** createComb function 
 * use userParams to create comb template
 *
 * function creates template of unit amplitude spikes at 1 zero spike +/- mm spikes on either side separated by 1/P 
 */
/*----------------------------------------------------------------*/
int createComb(ParamStruct *userParams, VectorStruct *comb, ExRegion *exr)
{
  static const char *fn = __func__;             /* store function name for log output */

  INT4 i, ind=0;				/* declare and initialise counters*/

  /* Allocate some memory for comb template vector */  
  if ( (comb->fvect=XLALCreateREAL8Vector(exr->fbins)) == NULL ){
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory to comb->fvect. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }

  /* zero all values first */
  for (i=0;i<exr->fbins;i++) { 
    comb->fvect->data[i]=0;						
  }
  
  /* Assign unity at spacings determined by template type for 1 zero spike + mm positive frequency spikes */
  for (i=0;i<=exr->mm;i++)	{
    ind=floor(0.5+ i/(userParams->orbitPeriod*exr->df));
    comb->fvect->data[ind]=1;						
  }
  
  /*  assign unity at spacings of 1/P for mm negative frequency spikes */
  for (i=1;i<=exr->mm;i++)	{
    ind=(exr->fbins) -floor(0.5 +i/(userParams->orbitPeriod*exr->df));
    comb->fvect->data[ind]=1;						
  }
 
  /* allocate rest of comb template structure */
  comb->df 		= exr->df;
  comb->dof 		= userParams->unitspikes;
  
  LogPrintf(LOG_DEBUG,"%s:  template with %d spikes created \n",fn, comb->dof);
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn);
  return XLAL_SUCCESS;

} /* createComb */


/*--------------------------------------------------------------- */
/** createTopHat function 
 * use userParams to create tophat template if tophat flag is raised
 *
 * function creates template of unit amplitude spikes for a band mm/P on either side of the central zero spike 
 */
/*----------------------------------------------------------------*/
int createTopHat(ParamStruct *userParams, VectorStruct *tophat, ExRegion *exr)
{
  static const char *fn = __func__;             /* store function name for log output */

  INT4 i=0, ind=0;				/* declare and initialise counters*/
  
  /* Allocate some memory for tophat template vector */  
  if ( (tophat->fvect=XLALCreateREAL8Vector(exr->fbins)) == NULL ){
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory to tophat->fvect. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }

  /* zero all values first */
  for (i=0;i<exr->fbins;i++) { 
    tophat->fvect->data[i]=0;						
  }
  
  /* Assign unity at spacings for 1 zero spike + dm = mm/P*df positive frequency spikes */
  for (i=0;i<=exr->dm;i++)	{
    tophat->fvect->data[i]=1;						
  }
  
  /*  assign unity for dm = mm/P*df negative frequency spikes */
  for (i=1;i<=exr->dm;i++)	{
    ind=(exr->fbins -i);
    tophat->fvect->data[ind]=1;						
  }
 
  /* allocate rest of tophat template structure */
  tophat->df 		= exr->df;
  tophat->dof		= userParams->unitspikes;
  
   
  LogPrintf(LOG_DEBUG,"%s:  template with %d spikes created. halfwidth(dm) = %d bins\n",fn, tophat->dof, exr->dm);
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn);
  return XLAL_SUCCESS;

} /* createTopHat */




/*--------------------------------------------------------------- */
/** ComputeCstat function 
 * uses comb and Fstat structures to calculate the Cstat 
 */
/*----------------------------------------------------------------*/
int ComputeCstat(VectorStruct *template, VectorStruct *Fstat, VectorStruct *Cstat, ExRegion *exr)
{
  static const char *fn = __func__;             /* store function name for log output */

  INT4 i;					/* initialise counters */
  INT4 N = exr->fbins;   			/* number of frequency bins in Fstat */

  /* Allocate memory for FFT plans and vectors */
  REAL8FFTPlan *pfwd = NULL;			
  REAL8FFTPlan *prev = NULL;
  COMPLEX16Vector *f_out = NULL;
  COMPLEX16Vector *t_out = NULL;
  COMPLEX16Vector *c_out = NULL;
  REAL8Vector *cstats = NULL;


  /* Create FFT plans and vectors */
  pfwd=XLALCreateREAL8FFTPlan(N, 1, 0);
  prev=XLALCreateREAL8FFTPlan(N, 0, 0);
  cstats=XLALCreateREAL8Vector(N);
  f_out=XLALCreateCOMPLEX16Vector(N/2 +1);
  t_out=XLALCreateCOMPLEX16Vector(N/2 +1);
  c_out=XLALCreateCOMPLEX16Vector(N/2 +1);

  /* Fourier transform fstat array */
  XLALREAL8ForwardFFT(f_out, Fstat->fvect, pfwd);
  
  /* Fourier transform template for convolution with fstat Fourier transform */
  XLALREAL8ForwardFFT(t_out, template->fvect, pfwd);			

  /* Perform convolution of fstat with template by multiplication in Fourier time domain */
  for (i=0;i<(N/2 +1); i++)	{
    c_out->data[i].real_FIXME = (creal(f_out->data[i]) * creal(t_out->data[i])) - (f_out->data[i].im * t_out->data[i].im); /* real part of c_out */
    c_out->data[i].im = (creal(f_out->data[i]) * t_out->data[i].im) + (f_out->data[i].im * creal(t_out->data[i])); /* imaginary part of c_out */
   }

  /* Inverse FFT back to frequency domain to retrieve Cstat */
  XLALREAL8ReverseFFT(cstats, c_out, prev);
  
  /* Allocate some memory for Cstat vector */  
  if ( (Cstat->fvect=XLALCreateREAL8Vector(exr->cbins)) == NULL) {
    LogPrintf (LOG_CRITICAL, "%s: Error allocating memory to Cstat->fvect. Error %d\n",fn,xlalErrno);
    return XLAL_ENOMEM;
  }
							
  /* Assign Cstat arrary values of Cstatistic */
  for (i=0; i<exr->cbins; i++) {
    Cstat->fvect->data[i] = cstats->data[i+exr->dm];  
  }
  
  /* fill in the rest of the Cstat struct */
  Cstat->dof 		= Fstat->dof*template->dof;
  Cstat->fmin 		= Fstat->fmin + Fstat->df*exr->dm;
  Cstat->df		= Fstat->df;
					
  /* Destroy FFT plans and free ComputeCstat function specific memory */
  XLALDestroyREAL8FFTPlan(pfwd);
  XLALDestroyREAL8FFTPlan(prev);
  XLALDestroyREAL8Vector(cstats);
  XLALDestroyCOMPLEX16Vector(f_out);
  XLALDestroyCOMPLEX16Vector(t_out);
  XLALDestroyCOMPLEX16Vector(c_out);
  
  LogPrintf(LOG_DEBUG,"%s: Cstat dof = %d \n",fn, Cstat->dof);
  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn);
  return XLAL_SUCCESS;

} /* ComputeCstat */


int OutputCstats(UserInput_t *uvar, ParamStruct *userParams, VectorStruct *Fstat, VectorStruct *Cstat, ExRegion *exr)
{
  static const char *fn = __func__;            	/* store function name for log output */
  
  FILE *Cstat_out = NULL;			/* output filename */
  INT4 i;					/* initialise counters */  

  /* open output file - also check it is good */
  if (LALUserVarWasSet(&uvar->outputCstat)) 	{
    if ((Cstat_out = fopen(uvar->outputCstat, "wb")) == NULL) 	{
      LogPrintf (LOG_CRITICAL, "%s: Error opening file '%s' for reading.. Error %d\n",fn,uvar->outputCstat,xlalErrno);
      return XLAL_EIO;
    }
    
    /* get full command line describing search */
    if ( (Cstat->comment = XLALUserVarGetLog(UVAR_LOGFMT_CMDLINE)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALUserVarGetLog failed with error = %d\n",fn,xlalErrno);
      return XLAL_EFAULT;
    }
    /* print Fstat header information and Cstat command line to output file */
    fprintf(Cstat_out, "%%%% *********************************************************************** \n");
    fprintf(Cstat_out, "%%%% ***************   Fstat Header   ****************** \n %s", Fstat->comment);
    fprintf(Cstat_out, "%%%% *********************************************************************** \n");
    fprintf(Cstat_out, "%%%% ***************   Cstat Header   ****************** \n");
    /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
     *  It should be modified to use git version information. */
    fprintf(Cstat_out, "%%%% %s\n%%%% cmdline: %s\n", "$Id$", Cstat->comment);
    
    /* Cstat header - output user input to file */
    fprintf(Cstat_out,"%%%% input: f0 = %f,\t asini = %lf, Porb = %lf,\t user fmin = %f,\n",userParams->f0,userParams->orbitasini,userParams->orbitPeriod,userParams->freq);
    fprintf(Cstat_out,"%%%% Fstat params: fmin = %lf \t df = %2.16e  \t fbins = %d \t dof = %d \n",Fstat->fmin, Fstat->df, Fstat->fvect->length, Fstat->dof);
    fprintf(Cstat_out,"%%%% Cstat params: fmin = %lf \t df = %2.16e  \t cbins = %d \t dof = %d \n",Cstat->fmin, Cstat->df, Cstat->fvect->length, Cstat->dof);
    fprintf(Cstat_out,"%%%% \n%%%% i \t frequency \t\t Fstat \t\t Cstat \n ");
    
    /* output fstat inputs and cstat output to file */
    for (i=0; i< exr->cbins; i++) {
      fprintf(Cstat_out,"%d\t%6.12f\t%lf\t%lf\n",i,Cstat->fmin+i*Cstat->df,Fstat->fvect->data[i+exr->dm],(Cstat->fvect->data[i]/exr->fbins));
    }
    fprintf(Cstat_out,"%% Done"); 
  }
  /* close file */
  fclose(Cstat_out);

  LogPrintf(LOG_DEBUG,"'%s' successfully completed. Leaving. \n",fn);
  return XLAL_SUCCESS;
  
} /* OutputCstats */

