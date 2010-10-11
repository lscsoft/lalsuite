/*
 * Copyright (C) 2010 Chris Messenger
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

/** \author C.Messenger
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * This code is designed to combine the results files from multiple analyses using
 * lalapps_SemiCoherentBinary into a single Bayesian results file.
 *
 * It takes as input a set of Bayesian output files divided up into contiguous frequency
 * chunks and produces new Bayes factors and posterior disributions through interpolation. 
 *
 */

/***********************************************************************************************/
/* includes */
#include <math.h>
#include <time.h>
#include <glob.h>
#include <stdio.h>
#include <gsl/gsl_interp.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_spline.h>        /* needed for the gsl interpolation */
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>
#include <lal/AVFactories.h>

/***********************************************************************************************/
/* some global constants */

#define STRINGLENGTH 256              /* the length of general string */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NFREQMAX 4                    /* the max dimensionality of the frequency derivitive grid */
#define NBINMAX 4                     /* the number of binary parameter dimensions */
#define WINGS_FACTOR 1.05             /* the safety factor in reading extra frequency from SFTs */
#define PCU_AREA 0.13                 /* the collecting area of a single PCU in square metres */
#define DEFAULT_SOURCE "SCOX1"        /* the default source name */
#define AMPVECLENGTH 25               /* the fixed number of amplitude values to sample */
#define NLOGLUT 64                    /* the number of elements in the log LUT */
#define NBESSELLUT 256                /* the number of elements in the bessel function LUT */
#define OVERRES 5                     /* the over-sampling used in the interpolation of pdfs */ 

/***********************************************************************************************/
/* some useful macros */

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/***********************************************************************************************/
/* define internal structures */

/** A structure that stores information about an individual results file
 */
typedef struct { 
  CHAR filename[STRINGLENGTH];      /**< the name of the file */
  UINT4 tstart;                     /**< file start time */
  REAL8 tspan;                      /**< file time span */
  REAL8 tobs;                       /**< file observationm time */
  REAL8 minfreq;                    /**< the minimum frequency */
  REAL8 band;                       /**< the frequency band */
  REAL8 mm;                         /**< the mismatch */
  UINT4 nseg;
  UINT4 ndim;
  REAL8 Bayes;
  REAL8 Bayes_fixed;
  /* REAL8 min_amp; */
/*   REAL8 max_amp; */
/*   REAL8 sig_amp; */
/*   REAL8 start_amp; */
/*   REAL8 delta_amp; */
/*   UINT4 length_amp; */
/*   CHAR prior_amp[STRINGLENGTH]; */
/*   REAL8Vector *logposterior_fixed_amp; */
  REAL8Vector *Bayes_perseg;
  UINT4Vector *Bayes_start;
  UINT4Vector *Bayes_end;
  CHAR **name;
  REAL8 min[NBINMAX+1];
  REAL8 max[NBINMAX+1];
  REAL8 sig[NBINMAX+1];
  REAL8 start[NBINMAX+1];
  REAL8 delta[NBINMAX+1];
  UINT4 length[NBINMAX+1];
  CHAR **prior;
  REAL8Vector *logposterior[NBINMAX+1];
  REAL8Vector *logposterior_fixed[NBINMAX+1];
} BayesianResultsFile;

/** A structure that stores information about a collection of Bayesian results files
 */
typedef struct { 
  UINT4 length;	                    /**< the number of channels */
  CHAR dir[STRINGLENGTH];           /**< the directory containing the files */
  BayesianResultsFile *file;        /**< a pointer to Bayesian Results file structures */
} BayesianResultsFileList;

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *inputdir;                   /**< the input directory */
  CHAR *outputdir;                  /**< the output directory */
  CHAR *source;                     /**< the name of the source */
  CHAR *obsid_pattern;              /**< the OBS ID substring */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */
RCSID( "$Id$");		        /* FIXME: use git-ID instead to set 'rcsid' */

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar, CHAR **clargs);
int XLALReadResultsDir(BayesianResultsFileList **resultsfiles, CHAR *inputdir, CHAR *pattern);
int XLALCombineBayesianResults(BayesianResultsFile **combinedresult,BayesianResultsFileList *resultsfiles);
REAL8 XLALLogSumExp(REAL8 logx,REAL8 logy);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;

/** The main function of semicoherentbinary.c
 *
 */
int main( int argc, char *argv[] )  {

  static const char *fn = __func__;               /* store function name for log output */
  UserInput_t uvar = empty_UserInput;             /* user input variables */
  CHAR *clargs = NULL;                            /* store the command line args */
  BayesianResultsFileList *resultsfiles = NULL;   /* the input results files */
  BayesianResultsFile *combinedresult = NULL;     /* the output combined result */
  UINT4 i,j;                                      /* counters */

  lalDebugLevel = 1;
  vrbflg = 1;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* setup LAL debug level */
  if (XLALGetDebugLevel(argc, argv, 'v')) {
    LogPrintf(LOG_CRITICAL,"%s : XLALGetDebugLevel() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  LogSetLevel(lalDebugLevel);

  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar,&clargs)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadUserVars() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",fn);

  /**********************************************************************************/
  /* READ FILES */
  /**********************************************************************************/

  /* get a list of results file names */
  if (XLALReadResultsDir(&resultsfiles,uvar.inputdir,uvar.obsid_pattern)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadResultsDir() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }

  /**********************************************************************************/
  /* COMBINE THE FILES */
  /**********************************************************************************/

  /* combine the results into a single result */
  if (XLALCombineBayesianResults(&combinedresult,resultsfiles)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCombineBayesianResults() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
 
  /**********************************************************************************/
  /* OUTPUT RESULTS TO FILE */
  /**********************************************************************************/

 
  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* free filelist params */
  for (i=0;i<resultsfiles->length;i++) {
    XLALDestroyREAL8Vector(resultsfiles->file[i].Bayes_perseg);
    XLALDestroyUINT4Vector(resultsfiles->file[i].Bayes_start);
    XLALDestroyUINT4Vector(resultsfiles->file[i].Bayes_end);
    for (j=0;j<NBINMAX+1;j++) {
      XLALFree(resultsfiles->file[i].name[j]);
      XLALFree(resultsfiles->file[i].prior[j]);
      XLALDestroyREAL8Vector(resultsfiles->file[i].logposterior[j]);
      XLALDestroyREAL8Vector(resultsfiles->file[i].logposterior_fixed[j]);
    }
    XLALFree(resultsfiles->file[i].name);
    XLALFree(resultsfiles->file[i].prior);
  }
  XLALFree(resultsfiles->file);
  XLALFree(resultsfiles);

  /* free combined results */
  XLALDestroyREAL8Vector(combinedresult->Bayes_perseg);
  XLALDestroyUINT4Vector(combinedresult->Bayes_start);
  XLALDestroyUINT4Vector(combinedresult->Bayes_end);
  XLALFree(combinedresult);

  /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();
  XLALFree(clargs);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",fn);

  LogPrintf(LOG_DEBUG,"%s : successfully completed.\n",fn);
  return 0;
  
} /* end of main */

/** Read in input user arguments
 *
 */
int XLALReadUserVars(int argc,            /**< [in] the command line argument counter */ 
		     char *argv[],        /**< [in] the command line arguments */
		     UserInput_t *uvar,   /**< [out] the user input structure */
		     CHAR **clargs        /**< [out] the command line args string */
		     )
{

  const CHAR *fn = __func__;   /* store function name for log output */
  CHAR *version_string;
  INT4 i;

  /* initialise user variables */
  uvar->obsid_pattern = NULL;

  /* initialise default source as SCOX1 */
  {
    UINT4 n = strlen(DEFAULT_SOURCE) + 1;
    uvar->source = XLALCalloc(n,sizeof(CHAR));
    snprintf(uvar->source,n,"%s",DEFAULT_SOURCE);
  }

  /* ---------- register all user-variables ---------- */
  XLALregBOOLUserStruct(help, 		        'h', UVAR_HELP,     "Print this message");
  XLALregSTRINGUserStruct(inputdir, 	        'i', UVAR_REQUIRED, "The input directory name"); 
  XLALregSTRINGUserStruct(outputdir, 	        'o', UVAR_REQUIRED, "The output directory name"); 
  XLALregSTRINGUserStruct(source, 	        'x', UVAR_OPTIONAL, "The source name (default SCOX1)"); 
  XLALregSTRINGUserStruct(obsid_pattern,        'O', UVAR_OPTIONAL, "The observation ID substring to match"); 
  XLALregBOOLUserStruct(version,                'V', UVAR_SPECIAL,  "Output code version");

  /* do ALL cmdline and cfgfile handling */
  if (XLALUserVarReadAllInput(argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput() failed with error = %d\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_EINVAL);
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

  /* put clargs into string */
  *clargs = XLALCalloc(1,sizeof(CHAR));
  for (i=0;i<argc;i++) {
    INT4 len = 2 + strlen(argv[i]) + strlen(*clargs);
    *clargs = XLALRealloc(*clargs,len*sizeof(CHAR));
    strcat(*clargs,argv[i]);
    strcat(*clargs," ");
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** Read in list of Bayesian results files
 */
int XLALReadResultsDir(BayesianResultsFileList **resultsfiles,     /**< [out] a structure containing a list of all input results files */
		       CHAR *inputdir,                             /**< [in] the input results directory */
		       CHAR *pattern                               /**< [in] the input APID pattern */
		       )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */
  INT4 nfiles;                    /* the number of returned files from glob */
  INT4 i,j,k;                     /* counters */
  glob_t pglob;
  CHAR glob_pattern[STRINGLENGTH];

  /* check input arguments */
  if ((*resultsfiles) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input results file list structure != NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (inputdir == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input directory string == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  

  /* if we have a filename patttern then set the global variable */
  if (pattern != NULL) snprintf(glob_pattern,STRINGLENGTH,"%s/*%s*.txt",inputdir,pattern); 
  else snprintf(glob_pattern,STRINGLENGTH,"%s/*.txt",inputdir); 
  LogPrintf(LOG_DEBUG,"%s : searching for file pattern %s\n",fn,glob_pattern);
  
  /* get the mode1 FITS file (apID 55) */
  if (glob(glob_pattern,0,NULL,&pglob)) {
    LogPrintf(LOG_CRITICAL,"%s : glob() failed to return a filelist.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  nfiles = pglob.gl_pathc;
  if (nfiles>0) {
    for (i=0;i<nfiles;i++) LogPrintf(LOG_DEBUG,"%s : found file %s\n",fn,pglob.gl_pathv[i]);
  }
  else if (nfiles == -1) {
    LogPrintf(LOG_CRITICAL,"%s : glob() failed.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  else if (nfiles == 0) {
    LogPrintf(LOG_CRITICAL,"%s : could not find any results files in directory %s matching pattern %s.\n",fn,inputdir,pattern);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  /* allocate memory for filelist structure */
  if (((*resultsfiles) = (BayesianResultsFileList *)XLALCalloc(1,sizeof(BayesianResultsFileList))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for results filelist structure.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*resultsfiles)->file = (BayesianResultsFile *)XLALCalloc(nfiles,sizeof(BayesianResultsFile))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for results file structures.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  snprintf((*resultsfiles)->dir,STRINGLENGTH,"%s",inputdir);
  (*resultsfiles)->length = nfiles;

  /* open each file and read in header information */
  for (i=0;i<nfiles;i++) {
    
    FILE *fp = NULL;
    CHAR line[STRINGLENGTH];
    CHAR *c = NULL;
    
    /* open the frame file */
    if ((fp = fopen(pglob.gl_pathv[i],"r")) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : unable to open file %s.\n",fn,pglob.gl_pathv[i]);
      XLAL_ERROR(fn,XLAL_EINVAL);
    } 
    LogPrintf(LOG_DEBUG,"%s : opened file %s.\n",fn,pglob.gl_pathv[i]);
      
    /* allocate memory for the parameter name and prior type */
    if (((*resultsfiles)->file[i].name = (CHAR **)XLALCalloc(NBINMAX+1,sizeof(CHAR*))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for name of parameter.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    if (((*resultsfiles)->file[i].prior = (CHAR **)XLALCalloc(NBINMAX+1,sizeof(CHAR*))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for type of prior.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    for (j=0;j<NBINMAX+1;j++) {
      if (((*resultsfiles)->file[i].name[j] = (CHAR *)XLALCalloc(STRINGLENGTH,sizeof(CHAR))) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for name of parameter.\n",fn);
	XLAL_ERROR(fn,XLAL_ENOMEM);
      }
      if (((*resultsfiles)->file[i].prior[j] = (CHAR *)XLALCalloc(STRINGLENGTH,sizeof(CHAR))) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for type of prior.\n",fn);
	XLAL_ERROR(fn,XLAL_ENOMEM);
      }
    }

    /* read in header information */
    while (fgets(line,sizeof(line),fp)) {   
      
      if ((c = strstr(line,"start time (GPS sec)"))) sscanf(strstr(c,"=")+2,"%d",&((*resultsfiles)->file[i].tstart));
      else if ((c = strstr(line,"observation span (sec)"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].tspan));
      else if ((c = strstr(line,"coherent time (sec)"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].tobs));
      else if ((c = strstr(line,"number of segments"))) sscanf(strstr(c,"=")+2,"%d",&((*resultsfiles)->file[i].nseg));
      else if ((c = strstr(line,"number of dimensions"))) sscanf(strstr(c,"=")+2,"%d",&((*resultsfiles)->file[i].ndim));
      else if ((c = strstr(line,"mismatch"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].mm));
      else if ((c = strstr(line,"log Bayes Factor (phase and amplitude marginalised per segment)	="))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].Bayes));
      else if ((c = strstr(line,"log Bayes Factor (phase marginalised per segment)"))) {
	sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].Bayes_fixed));
	
	/* allocate mem for the Bayesfactors per segment */
	if (((*resultsfiles)->file[i].Bayes_perseg = XLALCreateREAL8Vector((*resultsfiles)->file[i].nseg)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for Bayes factor vector.\n",fn);
	  XLAL_ERROR(fn,XLAL_ENOMEM);
	}
	if (((*resultsfiles)->file[i].Bayes_start = XLALCreateUINT4Vector((*resultsfiles)->file[i].nseg)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for Bayes factor vector.\n",fn);
	  XLAL_ERROR(fn,XLAL_ENOMEM);
	}
	if (((*resultsfiles)->file[i].Bayes_end = XLALCreateUINT4Vector((*resultsfiles)->file[i].nseg)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for Bayes factor vector.\n",fn);
	  XLAL_ERROR(fn,XLAL_ENOMEM);
	}

	/* if we've found the Bayes factor information line then the actual individual Bayes factors are 3 lines later */
	for (k=0;k<3;k++) fgets(line,sizeof(line),fp);
	for (k=0;k<(INT4)((*resultsfiles)->file[i].nseg);k++) {
	  fgets(line,sizeof(line),fp);
	  sscanf(line,"%d %d %lf",&((*resultsfiles)->file[i].Bayes_start->data[k]),
		 &((*resultsfiles)->file[i].Bayes_end->data[k]),
		 &((*resultsfiles)->file[i].Bayes_perseg->data[k]));
	  printf("read Bayes factors as %d %d %lf\n",((*resultsfiles)->file[i].Bayes_start->data[k]),
		 ((*resultsfiles)->file[i].Bayes_end->data[k]),
		 ((*resultsfiles)->file[i].Bayes_perseg->data[k]));
	}
	printf("***\n");

      }

     /*  else if ((c = strstr(line,"min_"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].min_amp)); */
/*       else if ((c = strstr(line,"max_amp"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].max_amp)); */
/*       else if ((c = strstr(line,"sig_amp"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].sig_amp)); */
/*       else if ((c = strstr(line,"start_amp"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].start_amp)); */
/*       else if ((c = strstr(line,"delta_amp"))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].delta_amp)); */
/*       else if ((c = strstr(line,"length_amp"))) { */
/* 	sscanf(strstr(c,"=")+2,"%d",&((*resultsfiles)->file[i].length_amp)); */
	
/* 	/\* allocate mem for the pdfs *\/ */
/* 	if (((*resultsfiles)->file[i].logposterior_fixed_amp = XLALCreateREAL8Vector((*resultsfiles)->file[i].length_amp)) == NULL) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for posterior vector.\n",fn); */
/* 	  XLAL_ERROR(fn,XLAL_ENOMEM); */
/* 	} */
	
/*       } */
/*       else if ((c = strstr(line,"prior_amp"))) { */
/* 	sscanf(strstr(c,"=")+2,"%s",((*resultsfiles)->file[i].prior_amp)); */
	
/* 	/\* if we've found the prior information line then the actual pdf is 3 lines later *\/ */
/* 	for (k=0;k<3;k++) fgets(line,sizeof(line),fp); */
/* 	for (k=0;k<(INT4)((*resultsfiles)->file[i].length_amp);k++) { */
/* 	  fgets(line,sizeof(line),fp); */
/* 	  sscanf(line,"%*e %le %*e %*e %*e %*e",&((*resultsfiles)->file[i].logposterior_fixed_amp->data[k])); */
/* 	  printf("read log posteriors as %lf\n",((*resultsfiles)->file[i].logposterior_fixed_amp->data[k])); */
/* 	} */
/* 	printf("***\n"); */
/*       } */
	
      /* loop over search parameters */
      for (j=0;j<NBINMAX+1;j++) {

	CHAR name_string[STRINGLENGTH];
	CHAR min_string[STRINGLENGTH];
	CHAR max_string[STRINGLENGTH];
	CHAR sig_string[STRINGLENGTH];
	CHAR start_string[STRINGLENGTH];
	CHAR delta_string[STRINGLENGTH];
	CHAR length_string[STRINGLENGTH];
	CHAR prior_string[STRINGLENGTH];
	
	sprintf(name_string,"name_%d",j);
	sprintf(min_string,"min_%d",j);
	sprintf(max_string,"max_%d",j);
	sprintf(sig_string,"sig_%d",j);
	sprintf(start_string,"start_%d",j);
	sprintf(delta_string,"delta_%d",j);
	sprintf(length_string,"length_%d",j);
	sprintf(prior_string,"prior_%d",j);	

	if ((c = strstr(line,name_string))) sscanf(strstr(c,"=")+2,"%s",((*resultsfiles)->file[i].name[j]));
	else if ((c = strstr(line,min_string))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].min[j]));
	else if ((c = strstr(line,max_string))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].max[j]));
	else if ((c = strstr(line,sig_string))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].sig[j]));
	else if ((c = strstr(line,start_string))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].start[j]));
	else if ((c = strstr(line,delta_string))) sscanf(strstr(c,"=")+2,"%lf",&((*resultsfiles)->file[i].delta[j]));
	else if ((c = strstr(line,length_string))) {
	  sscanf(strstr(c,"=")+2,"%d",&((*resultsfiles)->file[i].length[j]));
	  
	  /* allocate mem for the pdfs */
	  if (((*resultsfiles)->file[i].logposterior[j] = XLALCreateREAL8Vector((*resultsfiles)->file[i].length[j])) == NULL) {
	    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for posterior vector.\n",fn);
	    XLAL_ERROR(fn,XLAL_ENOMEM);
	  }
	  if (((*resultsfiles)->file[i].logposterior_fixed[j] = XLALCreateREAL8Vector((*resultsfiles)->file[i].length[j])) == NULL) {
	    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for posterior vector.\n",fn);
	    XLAL_ERROR(fn,XLAL_ENOMEM);
	  }
	  
	}
	else if ((c = strstr(line,prior_string))) {
	  sscanf(strstr(c,"=")+2,"%s",((*resultsfiles)->file[i].prior[j]));

	  /* if we've found the prior information line then the actual pdf is 3 lines later */
	  for (k=0;k<3;k++) fgets(line,sizeof(line),fp);
	  for (k=0;k<(INT4)((*resultsfiles)->file[i].length[j]);k++) {
	    fgets(line,sizeof(line),fp);
	    sscanf(line,"%*e %le %*e %le %*e %*e",&((*resultsfiles)->file[i].logposterior[j]->data[k]),&((*resultsfiles)->file[i].logposterior_fixed[j]->data[k]));
	    /* printf("read log posteriors as %lf %lf\n",((*resultsfiles)->file[i].logposterior[j]->data[k]),((*resultsfiles)->file[i].logposterior_fixed[j]->data[k])); */
	  }
	  /* printf("***\n"); */
	}
	
      } /* end loop over parameters */
      
    } /* end loop over lines */
 
    /* fill in extra info */
    snprintf((*resultsfiles)->file[i].filename,STRINGLENGTH,"%s",pglob.gl_pathv[i]);

    LogPrintf(LOG_DEBUG,"%s : read tstart as %d\n",fn,((*resultsfiles)->file[i].tstart));
    LogPrintf(LOG_DEBUG,"%s : read tspan as %f\n",fn,((*resultsfiles)->file[i].tspan));
    LogPrintf(LOG_DEBUG,"%s : read tobs as %f\n",fn,((*resultsfiles)->file[i].tobs));
    LogPrintf(LOG_DEBUG,"%s : read nseg as %d\n",fn,((*resultsfiles)->file[i].nseg));
    LogPrintf(LOG_DEBUG,"%s : read ndim as %d\n",fn,((*resultsfiles)->file[i].ndim));
    LogPrintf(LOG_DEBUG,"%s : read mismatch as %f\n",fn,((*resultsfiles)->file[i].mm));
    LogPrintf(LOG_DEBUG,"%s : read Bayes as %f\n",fn,((*resultsfiles)->file[i].Bayes));
    LogPrintf(LOG_DEBUG,"%s : read Bayes_fixed as %f\n",fn,((*resultsfiles)->file[i].Bayes_fixed));
  
   /*  LogPrintf(LOG_DEBUG,"%s : read amp min as %f\n",fn,((*resultsfiles)->file[i].min_amp)); */
/*     LogPrintf(LOG_DEBUG,"%s : read amp max as %f\n",fn,((*resultsfiles)->file[i].max_amp)); */
/*     LogPrintf(LOG_DEBUG,"%s : read amp sig as %f\n",fn,((*resultsfiles)->file[i].sig_amp)); */
/*     LogPrintf(LOG_DEBUG,"%s : read amp start as %f\n",fn,((*resultsfiles)->file[i].start_amp)); */
/*     LogPrintf(LOG_DEBUG,"%s : read amp delta as %f\n",fn,((*resultsfiles)->file[i].delta_amp)); */
/*     LogPrintf(LOG_DEBUG,"%s : read amp length as %d\n",fn,((*resultsfiles)->file[i].length_amp)); */
/*     LogPrintf(LOG_DEBUG,"%s : read amp prior as %s\n",fn,((*resultsfiles)->file[i].prior_amp)); */
    
    for (j=0;j<NBINMAX+1;j++) {
      LogPrintf(LOG_DEBUG,"%s : read name[%d] as %s\n",fn,j,((*resultsfiles)->file[i].name[j]));
      LogPrintf(LOG_DEBUG,"%s : read min[%d] as %f\n",fn,j,((*resultsfiles)->file[i].min[j]));
      LogPrintf(LOG_DEBUG,"%s : read max[%d] as %f\n",fn,j,((*resultsfiles)->file[i].max[j]));
      LogPrintf(LOG_DEBUG,"%s : read sig[%d] as %f\n",fn,j,((*resultsfiles)->file[i].sig[j]));
      LogPrintf(LOG_DEBUG,"%s : read start[%d] as %f\n",fn,j,((*resultsfiles)->file[i].start[j]));
      LogPrintf(LOG_DEBUG,"%s : read delta[%d] as %f\n",fn,j,((*resultsfiles)->file[i].delta[j]));
      LogPrintf(LOG_DEBUG,"%s : read length[%d] as %d\n",fn,j,((*resultsfiles)->file[i].length[j]));
      LogPrintf(LOG_DEBUG,"%s : read prior[%d] as %s\n",fn,j,((*resultsfiles)->file[i].prior[j]));
    }

    /* close the frame file */
    fclose(fp);
  
  }
  
  /* check for consistency between files e.g. all must have been performed */
  /* over the same range in all parameters (except frequency) and have used */
  /* the same data and mismatch.  Exit if any inconsistencies are found. */
  for (i=1;i<(INT4)(*resultsfiles)->length;i++) {
    
    if ((*resultsfiles)->file[i].tstart != (*resultsfiles)->file[0].tstart) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent start times for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    if ((*resultsfiles)->file[i].tspan != (*resultsfiles)->file[0].tspan) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent time spans for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    if ((*resultsfiles)->file[i].tobs != (*resultsfiles)->file[0].tobs) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent observation times for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    if ((*resultsfiles)->file[i].nseg != (*resultsfiles)->file[0].nseg) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent numbers of segments for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    if ((*resultsfiles)->file[i].ndim != (*resultsfiles)->file[0].ndim) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent number of search dimensions for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    if ((*resultsfiles)->file[i].mm != (*resultsfiles)->file[0].mm) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent mismatches for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }

    /* check parameter space consistency */
    for (j=0;j<NBINMAX+1;j++) {
      
      /* if not the frequency parameter */
      if (strcmp((*resultsfiles)->file[i].name[j],"nu")) {

	if ((*resultsfiles)->file[i].min[j] != (*resultsfiles)->file[0].min[j]) {
	  LogPrintf(LOG_CRITICAL,"%s : inconsistent min values for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}
	if ((*resultsfiles)->file[i].max[j] != (*resultsfiles)->file[0].max[j]) {
	  LogPrintf(LOG_CRITICAL,"%s : inconsistent max values for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}
	if ((*resultsfiles)->file[i].sig[j] != (*resultsfiles)->file[0].sig[j]) {
	  LogPrintf(LOG_CRITICAL,"%s : inconsistent sigma values for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}
	if (strcmp((*resultsfiles)->file[i].prior[j],(*resultsfiles)->file[0].prior[j])) {
	  LogPrintf(LOG_CRITICAL,"%s : inconsistent prior types for files %s and %s.\n",fn,(*resultsfiles)->file[0].filename,(*resultsfiles)->file[i].filename);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}

      }

    }
  }

  /* free original filelist */
  globfree(&pglob);
 
  /* read in the posteriors */
  exit(0);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** Combines the input results files into a single consistent result file
 *
 */
int XLALCombineBayesianResults(BayesianResultsFile **combinedresult,
			       BayesianResultsFileList *resultsfiles
			       )
{

  const CHAR *fn = __func__;      /* store function name for log output */
  UINT4 i,j,k;                    /* counters */
  REAL8 totalband = 0.0;          /* the total bandwidth */
  REAL8 *deltaband;               /* the bandwidths of each file */

  /* check input arguments */
  if ((*combinedresult) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input combinedresult structure != NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (resultsfiles == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, resultsfiles structure == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* allocate memory for the combined result */
  if (((*combinedresult) = (BayesianResultsFile *)XLALCalloc(1,sizeof(BayesianResultsFile))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output Bayesian results.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* first deal with header info */
  (*combinedresult)->nseg = 0;
  for (i=0;i<resultsfiles->length;i++) (*combinedresult)->nseg += resultsfiles->file[i].nseg;

  /* allocate mem for bandwidths */
  if ((deltaband = (REAL8 *)XLALCalloc(resultsfiles->length,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for delta band widths.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* compute entire bandwidth */
  totalband = 0.0;
  for (i=0;i<resultsfiles->length;i++) {   
    for (j=0;j<resultsfiles->file[i].ndim;j++) {
      if (strcmp(resultsfiles->file[i].name[j],"nu") == 0) {
	deltaband[i] = resultsfiles->file[i].max[j]-resultsfiles->file[i].min[j];
	printf("deltaband = %f\n",deltaband[i]);
	totalband += deltaband[i];
      }
    }
  }
  printf("totalband = %f\n",totalband);

  /* allocate memory for the combined result */
  if (((*combinedresult)->Bayes_perseg = XLALCreateREAL8Vector(resultsfiles->file[0].nseg)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output Bayesian segment results.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*combinedresult)->Bayes_start = XLALCreateUINT4Vector(resultsfiles->file[0].nseg)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output Bayesian segment results.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*combinedresult)->Bayes_end = XLALCreateUINT4Vector(resultsfiles->file[0].nseg)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output Bayesian segment results.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* combine Bayes factors per segment - remember that they are log(B) so have to be added carefully */
  /* they are also weighted by their bandwidth */
  for (j=0;j<(*combinedresult)->Bayes_perseg->length;j++) {

    (*combinedresult)->Bayes_perseg->data[j] = resultsfiles->file[0].Bayes_perseg->data[j] + log(deltaband[0]) - log(totalband);
    (*combinedresult)->Bayes_start->data[j] = resultsfiles->file[0].Bayes_start->data[j];
    (*combinedresult)->Bayes_end->data[j] = resultsfiles->file[0].Bayes_end->data[j];

    for (i=1;i<resultsfiles->length;i++) {
            
      REAL8 temp1 = (*combinedresult)->Bayes_perseg->data[j];
      REAL8 temp2 = resultsfiles->file[i].Bayes_perseg->data[j] + log(deltaband[i]) - log(totalband);
      
      (*combinedresult)->Bayes_perseg->data[j] = XLALLogSumExp(temp1,temp2);
    }
    LogPrintf(LOG_DEBUG,"%s : combined individual segment Bayes factors %d %d %lf\n",fn,(*combinedresult)->Bayes_start->data[j],
	   (*combinedresult)->Bayes_end->data[j],
	   (*combinedresult)->Bayes_perseg->data[j]);
    
  }

  /* combine overall Bayes factors - remember that they are log(B) so have to be added carefully */
  /* they are also weighted by their bandwidth */
  (*combinedresult)->Bayes = resultsfiles->file[0].Bayes + log(deltaband[0]) - log(totalband);
  (*combinedresult)->Bayes_fixed = resultsfiles->file[0].Bayes_fixed + log(deltaband[0]) - log(totalband);
  for (i=1;i<resultsfiles->length;i++) {
            
      REAL8 temp1 = (*combinedresult)->Bayes;
      REAL8 temp2 = resultsfiles->file[i].Bayes + log(deltaband[i]) - log(totalband);
      REAL8 temp3 = (*combinedresult)->Bayes_fixed;
      REAL8 temp4 = resultsfiles->file[i].Bayes_fixed + log(deltaband[i]) - log(totalband);
      
      (*combinedresult)->Bayes = XLALLogSumExp(temp1,temp2);
      (*combinedresult)->Bayes_fixed = XLALLogSumExp(temp3,temp4);
  }
  LogPrintf(LOG_DEBUG,"%s : combined Bayes factors %lf %lf\n",fn,(*combinedresult)->Bayes,(*combinedresult)->Bayes_fixed);

  /************************************************/
  /* combine posteriors */

  /* loop over the number of dimensions */

 /*  /\* first find out the pdf with the most samples *\/ */
/*   { */
/*     UINT4 maxlen = 0; */
/*     for (i=0;i<resultsfiles->length;i++) { */
/*       if (resultsfiles->file[i].logposterior_fixed_amp->length>maxlen) maxlen = resultsfiles->file[i].logposterior_fixed_amp->length; */
/*     } */
/*     printf("maxlen = %d\n",maxlen); */
    
/*     /\* we define an over-resolved interpolation scheme with OVERRES more points than the max value *\/ */
/*     if (((*combinedresult)->logposterior_fixed_amp = XLALCreateREAL8Vector((UINT4)(0.5 + OVERRES*maxlen))) == NULL) { */
/*       LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output amplitude posterior results.\n",fn); */
/*       XLAL_ERROR(fn,XLAL_ENOMEM); */
/*     } */

/*     /\* initialise results *\/ */
/*     for (j=0;j<(*combinedresult)->logposterior_fixed_amp->length;j++) (*combinedresult)->logposterior_fixed_amp->data[j] = -1e200; */

/*     /\* fill in some header info *\/ */
/*     (*combinedresult)->start_amp = resultsfiles->file[0].start_amp; */
/*     (*combinedresult)->length_amp = (*combinedresult)->logposterior_fixed_amp->length; */
/*     (*combinedresult)->delta_amp = ((resultsfiles->file[0].length_amp-1)*resultsfiles->file[0].delta_amp)/((*combinedresult)->length_amp - 1); */
/*     (*combinedresult)->ndim = resultsfiles->file[0].ndim; */
/*     (*combinedresult)->nseg = resultsfiles->file[0].nseg; */
/*     (*combinedresult)->band = totalband; */
    

   /*   CHAR filename[STRINGLENGTH];      /\**< the name of the file *\/ */
/*   UINT4 tstart;                     /\**< file start time *\/ */
/*   REAL8 tspan;                      /\**< file time span *\/ */
/*   REAL8 tobs;                       /\**< file observationm time *\/ */
/*   REAL8 minfreq;                    /\**< the minimum frequency *\/ */
/*   REAL8 band;                       /\**< the frequency band *\/ */
/*   REAL8 mm;                         /\**< the mismatch *\/ */
/*   UINT4 nseg; */
/*   UINT4 ndim; */
/*   REAL8 Bayes; */
/*   REAL8 Bayes_fixed; */
/*   REAL8 min_amp; */
/*   REAL8 max_amp; */
/*   REAL8 sig_amp; */
/*   REAL8 start_amp; */
/*   REAL8 delta_amp; */
/*   UINT4 length_amp; */
/*   CHAR prior_amp[STRINGLENGTH]; */
/*   REAL8Vector *logposterior_fixed_amp; */
/*   REAL8Vector *Bayes_perseg; */
/*   UINT4Vector *Bayes_start; */
/*   UINT4Vector *Bayes_end; */
/*   CHAR **name; */
/*   REAL8 min[NBINMAX]; */
/*   REAL8 max[NBINMAX]; */
/*   REAL8 sig[NBINMAX]; */
/*   REAL8 start[NBINMAX]; */
/*   REAL8 delta[NBINMAX]; */
/*   UINT4 length[NBINMAX]; */
/*   CHAR **prior; */
/*   REAL8Vector *logposterior[NBINMAX]; */
/*   REAL8Vector *logposterior_fixed[NBINMAX]; */


/*   } */

  /* perform interpolation on each amplitude posterior */
 /*  for (i=0;i<resultsfiles->length;i++) { */

/*     UINT4 N = resultsfiles->file[i].length_amp; */
/*     gsl_interp_accel *acc = gsl_interp_accel_alloc();        /\* gsl interpolation structures *\/ */
/*     gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,N); */
/*     REAL8 *x = NULL; */
/*     REAL8 *y = NULL; */
    
/*     /\* allocate memory for the temporary data *\/ */
/*     if ((x = (REAL8 *)XLALCalloc(N,sizeof(REAL8))) == NULL) { */
/*       LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary gsl input.\n",fn); */
/*       XLAL_ERROR(fn,XLAL_ENOMEM); */
/*     } */
/*     if ((y = (REAL8 *)XLALCalloc(N,sizeof(REAL8))) == NULL) { */
/*       LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary gsl input.\n",fn); */
/*       XLAL_ERROR(fn,XLAL_ENOMEM); */
/*     } */

/*     /\* fill in gsl interpolation input data *\/ */
/*     for (j=0;j<N;j++) { */
/*       x[j] = resultsfiles->file[i].start_amp + j*resultsfiles->file[i].delta_amp; */
/*       y[j] = resultsfiles->file[i].logposterior_fixed_amp->data[j]; */
/*       printf("x = %f y = %f\n",x[j],y[j]); */
/*     } */
/*     gsl_spline_init(spline,x,y,N); */

/*     /\* interpolate and add to result *\/ */
/*     for (j=0;j<(*combinedresult)->logposterior_fixed_amp->length;j++) { */
/*       REAL8 z = (*combinedresult)->start_amp + j*(*combinedresult)->delta_amp; */
/*       REAL8 temp1 = (*combinedresult)->logposterior_fixed_amp->data[j]; */
/*       REAL8 temp2 = gsl_spline_eval(spline,z,acc) + log(deltaband[i]) - log(totalband); */
/*       (*combinedresult)->logposterior_fixed_amp->data[j] = XLALLogSumExp(temp1,temp2); */
/*     } */

/*     /\* free memory *\/ */
/*     XLALFree(x); */
/*     XLALFree(y); */
/*     gsl_interp_accel_free(acc); */
/*     gsl_spline_free(spline); */

/*   } */

  /* output results to screen */
/*   for (j=0;j<(*combinedresult)->logposterior_fixed_amp->length;j++) { */
/*     REAL8 z = (*combinedresult)->start_amp + j*(*combinedresult)->delta_amp; */
/*     LogPrintf(LOG_DEBUG,"%s : combined amplitude posterior %f %f\n",fn,z,(*combinedresult)->logposterior_fixed_amp->data[j]); */
/*   } */

  /* loop over each search parameter */
  for (k=0;k<(*combinedresult)->ndim;k++) {
    
    UINT4 maxlen = 0;
    for (i=0;i<resultsfiles->length;i++) {
      if (resultsfiles->file[i].logposterior_fixed[k]->length>maxlen) maxlen = resultsfiles->file[i].logposterior_fixed[k]->length;
    }
    printf("maxlen = %d\n",maxlen);
    
    /* we define an over-resolved interpolation scheme with OVERRES more points than the max value */
    if (((*combinedresult)->logposterior[k] = XLALCreateREAL8Vector((UINT4)(0.5 + OVERRES*maxlen))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output amplitude posterior results.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    if (((*combinedresult)->logposterior_fixed[k] = XLALCreateREAL8Vector((UINT4)(0.5 + OVERRES*maxlen))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for the output amplitude posterior results.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /* initialise results */
    for (j=0;j<(*combinedresult)->logposterior[k]->length;j++) {
      (*combinedresult)->logposterior[k]->data[j] = -1e200;
      (*combinedresult)->logposterior_fixed[k]->data[j] = -1e200;
    }
    
    /* fill in some header info */
    (*combinedresult)->start[k] = resultsfiles->file[0].start[k];
    (*combinedresult)->length[k] = (*combinedresult)->logposterior[k]->length;
    (*combinedresult)->delta[k] = ((resultsfiles->file[0].length[k]-1)*resultsfiles->file[0].delta[k])/((*combinedresult)->length[k] - 1);
    
    /* perform interpolation on each file for this parameter posterior */
    for (i=0;i<resultsfiles->length;i++) {
      
      UINT4 N = resultsfiles->file[i].length[k];
      gsl_interp_accel *acc = gsl_interp_accel_alloc();        /* gsl interpolation structures */
      gsl_interp_accel *acc_fixed = gsl_interp_accel_alloc();
      gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,N);
      gsl_spline *spline_fixed = gsl_spline_alloc(gsl_interp_cspline,N);
      REAL8 *x = NULL;
      REAL8 *y = NULL;
      REAL8 *y_fixed = NULL;
      
      /* allocate memory for the temporary data */
      if ((x = (REAL8 *)XLALCalloc(N,sizeof(REAL8))) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary gsl input.\n",fn);
	XLAL_ERROR(fn,XLAL_ENOMEM);
      }
      if ((y = (REAL8 *)XLALCalloc(N,sizeof(REAL8))) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary gsl input.\n",fn);
	XLAL_ERROR(fn,XLAL_ENOMEM);
      }
      if ((y_fixed = (REAL8 *)XLALCalloc(N,sizeof(REAL8))) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary gsl input.\n",fn);
	XLAL_ERROR(fn,XLAL_ENOMEM);
      }
      
      /* fill in gsl interpolation input data */
      for (j=0;j<N;j++) {
	x[j] = resultsfiles->file[i].start[k] + j*resultsfiles->file[i].delta[k];
	y[j] = resultsfiles->file[i].logposterior[k]->data[j];
	y_fixed[j] = resultsfiles->file[i].logposterior_fixed[k]->data[j];
	printf("x = %f y = %f y_fixed = %f\n",x[j],y[j],y_fixed[j]);
      }
      gsl_spline_init(spline,x,y,N);
      gsl_spline_init(spline_fixed,x,y_fixed,N);
      
      /* interpolate and add to result */
      for (j=0;j<(*combinedresult)->logposterior_fixed[k]->length;j++) {
	REAL8 z = (*combinedresult)->start[k] + j*(*combinedresult)->delta[k];
	REAL8 temp1 = (*combinedresult)->logposterior[k]->data[j];
	REAL8 temp2 = gsl_spline_eval(spline,z,acc) + log(deltaband[i]) - log(totalband);
	REAL8 temp3 = (*combinedresult)->logposterior_fixed[k]->data[j];
	REAL8 temp4 = gsl_spline_eval(spline_fixed,z,acc_fixed) + log(deltaband[i]) - log(totalband);
	(*combinedresult)->logposterior[k]->data[j] = XLALLogSumExp(temp1,temp2);
	(*combinedresult)->logposterior_fixed[k]->data[j] = XLALLogSumExp(temp3,temp4);
      }
      
      /* free memory */
      XLALFree(x);
      XLALFree(y);
      XLALFree(y_fixed);
      gsl_interp_accel_free(acc);
      gsl_interp_accel_free(acc_fixed);
      gsl_spline_free(spline);
      gsl_spline_free(spline_fixed);
      
    }
   
    /* output results to screen */
    for (j=0;j<(*combinedresult)->logposterior[k]->length;j++) {
      REAL8 z = (*combinedresult)->start[k] + j*(*combinedresult)->delta[k];
      LogPrintf(LOG_DEBUG,"%s : combined posterior %f %f %f\n",fn,z,(*combinedresult)->logposterior[k]->data[j],(*combinedresult)->logposterior_fixed[k]->data[j]);
    }
    exit(0);
  }

  /* free memory*/
  XLALFree(deltaband);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
 
}


/* function to compute the log of the sum of the arguments of two logged quantities
 *
 * Eg. input log(x) and log(y) -> output log(x+y)
 *
 * If you do this by exponentiating first, then summing and then logging again you
 * can easily gat overflow errors.  We use a trick to avoid this.
 */
REAL8 XLALLogSumExp(REAL8 logx,      /**< [in] the log of x */  
		    REAL8 logy       /**< [in] the log of y */
		    )
{

  /* compute log(x + y) = logmax + log(1.0 + exp(logmin - logmax)) */
  /* this way the argument to the exponential is always negative */
  /* return logmax + log(1.0 + exp(logmin - logmax)); */
  if (logy>=logx) {
    return logy + log(1.0 + exp(logx - logy));
  }
  else { 
    return logx + log(1.0 + exp(logy - logx));
  }
  
}
