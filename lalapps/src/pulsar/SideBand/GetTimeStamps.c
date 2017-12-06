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
#include <lal/SFTfileIO.h>

/*---------- DEFINES ----------*/
#define MAXUINT4 2147483647

/*----- Error-codes -----*/
#define GETTIMESTAMPSC_ENULL 		1
#define GETTIMESTAMPSC_ESYS            	2
#define GETTIMESTAMPSC_EINPUT          	3
#define GETTIMESTAMPSC_EXLAL	        4

#define GETTIMESTAMPSC_MSGENULL        	"Arguments contained an unexpected null pointer"
#define GETTIMESTAMPSC_MSGESYS		"System call failed (probably file IO)"
#define GETTIMESTAMPSC_MSGEINPUT       	"Invalid input"
#define GETTIMESTAMPSC_MSGEXLAL		"XLALFunction-call failed"


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
CHAR* uvar_stampsfile; 
CHAR* uvar_sftdir;
INT4 uvar_tstart;
REAL8 uvar_duration;
BOOLEAN uvar_help;
CHAR* uvar_IFO;

REAL8 jumptemp;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *);

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
  FILE *fp = NULL;
  LALStatus status = blank_status;	/* initialize status */
  SFTConstraints XLAL_INIT_DECL(constraints);
  SFTCatalog *catalog = NULL;
  LIGOTimeGPS XLAL_INIT_DECL(start);
  LIGOTimeGPS XLAL_INIT_DECL(end);
  INT4 i;

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
  
  /* get an sft catalog */
  start.gpsSeconds = (INT4)uvar_tstart;
  end.gpsSeconds = (INT4)uvar_tstart + (INT4)uvar_duration;
  constraints.minStartTime = &start;
  constraints.maxStartTime = &end;
  LALSFTdataFind(&status,&catalog,uvar_sftdir,&constraints);
  
  /* output timestamps to file */
  if ((fp = fopen(uvar_stampsfile,"w"))==NULL) {
    XLALPrintError("\nError opening file '%s' for writing..\n\n",uvar_stampsfile);
    return (GETTIMESTAMPSC_ESYS);
  }
  
  for (i=0;i<(INT4)catalog->length;i++) {
    fprintf(fp,"%d %d\n",catalog->data[i].header.epoch.gpsSeconds,catalog->data[i].header.epoch.gpsNanoSeconds);
  }
  fclose(fp);

  if (lalDebugLevel) printf ("\nFinished outputting to file.\n");
  
  /* Free config-Variables and userInput stuff */
  LALDestroyUserVars(&status);
  
  /* clear memory */
  LALDestroySFTCatalog(&status,&catalog);
  
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
  uvar_stampsfile = NULL;
  uvar_sftdir = NULL;
  uvar_tstart = 0;
  uvar_duration = 0;
  uvar_IFO = NULL;

  /* register all our user-variables */
  LALregBOOLUserVar(status, 	help, 		                      'h', UVAR_HELP,     "Print this message"); 
  LALregSTRINGUserVar(status, 	stampsfile, 	                      'A', UVAR_OPTIONAL, "Name of timestamps file (default = NULL)");
  LALregSTRINGUserVar(status, 	sftdir, 	                      'B', UVAR_OPTIONAL, "Name of sft directory (default = NULL)");
  LALregINTUserVar(status,      tstart,                               'C', UVAR_REQUIRED, "GPS start time (default = 0)");
  LALregREALUserVar(status, 	duration,       	              'D', UVAR_REQUIRED, "Duration in seconds (default = 0)");
  LALregSTRINGUserVar(status, 	IFO,     	                      'E', UVAR_OPTIONAL, "Name of IFO (default = NULL)");

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */

