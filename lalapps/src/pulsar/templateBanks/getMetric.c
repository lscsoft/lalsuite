/*
*  Copyright (C) 2007 Reinhard Prix
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
 * \file getMetric.c
 * \ingroup pulsarApps
 * \author Reinhard Prix
 * \brief
 * Standalone code to calculated the metric in a given parameter-space point
 * using the MetricWrapper.
 */

/* ---------- includes ---------- */
#include <math.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/DopplerScan.h>

/* Error codes and messages */
#define GETMETRIC_ENORM 	0
#define GETMETRIC_ESUB  	1
#define GETMETRIC_EARG  	2
#define GETMETRIC_EBAD  	3
#define GETMETRIC_EFILE 	4
#define GETMETRIC_ENOARG 	5
#define GETMETRIC_EINPUT 	6
#define GETMETRIC_EDIM 		7


#define GETMETRIC_MSGENORM 	"Normal exit"
#define GETMETRIC_MSGESUB  	"Subroutine failed"
#define GETMETRIC_MSGEARG  	"Error parsing arguments"
#define GETMETRIC_MSGEBAD  	"Bad argument values"
#define GETMETRIC_MSGEFILE 	"File IO error"
#define GETMETRIC_MSGENOARG 	"Missing argument"
#define GETMETRIC_MSGEINPUT 	"Invalid user input"
#define GETMETRIC_MSGEDIM 	"Invalid metric dimension"

/*---------- local defines ---------- */
#define TRUE (1==1)
#define FALSE (1==0)

/* ---------- some local types ---------- */

/* User variables */
typedef struct {
  BOOLEAN help;
  CHAR* IFO;		/**< name of detector (LHO, LLO, GEO, VIRGO,..)*/
  REAL8 Alpha;		/**< skyposition Alpha: radians, equatorial coords. */
  REAL8 Delta;		/**< skyposition Delta: radians, equatorial coords. */
  REAL8 Freq;		/**< target-frequency */
  REAL8 f1dot;		/**< target 1. spindown-value df/dt */

  CHAR *ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *ephemSun;	/**< Sun ephemeris file to use */

  INT4  metricType;	/**< Metric function to use: ptole_analytic, ptole-numeric, ... */
  REAL8 startTime;	/**< GPS start time of observation */
  REAL8 duration;	/**< length of observation in seconds */
  BOOLEAN projectMetric;/**< project out frequency-component of metric */
} UserInput;

typedef struct {
  LALDetector *site;     	/**< detector of data to be searched */
  EphemerisData *ephemeris;/**< ephemeris data (from LALInitBarycenter()) */
  LIGOTimeGPS startTimeGPS;	/**< starttime of observation */
} ConfigVariables;

/*---------- empty structs for initializations ----------*/
static const ConfigVariables empty_ConfigVariables;
static const PtoleMetricIn empty_metricpar;
static const UserInput empty_UserInput;

/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *, UserInput *uvar, int argc, char *argv[] );
void initGeneral (LALStatus *status, ConfigVariables *cfg, const UserInput *uvar);
void printPulsarMetric(LALStatus *, const UserInput *uvar, const ConfigVariables *config);
void printFlatPulsarMetric (LALStatus *, const UserInput *uvar, const ConfigVariables *config);
void printMetric (LALStatus *status, const REAL8Vector *metric );

extern int vrbflg;

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = blank_status;
  ConfigVariables config = empty_ConfigVariables;
  UserInput uvar = empty_UserInput;

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (initUserVars (&status, &uvar, argc, argv), &status);	  

  if (uvar.help) 	/* help requested: we're done */
    exit (0);

  /* some general setup of the code */
  LAL_CALL (initGeneral (&status, &config, &uvar), &status);

  if ( uvar.metricType != 4) {
    LAL_CALL (printPulsarMetric(&status, &uvar, &config), &status);
  } else {
    LAL_CALL (printFlatPulsarMetric(&status, &uvar, &config), &status);
  }

  /* ----- Free memory ----- */
  LAL_CALL ( LALDestroyUserVars (&status), &status);


  XLALDestroyEphemerisData ( config.ephemeris );   /* Free ephemeris data */
  LALFree ( config.site );

  LALCheckMemoryLeaks(); 

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/** register all our "user-variables" */
void
initUserVars (LALStatus *status, UserInput *uvar, int argc, char *argv[])
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* ----- set defaults ----- */
  uvar->help = FALSE;
  uvar->IFO = NULL;

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->f1dot = 0.0;
  uvar->metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;

  uvar->startTime = 714180733;
  uvar->projectMetric = FALSE;

  /* ----- register all our user-variable ----- */

  XLALregBOOLUserStruct (  help,		'h', UVAR_HELP,     "Print this help/usage message" );
  XLALregSTRINGUserStruct (IFO,			'I', UVAR_REQUIRED, "Detector: H1, H2, L1, G1, ..." );

  XLALregREALUserStruct (  Alpha,		'a', UVAR_REQUIRED,  "skyposition Alpha in radians, equatorial coords." );
  XLALregREALUserStruct (  Delta, 		'd', UVAR_REQUIRED,  "skyposition Delta in radians, equatorial coords." );
  XLALregREALUserStruct (  Freq, 		'f', UVAR_REQUIRED,  "target frequency" );
  XLALregREALUserStruct (  f1dot, 		's', UVAR_OPTIONAL,  "first spindown-value df/dt" );
  XLALregINTUserStruct  (  metricType,   	'M', UVAR_OPTIONAL,  "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact, 4=FLAT[old]" );
  XLALregBOOLUserStruct (  projectMetric,   	 0, UVAR_OPTIONAL,  "Project metric onto frequency-surface" );
  XLALregREALUserStruct (  startTime,    	't', UVAR_OPTIONAL,  "GPS start time of observation" );
  XLALregREALUserStruct (  duration,		'T', UVAR_REQUIRED,  "Duration of observation in seconds" );

  XLALregSTRINGUserStruct (ephemEarth,   	 0,  UVAR_OPTIONAL,     "Earth ephemeris file to use");
  XLALregSTRINGUserStruct (ephemSun,     	 0,  UVAR_OPTIONAL,     "Sun ephemeris file to use");

  /* read cmdline & cfgfile  */
  TRY (LALUserVarReadAllInput (status->statusPtr, argc, argv), status);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

/*----------------------------------------------------------------------*/
/**
 * do some general initializations,
 * e.g. load ephemeris-files (if required), setup detector etc
 */
void
initGeneral (LALStatus *status, ConfigVariables *cfg, const UserInput *uvar)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  XLALGPSSetREAL8(&(cfg->startTimeGPS),  uvar->startTime);

  /* ---------- init ephemeris if needed ---------- */
  if ( uvar->metricType ==  LAL_PMETRIC_COH_EPHEM )
    {
      /* Init ephemerides */
      if ( (cfg->ephemeris = XLALInitBarycenter ( uvar->ephemEarth, uvar->ephemSun )) == NULL ) {
        ABORT (status, GETMETRIC_EFILE, GETMETRIC_MSGEFILE);
      }
    } /* end: init ephemeris data */

  /* ---------- initialize detector ---------- */
  if ( (cfg->site = XLALGetSiteInfo ( uvar->IFO )) == NULL ) {
    ABORT (status, GETMETRIC_EINPUT, GETMETRIC_MSGEINPUT);
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* initGeneral() */

/**
 * Call LALPulsarMetric(), which is in {f, alpha, delta, f1dot, ...} coordinates
 * and has non-constant coefficients.
 */
void
printPulsarMetric(LALStatus *status, const UserInput *uvar, const ConfigVariables *config)
{
  REAL8Vector *metric = NULL;
  PtoleMetricIn metricpar = empty_metricpar;
  UINT4 dim, a;
  static LALDetector this_detector;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  metricpar.position.system = COORDINATESYSTEM_EQUATORIAL;
  metricpar.position.longitude = uvar->Alpha;
  metricpar.position.latitude = uvar->Delta;
  
  TRY ( LALSCreateVector (status->statusPtr, &(metricpar.spindown), 1), status);
  metricpar.spindown->data[0] = uvar->f1dot / uvar->Freq; /* f1 = df/dt / f0 !!*/
  
  metricpar.epoch = config->startTimeGPS;
  metricpar.duration = (REAL4) uvar->duration;
  metricpar.maxFreq = uvar->Freq;
  
  this_detector = *(config->site);
  metricpar.site = &this_detector;
  metricpar.ephemeris = config->ephemeris;	/* needed for ephemeris-metrics */
  metricpar.metricType = uvar->metricType;
  
  TRY ( LALPulsarMetric(status->statusPtr, &metric, &metricpar), status);
  TRY ( LALSDestroyVector(status->statusPtr, &(metricpar.spindown)), status);
  
  if (uvar->projectMetric) {
    TRY ( LALProjectMetric( status->statusPtr, metric, 0 ), status);
  }
  
  /* convert f1 to f1dot: f1 = f1dot / f */
  dim = XLALFindMetricDim (metric );
  for ( a = 0; a < dim; a++ ) 
    {
      metric->data [ PMETRIC_INDEX(3,a) ] /= uvar->Freq;
      if ( a == 3 )
	metric->data [ PMETRIC_INDEX(3,a) ] /= uvar->Freq;
    }

  if ( lalDebugLevel )
    printf ("\n%%%% %s Metric: (f, alpha, delta, f1dot) \n",
	    uvar->projectMetric ? "Projected":"Unprojected");
  printf ("\n g_ij = \\\n");
  TRY ( printMetric(status->statusPtr, metric ), status);

  TRY( LALDDestroyVector (status->statusPtr, &metric), status);

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* printPulsarMetric() */


void
printFlatPulsarMetric (LALStatus *status, const UserInput *uvar, const ConfigVariables *config)
{
  REAL8Vector *metric = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* calculate the flat metric */
  TRY ( LALFlatPulsarMetric(status->statusPtr, 
			    &metric, config->startTimeGPS, uvar->duration, config->site), status);

  /* formatted output of metric coefficients */
  printf ("\n%%%% Flat Metric in dimensionless variables (kappaX, kappaY, w0, w1, w2, ...) \n");
  printf ("\ng_ij = \\\n");
  TRY ( printMetric(status->statusPtr, metric), status );

  /* free memory */
  TRY( LALDDestroyVector (status->statusPtr, &metric), status);

  DETATCHSTATUSPTR ( status );
  RETURN(status);

} /* printFlatPulsarMetric() */

/* Printout metric in 'standard' format of [freq, sky1, sky2, f1dot, f2dot, ... ]
 * while the input is in 'internal' format of [sky1, sky2, freq, f1dot, ... ]
 */
void
printMetric (LALStatus *status, const REAL8Vector *metric )
{
  UINT4 row, col, dim;

  INITSTATUS(status);

  ASSERT ( metric, status, GETMETRIC_EBAD, GETMETRIC_MSGEBAD);
  ASSERT ( metric->length, status, GETMETRIC_EBAD, GETMETRIC_MSGEBAD);

  dim = XLALFindMetricDim ( metric );

  /* print metric components (in octave/matlab format) */
  printf  ("[ ");
  for ( row = 0; row < dim; row ++ )
    {
      if ( row > 0 ) printf ("  ");
      for ( col = 0; col < dim; col ++ )
	{
	  printf (" %.16g", metric->data[ PMETRIC_INDEX(row,col) ] );
          if ( col < dim - 1 )
            printf (", ");
	}
      if ( row < dim - 1 )
	printf (";\n");
      else printf("];\n");
    }

  RETURN(status);

} /* printMetric() */
