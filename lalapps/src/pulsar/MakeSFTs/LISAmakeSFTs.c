/*
 * Copyright (C) 2006 Reinhard Prix
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
 * \author Reinhard Prix and John Whelan
 * \date 2006
 * \file 
 * \brief Read in MLDC timeseries-files and produce SFTs (v2) for them
 *
 * $Id$
 *
 */

/* ---------- includes ---------- */
#include <unistd.h>


#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdio.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/LISAspecifics.h>
#include <lal/LogPrintf.h>

/* lisaXML stuff */
#include "readxml.h"

RCSID ("$Id$");

/** \name Error codes */
/*@{*/
#define LISAMAKESFTS_ENORM 	0
#define LISAMAKESFTS_EINPUT  	1
#define LISAMAKESFTS_EMEM	2
#define LISAMAKESFTS_EFILE	3
#define LISAMAKESFTS_ENULL	4
#define LISAMAKESFTS_ENONULL	5

#define LISAMAKESFTS_MSGENORM 	"Normal exit"
#define LISAMAKESFTS_MSGEINPUT  "Bad argument values"
#define LISAMAKESFTS_MSGEMEM	"Out of memory"
#define LISAMAKESFTS_MSGEFILE	"File input/output error"
#define LISAMAKESFTS_MSGENULL	"Input contained illegal NULL"
#define LISAMAKESFTS_MSGENONULL	"Output container non-NULL"

/*@}*/

/*---------- DEFINES ----------*/
#define MAX_FILENAME_LEN   5012

#define TRUE    (1==1)
#define FALSE   (1==0)

/*---------- pseudo-LISA parameters ----------*/
#define LISA_ARM_LENGTH_SECONDS 16.6782

/*----- Macros ----- */
/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
static const LALUnit empty_LALUnit;
/*---------- Global variables ----------*/

/* User variables */
BOOLEAN uvar_help;
BOOLEAN uvar_lisasim;
CHAR *uvar_extraComment;
CHAR *uvar_outputDir;
CHAR *uvar_inputXML;
REAL8 uvar_Tsft;
CHAR *uvar_miscField;

/*---------- internal prototypes ----------*/
void initUserVars (LALStatus *status);
void ConvertLISAtimeseries2LAL ( LALStatus *status, MultiREAL4TimeSeries **lalTimeSeries, const TimeSeries *lisaTimeSeries );

CHAR *assembleDescription ( const CHAR *name, const CHAR *miscField );

/*==================== FUNCTION DEFINITIONS ====================*/

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */
  CHAR *add_comment = NULL;
  TimeSeries *lisaTimeSeries;		/* lisaXML timeseries-type */
  MultiREAL4TimeSeries *multiTs = NULL;	/* LAL-equivalent: hold 3 timeseries (X(t), Y(t), Z(t)) */
  SFTParams sftParams;
  UINT4 ifo, sidx, fidx;
  COMPLEX8FrequencySeries *sft = NULL;
  REAL4 fourpifL;
  COMPLEX8 ztmp;

  lalDebugLevel = 0;

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */
  
  /* set debug level */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);

  /* register all user-variables */
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* -- unfortunately getTDIdata() only works if xml-file is in local directory! ==> we need to cd there */
  {
    CHAR basedir[MAX_FILENAME_LEN], curdir[MAX_FILENAME_LEN];
    CHAR *tmp, *basename;
    /* remember current directory */
    if ( getcwd ( curdir, MAX_FILENAME_LEN ) == NULL ) {
      LogPrintf ( LOG_CRITICAL, "Current directory path '%s' exceeds maximal length %d\n", curdir, MAX_FILENAME_LEN );
      return LISAMAKESFTS_EINPUT;      
    }
    /* copy input-XML directory and split in basedir + basename */
    if ( strlen ( uvar_inputXML ) >= MAX_FILENAME_LEN ) {
      LogPrintf ( LOG_CRITICAL, "Input path '%s' exceeds maximal length %d \n", uvar_inputXML, MAX_FILENAME_LEN );
      return LISAMAKESFTS_EINPUT;      
    }
    strcpy ( basedir, uvar_inputXML );
    if ( (tmp = strrchr( basedir, '/')) == NULL ) {
      basedir[0] = '.'; basedir[1] = 0;
      basename = uvar_inputXML;
    }
    else {
      *tmp = 0;	/* just terminate string at last '/' */
      basename = tmp + 1;
    }
    /* change into input-XML directory for reading the time-series */
    if(chdir( basedir ) != 0)
      {
	LogPrintf (LOG_CRITICAL,  "Unable to change directory to xml-directory '%s'\n", basedir);
	return LISAMAKESFTS_EINPUT;
      }
    /* load xml-file and corresponding binary-data into lisaXML-type 'TimeSeries' */
    if ( (lisaTimeSeries = getTDIdata(basename)) == NULL ) {
      fprintf (stderr, "\nlisaXML::getTDIdata() failed for file '%s'\n\n",  uvar_inputXML );
      return LISAMAKESFTS_EFILE;
    }
    /* return to original directory */
    if(chdir( curdir ) != 0)
      {
	LogPrintf (LOG_CRITICAL,  "Unable to return to start-directory '%s'\n", curdir);
	return LISAMAKESFTS_EINPUT;
      }
  } /* chdir to xml-directory */

  /* convert lisaXML::TimeSeries -> LAL::REAL4TimeSeries */
  LAL_CALL ( ConvertLISAtimeseries2LAL ( &status, &multiTs, lisaTimeSeries ), &status );

  { /* build up full comment-string to be added to SFTs: 'generated by LISAmakeSFTs + RCSID + cmdline + user extraComment' */
    UINT4 len;
    CHAR *logstr = NULL;
    LAL_CALL ( LALUserVarGetLog ( &status, &logstr,  UVAR_LOGFMT_CMDLINE ), &status );
    len = 512 + strlen ( logstr );
    if ( ( add_comment = LALCalloc ( 1, len ) ) == NULL ) 
      return LISAMAKESFTS_EMEM;
    sprintf ( add_comment, "Generated by $Id$:\n%s", logstr );
    LALFree ( logstr );
  } /* assemble comment-string */

  /* convert timeseries -> SFTs and write them to disk */
  sftParams.Tsft = uvar_Tsft;
  sftParams.timestamps = NULL;
  sftParams.noiseSFTs = NULL;
  sftParams.make_v2SFTs = 1;

  for ( ifo=0; ifo < multiTs->length; ifo ++ )
    {
      SFTVector *SFTvect = NULL;
      CHAR *desc;
      if ( (desc = assembleDescription ( multiTs->data[ifo]->name, uvar_miscField )) == NULL )
	return -1;

      LAL_CALL ( LALSignalToSFTs (&status, &SFTvect, multiTs->data[ifo], &sftParams ), &status );
      /* Calibrate into "strain" using long-wavelength approx */
      for ( sidx=0; sidx < SFTvect->length; sidx++ )
	{
	  sft = &(SFTvect->data[sidx]);
	  for ( fidx=0; fidx < sft->data->length; fidx++ )
	    {
	      fourpifL = (2.0*LAL_TWOPI*LISA_ARM_LENGTH_SECONDS)
		* (sft->f0 + fidx * sft->deltaF);
	      if (uvar_lisasim)
		{
		  /* divide by    - 4 pi i f L */
		  ztmp = sft->data->data[fidx];
		  sft->data->data[fidx].re = - ztmp.im / fourpifL;
		  sft->data->data[fidx].im = ztmp.re / fourpifL;
		}
	      else
		{
		  sft->data->data[fidx].re /= (fourpifL * fourpifL);
		  sft->data->data[fidx].im /= (fourpifL * fourpifL);
		}
	    }
	}

      LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
      LALFree ( desc );
      LAL_CALL ( LALDestroySFTVector ( &status, &SFTvect ), &status );
    } /* for ifo < length */


  /* free memory */
  LALFree ( add_comment );
  for ( ifo = 0; ifo < multiTs->length ; ifo ++ )
    XLALDestroyREAL4TimeSeries ( multiTs->data[ifo] );
  LALFree ( multiTs->data );
  LALFree ( multiTs );
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALCheckMemoryLeaks(); 

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *status)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set defaults */
#define DEFAULT_OUTDIR  "./"
  uvar_outputDir = LALCalloc ( 1, strlen (DEFAULT_OUTDIR) + 1);
  strcpy ( uvar_outputDir, DEFAULT_OUTDIR );

  uvar_extraComment = NULL;
  uvar_miscField = NULL;
  uvar_Tsft = 604800;	

  /* now register all our user-variable */
  LALregSTRINGUserVar(status, inputXML,		'i', UVAR_REQUIRED, "XML file describing the LISA timeseries data");
  LALregREALUserVar(status,   Tsft,		'T', UVAR_OPTIONAL, "Length of SFTs to produce (in seconds)");
  LALregSTRINGUserVar(status, outputDir,	'o', UVAR_OPTIONAL, "Output directory for SFTs");
  LALregSTRINGUserVar(status, extraComment,	'C', UVAR_OPTIONAL, "Additional comment to be added to output-SFTs");
  LALregSTRINGUserVar(status, miscField,	'm', UVAR_OPTIONAL, "User-specifiable portion of the SFT-filename ('misc' field)");
  
  LALregBOOLUserVar(status,   lisasim,		's', UVAR_OPTIONAL, "TDI data are from LISA Simulator");

  LALregBOOLUserVar(status,   help,		'h', UVAR_HELP,     "Print this help/usage message");
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

/** Convert a lisaXML 'TimeSeries' into a LAL MultiREAL4TimeSeries */
void
ConvertLISAtimeseries2LAL ( LALStatus *status, MultiREAL4TimeSeries **lalTs, const TimeSeries *lisaTs )
{
  UINT4 nIFOs, i;
  MultiREAL4TimeSeries *ret = NULL;

  INITSTATUS( status, "ConvertLISAtimeseries2LAL", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( lalTs, status, LISAMAKESFTS_ENULL, LISAMAKESFTS_MSGENULL );
  ASSERT ( *lalTs == NULL, status, LISAMAKESFTS_ENONULL, LISAMAKESFTS_MSGENONULL );
  ASSERT ( lisaTs, status, LISAMAKESFTS_ENULL, LISAMAKESFTS_MSGENULL );

  nIFOs = (UINT4) lisaTs->Records;
  if ( nIFOs == 0 || nIFOs == 1 ) {
    ABORT ( status,  LISAMAKESFTS_EINPUT,  LISAMAKESFTS_MSGEINPUT );
  }
  nIFOs --; /* first lisa-timeseries is just timesteps */

  /* alloctate LAL container */
  if ( (ret = LALCalloc ( 1, sizeof ( MultiREAL4TimeSeries ) )) == NULL ) {
    ABORT ( status, LISAMAKESFTS_EMEM, LISAMAKESFTS_MSGEMEM );
  }
  if ( (ret->data = LALCalloc ( nIFOs, sizeof ( REAL4TimeSeries* ) )) == NULL ) {
    LALFree ( ret );
    ABORT ( status, LISAMAKESFTS_EMEM, LISAMAKESFTS_MSGEMEM );    
  }
  ret->length = nIFOs;

  /* allocate and convert individual timeseries X, Y, Z */
  for ( i=0; i < nIFOs; i ++ )
    {
      UINT4 l;
      CHAR name[LALNameLength];
      LIGOTimeGPS epoch = { 0, 0 };
      REAL8 f0 = 0;	/* no heterodyning */
      LALUnit units = empty_LALUnit;

      DataColumn *thisTs = lisaTs->Data[i+1];	/* skip first column: timesteps */
      REAL8 deltaT = thisTs->Cadence;
      size_t length = (size_t) thisTs->Length;

      /* Naming-convention: channel = {Z1, Z2, Z3} + ts-name + filename */
      LALSnprintf ( name, LALNameLength, "Z%d:%s_%s", i+1, thisTs->Name, lisaTs->FileName );
      name[LALNameLength-1] = 0; /* close string if it was truncated */

      epoch.gpsSeconds = thisTs->TimeOffset + LISA_TIME_ORIGIN;	/* offset for convenience of GPS-time ranges */

      if ( ( ret->data[i] = XLALCreateREAL4TimeSeries ( name, &epoch, f0, deltaT, &units, length )) == NULL )
	goto failed;

      /* now cast + copy all the data */
      for ( l=0; l < length; l ++ )
	ret->data[i]->data->data[l] = (REAL4) thisTs->data[l];

    } /* for i < nIFOs */

  /* ok: return final multiTimeSeries */
  (*lalTs) = ret;

  DETATCHSTATUSPTR (status);
  RETURN (status);

 failed:
  /* free all memory allocated */
  for ( i=0; i < nIFOs; i ++ )
    if ( ret->data[i] ) XLALDestroyREAL4TimeSeries ( ret->data[i] );
  if ( ret->data ) LALFree ( ret->data );
  if ( ret ) LALFree ( ret );

  ABORT ( status, LISAMAKESFTS_EMEM, LISAMAKESFTS_MSGEMEM );    

} /* ConvertLISAtimeseries2LAL() */


/* create a valid string for the 'description'-field in the SFT filename 
 * [according to SFTv2/Frame naming convention]
 */
CHAR *
assembleDescription ( const CHAR *name, const CHAR *miscField )
{
  CHAR *desc, *ptr;
  UINT4 len;
  const CHAR illegals[] = "-. ";

  if ( !name )
    return NULL;
  len = strlen ( name ) + 1;

  if ( miscField ) 
    len += strlen ( miscField ) + 1;

  if ( (desc = LALMalloc ( len * sizeof(CHAR) )) == NULL )
    return NULL;

  if ( (ptr = strchr( name, ':')) == NULL )
    return NULL;

  strcpy ( desc, ptr + 1);
  if ( uvar_miscField ) {
    strcat ( desc, "_" );
    strcat ( desc, uvar_miscField );
  }

  /* Now go through and replace all illegal characters ['-', '.', ' '] by '_' */
  ptr = desc;
  while ( *ptr != 0 )
    {
      if ( strchr ( illegals, *ptr ) != NULL )
	*ptr = '_';
      ptr ++;
    }

  return desc;

} /* assembleDescription() */
