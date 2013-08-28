/*
 * Copyright (C) 2007 Deepak Khurana
 * Copyright (C) 2006, 2007, 2008 Reinhard Prix, John T Whelan
 *
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
 * \author Reinhard Prix, John T Whelan, Deepak Khurana
 * \date 2006-2008
 * \file
 * \ingroup pulsarApps
 * \brief Read in MLDC timeseries-files and produce SFTs (v2) for them
 */

/* ---------- includes ---------- */
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdio.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/LISAspecifics.h>
#include <lal/LogPrintf.h>
#include <lal/Date.h>

#include <lalapps.h>

/* lisaXML stuff */
#include "readxml.h"

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

#define LAL_SQRT1_3   0.5773502691896257645091487805019575  /**< 1/sqrt(3) */

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
BOOLEAN uvar_raCorrections;
BOOLEAN uvar_makeX;
BOOLEAN uvar_makeY;
BOOLEAN uvar_makeZ;
BOOLEAN uvar_makeYminusZ;
BOOLEAN uvar_makeZminusX;
BOOLEAN uvar_makeXminusY;
BOOLEAN uvar_makeA;
BOOLEAN uvar_makeE;
BOOLEAN uvar_makeT;
CHAR *uvar_extraComment;
CHAR *uvar_outputDir;
CHAR *uvar_inputXML;
REAL8 uvar_Tsft;
CHAR *uvar_miscField;

/*---------- internal prototypes ----------*/
void initUserVars (LALStatus *status);
void ConvertLISAtimeseries2LAL ( LALStatus *status, MultiREAL4TimeSeries **lalTimeSeries, const TimeSeries *lisaTimeSeries );
CHAR *assembleDescription ( const CHAR *name, const CHAR *miscField );

void compute_R_f ( COMPLEX8 *R_f, REAL8 Freq, BOOLEAN useRAA, BOOLEAN isLISAsim );

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
  COMPLEX8 R_f,h_f;
  SFTVector **SFTVectList;
  SFTVector *SFTvect;
  BOOLEAN writeTDI[3];


  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (initUserVars (&status), &status);

  /* read cmdline & cfgfile  */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);

  /* make life easier in the loop */
  writeTDI[0] = uvar_makeX;
  writeTDI[1] = uvar_makeY;
  writeTDI[2] = uvar_makeZ;

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* ----- make sure output directory exists ---------- */
  {
    int ret;
    ret = mkdir ( uvar_outputDir, 0777);
    if ( (ret == -1) && ( errno != EEXIST ) )
      {
	int errsv = errno;
	LogPrintf (LOG_CRITICAL, "Failed to create directory '%s': %s\n", uvar_outputDir, strerror(errsv) );
	return -1;
      }
  }

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

  { /* build up full comment-string to be added to SFTs: 'generated by LISAmakeSFTs + VCSInfo + cmdline + user extraComment' */
    UINT4 len;
    CHAR *logstr = NULL;
    CHAR *VCSInfoString;
    LAL_CALL ( LALUserVarGetLog ( &status, &logstr,  UVAR_LOGFMT_CMDLINE ), &status );

    if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
      XLALPrintError("XLALGetVersionString(0) failed.\n");
      exit(1);
    }

    len = 1 + strlen( VCSInfoString ) + strlen ( logstr );
    if ( ( add_comment = LALCalloc ( 1, len ) ) == NULL )
      return LISAMAKESFTS_EMEM;
    snprintf ( add_comment, len, "%s\n%s", VCSInfoString, logstr );
    LALFree ( VCSInfoString );
    LALFree ( logstr );
  } /* assemble comment-string */

  /* convert timeseries -> SFTs and write them to disk */
  sftParams.Tsft = uvar_Tsft;
  sftParams.timestamps = NULL;
  sftParams.noiseSFTs = NULL;
  sftParams.window = NULL;

  if ( ( SFTVectList = LALCalloc( multiTs->length, sizeof(*SFTVectList) ) ) == NULL )
      return LISAMAKESFTS_EMEM;

  for ( ifo=0; ifo < multiTs->length; ifo ++ )
    {
      CHAR *desc;
      if ( (desc = assembleDescription ( multiTs->data[ifo]->name, uvar_miscField )) == NULL )
	return -1;

      SFTvect = NULL;

      LAL_CALL ( LALSignalToSFTs (&status, &SFTvect, multiTs->data[ifo], &sftParams ), &status );
      SFTVectList[ifo] = SFTvect;
      /* Calibrate into "strain" using long-wavelength approx */
      for ( sidx=0; sidx < SFTvect->length; sidx++ )
	{
	  sft = &(SFTvect->data[sidx]);
	  for ( fidx=0; fidx < sft->data->length; fidx++ )
	    {
	      REAL8 Freq = sft->f0 + fidx * sft->deltaF;

	      /* get complex "transfer-function" */
	      compute_R_f(&R_f, Freq, uvar_raCorrections, uvar_lisasim );

	      /* multiply SFT-data h_f by "transfer function" R_f */
	      h_f = sft->data->data[fidx];	/* copy original value */
	      sft->data->data[fidx].realf_FIXME = crealf(h_f) * crealf(R_f) - cimagf(h_f) * cimagf(R_f);	/*(a+ib)(c+id)=(ac-bd)+i(ad+bc)*/
	      sft->data->data[fidx].imagf_FIXME = crealf(h_f) * cimagf(R_f) + cimagf(h_f) * crealf(R_f);

	    } /* for fidx < sft->data->length */

	} /* for sidx < SFTvect->length */

      if (writeTDI[ifo]) {
	LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
      }
      LALFree ( desc );
    } /* for ifo < length */

  if ( uvar_makeYminusZ ) {
    CHAR *desc;
    CHAR comboname[MAX_FILENAME_LEN];
    UINT4 ifo1 = 1;
    UINT4 ifo2 = 2;

    snprintf ( comboname, MAX_FILENAME_LEN, "Z4:{%s}-{%s}", multiTs->data[ifo1]->name, multiTs->data[ifo2]->name );
    if ( (desc = assembleDescription ( comboname, uvar_miscField )) == NULL )
      return -1;
    if ( multiTs->length != 3 )
      {
	LogPrintf (LOG_CRITICAL,  "Need 3 input time series to make Y-Z; got %d\n", multiTs->length);
	return LISAMAKESFTS_EINPUT;
      }
    SFTvect = NULL;
    LAL_CALL ( LALSubtractSFTVectors (&status, &SFTvect, SFTVectList[ifo1], SFTVectList[ifo2]), &status );
    for ( sidx=0; sidx < SFTvect->length; sidx++ )
      {
	SFTvect->data[sidx].name[0] = 'Z';
	SFTvect->data[sidx].name[1] = '4';
      }
    LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
    LAL_CALL ( LALDestroySFTVector ( &status, &(SFTvect) ), &status );
    LALFree ( desc );

  } /* if uvar_makeYminusZ */

  if ( uvar_makeZminusX ) {
    CHAR *desc;
    CHAR comboname[MAX_FILENAME_LEN];
    UINT4 ifo1 = 2;
    UINT4 ifo2 = 0;

    snprintf ( comboname, MAX_FILENAME_LEN, "Z5:{%s}-{%s}", multiTs->data[ifo1]->name, multiTs->data[ifo2]->name );
    if ( (desc = assembleDescription ( comboname, uvar_miscField )) == NULL )
      return -1;
    if ( multiTs->length != 3 )
      {
	LogPrintf (LOG_CRITICAL,  "Need 3 input time series to make Z-X; got %d\n", multiTs->length);
	return LISAMAKESFTS_EINPUT;
      }
    SFTvect = NULL;
    LAL_CALL ( LALSubtractSFTVectors (&status, &SFTvect, SFTVectList[ifo1], SFTVectList[ifo2]), &status );
    for ( sidx=0; sidx < SFTvect->length; sidx++ )
      {
	SFTvect->data[sidx].name[0] = 'Z';
	SFTvect->data[sidx].name[1] = '5';
      }
    LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
    LAL_CALL ( LALDestroySFTVector ( &status, &(SFTvect) ), &status );
    LALFree ( desc );

  } /* if uvar_makeZminusX */

  if ( uvar_makeXminusY ) {
    CHAR *desc;
    CHAR comboname[MAX_FILENAME_LEN];
    UINT4 ifo1 = 0;
    UINT4 ifo2 = 1;

    snprintf ( comboname, MAX_FILENAME_LEN, "Z6:{%s}-{%s}", multiTs->data[ifo1]->name, multiTs->data[ifo2]->name );
    if ( (desc = assembleDescription ( comboname, uvar_miscField )) == NULL )
      return -1;
    if ( multiTs->length != 3 )
      {
	LogPrintf (LOG_CRITICAL,  "Need 3 input time series to make X-Y; got %d\n", multiTs->length);
	return LISAMAKESFTS_EINPUT;
      }
    SFTvect = NULL;
    LAL_CALL ( LALSubtractSFTVectors (&status, &SFTvect, SFTVectList[ifo1], SFTVectList[ifo2]), &status );
    for ( sidx=0; sidx < SFTvect->length; sidx++ )
      {
	SFTvect->data[sidx].name[0] = 'Z';
	SFTvect->data[sidx].name[1] = '6';
      }
    LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
    LAL_CALL ( LALDestroySFTVector ( &status, &(SFTvect) ), &status );
    LALFree ( desc );

  } /* if uvar_makeXminusY */

  if ( uvar_makeA ) {
    CHAR *desc;
    CHAR comboname[MAX_FILENAME_LEN];
    COMPLEX16Vector *weights = NULL;

    snprintf ( comboname, MAX_FILENAME_LEN, "%s", multiTs->data[0]->name );

    comboname[0] = 'Z';
    comboname[1] = '7';
    comboname[3] = 'A';

    if ( (desc = assembleDescription ( comboname, uvar_miscField )) == NULL )
      return -1;
    if ( multiTs->length != 3 )
      {
	LogPrintf (LOG_CRITICAL,
		   "Need 3 input time series to make A; got %d\n",
		   multiTs->length);
	return LISAMAKESFTS_EINPUT;
      }
    SFTvect = NULL;

    LAL_CALL ( LALZCreateVector ( &status, &weights, 3 ) , &status );
    weights->data[0].real_FIXME = (  2.0 / 3.0 );
    weights->data[1].real_FIXME = ( -1.0 / 3.0 );
    weights->data[2].real_FIXME = ( -1.0 / 3.0 );
    weights->data[0].imag_FIXME = 0.0;
    weights->data[1].imag_FIXME = 0.0;
    weights->data[2].imag_FIXME = 0.0;

    LAL_CALL ( LALLinearlyCombineSFTVectors (&status, &SFTvect,
					     SFTVectList, weights, comboname),
	       &status );

    LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
    LAL_CALL ( LALZDestroyVector ( &status, &weights ) , &status );
    LAL_CALL ( LALDestroySFTVector ( &status, &(SFTvect) ), &status );
    LALFree ( desc );

  } /* if uvar_makeA */

  if ( uvar_makeE ) {
    CHAR *desc;
    CHAR comboname[MAX_FILENAME_LEN];
    COMPLEX16Vector *weights = NULL;

    snprintf ( comboname, MAX_FILENAME_LEN, "%s", multiTs->data[0]->name );

    comboname[0] = 'Z';
    comboname[1] = '8';
    comboname[3] = 'E';

    if ( (desc = assembleDescription ( comboname, uvar_miscField )) == NULL )
      return -1;
    if ( multiTs->length != 3 )
      {
	LogPrintf (LOG_CRITICAL,
		   "Need 3 input time series to make E; got %d\n",
		   multiTs->length);
	return LISAMAKESFTS_EINPUT;
      }
    SFTvect = NULL;

    LAL_CALL ( LALZCreateVector ( &status, &weights, 3 ) , &status );
    weights->data[0].real_FIXME = 0;
    weights->data[1].real_FIXME = (-1.0 * LAL_SQRT1_3);
    weights->data[2].real_FIXME =         LAL_SQRT1_3;
    weights->data[0].imag_FIXME = 0.0;
    weights->data[1].imag_FIXME = 0.0;
    weights->data[2].imag_FIXME = 0.0;

    LAL_CALL ( LALLinearlyCombineSFTVectors (&status, &SFTvect,
					     SFTVectList, weights, comboname),
	       &status );

    LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
    LAL_CALL ( LALZDestroyVector ( &status, &weights ) , &status );
    LAL_CALL ( LALDestroySFTVector ( &status, &(SFTvect) ), &status );
    LALFree ( desc );

  } /* if uvar_makeE */

  if ( uvar_makeT ) {
    CHAR *desc;
    CHAR comboname[MAX_FILENAME_LEN];
    COMPLEX16Vector *weights = NULL;

    snprintf ( comboname, MAX_FILENAME_LEN, "%s", multiTs->data[0]->name );

    comboname[0] = 'Z';
    comboname[1] = '9';
    comboname[3] = 'T';

    if ( (desc = assembleDescription ( comboname, uvar_miscField )) == NULL )
      return -1;
    if ( multiTs->length != 3 )
      {
	LogPrintf (LOG_CRITICAL,
		   "Need 3 input time series to make T; got %d\n",
		   multiTs->length);
	return LISAMAKESFTS_EINPUT;
      }
    SFTvect = NULL;

    LAL_CALL ( LALZCreateVector ( &status, &weights, 3 ) , &status );
    weights->data[0].real_FIXME = ( LAL_SQRT2 / 3.0 );
    weights->data[1].real_FIXME = ( LAL_SQRT2 / 3.0 );
    weights->data[2].real_FIXME = ( LAL_SQRT2 / 3.0 );
    weights->data[0].imag_FIXME = 0.0;
    weights->data[1].imag_FIXME = 0.0;
    weights->data[2].imag_FIXME = 0.0;

    LAL_CALL ( LALLinearlyCombineSFTVectors (&status, &SFTvect,
					     SFTVectList, weights, comboname),
	       &status );
    for ( sidx=0; sidx < SFTvect->length; sidx++ )
      {
	SFTvect->data[sidx].name[0] = 'Z';
	SFTvect->data[sidx].name[1] = '9';
      }
    LAL_CALL ( LALWriteSFTVector2Dir (&status, SFTvect, uvar_outputDir, add_comment, desc ), &status );
    LAL_CALL ( LALZDestroyVector ( &status, &weights ) , &status );
    LAL_CALL ( LALDestroySFTVector ( &status, &(SFTvect) ), &status );
    LALFree ( desc );

  } /* if uvar_makeT */

  /* free memory */
  LALFree ( add_comment );
  for ( ifo = 0; ifo < multiTs->length ; ifo ++ )
    {
      LAL_CALL ( LALDestroySFTVector ( &status, &(SFTVectList[ifo]) ), &status );
      XLALDestroyREAL4TimeSeries ( multiTs->data[ifo] );
    }
  LALFree ( multiTs->data );
  LALFree ( multiTs );
  LALFree ( SFTVectList );
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set defaults */
#define DEFAULT_OUTDIR  "./"
  uvar_outputDir = LALCalloc ( 1, strlen (DEFAULT_OUTDIR) + 1);
  strcpy ( uvar_outputDir, DEFAULT_OUTDIR );

  uvar_extraComment = NULL;
  uvar_miscField = NULL;
  uvar_Tsft = 604800;
  uvar_makeX = 1;
  uvar_makeY = 0;
  uvar_makeZ = 0;
  uvar_makeYminusZ = 1;
  uvar_makeZminusX = 0;
  uvar_makeXminusY = 0;
  uvar_makeA = 0;
  uvar_makeE = 0;
  uvar_makeT = 0;

  uvar_lisasim = FALSE;
  uvar_raCorrections = FALSE;

  /* now register all our user-variable */
  LALregSTRINGUserVar(status, inputXML,		'i', UVAR_REQUIRED, "XML file describing the LISA timeseries data");
  LALregREALUserVar(status,   Tsft,		'T', UVAR_OPTIONAL, "Length of SFTs to produce (in seconds)");
  LALregSTRINGUserVar(status, outputDir,	'o', UVAR_OPTIONAL, "Output directory for SFTs");
  LALregSTRINGUserVar(status, extraComment,	'C', UVAR_OPTIONAL, "Additional comment to be added to output-SFTs");
  LALregSTRINGUserVar(status, miscField,	'm', UVAR_OPTIONAL, "User-specifiable portion of the SFT-filename ('misc' field)");

  LALregBOOLUserVar(status,   lisasim,		's', UVAR_OPTIONAL, "TDI data are from LISA Simulator");

  LALregBOOLUserVar(status,   raCorrections,	 0,  UVAR_OPTIONAL, "Use rigid adiabatic approximation");

  LALregBOOLUserVar(status,   makeX,		'X', UVAR_OPTIONAL, "Produce X");
  LALregBOOLUserVar(status,   makeY,		'Y', UVAR_OPTIONAL, "Produce Y");
  LALregBOOLUserVar(status,   makeZ,		'Z', UVAR_OPTIONAL, "Produce Z");

  LALregBOOLUserVar(status,   makeYminusZ,	'x', UVAR_OPTIONAL, "Produce Y-Z (combination independent of X)");
  LALregBOOLUserVar(status,   makeZminusX,	'y', UVAR_OPTIONAL, "Produce Z-X (combination independent of Y)");
  LALregBOOLUserVar(status,   makeXminusY,	'z', UVAR_OPTIONAL, "Produce X-Y (combination independent of Z)");

  LALregBOOLUserVar(status,   makeA,		'a', UVAR_OPTIONAL, "Produce (pseudo-)A");
  LALregBOOLUserVar(status,   makeE,		'e', UVAR_OPTIONAL, "Produce (pseudo-)E");
  LALregBOOLUserVar(status,   makeT,		't', UVAR_OPTIONAL, "Produce (pseudo-)T");

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

  INITSTATUS(status);
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
      snprintf ( name, LALNameLength, "Z%d:%s_%s", i+1, thisTs->Name, lisaTs->FileName );
      name[LALNameLength-1] = 0; /* close string if it was truncated */

      /* Workaround for LISAsim metadata error: read start time from t column of data */
      XLALGPSSetREAL8( &epoch, lisaTs->Data[0]->data[0] );

      epoch.gpsSeconds += LISA_TIME_ORIGIN;	/* offset for convenience of GPS-time ranges */

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
  const CHAR illegals[] = "-. {}()/";

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

  /* Now go through and replace all illegal characters ['-', '.', ' ', '{', '}', '(', ')', '/'] by '_' */
  ptr = desc;
  while ( *ptr != 0 )
    {
      if ( strchr ( illegals, *ptr ) != NULL )
	*ptr = '_';
      ptr ++;
    }

  return desc;

} /* assembleDescription() */

/**
 * Compute inverse of "transfer function" (i.e. direction-independent part of detector-response),
 * which we use to normalize SFTs with.
 *
 * Use long-wavelength-limit (useRAA=false) or the rigid-adiabatic approximation (useRAA=true),
 * for either synthLISA (isLISAsim=false) or isLISAsim=true
 */
void
compute_R_f ( COMPLEX8 *R_f, REAL8 Freq, BOOLEAN useRAA, BOOLEAN isLISAsim )
{

  REAL8 twopifL = LAL_TWOPI * LISA_ARM_LENGTH_SECONDS * Freq;
  REAL8 fourpifL = 2.0 * twopifL;

  if ( useRAA )
    {
      if ( isLISAsim ) 	/* RAA && LISAsim */
	{
	  R_f->realf_FIXME = - 0.5 * sin(fourpifL) / sin(twopifL);
	  R_f->imagf_FIXME =   0.5 * cos(fourpifL) / sin(twopifL);
	}
      else 		/* RAA && synthLISA */
	{
	  R_f->realf_FIXME = (0.5 / fourpifL) * cos(fourpifL) / sin(twopifL);
	  R_f->imagf_FIXME = (0.5 / fourpifL) * sin(fourpifL) / sin(twopifL);
	}
    } /* if useRAA */
  else
    {
      if ( isLISAsim )	/* LWL && LISAsim */
	{
	  R_f->realf_FIXME = 0;
	  R_f->imagf_FIXME = 1.0 / fourpifL;
	}
      else 		/* LWL && synthLISA */
	{
	  R_f->realf_FIXME = 1.0 / ( fourpifL * fourpifL );
	  R_f->imagf_FIXME = 0;
	}
    } /* if LWL */

  return;

} /* compute_R_f() */
