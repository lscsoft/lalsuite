/*
*  Copyright (C) 2009 Reinhard Prix
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
 * \file
 * \ingroup pulsarApps
 * \author R. Prix
 * \brief
 * compare two binary strain files, written from hardware-injection
 * streams, return maximal relative error.
 *
 * Note: the strain file need to start with magic REAL4 1234.5,
 * followed by an INT4 argument giving the number of elements.
 */

/* ---------- includes ---------- */
#include <stdio.h>

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/SFTutils.h>

#include <lalapps.h>

/* Error codes and messages */
#define MAKEFAKEDATAC_ENORM 	0
#define MAKEFAKEDATAC_ESUB  	1
#define MAKEFAKEDATAC_EARG  	2
#define MAKEFAKEDATAC_EBAD  	3
#define MAKEFAKEDATAC_EFILE 	4
#define MAKEFAKEDATAC_ENOARG 	5
#define MAKEFAKEDATAC_EINCOMPAT 6
#define MAKEFAKEDATAC_ENULL	7

#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "File IO error"
#define MAKEFAKEDATAC_MSGENOARG "Missing argument"
#define MAKEFAKEDATAC_MSGEINCOMPAT "Incompatible SFTs"
#define MAKEFAKEDATAC_MSGENULL	"Unexpected null pointer"

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

#define MYMAX(a,b) ( (a>b) ? a:b )
#define SQ(x) ( (x) * (x) )


/* User variables */
typedef struct {
  CHAR *infile1;
  CHAR *infile2;
  INT4 debug;
  BOOLEAN verbose;
  BOOLEAN help;
  BOOLEAN version;
  REAL8 relErrorMax;
} UserVar;


/* ----- local prototypes ----- */
REAL4 getMaxErrSFT (const SFTtype *sft1, const SFTtype *sft2);
REAL4 getMaxErrSFTVector (const SFTVector *sftvect1, const SFTVector *sftvect2);
void scalarProductSFT (LALStatus *stat, REAL4 *scalar, const SFTtype *sft1, const SFTtype *sft2);
void scalarProductSFTVector (LALStatus *stat, REAL4 *scalar, const SFTVector *sftvect1, const SFTVector *sftvect2);
void subtractSFTVectors (LALStatus *stat, SFTVector **ret, const SFTVector *sftvect1, const SFTVector *sftvect2);
void initUserVars (LALStatus *status, UserVar *uvar);
REAL4 XLALcompareREAL4Vectors ( REAL4Vector *ts1, REAL4Vector *ts2 );
REAL4Vector *XLALREAL4VectorFromFile ( const CHAR *fname );

extern int vrbflg;

/* Empty initializers */
static const LALStatus empty_LALStatus;
static const UserVar empty_UserVar;

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  LALStatus status = empty_LALStatus;	/* initialize status */
  UserVar uvar = empty_UserVar;
  REAL4Vector *ts1 = NULL, *ts2 = NULL;
  REAL4 maxd;


  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL ( initUserVars (&status, &uvar), &status );

  /* read cmdline & cfgfile  */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);

  if (uvar.help) 	/* help requested: we're done */
    exit (0);

  if (uvar.version)
    {
      CHAR *VCSInfoString;
      if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
        XLALPrintError("XLALGetVersionString(0) failed.\n");
        exit(1);
      }
      printf ("%s\n", VCSInfoString );
      XLALFree ( VCSInfoString );
      exit(0);
    }

  LogSetLevel ( lalDebugLevel );

  /* now read in the two timeseries */
  if ( (ts1 = XLALREAL4VectorFromFile ( uvar.infile1 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: failed to load timeseries file '%s'.\n", __func__, uvar.infile1 );
    return MAKEFAKEDATAC_EFILE;
  }

  if ( (ts2 = XLALREAL4VectorFromFile ( uvar.infile2 )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "%s: failed to load timeseries file '%s'.\n", __func__, uvar.infile2 );
    return MAKEFAKEDATAC_EFILE;
  }

  maxd = XLALcompareREAL4Vectors ( ts1, ts2 );
  if ( XLAL_IS_REAL8_FAIL_NAN ( maxd ) ) {
    LogPrintf (LOG_CRITICAL, "%s: XLALcompareTS() failed. xlalErrno = %d.\n", __func__, xlalErrno );
    return XLAL_EFUNC;
  }

  if ( maxd > uvar.relErrorMax )
    {
      LogPrintf (LOG_CRITICAL, "%s: FAILED. Maximal relative error %e exceeds tolerance %e.\n", __func__, maxd, uvar.relErrorMax );
      return XLAL_EFAILED;
    }

  LogPrintf (LOG_DEBUG, "%s: OK. Maximal relative error %e is within tolerance of %e.\n", __func__, maxd, uvar.relErrorMax );

  /* free memory */
  XLALDestroyREAL4Vector ( ts1 );
  XLALDestroyREAL4Vector ( ts2 );

  LALDestroyUserVars (&status);

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *status, UserVar *uvar)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set some defaults */
  uvar->debug = lalDebugLevel;
  uvar->verbose = FALSE;
  uvar->relErrorMax = 1e-4;

  /* now register all our user-variable */
  LALregBOOLUserStruct(status,   help,	        'h', UVAR_HELP,     "Print this help/usage message");
  LALregSTRINGUserStruct(status, infile1,	'1', UVAR_REQUIRED, "First timeseries input file");
  LALregSTRINGUserStruct(status, infile2,	'2', UVAR_REQUIRED, "Second timeseries input file");
  LALregBOOLUserStruct(status,   verbose,	'v', UVAR_OPTIONAL, "Verbose output of differences");
  LALregREALUserStruct(status,   relErrorMax,   'e', UVAR_OPTIONAL, "Maximal relative error acceptable to 'pass' comparison");
  LALregBOOLUserStruct(status,   version,	'V', UVAR_SPECIAL,  "Output version information");


  DETATCHSTATUSPTR (status );
  RETURN (status );

} /* initUserVars() */


/**
 * Compare two REAL8 vectors, returns a measure of the difference.
 */
REAL4
XLALcompareREAL4Vectors ( REAL4Vector *ts1, REAL4Vector *ts2 )
{
  UINT4 i, numSteps;
  REAL8 total_power, maxdiff, sumdiff, maxpower;
  REAL8 reldiff, reldiff_max, reldiff_avg;

  if ( !ts1 || !ts2 || !ts1->data || !ts2->data ) {
    XLALPrintError ("%s: illegal NULL input.\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EINVAL );
  }

  numSteps = ts1->length;
  if ( ts2->length != numSteps ) {
    XLALPrintError ("%s: number of timesteps of ts1 (%d) differs from ts2 (%d).\n", __func__, numSteps, ts2->length );
    XLAL_ERROR_REAL8 ( XLAL_EBADLEN );
  }

  sumdiff = 0;
  maxdiff = 0;
  total_power = 0;
  maxpower = 0;
  for ( i=0; i < numSteps; i ++ )
    {
      REAL8 power, diff;

      diff = fabs ( ts1->data[i] - ts2->data[i] );

      if ( diff > maxdiff )
	maxdiff = diff;

      sumdiff += diff;

      power = fabs ( ts1->data[i] );
      if ( power > maxpower )
	maxpower = power;

      total_power += power;

    } /* for i < numSteps */

  reldiff_max = maxdiff / maxpower;
  reldiff_avg = sumdiff / total_power;

  LogPrintf (LOG_DEBUG, "%s: maximal difference = %g, maximal amplitude = %g ==> relative error %g\n", __func__, maxdiff, maxpower, reldiff_max );
  LogPrintf (LOG_DEBUG, "%s: total difference = %g, total summed amplitude = %g ==> relative avg error %g\n", __func__, sumdiff, total_power, reldiff_avg );

  reldiff = MYMAX ( reldiff_max, reldiff_avg );

  return (REAL4)reldiff;

  } /* XLALcompareREAL4Vectors() */


/**
 * Load timeseries from binary input file into REAL4Vector
 */
REAL4Vector *
XLALREAL4VectorFromFile ( const CHAR *fname )
{
  REAL4Vector *vect;
  FILE *fp;
  REAL4 magic;
  UINT4 len;

  if ( !fname ) {
    XLALPrintError ("%s: filename is NULL.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  if ( (fp = fopen ( fname, "rb")) == NULL ) {
    XLALPrintError ("%s: failed to open file '%s' for reading.\n", __func__, fname );
    XLAL_ERROR_NULL ( XLAL_ESYS );
  }

  if ( (fread( &magic, sizeof(magic), 1, fp ) != 1) || (magic != 1234.5) ) {
    XLALPrintError ("%s: file '%s' has wrong magic byte (%f) != 1234.5.\n", __func__, fname, magic );
    fclose(fp);
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  if ( fread ( &len, sizeof(len), 1, fp ) != 1 ) {
    XLALPrintError ("%s: failed to read UINT4 length from file '%s'\n", __func__, fname );
    fclose(fp);
    XLAL_ERROR_NULL ( XLAL_ESYS );
  }

  if ( ( vect = XLALCreateREAL4Vector ( len )) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, len );
    fclose(fp);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  if ( len != fread ( vect->data, sizeof(*vect->data), len, fp ) ) {
    XLALPrintError ("%s: failed to read %d REAL4s from input file '%s'.\n", __func__, len, fname );
    fclose (fp);
    XLALDestroyREAL4Vector ( vect );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  fclose (fp );

  return vect;

} /* XLALREAL4VectorFromFile() */
