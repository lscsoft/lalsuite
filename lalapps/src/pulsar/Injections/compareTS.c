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
 * \ingroup lalapps_pulsar_Tools
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

/* User variables */
typedef struct {
  CHAR *infile1;
  CHAR *infile2;
  INT4 debug;
  BOOLEAN verbose;
  BOOLEAN version;
  REAL8 relErrorMax;
} UserVar;


/* ----- local prototypes ----- */
REAL4 getMaxErrSFT (const SFTtype *sft1, const SFTtype *sft2);
REAL4 getMaxErrSFTVector (const SFTVector *sftvect1, const SFTVector *sftvect2);
int initUserVars ( UserVar *uvar );
REAL4 XLALcompareREAL4Vectors ( REAL4Vector *ts1, REAL4Vector *ts2 );
REAL4Vector *XLALREAL4VectorFromFile ( const CHAR *fname );

extern int vrbflg;

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  /* register all user-variables */
  UserVar XLAL_INIT_DECL(uvar);
  XLAL_CHECK_MAIN ( initUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  if (uvar.version)
    {
      CHAR *VCSInfoString;
      XLAL_CHECK_MAIN ( (VCSInfoString = XLALGetVersionString(0)) != NULL, XLAL_EFUNC );
      printf ("%s\n", VCSInfoString );
      XLALFree ( VCSInfoString );
      exit(0);
    }

  /* now read in the two timeseries */
  REAL4Vector *ts1, *ts2;
  XLAL_CHECK_MAIN ( (ts1 = XLALREAL4VectorFromFile ( uvar.infile1 )) != NULL, XLAL_EFUNC, "Failed to load timeseries file '%s'.\n", uvar.infile1 );
  XLAL_CHECK_MAIN ( (ts2 = XLALREAL4VectorFromFile ( uvar.infile2 )) != NULL, XLAL_EFUNC, "Failed to load timeseries file '%s'.\n", uvar.infile2 );

  REAL4 maxd = XLALcompareREAL4Vectors ( ts1, ts2 );
  XLAL_CHECK_MAIN ( ! XLAL_IS_REAL8_FAIL_NAN ( maxd ), XLAL_EFUNC );

  XLAL_CHECK_MAIN ( maxd <= uvar.relErrorMax, XLAL_ETOL, "FAILED. Maximal relative error %e exceeds tolerance %e.\n", maxd, uvar.relErrorMax );

  XLALPrintInfo ("OK. Maximal relative error %e is within tolerance of %e.\n", maxd, uvar.relErrorMax );

  /* free memory */
  XLALDestroyREAL4Vector ( ts1 );
  XLALDestroyREAL4Vector ( ts2 );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
int
initUserVars ( UserVar *uvar )
{
  /* set some defaults */
  uvar->debug = lalDebugLevel;
  uvar->verbose = 0;
  uvar->relErrorMax = 1e-4;

  /* now register all our user-variable */
  XLALRegisterUvarMember(	infile1,	STRING, '1', REQUIRED, "First timeseries input file");
  XLALRegisterUvarMember( 	infile2,	STRING, '2', REQUIRED, "Second timeseries input file");
  XLALRegisterUvarMember( 	verbose,	BOOLEAN, 'v', OPTIONAL, "Verbose output of differences");
  XLALRegisterUvarMember( 	relErrorMax,   	REAL8, 'e', OPTIONAL, "Maximal relative error acceptable to 'pass' comparison");
  XLALRegisterUvarMember( 	version,	BOOLEAN, 'V', SPECIAL,  "Output version information");

  return XLAL_SUCCESS;
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

  reldiff = fmax ( reldiff_max, reldiff_avg );

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
