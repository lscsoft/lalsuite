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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup lalpulsar_bin_SFTTools
 * \author R. Prix
 */

/* ---------- includes ---------- */
#include "config.h"

#include <lal/Date.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LALPulsarVCSInfo.h>

/* User variables */
typedef struct
{
  CHAR *sftBname1;
  CHAR *sftBname2;
  INT4 debug;
  BOOLEAN verbose;
  BOOLEAN quiet;
  REAL8 relErrorMax;
} UserInput_t;

/* local prototypes */
int initUserVars ( UserInput_t *uvar );
REAL4 getMaxErrSFT (const SFTtype *sft1, const SFTtype *sft2);
REAL4 getMaxErrSFTVector (const SFTVector *sftvect1, const SFTVector *sftvect2);
REAL4 scalarProductSFT ( const SFTtype *sft1, const SFTtype *sft2 );
REAL4 scalarProductSFTVector ( const SFTVector *sftvect1, const SFTVector *sftvect2 );
SFTVector *subtractSFTVectors ( const SFTVector *sftvect1, const SFTVector *sftvect2 );

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  SFTVector *diffs = NULL;
  REAL8 maxd = 0;

  /* register all user-variables */
  UserInput_t XLAL_INIT_DECL(uvar);
  XLAL_CHECK_MAIN ( initUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  SFTCatalog *catalog;
  XLAL_CHECK_MAIN ( (catalog = XLALSFTdataFind ( uvar.sftBname1, NULL )) != NULL, XLAL_EFUNC );
  SFTVector *SFTs1;
  XLAL_CHECK_MAIN ( (SFTs1 = XLALLoadSFTs( catalog, -1, -1 )) != NULL, XLAL_EFUNC );
  XLALDestroySFTCatalog ( catalog );

  XLAL_CHECK_MAIN ( (catalog = XLALSFTdataFind ( uvar.sftBname2, NULL )) != NULL, XLAL_EFUNC );
  SFTVector *SFTs2;
  XLAL_CHECK_MAIN ( (SFTs2 = XLALLoadSFTs( catalog, -1, -1 )) != NULL, XLAL_EFUNC );
  XLALDestroySFTCatalog ( catalog );

  /* ---------- do some sanity checks of consistency of SFTs ----------*/
  XLAL_CHECK_MAIN ( SFTs1->length == SFTs2->length, XLAL_EINVAL, "Number of SFTs differ for SFTbname1 and SFTbname2!\n");

  for ( UINT4 i=0; i < SFTs1->length; i++ )
    {
      SFTtype *sft1 = &(SFTs1->data[i]);
      SFTtype *sft2 = &(SFTs2->data[i]);

      XLAL_CHECK_MAIN( strcmp( sft1->name, sft2->name ) == 0, XLAL_EINVAL, "\nERROR SFT %d: detector-prefix differ! '%s' != '%s'\n", i, sft1->name, sft2->name );

      XLAL_CHECK_MAIN ( sft1->data->length == sft2->data->length, XLAL_EINVAL, "\nERROR SFT %d: lengths differ! %d != %d\n", i, sft1->data->length, sft2->data->length );

      REAL8 Tdiff = XLALGPSDiff(&(sft1->epoch), &(sft2->epoch));
      CHAR buf1[32], buf2[32];;
      XLAL_CHECK_MAIN ( Tdiff == 0.0, XLAL_EINVAL, "SFT %d: epochs differ: %s vs %s\n", i, XLALGPSToStr(buf1,&sft1->epoch), XLALGPSToStr(buf2,&sft2->epoch) );
      XLAL_CHECK_MAIN ( sft1->f0 == sft2->f0, XLAL_EINVAL, "ERROR SFT %d: fmin differ: %fHz vs %fHz\n", i, sft1->f0, sft2->f0 );
      XLAL_CHECK_MAIN ( sft1->deltaF == sft2->deltaF, XLAL_EINVAL, "ERROR SFT %d: deltaF differs: %fHz vs %fHz\n", i, sft1->deltaF, sft2->deltaF );
    } /* for i < numSFTs */

  /*---------- now do some actual comparisons ----------*/
  XLAL_CHECK_MAIN ( (diffs = subtractSFTVectors ( SFTs1, SFTs2)) != NULL, XLAL_EFUNC );

  if ( uvar.verbose)
    {
      for ( UINT4 i=0; i < SFTs1->length; i++)
	{
	  SFTtype *sft1 = &(SFTs1->data[i]);
	  SFTtype *sft2 = &(SFTs2->data[i]);

	  XLALPrintInfo ("i=%02d: ", i);
	  REAL4 norm1 = scalarProductSFT ( sft1, sft1 );
	  norm1 = sqrt(norm1);
	  REAL4 norm2 = scalarProductSFT ( sft2, sft2 );
	  norm2 = sqrt(norm2);
	  REAL4 scalar = scalarProductSFT ( sft1, sft2 );

	  REAL4 normdiff = scalarProductSFT ( &(diffs->data[i]), &(diffs->data[i]) );

	  REAL4 d1 = (norm1 - norm2)/norm1;
	  REAL4 d2 = 1.0 - scalar / (norm1*norm2);
	  REAL4 d3 = normdiff / (norm1*norm1 + norm2*norm2 );
	  REAL4 d4 = getMaxErrSFT (sft1, sft2);
	  XLALPrintInfo ("SFT #%d: (|x|-|y|)/|x|=%10.3e, 1-x.y/(|x||y|)=%10.3e, |x-y|^2/(|x|^2+|y|^2))=%10.3e, maxErrSFT=%10.3e\n", i, d1, d2, d3, d4);
	} /* for i < SFTs->length */
    } /* if verbose */

  /* ---------- COMBINED measures ---------- */
  {
    REAL4 ret;
    ret = scalarProductSFTVector ( SFTs1, SFTs1 );
    REAL8 norm1 = sqrt( (REAL8)ret );

    ret = scalarProductSFTVector ( SFTs2, SFTs2 );
    REAL8 norm2 = sqrt( (REAL8)ret );

    ret = scalarProductSFTVector ( SFTs1, SFTs2 );
    REAL8 scalar = (REAL8) ret;

    ret = scalarProductSFTVector ( diffs, diffs );
    REAL8 normdiff = (REAL8) ret;

    REAL8 d1 = (norm1 - norm2)/norm1;
    maxd = fmax ( maxd, d1 );
    REAL8 d2 = 1.0 - scalar / (norm1*norm2);
    maxd = fmax ( maxd, d2 );
    REAL8 d3 = normdiff / ( norm1*norm1 + norm2*norm2);
    maxd = fmax ( maxd, d3 );
    REAL8 d4 = getMaxErrSFTVector (SFTs1, SFTs2);
    maxd = fmax ( maxd, d4 );

    if ( uvar.verbose ) {
      printf ("TOTAL: (|x|-|y|)/|x|=%10.3e, 1-x.y/(|x||y|)=%10.3e, |x-y|^2/(|x|^2+|y|^2)=%10.3e, maxErrSFT=%10.3e\n", d1, d2, d3, d4);
      printf ("COMPARE: maxd=%10.3e, relErrMax=%10.3e\n", maxd, uvar.relErrorMax);
    }

  } /* combined total measures */


  /* free memory */
  XLALDestroySFTVector ( SFTs1 );
  XLALDestroySFTVector ( SFTs2 );
  XLALDestroySFTVector ( diffs );
  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  if ( maxd <= uvar.relErrorMax ) {
    return 0;
  }
  else {
    if ( !uvar.quiet ) {
      XLALPrintError("Tolerance exceeded! maxd=%10.3e, relErrMax=%10.3e\n", maxd, uvar.relErrorMax);
    }
    return 1;
  }

} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
int
initUserVars ( UserInput_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  /* set some defaults */
  uvar->debug = lalDebugLevel;
  uvar->verbose = 0;
  uvar->quiet = 0;
  uvar->relErrorMax = 1e-4;

  /* now register all our user-variable */

  XLAL_CHECK ( XLALRegisterUvarMember( sftBname1,       STRING,  '1', REQUIRED, "Path and basefilename for SFTs1. Possibilities are:\n"
                                       " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( sftBname2,       STRING,  '2', REQUIRED, "Path and basefilename for SFTs2. Possibilities are:\n"
                                       " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( verbose,         BOOLEAN, 'V', OPTIONAL, "Verbose output of differences") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( quiet,           BOOLEAN, 'q', OPTIONAL, "No output of differences even on failure") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( relErrorMax,     REAL8,   'e', OPTIONAL, "Maximal relative error acceptable to 'pass' comparison") == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} /* initUserVars() */

/* for two SFTs: get maximal value of |X_k - Y_k|^2 / max(|X_k|^2,|Y_k|^2) */
REAL4
getMaxErrSFT (const SFTtype *sft1, const SFTtype *sft2)
{
  UINT4 i;
  REAL8 maxDiff, maxAmpl;
  REAL4 maxErr;

  maxDiff = 0;
  maxAmpl = 0;
  for (i=0; i < sft1->data->length; i++)
    {
      REAL8 diff, A1, A2, Ampl;
      REAL8 re1, re2, im1, im2;
      re1 = crealf(sft1->data->data[i]);
      im1 = cimagf(sft1->data->data[i]);
      re2 = crealf(sft2->data->data[i]);
      im2 = cimagf(sft2->data->data[i]);

      diff = (re1 - re2)*(re1 - re2) + (im1 - im2)*(im1 - im2);
      A1 = re1*re1 + im1*im1;
      A2 = re2*re2 + im2*im2;
      Ampl = fmax(A1, A2);

      maxDiff = fmax(maxDiff, diff);
      maxAmpl = fmax(maxAmpl, Ampl);

    } /* for i */

  maxErr = maxDiff / maxAmpl;

  return(maxErr);

} /* getMaxErrSFT() */

REAL4
getMaxErrSFTVector (const SFTVector *sftvect1, const SFTVector *sftvect2)
{
  UINT4 i;
  REAL4 maxErr, thisErr;

  maxErr = 0;

  for (i=0; i<sftvect1->length; i++)
    {
      thisErr = getMaxErrSFT ( &(sftvect1->data[i]), &(sftvect2->data[i]));
      maxErr = fmax (maxErr, thisErr);
    }

  return maxErr;

} /* getMaxErrSFTVector() */


/*--------------------------------------------------
 * implements a straightforward L2 scalar product of
 * two time-series x_i and y_i : x*y = sum_i x_i y_i
 * in Fourier-space, which is 2/N * Re( sum_i X_i Y*_i)
 *--------------------------------------------------*/
REAL4
scalarProductSFT ( const SFTtype *sft1, const SFTtype *sft2 )
{
  XLAL_CHECK_REAL4 ( (sft1 != NULL) && (sft2 != NULL), XLAL_EINVAL );
  XLAL_CHECK_REAL4 ( sft1->data->length == sft2->data->length, XLAL_EINVAL );

  /* we do the calculation in REAL8 to avoid accumulating roundoff-errors */
  REAL8 prod = 0;
  for ( UINT4 i=0; i < sft1->data->length; i++ )
    {
      REAL8 xre = (REAL8)crealf(sft1->data->data[i]);
      REAL8 xim = (REAL8)cimagf(sft1->data->data[i]);
      REAL8 yre = (REAL8)crealf(sft2->data->data[i]);
      REAL8 yim = (REAL8)cimagf(sft2->data->data[i]);

      prod +=  xre * yre + xim * yim;

    } /* for i < SFT-length */

  prod *= 2.0 / ((REAL8)sft1->data->length);

  return (REAL4) prod;

} /* scalarProductSFT() */

/*--------------------------------------------------
 * extend the previous definition to an SFT-vector
 * this is simply the sum of individual SFT-products
 *--------------------------------------------------*/
REAL4
scalarProductSFTVector ( const SFTVector *sftvect1, const SFTVector *sftvect2 )
{
  XLAL_CHECK_REAL4 ( (sftvect1 != NULL) && (sftvect2 != NULL), XLAL_EINVAL );
  XLAL_CHECK_REAL4 ( sftvect1->length == sftvect2->length, XLAL_EINVAL );

  REAL8 prod = 0;
  for ( UINT4 i=0; i < sftvect1->length; i++ )
    {
      REAL4 xy = scalarProductSFT ( &(sftvect1->data[i]), &(sftvect2->data[i]) );
      prod += (REAL8) xy;
    }

  return (REAL4)prod;

} /* scalarProductSFTVector() */


/*--------------------------------------------------
 * calculate the difference of two SFT-vectors
 *--------------------------------------------------*/
SFTVector *
subtractSFTVectors ( const SFTVector *sftvect1, const SFTVector *sftvect2 )
{
  XLAL_CHECK_NULL ( (sftvect1 != NULL) && (sftvect2 != NULL), XLAL_EINVAL );
  XLAL_CHECK_NULL ( sftvect1->length == sftvect2->length, XLAL_EINVAL );

  UINT4 numSFTs = sftvect1->length;
  UINT4 numBins = sftvect1->data[0].data->length;

  SFTVector *vect;
  XLAL_CHECK_NULL ( (vect = XLALCreateSFTVector ( numSFTs, numBins )) != NULL, XLAL_EFUNC );

  for ( UINT4 alpha = 0; alpha < numSFTs; alpha ++ )
    {
      for ( UINT4 j=0; j < numBins; j++ )
	{
	  vect->data[alpha].data->data[j] = sftvect1->data[alpha].data->data[j] - sftvect2->data[alpha].data->data[j];
	}
    } /* for alpha < numSFTs */

  return vect;

} /* subtractSFTVectors() */
