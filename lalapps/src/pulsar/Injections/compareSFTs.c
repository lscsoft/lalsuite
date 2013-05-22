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
 * \file
 * \ingroup pulsarApps
 * \author R. Prix
 */

/* ---------- includes ---------- */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>

/* Error codes and messages */
/**\name Error Codes */ /*@{*/
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

/*@}*/

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

#define mymax(a,b) ( (a>b) ? a:b )
/* local prototypes */
/* Prototypes for the functions defined in this file */
void initUserVars (LALStatus *stat);
REAL4 getMaxErrSFT (const SFTtype *sft1, const SFTtype *sft2);
REAL4 getMaxErrSFTVector (const SFTVector *sftvect1, const SFTVector *sftvect2);
void scalarProductSFT (LALStatus *stat, REAL4 *scalar, const SFTtype *sft1, const SFTtype *sft2);
void scalarProductSFTVector (LALStatus *stat, REAL4 *scalar, const SFTVector *sftvect1, const SFTVector *sftvect2);
void subtractSFTVectors (LALStatus *stat, SFTVector **ret, const SFTVector *sftvect1, const SFTVector *sftvect2);

extern int vrbflg;

/*----------------------------------------------------------------------*/
static const LALStatus empty_status;
static const SFTConstraints empty_constraints;
/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_sftBname1;
CHAR *uvar_sftBname2;
INT4 uvar_debug;
BOOLEAN uvar_verbose;
BOOLEAN uvar_help;
REAL8 uvar_relErrorMax;

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  LALStatus status = empty_status;	/* initialize status */
  SFTConstraints constraints = empty_constraints;
  CHAR detector[2] = "??";	/* allow reading v1-SFTs without detector-info */
  SFTVector *SFTs1 = NULL, *SFTs2 = NULL;
  SFTVector *diffs = NULL;
  SFTCatalog *catalog = NULL;
  UINT4 i;
  REAL8 maxd = 0;

  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (initUserVars (&status), &status);

  /* read cmdline & cfgfile  */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* now read in the two complete sft-vectors */
  constraints.detector = detector;

  LAL_CALL (LALSFTdataFind (&status, &catalog, uvar_sftBname1, &constraints ), &status );
  LAL_CALL (LALLoadSFTs(&status, &SFTs1, catalog, -1, -1 ), &status );
  LAL_CALL (LALDestroySFTCatalog(&status, &catalog), &status );
  catalog = NULL;

  LAL_CALL (LALSFTdataFind (&status, &catalog, uvar_sftBname2, &constraints ), &status );
  LAL_CALL (LALLoadSFTs(&status, &SFTs2, catalog, -1, -1 ), &status );
  LAL_CALL (LALDestroySFTCatalog(&status, &catalog), &status );

  /* ---------- do some sanity checks of consistency of SFTs ----------*/
  if (SFTs1->length != SFTs2->length) {
    XLALPrintError ("Warning: number of SFTs differ for SFTbname1 and SFTbname2!\n");
    exit(1);
  }
  for (i=0; i < SFTs1->length; i++)
    {
      REAL8 Tdiff;
      SFTtype *sft1, *sft2;
      sft1 = &(SFTs1->data[i]);
      sft2 = &(SFTs2->data[i]);

      if( strcmp( sft1->name, sft2->name ) )
	{
	  if ( lalDebugLevel ) XLALPrintError("WARNING SFT %d: detector-prefix differ! '%s' != '%s'\n", i, sft1->name, sft2->name );
	  /* exit (1); */  /* can't be too strict here, as we also allow v1-SFTs, which don't have detector-name */
	}

      if (sft1->data->length != sft2->data->length)
	{
	  XLALPrintError ("\nERROR SFT %d: lengths differ! %d != %d\n", i, sft1->data->length, sft2->data->length);
	  exit(1);
	}
      Tdiff = XLALGPSDiff(&(sft1->epoch), &(sft2->epoch));
      if ( Tdiff != 0.0 )
	XLALPrintError ("WARNING SFT %d: epochs differ: (%d s, %d ns)  vs (%d s, %d ns)\n", i,
		       sft1->epoch.gpsSeconds, sft1->epoch.gpsNanoSeconds, sft2->epoch.gpsSeconds, sft2->epoch.gpsNanoSeconds);

      if ( sft1->f0 != sft2->f0)
	{
	  XLALPrintError ("ERROR SFT %d: fmin differ: %fHz vs %fHz\n", i, sft1->f0, sft2->f0);
	  exit(1);
	}

      if ( sft1->deltaF != sft2->deltaF )
	{
	  XLALPrintError ("ERROR SFT %d: deltaF differs: %fHz vs %fHz\n", i, sft1->deltaF, sft2->deltaF);
	  exit(1);
	}
    } /* for i < numSFTs */

  /*---------- now do some actual comparisons ----------*/
  LAL_CALL (subtractSFTVectors (&status, &diffs, SFTs1, SFTs2), &status);

  if ( uvar_verbose)
    {
      for (i=0; i < SFTs1->length; i++)
	{
	  SFTtype *sft1 = &(SFTs1->data[i]);
	  SFTtype *sft2 = &(SFTs2->data[i]);

	  REAL4 d1, d2, d3, d4;
	  REAL4 scalar, norm1, norm2, normdiff;
	  printf ("i=%02d: ", i);
	  LAL_CALL( scalarProductSFT(&status, &norm1, sft1, sft1 ), &status);
	  norm1 = sqrt(norm1);
	  LAL_CALL( scalarProductSFT(&status, &norm2, sft2, sft2 ), &status);
	  norm2 = sqrt(norm2);
	  LAL_CALL( scalarProductSFT(&status, &scalar, sft1, sft2 ), &status);

	  LAL_CALL( scalarProductSFT(&status, &normdiff, &(diffs->data[i]), &(diffs->data[i]) ), &status);

	  d1 = (norm1 - norm2)/norm1;
	  d2 = 1.0 - scalar / (norm1*norm2);
	  d3 = normdiff / (norm1*norm1 + norm2*norm2 );
	  d4 = getMaxErrSFT (sft1, sft2);
	  printf ("(|x|-|y|)/|x|=%10.3e, 1-x.y/(|x||y|)=%10.3e, |x-y|^2/(|x|^2+|y|^2))=%10.3e, maxErr=%10.3e\n", d1, d2, d3, d4);
	} /* for i < SFTs->length */
    } /* if verbose */

  /* ---------- COMBINED measures ---------- */
  {
    REAL4 ret;
    REAL8 norm1, norm2, normdiff, scalar;
    REAL8 d1, d2, d3,d4;

    LAL_CALL( scalarProductSFTVector (&status, &ret, SFTs1, SFTs1 ), &status);
    norm1 = sqrt( (REAL8)ret );

    LAL_CALL( scalarProductSFTVector (&status, &ret, SFTs2, SFTs2 ), &status);
    norm2 = sqrt( (REAL8)ret );

    LAL_CALL( scalarProductSFTVector (&status, &ret, SFTs1, SFTs2 ), &status);
    scalar = (REAL8) ret;

    LAL_CALL( scalarProductSFTVector(&status, &ret, diffs, diffs ), &status);
    normdiff = (REAL8) ret;

    d1 = (norm1 - norm2)/norm1;
    maxd = mymax ( maxd, d1 );
    d2 = 1.0 - scalar / (norm1*norm2);
    maxd = mymax ( maxd, d2 );
    d3 = normdiff / ( norm1*norm1 + norm2*norm2);
    maxd = mymax ( maxd, d3 );
    d4 = getMaxErrSFTVector (SFTs1, SFTs2);
    maxd = mymax ( maxd, d4 );

    if ( uvar_verbose )
      printf ("\nTOTAL:(|x|-|y|)/|x|=%10.3e, 1-x.y/(|x||y|)=%10.3e, |x-y|^2/(|x|^2+|y|^2)=%10.3e, maxErr=%10.3e\n", d1, d2, d3, d4);
    else
      printf ("%10.3e  %10.3e\n", d3, d4);

  } /* combined total measures */


  /* free memory */
  LAL_CALL (LALDestroySFTVector(&status, &SFTs1), &status);
  LAL_CALL (LALDestroySFTVector(&status, &SFTs2), &status);
  LAL_CALL (LALDestroySFTVector(&status, &diffs), &status);
  LAL_CALL (LALDestroyUserVars (&status), &status);


  LALCheckMemoryLeaks();

  if ( maxd <= uvar_relErrorMax )
    return 0;
  else
    return 1;

} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *stat)
{
  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  /* set some defaults */
  uvar_debug = lalDebugLevel;
  uvar_verbose = FALSE;
  uvar_relErrorMax = 1e-4;

  /* now register all our user-variable */

  LALregBOOLUserVar(stat,   help,	'h', UVAR_HELP,     "Print this help/usage message");
  LALregSTRINGUserVar(stat, sftBname1,	'1', UVAR_REQUIRED, "Path and basefilename for SFTs1");
  LALregSTRINGUserVar(stat, sftBname2,	'2', UVAR_REQUIRED, "Path and basefilename for SFTs2");
  LALregBOOLUserVar(stat,   verbose,	'v', UVAR_OPTIONAL, "Verbose output of differences");
  LALregREALUserVar(stat,   relErrorMax,'e', UVAR_OPTIONAL, "Maximal relative error acceptable to 'pass' comparison");


  DETATCHSTATUSPTR (stat);
  RETURN (stat);

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
      Ampl = mymax(A1, A2);

      maxDiff = mymax(maxDiff, diff);
      maxAmpl = mymax(maxAmpl, Ampl);

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
      maxErr = mymax (maxErr, thisErr);
    }

  return(maxErr);

} /* getMaxErrSFTVector() */


/*--------------------------------------------------
 * implements a straightforward L2 scalar product of
 * two time-series x_i and y_i : x*y = sum_i x_i y_i
 * in Fourier-space, which is 2/N * Re( sum_i X_i Y*_i)
 *--------------------------------------------------*/
void
scalarProductSFT (LALStatus *stat, REAL4 *scalar, const SFTtype *sft1, const SFTtype *sft2)
{
  UINT4 i;
  REAL8 prod;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  ASSERT ( scalar, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sft1, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sft2, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sft1->data->length == sft2->data->length, stat, MAKEFAKEDATAC_EINCOMPAT, MAKEFAKEDATAC_MSGEINCOMPAT);

  /* we do the calculation in REAL8 to avoid accumulating roundoff-errors */
  prod = 0;
  for (i=0; i < sft1->data->length; i++)
    {
      REAL8 xre, xim, yre, yim;
      xre = (REAL8)crealf(sft1->data->data[i]);
      xim = (REAL8)cimagf(sft1->data->data[i]);
      yre = (REAL8)crealf(sft2->data->data[i]);
      yim = (REAL8)cimagf(sft2->data->data[i]);

      prod +=  xre * yre + xim * yim;

    } /* for i < SFT-length */

  prod *= 2.0 / ((REAL8)sft1->data->length);

  /* cast back to REAL4 */
  (*scalar) = (REAL4) prod;


  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* scalarProductSFT() */

/*--------------------------------------------------
 * extend the previous definition to an SFT-vector
 * this is simply the sum of individual SFT-products
 *--------------------------------------------------*/
void
scalarProductSFTVector (LALStatus *stat, REAL4 *scalar, const SFTVector *sftvect1, const SFTVector *sftvect2)
{
  UINT4 i;
  REAL8 prod;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  ASSERT ( scalar, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sftvect1, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sftvect2, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sftvect1->length == sftvect2->length, stat, MAKEFAKEDATAC_EINCOMPAT, MAKEFAKEDATAC_MSGEINCOMPAT);

  prod = 0;
  for (i=0; i < sftvect1->length; i++)
    {
      REAL4 xy;
      TRY ( scalarProductSFT (stat->statusPtr, &xy, &(sftvect1->data[i]), &(sftvect2->data[i]) ), stat);
      prod += (REAL8) xy;
    }

  (*scalar) = (REAL4) prod;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* scalarProductSFTVector() */


/*--------------------------------------------------
 * calculate the difference of two SFT-vectors
 *--------------------------------------------------*/
void
subtractSFTVectors (LALStatus *stat, SFTVector **ret, const SFTVector *sftvect1, const SFTVector *sftvect2)
{
  UINT4 alpha, j;
  UINT4 M, N;
  SFTVector *vect = NULL;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  ASSERT ( ret, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( *ret == NULL, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sftvect1, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sftvect2, stat, MAKEFAKEDATAC_ENULL, MAKEFAKEDATAC_MSGENULL);
  ASSERT ( sftvect1->length == sftvect2->length, stat, MAKEFAKEDATAC_EINCOMPAT, MAKEFAKEDATAC_MSGEINCOMPAT);

  M = sftvect1->length;
  N = sftvect1->data[0].data->length;

  TRY ( LALCreateSFTVector (stat->statusPtr, &vect, M, N), stat);

  for (alpha = 0; alpha < M; alpha ++)
    {
      for (j=0; j < N; j++)
	{
	  vect->data[alpha].data->data[j].realf_FIXME = crealf(sftvect1->data[alpha].data->data[j]) - crealf(sftvect2->data[alpha].data->data[j]);
	  vect->data[alpha].data->data[j].imagf_FIXME = cimagf(sftvect1->data[alpha].data->data[j]) - cimagf(sftvect2->data[alpha].data->data[j]);
	} /* for j < N */

    } /* for alpha < M */

  (*ret) = vect;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* subtractSFTVectors() */
