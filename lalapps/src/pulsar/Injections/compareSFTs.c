 /*-----------------------------------------------------------------------
 *
 * File Name: compareSFTs.c
 *
 * Authors: Prix, R.
 *
 * Revision: $Id$
 *           
 *-----------------------------------------------------------------------
 */

/* ---------- includes ---------- */
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>

RCSID ("$Id");

/* Error codes and messages */
/* <lalErrTable file="MAKEFAKEDATACErrorTable"> */
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

/* </lalErrTable> */

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

/* local prototypes */
/* Prototypes for the functions defined in this file */
void initUserVars (LALStatus *stat);
REAL4 mymax (REAL4 x, REAL4 y);
void compare_SFTs (const SFTtype *sft1, const SFTtype *sft2);
void scalarProductSFT (LALStatus *stat, REAL4 *scalar, const SFTtype *sft1, const SFTtype *sft2);
void scalarProductSFTVector (LALStatus *stat, REAL4 *scalar, const SFTVector *sftvect1, const SFTVector *sftvect2);
void subtractSFTVectors (LALStatus *stat, SFTVector **ret, const SFTVector *sftvect1, const SFTVector *sftvect2);

extern int vrbflg;

/*----------------------------------------------------------------------*/
static const LALStatus empty_status;
/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_sftBname1;
CHAR *uvar_sftBname2;
INT4 uvar_debug;
BOOLEAN uvar_verbose;
BOOLEAN uvar_help;

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  SFTVector *SFTs1 = NULL, *SFTs2 = NULL;
  SFTVector *diffs = NULL;
  UINT4 i;

  lalDebugLevel = 0;
  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'd'), &status);
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* now read in the two complete sft-vectors */
  LAL_CALL (LALReadSFTfiles (&status, &SFTs1, 0, 0, 0, uvar_sftBname1), &status);
  LAL_CALL (LALReadSFTfiles (&status, &SFTs2, 0, 0, 0, uvar_sftBname2), &status);

  if (SFTs1->length != SFTs2->length) {
    LALPrintError ("Warning: number of SFTs differ for SFTbname1 and SFTbname2!\n");
  }

#if 0 
  /* OBSOLETE STUFF */
  printf ("\n******************** OLD METHOD ********************\n");

  printf ("\nCompareSFTs: Legend: 'max' := maxERR / max, 'mean' := meanERR / mean \n");
  for (i=0; i < mymax(SFTs1->length, SFTs2->length); i++)
    {
      printf ("i=%02d: ", i);
      compare_SFTs ( &(SFTs1->data[i]), &(SFTs2->data[i]) );
    }
  printf ("\n******************** NEW METHOD ********************\n");
#endif


  LAL_CALL (subtractSFTVectors (&status, &diffs, SFTs1, SFTs2), &status);

  if ( uvar_verbose)
    {
      for (i=0; i < mymax(SFTs1->length, SFTs2->length); i++)
	{
	  REAL4 d1, d2, d3;
	  REAL4 scalar, norm1, norm2, normdiff;
	  printf ("i=%02d: ", i);
	  LAL_CALL( scalarProductSFT(&status, &norm1, &(SFTs1->data[i]), &(SFTs1->data[i]) ), &status);
	  norm1 = sqrt(norm1);
	  LAL_CALL( scalarProductSFT(&status, &norm2, &(SFTs2->data[i]), &(SFTs2->data[i]) ), &status);
	  norm2 = sqrt(norm2);
	  LAL_CALL( scalarProductSFT(&status, &scalar, &(SFTs1->data[i]), &(SFTs2->data[i]) ), &status);
	  
	  LAL_CALL( scalarProductSFT(&status, &normdiff, &(diffs->data[i]), &(diffs->data[i]) ), &status);
	  
	  d1 = (norm1 - norm2)/norm1;
	  d2 = 1.0 - scalar / (norm1*norm2);
	  d3 = normdiff / (norm1*norm1 + norm2*norm2 );
	  
	  printf (" (|x|-|y|)/|x| = %10.3e, 1 - x.y/(|x| |y|) = %10.3e, |x - y|^2/(|x|^2+|y|^2)) = %10.3e\n", d1, d2, d3);
	} /* for i < SFTs->length */
    } /* if verbose */

  /* COMBINED measures */
  {
    REAL4 ret;
    REAL8 norm1, norm2, normdiff, scalar;
    REAL8 d1, d2, d3;


    LAL_CALL( scalarProductSFTVector (&status, &ret, SFTs1, SFTs1 ), &status);
    norm1 = sqrt( (REAL8)ret );

    LAL_CALL( scalarProductSFTVector (&status, &ret, SFTs2, SFTs2 ), &status);
    norm2 = sqrt( (REAL8)ret );

    LAL_CALL( scalarProductSFTVector (&status, &ret, SFTs1, SFTs2 ), &status);
    scalar = (REAL8) ret;

    LAL_CALL( scalarProductSFTVector(&status, &ret, diffs, diffs ), &status);
    normdiff = (REAL8) ret;

    d1 = (norm1 - norm2)/norm1;
    d2 = 1.0 - scalar / (norm1*norm2);
    d3 = normdiff / ( norm1*norm1 + norm2*norm2);

    if ( uvar_verbose )
      printf ("\nTOTAL: (|x|-|y|)/|x| = %10.3e, 1 - x.y/(|x| |y|) = %10.3e, |x - y|^2/(|x|^2+|y|^2) = %10.3e\n", d1, d2, d3);
    else
      printf ("%10.3e", d3);

  } /* combined total measures */


  /* free memory */
  LAL_CALL (LALDestroySFTVector(&status, &SFTs1), &status);
  LAL_CALL (LALDestroySFTVector(&status, &SFTs2), &status);
  LAL_CALL (LALDestroyUserVars (&status), &status);


  LALCheckMemoryLeaks(); 

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *stat)
{
  INITSTATUS( stat, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* set some defaults */
  uvar_debug = lalDebugLevel;
  uvar_verbose = FALSE;

  /* now register all our user-variable */

  LALregSTRINGUserVar(stat, sftBname1,	'1', UVAR_REQUIRED, "Path and basefilename for SFTs1");
  LALregSTRINGUserVar(stat, sftBname2,	'2', UVAR_REQUIRED, "Path and basefilename for SFTs2");
  LALregBOOLUserVar(stat,   verbose,	'v', UVAR_OPTIONAL, "Verbose output of differences");
  LALregBOOLUserVar(stat,   help,	'h', UVAR_HELP,     "Print this help/usage message");

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* initUserVars() */


REAL4 mymax (REAL4 x, REAL4 y)
{
  return (x > y ? x : y);
}
/* little debug-function: compare two sft's and print relative errors */
void
compare_SFTs (const SFTtype *sft1, const SFTtype *sft2)
{
  static LALStatus status;
  UINT4 i;
  REAL4 errpow= 0, errph = 0;
  REAL4 meanerr = 0, meanph = 0;
  REAL4 re1, re2, im1, im2, pow1, pow2, phase1, phase2;
  REAL8 Tdiff, tmp;
  INT4 ind_pow = 0, ind_phase = 0;
  REAL4 mean = 0, sftmax=0;

  if (sft1->data->length != sft2->data->length) 
    {
      LALPrintError ("\ncompare_SFTs(): lengths differ! %d != %d\n", sft1->data->length, sft2->data->length);
      return;
    }
  LALDeltaFloatGPS (&status, &Tdiff, &(sft1->epoch), &(sft2->epoch));
  if ( Tdiff != 0.0 ) 
    LALPrintError ("epochs differ: (%d s, %d ns)  vs (%d s, %d ns)\n", 
	    sft1->epoch.gpsSeconds, sft1->epoch.gpsNanoSeconds, sft2->epoch.gpsSeconds, sft2->epoch.gpsNanoSeconds);

  if ( sft1->f0 != sft2->f0)
    LALPrintError ("fmin differ: %fHz vs %fHz\n", sft1->f0, sft2->f0);

  if ( sft1->deltaF != sft2->deltaF )
    LALPrintError ("deltaF differs: %fHz vs %fHz\n", sft1->deltaF, sft2->deltaF);

  for (i=0; i < sft1->data->length; i++)
    {
      re1 = sft1->data->data[i].re;
      im1 = sft1->data->data[i].im;
      re2 = sft2->data->data[i].re;
      im2 = sft2->data->data[i].im;

      pow1 = sqrt(re1*re1 + im1*im1);
      pow2 = sqrt(re2*re2 + im2*im2);

      mean += pow1;
      sftmax = mymax (pow1, sftmax);

      phase1 = atan2 (im1, re1);
      phase2 = atan2 (im2, re2);
      /* error in power */
      tmp = fabs(pow1-pow2);
      meanerr += tmp;
      if (tmp > errpow) 
	{
	  errpow = tmp;
	  ind_pow = i;
	}
      /* error in phase */
      tmp = fabs(phase1 - phase2);
      if (tmp > LAL_PI)
	tmp = LAL_TWOPI - tmp;

      meanph += tmp;
      if (tmp > errph)
	{
	  errph = tmp;
	  ind_phase = i;
	}

    } /* for i */

  mean /= sft1->data->length;
  meanerr /= sft1->data->length;
  meanph /= sft1->data->length;

  printf ("POWER: max = %e @ %6.2fHz  mean = %e, PHASE: max = %e @ %6.2fHz, mean = %e \n", 
	  errpow/sftmax, sft1->f0 + ind_pow * sft1->deltaF, meanerr/mean, 
	  (REAL4)(errph), sft2->f0 + ind_phase * sft1->deltaF, (REAL4)(meanph));

  return;

} /* compare_SFTs() */

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

  INITSTATUS( stat, "scalarProductSFT", rcsid );
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
      xre = (REAL8)sft1->data->data[i].re;
      xim = (REAL8)sft1->data->data[i].im;
      yre = (REAL8)sft2->data->data[i].re;
      yim = (REAL8)sft2->data->data[i].im;

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

  INITSTATUS( stat, "scalarProductSFT", rcsid );
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

  INITSTATUS( stat, "subtractSFTVector", rcsid );
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
	  vect->data[alpha].data->data[j].re = sftvect1->data[alpha].data->data[j].re - sftvect2->data[alpha].data->data[j].re;
	  vect->data[alpha].data->data[j].im = sftvect1->data[alpha].data->data[j].im - sftvect2->data[alpha].data->data[j].im;
	} /* for j < N */

    } /* for alpha < M */

  (*ret) = vect;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* subtractSFTVectors() */
