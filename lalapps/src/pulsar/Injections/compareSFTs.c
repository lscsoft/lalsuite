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


#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "File IO error"
#define MAKEFAKEDATAC_MSGENOARG "Missing argument"


/* </lalErrTable> */

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

/* local prototypes */
/* Prototypes for the functions defined in this file */
void initUserVars (LALStatus *stat);
REAL4 mymax (REAL4 x, REAL4 y);
void compare_SFTs (const SFTtype *sft1, const SFTtype *sft2);

extern int vrbflg;

/*----------------------------------------------------------------------*/
static const LALStatus empty_status;
/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_sftBname1;
CHAR *uvar_sftBname2;
INT4 uvar_debug;
BOOLEAN uvar_help;

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  SFTVector *SFTs1 = NULL, *SFTs2 = NULL;
  UINT4 i;

  lalDebugLevel = 0;
  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);
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

  printf ("\nCompareSFTs: Legend: 'max' := maxERR / max, 'mean' := meanERR / mean \n");
  for (i=0; i < mymax(SFTs1->length, SFTs2->length); i++)
    {
      printf ("i=%02d: ", i);
      compare_SFTs ( &(SFTs1->data[i]), &(SFTs2->data[i]) );
    }


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

  /* now register all our user-variable */

  LALregSTRINGUserVar(stat, sftBname1,	'1', UVAR_REQUIRED, "Path and basefilename for SFTs1");
  LALregSTRINGUserVar(stat, sftBname2,	'2', UVAR_REQUIRED, "Path and basefilename for SFTs2");
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
      printf ("\ncompare_SFTs(): lengths differ! %d != %d\n", sft1->data->length, sft2->data->length);
      return;
    }
  LALDeltaFloatGPS (&status, &Tdiff, &(sft1->epoch), &(sft2->epoch));
  if ( Tdiff != 0.0 ) 
    printf ("epochs differ: (%d s, %d ns)  vs (%d s, %d ns)\n", 
	    sft1->epoch.gpsSeconds, sft1->epoch.gpsNanoSeconds, sft2->epoch.gpsSeconds, sft2->epoch.gpsNanoSeconds);

  if ( sft1->f0 != sft2->f0)
    printf ("fmin differ: %fHz vs %fHz\n", sft1->f0, sft2->f0);

  if ( sft1->deltaF != sft2->deltaF )
    printf ("deltaF differs: %fHz vs %fHz\n", sft1->deltaF, sft2->deltaF);

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

