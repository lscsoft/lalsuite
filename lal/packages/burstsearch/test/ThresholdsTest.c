/*-----------------------------------------------------------------------
 *
 * File Name: ThresholdsTest.c
 *
 * Author: Eanna Flanagan
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * main()
 *
 * SYNOPSIS
 *
 * DESCRIPTION
 * Test suite for functions in Thresholds.c
 *
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 *
 * CALLS
 * LALOverlap()
 * LALSCreateVector()
 * LALSDestroyVector()
 * FindRoot()
 *
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/Thresholds.h>


#define _CODES(x) #x
#define CODES(x) _CODES(x)


NRCSID (MAIN, "$Id$");


extern char *optarg;
extern int   optind;

INT4 lalDebugLevel = 0;   /* set to 2 to get full status information for tests */
INT4 verbose    = 1;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);



int
main (int argc, char *argv[])
{
  static LALStatus    status;
  ChisqCdfIn       input1;
  Chi2ThresholdIn  input2;
  RhoThresholdIn   input3;

  REAL8 chi2;
  REAL8 dof;
  REAL8 rho;
  REAL8 alpha;
  REAL8 alpha1;
  REAL8 alpha2;
  REAL8 beta;
  REAL8 temp1;
  REAL8 temp2;

  /*
   *
   * Parse the command line options
   *
   */


  ParseOptions (argc, argv);

  /* 
   *  Check to make sure the functions return the correct values.  
   *  First time around 
   */

  chi2 = 2.3l;
  dof = 8.0;
  rho = 2.2;
  input1.chi2 = chi2;
  input1.dof = dof;
  input1.nonCentral = rho*rho;


  LALChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(0), 1);
  LALOneMinusChisqCdf (&status, &alpha2, &input1);
  TestStatus (&status, CODES(0), 1);
  LALNoncChisqCdf( &status, &beta, &input1);
  TestStatus (&status, CODES(0), 1);

  alpha = 1.0 - alpha1;
  if(verbose)
    {
      printf("-- Test 1 of LALChisqCdf(), LALNoncChisqCdf(), LALChi2Threshold() and LALRhoThreshold() --\n");
      printf("chi2 is %f and dof is %f\n",chi2,dof);
      printf("rho is %f\n",rho);
      printf("alpha is %f, should be 0.970406\n",alpha);
      printf("alpha is %f, should be 0.970406\n",alpha2);
      printf("beta is %f, should be 0.00439452\n",beta);
    }

  /* test were the correct values obtained */
  if( (fabs(alpha - 0.970406) > 1e-5) || (fabs(beta - 0.00439452) > 1e-7))
    {
      fprintf( stderr, "Incorrect values returned\n");
      return 1;
    }
      
  input2.falseAlarm = alpha;
  input2.dof = dof;
  LALChi2Threshold( &status, &temp1, &input2);
  TestStatus (&status, CODES(0), 1);
  input3.chi2 = chi2;
  input3.dof = dof;
  input3.falseDismissal=beta;
  LALRhoThreshold( &status, &temp2, &input3);
  TestStatus (&status, CODES(0), 1);

  if(verbose)
    {
      printf("Working backwards, chi2 is %f\n",temp1);
      printf("Working backwards, rho is %f\n",temp2);      
    }
  
  
  /* test were the correct values obtained */
  if( (fabs(temp1-chi2) > 1e-5*chi2) || (fabs(temp2 - rho) > 1e-5*rho))
    {
      fprintf( stderr, "Incorrect values returned\n");
      return 1;
    }

  /* 
   *  Check to make sure the functions return the correct values.  
   *  Second time around 
   */
  chi2 = 12.3l;
  dof = 3.1;
  rho = 2.2;
  input1.chi2 = chi2;
  input1.dof = dof;
  input1.nonCentral = rho*rho;

  LALChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(0), 1);
  LALNoncChisqCdf( &status, &beta, &input1);
  TestStatus (&status, CODES(0), 1);

  alpha = 1.0 - alpha1;
  if(verbose)
    {
      printf("-- Test 2 of LALChisqCdf(), LALNoncChisqCdf(), LALChi2Threshold() and LALRhoThreshold() --\n");
      printf("chi2 is %f and dof is %f\n",chi2,dof);
      printf("rho is %f\n",rho);
      printf("alpha is %f, should be 0.007066\n",alpha);
      printf("beta is %f, should be 0.822575\n",beta);
    }

  /* test were the correct values obtained */
  if( (fabs(alpha - 0.007066) > 1e-7) || (fabs(beta - 0.822575) > 1e-6))
    {
      fprintf( stderr, "Incorrect values returned\n");
      return 1;
    }
      
  input2.falseAlarm = alpha;
  input2.dof = dof;
  LALChi2Threshold( &status, &temp1, &input2);
  TestStatus (&status, CODES(0), 1);
  input3.chi2 = chi2;
  input3.dof = dof;
  input3.falseDismissal=beta;
  LALRhoThreshold( &status, &temp2, &input3);
  TestStatus (&status, CODES(0), 1);

  if(verbose)
    {
      printf("Working backwards, chi2 is %f\n",temp1);
      printf("Working backwards, rho is %f\n",temp2);      
    }
  
  /* test were the correct values obtained */
  if( (fabs(temp1-chi2) > 1e-5*chi2) || (fabs(temp2 - rho) > 1e-5*rho))
    {
      fprintf( stderr, "Incorrect values returned\n");
      return 1;
    }



  /*
   *
   * Check to make sure that correct error codes are generated.
   *
   */

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {

  if (verbose || lalDebugLevel)
  {
    printf ("\n===== Check Errors =====\n");
  }


  /* one of the arguments is a null pointer */

  if (verbose)
  {
    printf ("\n----- Null Pointer Error: Code 1 (8 times) \n");
  }

  LALChisqCdf (&status, &alpha1, NULL);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALChisqCdf (&status, NULL, &input1);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALNoncChisqCdf (&status, &alpha1, NULL);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALNoncChisqCdf (&status, NULL, &input1);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALChi2Threshold( &status, &temp1, NULL);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALChi2Threshold( &status, NULL, &input2);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALRhoThreshold( &status, &temp2, NULL);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);

  LALRhoThreshold( &status, NULL, &input3);
  TestStatus (&status, CODES(LAL_NULL_ERR), 1);




  /* Arguments are non-postive */

  if (verbose)
  {
    printf ("\n----- Negative arguments Error: Code 2 (8 times) \n");
  }

  input1.chi2 *= -1.0;
  LALChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  LALNoncChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input1.chi2 *= -1.0;  /* set it back to positive for remaining tests */

  input1.dof *= -1.0;
  LALChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  LALNoncChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input1.dof *= -1.0;

  input1.nonCentral *= -1.0;
  LALNoncChisqCdf (&status, &alpha1, &input1);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input1.nonCentral *= -1.0;

  input2.dof *= -1;
  LALChi2Threshold( &status, &temp1, &input2);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input2.dof *= -1;
  
  input3.dof *= -1;
  LALRhoThreshold( &status, &temp2, &input3);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input3.dof *= -1;

  input3.chi2 *= -1;
  LALRhoThreshold( &status, &temp2, &input3);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input3.chi2 *= -1;




  /*
   * 
   *  Maximum iterations exceeded error: this test no longer performed
   *  since GSL is used to compute these distributions
   *
   */



  /* 
   *  Test LALNoncChisqCdf for a recursive error 
   *
   *  This test currently commented out since it causes a memory
   *  leak and causes the LALCheckMemoryLeaks() test below to fail.
   *  The memory leak is due to the way the LAL LALStatus macros are
   *  currently written.  
   *
   */

  /*  LALNoncChisqCdf (&status, &alpha1, &input1);
      TestStatus (&status, CODES(-1), 1); */




  /* reset parameters to original values */
  input1.dof = dof;
  input1.chi2 = chi2;
  
  /* 
   *  There is no test here for exceeding the maximum number of
   *  iterations in LALNoncChisqCdf() since I could not find a set
   *   of parameters which caused this condition to occur.
   */






  /* Supplied probabilities must lie between 0 and 1 */

  if (verbose)
  {
    printf ("\n----- Probability out of range Error: Code 8 (4 times) \n");
  }


  input2.falseAlarm = -1.0;
  LALChi2Threshold (&status, &temp1, &input2);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input2.falseAlarm = 2.0;
  LALChi2Threshold (&status, &alpha1, &input2);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  /* set it back to original value for remaining tests */
  input2.falseAlarm= alpha;  



  input3.falseDismissal = -1.0;
  LALRhoThreshold (&status, &temp2, &input3);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  input3.falseDismissal = 2.0;
  LALRhoThreshold (&status, &temp2, &input3);
  TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
  /* set it back to original value for remaining tests */
  input3.falseDismissal= beta;  

  }
#endif

  LALCheckMemoryLeaks ();

  if(verbose)  printf("PASS: all tests\n");

  return 0;
}




/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (LALStatus *status, const char *ignored, int exitcode)
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    /* REPORTSTATUS (status);*/
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -q         quiet: run silently\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  fprintf (stderr, "  -d level   set lalDebugLevel to level\n");
  exit (exitcode);
}


/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}
