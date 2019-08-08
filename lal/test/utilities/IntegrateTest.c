/*
*  Copyright (C) 2007 Jolien Creighton
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
 * \ingroup Integrate_h
 *
 * \brief Tests the routines in \ref Integrate_h by performing a suite of numerical
 * integrations and checking the accuracy of the results.
 *
 * ### Usage ###
 *
 * \code
 * IntegrateTest [options]
 * Options:
 * -h         print this message
 * -q         quiet: run silently
 * -v         verbose: print extra information
 * -d level   set lalDebugLevel to level
 * \endcode
 *
 * ### Exit codes ###
 *
 * <table><tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 */

/** \cond DONT_DOXYGEN */
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/Integrate.h>
#include <lal/LALString.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define CODES_(x) #x
#define CODES(x) CODES_(x)

/*
 *
 * These functions are used as integrands for the test integrations.
 *
 */


static REAL8 xff1 (REAL8 x, void *p)
{
  REAL8  y;
  REAL8  x2 = x*x;
  REAL8  x4 = x2*x2;
  INT4  *n;
  if (p == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  ++(*(n = (INT4 *)p));
  y = x4*log(x + sqrt(x2 + 1));
  return y;
}

static REAL8 xff2 (REAL8 x, void *p)
{
  REAL8  y;
  INT4  *n;
  if (p == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  ++(*(n = (INT4 *)p));
  y = 1/(x*x*x);
  return y;
}

static REAL8 xff3 (REAL8 x, void *p)
{
  REAL8  y;
  INT4  *n;
  if (p == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  ++(*(n = (INT4 *)p));
  y = exp(-x*x/2);
  return y;
}

static REAL8 xff4 (REAL8 x, void *p)
{
  REAL8  y;
  INT4  *n;
  if (p == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  ++(*(n = (INT4 *)p));
  y = 1/sqrt(x);
  return y;
}

static REAL8 xff5 (REAL8 x, void *p)
{
  REAL8  y;
  INT4  *n;
  if (p == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  ++(*(n = (INT4 *)p));
  y = x + 1/sqrt(5 - x);
  return y;
}

#if defined(NDEBUG) || defined(LAL_NDEBUG)
/* debugging is turned off */
#else
/*
 *
 * This function produces random numbers.  Integration of this function should
 * not converge.  Make this routine fast... no status handling!
 *
 */
static REAL8 xbbad (REAL8 UNUSED x, void *p)
{
  return (REAL8)(++(*(INT4 *)p));
}
#endif



int   verbose    = 0;

static void Usage (const char *program, int exitflag);

static void ParseOptions (int argc, char *argv[]);

int main (int argc, char *argv[])
{
  const REAL8   depsilon = 1e-13; /* not as good as expected (1e-15) */
  REAL8         dresult;
  long double   expect;
  INT4          count;

  ParseOptions (argc, argv);


  /*
   *
   * Test 1: Integrate a regular function over a closed interval.
   *
   */


  if ( verbose )
    printf ("Test 1:"
        " Integrate a regular function over a closed interval.\n");
  count = 0;
  expect = 8.153364119811650205L;
  dresult = XLALREAL8RombergIntegrate (&xff1, &count, 0, 2, ClosedInterval);
  if (xlalErrno)
    abort();
  else if (verbose)
    printf("\nXLALREAL8RombergIntegrate exitted with xlalErrno: %d\n", xlalErrno);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", dresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(dresult - expect) > depsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 2: Integrate to infinity a function with power-law fall-off.
   *
   */


  if ( verbose )
    printf ("\nTest 2:"
        " Integrate to infinity a function with power-law fall-off.\n");
  count = 0;
  expect = 1.0L/200.0L;
  dresult = XLALREAL8RombergIntegrate (&xff2, &count, 10, 1e300, InfiniteDomainPow);
  if (xlalErrno)
    abort();
  else if (verbose)
    printf("\nXLALREAL8RombergIntegrate exitted with xlalErrno: %d\n", xlalErrno);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", dresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(dresult - expect) > depsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 3: Integrate to infinity a function that falls off exponentially.
   *
   */


  if ( verbose )
    printf ("\nTest 3:"
        " Integrate to infinity a function that falls off exponentially.");
  count  = 0;
  expect = 0.0570261239928920483L;
  dresult = XLALREAL8RombergIntegrate (&xff3, &count, 2, 1e300, InfiniteDomainExp);
  if (xlalErrno)
    abort();
  else if (verbose)
    printf("\nXLALREAL8RombergIntegrate exitted with xlalErrno: %d\n", xlalErrno);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", dresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(dresult - expect) > depsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 4: Integrate an integrable singularity at the lower limit.
   *
   */


  if ( verbose )
    printf ("\nTest 4:"
        " Integrate an integrable singularity at the lower limit.\n");
  count  = 0;
  expect = 2.0L;
  dresult = XLALREAL8RombergIntegrate (&xff4, &count, 0, 1, SingularLowerLimit);
  if (xlalErrno)
    abort();
  else if (verbose)
    printf("\nXLALREAL8RombergIntegrate exitted with xlalErrno: %d\n", xlalErrno);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", dresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(dresult - expect) > depsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 5: Integrate an integrable singularity at the upper limit.
   *
   */


  if ( verbose )
    printf ("\nTest 5:"
        " Integrate an integrable singularity at the upper limit.\n");
  count  = 0;
  expect = 6.5L;
  dresult = XLALREAL8RombergIntegrate (&xff5, &count, 4, 5, SingularUpperLimit);
  if (xlalErrno)
    abort();
  else if (verbose)
    printf("\nXLALREAL8RombergIntegrate exitted with xlalErrno: %d\n", xlalErrno);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", dresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(dresult - expect) > depsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  LALCheckMemoryLeaks ();


  /*
   *
   * Check error conditions.
   *
   */

#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {
    if ( verbose )
      printf ("\nChecking error conditions:\n");

    if ( verbose )
      printf ("\nUnknown integral type:\r");
    dresult = XLALREAL8RombergIntegrate (&xff1, &count, 0, 2, 999);
    if (xlalErrno == XLAL_EINVAL)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Unknown integral type check passed.\n");

    if ( verbose )
      printf ("\nMaximum iterations exceeded:\r");
    dresult = XLALREAL8RombergIntegrate (&xbbad, &count, 0, 2, ClosedInterval);
    if (xlalErrno == XLAL_EMAXITER)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Maximum iterations exceeded check passed.\n");

    if ( verbose )
      printf ("\nRecursive error:\r");
    dresult = XLALREAL8RombergIntegrate (&xff1, NULL, 0, 2, ClosedInterval);
    if (xlalErrno == XLAL_EFUNC + XLAL_EFAULT)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Recursive error check passed.\n");
  }

#endif

  return 0;
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
  FILE *fp;

  while (1)
  {
    int c = -1;

    c = LALgetopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently */
        fp = freopen ("/dev/null", "w", stderr);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        fp = freopen ("/dev/null", "w", stdout);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (LALoptind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}

/** \endcond */
