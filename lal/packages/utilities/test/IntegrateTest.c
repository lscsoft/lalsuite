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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/Integrate.h>

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


static void f1 (LALStatus *s, REAL4 *y, REAL4 x, void *p)
{
  REAL4  x2 = x*x;
  REAL4  x4 = x2*x2;
  INT4  *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x4*log(x + sqrt(x2 + 1));
  RETURN (s);
}

static void ff1 (LALStatus *s, REAL8 *y, REAL8 x, void *p)
{
  REAL8  x2 = x*x;
  REAL8  x4 = x2*x2;
  INT4  *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x4*log(x + sqrt(x2 + 1));
  RETURN (s);
}

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

static void f2 (LALStatus *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/(x*x*x);
  RETURN (s);
}

static void ff2 (LALStatus *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/(x*x*x);
  RETURN (s);
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

static void f3 (LALStatus *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = exp(-x*x/2);
  RETURN (s);
}

static void ff3 (LALStatus *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = exp(-x*x/2);
  RETURN (s);
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

static void f4 (LALStatus *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/sqrt(x);
  RETURN (s);
}

static void ff4 (LALStatus *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/sqrt(x);
  RETURN (s);
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

static void f5 (LALStatus *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x + 1/sqrt(5 - x);
  RETURN (s);
}

static void ff5 (LALStatus *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x + 1/sqrt(5 - x);
  RETURN (s);
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

static void g (LALStatus *s, REAL4 *z, REAL4 x, void *p)
{
  REAL4 y;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  y  = *((REAL4 *)p);
  *z = exp(-(x*x + y*y)/2);
  RETURN (s);
}

static void gg (LALStatus *s, REAL8 *z, REAL8 x, void *p)
{
  REAL8 y;
  INITSTATUS(s);
  ASSERT (p, s, 1, "Null pointer");
  y  = *((REAL8 *)p);
  *z = exp(-(x*x + y*y)/2);
  RETURN (s);
}

static void h (LALStatus *s, REAL4 *z, REAL4 y, void *p)
{
  SIntegrateIn intinp;
  INITSTATUS(s);
  ATTATCHSTATUSPTR (s);
  if ( p )
    ABORT ( s, 2, "Non-null pointer");
  p = NULL;
  intinp.function = g;
  intinp.xmin     = 0;
  intinp.xmax     = sqrt(y);
  intinp.type     = ClosedInterval;
  LALSRombergIntegrate (s->statusPtr, z, &intinp, &y);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
}

static void hh (LALStatus *s, REAL8 *z, REAL8 y, void *p)
{
  DIntegrateIn intinp;
  INITSTATUS(s);
  ATTATCHSTATUSPTR (s);
  if ( p )
    ABORT ( s, 2, "Non-null pointer");
  p = NULL;
  intinp.function = gg;
  intinp.xmin     = 0;
  intinp.xmax     = sqrt(y);
  intinp.type     = ClosedInterval;
  LALDRombergIntegrate (s->statusPtr, z, &intinp, &y);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
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
static void bad (LALStatus UNUSED *s, REAL4 *y, REAL4 UNUSED x, void *p)
{
  INT4 *n = (INT4 *)p;
  *y = *n = 1664525L*(*n) + 1013904223L;
}

static void bbad (LALStatus UNUSED *s, REAL8 *y, REAL8 UNUSED x, void *p)
{
  *y = (REAL8)(++(*(INT4 *)p));
}

static REAL8 xbbad (REAL8 UNUSED x, void *p)
{
  return (REAL8)(++(*(INT4 *)p));
}
#endif



extern char *optarg;
extern int   optind;

int   verbose    = 0;

static void Usage (const char *program, int exitflag);

static void ParseOptions (int argc, char *argv[]);

static void TestStatus (LALStatus *status, const char *expectCodes, int exitCode);

#if defined(NDEBUG) || defined(LAL_NDEBUG)
/* debugging is turned off */
#else
static void ClearStatus (LALStatus *status);
#endif

int main (int argc, char *argv[])
{
  const REAL4   sepsilon = 1e-6;
  const REAL8   depsilon = 1e-13; /* not as good as expected (1e-15) */
  static LALStatus status;
  SIntegrateIn  sintinp;
  DIntegrateIn  dintinp;
  REAL4         sresult;
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

  sintinp.function = f1;
  sintinp.xmin     = 0;
  sintinp.xmax     = 2;
  sintinp.type     = ClosedInterval;
  dintinp.function = ff1;
  dintinp.xmin     = 0;
  dintinp.xmax     = 2;
  dintinp.type     = ClosedInterval;

  count  = 0;
  expect = 8.153364119811650205L;
  LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", sresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(sresult - expect) > sepsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count = 0;
  LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
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

  count = 0;
  dresult = XLALREAL8RombergIntegrate (&xff1, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
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
  sintinp.function = f2;
  sintinp.xmin     = 10;
  sintinp.xmax     = 1e30;
  sintinp.type     = InfiniteDomainPow;
  dintinp.function = ff2;
  dintinp.xmin     = 10;
  dintinp.xmax     = 1e300;
  dintinp.type     = InfiniteDomainPow;

  count  = 0;
  expect = 1.0L/200.0L;
  LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", sresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(sresult - expect) > sepsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count = 0;
  LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
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

  count = 0;
  dresult = XLALREAL8RombergIntegrate (&xff2, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
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
  sintinp.function = f3;
  sintinp.xmin     = 2;
  sintinp.xmax     = 1e30;
  sintinp.type     = InfiniteDomainExp;
  dintinp.function = ff3;
  dintinp.xmin     = 2;
  dintinp.xmax     = 1e300;
  dintinp.type     = InfiniteDomainExp;

  count  = 0;
  expect = 0.0570261239928920483L;
  LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", sresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(sresult - expect) > sepsilon*fabsl(expect))
  {
    if ( verbose )
      fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count = 0;
  LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
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

  count = 0;
  dresult = XLALREAL8RombergIntegrate (&xff3, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
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
  sintinp.function = f4;
  sintinp.xmin     = 0;
  sintinp.xmax     = 1;
  sintinp.type     = SingularLowerLimit;
  dintinp.function = ff4;
  dintinp.xmin     = 0;
  dintinp.xmax     = 1;
  dintinp.type     = SingularLowerLimit;

  count  = 0;
  expect = 2.0L;
  LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", sresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  if (fabsl(sresult - expect) > sepsilon*fabsl(expect))
  {
  if ( verbose )
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count  = 0;
  LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
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

  count = 0;
  dresult = XLALREAL8RombergIntegrate (&xff4, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
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
  sintinp.function = f5;
  sintinp.xmin     = 4;
  sintinp.xmax     = 5;
  sintinp.type     = SingularUpperLimit;
  dintinp.function = ff5;
  dintinp.xmin     = 4;
  dintinp.xmax     = 5;
  dintinp.type     = SingularUpperLimit;

  count  = 0;
  expect = 6.5L;
  LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", sresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  /* this doesn't work so well: multiply tolerance by factor of three */
  if (fabsl(sresult - expect) > 3*sepsilon*fabsl(expect))
  {
    if ( verbose )
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count  = 0;
  LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("number of function calls: %d\n", count);
  if ( verbose )
    printf ("result: %.15f\n", dresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  /* this doesn't work so well: multiply tolerance by factor of three */
  if (fabsl(dresult - expect) > 3*depsilon*fabsl(expect))
  {
    if ( verbose )
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count = 0;
  dresult = XLALREAL8RombergIntegrate (&xff5, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
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
   * Test 6: Two-dimensional integral.
   *
   */


  if ( verbose )
  printf ("\nTest 6: Two-dimensional integral.\n");
  sintinp.function = h;
  sintinp.xmin     = 0;
  sintinp.xmax     = 10;
  sintinp.type     = OpenInterval;
  dintinp.function = hh;
  dintinp.xmin     = 0;
  dintinp.xmax     = 10;
  dintinp.type     = OpenInterval;

  expect = 0.88274109326014810823L;
  LALSRombergIntegrate (&status, &sresult, &sintinp, NULL);
  TestStatus (&status, CODES(0), 1);
  if ( verbose )
    printf ("result: %.15f\n", sresult);
  if ( verbose )
    printf ("expect: %.15Lf\n", expect);
  /* integral isn't very accurate because we needed to use an open interval */
  if ( verbose )
    printf ("error:  %.2f%%\n", 100*fabsl(sresult - expect)/fabsl(expect));
  /*
   * don't do 2d double-precision: it takes too long!
   *
   * LALDRombergIntegrate (&status, &dresult, &dintinp, NULL);
   * TestStatus (&status, CODES(0), 1);
   * if ( verbose )
   * printf ("result: %.15f\n", dresult);
   * if ( verbose )
   * printf ("expect: %.15Lf\n", expect);
   * if ( verbose )
   * printf ("error:  %.2f%%\n", 100*fabs(dresult - expect)/fabs(expect));
   *
   */


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
      printf ("\nNull pointer:\r");
    LALSRombergIntegrate (&status, NULL, &sintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_ENULL), 1);
    LALDRombergIntegrate (&status, NULL, &dintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_ENULL), 1);
    if ( verbose )
      printf ("Null pointer check passed.\n");

    if ( verbose )
      printf ("\nNull pointer:\r");
    LALSRombergIntegrate (&status, &sresult, NULL, &count);
    TestStatus (&status, CODES(INTEGRATEH_ENULL), 1);
    LALDRombergIntegrate (&status, &dresult, NULL, &count);
    TestStatus (&status, CODES(INTEGRATEH_ENULL), 1);
    if ( verbose )
      printf ("Null pointer check passed.\n");

    if ( verbose )
      printf ("\nNull pointer:\r");
    sintinp.function = NULL;
    sintinp.xmin     = 0;
    sintinp.xmax     = 2;
    sintinp.type     = ClosedInterval;
    dintinp.function = NULL;
    dintinp.xmin     = 0;
    dintinp.xmax     = 2;
    dintinp.type     = ClosedInterval;
    LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_ENULL), 1);
    LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_ENULL), 1);
    dresult = XLALREAL8RombergIntegrate (NULL, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
    if (xlalErrno == XLAL_EFAULT)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Null pointer check passed.\n");

    if ( verbose )
      printf ("\nInvalid domain:\r");
    sintinp.function = f1;
    sintinp.xmin     = 0;
    sintinp.xmax     = 0;
    sintinp.type     = ClosedInterval;
    dintinp.function = ff1;
    dintinp.xmin     = 0;
    dintinp.xmax     = 0;
    dintinp.type     = ClosedInterval;
    LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_EIDOM), 1);
    LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_EIDOM), 1);
    dresult = XLALREAL8RombergIntegrate (&xff1, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
    if (xlalErrno == XLAL_EDOM)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Invalid domain check passed.\n");

    if ( verbose )
      printf ("\nUnknown integral type:\r");
    sintinp.function = f1;
    sintinp.xmin     = 0;
    sintinp.xmax     = 2;
    sintinp.type     = 999;
    dintinp.function = ff1;
    dintinp.xmin     = 0;
    dintinp.xmax     = 2;
    dintinp.type     = 999;
    LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_ETYPE), 1);
    LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_ETYPE), 1);
    dresult = XLALREAL8RombergIntegrate (&xff1, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
    if (xlalErrno == XLAL_EINVAL)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Unknown integral type check passed.\n");

    if ( verbose )
      printf ("\nMaximum iterations exceeded:\r");
    sintinp.function = bad;  /* bad is a quick random number generator */
    sintinp.xmin     = 0;
    sintinp.xmax     = 2;
    sintinp.type     = ClosedInterval;
    dintinp.function = bbad;  /* bbad is a quick random number generator */
    dintinp.xmin     = 0;
    dintinp.xmax     = 2;
    dintinp.type     = ClosedInterval;
    count            = 13;   /* count is now used as a random number seed */
    LALSRombergIntegrate (&status, &sresult, &sintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_EMXIT), 1);
    count = 1;
    LALDRombergIntegrate (&status, &dresult, &dintinp, &count);
    TestStatus (&status, CODES(INTEGRATEH_EMXIT), 1);
    dresult = XLALREAL8RombergIntegrate (&xbbad, &count, dintinp.xmin, dintinp.xmax, dintinp.type);
    if (xlalErrno == XLAL_EMAXITER)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Maximum iterations exceeded check passed.\n");

    if ( verbose )
      printf ("\nRecursive error:\r");
    sintinp.function = f1;
    sintinp.xmin     = 0;
    sintinp.xmax     = 2;
    sintinp.type     = ClosedInterval;
    dintinp.function = ff1;
    dintinp.xmin     = 0;
    dintinp.xmax     = 2;
    dintinp.type     = ClosedInterval;
    LALSRombergIntegrate (&status, &sresult, &sintinp, NULL);
    TestStatus (&status, CODES(-1), 1);
    ClearStatus (&status);
    LALDRombergIntegrate (&status, &dresult, &dintinp, NULL);
    TestStatus (&status, CODES(-1), 1);
    dresult = XLALREAL8RombergIntegrate (&xff1, NULL, dintinp.xmin, dintinp.xmax, dintinp.type);
    if (xlalErrno == XLAL_EFUNC + XLAL_EFAULT)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Recursive error check passed.\n");
    ClearStatus (&status);
  }

#endif

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
    REPORTSTATUS (status);
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


#if defined(NDEBUG) || defined(LAL_NDEBUG)
/* debugging is turned off */
#else
/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
static void
ClearStatus (LALStatus *status)
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}
#endif

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

    c = getopt (argc, argv, "hqvd:");
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

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}

/** \endcond */
