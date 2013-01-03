/*
 * Copyright (C) 2009 Reinhard Prix
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

#include <math.h>
#include <sys/times.h>

#include <lal/LALMalloc.h>
#include <lal/LALInitBarycenter.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALConstants.h>

#include <lal/UniversalDopplerMetric.h>

/** \author Reinhard Prix
 * \file
 * \ingroup UniversalDopplerMetric_h
 * \brief Tests for exported functions in UniversalDopplerMetric
 *
 */

// ---------- defines --------------------
#define TEST_OK		0
#define TEST_FAIL	1
#define TEST_VAL	2

// ---------- Macros --------------------

// copy 3 components of Euklidean vector
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)
// compute norm of Euklidean 3-vector
#define NORM(v) ( sqrt ( (v)[0]*(v)[0] + (v)[1]*(v)[1] + (v)[2]*(v)[2] )  )
// compute scalar division of Euklidean 3-vector by div
#define DIV_VECT(dst,src,div) do { (dst)[0] = (src)[0]/(div); (dst)[1] = (src)[1]/(div); (dst)[2] = (src)[2]/(div); } while(0)
#define SUB_VECT(dst,src) do { (dst)[0] -= (src)[0]; (dst)[1] -= (src)[1]; (dst)[2] -= (src)[2]; } while(0)
#define MULT_VECT(v,lam) do{ (v)[0] *= (lam); (v)[1] *= (lam); (v)[2] *= (lam); } while(0)

// ---------- global variables --------------------
static LALStatus empty_status;

extern int lalDebugLevel;

// ---------- local prototypes
static int test_XLALComputeOrbitalDerivatives(void);

// ---------- function definitions --------------------
/** MAIN function: calls a number of unit-tests
 */
int main( void )
{
  const char *fn = __func__;
  INT4 ret;

  lalDebugLevel = 1;

  ret = test_XLALComputeOrbitalDerivatives();

  if ( ret != TEST_OK )
    {
      XLALPrintError ("%s: test_XLALComputeOrbitalDerivatives() failed.\n", fn );
      return ret;
    }

  /* check for memory leaks */
  LALCheckMemoryLeaks();


  /* all tests passed */
  return TEST_OK;

} /* main() */


/** Unit test function for XLALComputeOrbitalDerivatives()
 */
static int
test_XLALComputeOrbitalDerivatives ( void )
{
  const char *fn = __func__;

  // ----- load an example ephemeris, describing a pure cicular 2D
  // orbit w period of one year
  EphemerisData edat = empty_EphemerisData;
  LALStatus status = empty_status;
  CHAR earthfname[] = DATADIR "circularEphem.dat";
  CHAR sunfname[] = DATADIR "sun00-04.dat";

  printf ("\nEntering test_XLALComputeOrbitalDerivatives() ... ");

  edat.ephiles.earthEphemeris = earthfname;
  edat.ephiles.sunEphemeris   = sunfname;

  LALInitBarycenter( &status, &edat);
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALInitBarycenter() failed!\n\n", fn);
    return TEST_FAIL;
  }


  /* Unit test for XLALComputeOrbitalDerivatives() */
  // ----------------------------------------
  // the tests consists in using a purely 2D circular orbit to
  // compute the derivatives, so we can analytically check them!

  // one line from the middle of cicularEphem.dat:
  // 700244800    498.41220200278  24.31154155846971  0    -4.84039532365306e-06  9.923320107134072e-05  0    -1.975719726623028e-11 -9.63716218195963e-13 0

  // we can use this to check derivatives '0-'2'
  vect3Dlist_t *rOrb_n = NULL;
  LIGOTimeGPS t0 = { 700244800, 0 };
  vect3D_t r_0 = {498.41220200278, 24.31154155846971,  0 };
  vect3D_t r_1 = {-4.84039532365306e-06,   9.923320107134072e-05,  0 };
  vect3D_t r_2 = { -1.975719726623028e-11, -9.63716218195963e-13,  0 };
  UINT4 maxorder = 4;

  if ( ( rOrb_n = XLALComputeOrbitalDerivatives ( maxorder, &t0, &edat )) == NULL ) {
    XLALPrintError ("%s: call to XLALComputeOrbitalDerivatives() failed!.\n", fn );
    return TEST_FAIL;
  }

  // ----- first check differences to known first derivatives 0 - 2:
  vect3D_t diff;
  REAL8 reldiff;
  REAL8 reltol = 1e-6;
  UINT4 i;

  /* order = 0 */
  for (i=0; i<3; i++) diff[i] = r_0[i] - rOrb_n->data[0][i];
  reldiff = sqrt ( NORM(diff) / NORM(r_0) );
  if ( reldiff > reltol ) {
    XLALPrintError ("%s: relative error %g on 0th order r_0 exceeds tolerance of %g.\n", fn, reldiff, reltol );
    return TEST_VAL;
  }

  /* order = 1 */
  for (i=0; i<3; i++) diff[i] = r_1[i] - rOrb_n->data[1][i];
  reldiff = sqrt ( NORM(diff) / NORM(r_0) );
  if ( reldiff > reltol ) {
    XLALPrintError ("%s: relative error %g on 1st order r_1 exceeds tolerance of %g.\n", fn, reldiff, reltol );
    return TEST_VAL;
  }

  /* order = 1 */
  for (i=0; i<3; i++) diff[i] = r_2[i] - rOrb_n->data[2][i];
  reldiff = sqrt ( NORM(diff) / NORM(r_0) );
  if ( reldiff > reltol ) {
    XLALPrintError ("%s: relative error %g on 2n order r_2 exceeds tolerance of %g.\n", fn, reldiff, reltol );
    return TEST_VAL;
  }

  // ----- second check derivatives against known analytic results
  REAL8 RorbC = LAL_AU_SI / LAL_C_SI;
  REAL8 Om = LAL_TWOPI / LAL_YRSID_SI;
  vect3D_t nR, nV;
  REAL8 norm;
  vect3D_t rTest_n[maxorder+1];

  norm = NORM( r_0 );
  DIV_VECT (nR, r_0, norm );
  norm = NORM( r_1 );
  DIV_VECT (nV, r_1, norm );

  /* r_0 */
  COPY_VECT(rTest_n[0], nR );
  MULT_VECT(rTest_n[0], RorbC );
  /* r_1 */
  COPY_VECT(rTest_n[1], nV );
  MULT_VECT(rTest_n[1], RorbC*Om );
  /* r_2 */
  COPY_VECT(rTest_n[2], nR );
  MULT_VECT(rTest_n[2], -RorbC*Om*Om );
  /* r_3 */
  COPY_VECT(rTest_n[3], nV );
  MULT_VECT(rTest_n[3], -RorbC*Om*Om*Om );
  /* r_4 */
  COPY_VECT(rTest_n[4], nR );
  MULT_VECT(rTest_n[4], RorbC*Om*Om*Om*Om );


  UINT4 n;
  reltol = 1e-2;	/* 1% error should be enough to enforce order 0-2 are ~1e-5 - 1e-8, order 3-4 are ~ 5e-3*/
  for ( n=0; n <= maxorder; n ++ )
    {
      for (i=0; i<3; i++) diff[i] = rTest_n[n][i] - rOrb_n->data[n][i];
      reldiff = sqrt ( NORM(diff) / NORM(rTest_n[n]) );
      XLALPrintInfo ("order %d: relative difference = %g\n", n, reldiff );
      if ( reldiff > reltol ) {
        XLALPrintError ("%s: relative error %g on r_%d exceeds tolerance of %g.\n", fn, reldiff, n, reltol );
        XLALPrintError ("rd[%d] = {%g, %g, %g},  rTest[%d] = {%g, %g, %g}\n",
                        n, rOrb_n->data[n][0], rOrb_n->data[n][1], rOrb_n->data[n][2],
                        n, rTest_n[n][0], rTest_n[n][1], rTest_n[n][2] );
        return TEST_VAL;
      }
    } /* for n <= maxorder */

  /* free memory */
  XLALDestroyVect3Dlist ( rOrb_n );
  XLALFree ( edat.ephemE );
  XLALFree ( edat.ephemS );

  printf ("OK.\n\n");

  return TEST_OK;

} /* test_XLALComputeOrbitalDerivatives() */
