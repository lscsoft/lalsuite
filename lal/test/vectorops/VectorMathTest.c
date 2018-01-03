/*
*  Copyright (C) 2015 Reinhard Prix, Karl Wette
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
#include <stdlib.h>

#include <config.h>

#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>	 // for timing function XLALGetCPUTime()
#include <lal/UserInput.h>

#include <lal/VectorMath.h>

/* for access to internal prototypes for FPU functions, for reference results */
#include "../../src/vectorops/VectorMath_internal.h"

// ---------- Macros ----------
#define frand() (rand() / (REAL4)RAND_MAX)
#define Relerr(dx,x) (fabsf(x)>0 ? fabsf((dx)/(x)) : fabsf(dx) )

// ----- test and benchmark operators with 1 REAL4 vector input and 1 REAL4 vector output (S2S) ----------
#define TESTBENCH_VECTORMATH_S2S(name,in)                               \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_FPU( xOutRef, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                           \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( xOut, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                           \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL4 err = fabsf ( xOut[i] - xOutRef[i] );                      \
      REAL4 relerr = Relerr ( err, xOutRef[i] );                       \
      maxErr    = fmaxf ( err, maxErr );                                \
      maxRelerr = fmaxf ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL4_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxRelerr, reltol ); \
  }

#define TESTBENCH_VECTORMATH_S2SS(name,in)                              \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_FPU( xOutRef, xOutRef2, xIn, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                               \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( xOut, xOut2, xIn, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                           \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ ) {                            \
      REAL4 err1 = fabsf ( xOut[i] - xOutRef[i] );                     \
      REAL4 err2 = fabsf ( xOut2[i] - xOutRef2[i] );                   \
      REAL4 relerr1 = Relerr ( err1, xOutRef[i] );                      \
      REAL4 relerr2 = Relerr ( err2, xOutRef2[i] );                    \
      maxErr    = fmaxf ( err1, maxErr );                                \
      maxErr    = fmaxf ( err2, maxErr );                               \
      maxRelerr = fmaxf ( relerr1, maxRelerr );                          \
      maxRelerr = fmaxf ( relerr2, maxRelerr );                         \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL4_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, reltol ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxRelerr, reltol ); \
  }


#define TESTBENCH_VECTORMATH_SS2S(name,in1,in2)                         \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_FPU( xOutRef, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( xOut, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL4 err = fabsf ( xOut[i] - xOutRef[i] );                      \
      REAL4 relerr = Relerr ( err, xOutRef[i] );                       \
      maxErr    = fmaxf ( err, maxErr );                                \
      maxRelerr = fmaxf ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL4_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxRelerr, reltol ); \
  }


// local types
typedef struct
{
  INT4 randSeed;	/**< allow user to specify random-number seed for reproducible noise-realizations */
  INT4 Nruns;		// number of repated timing 'runs' to average over in order to improve variance of result
  INT4 inAlign;		// alignment of input vectors; default is sizeof(void*), i.e. no particular alignment
  INT4 outAlign;	// alignment of output vectors; default is sizeof(void*), i.e. no particular alignment
} UserInput_t;


// ---------- main ----------
int
main ( int argc, char *argv[] )
{
  UserInput_t XLAL_INIT_DECL(uvar_s);
  UserInput_t *uvar = &uvar_s;

  uvar->randSeed = 1;
  uvar->Nruns = 1;
  uvar->inAlign = uvar->outAlign = sizeof(void*);
  // ---------- register user-variable ----------
  XLALRegisterUvarMember(  randSeed,            INT4, 's', OPTIONAL, "Random-number seed");
  XLALRegisterUvarMember(  Nruns,               INT4, 'r', OPTIONAL, "Number of repeated timing 'runs' to average over (=improves variance)" );
  XLALRegisterUvarMember(  inAlign,             INT4, 'a', OPTIONAL, "Alignment of input vectors; default is sizeof(void*), i.e. no particular alignment" );
  XLALRegisterUvarMember(  outAlign,            INT4, 'b', OPTIONAL, "Alignment of output vectors; default is sizeof(void*), i.e. no particular alignment" );

  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  srand ( uvar->randSeed );
  XLAL_CHECK ( uvar->Nruns >= 1, XLAL_EDOM );
  UINT4 Nruns = (UINT4)uvar->Nruns;

  UINT4 Ntrials = 1000000 + 7;
  REAL4VectorAligned *xIn_a, *xIn2_a, *xOut_a, *xOut2_a;
  XLAL_CHECK ( ( xIn_a   = XLALCreateREAL4VectorAligned ( Ntrials, uvar->inAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xIn2_a  = XLALCreateREAL4VectorAligned ( Ntrials, uvar->inAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOut_a  = XLALCreateREAL4VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOut2_a = XLALCreateREAL4VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );
  REAL4VectorAligned *xOutRef_a, *xOutRef2_a;
  XLAL_CHECK ( (xOutRef_a  = XLALCreateREAL4VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (xOutRef2_a = XLALCreateREAL4VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );

  // extract aligned REAL4 vectors from these
  REAL4 *xIn      = xIn_a->data;
  REAL4 *xIn2     = xIn2_a->data;
  REAL4 *xOut     = xOut_a->data;
  REAL4 *xOut2    = xOut2_a->data;
  REAL4 *xOutRef  = xOutRef_a->data;
  REAL4 *xOutRef2 = xOutRef2_a->data;

  REAL8 tic, toc;
  REAL4 maxErr = 0, maxRelerr = 0;
  REAL4 abstol, reltol;

  XLALPrintInfo ("Testing sin(x), cos(x) for x in [-1000, 1000]\n");
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i] = 2000 * ( frand() - 0.5 );
  }
  abstol = 2e-7, reltol = 1e-5;
  // ==================== SIN() ====================
  TESTBENCH_VECTORMATH_S2S(Sin,xIn);

  // ==================== COS() ====================
  TESTBENCH_VECTORMATH_S2S(Cos,xIn);

  // ==================== SINCOS() ====================
  TESTBENCH_VECTORMATH_S2SS(SinCos,xIn);

  // ==================== SINCOS(2PI*x) ====================
  TESTBENCH_VECTORMATH_S2SS(SinCos2Pi,xIn);

  // ==================== EXP() ====================
  XLALPrintInfo ("\nTesting exp(x) for x in [-10, 10]\n");
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i] = 20 * ( frand() - 0.5 );
  }

  abstol = 3e-3, reltol = 2e-7;
  TESTBENCH_VECTORMATH_S2S(Exp,xIn);

  // ==================== LOG() ====================
  XLALPrintInfo ("\nTesting log(x) for x in (0, 10000]\n");
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i] = 10000.0f * frand() + 1e-6;
  } // for i < Ntrials
  abstol = 2e-6, reltol = 2e-7;

  TESTBENCH_VECTORMATH_S2S(Log,xIn);

  // ==================== ADD,MUL ====================
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i]  = -10000.0f + 20000.0f * frand() + 1e-6;
    xIn2[i] = -10000.0f + 20000.0f * frand() + 1e-6;
  } // for i < Ntrials
  abstol = 2e-7, reltol = 2e-7;

  XLALPrintInfo ("\nTesting add,multiply,shift,scale(x,y) for x,y in (-10000, 10000]\n");
  TESTBENCH_VECTORMATH_SS2S(Add,xIn,xIn2);
  TESTBENCH_VECTORMATH_SS2S(Multiply,xIn,xIn2);

  TESTBENCH_VECTORMATH_SS2S(Shift,xIn[0],xIn2);
  TESTBENCH_VECTORMATH_SS2S(Scale,xIn[0],xIn2);

  XLALPrintInfo ("\n");

  // ---------- clean up memory ----------
  XLALDestroyREAL4VectorAligned ( xIn_a );
  XLALDestroyREAL4VectorAligned ( xIn2_a );
  XLALDestroyREAL4VectorAligned ( xOut_a );
  XLALDestroyREAL4VectorAligned ( xOut2_a );

  XLALDestroyREAL4VectorAligned ( xOutRef_a );
  XLALDestroyREAL4VectorAligned ( xOutRef2_a );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()
