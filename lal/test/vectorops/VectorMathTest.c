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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/
#include <math.h>
#include <stdlib.h>

#include <config.h>

#include <lal/LALVCSInfo.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>	 // for timing function XLALGetCPUTime()
#include <lal/UserInput.h>

#include <lal/VectorMath.h>

/* for access to internal prototypes for generic (GEN) functions, for reference results */
#include <vectorops/VectorMath_internal.h>

// ---------- Macros ----------
#define frand() (rand() / (REAL4)RAND_MAX)
#define Relerr(dx,x) (fabsf(x)>0 ? fabsf((dx)/(x)) : fabsf(dx) )
#define Relerrd(dx,x) (fabs(x)>0 ? fabs((dx)/(x)) : fabs(dx) )
#define cRelerr(dx,x) (cabsf(x)>0 ? cabsf((dx)/(x)) : fabsf(dx) )

// ----- test and benchmark operators with 1 REAL4 vector input and 1 INT4 vector output (S2I) ----------
#define TESTBENCH_VECTORMATH_S2I(name,in)                               \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( xOutRefI4->data, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( xOutI4->data, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      XLAL_CHECK ( xOutI4->data[i] == xOutRefI4->data[i], XLAL_ETOL, "%s: found element #%u (%i) differs from reference (%i)", #name, i, xOutI4->data[i], xOutRefI4->data[i] ); \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec\n", XLALVector##name##REAL4_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6 ); \
  }

// ----- test and benchmark operators with 1 REAL4 vector input and 1 REAL4 scalar output (S2s) ----------
#define TESTBENCH_VECTORMATH_S2s(name,in)                               \
  {                                                                     \
    REAL4 sOutRef, sOut;                                                \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( &sOutRef, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                           \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( &sOut, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                           \
    maxErr = fabsf( sOut - sOutRef );                                   \
    maxRelerr = Relerr( maxErr, sOutRef );                                \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL4_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL4", maxRelerr, reltol ); \
  }

// ----- test and benchmark operators with 1 REAL4 vector input and 1 REAL4 vector output (S2S) ----------
#define TESTBENCH_VECTORMATH_S2S(name,in)                               \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( xOutRef, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
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

// ----- test and benchmark operators with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define TESTBENCH_VECTORMATH_S2SS(name,in)                              \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( xOutRef, xOutRef2, xIn, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
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

// ----- test and benchmark operators with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define TESTBENCH_VECTORMATH_SS2S(name,in1,in2)                         \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( xOutRef, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
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

// ----- test and benchmark operators with 3 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define TESTBENCH_VECTORMATH_SSS2S(name,in1,in2,in3)                    \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( xOutRef, in1, in2, in3, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( xOut, in1, in2, in3, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
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

// ----- test and benchmark operators with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
#define TESTBENCH_VECTORMATH_SS2uU(name,in1,in2)                        \
  {                                                                     \
    UINT4 xCount = 0, xCountRef = 0;                                    \
    XLAL_CHECK ( XLALVector##name##REAL4_GEN( &xCountRef, xOutRefU4->data, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL4( &xCount, xOutU4->data, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
      XLAL_CHECK ( xCount == xCountRef, XLAL_ETOL, "%s: count of found elements (%u) differs from reference (%u)", #name, xCount, xCountRef ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    for ( UINT4 i = 0; i < xCount; i ++ )                               \
    {                                                                   \
      XLAL_CHECK ( xOutU4->data[i] == xOutRefU4->data[i], XLAL_ETOL, "%s: found element #%u (%u) differs from reference (%u)", #name, i, xOutU4->data[i], xOutRefU4->data[i] ); \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec\n", XLALVector##name##REAL4_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6 ); \
  }

// ----- test and benchmark operators with 1 REAL8 vector input and 1 REAL8 scalar output (D2d) ----------
#define TESTBENCH_VECTORMATH_D2d(name,in)                               \
  {                                                                     \
    REAL8 sOutRef, sOut;                                                \
    XLAL_CHECK ( XLALVector##name##REAL8_GEN( &sOutRef, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                           \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL8( &sOut, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                           \
    maxErr = fabs( sOut - sOutRef );                                   \
    maxRelerr = Relerrd( maxErr, sOutRef );                                \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL8_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxRelerr, reltol ); \
  }

// ----- test and benchmark operators with 2 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
#define TESTBENCH_VECTORMATH_DD2D(name,in1,in2)                         \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL8_GEN( xOutRefD, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL8( xOutD, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL8 err = fabs ( xOutD[i] - xOutRefD[i] );                      \
      REAL8 relerr = Relerrd ( err, xOutRefD[i] );                       \
      maxErr    = fmax ( err, maxErr );                                \
      maxRelerr = fmax ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL8_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxRelerr, reltol ); \
  }

// ----- test and benchmark operators with 3 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
#define TESTBENCH_VECTORMATH_DDD2D(name,in1,in2,in3)                    \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL8_GEN( xOutRefD, in1, in2, in3, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL8( xOutD, in1, in2, in3, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL8 err = fabs ( xOutD[i] - xOutRefD[i] );                      \
      REAL8 relerr = Relerrd ( err, xOutRefD[i] );                       \
      maxErr    = fmax ( err, maxErr );                                \
      maxRelerr = fmax ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL8_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxRelerr, reltol ); \
  }

// ----- test and benchmark operators with 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
#define TESTBENCH_VECTORMATH_CC2C(name,in1,in2)                         \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##COMPLEX8_GEN( xOutRefC, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##COMPLEX8( xOutC, in1, in2, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL4 err = cabsf ( xOutC[i] - xOutRefC[i] );                      \
      REAL4 relerr = cRelerr ( err, xOutRefC[i] );                       \
      maxErr    = fmaxf ( err, maxErr );                                \
      maxRelerr = fmaxf ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##COMPLEX8_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "COMPLEX8", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "COMPLEX8", maxRelerr, reltol ); \
  }

// ----- test and benchmark operators with 3 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
#define TESTBENCH_VECTORMATH_CCC2C(name,in1,in2,in3)                    \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##COMPLEX8_GEN( xOutRefC, in1, in2, in3, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                             \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##COMPLEX8( xOutC, in1, in2, in3, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                             \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL4 err = cabsf ( xOutC[i] - xOutRefC[i] );                      \
      REAL4 relerr = cRelerr ( err, xOutRefC[i] );                       \
      maxErr    = fmaxf ( err, maxErr );                                \
      maxRelerr = fmaxf ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##COMPLEX8_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "COMPLEX8", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "COMPLEX8", maxRelerr, reltol ); \
  }

// ----- test and benchmark operators with 1 REAL8 vector input and 1 REAL8 vector output (D2D) ----------
#define TESTBENCH_VECTORMATH_D2D(name,in)                               \
  {                                                                     \
    XLAL_CHECK ( XLALVector##name##REAL8_GEN( xOutRefD, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    tic = XLALGetCPUTime();                                           \
    for (UINT4 l=0; l < Nruns; l ++ ) {                                 \
      XLAL_CHECK ( XLALVector##name##REAL8( xOutD, in, Ntrials ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    }                                                                   \
    toc = XLALGetCPUTime();                                           \
    maxErr = maxRelerr = 0;                                             \
    for ( UINT4 i = 0; i < Ntrials; i ++ )                              \
    {                                                                   \
      REAL8 err = fabs ( xOut[i] - xOutRef[i] );                      \
      REAL8 relerr = Relerrd ( err, xOutRef[i] );                       \
      maxErr    = fmax ( err, maxErr );                                \
      maxRelerr = fmax ( relerr, maxRelerr );                          \
    }                                                                   \
    XLALPrintInfo ( "%-32s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                    XLALVector##name##REAL8_name, (REAL8)Ntrials * Nruns / (toc - tic)/1e6, maxErr, (abstol), maxRelerr, (reltol) ); \
    XLAL_CHECK ( (maxErr <= (abstol)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxErr, abstol ); \
    XLAL_CHECK ( (maxRelerr <= (reltol)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", #name "REAL8", maxRelerr, reltol ); \
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
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  srand ( uvar->randSeed );
  XLAL_CHECK ( uvar->Nruns >= 1, XLAL_EDOM );
  UINT4 Nruns = (UINT4)uvar->Nruns;

  UINT4 Ntrials = 1000000 + 7;

  INT4Vector *xOutI4;
  INT4Vector *xOutRefI4;
  XLAL_CHECK ( ( xOutI4 = XLALCreateINT4Vector ( Ntrials )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOutRefI4 = XLALCreateINT4Vector ( Ntrials )) != NULL, XLAL_EFUNC );

  UINT4Vector *xOutU4;
  UINT4Vector *xOutRefU4;
  XLAL_CHECK ( ( xOutU4 = XLALCreateUINT4Vector ( Ntrials )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOutRefU4 = XLALCreateUINT4Vector ( Ntrials )) != NULL, XLAL_EFUNC );

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

  REAL8VectorAligned *xInD_a, *xIn2D_a, *xOutD_a, *xOutRefD_a;
  XLAL_CHECK ( ( xInD_a   = XLALCreateREAL8VectorAligned ( Ntrials, uvar->inAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xIn2D_a  = XLALCreateREAL8VectorAligned ( Ntrials, uvar->inAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOutD_a  = XLALCreateREAL8VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (xOutRefD_a= XLALCreateREAL8VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );

  // extract aligned REAL8 vectors from these
  REAL8 *xInD      = xInD_a->data;
  REAL8 *xIn2D     = xIn2D_a->data;
  REAL8 *xOutD     = xOutD_a->data;
  REAL8 *xOutRefD  = xOutRefD_a->data;

  COMPLEX8VectorAligned *xInC_a, *xIn2C_a, *xOutC_a, *xOutRefC_a;
  XLAL_CHECK ( ( xInC_a   = XLALCreateCOMPLEX8VectorAligned ( Ntrials, uvar->inAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xIn2C_a  = XLALCreateCOMPLEX8VectorAligned ( Ntrials, uvar->inAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOutC_a  = XLALCreateCOMPLEX8VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (xOutRefC_a  = XLALCreateCOMPLEX8VectorAligned ( Ntrials, uvar->outAlign )) != NULL, XLAL_EFUNC );

  // extract aligned COMPLEX8 vectors from these
  COMPLEX8 *xInC      = xInC_a->data;
  COMPLEX8 *xIn2C     = xIn2C_a->data;
  COMPLEX8 *xOutC     = xOutC_a->data;
  COMPLEX8 *xOutRefC  = xOutRefC_a->data;

  REAL8 tic, toc;
  REAL4 maxErr = 0, maxRelerr = 0;
  REAL4 abstol, reltol;

  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i] = 2000 * ( frand() - 0.5 );
    xInD[i] = xIn[i];
  }
  abstol = 2e-7, reltol = 1e-5;

  XLALPrintInfo ("Testing (INT4) for x in [-1000, 1000]\n");
  // ==================== (INT4) ====================
  TESTBENCH_VECTORMATH_S2I(INT4From,xIn);

  // ==================== SCALAR MAX ====================
  XLALPrintInfo ("Testing scalar max for x in [-1000, 1000]\n");
  TESTBENCH_VECTORMATH_S2s(ScalarMax,xIn);
  TESTBENCH_VECTORMATH_D2d(ScalarMax,xInD);

  XLALPrintInfo ("Testing sin(x), cos(x) for x in [-1000, 1000]\n");
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

  abstol = 4e-3, reltol = 3e-7;
  TESTBENCH_VECTORMATH_S2S(Exp,xIn);

  // ==================== LOG() ====================
  XLALPrintInfo ("\nTesting log(x) for x in (0, 10000]\n");
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i] = 10000.0f * frand() + 1e-6;
  } // for i < Ntrials
  abstol = 2e-6, reltol = 4e-7;

  TESTBENCH_VECTORMATH_S2S(Log,xIn);

  // ==================== ADD,MUL,ROUND ====================
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i]  = -10000.0f + 20000.0f * frand() + 1e-6;
    xIn2[i] = -10000.0f + 20000.0f * frand() + 1e-6;
    xInD[i] = -100000.0 + 200000.0 * frand() + 1e-6;
    xIn2D[i]= -100000.0 + 200000.0 * frand() + 1e-6;
    xInC[i] = -10000.0f + 20000.0f * frand() + 1e-6 + ( -10000.0f + 20000.0f * frand() + 1e-6 ) * _Complex_I;
    xIn2C[i]= -10000.0f + 20000.0f * frand() + 1e-6 + ( -10000.0f + 20000.0f * frand() + 1e-6 ) * _Complex_I;
  } // for i < Ntrials
  abstol = 2e-7, reltol = 2e-7;

  XLALPrintInfo ("\nTesting round(x) for x in (-10000, 10000]\n");
  TESTBENCH_VECTORMATH_S2S(Round,xIn);
  TESTBENCH_VECTORMATH_D2D(Round,xInD);

  XLALPrintInfo ("\nTesting add,multiply,shift,scale(x,y) for x,y in (-10000, 10000]\n");
  TESTBENCH_VECTORMATH_SS2S(Add,xIn,xIn2);
  TESTBENCH_VECTORMATH_SS2S(Sub,xIn,xIn2);
  TESTBENCH_VECTORMATH_SS2S(Multiply,xIn,xIn2);
  TESTBENCH_VECTORMATH_SS2S(Max,xIn,xIn2);

  TESTBENCH_VECTORMATH_SS2S(Shift,xIn[0],xIn2);
  TESTBENCH_VECTORMATH_SS2S(Scale,xIn[0],xIn2);
  TESTBENCH_VECTORMATH_SSS2S(ScaleAdd,xIn[0],xIn,xIn2);

  TESTBENCH_VECTORMATH_DD2D(Add,xInD,xIn2D);
  TESTBENCH_VECTORMATH_DD2D(Sub,xInD,xIn2D);
  TESTBENCH_VECTORMATH_DD2D(Multiply,xInD,xIn2D);
  TESTBENCH_VECTORMATH_DD2D(Max,xInD,xIn2D);

  TESTBENCH_VECTORMATH_DD2D(Shift,xInD[0],xIn2D);
  TESTBENCH_VECTORMATH_DD2D(Scale,xInD[0],xIn2D);
  TESTBENCH_VECTORMATH_DDD2D(ScaleAdd,xInD[0],xInD,xIn2D);

  TESTBENCH_VECTORMATH_CC2C(Multiply,xInC,xIn2C);
  TESTBENCH_VECTORMATH_CC2C(Add,xInC,xIn2C);

  TESTBENCH_VECTORMATH_CC2C(Scale,xInC[0],xIn2C);
  TESTBENCH_VECTORMATH_CC2C(Shift,xInC[0],xIn2C);
  TESTBENCH_VECTORMATH_CCC2C(ScaleAdd,xIn[0],xInC,xIn2C);

  // ==================== FIND ====================
  for ( UINT4 i = 0; i < Ntrials; i ++ ) {
    xIn[i]  = -10000.0f + 20000.0f * frand() + 1e-6;
    xIn2[i] = -10000.0f + 20000.0f * frand() + 1e-6;
  } // for i < Ntrials

  XLALPrintInfo ("\nTesting find for x,y in (-10000, 10000]\n");
  TESTBENCH_VECTORMATH_SS2uU(FindVectorLessEqual,xIn,xIn2);

  TESTBENCH_VECTORMATH_SS2uU(FindScalarLessEqual,xIn[0],xIn2);

  XLALPrintInfo ("\n");

  // ---------- clean up memory ----------
  XLALDestroyINT4Vector ( xOutI4 );
  XLALDestroyINT4Vector ( xOutRefI4 );

  XLALDestroyUINT4Vector ( xOutU4 );
  XLALDestroyUINT4Vector ( xOutRefU4 );

  XLALDestroyREAL4VectorAligned ( xIn_a );
  XLALDestroyREAL4VectorAligned ( xIn2_a );
  XLALDestroyREAL4VectorAligned ( xOut_a );
  XLALDestroyREAL4VectorAligned ( xOut2_a );

  XLALDestroyREAL4VectorAligned ( xOutRef_a );
  XLALDestroyREAL4VectorAligned ( xOutRef2_a );

  XLALDestroyREAL8VectorAligned ( xInD_a );
  XLALDestroyREAL8VectorAligned ( xIn2D_a );
  XLALDestroyREAL8VectorAligned ( xOutD_a );
  XLALDestroyREAL8VectorAligned ( xOutRefD_a );

  XLALDestroyCOMPLEX8VectorAligned ( xInC_a );
  XLALDestroyCOMPLEX8VectorAligned ( xIn2C_a );
  XLALDestroyCOMPLEX8VectorAligned ( xOutC_a );
  XLALDestroyCOMPLEX8VectorAligned ( xOutRefC_a );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()
