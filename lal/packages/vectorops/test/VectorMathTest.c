/*
*  Copyright (C) 2015 Reinhard Prix
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
#include <lal/LogPrintf.h>	 // for timing function XLALGetTimeOfDay()
#include <lal/UserInput.h>

#include <lal/VectorMath.h>

// ---------- Macros ----------
#define frand() (rand() / (REAL4)RAND_MAX)
#define Relerr(dx,x) (fabsf(x)>0 ? fabsf((dx)/(x)) : fabsf(dx) )

#define BENCH_FUNCF_1OUT(funcf,abstol0,reltol0)                         \
  {                                                                     \
    XLAL_CHECK ( XLALVectorDeviceSet ( VECTORDEVICE_FPU ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    XLAL_CHECK ( XLALVector##funcf( xOut_Ref, xIn ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    for ( UINT4 dev = VECTORDEVICE_START+1; dev < VECTORDEVICE_END; dev ++ )  \
      {                                                                 \
        if ( !XLALVectorDeviceIsAvailable ( dev ) ) { continue; };      \
        XLAL_CHECK ( XLALVectorDeviceSet ( dev ) == XLAL_SUCCESS, XLAL_EFUNC ); \
        tic = XLALGetTimeOfDay();                                       \
        for (UINT4 l=0; l < Nruns; l ++ ) {                             \
          XLAL_CHECK ( XLALVector##funcf( xOut, xIn ) == XLAL_SUCCESS, XLAL_EFUNC ); \
        }                                                               \
        toc = XLALGetTimeOfDay();                                       \
        maxErr = maxRelerr = 0;                                         \
        for ( UINT4 i = 0; i < xIn->length; i ++ )                      \
          {                                                             \
            REAL4 err = fabsf ( xOut->data[i] - xOut_Ref->data[i] );    \
            REAL4 relerr = Relerr ( err, xOut_Ref->data[i] );           \
            maxErr    = fmaxf ( err, maxErr );                          \
            maxRelerr = fmaxf ( relerr, maxRelerr );                    \
          }                                                             \
        const CHAR *devName = XLALVectorDeviceName(dev);                \
        char buf[256];                                                  \
        sprintf (buf, "%s() %s", #funcf, devName );                     \
        XLALPrintInfo ( "%-16s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                        buf, (REAL8)xIn->length * Nruns / (toc - tic)/1e6, maxErr, (abstol0), maxRelerr, (reltol0) ); \
        XLAL_CHECK ( (maxErr <= (abstol0)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", devName, maxErr, abstol0 ); \
        XLAL_CHECK ( (maxRelerr <= (reltol0)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", devName, maxRelerr, reltol0 ); \
      }                                                                 \
  }

#define BENCH_FUNCF_2OUT(funcf,abstol0,reltol0)                         \
  {                                                                     \
    XLAL_CHECK ( XLALVectorDeviceSet ( VECTORDEVICE_FPU ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    XLAL_CHECK ( XLALVector##funcf( xOut_Ref, xOut2_Ref, xIn ) == XLAL_SUCCESS, XLAL_EFUNC ); \
    for ( UINT4 dev = VECTORDEVICE_START+1; dev < VECTORDEVICE_END; dev ++ )  \
      {                                                                 \
        if ( !XLALVectorDeviceIsAvailable ( dev ) ) { continue; };      \
        XLAL_CHECK ( XLALVectorDeviceSet ( dev ) == XLAL_SUCCESS, XLAL_EFUNC ); \
        tic = XLALGetTimeOfDay();                                       \
        for (UINT4 l=0; l < Nruns; l ++ ) {                             \
          XLAL_CHECK ( XLALVector##funcf( xOut, xOut2, xIn ) == XLAL_SUCCESS, XLAL_EFUNC ); \
        }                                                               \
        toc = XLALGetTimeOfDay();                                       \
        maxErr = maxRelerr = 0;                                         \
        for ( UINT4 i = 0; i < xIn->length; i ++ ) {                    \
          REAL4 err  = fabsf ( xOut->data[i] - xOut_Ref->data[i] );     \
          REAL4 err2 = fabsf ( xOut2->data[i] - xOut2_Ref->data[i] );   \
          REAL4 relerr  = Relerr ( err, xOut_Ref->data[i] );            \
          REAL4 relerr2 = Relerr ( err2, xOut2_Ref->data[i] );          \
          maxErr    = fmaxf ( err, maxErr );                            \
          maxErr    = fmaxf ( err2, maxErr );                           \
          maxRelerr = fmaxf ( relerr, maxRelerr );                      \
          maxRelerr = fmaxf ( relerr2, maxRelerr );                     \
        }                                                               \
        const CHAR *devName = XLALVectorDeviceName(dev);                \
        char buf[256];                                                  \
        sprintf (buf, "%s() %s", #funcf, devName );                     \
        XLALPrintInfo ( "%-16s: %4.0f Mops/sec [maxErr = %7.2g (tol=%7.2g), maxRelerr = %7.2g (tol=%7.2g)]\n", \
                        buf, (REAL8)xIn->length * Nruns / (toc - tic)/1e6, maxErr, (abstol0), maxRelerr, (reltol0) ); \
        XLAL_CHECK ( (maxErr <= (abstol0)), XLAL_ETOL, "%s: absolute error (%g) exceeds tolerance (%g)\n", devName, maxErr, abstol0 ); \
        XLAL_CHECK ( (maxRelerr <= (reltol0)), XLAL_ETOL, "%s: relative error (%g) exceeds tolerance (%g)\n", devName, maxRelerr, reltol0 ); \
      }                                                                 \
  }

// local types
typedef struct
{
  BOOLEAN help;		/**< Print this help/usage message */
  INT4 randSeed;	/**< allow user to specify random-number seed for reproducible noise-realizations */
  INT4 Nruns;		// number of repeated timing 'runs' to average over in order to improve variance of result
} UserInput_t;


// ---------- main ----------
int
main ( int argc, char *argv[] )
{
  UserInput_t XLAL_INIT_DECL(uvar_s);
  UserInput_t *uvar = &uvar_s;

  uvar->randSeed = 1;
  uvar->Nruns = 1;
  // ---------- register user-variable ----------
  XLALregBOOLUserStruct (  help,                'h', UVAR_HELP    , "Print this help/usage message");
  XLALregINTUserStruct  (  randSeed,            's', UVAR_OPTIONAL, "Random-number seed");
  XLALregINTUserStruct  (  Nruns,               'r', UVAR_OPTIONAL, "Number of repeated timing 'runs' to average over (=improves variance)" );

  XLAL_CHECK ( XLALUserVarReadAllInput ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( uvar->help ) {	/* if help was requested, we're done */
    exit (0);
  }

  srand ( uvar->randSeed );
  XLAL_CHECK ( uvar->Nruns >= 1, XLAL_EDOM );
  UINT4 Nruns = (UINT4)uvar->Nruns;

  UINT4 Ntrials = 1000000 + 7;
  REAL4VectorAligned *xIn_a, *xOut_a, *xOut2_a;
  XLAL_CHECK ( ( xIn_a  =  XLALCreateREAL4VectorAligned ( Ntrials, 32 )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOut_a =  XLALCreateREAL4VectorAligned ( Ntrials, 32 )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ( xOut2_a = XLALCreateREAL4VectorAligned ( Ntrials, 32 )) != NULL, XLAL_EFUNC );
  // alias these into standard REAL4Vectors for convenience
  REAL4Vector *xIn   = (REAL4Vector*)xIn_a;
  REAL4Vector *xOut  = (REAL4Vector*)xOut_a;
  REAL4Vector *xOut2 = (REAL4Vector*)xOut2_a;
  REAL4Vector *xOut_Ref, *xOut2_Ref;
  XLAL_CHECK ( (xOut_Ref = XLALCreateREAL4Vector ( Ntrials )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (xOut2_Ref = XLALCreateREAL4Vector ( Ntrials )) != NULL, XLAL_EFUNC );

  REAL8 tic, toc;
  REAL4 maxErr = 0, maxRelerr = 0;
  REAL4 abstol, reltol;

  XLALPrintInfo ("\n%s\n\n", XLALVectorDeviceHelpString() );

  // ---------- input data x in [-1000, 1000] for sin(),cos() ----------
  XLALPrintInfo ("Testing sinf(x), cosf(x) for x in [-1000, 1000]\n");
  for ( UINT4 i = 0; i < xIn->length; i ++ ) {
    xIn->data[i] = 2000 * ( frand() - 0.5 );
  }
  abstol = 2e-7, reltol = 1e-5;
  // ==================== SINF() ====================
  BENCH_FUNCF_1OUT(Sinf,abstol,reltol);
  XLALPrintInfo ("\n");

  // ==================== COSF() ====================
  BENCH_FUNCF_1OUT(Cosf,abstol,reltol);
  XLALPrintInfo ("\n");

  // ==================== SINCOSF() ====================
  BENCH_FUNCF_2OUT(SinCosf,abstol,reltol);
  XLALPrintInfo ("\n");

  // ==================== SINCOSF(2PI*x) ====================
  BENCH_FUNCF_2OUT(SinCosf2PI,abstol,reltol);
  XLALPrintInfo ("\n");

  // ==================== EXPF() ====================
  // ---------- input data x in [-10, 10] for sin(),cos() ----------
  XLALPrintInfo ("Testing expf(x) for x in [-10, 10]\n");
  for ( UINT4 i = 0; i < xIn->length; i ++ ) {
    xIn->data[i] = 20 * ( frand() - 0.5 );
  }
  abstol = 3e-3, reltol = 2e-7;

  BENCH_FUNCF_1OUT(Expf,abstol,reltol);
  XLALPrintInfo ("\n");

  // ==================== LOGF() ====================
  // ---------- input data x in [0, 10000] for logf(x) ----------
  XLALPrintInfo ("Testing logf(x) for x in (0, 10000]\n");
  for ( UINT4 i = 0; i < xIn->length; i ++ ) {
    xIn->data[i] = 10000.0f * frand() + 1e-6;
  } // for i < Ntrials
  abstol = 2e-6, reltol = 2e-7;

  BENCH_FUNCF_1OUT(Logf,abstol,reltol);
  XLALPrintInfo ("\n");

  // ---------- clean up memory ----------
  XLALDestroyREAL4VectorAligned ( xIn_a );
  XLALDestroyREAL4VectorAligned ( xOut_a );
  XLALDestroyREAL4VectorAligned ( xOut2_a );

  XLALDestroyREAL4Vector ( xOut_Ref );
  XLALDestroyREAL4Vector ( xOut2_Ref );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()
