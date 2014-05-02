/*
 *
 * Copyright (C) 2005 Reinhard Prix
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

/*********************************************************************************/
/**
 * \author Reinhard Prix
 * \file
 * \brief Test for ExtrapolatePulsarSpins().
 *
 */
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ExtrapolatePulsarSpins.h>

/** \name Error codes */
/*@{*/
#define PULSARSPINTESTC_ENORM 	0
#define PULSARSPINTESTC_ESUB  	1

#define PULSARSPINTESTC_MSGENORM    "Normal exit"
#define PULSARSPINTESTC_MSGESUB     "Subroutine failed"
/*@}*/


#define RELERROR(x, y) fabs( 2.0 * ((x) - (y)) / ( (x) + (y) ) )


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( PULSARSPINTESTC_ESUB, PULSARSPINTESTC_MSGESUB,      	     \
           "Function call \"" #func "\" failed:" );                  \
    return PULSARSPINTESTC_ESUB;                                     \
  }                                                                  \
} while (0)


/**
 * Very simple test: given spin-params at \f$\tau_0\f$, extrapolate them to
 * \f$\tau_1\f$ and compare to reference-result...
 */
int main(int argc, char *argv[])
{
  LALStatus XLAL_INIT_DECL(status);
  PulsarSpins result;
  PulsarSpins fkdot0, fkdot1;
  REAL8 tolerance = 1.0e-12;
  REAL8 tolerancePhi = 1.0e-6;

  PulsarSpinRange range0, range2;
  PulsarSpinRange rangeResult;

  LIGOTimeGPS epoch0 = {714180733, 0};
  LIGOTimeGPS epoch1 = {714180733 + 94608000, 0};	/* 3 years later */
  LIGOTimeGPS epoch2 = {714180733 - 94608000, 0};	/* 3 years earlier */

  if ( argc == 1 )
    argc = 1;	/* avoid warning */

  /* set up initial spin-values */
  XLAL_INIT_MEM( fkdot0 );
  fkdot0[0] = 300.0;
  fkdot0[1] = -1.e-7;
  fkdot0[2] = 1e-15;
  fkdot0[3] = -1e-22;

  /* Result: (from ExtrPulsarSpins.m) */
  XLAL_INIT_MEM ( result );
  result[0] =  2.809011145986047e+02;
  result[1] = -4.529256832000000e-07;
  result[2] = -8.460800000000001e-15;
  result[3] = -1.000000000000000e-22;

  /* first propagate single spin-vector */
  printf(" \n ----- Test1: LALExtrapolatePulsarSpins() ----- \n");
  printf("Input @ tau0 = %d.%09d : [%.10g, %.10g, %.10g, %.10g ]\n",
	 epoch0.gpsSeconds, epoch0.gpsNanoSeconds,
	 fkdot0[0], fkdot0[1], fkdot0[2], fkdot0[3] );

  SUB ( LALExtrapolatePulsarSpins (&status, fkdot1, epoch1, fkdot0, epoch0), &status );

  printf("Output2@ tau1 = %d.%09d : [%.10g, %.10g, %.10g, %.10g]\n",
	 epoch1.gpsSeconds, epoch1.gpsNanoSeconds,
	 fkdot1[0], fkdot1[1], fkdot1[2], fkdot1[3] );

  printf ("Reference-result:                  : [%.10g, %.10g, %.10g, %.10g]\n",
	  result[0], result[1], result[2], result[3] );

  if ( (RELERROR(fkdot1[0], result[0]) > tolerance) ||
       (RELERROR(fkdot1[1], result[1]) > tolerance) ||
       (RELERROR(fkdot1[2], result[2]) > tolerance) ||
       (RELERROR(fkdot1[3], result[3]) > tolerance) )
    {
      XLALPrintError ( "\nRelative error of XLALExtrapolatePulsarSpins() exceeds tolerance of %g \n\n", tolerance);
      return -1;
    }
  else
    printf ("\n ==> OK. LALExtrapolatePulsarSpins() lies within %g of the reference-result!\n", tolerance);

  /* ----- propagate phase from epoch1 --> epoch0, given fkdot0 ----- */
  {
    REAL8 phi1Result = 3.951107892803490;	/* from ExtrapolatePulsarSpins.m */
    REAL8 phi0, phi1;

    phi0 = 1;
    if ( XLALExtrapolatePulsarPhase ( &phi1, fkdot1, epoch1, phi0, epoch0 ) != XLAL_SUCCESS ) {
      ERROR( PULSARSPINTESTC_ESUB, PULSARSPINTESTC_MSGESUB, "Function call \"XLALExtrapolatePulsarPhase\" failed:" );
      return PULSARSPINTESTC_ESUB;
    }

    printf ("\nExtrapolated phase phi1 = %.16f, Reference-result = %.16f\n", phi1, phi1Result );
    if ( RELERROR(phi1, phi1Result) > tolerancePhi )
      {
	XLALPrintError ( "\nRelative error of LALExtrapolatePulsarPhase() exceeds tolerance of %g \n\n", tolerancePhi);
	return -1;
      }
    else
      printf ("\n ==> OK. LALExtrapolatePulsarPhase() lies within %g of the reference-result!\n", tolerancePhi);
  }

  /* ----- now test LALExtrapolatePulsarSpinRange() ----- */
  /* set up initial spin-range */
  range0.refTime = epoch0;
  XLAL_INIT_MEM ( range0.fkdot );
  range0.fkdot[0] = fkdot0[0];
  range0.fkdot[1] = fkdot0[1];
  range0.fkdot[2] = fkdot0[2];
  range0.fkdot[3] = fkdot0[3];

  XLAL_INIT_MEM ( range0.fkdotBand );
  range0.fkdotBand[0] = 0;
  range0.fkdotBand[1] = -1.0e-7;
  range0.fkdotBand[2] =  1.0e-15;
  range0.fkdotBand[3] = -1.0e-22;

  /* set up the range reference-result (from ExtrapolatePulsarSpins.m) */
  rangeResult.refTime.gpsSeconds = 619572733;
  rangeResult.refTime.gpsNanoSeconds = 0;

  XLAL_INIT_MEM ( rangeResult.fkdot );
  rangeResult.fkdot[0] =  3.280495590653952e+02;
  rangeResult.fkdot[1] = -1.284283366400000e-06;
  rangeResult.fkdot[2] =  1.046080000000000e-14;
  rangeResult.fkdot[3] = -2.000000000000000e-22;

  XLAL_INIT_MEM ( rangeResult.fkdotBand );
  rangeResult.fkdotBand[0] =  2.804955906539521e+01;
  rangeResult.fkdotBand[1] =  6.421416832000000e-07;
  rangeResult.fkdotBand[2] =  1.046080000000000e-14;
  rangeResult.fkdotBand[3] =  1.000000000000000e-22;

  printf(" \n ----- Test2: LALExtrapolatePulsarSpinRange() ----- \n");
  printf ("Input: epoch = %d.%09d \n", range0.refTime.gpsSeconds, range0.refTime.gpsNanoSeconds );
  printf ("Input: fkdot     = [%.10g, %.10g, %.10g, %.10g ]\n",
	  range0.fkdot[0], range0.fkdot[1], range0.fkdot[2],  range0.fkdot[3] );
  printf ("Input: fkdotBand = [%.10g, %.10g, %.10g, %.10g ]\n",
	  range0.fkdotBand[0], range0.fkdotBand[1], range0.fkdotBand[2], range0.fkdotBand[3] );

  SUB ( LALExtrapolatePulsarSpinRange (&status, &range2, epoch2, &range0 ), &status );

  printf ("\n");
  printf ("Output: epoch = %d.%09d \n", range2.refTime.gpsSeconds, range2.refTime.gpsNanoSeconds );
  printf ("Output: fkdot     = [%.10g, %.10g, %.10g, %.10g ]\n",
	  range2.fkdot[0], range2.fkdot[1], range2.fkdot[2], range2.fkdot[3] );
  printf ("Output: fkdotBand = [%.10g, %.10g, %.10g, %.10g ]\n",
	  range2.fkdotBand[0], range2.fkdotBand[1], range2.fkdotBand[2], range2.fkdotBand[3] );

  printf ("\n");
  printf ("Octave : fkdot     = [%.10g, %.10g, %.10g, %.10g ]\n",
	  rangeResult.fkdot[0], rangeResult.fkdot[1], rangeResult.fkdot[2], rangeResult.fkdot[3] );
  printf ("Octave : fkdotBand = [%.10g, %.10g, %.10g, %.10g ]\n",
	  rangeResult.fkdotBand[0], rangeResult.fkdotBand[1], rangeResult.fkdotBand[2], rangeResult.fkdotBand[3] );

  if ( (range2.refTime.gpsSeconds != rangeResult.refTime.gpsSeconds)
       || ( range2.refTime.gpsNanoSeconds != rangeResult.refTime.gpsNanoSeconds ) )
    {
      XLALPrintError ("\nOutput-range has wrong epoch\n");
      return -1;
    }

  if ( (RELERROR(range2.fkdot[0], rangeResult.fkdot[0]) > tolerance) ||
       (RELERROR(range2.fkdot[1], rangeResult.fkdot[1]) > tolerance) ||
       (RELERROR(range2.fkdot[2], rangeResult.fkdot[2]) > tolerance) ||
       (RELERROR(range2.fkdot[3], rangeResult.fkdot[3]) > tolerance) ||

       (RELERROR(range2.fkdotBand[0], rangeResult.fkdotBand[0]) > tolerance) ||
       (RELERROR(range2.fkdotBand[1], rangeResult.fkdotBand[1]) > tolerance) ||
       (RELERROR(range2.fkdotBand[2], rangeResult.fkdotBand[2]) > tolerance) ||
       (RELERROR(range2.fkdotBand[3], rangeResult.fkdotBand[3]) > tolerance)
       )
    {
      XLALPrintError ( "\nRelative error of LALExtrapolatePulsarSpinRange() exceeds tolerance of %g \n", tolerance );
      return -1;
    }
  else
    printf ("\n ==> OK. LALExtrapolatePulsarSpinRange() lies within %g of the reference-result!\n\n", tolerance);

  LALCheckMemoryLeaks();

  return 0;	/* OK */

} /* main() */
