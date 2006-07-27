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
/** \author Reinhard Prix
 * \file 
 * \brief Test for ExtrapolatePulsarSpins().
 *                                                                          
 *********************************************************************************/
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ExtrapolatePulsarSpins.h>

NRCSID (PULSARSPINTESTC, "$Id$");

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
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, PULSARSPINTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              PULSARSPINTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( PULSARSPINTESTC_ESUB, PULSARSPINTESTC_MSGESUB,      	     \
           "Function call \"" #func "\" failed:" );                  \
    return PULSARSPINTESTC_ESUB;                                     \
  }                                                                  \
} while (0)


static LALStatus empty_status;
/** Very simple test: given spin-params at \f$\tau_0\f$, extrapolate them to
 *  \f$\tau_1\f$ and compare to reference-result...
 */
int main(int argc, char *argv[])
{
  LALStatus status = empty_status;
  PulsarSpins result;
  REAL8Vector *fkdot0 = NULL; 
  REAL8Vector *fkdot1 = NULL;
  REAL8 tolerance = 1.0e-12;
  PulsarSpins spins0, spins1;

  LALPulsarSpinRange *range0, *range2;
  LALPulsarSpinRange *rangeResult;

  LIGOTimeGPS epoch0 = {714180733, 0};
  LIGOTimeGPS epoch1 = {714180733 + 94608000, 0};	/* 3 years later */
  LIGOTimeGPS epoch2 = {714180733 - 94608000, 0};	/* 3 years earlier */

  if ( argc == 1 )
    argc = 1;	/* avoid warning */

  fkdot0 = XLALCreateREAL8Vector ( 4 );	/* 3 spindowns */
  fkdot1 = XLALCreateREAL8Vector ( 4 );

  range0 = XLALCreatePulsarSpinRange ( 5 );	/* 4 spindowns */
  range2 = XLALCreatePulsarSpinRange ( 5 );
  rangeResult = XLALCreatePulsarSpinRange ( 5 );
  /* set up initial spin-values */
  fkdot0->data[0] = 300.0;
  fkdot0->data[1] = -1.e-7;
  fkdot0->data[2] = 1e-15;
  fkdot0->data[3] = -1e-22;

  /* do the same extrapolations with the 'PulsarSpins' type */
  spins0[0] = fkdot0->data[0];
  spins0[1] = fkdot0->data[1];
  spins0[2] = fkdot0->data[2];
  spins0[3] = fkdot0->data[3];

  /* Result: (from ExtrPulsarSpins.m) */
  result[0] =  2.809011145986047e+02;
  result[1] = -4.529256832000000e-07;
  result[2] = -8.460800000000001e-15;
  result[3] = -1.000000000000000e-22;

  /* first propagate single spin-vector */
  printf(" \n ----- Test1: LALExtrapolatePulsarSpins() ----- \n");
  printf("Input @ tau0 = %d.%09d : [%.10g, %.10g, %.10g, %.10g ]\n", 
	 epoch0.gpsSeconds, epoch0.gpsNanoSeconds,
	 fkdot0->data[0], fkdot0->data[1], fkdot0->data[2], fkdot0->data[3] );

  SUB ( LALExtrapolatePulsarSpins (&status, fkdot1, epoch1, fkdot0, epoch0), &status );

  printf("Output@ tau1 = %d.%09d : [%.10g, %.10g, %.10g, %.10g]\n",
	 epoch1.gpsSeconds, epoch1.gpsNanoSeconds,
	 fkdot1->data[0], fkdot1->data[1], fkdot1->data[2], fkdot1->data[3] );
  
  /* do the same with LALExtrapolatePulsarSpins2() */
  SUB ( LALExtrapolatePulsarSpins2 (&status, spins1, epoch1, spins0, epoch0), &status );

  printf("Output2@ tau1 = %d.%09d : [%.10g, %.10g, %.10g, %.10g]\n",
	 epoch1.gpsSeconds, epoch1.gpsNanoSeconds,
	 spins1[0], spins1[1], spins1[2], spins1[3] );

  printf ("Reference-result:                  : [%.10g, %.10g, %.10g, %.10g]\n",
	  result[0], result[1], result[2], result[3] );

  if ( (RELERROR(fkdot1->data[0], result[0]) > tolerance) ||
       (RELERROR(fkdot1->data[1], result[1]) > tolerance) ||
       (RELERROR(fkdot1->data[2], result[2]) > tolerance) ||
       (RELERROR(fkdot1->data[3], result[3]) > tolerance) )
    {
      LALPrintError ( "\nRelative error of XLALExtrapolatePulsarSpins() exceeds tolerance of %g \n\n",
		      tolerance);
      return -1;
    }
  else
    printf ("\n ==> OK. LALExtrapolatePulsarSpins() lies within %g of the reference-result!\n",
	    tolerance);

  if ( (RELERROR(spins1[0], result[0]) > tolerance) ||
       (RELERROR(spins1[1], result[1]) > tolerance) ||
       (RELERROR(spins1[2], result[2]) > tolerance) ||
       (RELERROR(spins1[3], result[3]) > tolerance) )
    {
      LALPrintError ( "\nRelative error of LALExtrapolatePulsarSpins2() exceeds tolerance of %g \n\n",
		      tolerance);
      return -1;
    }
  else
    printf ("\n ==> OK. LALExtrapolatePulsarSpins2() lies within %g of the reference-result!\n",
	    tolerance);

  /* ----- propagate phase from epoch1 --> epoch0, given fkdot0 ----- */
  {

    REAL8 phi1Result = 3.951107892803490;	/* from ExtrapolatePulsarSpins.m */
    REAL8 phi0, phi1;

    phi0 = 1;
    SUB ( LALExtrapolatePulsarPhase ( &status, &phi1, spins1, epoch1, phi0, epoch0 ), &status );

    printf ("\nExtrapolated phase phi1 = %.16f, Reference-result = %.16f\n", phi1, phi1Result );
    if ( RELERROR(phi1, phi1Result) > tolerance ) 
      {
	LALPrintError ( "\nRelative error of LALExtrapolatePulsarPhase() exceeds tolerance of %g \n\n",
			tolerance);
	return -1;
      }
    else
      printf ("\n ==> OK. LALExtrapolatePulsarPhase() lies within %g of the reference-result!\n",
	      tolerance);
  }

  /* ----- now test LALExtrapolatePulsarSpinRange() ----- */
  /* set up initial spin-range */
  range0->epoch = epoch0;
  range0->fkdot->data[0] = fkdot0->data[0];
  range0->fkdot->data[1] = fkdot0->data[1];
  range0->fkdot->data[2] = fkdot0->data[2];
  range0->fkdot->data[3] = fkdot0->data[3];
  range0->fkdot->data[4] = 2.0000e-29;

  range0->fkdotBand->data[0] = 0;
  range0->fkdotBand->data[1] = -1.0e-7;
  range0->fkdotBand->data[2] =  1.0e-15;
  range0->fkdotBand->data[3] = -1.0e-22;
  range0->fkdotBand->data[4] = -1.0e-30;

  /* set up the range reference-result (from ExtrapolatePulsarSpins.m) */
  rangeResult->epoch.gpsSeconds = 619572733;
  rangeResult->epoch.gpsNanoSeconds = 0;
  rangeResult->fkdot->data[0] = 3.914735849716052e+02;
  rangeResult->fkdot->data[1] = -4.106967813079040e-06;
  rangeResult->fkdot->data[2] =  9.549219980800000e-14;
  rangeResult->fkdot->data[3] =  -2.092160000000000e-21;
  rangeResult->fkdot->data[4] =   1.900000000000000e-29;

  rangeResult->fkdotBand->data[0] =  3.138766569203784e+01;
  rangeResult->fkdotBand->data[1] =  7.832759055339521e-07;
  rangeResult->fkdotBand->data[2] =  1.493613683200000e-14;
  rangeResult->fkdotBand->data[3] =  1.946080000000001e-22;
  rangeResult->fkdotBand->data[4] =  1.000000000000000e-30;


  printf(" \n ----- Test2: LALExtrapolatePulsarSpinRange() ----- \n");
  printf ("Input: epoch = %d.%09d \n", range0->epoch.gpsSeconds, range0->epoch.gpsNanoSeconds );
  printf ("Input: fkdot     = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range0->fkdot->data[0], range0->fkdot->data[1], range0->fkdot->data[2], 
		 range0->fkdot->data[3], range0->fkdot->data[4] ); 
  printf ("Input: fkdotBand = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range0->fkdotBand->data[0], range0->fkdotBand->data[1], range0->fkdotBand->data[2], 
		 range0->fkdotBand->data[3], range0->fkdotBand->data[4] ); 

  SUB ( LALExtrapolatePulsarSpinRange (&status, range2, epoch2, range0 ), &status ); 

  printf ("\n");
  printf ("Output: epoch = %d.%09d \n", range2->epoch.gpsSeconds, range2->epoch.gpsNanoSeconds );
  printf ("Output: fkdot     = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range2->fkdot->data[0], range2->fkdot->data[1], range2->fkdot->data[2], 
		 range2->fkdot->data[3], range2->fkdot->data[4] ); 
  printf ("Output: fkdotBand = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range2->fkdotBand->data[0], range2->fkdotBand->data[1], range2->fkdotBand->data[2], 
		 range2->fkdotBand->data[3], range2->fkdotBand->data[4] ); 

  printf ("\n");
  printf ("Octave : fkdot     = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 rangeResult->fkdot->data[0], rangeResult->fkdot->data[1], 
		 rangeResult->fkdot->data[2], rangeResult->fkdot->data[3], 
		 rangeResult->fkdot->data[4] ); 
  printf ("Octave : fkdotBand = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 rangeResult->fkdotBand->data[0], rangeResult->fkdotBand->data[1], 
		 rangeResult->fkdotBand->data[2], rangeResult->fkdotBand->data[3], 
		 rangeResult->fkdotBand->data[4] ); 
  
  if ( (range2->epoch.gpsSeconds != rangeResult->epoch.gpsSeconds) 
       || ( range2->epoch.gpsNanoSeconds != rangeResult->epoch.gpsNanoSeconds ) )
    {
      LALPrintError ("\nOutput-range has wrong epoch\n");
      return -1;
    }

  if ( (RELERROR(range2->fkdot->data[0], rangeResult->fkdot->data[0]) > tolerance) ||
       (RELERROR(range2->fkdot->data[1], rangeResult->fkdot->data[1]) > tolerance) ||
       (RELERROR(range2->fkdot->data[2], rangeResult->fkdot->data[2]) > tolerance) ||
       (RELERROR(range2->fkdot->data[3], rangeResult->fkdot->data[3]) > tolerance) ||
       (RELERROR(range2->fkdot->data[4], rangeResult->fkdot->data[4]) > tolerance) ||

       (RELERROR(range2->fkdotBand->data[0], rangeResult->fkdotBand->data[0]) > tolerance) ||
       (RELERROR(range2->fkdotBand->data[1], rangeResult->fkdotBand->data[1]) > tolerance) ||
       (RELERROR(range2->fkdotBand->data[2], rangeResult->fkdotBand->data[2]) > tolerance) ||
       (RELERROR(range2->fkdotBand->data[3], rangeResult->fkdotBand->data[3]) > tolerance) ||
       (RELERROR(range2->fkdotBand->data[4], rangeResult->fkdotBand->data[4]) > tolerance) 
       )
    {
      LALPrintError ( "\nRelative error of LALExtrapolatePulsarSpinRange() exceeds tolerance of %g \n",
		      tolerance );
      return -1;
    }
  else
    printf ("\n ==> OK. LALExtrapolatePulsarSpinRange() lies within %g of the reference-result!\n\n",
	    tolerance);

  LALDDestroyVector (&status, &fkdot0);
  LALDDestroyVector (&status, &fkdot1);

  XLALDestroyPulsarSpinRange(range0);
  XLALDestroyPulsarSpinRange(range2);
  XLALDestroyPulsarSpinRange(rangeResult);

  LALCheckMemoryLeaks();

  return 0;	/* OK */

} /* main() */
