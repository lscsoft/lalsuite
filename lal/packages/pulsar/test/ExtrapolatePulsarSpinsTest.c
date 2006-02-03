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
  REAL8 result[]= {347.6632471314574, 2.36975876347904E-6, 8.104593663999999E-14, 1.79216E-21, 2.0E-29};
  REAL8Vector *fkdot0 = NULL; 
  REAL8Vector *fkdot1 = NULL;
  REAL8 tolerance = 1.0e-12;
  LALPulsarSpinRange *range0, *range1;
  LALPulsarSpinRange *rangeResult;

  LIGOTimeGPS epoch0 = {714180733, 0};
  LIGOTimeGPS epoch1 = {714180733 + 94608000, 0};	/* 3 years later */

  if ( argc == 1 )
    argc = 1;	/* avoid warning */

  fkdot0 = XLALCreateREAL8Vector ( 5 );	/* 4 spindowns */
  fkdot1 = XLALCreateREAL8Vector ( 5 );
  range0 = XLALCreatePulsarSpinRange ( 5 );
  range1 = XLALCreatePulsarSpinRange ( 5 );
  rangeResult = XLALCreatePulsarSpinRange ( 5 );
  /* set up initial spin-values */
  fkdot0->data[0] = 300.0;
  fkdot0->data[1] = -1.e-7;
  fkdot0->data[2] = 1e-15;
  fkdot0->data[3] = -1e-22;
  fkdot0->data[4] = 2e-29;

  /* set up initial spin-range */
  range0->epoch = epoch0;
  range0->fkdot->data[0] = fkdot0->data[0];
  range0->fkdot->data[1] = fkdot0->data[1];
  range0->fkdot->data[2] = fkdot0->data[2];
  range0->fkdot->data[3] = fkdot0->data[3];
  range0->fkdot->data[4] = fkdot0->data[4];

  range0->fkdotBand->data[0] = 0;
  range0->fkdotBand->data[1] = -1.0e-7;
  range0->fkdotBand->data[2] =  1.0e-15;
  range0->fkdotBand->data[3] = -1.0e-22;
  range0->fkdotBand->data[4] = -1.0e-30;

  /* set up the range reference-result (from ExtrapolatePulsarSpins.m) */
  rangeResult->epoch.gpsSeconds = 808788733;
  rangeResult->epoch.gpsNanoSeconds = 0;
  rangeResult->fkdot->data[0] = 3.207509182714196e+02;
  rangeResult->fkdot->data[1] = 1.681090857945088e-06;
  rangeResult->fkdot->data[2] = 6.710979980799999e-14;
  rangeResult->fkdot->data[3] = 1.597552000000000e-21;
  rangeResult->fkdot->data[4] = 1.900000000000000e-29;

  rangeResult->fkdotBand->data[0] = 3.138766569203784e+01;
  rangeResult->fkdotBand->data[1] = 7.832759055339521e-07;
  rangeResult->fkdotBand->data[2] = 1.493613683200000e-14;
  rangeResult->fkdotBand->data[3] = 1.946080000000003e-22;
  rangeResult->fkdotBand->data[4] = 1.000000000000000e-30;

  /* first propagate single spin-vector */
  SUB ( LALExtrapolatePulsarSpins (&status, fkdot1, epoch1, fkdot0, epoch0), &status );

  printf("Input spin-params: (%.10g, %.10g, %.10g, %.10g, %.10g)\n", 
	 fkdot0->data[0], fkdot0->data[1], fkdot0->data[2], fkdot0->data[3], fkdot0->data[4]);
  printf("Output spin-params: (%.10g, %.10g, %.10g, %.10g, %.10g)\n",
	 fkdot1->data[0], fkdot1->data[1], fkdot1->data[2], fkdot1->data[3], fkdot1->data[4]);
  
  /* this is a maxima-script to calculate the result:
     f(t) := f0 + f1 * t/1! + f2 * t^2 / 2! + f3 * t^3/3! + f4 * t^4/4!;
     f1dot(t) := at ( diff( f(x), x, 1), x = t);
     f2dot(t) := at ( diff( f(x), x, 2), x = t);
     f3dot(t) := at ( diff( f(x), x, 3), x = t);
     f4dot(t) := at ( diff( f(x), x, 4), x = t);
     f0:300.0; f1:-1.e-7; f2:1e-15; f3:-1e-22; f4:2e-29$
     tau:94608000;
     [ f(tau), f1dot(tau), f2dot(tau), f3dot(tau), f4dot(tau)]; 

     to exectute, simply paste into a maxima-buffer. The result is:
     [347.6632471314574, 2.36975876347904E-6, 8.104593663999999E-14, 1.79216E-21, 2.0E-29]
     and we use this to check the output of ExtrapolatePulsarSpins():
  */
  printf ("Reference-result is: (%.10g, %.10g, %.10g, %.10g, %.10g)\n",
	  result[0], result[1], result[2], result[3], result[4]);

  if ( (RELERROR(fkdot1->data[0], result[0]) > tolerance) ||
       (RELERROR(fkdot1->data[1], result[1]) > tolerance) ||
       (RELERROR(fkdot1->data[2], result[2]) > tolerance) ||
       (RELERROR(fkdot1->data[3], result[3]) > tolerance) ||
       (RELERROR(fkdot1->data[4], result[4]) > tolerance) )
    {
      LALPrintError ( "\nRelative error of XLALExtrapolatePulsarSpins() exceeds tolerance of %g ",
		      tolerance);
      return -1;
    }
  else
    printf ("\nOK. LALExtrapolatePulsarSpins() lies within %g of the reference-result!\n",
	    tolerance);
  
  LALDDestroyVector (&status, &fkdot0);
  LALDDestroyVector (&status, &fkdot1);


  /* ----- now test LALExtrapolatePulsarSpinRange() ----- */
  LALPrintError ("Input: epoch = %d.%09d \n", range0->epoch.gpsSeconds, range0->epoch.gpsNanoSeconds );
  LALPrintError ("Input: fkdot     = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range0->fkdot->data[0], range0->fkdot->data[1], range0->fkdot->data[2], 
		 range0->fkdot->data[3], range0->fkdot->data[4] ); 
  LALPrintError ("Input: fkdotBand = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range0->fkdotBand->data[0], range0->fkdotBand->data[1], range0->fkdotBand->data[2], 
		 range0->fkdotBand->data[3], range0->fkdotBand->data[4] ); 

  SUB ( LALExtrapolatePulsarSpinRange (&status, range1, epoch1, range0 ), &status ); 

  LALPrintError ("\n");
  LALPrintError ("Output: epoch = %d.%09d \n", range1->epoch.gpsSeconds, range1->epoch.gpsNanoSeconds );
  LALPrintError ("Output: fkdot     = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range1->fkdot->data[0], range1->fkdot->data[1], range1->fkdot->data[2], 
		 range1->fkdot->data[3], range1->fkdot->data[4] ); 
  LALPrintError ("Output: fkdotBand = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 range1->fkdotBand->data[0], range1->fkdotBand->data[1], range1->fkdotBand->data[2], 
		 range1->fkdotBand->data[3], range1->fkdotBand->data[4] ); 

  LALPrintError ("\n");
  LALPrintError ("Octave : fkdot     = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 rangeResult->fkdot->data[0], rangeResult->fkdot->data[1], 
		 rangeResult->fkdot->data[2], rangeResult->fkdot->data[3], 
		 rangeResult->fkdot->data[4] ); 
  LALPrintError ("Octave : fkdotBand = [%.10g, %.10g, %.10g, %.10g, %.10g ]\n",
		 rangeResult->fkdotBand->data[0], rangeResult->fkdotBand->data[1], 
		 rangeResult->fkdotBand->data[2], rangeResult->fkdotBand->data[3], 
		 rangeResult->fkdotBand->data[4] ); 
  
  if ( (range1->epoch.gpsSeconds != rangeResult->epoch.gpsSeconds) 
       || ( range1->epoch.gpsNanoSeconds != rangeResult->epoch.gpsNanoSeconds ) )
    {
      LALPrintError ("\nOutput-range has wrong epoch\n");
      return -1;
    }

  if ( (RELERROR(range1->fkdot->data[0], rangeResult->fkdot->data[0]) > tolerance) ||
       (RELERROR(range1->fkdot->data[1], rangeResult->fkdot->data[1]) > tolerance) ||
       (RELERROR(range1->fkdot->data[2], rangeResult->fkdot->data[2]) > tolerance) ||
       (RELERROR(range1->fkdot->data[3], rangeResult->fkdot->data[3]) > tolerance) ||
       (RELERROR(range1->fkdot->data[4], rangeResult->fkdot->data[4]) > tolerance) ||

       (RELERROR(range1->fkdotBand->data[0], rangeResult->fkdotBand->data[0]) > tolerance) ||
       (RELERROR(range1->fkdotBand->data[1], rangeResult->fkdotBand->data[1]) > tolerance) ||
       (RELERROR(range1->fkdotBand->data[2], rangeResult->fkdotBand->data[2]) > tolerance) ||
       (RELERROR(range1->fkdotBand->data[3], rangeResult->fkdotBand->data[3]) > tolerance) ||
       (RELERROR(range1->fkdotBand->data[4], rangeResult->fkdotBand->data[4]) > tolerance) 
       )
    {
      LALPrintError ( "\nRelative error of LALExtrapolatePulsarSpinRange() exceeds tolerance of %g \n",
		      tolerance );
      return -1;
    }
  else
    printf ("\nOK. LALExtrapolatePulsarSpinRange() lies within %g of the reference-result!\n",
	    tolerance);


  XLALDestroyPulsarSpinRange(range0);
  XLALDestroyPulsarSpinRange(range1);
  XLALDestroyPulsarSpinRange(rangeResult);

  LALCheckMemoryLeaks();

  return 0;	/* OK */

} /* main() */
