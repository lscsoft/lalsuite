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

#define RELERROR(x, y) fabs( 2.0 * ((x) - (y)) / ( (x) + (y) ) )

static LALStatus empty_status;
/** Very simple test: given spin-params at \f$\tau_0\f$, extrapolate them to
 *  \f$\tau_1\f$ and compare to reference-result...
 */
int main(void)
{
  LALStatus status = empty_status;
  REAL8 result[] = {347.6632471314574, 2.36975876347904E-6, 8.104593663999999E-14, 1.79216E-21, 2.0E-29};
  REAL8Vector *fkdot0 = NULL; 
  REAL8Vector *fkdot1 = NULL;
  REAL8 tolerance = 1.0e-12;

  LIGOTimeGPS epoch0 = {714180733, 0};
  LIGOTimeGPS epoch1 = {714180733 + 94608000, 0};	/* 3 years later */

  LALDCreateVector (&status, &fkdot0, 5 );	/* 4 spindowns */
  LALDCreateVector (&status, &fkdot1, 5 );	/* 4 spindowns */

  fkdot0->data[0] = 300.0;
  fkdot0->data[1] = -1.e-7;
  fkdot0->data[2] = 1e-15;
  fkdot0->data[3] = -1e-22;
  fkdot0->data[4] = 2e-29;

  if ( XLALExtrapolatePulsarSpins (fkdot1, epoch1, fkdot0, epoch0) )
    {
      int code = xlalErrno;
      LALPrintError ("\nXLALExtrapolatePulsarSpins() failed (xlalErrno = %d)!\n\n", code);
      return -1;
    }

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
    printf ("\nOK. XLALExtrapolatePulsarSpins() lies within %g of the reference-result!\n",
	    tolerance);
  
  LALDDestroyVector (&status, &fkdot0);
  LALDDestroyVector (&status, &fkdot1);

  return 0;	/* OK */

} /* main() */
