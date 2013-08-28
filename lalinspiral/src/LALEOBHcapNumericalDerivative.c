/*
*  Copyright (C) 2010 Craig Robinson 
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
 * \author Craig Robinson
 *
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do just that.
 *
 */

#include <lal/LALInspiral.h>

#include <gsl/gsl_deriv.h>


/* We need to encapsulate the data for the GSL derivative function */
struct HcapDerivParams
{
   REAL8Vector           *values;
   REAL8Vector           *tmpVec;
   InspiralDerivativesIn *ak;
   UINT4                 varyParam;
};

static double GSLHamiltonianWrapper( double x, void *params )

int XLALHcapNumericalDerivative(
                          REAL8Vector           *values,
                          REAL8Vector           *dvalues,
                          InspiralDerivativesIn *ak
                               )
{

  static const char func[] = "XLALHcapNumericalDerivative";

  static const REAL8 STEP_SIZE = 1.0e-4;

  HcapDerivParams params;

  gsl_function F;
  INT4         gslStatus;

  UINT4 i;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */ 
  params.values  = values;
  params.ak      = ak;
  params.tmpVec  = XLALCreateREAL8Vector( values->length );
  if ( !params->tmpVec )
  {
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  F.function = &GSLHamiltonianWrapper;
  F.params   = params;

  /* Now calculate derivatives w.r.t. each parameter */
  for ( i = 0; i < values->length; i++ )
  {
    params->varyParam = i;
    XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values->data[i], 
                      STEP_SIZE, dvalues->data[i], absErr ) );

    if ( gslStatus != GSL_SUCCESS )
    {
      XLALDestroyREAL8Vector( params.tmpVec );
      XLALPrintError( "Failure in GSL function\n" );
      XLAL_ERROR( func, XLAL_EFUNC );
    }
  }

  return XLAL_SUCCESS;
}


/* Wrapper for GSL to call the Hamiltonian function */
static double GSLHamiltonianWrapper( double x, void *params )
{
  struct HcapDerivParams *dParams = (struct HcapDerivParams *)params;

  REAL8Vector *tmpVec = dParams->tmpVec;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec->data, dParams->values->data, 
              tmpVec->length * sizeof(REAL8) );

  /* Set the relevant entry in the vector to the correct value */
  tmpVec->data[varyParam] = x;

  return XLALEffectiveHamiltonian( tmpVec, dParams->ak );
}
