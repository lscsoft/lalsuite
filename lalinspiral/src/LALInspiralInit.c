/*
*  Copyright (C) 2007 Thomas Cokelaer, Drew Keppel
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
 * \author Cokelaer T.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Module to initialize some parameters for waveform generation.
 *
 * \heading{Description}
 * The input parameters is an InspiralTemplate structure which provides the waveform parameters
 * such as masses, lower frequency ... . The function \c LALInspiralInit calls the
 * \c LALInspiralParameterCalc function in order to  compute all the mass parameters. Then,
 * \c LALInspiralRestrictedAmplitude function is called to get the restricted newtonian
 * amplitude. LALInspiralWavelength, LALInspiralSetup and LALInspiralChooseModel are also called
 * in order to estimate the waveform length which is stored in an output structure called
 * \c InspiralInit. We also stored Energy, flux and evolution function of flux and energy in
 * that structure.
 *
 * The  \c LALInspiralChooseModel function might failed or send a non zero status code.
 * That function force it to be zero therefore the codes which  use LALInspiralInit (mainly
 * injection code right now) won't stopped. Of course, if status code is non zero, we have to keep
 * trace of it. Thus, the length of the waveform is fixed to zero in case of problems such as
 * negative length, cutoff frequency lower than the lower cutoff frequency ... .
 *
 * \heading{Uses}
 * \code
 * LALInspiralParameterCalc
 * LALInspiralRestrictedAmplitude
 * LALInspiralWaveLength
 * LALInspiralChooseModel
 * LALInspiralSetup
 * \endcode
 *
 * \heading{Notes}
 * There is only one assert on the InspiralTemplate variable since  all relevant asserts
 * are already included in the different functions which are called throughout the LALInspiralInit
 * function.
 *
 */


#include <lal/LALInspiral.h>
#define  LALINSPIRALINIT_LENGTHOVERESTIMATION  0.1       /* 10 % */

void
LALInspiralInit (LALStatus        *status,
		 InspiralTemplate *params,
		 InspiralInit     *paramsInit)
{
  XLALPrintDeprecationWarning("LALInspiralInit", "XLALInspiralInit");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if ( XLALInspiralInit(params, paramsInit) == XLAL_FAILURE )
  {
    ABORTXLAL( status );
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int
XLALInspiralInit (InspiralTemplate *params,
		  InspiralInit     *paramsInit)
{
  UINT4 ndx;
  REAL8 x;


  if (params == NULL)
    XLAL_ERROR(XLAL_EFAULT);

  if ( XLALInspiralParameterCalc(params) == XLAL_FAILURE )
  {
    XLAL_ERROR(XLAL_EFUNC);
  }

  if ( XLALInspiralRestrictedAmplitude(params) == XLAL_FAILURE )
  {
    XLAL_ERROR(XLAL_EFUNC);
  }

  if ( XLALInspiralSetup(&(paramsInit->ak), params) == XLAL_FAILURE )
  {
    XLAL_ERROR(XLAL_EFUNC);
  }
  if ( XLALInspiralChooseModel(&(paramsInit->func), &(paramsInit->ak), params) 
       == XLAL_FAILURE )
  {
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* The parameters have been initialized now. However, we can have some problems
     with the LALInspiralChooseModel related to bad estimation of the length.

     We first need to check that the length is not wrong.
     Then, to check that flso is > fLow

     keep a security length higher than the one given by ChooseModel
  */

  if( params->fCutoff < params->fLower){
    XLALPrintWarning(LALINSPIRALH_MSGEFLOWER);
    paramsInit->nbins = 0;

    XLALPrintInfo("XLAL Info - %s: Estimated Length (seconds) requested = %f | fCutoff = %f\n", __func__, paramsInit->ak.tn, params->fCutoff);

    return XLAL_SUCCESS;
  }

  if( paramsInit->ak.tn <=0 || params->tC <= 0){
    XLALPrintWarning(LALINSPIRALH_MSGESIZE);
    paramsInit->nbins = 0;
    XLALPrintInfo("XLAL Info - %s: Estimated Length (seconds) requested = %f\n", __func__, paramsInit->ak.tn);

    return XLAL_SUCCESS;
  }

  /* if everything is fine is ChooseModel then we
   * estimate the waveform length. */
  /* we add a minimal value and 10 % of overestimation */

  x = (1.+ LALINSPIRALINIT_LENGTHOVERESTIMATION)
    * (paramsInit->ak.tn + 1 ) * params->tSampling
    + params->nStartPad + params->nEndPad ;
  ndx = ceil(log10(x)/log10(2.));
  paramsInit->nbins = pow(2, ndx) ;

  XLALPrintInfo("XLAL Info - %s: Estimated Length (seconds)= %f | Allocated length (bins) = %d\n", __func__, paramsInit->ak.tn, paramsInit->nbins);


  return XLAL_SUCCESS;
}
