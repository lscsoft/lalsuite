/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash
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

#include <lal/Interpolate.h>
#include <lal/LALInspiralBank.h>


/**
 * \ingroup LALInspiralBank_h
 * \brief NO BRIEF DESCRIPTION
 * \author Brown, D. A., and Sathyaprakash, B. S.
 *
 * The moments of the noise curve are defined as
 * \f{equation}{
 * I(q)  \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}}
 * \frac{x^{-q/3}}{S_{h}(x)} \, dx \,.
 * \f}
 * This function calculates the integrand of this integral, i.e.\ for a given \f$x\f$
 * it calculates
 * \f{equation}{
 * \frac{x^{-q/3}}{S_{h}(x)} \,\,.
 * \f}
 * by interpolating the frequency series containing \f$S_h(f)\f$.
 */
void
LALInspiralMomentsIntegrand(
    LALStatus  *status,		/**< LAL status pointer */
    REAL8      *integrand,	/**< [out] the value of the integrand */
    REAL8       x,		/**< [in] the point where the integrand is required */
    void       *params		/**< [in] of type \c InspiralMomentsIn containing the details required in moments calculation */
    )

{
   InspiralMomentsIn   *integrandParams;

   DInterpolateOut      interpOutput;
   DInterpolatePar      interpParams;

   UINT4                numInterpPts = 4;
   REAL8                f[4];
   REAL8                fMin;
   REAL8                fMax;
   REAL8                deltaF;
   UINT8                freqIndex;

   INITSTATUS(status);
   ATTATCHSTATUSPTR( status );

   ASSERT( params, status,
       LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );

   integrandParams = (InspiralMomentsIn *) params;

   /* check that we have a pointer to a frequency series and it has data */
   ASSERT( integrandParams->shf, status,
       LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
   ASSERT( integrandParams->shf->data, status,
       LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
   ASSERT( integrandParams->shf->data->data, status,
       LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );

   /* the minimum and maximum frequency where we have four points */
   deltaF = integrandParams->shf->deltaF;
   fMin = integrandParams->shf->f0 + deltaF;
   fMax = integrandParams->shf->f0 +
     ((REAL8) integrandParams->shf->data->length - 2 ) * deltaF;

   /* locate the nearest point in the frequency series to the desired point */
   if ( x <= fMin )
   {
     freqIndex = 1;
   }
   else if ( x >= fMax )
   {
     freqIndex = integrandParams->shf->data->length - 3;
   }
   else
   {
     freqIndex = (UINT8) floor( (x - integrandParams->shf->f0) / deltaF );
   }

   /* set up the frequency values for interpolation */
   f[0] = (REAL8)(freqIndex - 1) * deltaF;
   f[1] = (REAL8)(freqIndex) * deltaF;
   f[2] = (REAL8)(freqIndex + 1) * deltaF;
   f[3] = (REAL8)(freqIndex + 2) * deltaF;

   /* set up the interpolation parameters */
   interpParams.n = numInterpPts;
   interpParams.x = f;
   interpParams.y = integrandParams->shf->data->data + freqIndex - 1;

   /* perform the interpolation... */
   LALDPolynomialInterpolation( status->statusPtr, &interpOutput, x,
       &interpParams );
   CHECKSTATUSPTR( status );

   /* ...and check for a negative value of shf */
   if ( interpOutput.y < 0 )
   {
     ABORT( status, 999, "interpolation output is negative" );
   }

   /* now compute the integrand using shf */
   *integrand = pow( x, -(integrandParams->ndx) ) / interpOutput.y;

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
