#ifndef _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_H
#define _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_H

#include <unistd.h>
#include <math.h>
#include <gsl/gsl_deriv.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficients.c"
#include "LALSimIMREOBFactorizedWaveform.c"

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


static
REAL8 GSLSpinHamiltonianWrapper( double x, void *params );

static
     int XLALSpinHcapNumericalDerivative(
                 double UNUSED     t,         /**<< UNUSED */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               );

static
     REAL8 XLALSpinHcapNumDerivWRTParam(
                 const INT4 paramIdx,      /**<< Index of the parameters */
                 const REAL8 values[],     /**<< Dynamical variables */
                 SpinEOBParams *funcParams /**<< EOB Parameters */
                 );

static
     int XLALSpinHcapNumericalDerivativeNoFlux(
                 double UNUSED     t,         /**<< UNUSED */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               );




#endif 
