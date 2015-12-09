/*
 * Copyright (C) 2015 S. Babak
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

#ifndef _LALSIMFINDATTACHTIME_H
#define _LALSIMFINDATTACHTIME_H

#include <math.h>

#include <lal/LALGSL.h>

#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>


#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

//#include "LALSimIMRSpinEOB.h"


#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>




#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

double  XLALSimLocateOmegaTime(
    REAL8Array *dynamicsHi,
    unsigned int numdynvars,
    unsigned int retLenHi,
    SpinEOBParams   seobParams,
    SpinEOBHCoeffs  seobCoeffs,
    REAL8 m1,
    REAL8 m2,
    REAL8Vector *radiusVec,
    int *found,
    REAL8* tMaxOmega
);


double XLALSimLocateAmplTime(
    REAL8Vector *timeHi,
    COMPLEX16Vector *hP22,
    REAL8Vector *radiusVec,
    int *found,
    REAL8* tMaxAmp
);

INT4 XLALSimCheckRDattachment(
    REAL8Vector * signal1,
    REAL8Vector * signal2,
    REAL8* ratio,
    const REAL8 tAtt,
    const INT4 l,
    const INT4 m,
    const REAL8 dt,
    const REAL8 mass1,
    const REAL8 mass2,
    const REAL8 spin1x,
    const REAL8 spin1y,
    const REAL8 spin1z,
    const REAL8 spin2x,
    const REAL8 spin2y,
    const REAL8 spin2z,
    REAL8Vector * timeVec,
    REAL8Vector * matchrange,
    Approximant approximant,
    const REAL8 JLN,
    REAL8 * timediff
);

int XLALSimAdjustRDattachmentTime(
    REAL8Vector * signal1,
    REAL8Vector * signal2,
    COMPLEX16TimeSeries* h22,
    COMPLEX16TimeSeries* h2m2,
    REAL8* ratio22,
    REAL8* ratio2m2,
    REAL8* tAtt,
    const REAL8 thr,
    const REAL8 dt,
    const REAL8 mass1,
    const REAL8 mass2,
    const REAL8 spin1x,
    const REAL8 spin1y,
    const REAL8 spin1z,
    const REAL8 spin2x,
    const REAL8 spin2y,
    const REAL8 spin2z,
    REAL8Vector * timeVec,
    REAL8Vector * matchrange,
    Approximant approximant,
    const REAL8 JLN,
    const REAL8 combsize,
    const REAL8 tMaxOmega,
    const REAL8 tMaxAmp
);



#endif
