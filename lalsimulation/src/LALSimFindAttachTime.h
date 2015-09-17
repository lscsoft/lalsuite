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
//#include <lal/LALSimInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>


#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

//#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimIMREOBHybridRingdown.c"


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
    REAL8 *radiusVec,
    int *found,
    REAL8* tMaxOmega
);


double XLALSimLocateAmplTime(
    REAL8Vector *timeHi, 
    COMPLEX16Vector *hP22,
    REAL8 *radiusVec,
    int *found,
    REAL8* tMaxAmp
);
   
INT4 XLALSimCheckRDattachment(
    REAL8Vector * signal1,	/**<< Real of inspiral waveform to which we attach ringdown */
    REAL8Vector * signal2,	/**<< Imag of inspiral waveform to which we attach ringdown */
    REAL8* ratio,           /**<< output ratio  */
    const REAL8 tAtt,       /**<< time of RD attachment */
    const INT4 l,	/**<< Current mode l */
    const INT4 m,	/**<< Current mode m */
    const REAL8 dt,	/**<< Sample time step (in seconds) */
    const REAL8 mass1,	/**<< First component mass (in Solar masses) */
    const REAL8 mass2,	/**<< Second component mass (in Solar masses) */
    const REAL8 spin1x,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1y,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1z,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin2x,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2y,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2z,	/**<<The spin of the second object; only needed for spin waveforms */
    REAL8Vector * timeVec,	/**<< Vector containing the time values */
    REAL8Vector * matchrange,	/**<< Time values chosen as points for performing comb matching */
    Approximant approximant,	/**<<The waveform approximant being used */
    const REAL8 JLN           /**<< cosine of the angle between J and LN at the light ring */   
);

int XLALSimAdjustRDattachmentTime( 
    REAL8Vector * signal1,	/**<< Output Real of inspiral waveform to which we attach ringdown */
    REAL8Vector * signal2,	/**<< Output Imag of inspiral waveform to which we attach ringdown */
    COMPLEX16TimeSeries* h22,   /**<< input time series (inspiral) */
    COMPLEX16TimeSeries* h2m2,  /**<< input time series (inspiral) */
    REAL8* ratio22,      /**<< output ratio for 2,2 mode */
    REAL8* ratio2m2,     /**<< output ratio  for 2,-2 mode*/
    REAL8* tAtt,       /**<< output/input time of RD attachment */
    const REAL8 thr,        /**<< threshold on the ratio */
    const REAL8 dt,	/**<< Sample time step (in seconds) */
    const REAL8 mass1,	/**<< First component mass (in Solar masses) */
    const REAL8 mass2,	/**<< Second component mass (in Solar masses) */
    const REAL8 spin1x,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1y,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1z,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin2x,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2y,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2z,	/**<<The spin of the second object; only needed for spin waveforms */
    REAL8Vector * timeVec,	/**<< Vector containing the time values */
    REAL8Vector * matchrange,	/**<< Time values chosen as points for performing comb matching */
    Approximant approximant,	/**<<The waveform approximant being used */
    const REAL8 JLN,           /**<< cosine of the angle between J and LN at the light ring */
    const REAL8 combsize,        /**<< combsize for RD attachment */
    const REAL8 tMaxOmega,
    const REAL8 tMaxAmp
);
 

    
#endif    
