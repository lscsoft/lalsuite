/*
 *  LALInferenceCalibrationErrors.c:  Bayesian Followup Calibration Error routines.
 *
 *  Copyright (C) 2014 Salvatore Vitale, Ryan Lynch
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
 
#ifndef _LALINFERENCECALIBRATIONERRORS_H  /* Double-include protection. */
#define _LALINFERENCECALIBRATIONERRORS_H
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/TimeSeries.h>
#include <lal/LALInference.h>
#include <lal/LALDatatypes.h>
#include <math.h>

// GSL PACKAGES
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define LALCALIBRATIONERRORSH_ENULL 1
#define LALCALIBRATIONERRORSH_EDIV0 2
#define LALCALIBRATIONERRORSH_MSGENULL "Null pointer"
#define LALCALIBRATIONERRORSH_MSGEDIV0 "Division by zero"

void LALInferenceApplyCalibrationErrors(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine);
#endif

