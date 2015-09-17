/*
 *  Copyright (C) 2014 Sylvain Marsat
 *  Reduced Order Model for EOBNRv2HM
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
 * \author Sylvain Marsat
 *
 * \file
 *
 * \brief C code headers for EOBNRv2HM reduced order model (non-spinning version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details on the reduced order method.
 * See arXiv:1106.1021 for the EOBNRv2HM model.
 *
 * Borrows from the SEOBNR ROM LAL code written by Michael Puerrer and John Veitch.
 *
 * The binary data files are available at [TBD].
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * Parameter ranges:
 *   q = 1-6
 *   No spin
 *   Mtot >= 20Msun for fstart=9Hz
 *
 */

#ifndef _LALSIMIMREOBNRV2HMROM_H
#define _LALSIMIMREOBNRV2HMROM_H

#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>

#include "LALSimIMREOBNRv2HMROMstruct.h"

/**************************************************/
/**************** Prototypes **********************/

/* Functions to load, initalize and cleanup data */
INT4 EOBNRv2HMROM_Init_LALDATA(void);
INT4 EOBNRv2HMROM_Init(const char dir[]);

void EOBNRHMROMdata_Init(EOBNRHMROMdata **data);
void EOBNRHMROMdata_interp_Init(EOBNRHMROMdata_interp **data_interp);
void EOBNRHMROMdata_coeff_Init(EOBNRHMROMdata_coeff **data_coeff);

void EOBNRHMROMdata_Cleanup(EOBNRHMROMdata *data);
void EOBNRHMROMdata_interp_Cleanup(EOBNRHMROMdata_interp *data_interp);
void EOBNRHMROMdata_coeff_Cleanup(EOBNRHMROMdata_coeff *data_coeff);

/* Functions to read data */
void Err_Handler(const char *reason, const char *file, int line, int gsl_errno);
INT4 Read_Vector(const char dir[], const char fname[], gsl_vector *v);
INT4 Read_Matrix(const char dir[], const char fname[], gsl_matrix *m);
INT4 Read_Data_Mode(const char dir[], const INT4 mode[2], EOBNRHMROMdata *data);

/* Functions to interpolate the data and to evaluate the interpolated data for a given q */

INT4 Evaluate_Spline_Data(
  const REAL8 q,                           /* Input: q-value for which projection coefficients should be evaluated */
  const EOBNRHMROMdata_interp* data_interp,  /* Input: data in interpolated form */
  EOBNRHMROMdata_coeff* data_coeff           /* Output: vectors of projection coefficients and shifts in time and phase */
);

INT4 Interpolate_Spline_Data(
  const EOBNRHMROMdata *data,           /* Input: data in vector/matrix form to interpolate */
  EOBNRHMROMdata_interp *data_interp    /* Output: interpolated data */
);

/* Functions for waveform reconstruction */

INT4 EOBNRv2HMROMCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  REAL8 phiRef,
  REAL8 deltaF,
  REAL8 fLow,
  REAL8 fHigh,
  REAL8 fRef,
  REAL8 distance,
  REAL8 inclination,
  REAL8 Mtot_sec,
  REAL8 q);

/* Compute waveform in LAL format */
/*INT4 XLALSimIMREOBNRv2HMROM(
  struct tagCOMPLEX16FrequencySeries **hptilde,
  struct tagCOMPLEX16FrequencySeries **hctilde,
  REAL8 phiRef,
  REAL8 deltaF,
  REAL8 fLow,
  REAL8 fHigh,
  REAL8 fRef,
  REAL8 distance,
  REAL8 inclination,
  REAL8 m1SI,
  REAL8 m2SI,
  const int higherModesFlag);*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMREOBNRV2HMROM_H */
