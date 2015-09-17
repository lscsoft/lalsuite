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
 * \brief C code headers for structures for EOBNRv2HM reduced order model (non-spinning version).
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

#ifndef _LALSIMIMREOBNRV2HMROMSTRUCT_H
#define _LALSIMIMREOBNRV2HMROMSTRUCT_H

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


/***************************************************/
/*************** Type definitions ******************/

typedef struct tagSplineList {
    gsl_spline*            spline; /* The gsl spline */
    gsl_interp_accel*      accel; /* The gsl accelerator */
    UINT4                  i; /* Index in the list  */
    struct tagSplineList*  next; /* Pointer to the next list element */
} SplineList;

typedef struct tagEOBNRHMROMdata
{
  gsl_vector* q;
  gsl_vector* freq;
  gsl_matrix* Camp;
  gsl_matrix* Cphi;
  gsl_matrix* Bamp;
  gsl_matrix* Bphi;
  gsl_vector* shifttime;
  gsl_vector* shiftphase;
} EOBNRHMROMdata;

typedef struct tagEOBNRHMROMdata_interp
{
  SplineList* Camp_interp; /* List of splines for the amp coefficients - SplineList, index of reduced basis */
  SplineList* Cphi_interp; /* List of splines for the phase coefficients - SplineList, index of reduced basis */
  SplineList* shifttime_interp; /* interpolated shift in time - SplineList with one element */
  SplineList* shiftphase_interp; /* interpolated shift in phase - SplineList with one element */
} EOBNRHMROMdata_interp;

typedef struct tagEOBNRHMROMdata_coeff
{
  gsl_vector* Camp_coeff;
  gsl_vector* Cphi_coeff;
  double*     shifttime_coeff;
  double*     shiftphase_coeff;
} EOBNRHMROMdata_coeff;

typedef struct tagListmodesEOBNRHMROMdata
{
    EOBNRHMROMdata*                     data; /* The ROM data. */
    UINT4                               l; /* Node mode l  */
    INT4                                m; /* Node submode m  */
    struct tagListmodesEOBNRHMROMdata*  next; /* next pointer */
} ListmodesEOBNRHMROMdata;

typedef struct tagListmodesEOBNRHMROMdata_interp
{
    EOBNRHMROMdata_interp*                     data_interp; /* The splines built from the coefficients. */
    UINT4                                      l; /* Node mode l  */
    INT4                                       m; /* Node submode m  */
    struct tagListmodesEOBNRHMROMdata_interp*  next; /* next pointer */
} ListmodesEOBNRHMROMdata_interp;

typedef struct tagListmodesEOBNRHMROMdata_coeff
{
    EOBNRHMROMdata_coeff*                     data_coeff; /* The data of coefficients. */
    UINT4                                     l; /* Node mode l  */
    INT4                                      m; /* Node submode m  */
    struct tagListmodesEOBNRHMROMdata_coeff*  next; /* next pointer */
} ListmodesEOBNRHMROMdata_coeff;


/**********************************************************/
/**************** Internal functions **********************/

/* Functions associated to list manipulations */
SplineList* SplineList_AddElementNoCopy(
	   SplineList* appended,  /* List structure to prepend to */
	   gsl_spline* spline,  /* spline to contain */
           gsl_interp_accel* accel,  /* accelerator to contain */
	   UINT4 i /* index in the list */
);
SplineList* SplineList_GetElement(
	   SplineList* const splinelist,  /* List structure to get a particular mode from */
	   const UINT4 i /* index in the list */
);
void SplineList_Destroy(
	   SplineList* list  /* List structure to destroy; notice that the content is destroyed too */
);
ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata* indata,  /* data to contain */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_GetMode(
	   ListmodesEOBNRHMROMdata* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
void ListmodesEOBNRHMROMdata_Destroy(
	   ListmodesEOBNRHMROMdata* list  /* List structure to destroy; notice that the data is destroyed too */
);
ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata_interp* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata_interp* data,  /* data to contain */
	   UINT4 l, /* major mode number */
	   INT4 m  /* minor mode number */
);
ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_GetMode(
	   ListmodesEOBNRHMROMdata_interp* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
void ListmodesEOBNRHMROMdata_interp_Destroy(
	   ListmodesEOBNRHMROMdata_interp* list  /* List structure to destroy; notice that the data is destroyed too */
);
ListmodesEOBNRHMROMdata_coeff* ListmodesEOBNRHMROMdata_coeff_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata_coeff* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata_coeff* data,  /* data to contain */
	   UINT4 l, /* major mode number */
	   INT4 m  /* minor mode number */
);
ListmodesEOBNRHMROMdata_coeff* ListmodesEOBNRHMROMdata_coeff_GetMode(
	   ListmodesEOBNRHMROMdata_coeff* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
void ListmodesEOBNRHMROMdata_coeff_Destroy(
	   ListmodesEOBNRHMROMdata_coeff* list  /* List structure to destroy; notice that the data is destroyed too */
);

/* Additional function reproducing XLALSpinWeightedSphericalHarmonic (which disappeared ?) and XLALSimAddMode for frequency-domain structures */
COMPLEX16 SpinWeightedSphericalHarmonic(REAL8 theta, REAL8 phi, INT4 s, INT4 l, INT4 m); /* Currently only supports s=-2, l=2,3,4,5 modes */
INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMREOBNRV2HMROMSTRUCT_H */
