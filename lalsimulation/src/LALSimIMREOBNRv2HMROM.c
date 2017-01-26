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
 * \brief C code for EOBNRv2HM reduced order model (non-spinning version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details on the reduced order method.
 * See arXiv:1106.1021 for the EOBNRv2HM model.
 *
 * Borrows from the SEOBNR ROM LAL code written by Michael Puerrer and John Veitch.
 *
 * The binary data files are available at [TBD].
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * Parameter ranges:
 *   q = 1-11.98
 *   No spin
 *   Mtot >= 8Msun for fstart=10Hz
 *
 */

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
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>

#include "LALSimIMRSEOBNRROMUtilities.c"
#include "LALSimIMREOBNRv2HMROMUtilities.c"



/*****************************************************************************/
/**************************** General definitions ****************************/

/* NOTE: some integer constants #define statements are in LALSimIMREOBNRv2HMROMUtilities.c */

/* Number and list of modes to generate - to be modified ultimately to allow for a selection of the desired mode(s) */
/* By convention the first mode of the list, assumed to be the 22, is used to set phiRef */
#define EOBNRV2_ROM_NUM_MODES_MAX 5
static INT4 nbmode = EOBNRV2_ROM_NUM_MODES_MAX;
static const INT4 listmode[EOBNRV2_ROM_NUM_MODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {4,4}, {5,5} };

/* Maximal mass ratio covered by the model - q=12 (almost) for now */
static const REAL8 q_max = 11.9894197212;
/* Minimal geometric frequency covered by the model - f=8Hz for M=10Msol for now */
static const REAL8 Mf_ROM_min = 0.0003940393857519091;
/* Maximal geometric frequency covered by the model - Mf=0.285 for the 55 mode - used as default */
static const REAL8 Mf_ROM_max = 0.285;
/* Total mass (in units of solar mass) used to generate the ROM - used to restore the correct amplitude (global factor M) when coming back to physical units */
static const REAL8 M_ROM = 10.;


/******************************************************************/
/********* Definitions for the persistent structures **************/

/* SEOBNR-ROM structures are generalized to lists */
static ListmodesEOBNRHMROMdata* __lalsim_EOBNRv2HMROM_data_init = NULL; /* for initialization only */
static ListmodesEOBNRHMROMdata** const __lalsim_EOBNRv2HMROM_data = &__lalsim_EOBNRv2HMROM_data_init;
static ListmodesEOBNRHMROMdata_interp* __lalsim_EOBNRv2HMROM_interp_init = NULL; /* for initialization only */
static ListmodesEOBNRHMROMdata_interp** const __lalsim_EOBNRv2HMROM_interp = &__lalsim_EOBNRv2HMROM_interp_init;
/* Tag indicating whether the data has been loaded and interpolated */
static INT4 __lalsim_EOBNRv2HMROM_setup = XLAL_FAILURE; /* To be set to XLAL_SUCCESS after initialization*/


/****************************************************************************/
/******************************** Prototypes ********************************/

/* Functions to interpolate the data and to evaluate the interpolated data for a given q */
static INT4 Evaluate_Spline_Data(
  const REAL8 q,                           /* Input: q-value for which projection coefficients should be evaluated */
  const EOBNRHMROMdata_interp* data_interp,  /* Input: data in interpolated form */
  EOBNRHMROMdata_coeff* data_coeff           /* Output: vectors of projection coefficients and shifts in time and phase */
);
static INT4 Interpolate_Spline_Data(
  const EOBNRHMROMdata *data,           /* Input: data in vector/matrix form to interpolate */
  EOBNRHMROMdata_interp *data_interp    /* Output: interpolated data */
);

/* Functions to load and initalize ROM data */
static INT4 EOBNRv2HMROM_Init_LALDATA(void);
static INT4 EOBNRv2HMROM_Init(const char dir[]);

/* Core functions for waveform reconstruction */
static INT4 EOBNRv2HMROMCore(
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


/************************************************************************************/
/******************************** Internal functions ********************************/

/* Return the closest higher power of 2 */
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/* Arbitrary tuned q-dependent functions by which the frequencies for the 44 and 55 modes have been multiplied (to put the ringdowns at the same Mf). The frequencies stored in the data for the 44 and 55 modes are the rescaled ones, not the original ones. */
static REAL8 Scaling44(const REAL8 q) {
  return 1.-1./4.*exp(-(q-1.)/5.);
}
static REAL8 Scaling55(const REAL8 q) {
  return 1.-1./3.*exp(-(q-1.)/5.);
}

/* Function evaluating eta as a function of q */
static REAL8 EtaOfq(const REAL8 q) {
  return q/(1.+q)/(1.+q);
}
/* Function evaluating delta m/m = (m1-m2)/(m1+m2) as a function of q*/
static REAL8 DeltaOfq(const REAL8 q) {
  return( (q-1.)/(q+1.) );
}

/* Fit of the frequency of the 22 mode at the peak amplitude - from table III in the EOBNRv2HM paper, Pan&al 1106 */
static double omega22peakOfq(const double q) {
  double eta = EtaOfq(q);
  return 0.2733 + 0.2316*eta + 0.4463*eta*eta;
}

/* Amplitude factors scaled out for each mode; this ddoes not include the global amplitude factor scaled out from all modes. */
static REAL8 ModeAmpFactor(const INT4 l, const INT4 m, const REAL8 q) {
  REAL8 eta = EtaOfq(q);
  if( l==2 && m==2 ) return(sqrt(eta));
  else if( l==2 && m==1 ) return( sqrt(eta)*DeltaOfq(q) );
  else if( l==3 && m==3 ) return( sqrt(eta)*DeltaOfq(q) );
  else if( l==4 && m==4 ) return( sqrt(Scaling44(q))*sqrt(eta)*(1.-3.*eta) );
  else if( l==5 && m==5 ) return( pow(Scaling55(q), 1./6)*sqrt(eta)*DeltaOfq(q)*(1.-2.*eta) );
  else {
    fprintf(stderr, "Unknown mode (%d,%d) for the amplitude factor.\n", l, m);
    return(XLAL_FAILURE);
  }
}

/* Function interpolating the data in matrix/vector form produces an interpolated data in the form of SplineLists */
static INT4 Interpolate_Spline_Data(const EOBNRHMROMdata *data, EOBNRHMROMdata_interp *data_interp) {

  SplineList* splinelist;
  gsl_spline* spline;
  gsl_interp_accel* accel;
  gsl_vector* matrixline;
  gsl_vector* vector;
  INT4 j;

  /* Interpolating Camp */
  splinelist = data_interp->Camp_interp;
  for (j = 0; j<nk_amp; j++) {
    matrixline = gsl_vector_alloc(nbwf);
    gsl_matrix_get_row(matrixline, data->Camp, j);

    accel = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
    gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(matrixline, 0), nbwf);

    splinelist = SplineList_AddElementNoCopy(splinelist, spline,  accel, j);
    gsl_vector_free(matrixline);
  }
  data_interp->Camp_interp = splinelist;

  /* Interpolating Cphi */
  splinelist = data_interp->Cphi_interp;
  for (j = 0; j<nk_phi; j++) {
    matrixline = gsl_vector_alloc(nbwf);
    gsl_matrix_get_row(matrixline, data->Cphi, j);

    accel = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
    gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(matrixline, 0), nbwf);

    splinelist = SplineList_AddElementNoCopy(splinelist, spline,  accel, j);
    gsl_vector_free(matrixline);
  }
  data_interp->Cphi_interp = splinelist;

  /* Interpolating shifttime */
  splinelist = data_interp->shifttime_interp;
  vector = data->shifttime;

  accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
  gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(vector, 0), nbwf);

  splinelist = SplineList_AddElementNoCopy(NULL, spline,  accel, 0); /* This SplineList has only 1 element */
  data_interp->shifttime_interp = splinelist;

  /* Interpolating shiftphase */
  splinelist = data_interp->shiftphase_interp;
  vector = data->shiftphase;

  accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
  gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(vector, 0), nbwf);

  splinelist = SplineList_AddElementNoCopy(NULL, spline,  accel, 0); /* This SplineList has only 1 element */
  data_interp->shiftphase_interp = splinelist;

  return XLAL_SUCCESS;
}

/* Function taking as input interpolated data in the form of SplineLists
 * evaluates for a given q the projection coefficients and shifts in time and phase
*/
static INT4 Evaluate_Spline_Data(const REAL8 q, const EOBNRHMROMdata_interp* data_interp, EOBNRHMROMdata_coeff* data_coeff){
  INT4 j;
  SplineList* splinelist;

  /* Evaluating the vector of projection coefficients for the amplitude */
  for (j=0; j<nk_amp; j++) {
    splinelist = SplineList_GetElement(data_interp->Camp_interp, j);
    gsl_vector_set(data_coeff->Camp_coeff, j, gsl_spline_eval(splinelist->spline, q, splinelist->accel));
  }
  /* Evaluating the vector of projection coefficients for the phase */
  for (j=0; j<nk_phi; j++) {
    splinelist = SplineList_GetElement(data_interp->Cphi_interp, j);
    gsl_vector_set(data_coeff->Cphi_coeff, j, gsl_spline_eval(splinelist->spline, q, splinelist->accel));
  }
  /* Evaluating the shift in time */
  splinelist = SplineList_GetElement(data_interp->shifttime_interp, 0); /* This SplineList has only one element */
  *(data_coeff->shifttime_coeff) = gsl_spline_eval(splinelist->spline, q, splinelist->accel);
  /* Evaluating the shift in phase */
  splinelist = SplineList_GetElement(data_interp->shiftphase_interp, 0); /* This SplineList has only one element */
  *(data_coeff->shiftphase_coeff) = gsl_spline_eval(splinelist->spline, q, splinelist->accel);

  return XLAL_SUCCESS;
}

/*
 * Core function for computing the ROM waveform.
 * Evaluates projection coefficients and shifts in time and phase at desired q.
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static INT4 EOBNRv2HMROMCore(
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
  REAL8 q)
{
  INT4 ret = XLAL_SUCCESS;
  INT4 i;
  INT4 j;
  double tpeak22estimate = 0.;
  /* Check output arrays */
  if(!hptilde || !hctilde) XLAL_ERROR(XLAL_EFAULT);
  if(*hptilde || *hctilde)
  {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p\n",(*hptilde),(*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }

  /* Check if the data has been set up */
  if(__lalsim_EOBNRv2HMROM_setup) {
    XLALPrintError("Error: the ROM data has not been set up\n");
    XLAL_ERROR(XLAL_EFAULT);
  }
  /* Set the global pointers to data */
  ListmodesEOBNRHMROMdata* listdata = *__lalsim_EOBNRv2HMROM_data;
  ListmodesEOBNRHMROMdata_interp* listdata_interp = *__lalsim_EOBNRv2HMROM_interp;

  /* Global amplitude prefactor - includes total mass scaling, Fourier scaling, distance scaling, and undoing an additional arbitrary scaling */
  REAL8 Mtot_msol = Mtot_sec / LAL_MTSUN_SI; /* Mtot_msol and M_ROM in units of solar mass */
  REAL8 amp0 = (Mtot_msol/M_ROM) * Mtot_sec * 1.E-16 * 1.E6 * LAL_PC_SI / distance;

  /* Highest allowed geometric frequency for the first mode of listmode in the ROM - used for fRef
   * by convention, we use the first mode of listmode (presumably the 22) for phiref */
  ListmodesEOBNRHMROMdata* listdata_ref = ListmodesEOBNRHMROMdata_GetMode(listdata, listmode[0][0], listmode[0][1]);
  EOBNRHMROMdata* data_ref = listdata_ref->data;
  REAL8 Mf_ROM_max_ref = gsl_vector_get(data_ref->freq, nbfreq-1);
  /* Convert to geometric units for frequency */
  REAL8 fLow_geom = fLow * Mtot_sec;
  REAL8 fHigh_geom = fHigh * Mtot_sec;
  REAL8 fRef_geom = fRef * Mtot_sec;
  REAL8 deltaF_geom = deltaF * Mtot_sec;

  /* Enforce allowed geometric frequency range */
  if (fLow_geom < Mf_ROM_min) { /* Enforce minimal frequency */
    XLALPrintWarning("Starting frequency Mflow=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
    fLow_geom = Mf_ROM_min;
  }
  /* Default highest frequency */
  if (fHigh == 0)
    fHigh_geom = Mf_ROM_max;
  /* In case the user asks for a frequency higher than covered by the ROM, we keep it that way as we will just 0-pad the waveform (and do it anyway for some modes) */
  if (fRef_geom > Mf_ROM_max_ref || fRef_geom == 0)
    fRef_geom = Mf_ROM_max_ref; /* If fRef > fhigh or 0 we reset fRef to default value of cutoff frequency for the first mode of the list (presumably the 22 mode) */
  if (0 < fRef_geom && fRef_geom < Mf_ROM_min) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Setting it to the lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }
  /* Set up output array with size closest power of 2 - fHigh is the upper frequency specified by the user */
  size_t nbpt = NextPow2(fHigh_geom / deltaF_geom) + 1;

  /* Internal storage for the projection coefficients and shifts in time and phase */
  /* Initialized only once, and reused for the different modes */
  EOBNRHMROMdata_coeff *data_coeff = NULL;
  EOBNRHMROMdata_coeff_Init(&data_coeff);
  /* Create spherical harmonic frequency series that will contain the hlm's */
  SphHarmFrequencySeries** hlmsphharmfreqseries = XLALMalloc(sizeof(SphHarmFrequencySeries));
  *hlmsphharmfreqseries = NULL;

  /* GPS time definition - common to all modes */
  LIGOTimeGPS tC;
  XLALGPSAdd(&tC, -1. / deltaF);  /* coalesce at t=0 */

  /* The phase change imposed by phiref, from the phase of the first mode in the list - to be set in the first step of the loop on the modes */
  REAL8 phase_change_ref = 0;

  /* Main loop over the modes */
  for( i=0; i<nbmode; i++ ){
    UINT4 l = listmode[i][0];
    INT4 m = listmode[i][1];

    /* Getting the relevant modes in the lists of data */
    ListmodesEOBNRHMROMdata* listdata_mode = ListmodesEOBNRHMROMdata_GetMode(listdata, l, m);
    ListmodesEOBNRHMROMdata_interp* listdata_interp_mode = ListmodesEOBNRHMROMdata_interp_GetMode(listdata_interp, l, m);

    /* Evaluating the projection coefficients and shift in time and phase */
    ret |= Evaluate_Spline_Data(q, listdata_interp_mode->data_interp, data_coeff);

    /* Evaluating the unnormalized amplitude and unshifted phase vectors for the mode */
    /* Notice a change in convention: B matrices are transposed with respect to the B matrices in SEOBNRROM */
    /* amp_pts = Bamp . Camp_coeff */
    /* phi_pts = Bphi . Cphi_coeff */
    gsl_vector* amp_f = gsl_vector_alloc(nbfreq);
    gsl_vector* phi_f = gsl_vector_alloc(nbfreq);
    gsl_blas_dgemv(CblasNoTrans, 1.0, listdata_mode->data->Bamp, data_coeff->Camp_coeff, 0.0, amp_f);
    gsl_blas_dgemv(CblasNoTrans, 1.0, listdata_mode->data->Bphi, data_coeff->Cphi_coeff, 0.0, phi_f);

    /* The downsampled frequencies for the mode - we undo the rescaling of the frequency for the 44 and 55 modes */
    gsl_vector* freq_ds = gsl_vector_alloc(nbfreq);
    gsl_vector_memcpy(freq_ds, listdata_mode->data->freq);
    if ( l==4 && m==4) gsl_vector_scale( freq_ds, 1./Scaling44(q));
    if ( l==5 && m==5) gsl_vector_scale( freq_ds, 1./Scaling55(q));

    /* Evaluating the shifts in time and phase - conditional scaling for the 44 and 55 modes */
    /* Note: the stored values of 'shifttime' correspond actually to 2pi*Deltat */
    SplineList* shifttime_splinelist = listdata_interp_mode->data_interp->shifttime_interp;
    SplineList* shiftphase_splinelist = listdata_interp_mode->data_interp->shiftphase_interp;
    REAL8 twopishifttime;
    if( l==4 && m==4) {
      twopishifttime = gsl_spline_eval(shifttime_splinelist->spline, q, shifttime_splinelist->accel) * Scaling44(q);
    }
    else if( l==5 && m==5) {
      twopishifttime = gsl_spline_eval(shifttime_splinelist->spline, q, shifttime_splinelist->accel) * Scaling55(q);
    }
    else {
      twopishifttime = gsl_spline_eval(shifttime_splinelist->spline, q, shifttime_splinelist->accel);
    }
    REAL8 shiftphase = gsl_spline_eval(shiftphase_splinelist->spline, q, shiftphase_splinelist->accel);

    /* If first mode in the list, assumed to be the 22 mode, set totalshifttime and phase_change_ref */
    if( i==0 ) {
      if(l==2 && m==2) {
      /* Setup 1d cubic spline for the phase of the 22 mode */
      gsl_interp_accel* accel_phi22 = gsl_interp_accel_alloc();
      gsl_spline* spline_phi22 = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
      gsl_spline_init(spline_phi22, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(phi_f,0), nbfreq);
      /* Compute the shift in time needed to set the peak of the 22 mode roughly at t=0 */
      /* We use the SPA formula tf = -(1/2pi)*dPsi/df to estimate the correspondence between frequency and time */
      /* The frequency corresponding to the 22 peak is omega22peak/2pi, with omega22peak taken from the fit to NR in Pan&al 1106 EOBNRv2HM paper */
      double f22peak = fmin(omega22peakOfq(q)/(2*LAL_PI), Mf_ROM_max_ref); /* We ensure we evaluate the spline within its range */
      tpeak22estimate = -1./(2*LAL_PI) * gsl_spline_eval_deriv(spline_phi22, f22peak, accel_phi22);
      /* Determine the change in phase (to be propagated to all modes) required to have phi22(fRef) = 2*phiRef */
      phase_change_ref = 2*phiRef + (gsl_spline_eval(spline_phi22, fRef_geom, accel_phi22) - (twopishifttime - 2*LAL_PI*tpeak22estimate) * fRef_geom - shiftphase);

      gsl_spline_free(spline_phi22);
      gsl_interp_accel_free(accel_phi22);
      }
      else {
	XLALPrintError("Error: the first mode in listmode must be the 22 mode to set the changes in phase and time \n");
	XLAL_ERROR(XLAL_EFAILED);
      }
    }
    /* Total shift in time, and total change in phase for this mode */
    double totaltwopishifttime = twopishifttime - 2*LAL_PI*tpeak22estimate;
    double constphaseshift = (double) m/listmode[0][1] * phase_change_ref + shiftphase;

    /* Initialize the complex series for the mode - notice that metadata used here is useless, only the one for the final output will matter */
    COMPLEX16FrequencySeries* mode = XLALCreateCOMPLEX16FrequencySeries("mode hlm", &tC, 0.0, deltaF, &lalStrainUnit, nbpt);
    memset(mode->data->data, 0, nbpt * sizeof(COMPLEX16));
    /* Setup 1d cubic spline for the phase and amplitude of the mode */
    gsl_interp_accel* accel_phi = gsl_interp_accel_alloc();
    gsl_interp_accel* accel_amp = gsl_interp_accel_alloc();
    gsl_spline* spline_phi = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
    gsl_spline* spline_amp = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
    gsl_spline_init(spline_phi, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(phi_f,0), nbfreq);
    gsl_spline_init(spline_amp, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(amp_f,0), nbfreq);
    /* Interval in frequency covered by the ROM */
    REAL8 fLow_geom_mode = gsl_vector_get(freq_ds, 0);
    REAL8 fHigh_geom_mode = fmin(gsl_vector_get(freq_ds, nbfreq-1), fHigh_geom);
    /* Initialize the loop - values outside this range in j are 0 by default */
    INT4 jStart = (UINT4) ceil(fLow_geom_mode / deltaF_geom);
    INT4 jStop = (UINT4) ceil(fHigh_geom_mode / deltaF_geom);
    COMPLEX16 *modedata = mode->data->data;
    /* Mode-dependent complete amplitude prefactor */
    REAL8 amp_pre = amp0 * ModeAmpFactor( l, m, q);
    /* Loop on the frequency samples chosen to evaluate the waveform */
    /* We set apart the first and last step to avoid falling outside of the range of the splines by numerical errors */
    REAL8 f, A, phase;

    f = fmax(fLow_geom_mode, jStart*deltaF_geom);
    A = gsl_spline_eval(spline_amp, f, accel_amp);
    phase = -gsl_spline_eval(spline_phi, f, accel_phi) + totaltwopishifttime * f + constphaseshift; /* Minus sign put here, in the internals of the ROM model \Psi = -phase */
    modedata[jStart] = amp_pre * A * cexp(I*phase);

    for (j=jStart+1; j<jStop-1; j++) {
      f = j*deltaF_geom;
      A = gsl_spline_eval(spline_amp, f, accel_amp);
      phase = -gsl_spline_eval(spline_phi, f, accel_phi) + totaltwopishifttime * f + constphaseshift; /* Minus sign put here, in the internals of the ROM model \Psi = -phase */
      modedata[j] = amp_pre * A * cexp(I*phase);
    }

    f = fmin(fHigh_geom_mode, (jStop-1)*deltaF_geom);
    A = gsl_spline_eval(spline_amp, f, accel_amp);
    phase = -gsl_spline_eval(spline_phi, f, accel_phi) + totaltwopishifttime * f + constphaseshift; /* Minus sign put here, in the internals of the ROM model \Psi = -phase */
    modedata[jStop-1] = amp_pre * A * cexp(I*phase);

    /* Add the computed mode to the SphHarmFrequencySeries structure */
    *hlmsphharmfreqseries = XLALSphHarmFrequencySeriesAddMode(*hlmsphharmfreqseries, mode, l, m);

    /* Cleanup for the mode */
    gsl_spline_free(spline_amp);
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(accel_amp);
    gsl_interp_accel_free(accel_phi);
    gsl_vector_free(amp_f);
    gsl_vector_free(phi_f);
    gsl_vector_free(freq_ds);
    XLALDestroyCOMPLEX16FrequencySeries(mode);

  }
  /* Cleanup of the coefficients data structure */
  EOBNRHMROMdata_coeff_Cleanup(data_coeff);

  /* Combining the modes for a hplus, hcross output */
  /* Initialize the complex series hplus, hcross */
  *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, nbpt);
  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, nbpt);

  if (!(hptilde) || !(*hctilde)) XLAL_ERROR(XLAL_EFUNC);
  memset((*hptilde)->data->data, 0, nbpt * sizeof(COMPLEX16));
  memset((*hctilde)->data->data, 0, nbpt * sizeof(COMPLEX16));

  XLALUnitDivide(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
  XLALUnitDivide(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  /* Adding the modes to form hplus, hcross
   * - use of a function that copies XLALSimAddMode but for Fourier domain structures */
  INT4 sym; /* sym will decide whether to add the -m mode (when equatorial symmetry is present) */
  for( i=0; i<nbmode; i++){
    INT4 l = listmode[i][0];
    INT4 m = listmode[i][1];
    COMPLEX16FrequencySeries* mode = XLALSphHarmFrequencySeriesGetMode(*hlmsphharmfreqseries, l, m);
    if ( m==0 ) sym = 0; /* We test for hypothetical m=0 modes */
    else sym = 1;
    FDAddMode( *hptilde, *hctilde, mode, inclination, 0., l, m, sym); /* The phase \Phi is set to 0 - assumes phiRef is defined as half the phase of the 22 mode h22 (or the first mode in the list), not for h = hplus-I hcross */
  }

  /* Destroying the list of frequency series for the modes, including the COMPLEX16FrequencySeries that it contains */
  XLALDestroySphHarmFrequencySeries(*hlmsphharmfreqseries);
  XLALFree(hlmsphharmfreqseries);

  /* Additional complex conjugation of hptilde, hctilde - due to the difference in convention for the Fourier transform between LAL and the ROM internals */
  COMPLEX16* datap = (*hptilde)->data->data;
  COMPLEX16* datac = (*hctilde)->data->data;
  for ( j = 0; j < (INT4) (*hptilde)->data->length; ++j ) {
    datap[j] = conj(datap[j]);
  }
  for ( j = 0; j < (INT4) (*hctilde)->data->length; ++j ) {
    datac[j] = conj(datac[j]);
  }

  return(XLAL_SUCCESS);
}

/* Setup EOBNRv2HMROM model using data files installed in $LAL_DATA_PATH */
static INT4 EOBNRv2HMROM_Init_LALDATA(void) {
  if (!__lalsim_EOBNRv2HMROM_setup) return XLAL_SUCCESS;

  INT4 ret=XLAL_FAILURE;
  char *envpath=NULL;
  char path[32768];
  char *brkt,*word;
  envpath=getenv("LAL_DATA_PATH");
  if(!envpath) {
    XLALPrintError("XLAL Error: the environment variable LAL_DATA_PATH, giving the path to the ROM data, seems undefined\n");
    return(XLAL_FAILURE);
  }
  strncpy(path,envpath,sizeof(path));

  for(word=strtok_r(path,":",&brkt); word; word=strtok_r(NULL,":",&brkt))
  {
    ret = EOBNRv2HMROM_Init(word);
    if(ret == XLAL_SUCCESS) break;
  }
  if(ret!=XLAL_SUCCESS) {
    XLALPrintError("Unable to find EOBNRv2HMROM data files in $LAL_DATA_PATH\n");
    exit(XLAL_FAILURE);
  }
  __lalsim_EOBNRv2HMROM_setup = ret;
  return(ret);
}

/* Setup EOBNRv2HMROM model using data files installed in dir */
static INT4 EOBNRv2HMROM_Init(const char dir[]) {
  if(!__lalsim_EOBNRv2HMROM_setup) {
    XLALPrintError("Error: EOBNRHMROMdata was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  INT4 ret = XLAL_SUCCESS;
  ListmodesEOBNRHMROMdata* listdata = *__lalsim_EOBNRv2HMROM_data;
  ListmodesEOBNRHMROMdata_interp* listdata_interp = *__lalsim_EOBNRv2HMROM_interp;

  for (UINT4 j=0; j<EOBNRV2_ROM_NUM_MODES_MAX; j++) {

    EOBNRHMROMdata* data = NULL;
    EOBNRHMROMdata_Init(&data);
    ret |= Read_Data_Mode( dir, listmode[j], data);
    if (!ret) {
      listdata = ListmodesEOBNRHMROMdata_AddModeNoCopy( listdata, data, listmode[j][0], listmode[j][1]);

      EOBNRHMROMdata_interp* data_interp = NULL;
      EOBNRHMROMdata_interp_Init(&data_interp);
      ret |= Interpolate_Spline_Data(data, data_interp);
      if (!ret) listdata_interp = ListmodesEOBNRHMROMdata_interp_AddModeNoCopy( listdata_interp, data_interp, listmode[j][0], listmode[j][1]);
    }
  }

  __lalsim_EOBNRv2HMROM_setup = ret;
  if (!ret) {
    *__lalsim_EOBNRv2HMROM_data = listdata;
    *__lalsim_EOBNRv2HMROM_interp = listdata_interp;
  }
  return(ret);
}


/******************************************************************************************/
/****************************** Waveform generation function ******************************/

/* Compute waveform in LAL format */
INT4 XLALSimIMREOBNRv2HMROM(
  struct tagCOMPLEX16FrequencySeries **hptilde, /* Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /* Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /* Half-phase of the 22 mode (~orbital phase) at reference frequency */
  REAL8 deltaF,                                 /* Sampling frequency (Hz) */
  REAL8 fLow,                                   /* Start frequency (Hz) - 0 defaults to the lowest Mf of the lowest m>=2 mode */
  REAL8 fHigh,                                  /* End frequency (Hz) - 0 defaults to the highest Mf of the highest m mode */
  REAL8 fRef,                                   /* Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /* Distance of source (m) */
  REAL8 inclination,                            /* Inclination of source (rad) */
  REAL8 m1SI,                                   /* Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /* Mass of companion 2 (kg) */
  const int higherModesFlag)                    /* flag to indicate use of higher modes */
{

  /* Set the number of modes depending on whether the user wants higher order modes */
  if ( higherModesFlag == 0 )
  {
    nbmode = 1;
  }
  else if ( higherModesFlag == 1 )
  {
    nbmode = EOBNRV2_ROM_NUM_MODES_MAX;
  }
  else
  {
    XLALPrintError( "Higher mode flag appears to be uninitialised "
        "(expected 0 or 1, but got %d\n)", higherModesFlag );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Get masses in terms of solar mass */
  REAL8 mass1 = m1SI / LAL_MSUN_SI;
  REAL8 mass2 = m2SI / LAL_MSUN_SI;
  REAL8 Mtot = mass1 + mass2;
  REAL8 q = fmax(mass1/mass2, mass2/mass1);    /* Mass-ratio >1 by convention*/
  REAL8 Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  if ( q > q_max ) {
    XLALPrintError( "XLAL Error - %s: q out of range!\nEOBNRv2HMROM is only available for a mass ratio in the range q <= %g.\n", __func__, q_max);
    XLAL_ERROR( XLAL_EDOM );
  }

  /* Set up (load and build interpolation) ROM data if not setup already */
  EOBNRv2HMROM_Init_LALDATA();

  /* Call the Core function */
  INT4 retcode = EOBNRv2HMROMCore(hptilde,hctilde,
            phiRef, deltaF, fLow, fHigh, fRef, distance, inclination, Mtot_sec, q);

  return(retcode);
}
