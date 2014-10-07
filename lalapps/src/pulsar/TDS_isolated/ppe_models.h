/*******************************************************************************
  Matt Pitkin, Colin Gill, John Veitch - 2011

  ppe_models.h

  Header file for ppe_models.c

*******************************************************************************/

/*
  Author:
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Header file for the signal models functions used in parameter
 * estimation code for known pulsar searches using the nested sampling
 * algorithm.
 */

#ifndef _PPE_MODELS_H
#define _PPE_MODELS_H

#include "pulsar_parameter_estimation_nested.h"

#ifdef __cplusplus
extern "C" {
#endif

/* global variables */

/** The inverse of the factorials of the numbers 0 to 6. */
static const REAL8 inv_fact[7] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* model functions */
void get_pulsar_model( LALInferenceModel *model, LALInferenceIFOData *data );

REAL8 rescale_parameter( LALInferenceModel *model, LALInferenceIFOData *data, const CHAR *parname );

void pulsar_model( BinaryPulsarParams params, LALInferenceIFOData *data );

REAL8Vector *get_phase_model( BinaryPulsarParams params, LALInferenceIFOData *data, REAL8 freqFactor,
                              UINT4 downsampled );

REAL8Vector *get_ssb_delay( BinaryPulsarParams pars, LIGOTimeGPSVector *datatimes, EphemerisData *ephem,
                            TimeCorrectionData *tdat, TimeCorrectionType ttype, LALDetector *detector,
                            REAL8 interptime );

REAL8Vector *get_bsb_delay( BinaryPulsarParams pars, LIGOTimeGPSVector *datatimes, REAL8Vector *dts,
                            EphemerisData *ephem );

void get_triaxial_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data );

void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data );

void get_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data );

REAL8 noise_only_model( LALInferenceRunState *runState );

REAL8 get_phase_mismatch( REAL8Vector *phi1, REAL8Vector *phi2, LIGOTimeGPSVector *ts );

void get_earth_pos_vel( EarthState *earth, EphemerisData *ephem, LIGOTimeGPS *t );

/* functions to convert between parameters */
void phi0_psi_transform( REAL8 phi0, REAL8 psi, REAL8 *phi0prime, REAL8 *psiprime );

void inverse_phi0_psi_transform( REAL8 phi0prime, REAL8 psiprime, REAL8 *phi0, REAL8 *psi );

void invert_source_params( BinaryPulsarParams *params );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_MODELS_H */
