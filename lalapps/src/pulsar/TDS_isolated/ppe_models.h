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
 * \ingroup lalapps_pulsar
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

/* model functions */
void get_pulsar_model( LALInferenceModel *model );

REAL8 rescale_parameter( LALInferenceModel *model, LALInferenceIFOModel *ifo, const CHAR *parname );

void pulsar_model( BinaryPulsarParams params, LALInferenceIFOModel *ifo );

REAL8Vector *get_phase_model( BinaryPulsarParams params, LALInferenceIFOModel *ifo, REAL8 freqFactor );

REAL8Vector *get_ssb_delay( BinaryPulsarParams pars, LIGOTimeGPSVector *datatimes, EphemerisData *ephem,
                            TimeCorrectionData *tdat, TimeCorrectionType ttype, LALDetector *detector,
                            REAL8 interptime );

REAL8Vector *get_bsb_delay( BinaryPulsarParams pars, LIGOTimeGPSVector *datatimes, REAL8Vector *dts,
                            EphemerisData *ephem );

void get_triaxial_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo );

void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo );

void get_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo );

REAL8 get_phase_mismatch( REAL8Vector *phi1, REAL8Vector *phi2, LIGOTimeGPSVector *ts );

void get_earth_pos_vel( EarthState *earth, EphemerisData *ephem, LIGOTimeGPS *t );

void response_lookup_table( REAL8 t0, LALDetAndSource detNSource, INT4 timeSteps, REAL8Vector *a1,
                            REAL8Vector *b1, REAL8Vector *a2, REAL8Vector *b2, REAL8Vector *a3,
                            REAL8Vector *b3 );

/* functions to convert between parameters */
void phi0_psi_transform( REAL8 phi0, REAL8 psi, REAL8 *phi0prime, REAL8 *psiprime );

void inverse_phi0_psi_transform( REAL8 phi0prime, REAL8 psiprime, REAL8 *phi0, REAL8 *psi );

void invert_source_params( BinaryPulsarParams *params );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_MODELS_H */
