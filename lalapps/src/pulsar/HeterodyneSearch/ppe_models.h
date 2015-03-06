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
 * \ingroup lalapps_pulsar_HeterodyneSearch
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

void add_pulsar_parameter( LALInferenceVariables *var, PulsarParameters *params, const CHAR *parname );

void add_variable_parameter( PulsarParameters *params, LALInferenceVariables *var, const CHAR *parname, LALInferenceParamVaryType vary );

void pulsar_model( PulsarParameters *params, LALInferenceIFOModel *ifo );

void set_nonGR_model_parameters( PulsarParameters *pars, char* nonGRmodel );

REAL8Vector *get_phase_model( PulsarParameters *params, LALInferenceIFOModel *ifo, REAL8 freqFactor );

REAL8Vector *get_ssb_delay( PulsarParameters *pars, LIGOTimeGPSVector *datatimes, EphemerisData *ephem,
                            TimeCorrectionData *tdat, TimeCorrectionType ttype, LALDetector *detector,
                            REAL8 interptime );

REAL8Vector *get_bsb_delay( PulsarParameters *pars, LIGOTimeGPSVector *datatimes, REAL8Vector *dts,
                            EphemerisData *ephem );

void get_triaxial_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo );

void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo );

void get_amplitude_model( PulsarParameters *pars, LALInferenceIFOModel *ifo );

REAL8 get_phase_mismatch( REAL8Vector *phi1, REAL8Vector *phi2, LIGOTimeGPSVector *ts );

void get_earth_pos_vel( EarthState *earth, EphemerisData *ephem, LIGOTimeGPS *t );

void response_lookup_table( REAL8 t0, LALDetAndSource detNSource, INT4 timeSteps, REAL8 avedt,
                            REAL8Vector *a1, REAL8Vector *b1, REAL8Vector *a2, REAL8Vector *b2,
                            REAL8Vector *a3, REAL8Vector *b3 );

/* functions to convert between parameters */

void invert_source_params( PulsarParameters *params );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_MODELS_H */
