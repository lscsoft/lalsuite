/*
*  Copyright (C) 2014 Matthew Pitkin
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
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Header file for the helper functions for the parameter estimation code for known pulsar
 * searches using the nested sampling algorithm.
 */

#ifndef _PPE_UTILS_H
#define _PPE_UTILS_H

#include "pulsar_parameter_estimation_nested.h"
#include "ppe_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

void compute_variance( LALInferenceIFOData *data, LALInferenceIFOModel *model );
COMPLEX16Vector *subtract_running_median( COMPLEX16Vector *data );

/* functions for finding change points in the data */
UINT4Vector *get_chunk_lengths( LALInferenceIFOModel *ifo, INT4 chunkMax );
UINT4Vector *chop_n_merge( LALInferenceIFOData *data, INT4 chunkMin, INT4 chunkMax );
UINT4Vector *chop_data( gsl_vector_complex *data, INT4 chunkMin );
UINT4 find_change_point( gsl_vector_complex *data, REAL8 *logodds, INT4 chunkMin );
void rechop_data( UINT4Vector *segs, INT4 chunkMax, INT4 chunkMin );
void merge_data( COMPLEX16Vector *data, UINT4Vector *segs );

void gzip_output( LALInferenceRunState *runState );
INT4 count_csv( CHAR *csvline );
INT4 recognised_parameter( CHAR *parname );
void check_and_add_fixed_variable( LALInferenceVariables *vars, const char *name, void *value, LALInferenceVariableType type );

TimeCorrectionType XLALAutoSetEphemerisFiles( CHAR *efile, CHAR *sfile,
                                              CHAR *tfile,
                                              BinaryPulsarParams pulsar,
                                              INT4 gpsstart, INT4 gpsend );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_UTILS_H */
