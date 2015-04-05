/*
*  Copyright (C) 2010 -- 2014 Evan Goetz
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

#ifndef __SFTFUNCTIONS_H__
#define __SFTFUNCTIONS_H__

#include <lal/SSBtimes.h>
#include <lal/LALComputeAM.h>
#include <lal/LALInitBarycenter.h>

#include <gsl/gsl_sf_trig.h>

#include "TwoSpectTypes.h"
#include "vectormath.h"

SFTCatalog * findSFTdata(UserInput_t *params);
MultiSFTVector * extractSFTband(SFTCatalog *catalog, REAL8 minfbin, REAL8 maxfbin);
MultiSFTVector * getMultiSFTVector(UserInput_t *params, REAL8 minfbin, REAL8 maxfbin);
MultiSFTVector * generateSFTdata(UserInput_t *uvar, MultiLALDetector *detectors, EphemerisData *edat, INT4 maxbinshift, gsl_rng *rng);
REAL4Vector * coherentlyAddSFTs(MultiSFTVector *multiSFTs, MultiSSBtimes *multissb, MultiAMCoeffs *multiAMcoefficients, INT4 cosiSign, REAL8 *assumeNScosi, REAL8 *assumeNSpsi, UserInput_t *params);
REAL4Vector * convertSFTdataToPowers(SFTVector *sfts, UserInput_t *params, REAL8 normalization);
REAL4Vector * readInSFTs(UserInput_t *params, REAL8 normalization, REAL8 minfbin, REAL8 maxfbin);

MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTCatalog(SFTCatalog *catalog);
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTs(MultiSFTVector *multiSFTvector, UserInput_t *params);
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSegmentsFile(LALStringVector *filenames, REAL8 t0, REAL8 Tsft, REAL8 SFToverlap, REAL8 dur);

INT4Vector * markBadSFTs(REAL4Vector *tfdata, UserInput_t *params);
INT4Vector * existingSFTs(REAL4Vector *tfdata, UINT4 numffts);
void removeBadSFTs(REAL4Vector *tfdata, INT4Vector *badsfts);

INT4 slideTFdata(REAL4Vector *output, UserInput_t *params, REAL4Vector *tfdata, INT4Vector *binshifts);
INT4 tfRngMeans(REAL4Vector *output, REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize);
INT4 tfMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, INT4 numffts, INT4 numfbins);
INT4 tfWeight(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4VectorAligned *antPatternWeights, INT4Vector *indexValuesOfExistingSFTs, UserInput_t *params);

REAL8 determineSumOfWeights(REAL4Vector *antweightssq, REAL4Vector *rngMeanssq);

#endif
