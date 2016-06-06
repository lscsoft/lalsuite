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

typedef struct {
   SkyPosition assumeNSpos;
   REAL8 *assumeNScosi;
   REAL8 *assumeNSpsi;
   REAL8 *assumeNSGWfreq;
   REAL8 *assumeNSorbitP;
   REAL8 *assumeNSasini;
   LIGOTimeGPS *assumeNSorbitTp;
   LIGOTimeGPS *assumeNSrefTime;
} assumeNSparams;

SFTCatalog * findSFTdata(const UserInput_t *params);
MultiSFTVector * extractSFTband(const SFTCatalog *catalog, const REAL8 minfbin, const REAL8 maxfbin);
MultiSFTVector * getMultiSFTVector(UserInput_t *params, const REAL8 minfbin, const REAL8 maxfbin);
MultiSFTVector * generateSFTdata(UserInput_t *uvar, const MultiLALDetector *detectors, const EphemerisData *edat, const INT4 maxbinshift, const gsl_rng *rng);
REAL4VectorAligned * coherentlyAddSFTs(const MultiSFTVector *multiSFTvector, const MultiSSBtimes *multissb, const MultiAMCoeffs *multiAMcoeffs, const LIGOTimeGPSVector *jointTimestamps, const REAL4VectorAlignedArray *backgroundRatio, const INT4 cosiSign, const assumeNSparams *NSparams, const UserInput_t *params, REAL4VectorAligned *backgroundScaling);
REAL4VectorAligned * convertSFTdataToPowers(const SFTVector *sfts, const UserInput_t *params, const REAL8 normalization);
REAL4VectorAligned * readInSFTs(UserInput_t *params, const REAL8 normalization, const REAL8 minfbin, const REAL8 maxfbin);

REAL4TimeSeries * computeNSfreqTS(const PulsarParams *pulsarParams, LIGOTimeGPS tepoch, REAL8 Tsft, REAL8 SFToverlap, REAL8 duration);

MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTCatalog(const SFTCatalog *catalog);
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTs(const MultiSFTVector *multiSFTvector, const UserInput_t *params);
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSegmentsFile(const LALStringVector *filenames, const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 dur);

LIGOTimeGPSVector * jointTimestampsFromMultiTimestamps(const MultiLIGOTimeGPSVector *multiTimestamps);
MultiLIGOTimeGPSVector * squeezeMultiTimestamps(const MultiLIGOTimeGPSVector *multiTS);

INT4Vector * markBadSFTs(const REAL4VectorAligned *tfdata, const UserInput_t *params);
INT4Vector * existingSFTs(const REAL4VectorAligned *tfdata, const UINT4 numffts);
void removeBadSFTs(REAL4VectorAligned *tfdata, const INT4Vector *badsfts);

INT4 slideTFdata(REAL4VectorAligned *output, const UserInput_t *params, const REAL4VectorAligned *tfdata, const INT4Vector *binshifts);
INT4 tfRngMeans(REAL4VectorAligned *output, const REAL4VectorAligned *tfdata, const UINT4 numffts, const UINT4 numfbins, const UINT4 blksize);
INT4 tfMeanSubtract(REAL4VectorAligned *tfdata, const REAL4VectorAligned *rngMeans, const REAL4VectorAligned *backgroundScaling, const UINT4 numffts, const UINT4 numfbins);
INT4 tfWeight(REAL4VectorAligned *output, const REAL4VectorAligned *tfdata, REAL4VectorAligned *rngMeans, REAL4VectorAligned *antPatternWeights, const REAL4VectorAligned *backgroundScaling, const INT4Vector *indexValuesOfExistingSFTs, const UserInput_t *params);
INT4 replaceTFdataWithSubsequentTFdata(REAL4VectorAlignedArray *tfdataarray, const UINT4 numffts);
INT4 checkBackgroundScaling(const REAL4VectorAligned *background, const REAL4VectorAligned *backgroundScaling, const INT4Vector *SFTexistVector);

REAL8 determineSumOfWeights(const REAL4VectorAligned *antweightssq, const REAL4VectorAligned *rngMeanssq);

INT4 printSFTtimestamps2File(const MultiSFTVector *multiSFTvector, const UserInput_t *params);


#endif
