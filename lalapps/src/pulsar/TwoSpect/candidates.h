/*
*  Copyright (C) 2010, 2011 Evan Goetz
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

#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

#include <lal/RealFFT.h>
#include "TwoSpectTypes.h"

candidateVector * createcandidateVector(const UINT4 length);
candidateVector * resizecandidateVector(candidateVector *vector, const UINT4 length);
candidateVector * keepMostSignificantCandidates(const candidateVector *input, const UserInput_t *params);
void destroycandidateVector(candidateVector *vector);

void loadCandidateData(candidate* output,
                       const REAL8 fsig,
                       const REAL8 period,
                       const REAL8 moddepth,
                       const REAL4 ra,
                       const REAL4 dec,
                       const REAL8 statval,
                       const REAL8 h0,
                       const REAL8 prob,
                       const INT4 proberrcode,
                       const REAL8 normalization,
                       const INT4 templateVectorIndex,
                       const BOOLEAN lineContamination);

INT4 analyzeOneTemplate(candidate *output,
                        const candidate *input,
                        const ffdataStruct *ffdata,
                        const REAL4VectorAligned *aveNoise,
                        const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                        const UserInput_t *params,
                        const REAL4FFTPlan *plan,
                        const gsl_rng *rng,
                        const BOOLEAN exactflag);
INT4 bruteForceTemplateSearch(candidate *output,
                              const candidate input,
                              const TwoSpectParamSpaceSearchVals *paramspace,
                              const UserInput_t *params,
                              const REAL4VectorAligned *ffdata,
                              const REAL4VectorAligned *aveNoise,
                              const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                              const REAL4FFTPlan *secondFFTplan,
                              const gsl_rng *rng,
                              const BOOLEAN useExactTemplates);
INT4 bruteForceTemplateTest(candidateVector **output,
                            const candidate input,
                            const TwoSpectParamSpaceSearchVals *paramspace,
                            const UserInput_t *params,
                            const REAL4VectorAligned *ffdata,
                            const REAL4VectorAligned *aveNoise,
                            const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                            const REAL4FFTPlan *secondFFTplan,
                            const gsl_rng *rng,
                            const BOOLEAN useExactTemplates);
INT4 templateSearch_scox1Style(candidateVector **output,
                               const REAL8 fminimum,
                               const REAL8 fspan,
                               const REAL8 period,
                               const REAL8 asini,
                               const REAL8 asinisigma,
                               const SkyPosition skypos,
                               const UserInput_t *params,
                               const REAL4VectorAligned *ffdata,
                               const REAL4VectorAligned *aveNoise,
                               const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                               const REAL4VectorSequence *trackedlines,
                               const REAL4FFTPlan *secondFFTplan,
                               const gsl_rng *rng,
                               const BOOLEAN useExactTemplates);
INT4 templateSearch_fixedDf(candidateVector **output,
                               const LALStringVector *dffixed,
                               const REAL8 fminimum,
                               const REAL8 fspan,
                               const REAL8 period,
                               const SkyPosition skypos,
                               const UserInput_t *params,
                               const REAL4VectorAligned *ffdata,
                               const REAL4VectorAligned *aveNoise,
                               const REAL4VectorAligned *aveTFnoisePerFbinRatio,
			       const REAL4VectorSequence *trackedlines,
                               const REAL4FFTPlan *secondFFTplan,
                               const gsl_rng *rng,
                               const BOOLEAN useExactTemplates);
INT4 clusterCandidates(candidateVector **output,
                       const candidateVector *input,
                       const ffdataStruct *ffdata,
                       const UserInput_t *params,
                       const REAL4VectorAligned *ffplanenoise,
                       const REAL4VectorAligned *fbinaveratios,
                       const gsl_rng *rng,
                       const BOOLEAN exactflag);
INT4 testIHScandidates(candidateVector **output,
                       const candidateVector *ihsCandidates,
                       const ffdataStruct *ffdata,
                       const REAL4VectorAligned *aveNoise,
                       const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                       const SkyPosition pos,
                       const UserInput_t *params,
                       const gsl_rng *rng);
INT4 testTwoSpectTemplateVector(candidateVector *output,
                                const TwoSpectTemplateVector *templateVec,
                                const ffdataStruct *ffdata,
                                const REAL4VectorAligned *aveNoise,
                                const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                                const SkyPosition skypos,
                                const UserInput_t *params,
                                const gsl_rng *rng,
                                const UINT4 templateLen);
INT4 analyzeCandidatesTemplateFromVector(candidateVector *output,
                                         const candidateVector *input,
                                         const TwoSpectTemplateVector *vector,
                                         const ffdataStruct *ffdata,
                                         const REAL4VectorAligned *aveNoise,
                                         const REAL4VectorAligned *aveTFnoisePerFbinRatio,
                                         const UserInput_t *params,
                                         const gsl_rng *rng,
                                         const UINT4 templateLen);

INT4 writeCandidateVector2File(const CHAR *outputfile, const candidateVector *input);

REAL8 maxModDepth(const REAL8 period, const REAL8 cohtime);
REAL8 minPeriod(const REAL8 moddepth, const REAL8 cohtime);
REAL8 calculateR(const REAL4VectorAligned *ffdata, const TwoSpectTemplate *template, const REAL4VectorAligned *noise, const REAL4VectorAligned *fbinaveratios);

#endif


