/** \defgroup InspiralInjectionParams
 * \ingroup inject
 * \author D. Brown, J. Creighton, S. Fairhurst, G. Jones, E. Messaritaki
 *
 * \brief Module for generating randomly distributed inspiral parameters
 *

 *
 */

/** \file InspiralInjecionParams.h
 *  \ingroup InspiralInjectionParams
 * \date $Date$
 *
 *
 */

/** enum containing the different ways in which the distance to
    injections can be distributed */
typedef enum
{
  unknownDistanceDist,
  distFromSourceFile,
  uniformDistance,
  uniformLogDistance,
  uniformVolume
}
DistanceDistribution;

/** enum containing the different ways in which the sky location of
    injections can be distributed */
typedef enum
{
  unknownLocationDist,
  locationFromSourceFile,
  locationFromExttrigFile,
  uniformSkyLocation,
  fixedSkyLocation
}
SkyLocationDistribution;

/** enum containing the different ways in which the masses of
    injections can be distributed */
typedef enum
{
  unknownMassDist,
  massFromSourceFile,
  massFromNRFile,
  uniformTotalMass,
  uniformComponentMass,
  logComponentMass,
  gaussianMassDist,
  uniformTotalMassRatio,
  logMassUniformTotalMassRatio,
}
MassDistribution;

/** enum containing the different ways in which the inclinations of
    injections can be distributed */
typedef enum
{
  unknownInclDist,
  uniformInclDist,
  gaussianInclDist,
  fixedInclDist
}
InclDistribution;

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( NRWAVEIOC, "$Id$");

SimInspiralTable* XLALRandomInspiralTime( SimInspiralTable *inj,
  RandomParams *randParams,
  LIGOTimeGPS startTime,
  REAL4 timeWindow );

SimInspiralTable* XLALRandomInspiralDistance( SimInspiralTable *inj,
    RandomParams *randParams,
    DistanceDistribution dDist,
    REAL4  distMin,
    REAL4  distMax );

SimInspiralTable* XLALRandomInspiralSkyLocation( SimInspiralTable *inj,
    RandomParams *randParams);

void XLALRandomInspiralMilkywayLocation(
    REAL8 *rightAscension,
    REAL8 *declination,
    REAL8 *distance,
    RandomParams *randParams);

SimInspiralTable* XLALRandomInspiralOrientation( SimInspiralTable *inj,
    RandomParams *randParams,
    InclDistribution iDist,
    REAL4   inclinationPeak  );

SimInspiralTable* XLALRandomInspiralMasses( SimInspiralTable *inj,
    RandomParams *randParams,
    MassDistribution mDistr,
    REAL4  mass1Min,
    REAL4  mass1Max,
    REAL4  mass2Min,
    REAL4  mass2Max,
    REAL4  minTotalMass,
    REAL4  maxTotalMass  );

SimInspiralTable* XLALGaussianInspiralMasses( SimInspiralTable *inj,
    RandomParams *randParams,
    REAL4  mass1Min,
    REAL4  mass1Max,
    REAL4  mass1Mean,
    REAL4  mass1Std,
    REAL4  mass2Min,
    REAL4  mass2Max,
    REAL4  mass2Mean,
    REAL4  mass2Std);

SimInspiralTable* XLALRandomInspiralTotalMassRatio( SimInspiralTable *inj,
    RandomParams *randParams,
    MassDistribution mDist,
    REAL4  minTotalMass,
    REAL4  maxTotalMass,
    REAL4  minMassRatio,
    REAL4  maxMassRatio);

SimInspiralTable* XLALRandomInspiralSpins( SimInspiralTable *inj,
    RandomParams *randParams,
    REAL4  spin1Min,
    REAL4  spin1Max,
    REAL4  spin2Min,
    REAL4  spin2Max,
    REAL4  kappa1Min,
    REAL4  kappa1Max,
    REAL4  abskappa1Min,
    REAL4  abskappa1Max);

SimInspiralTable* XLALRandomNRInjectTotalMass(
    SimInspiralTable *inj,
    RandomParams *randParams,
    REAL4  minTotalMass,
    REAL4  maxTotalMass,
    SimInspiralTable *nrInjParams
    );

SimInspiralTable *XLALInspiralSiteTimeAndDist(
    SimInspiralTable  *inj,
    LALDetector       *detector,
    LIGOTimeGPS       *endTime,
    REAL4             *effDist);

SimInspiralTable *XLALPopulateSimInspiralSiteInfo(
    SimInspiralTable           *inj );


COMPLEX8FrequencySeries *generateActuation(
    COMPLEX8FrequencySeries *resp,
    REAL4                    ETMcal,
    REAL4                    pendF,
    REAL4                    pendQ );


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

