/********************************** <lalVerbatim file="ExcessPowerHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _EXCESSPOWER_H
#define _EXCESSPOWER_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/TFTransform.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (EXCESSPOWERH, "$Id$");

/******** <lalErrTable file="ExcessPowerHErrTab"> ********/
#define EXCESSPOWERH_ENULLP       1
#define EXCESSPOWERH_EPOSARG      2
#define EXCESSPOWERH_EPOW2        4
#define EXCESSPOWERH_EMALLOC      8
#define EXCESSPOWERH_EINCOMP      16
#define EXCESSPOWERH_EORDER       32
#define EXCESSPOWERH_ENONNULL     64
#define EXCESSPOWERH_ETILES       65
#define EXCESSPOWERH_EDELF        128


#define EXCESSPOWERH_MSGENULLP    "Null pointer"
#define EXCESSPOWERH_MSGEPOSARG   "Arguments must be non-negative"
#define EXCESSPOWERH_MSGEPOW2     "Length of supplied data must be a power of 2"
#define EXCESSPOWERH_MSGEMALLOC   "Malloc failure"
#define EXCESSPOWERH_MSGEINCOMP   "Incompatible arguments"
#define EXCESSPOWERH_MSGEORDER    "Routines called in illegal order"
#define EXCESSPOWERH_MSGENONNULL  "Null pointer expected"
#define EXCESSPOWERH_MSGETILES    "Malloc failed while assigning memory for a tile"
#define EXCESSPOWERH_MSGEDELF     "Inconsistent deltaF in spectrum and data"
/******** </lalErrTable> ********/



typedef struct
tagTFTile
{
  INT4                             fstart;
  INT4                             fend;
  INT4                             tstart;
  INT4                             tend;
  INT4                             whichPlane;
  REAL8                            deltaT;      /* deltaF will always be 1/deltaT     */
  REAL8                            excessPower;
  REAL8                            alpha;
  REAL8                            weight;
  BOOLEAN                          firstCutFlag;  
  struct tagTFTile                 *nextTile;
}
TFTile;


typedef struct
tagTFTiling
{
  TFTile                           *firstTile;  /* linked list of Tiles */
  INT4                             numTiles;
  COMPLEX8TimeFrequencyPlane       **tfp;       /* Vector of pointers */
  ComplexDFTParams                 **dftParams; /* Vector of pointers */
  INT4                             numPlanes;
  BOOLEAN                          planesComputed;      
  BOOLEAN                          excessPowerComputed;
  BOOLEAN                          tilesSorted;
}
TFTiling;


typedef struct
tagCreateTFTilingIn
{
  INT4                             overlapFactor;
  INT4                             minFreqBins;
  INT4                             minTimeBins;
  REAL8                            flow;    /* lowest freq to search       */
  REAL8                            deltaF;  
  INT4                             length;  
}
CreateTFTilingIn;


typedef struct
tagComputeExcessPowerIn
{
  REAL8                            numSigmaMin;
  REAL8                            alphaDefault;                             
}
ComputeExcessPowerIn;


void
LALAddWhiteNoise (
               LALStatus                               *status,
               COMPLEX8Vector                       *v,
               REAL8                                noiseLevel
               );


void
LALCreateTFTiling (
                LALStatus                              *status,
                TFTiling                            **tfTiling,
                CreateTFTilingIn                    *input
                );


void
LALDestroyTFTiling (
                 LALStatus                             *status,
                 TFTiling                           **tfTiling
                 );


void
LALComputeTFPlanes (
                 LALStatus                             *status,
                 TFTiling                           *tfTiling,
                 COMPLEX8FrequencySeries            *freqSeries
                 );


void
LALComputeExcessPower (
                    LALStatus                          *status,
                    TFTiling                        *tfTiling,
                    ComputeExcessPowerIn            *input
                    );


void
LALSortTFTiling (
              LALStatus                                *status,
              TFTiling                              *tfTiling
              );


void
LALCountEPEvents (
               LALStatus                               *status,
               INT4                                 *numEvents,
               TFTiling                             *tfTiling,
               REAL8                                alphaThreshold
               );


void
LALComputeLikelihood (
                   LALStatus                           *status,
                   REAL8                            *lambda,
                   TFTiling                         *tfTiling
                   );


void 
LALPrintTFTileList (
                 LALStatus                             *status,
                 FILE                               *fp,
                 TFTiling                           *tfTiling,
                 INT4                               maxTiles
                 );

void 
PrintTFTileList1 (
                 LALStatus                             *status,
                 FILE                               *fp,
                 TFTiling                           *tfTiling,
                 INT4                               maxTiles
                 );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif




