/*-----------------------------------------------------------------------
 *
 * File Name: ExcessPower.h
 *
 * Author: Eanna Flanagan
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * ExcessPower.h
 *
 * SYNOPSIS
 * #include "ExcessPower.h"
 *
 * DESCRIPTION
 * Error codes, typedefs, and protypes for the functions chisqCdf() and
 * related functions
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */


#ifndef _EXCESSPOWER_H
#define _EXCESSPOWER_H


#include "LALStdlib.h"
#include "LALDatatypes.h"
#include "TFTransform.h"
#include "LALRCSID.h"

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (EXCESSPOWERH, "$Id$");


#define EXCESSPOWER_ENULLP       1
#define EXCESSPOWER_EPOSARG      2
#define EXCESSPOWER_EPOW2        4
#define EXCESSPOWER_EMALLOC      8
#define EXCESSPOWER_EINCOMP      16
#define EXCESSPOWER_EORDER       32
#define EXCESSPOWER_ENONNULL     64


#define EXCESSPOWER_MSGENULLP    "Null pointer"
#define EXCESSPOWER_MSGEPOSARG   "Arguments must be non-negative"
#define EXCESSPOWER_MSGEPOW2     "Length of supplied data must be a power of 2"
#define EXCESSPOWER_MSGEMALLOC   "Malloc failure"
#define EXCESSPOWER_MSGEINCOMP   "Incompatible arguments"
#define EXCESSPOWER_MSGEORDER    "Routines called in illegal order"
#define EXCESSPOWER_MSGENONNULL  "Null pointer expected"




typedef struct
tagTFTile
{
  INT4                             fstart;
  INT4                             fend;
  INT4                             tstart;
  INT4                             tend;
  INT4                             whichPlane;
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
AddWhiteNoise (
               Status                               *status,
               COMPLEX8Vector                       *v,
               REAL8                                noiseLevel
               );


void
CreateTFTiling (
                Status                              *status,
                TFTiling                            **tfTiling,
                CreateTFTilingIn                    *input
                );


void
DestroyTFTiling (
                 Status                             *status,
                 TFTiling                           **tfTiling
                 );


void
ComputeTFPlanes (
                 Status                             *status,
                 TFTiling                           *tfTiling,
                 COMPLEX8FrequencySeries            *freqSeries
                 );


void
ComputeExcessPower (
                    Status                          *status,
                    TFTiling                        *tfTiling,
                    ComputeExcessPowerIn            *input
                    );


void
SortTFTiling (
              Status                                *status,
              TFTiling                              *tfTiling
              );


void
CountEPEvents (
               Status                               *status,
               INT4                                 *numEvents,
               TFTiling                             *tfTiling,
               REAL8                                alphaThreshold
               );


void
ComputeLikelihood (
                   Status                           *status,
                   REAL8                            *lambda,
                   TFTiling                         *tfTiling
                   );


void 
PrintTFTileList (
                 Status                             *status,
                 FILE                               *fp,
                 TFTiling                           *tfTiling,
                 INT4                               maxTiles
                 );

void 
PrintTFTileList1 (
                 Status                             *status,
                 FILE                               *fp,
                 TFTiling                           *tfTiling,
                 INT4                               maxTiles
                 );


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif




