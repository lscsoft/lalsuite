/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpExch.h
 *
 * Author: Allen, B. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPEXCH_H
#define _FINDCHIRPEXCH_H

#include <lal/LALDatatypes.h>
#include <lal/Comm.h>
#include <lal/DataBuffer.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FINDCHIRPEXCHH, "$Id$");

#define FINDCHIRPEXCH_ENULL 1
#define FINDCHIRPEXCH_ENNUL 2
#define FINDCHIRPEXCH_ENOBJ 4
#define FINDCHIRPEXCH_EHAND 8

#define FINDCHIRPEXCH_MSGENULL "Null pointer"
#define FINDCHIRPEXCH_MSGENNUL "Non-null pointer"
#define FINDCHIRPEXCH_MSGENOBJ "Invalid number of objects"
#define FINDCHIRPEXCH_MSGEHAND "Wrong handshake"

typedef enum
{
  ExchDataSegment,
  ExchInspiralBankIn,
  ExchInspiralTemplate,
  ExchInspiralEvent,
  ExchFinished
}
ExchObjectType;

typedef struct
tagExchParams
{
  ExchObjectType exchObjectType;
  INT4           send;
  INT4           numObjects;
  INT4           partnerProcNum;
}
ExchParams;

typedef enum
{
  M1AndM2,
  TotalMassAndEta,
  TotalMassAndMu
}
InputMasses;


typedef enum
{
  Best,
  TaylorTime20,
  TaylorFreq20,
  Papprox20
}
InspiralMethod;


typedef struct
tagInspiralTemplate
{
  INT4           number;
  REAL4          mass1; 
  REAL4          mass2;
  REAL4          spin1[3];
  REAL4          spin2[3];
  REAL4          inclination;
  REAL4          eccentricity;
  REAL8          totalMass; 
  REAL8          mu; 
  REAL8          eta;
  REAL8          fLower;
  REAL8          fCutoff;
  REAL8          tSampling;
  REAL8          phaseShift;
  INT4           nStartPad;
  INT4           nEndPad;
  InputMasses    massChoice;
  InspiralMethod method;
 }
InspiralTemplate;

/* the various interferometer codes */
typedef enum
{
  Caltech40m, Hanford4km, Hanford2km, Livingston4km, GEO600m, TAMA300m
}
Detector;

/* input for specifying a template bank */
typedef struct
tagInspiralBankIn
{
  REAL4          mMin;
  REAL4          mMax;
  REAL4          ffCoarse;
  REAL4          ffFine;
  Detector       detector;
  InspiralMethod method;
}
InspiralBankIn;


typedef struct
tagInspiralEvent
{
  LIGOTimeGPS      time;
  InspiralTemplate tmplt;
  REAL4            snrsq;
  REAL4            chisq;
  REAL4            sigma;
}
InspiralEvent;




void
LALInitializeExchange (
    LALStatus      *status,
    ExchParams **exchParamsOut,
    ExchParams  *exchParamsInp,
    INT4         myProcNum
    );

void
LALFinalizeExchange (
    LALStatus      *status,
    ExchParams **exchParams
    );

void
LALExchangeDataSegment (
    LALStatus      *status,
    DataSegment *segment,
    ExchParams  *exchParams
    );

void
LALExchangeInspiralBankIn (
    LALStatus         *status,
    InspiralBankIn *bankIn,
    ExchParams     *exchParams
    );

void
LALExchangeInspiralTemplate (
    LALStatus           *status,
    InspiralTemplate *tmplt,
    ExchParams       *exchParams
    );

void
LALExchangeInspiralEvent (
    LALStatus        *status,
    InspiralEvent *event,
    ExchParams    *exchParams
    );


#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPEXCH_H */
