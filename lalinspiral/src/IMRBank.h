#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOMetadataTables.h>
#include <gsl/gsl_rng.h>

#ifndef _IMRBANK_H
#define _IMRBANK_H


/** New structure to hold moments full arrays!! */
typedef struct
  {
  REAL8Vector * minus3[25];
  REAL8Vector * plus3[25];
  REAL8Vector * logplus3[25];
  REAL8Vector * logminus3[25];
  REAL8Vector * logsqplus3[25];
  REAL8Vector * logsqminus3[25];
  REAL8FrequencySeries *psd;
  UINT4 length;
  REAL8 deltaF;
  REAL8 flow;
  LIGOTimeGPS epoch;
  REAL8 f0;
  LALUnit sampleUnits;
  } IMRBankCumulativeNoiseMoments;

/** A more convenient metric type */
typedef struct
  {
  REAL8 data[3][3];
  REAL8 m1;
  REAL8 m2;
  REAL8 M;
  REAL8 eta;
  REAL8 tau0;
  REAL8 tau3;
  } IMRBankMetric;

/** This holds a box in m1,m2 */
typedef struct tagIMRBankMassRegion
  {
  REAL8 mbox[3];
  struct tagIMRBankMassRegion *next;
  } IMRBankMassRegion;

int XLALTileIMRBankMassRegion(InspiralCoarseBankIn *in, SnglInspiralTable **first);

#endif
