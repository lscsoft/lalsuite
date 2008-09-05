#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOMetadataTables.h>
#include <gsl/gsl_rng.h>

/* New structure to hold moments full arrays!! */
typedef struct tagLALIMRBankCumulativeNoiseMoments
  {
  REAL8Vector * minus3[15];
  REAL8Vector * plus3[15];
  REAL8Vector * logplus3[15];
  REAL8Vector * logminus3[15];
  REAL8Vector * logsqplus3[15];
  REAL8Vector * logsqminus3[15];
  REAL8FrequencySeries *psd;
  UINT4 length;
  REAL8 deltaF;
  REAL8 flow;
  LIGOTimeGPS epoch;
  REAL8 f0;
  LALUnit sampleUnits;
  } LALIMRBankCumulativeNoiseMoments;

typedef struct tagLALIMRBankMetric
  {
  REAL8 data[2][2];
  } LALIMRBankMetric;

typedef struct tagLALIMRBankMassRegion
  {
  REAL8 mbox[3];
  struct tagLALIMRBankMassRegion *next;
  } LALIMRBankMassRegion;

int XLALCreateLALIMRBankCumulativeNoiseMoments(LALIMRBankCumulativeNoiseMoments *moments, REAL8FrequencySeries *psd, REAL8 flow);


int XLALComputeIMRBankCumulativeNoiseMoment(LALIMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag);


int XLALDestroyLALIMRBankCumulativeNoiseMoments(LALIMRBankCumulativeNoiseMoments *moments);

REAL8 XLALComputeLALIMRBankMetricMassMass(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, LALIMRBankCumulativeNoiseMoments *I);

REAL8 XLALComputeLALIMRBankMetricEtaEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, LALIMRBankCumulativeNoiseMoments *I);

REAL8 XLALComputeLALIMRBankMetricMassEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, LALIMRBankCumulativeNoiseMoments *I);

int XLALComputeLALIMRBankMetric(REAL8 mass1, REAL8 mass2, LALIMRBankCumulativeNoiseMoments *moments, LALIMRBankMetric *metric);

INT4 XLALTileLALIMRBankMassRegion(InspiralCoarseBankIn *in, SnglInspiralTable **first);
