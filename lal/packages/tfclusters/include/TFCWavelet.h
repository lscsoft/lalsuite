#ifndef _TFCWAVELET_H_
#define _TFCWAVELET_H_

/************************************ <lalVerbatim file="TFCWaveletH">
Author: Sylvestre, J.
$Id$
************************************* </lalVerbatim> */

#include <lal/LALDatatypes.h>
#include <lal/LALStatusMacros.h>
#include <lal/TFClusters.h>

#define TFCWAVELETH_ENULLP 1000
#define TFCWAVELETH_ELEN 1001
#define TFCWAVELETH_ENNULLP 1002
#define TFCWAVELETH_EMALLOC 1003
#define TFCWAVELETH_EINCOMP 1004

#define TFCWAVELETH_MSGENULLP "Null pointer"
#define TFCWAVELETH_MSGELEN   "Invalid input size"
#define TFCWAVELETH_MSGENNULLP "Non-null pointer"
#define TFCWAVELETH_MSGEMALLOC "Memory allocation error"
#define TFCWAVELETH_MSGEINCOMP "Incompatible time series and TF parameters"

/* Daubechies coefficients for size 4, 12 & 20 */
#if 0
static const REAL4 daub2[4]={0.4829629131445341,0.8365163037378079,
                    0.2241438680420134,-0.1294095225512604};
static const REAL4 daub6[12]={0.111540743350, 0.494623890398, 0.751133908021,
                      0.315250351709,-0.226264693965,-0.129766867567,
                      0.097501605587, 0.027522865530,-0.031582039318,
                      0.000553842201, 0.004777257511,-0.001077301085};
static const REAL4 daub10[20]={0.026670057901, 0.188176800078, 0.527201188932,
                      0.688459039454, 0.281172343661,-0.249846424327,
                      -0.195946274377, 0.127369340336, 0.093057364604,
                      -0.071394147166,-0.029457536822, 0.033212674059,
                      0.003606553567,-0.010733175483, 0.001395351747,
                      0.001992405295,-0.000685856695,-0.000116466855,
                      0.000093588670,-0.000013264203};
#endif

typedef struct tagTFCWavelet {

  REAL4Vector *smooth;
  REAL4Vector *detail;

  REAL4Vector *input;

  REAL4Vector *wavelet;
} TFCWavelet;

typedef struct tagTFCWParams {

  INT4 timeBins;
  INT4 freqBins;

  REAL8 deltaT;

  REAL8 minScale;
  REAL8 maxScale;

  REAL4Vector *wavelet;
} TFCWParams;

void 
LALWaveletFilter (
		  LALStatus *status,
		  TFCWavelet *wave
		  );


void 
LALComputeWaveletSpectrogram (
			      LALStatus *status, 
			      Spectrogram *out, 
			      TFCWParams *tspec, 
			      REAL4TimeSeries *tseries
			      );


#endif
