/********************************** <lalVerbatim file="ExcessPowerHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _EPSEARCH_H
#define _EPSEARCH_H

#include <lal/ExcessPower.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TFTransform.h>
#include <lal/Window.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID(EPSEARCHH, "$Id$");


typedef struct
tagEPSearchParams {
	CHAR                 *channelName;
	BOOLEAN               printSpectrum;
	BOOLEAN               useOverWhitening;	
	size_t                eventLimit;
	UINT4                 windowLength;
	UINT4                 windowShift;
	REAL8                 lnalphaThreshold;
	AvgSpecMethod         method;
	CreateTFTilingIn      tfTilingInput;
	TFPlaneParams         tfPlaneParams;
	WindowType            windowType;
	INT4                  minFreqBins;
} EPSearchParams;

int
XLALEPSearch(
	SnglBurstTable  **burstEvent,
	const REAL4TimeSeries  *tseries,
	EPSearchParams   *params
);

void
LALEPConditionData(
	LALStatus        *status,
	REAL4TimeSeries  *series,
	REAL8             flow,
	REAL8             resampledeltaT,
	ResampleTSFilter  resampleFiltType,
	INT4              corruption
);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
