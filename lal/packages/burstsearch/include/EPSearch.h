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
#include <lal/Window.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID(EPSEARCHH, "$Id$");


typedef struct
tagEPSearchParams {
	CHAR                 *channelName;
	BOOLEAN               printSpectrum;
	INT4                  eventLimit;
	UINT4                 windowLength;
	UINT4                 windowShift;
	REAL8                 alphaThreshold;
	AvgSpecMethod         method;
        CreateTFTilingIn      tfTilingInput;
        TFPlaneParams         tfPlaneParams;
	ComputeExcessPowerIn  compEPInput;
	WindowType            windowType;
} EPSearchParams;

void
EPSearch(
	LALStatus        *status,
	REAL4TimeSeries  *tseries,
	EPSearchParams   *params,
	SnglBurstTable  **burstEvent
);

void
EPConditionData(
	LALStatus        *status,
	REAL4TimeSeries  *series,
	REAL4             flow,
	REAL8             resampledeltaT,
	ResampleTSFilter  resampleFiltType,
	INT4              corruption
);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
