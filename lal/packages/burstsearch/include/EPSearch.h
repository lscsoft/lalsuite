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
	REAL4Window          *window;
	UINT4                 windowShift;
	REAL8                 lnalphaThreshold;
	AvgSpecMethod         method;
	CreateTFTilingIn      tfTilingInput;
	TFPlaneParams         tfPlaneParams;
} EPSearchParams;

int
XLALEPSearch(
	SnglBurstTable  **burstEvent,
	const REAL4TimeSeries  *tseries,
	EPSearchParams   *params
);

int
XLALEPConditionData(
	REAL4TimeSeries  *series,
	REAL8             flow,
	REAL8             resampledeltaT,
	INT4              corruption
);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
