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
#include <lal/LALInspiral.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID(EPSEARCHH, "$Id$");

typedef struct
tagMergerDurationLimits {
	REAL4      high;
	REAL4      low;
} MergerDurationLimits;


REAL4
XLALMergerFrequency( 
    InspiralTemplate *params  
    );

REAL4
XLALQnrFrequency( 
    InspiralTemplate *params  
    );

REAL4
XLALMergerEnergyFraction(
	InspiralTemplate *params
	);

void 
XLALMergerDuration(
	InspiralTemplate *params,
	MergerDurationLimits *dur
	);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
