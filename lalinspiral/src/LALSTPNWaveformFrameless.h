#include <math.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>

#include <lal/LALAdaptiveRungeKutta4.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

void LALSTPNFramelessAdaptiveWaveformEngine(LALStatus *status,
    REAL4Vector *signalvec1, REAL4Vector *signalvec2, UINT4 *countback,
    InspiralTemplate *params,InspiralInit *paramsInit
    );

int XLALSTPNFramelessAdaptiveWaveformEngine(REAL4Vector *signalvec1,
    REAL4Vector *signalvec2, UINT4 *countback, InspiralTemplate *params, 
    InspiralInit *paramsInit
    );

#ifdef  __cplusplus
#pragma {
}
#endif

