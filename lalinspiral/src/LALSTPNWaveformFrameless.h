#include <math.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>

#include <lal/LALAdaptiveRungeKuttaIntegrator.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

void LALSTPNFramelessAdaptiveWaveformEngine(LALStatus *status,
    REAL4Vector *signalvec1, REAL4Vector *signalvec2, UINT4 *countback,
    InspiralTemplate *params,InspiralInit *paramsInit
    );

int XLALSTPNFramelessAdaptiveWaveformEngine(REAL4Vector *signalvec1,
    REAL4Vector *signalvec2, UINT4 *countback, InspiralTemplate *params, 
    InspiralInit *paramsInit
    );

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
