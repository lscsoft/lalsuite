#ifndef _LALSTPNWAVEFORM2_H
#define _LALSTPNWAVEFORM2_H

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

extern int newswitch;

void
LALSTPNAdaptiveWaveformEngine( LALStatus *status,
                							 REAL4Vector *signalvec1,REAL4Vector *signalvec2,
                							 REAL4Vector *a,REAL4Vector *ff,REAL8Vector *phi,REAL4Vector *shift,
                							 UINT4 *countback,
                							 InspiralTemplate *params,InspiralInit *paramsInit
														 );

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSTPNWAVEFORM2_H */
