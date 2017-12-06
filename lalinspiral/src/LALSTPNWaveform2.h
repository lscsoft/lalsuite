#ifndef _LALSTPNWAVEFORM2_H
#define _LALSTPNWAVEFORM2_H

#include <math.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>

#include <lal/LALAdaptiveRungeKutta4.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

extern int newswitch;

void
LALSTPNAdaptiveWaveformEngine( LALStatus *status,
                							 REAL4Vector *signalvec1,REAL4Vector *signalvec2,
                							 REAL4Vector *a,REAL4Vector *ff,REAL8Vector *phi,REAL4Vector *shift,
                							 UINT4 *countback,
                							 InspiralTemplate *params,InspiralInit *paramsInit
														 );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LALSTPNWAVEFORM2_H */
