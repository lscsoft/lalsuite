#include <math.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>

#include <LALAdaptiveRungeKutta4.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSTPN_TEST_ENERGY			1025
#define LALSTPN_TEST_OMEGADOT		1026
#define LALSTPN_TEST_COORDINATE	1027
#define LALSTPN_TEST_OMEGANAN		1028

#define LALSTPN_DERIVATIVE_OMEGANONPOS	1030
#define LALSTPN_DERIVATIVE_COORDINATE		1031

void
LALSTPNAdaptiveWaveformEngineFrameless( LALStatus *status,
    REAL4Vector *signalvec1,REAL4Vector *signalvec2,
    UINT4 *countback,
    InspiralTemplate *params,InspiralInit *paramsInit
    );

#ifdef  __cplusplus
#pragma {
}
#endif

