#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <lal/FrequencySeries.h>
#include <sys/types.h>

typedef struct tagWS {
    size_t n;
    fftwf_plan plan;
    COMPLEX8 *zf;
    COMPLEX8 *zt;
} WS;

WS *XLALCreateSBankWorkspaceCache(void);
void XLALDestroySBankWorkspaceCache(WS *workspace_cache);
REAL8 XLALInspiralSBankComputeMatch(const COMPLEX8FrequencySeries *inj, const COMPLEX8FrequencySeries *tmplt, WS *workspace_cache);
