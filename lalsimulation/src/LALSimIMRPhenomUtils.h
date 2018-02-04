#ifndef _LALSIM_IMR_PHENOMUTILS_H
#define _LALSIM_IMR_PHENOMUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>

// void XLALSimPhenomUtilsTest(void);

double XLALSimPhenomUtilsMftoHz(REAL8 Mf, REAL8 Mtot_Msun);
double XLALSimPhenomUtilsHztoMf(REAL8 fHz, REAL8 Mtot_Msun);

double XLALSimPhenomUtilsFDamp0(REAL8 Mtot_Msun, REAL8 distance);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMUTILS_H */
