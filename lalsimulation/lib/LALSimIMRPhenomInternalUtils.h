#ifndef _LALSIM_IMR_PHENOMINTERNALUTILS_H
#define _LALSIM_IMR_PHENOMINTERNALUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

#include <stdbool.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/XLALError.h>
#include <gsl/gsl_math.h>

// UNUSED void PhenomInternal_UtilsTest(void);

UNUSED bool PhenomInternal_approx_equal(REAL8 x, REAL8 y, REAL8 epsilon);

UNUSED void PhenomInternal_nudge(REAL8 *x, REAL8 X, REAL8 epsilon);

UNUSED size_t PhenomInternal_NextPow2(const size_t n);

UNUSED int PhenomInternal_AlignedSpinEnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1z, REAL8 *chi2z);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMINTERNALUTILS_H */
