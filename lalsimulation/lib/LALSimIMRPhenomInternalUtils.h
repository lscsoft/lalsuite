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
#include <lal/Sequence.h>
#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>

// UNUSED void PhenomInternal_UtilsTest(void);

/**
 * Tolerance used below which numbers are treated as zero for the calculation of atan2
 */
#define PhenomInternal_MAX_TOL_ATAN 1.0e-15

UNUSED REAL8 PhenomInternal_atan2tol(REAL8 x, REAL8 y, REAL8 tol);

UNUSED bool PhenomInternal_approx_equal(REAL8 x, REAL8 y, REAL8 epsilon);

UNUSED void PhenomInternal_nudge(REAL8 *x, REAL8 X, REAL8 epsilon);

UNUSED size_t PhenomInternal_NextPow2(const size_t n);

UNUSED int PhenomInternal_AlignedSpinEnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1z, REAL8 *chi2z);

UNUSED int PhenomInternal_PrecessingSpinEnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1x, REAL8 *chi1y, REAL8 *chi1z, REAL8 *chi2x, REAL8 *chi2y, REAL8 *chi2z);

UNUSED void PhenomInternal_ComputeIMRPhenomPv3CartesianToPolar(REAL8 *polar, REAL8 *azimuthal, REAL8 *magnitude, REAL8 x, REAL8 y, REAL8 z);

UNUSED double PhenomInternal_EradRational0815_s(double eta, double s);
UNUSED double PhenomInternal_EradRational0815(double eta, double chi1, double chi2);

UNUSED double PhenomInternal_OrbAngMom3PN(
    const double f_orb_hz,
    const double m1,
    const double m2,
    const double s1x,
    const double s1y,
    const double s1z,
    const double s2x,
    const double s2y,
    const double s2z,
    const double f_0,
    const int ExpansionOrder);

UNUSED int PhenomInternal_IMRPhenomHMFDAddMode(
    COMPLEX16FrequencySeries *hptilde,
    COMPLEX16FrequencySeries *hctilde,
    COMPLEX16FrequencySeries *hlmtilde,
    REAL8 theta,
    REAL8 phi,
    INT4 l,
    INT4 m,
    INT4 sym);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMINTERNALUTILS_H */
