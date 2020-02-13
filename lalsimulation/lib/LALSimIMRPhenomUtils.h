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

// the following definitions are used in XLALSimIMRPhenomPv3HMComputeWignerdElements
#define FIVE_OVER_16 0.3125
#define THREE_OVER_EIGHT 0.375
#define ONE_OVER_TWO 0.5
#define ONE_OVER_EIGHT 0.125
#define ONE_OVER_FOUR 0.25
#define SQRT_6 2.44948974278317788
#define FIFTEEN_OVER_32 0.46875
#define THREE_OVER_16 0.1875
#define ONE_OVER_32 0.03125
#define mSQRT_6_OVER_32 -0.07654655446197431
#define SQRT_15_OVER_32 0.12103072956898178
#define mSQRT_15_OVER_32 -0.12103072956898178
#define SQRT_5_OVER_16 0.13975424859373686
#define SQRT_6_OVER_32 0.07654655446197431
#define SQRT_10_OVER_32 0.09882117688026186
#define mSQRT_10_OVER_32 -0.09882117688026186
#define SQRT_30_OVER_16 0.3423265984407288

#define SEVEN_OVER_16 0.4375
#define SEVEN_OVER_32 0.21875
#define ONE_OVER_128 0.0078125
#define THIRTYFIVE_OVER_128 0.2734375
#define ONE_OVER_16 0.0625
#define SQRT_2_OVER_64 0.02209708691207961
#define mSQRT_2_OVER_64 -0.02209708691207961
#define SQRT_7_OVER_16 0.16535945694153692
#define SQRT_14_OVER_32 0.11692679333668567
#define SQRT_70_16 0.5229125165837972
#define NINE_OVER_32 0.28125
#define SQRT_7_OVER_32 0.08267972847076846
#define SQRT_35_OVER_4 1.479019945774904

// void XLALSimPhenomUtilsTest(void);

double XLALSimPhenomUtilsMftoHz(REAL8 Mf, REAL8 Mtot_Msun);
double XLALSimPhenomUtilsHztoMf(REAL8 fHz, REAL8 Mtot_Msun);

double XLALSimPhenomUtilsFDamp0(REAL8 Mtot_Msun, REAL8 distance);

REAL8 XLALSimPhenomUtilsPhenomPv2FinalSpin(const REAL8 m1, const REAL8 m2, const REAL8 chi1_l, const REAL8 chi2_l, const REAL8 chip);

REAL8 XLALSimPhenomUtilsPhenomPv3HMFinalSpin(const REAL8 m1, const REAL8 m2, const REAL8 chi1x, const REAL8 chi1y, const REAL8 chi1z, const REAL8 chi2x, const REAL8 chi2y, const REAL8 chi2z);

double XLALSimPhenomUtilsIMRPhenomDFinalMass(REAL8 m1, REAL8 m2, REAL8 chi1z, REAL8 chi2z);

REAL8 XLALSimPhenomUtilsChiP(const REAL8 m1, const REAL8 m2, const REAL8 s1x, const REAL8 s1y, const REAL8 s2x, const REAL8 s2y);

int XLALSimPhenomUtilsPhenomPv3HMWignerdElement(REAL8 *wig_d, UINT4 ell, INT4 mprime, INT4 mm, REAL8 b);

/**
 * a strcut to keep the wigner-d matrix elements
 */
typedef struct tagIMRPhenomPv3HMWignderStruct
{
    REAL8 d22[2][5]; /**< wigner-d matrix elements for ell=2, mprime=+/-2 and positive[0] and negative[1] mprime */
    REAL8 d21[2][5]; /**< wigner-d matrix elements for ell=2, mprime=+/-1 and positive[0] and negative[1] mprime */
    REAL8 d33[2][7]; /**< wigner-d matrix elements for ell=3, mprime=+/-3 and positive[0] and negative[1] mprime */
    REAL8 d32[2][7]; /**< wigner-d matrix elements for ell=3, mprime=+/-2 and positive[0] and negative[1] mprime */
    REAL8 d44[2][9]; /**< wigner-d matrix elements for ell=4, mprime=+/-4 and positive[0] and negative[1] mprime */
    REAL8 d43[2][9]; /**< wigner-d matrix elements for ell=4, mprime=+/-3 and positive[0] and negative[1] mprime */
} IMRPhenomPv3HMWignderStruct;

/**
 * a strcut to keep the complex exponential terms of the alpha precession angle
 */
typedef struct tagIMRPhenomPv3HMAlphaStruct
{
    COMPLEX16 alpha2[5]; /**< optimised elements for complex exponential of the alpha angle for ell = 2*/
    COMPLEX16 alpha3[7]; /**< optimised elements for complex exponential of the alpha angle for ell = 3*/
    COMPLEX16 alpha4[9]; /**< optimised elements for complex exponential of the alpha angle for ell = 4*/
} IMRPhenomPv3HMAlphaStruct;

/**
 * a strcut to keep the spherical harmonic terms
 */
typedef struct tagIMRPhenomPv3HMYlmStruct
{
    COMPLEX16 Ylm2[2][5]; /**< optimised elements Ylms for ell = 2. [0] for Ylm and [1] for conj(Ylm) */
    COMPLEX16 Ylm3[2][7]; /**< optimised elements Ylms for ell = 3. [0] for Ylm and [1] for conj(Ylm) */
    COMPLEX16 Ylm4[2][9]; /**< optimised elements Ylms for ell = 4. [0] for Ylm and [1] for conj(Ylm) */
} IMRPhenomPv3HMYlmStruct;

int XLALSimIMRPhenomPv3HMComputeWignerdElements(IMRPhenomPv3HMWignderStruct **p, UNUSED UINT4 ell, UNUSED INT4 mprime, UNUSED REAL8 b);

int XLALSimIMRPhenomPv3HMComputeAlphaElements(IMRPhenomPv3HMAlphaStruct **p, UINT4 ell, REAL8 alpha);

IMRPhenomPv3HMYlmStruct *XLALSimIMRPhenomPv3HMComputeYlmElements(REAL8 theta, REAL8 phi, INT4 ell);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMUTILS_H */
