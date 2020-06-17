/*
 * LALSimInspiralEOBPostAdiabatic.h
 */

#include "LALSimIMRSpinEOB.h"

#ifndef LALSimInspiralEOBPostAdiabatic_h
#define LALSimInspiralEOBPostAdiabatic_h

int
XLALSimInspiralEOBPostAdiabatic(
    REAL8Vector ** dynamics,
    REAL8 deltaT,
    const REAL8 m1SI,
    const REAL8 m2SI,
    const REAL8 spin1z,
    const REAL8 spin2z,
    UINT4 SpinAlignedEOBversion,
    SpinEOBParams *seobParams,
    EOBNonQCCoeffs *nqcCoeffs
);

double
XLALSimInspiralEOBPACalculateSymmetricMassRatio(
    const REAL8 q
);

double
XLALSimInspiralEOBPostAdiabaticX1(
    const REAL8 nu
);

double
XLALSimInspiralEOBPostAdiabaticX2(
    const REAL8 nu
);

double
XLALSimInspiralEOBPostAdiabaticTimeUnitsFactor(
    REAL8 M
);

double
XLALSimInspiralEOBPostAdiabaticDynr0Kepler(
    REAL8 f0
);

REAL8
XLALSimInspiralEOBPostAdiabaticTotalSpin(
    REAL8 q,
    REAL8 a1,
    REAL8 a2
);

REAL8
XLALSimInspiralEOBPostAdiabaticFinalRadius(
    REAL8 q,
    REAL8 a1,
    REAL8 a2
);

double
XLALSimInspiralEOBPostAdiabaticFitGlobalc3(
    REAL8 nu,
    REAL8 a1,
    REAL8 a2
);

double
XLALSimInspiralEOBPostAdiabaticz3(
    const REAL8 nu
);

int
XLALSimInspiralEOBPACalculateRadialGrid(
    REAL8Vector *rVec,
    LALDict *LALParams
);

double
XLALSimInspiralEOBPostAdiabaticdprstarFunc(
    REAL8 prstar_sol,
    void *params
);

double
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
    REAL8 pphi_sol,
    void *params
);

double
XLALSimInspiralEOBPostAdiabaticj0Func(
    REAL8 j0_sol,
    void *params
);

struct PostAdiabaticRootSolveParams
{
    REAL8 r;
    REAL8 rc;
    REAL8 drcBydr;
    REAL8 uc2;
    REAL8 ducBydr;
    REAL8 prstar;
    REAL8 dprstarBydr;
    REAL8 phi;
    REAL8 pphi;
    REAL8 dpphiBydr;
    REAL8 A;
    REAL8 dA;
    REAL8 B;
    REAL8 dAuc2Bydr;
    REAL8 HeffOrb;
    REAL8 dGBydr;
    REAL8 dGBydprstar;
    REAL8 dr;
    REAL8 omega;
    REAL8 csi;
    SpinEOBParams *seobParams;
    EOBNonQCCoeffs *nqcCoeffs;
    LALDict *LALParams;
};

struct PostAdiabaticRoot
{
    REAL8 root;
    INT4 status;
    INT4 nIter;
};

int
XLALSimInspiralEOBPostAdiabaticRootFinder(
    struct PostAdiabaticRoot *result,
    double (*Func)(REAL8, void *),
    struct PostAdiabaticRootSolveParams *params,
    REAL8 x_lower,
    REAL8 x_upper,
    REAL8 absTol,
    REAL8 relTol
);

REAL8Vector
XLALReverseREAL8Vector(
    REAL8Vector *Vec
);

REAL8Vector
XLALOffsetREAL8Vector(
    REAL8Vector *Vec,
    REAL8 offset
);

REAL8Vector
XLALPostAdiabaticSplineDerivative(
    REAL8Vector *VecX,
    REAL8Vector *VecY
);

REAL8Vector
XLALFDDerivative1Order2(
    REAL8Vector *XVec,
    REAL8Vector *YVec
);

REAL8Vector
XLALFDDerivative1Order4(
    REAL8Vector *XVec,
    REAL8Vector *YVec
);

REAL8Vector
XLALFDDerivative1Order6(
    REAL8Vector *XVec,
    REAL8Vector *YVec
);

REAL8Vector
XLALFDDerivative1Order8(
    REAL8Vector *XVec,
    REAL8Vector *YVec
);

REAL8Vector
XLALCumulativeIntegral3(
    REAL8Vector *XVec,
    REAL8Vector *YVec
);

REAL8
XLALSimInspiralEOBPAHamiltonianWrapper(
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    SpinEOBHCoeffs *seobCoeffs,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPAHamiltonianDerivative(
    REAL8 h,
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    SpinEOBHCoeffs *seobCoeffs,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
    REAL8 h,
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    SpinEOBHCoeffs *seobCoeffs,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPAFluxWrapper(
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    REAL8 omega,
    REAL8 H,
    SpinEOBParams *seobParams,
    EOBNonQCCoeffs *nqcCoeffs,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPANewtonianj0(
    REAL8 r
);

#endif