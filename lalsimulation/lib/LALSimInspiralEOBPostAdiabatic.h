/*
 * LALSimInspiralEOBPostAdiabatic.h
 */

#ifndef LALSimInspiralEOBPostAdiabatic_h
#define LALSimInspiralEOBPostAdiabatic_h

int
XLALSimInspiralEOBPostAdiabatic(
    REAL8TimeSeries ** dynamics,
    REAL8 deltaT,
    const REAL8 m1SI,
    const REAL8 m2SI
);

double
XLALSimInspiralEOBPostAdiabaticSymmetricMassRatio(
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
XLALSimInspiralEOBPostAdiabaticMetricS(
    REAL8 *A,
    REAL8 *B,
    REAL8 *dA,
    REAL8 *d2A,
    REAL8 *dB,
    REAL8 r,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadius(
    REAL8 *rc,
    REAL8 *drcBydr,
    REAL8 *d2rcBydr2,
    REAL8 r,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadiusLO(
    REAL8 *rc,
    REAL8 *drcBydr,
    REAL8 *d2rcBydr2,
    REAL8 r,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticMetricA5PNlog(
    REAL8 *Aorb,
    REAL8 *dAorbBydu,
    REAL8 *d2AorbBydu2,
    REAL8 r,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticsGSDynamics(
    REAL8 *ggm,
    REAL8 r,
    REAL8 rc,
    REAL8 drcBydr,
    REAL8 prstar,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticHamiltonianS(
    REAL8 *H,
    REAL8 *Heff,
    REAL8 *Heff_orb,
    REAL8 *dHeff_dr,
    REAL8 *dHeff_dprstar,
    REAL8 *dHeff_dpphi,
    REAL8 *d2Heff_dprstar20,
    REAL8 r,
    REAL8 rc,
    REAL8 drc_dr,
    REAL8 pphi,
    REAL8 prstar,
    REAL8 S,
    REAL8 Sstar,
    REAL8 A,
    REAL8 dA,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticFluxS(
    REAL8 *Flux,
    REAL8 x,
    REAL8 Omega,
    REAL8 r_omega,
    REAL8 E,
    REAL8 Heff,
    REAL8 jhat,
    REAL8 r,
    REAL8 pr_star,
    REAL8 ddotr,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticFlmNewtonian(
    REAL8 *Nlm,
    REAL8 x,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticTlmFlux(
    REAL8 *MTlm,
    REAL8 w
);

int
XLALSimInspiralEOBPostAdiabaticGetWavFlmS(
    REAL8 *rholm,
    REAL8 *flm,
    REAL8 x,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticGetWavFlmSSLO(
    REAL8 *rholm,
    REAL8 *flm,
    REAL8 x,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticWavFlm(
    REAL8 *rholm,
    REAL8 *flm,
    REAL8 x,
    LALDict *LALParams
);

double
XLALSimInspiralEOBPostAdiabaticEulerLog(
    REAL8 x,
    INT4 m
);

int
XLALSimInspiralEOBPostAdiabaticHorizonFluxS(
    REAL8 *hatFH,
    REAL8 x,
    REAL8 Heff,
    REAL8 jhat,
    LALDict *LALParams
);

int
XLALSimInspiralEOBPostAdiabaticHorizonFlux(
    REAL8 *hatFH,
    REAL8 x,
    REAL8 Heff,
    REAL8 jhat,
    LALDict *LALParams
);

double
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
    REAL8 prstar_sol,
    void *params
);

double
XLALSimInspiralEOBPostAdiabaticdprstarFunc(
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
    REAL8 pphi;
    REAL8 dpphiBydr;
    REAL8 A;
    REAL8 dA;
    REAL8 B;
    REAL8 dAuc2Bydr;
    REAL8 HeffOrb;
    REAL8 dGBydr;
    REAL8 dGBydprstar;
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
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPAHamiltonianDerivative(
    REAL8 h,
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
    REAL8 h,
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPAFluxWrapper(
    REAL8 r,
    REAL8 prstar,
    REAL8 pphi,
    REAL8 omega,
    REAL8 H,
    LALDict *LALParams
);

REAL8
XLALSimInspiralEOBPANewtonianj0(
    REAL8 r
);

#endif