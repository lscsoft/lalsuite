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




















