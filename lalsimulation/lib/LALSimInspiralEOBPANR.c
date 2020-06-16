/*
 * EOB Metric potentials A(r), B(r), and their derivatives, spin version
 */
int
XLALSimInspiralEOBPostAdiabaticMetricS(
	REAL8 *A,
	REAL8 *B,
	REAL8 *dA,
	REAL8 *d2A,
	REAL8 *dB,
	REAL8 r,
	LALDict *LALParams)
{
    UNUSED REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
    INT4 useTidal = XLALDictLookupINT4Value(LALParams, "useTidal");

    REAL8 u = 1. / r;
    REAL8 u2 = u * u;
    REAL8 u3 = u2 * u;
    REAL8 u4 = u2 * u2;
     
    REAL8 rc, drcBydr, d2rcBydr2;
    XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadius(&rc, &drcBydr, &d2rcBydr2, r, LALParams);

    /* A potential and derivative with respect to u */
    REAL8 Aorb, dAorbBydu, d2AorbBydu2;
    XLALSimInspiralEOBPostAdiabaticMetricA5PNlog(&Aorb, &dAorbBydu, &d2AorbBydu2, rc, LALParams);
    
    /* Add here tides if needed */
    if (useTidal)
    {
    	// Include tides code if needed
        ;
    }
    
    /* A potential and derivative with respect to r */  
    REAL8 uc  = 1. / rc;
    REAL8 uc2 = uc * uc;
    REAL8 uc3 = uc2 * uc;
    REAL8 uc4 = uc2 * uc2;
    
    REAL8 dAorb  = -dAorbBydu * uc2;
    REAL8 d2Aorb = 2.*dAorbBydu*uc3 + d2AorbBydu2*uc4;
    
    /* Correct A for spin */
    REAL8 AKerr_Multipole = (1.+2.*uc) / (1.+2.*u);
    REAL8 fss = 1.;
         
    *A = Aorb * AKerr_Multipole * fss;    
    *dA = dAorb*drcBydr*(1.+2.*uc)/(1.+2.*u) - 2.*Aorb*drcBydr*uc2/(1.+2.*u) + 2.*Aorb*(1.+2.*uc)*u2/((1.+2.*u)*(1.+2.*u));
    *d2A = d2Aorb*(1.+2.*uc)/(1.+2.*u) + 4.*dAorb*( u2*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) - uc2/(1.+2.*u)*drcBydr) + Aorb*(-4.*u3*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) + 8.*u4*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)*(1.+2.*u))+4.*uc3*(1.+2.*u)*drcBydr*drcBydr - 2.*uc2/(1.+2.*u)*d2rcBydr2);
    
    // /* D potential and derivative with respect to r */
    REAL8 Dp = 1.0 + 6.*nu*uc2 - 2.*(3.0*nu-26.0)*nu*uc3;
    REAL8 D  = 1./Dp;
    REAL8 dD = 6.*uc2*(2.*nu*uc-(3.*nu-26.)*nu*uc2)*D*D;
    
    // /* B potential and derivative with respect to r */
    *B = r * r * uc2 * D / (*A);
    *dB = (dD*(*A) - D*(*dA))/((*A)*(*A));

    return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadius(
    REAL8 *rc,
    REAL8 *drcBydr,
    REAL8 *d2rcBydr2,
    REAL8 r,
    LALDict *LALParams)
{
    const char *centrifugalRadius = XLALDictLookupStringValue(LALParams, "centrifugalRadius");
    // const char centrifugalRadius[] = "LO";

    if (strcmp(centrifugalRadius, "LO") == 0)
    {
        XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadiusLO(rc, drcBydr, d2rcBydr2, r, LALParams);
    }
    else if (strcmp(centrifugalRadius, "NLO") == 0)
    {
        XLALPrintError("The option NLO is not available yet.\n");
    }
    else if (strcmp(centrifugalRadius, "NNLO") == 0)
    {
        XLALPrintError("The option NNLO is not available yet.\n");
    }
    else if (strcmp(centrifugalRadius, "NNLOS4") == 0)
    {
        XLALPrintError("The option NNLOS4 is not available yet.\n");
    }
    else if (strcmp(centrifugalRadius, "NOSPIN") == 0)
    {
        XLALPrintError("The option NOSPIN is not available yet.\n");
    }
    else if (strcmp(centrifugalRadius, "NOTIDES") == 0)
    {
        XLALPrintError("The option NOTIDES is not available yet.\n");
    }
    else
    {
        XLALPrintError("Wrong option for centrifugalRadius parameter, use one of \"LO\", \"NLO\", \"NNLO\", \"NNLOS4\", \"NOSPIN\", and \"NOTIDES\".\n");
    }

    return XLAL_SUCCESS;
}

double
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
    REAL8 prstar_sol,
    void *params)
{
    struct PostAdiabaticRootSolveParams *prstarParams = (struct PostAdiabaticRootSolveParams *)params;

    const REAL8 nu = XLALDictLookupREAL8Value(prstarParams->LALParams, "nu");
    const REAL8 z3 = XLALDictLookupREAL8Value(prstarParams->LALParams, "z3");
    const REAL8 S = XLALDictLookupREAL8Value(prstarParams->LALParams, "S");
    const REAL8 Sstar = XLALDictLookupREAL8Value(prstarParams->LALParams, "Sstar");

    REAL8 r = prstarParams->r;
    REAL8 rc = prstarParams->rc;
    REAL8 drc_dr = prstarParams->drcBydr;
    REAL8 uc2 = prstarParams->uc2;
    REAL8 duc_dr = prstarParams->ducBydr;
    REAL8 pphi = prstarParams->pphi;
    REAL8 dpphi_dr = prstarParams->dpphiBydr;
    REAL8 A = prstarParams->A;
    REAL8 B = prstarParams->B;
    REAL8 dA = prstarParams->dA;
    LALDict *LALParams = prstarParams->LALParams;

    REAL8 ggm[14];

    UNUSED REAL8 G;
    UNUSED REAL8 dG_dr;
    UNUSED REAL8 dG_dprstar;
    UNUSED REAL8 dG_dprstarbyprstar;

    REAL8 H;
    REAL8 Heff;
    REAL8 Heff_orb;
    REAL8 dHeff_dr;
    REAL8 dHeff_dprstar;
    REAL8 dHeff_dpphi;
    REAL8 d2Heff_dprstar20;

    REAL8 E;
    REAL8 Omg;

    REAL8 Heff_orb_f;
    REAL8 Heff_f;
    REAL8 E_f;

    REAL8 psi;
    REAL8 r_omg2;
    REAL8 v_phi;
    REAL8 x;
    REAL8 jhat;
    REAL8 flux;

    REAL8 prefactor;
    REAL8 main_sol;

    REAL8 total;

    XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rc, drc_dr, prstar_sol, LALParams);

    G = ggm[2]*S + ggm[3]*Sstar;
    dG_dr = ggm[6]*S + ggm[7]*Sstar;
    dG_dprstar = ggm[4]*S + ggm[5]*Sstar;
    dG_dprstarbyprstar = ggm[10]*S + ggm[11]*Sstar;

    XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, r, rc, drc_dr, pphi, prstar_sol, S, Sstar, A, dA, LALParams);

    E = nu * H;
    Omg = dHeff_dpphi / E;

    Heff_orb_f = sqrt(A * (1.+pphi*pphi*uc2));
    Heff_f = G*pphi + Heff_orb_f;
    E_f = sqrt(1 + 2*nu*(Heff_f-1));

    psi = (duc_dr + dG_dr*rc*sqrt(A/(pphi*pphi) + A*uc2)/A)/(-0.5*dA);
    r_omg2 = pow(((1./sqrt(rc*rc*rc*psi))+G)/(E_f), -2./3.);
    v_phi = r_omg2 * Omg;
    x = v_phi * v_phi;
    jhat = pphi / (r_omg2*v_phi);
    XLALSimInspiralEOBPostAdiabaticFluxS(&flux, x, Omg, r_omg2, E, Heff, jhat, r, prstar_sol, 0.0, LALParams);

    prefactor = sqrt(A/B) * 1. / (nu*H*Heff_orb);
    main_sol = prstar_sol*(1+2*z3*A/(rc*rc)*prstar_sol*prstar_sol) + Heff_orb*pphi*dG_dprstar;

    total = prefactor * main_sol;

    return dpphi_dr*total - flux;
}

double
XLALSimInspiralEOBPostAdiabaticdprstarFunc(
    REAL8 pphi_sol,
    void *params)
{
    struct PostAdiabaticRootSolveParams *pphiParams = (struct PostAdiabaticRootSolveParams *)params;

    const REAL8 z3 = XLALDictLookupREAL8Value(pphiParams->LALParams, "z3");
    const REAL8 S = XLALDictLookupREAL8Value(pphiParams->LALParams, "S");
    const REAL8 Sstar = XLALDictLookupREAL8Value(pphiParams->LALParams, "Sstar");

    REAL8 r = pphiParams->r;
    REAL8 rc = pphiParams->rc;
    REAL8 drcBydr = pphiParams->drcBydr;
    REAL8 dAuc2Bydr = pphiParams->dAuc2Bydr;
    REAL8 dprstarBydr = pphiParams->dprstarBydr;
    REAL8 dA = pphiParams->dA;
    REAL8 prstar = pphiParams->prstar;
    REAL8 A = pphiParams->A;
    REAL8 uc2 = pphiParams->uc2;
    LALDict *LALParams = pphiParams->LALParams;

    REAL8 ggm[14];

    UNUSED REAL8 G;
    UNUSED REAL8 dGBydr;
    UNUSED REAL8 dGBydprstar;
    UNUSED REAL8 dG_dprstarbyprstar;

    XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rc, drcBydr, prstar, LALParams);

    G = ggm[2]*S + ggm[3]*Sstar;
    dGBydr = ggm[6]*S + ggm[7]*Sstar;
    dGBydprstar = ggm[4]*S + ggm[5]*Sstar;

    REAL8 H;
    REAL8 Heff;
    REAL8 HeffOrb;
    REAL8 dHeff_dr;
    REAL8 dHeff_dprstar;
    REAL8 dHeff_dpphi;
    REAL8 d2Heff_dprstar20;

    XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &Heff, &HeffOrb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, r, rc, drcBydr, pphi_sol, prstar, S, Sstar, A, dA, LALParams);

    REAL8 aCoeff;
    REAL8 bCoeff;
    REAL8 cCoeff;

    REAL8 result;

    aCoeff = dAuc2Bydr;
    bCoeff = 2 * HeffOrb * (dGBydr+dGBydprstar*dprstarBydr);
    cCoeff = dA + 2*prstar*dprstarBydr*(1+2*z3*A*uc2*prstar*prstar) + z3*dAuc2Bydr*gsl_pow_int(prstar, 4);

    result = aCoeff*pphi_sol*pphi_sol + bCoeff*pphi_sol + cCoeff;

    return result;
}