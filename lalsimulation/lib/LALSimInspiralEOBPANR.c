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