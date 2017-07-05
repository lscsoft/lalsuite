#ifndef _LALSIMIMRSPINPRECEOBWFGEN_C
#define _LALSIMIMRSPINPRECEOBWFGEN_C

/**************************************
 **************************************
 ** OPTV3:
 ** This function outputs hp2,  hp1,  h0,  hm1, and hm2
 ** which are the h modes for l=2; m = 2, 1, ..., -2
 ** It is Organized into 7 steps:
 **       (1) the Hamiltonian
 **       (2) seobCoeffs
 **       (3) dvalues
 **       (4) omega
 **       (5) hCoeffs
 **       (6) polarDynamics
 **       (7) hLM
 **************************************
 **************************************/
static void XLALEOBSpinPrecCalcHModes(
                  COMPLEX16 * hp2, /**<< OUTPUT: h; m = 2 mode */
                  COMPLEX16 * hp1, /**<< OUTPUT: h; m = 1 mode */
                  COMPLEX16 * h0,  /**<< OUTPUT: h; m = 0 mode */
                  COMPLEX16 * hm1, /**<< OUTPUT: h; m = -1 mode */
                  COMPLEX16 * hm2, /**<< OUTPUT: h; m = -2 mode */
                  SpinEOBParams * seobParams, /**<< EOB Parameters */
                  REAL8Vector * values /**<< Dynamical variables */
                  ){
    COMPLEX16 hLM, hNQC;
    REAL8 ham, omega, v, magR2, s1dotLN, s2dotLN, chiS, chiA, tplspin;
    REAL8 rcrossrdot[3], rcrossp[3];
    REAL8 a2 = 0.0;
    REAL8 m1 = seobParams->eobParams->m1;
    REAL8 m2 = seobParams->eobParams->m2;

    REAL8 dvaluesData[14] = {0.};
    REAL8 polarDynamicsData[4] = {0.};
    REAL8Vector dvalues;
    REAL8Vector polarDynamics;
    dvalues.length = 14;
    polarDynamics.length = 4;
    polarDynamics.data = polarDynamicsData;
    dvalues.data = dvaluesData;

    /**  OPTV3: (1): the Hamiltonian **/
    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3099-3102*/
    for(int i=0; i<3; i++){
        seobParams->sigmaKerr->data[i] = values->data[i+6] + values->data[i+9];
        seobParams->sigmaStar->data[i] = (m2/m1)*values->data[i+6] + (m1/m2)*values->data[i+9];
        a2 += seobParams->sigmaKerr->data[i]*seobParams->sigmaKerr->data[i];
    }

    /* Calculate the value of the Hamiltonian */
    {

      REAL8 rCartData[3], pCartData[3], s1CartData[3], s2CartData[3];

      REAL8Vector rCartVec;
      REAL8Vector pCartVec;
      REAL8Vector s1CartVec;
      REAL8Vector s2CartVec;

      rCartVec.length = pCartVec.length = s1CartVec.length = s2CartVec.length = 3;

      rCartVec.data = rCartData;
      pCartVec.data = pCartData;
      s1CartVec.data = s1CartData;
      s2CartVec.data = s2CartData;

      memcpy(rCartVec.data,    values->data,     3*sizeof(REAL8));
      memcpy(pCartVec.data,    values->data+3,   3*sizeof(REAL8));
      memcpy(s1CartVec.data,   values->data+6,   3*sizeof(REAL8));
      memcpy(s2CartVec.data,   values->data+9,   3*sizeof(REAL8));

      REAL8Vector * rCart  = &rCartVec;
      REAL8Vector * pCart  = &pCartVec;
      REAL8Vector * s1Cart = &s1CartVec;
      REAL8Vector * s2Cart = &s2CartVec;

      ham = XLALSimIMRSpinPrecEOBHamiltonian( seobParams->eobParams->eta, rCart, pCart, s1Cart, s2Cart, seobParams->sigmaKerr, seobParams->sigmaStar, seobParams->tortoise, seobParams->seobCoeffs );
    }
    /** OPTV3: (2): hCoeffs**/
    seobParams->a = sqrt(a2); /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3102*/

    // OPTV3: SpinAlignedEOBversion = 2
    XLALSimIMRCalculateSpinEOBHCoeffs( seobParams->seobCoeffs, seobParams->eobParams->eta, seobParams->a, 2); /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3127*/

    /** OPTV3: (3): dValues**/
    memset( dvalues.data, 0, 14*sizeof(REAL8));
    XLALSpinPrecHcapRvecDerivative_exact(0.0, values->data, dvalues.data, seobParams); /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3006*/

    /** OPTV3: (4): omega **/

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3017-3020*/
    /* Calculate omega */
    cross_product(values->data,dvalues.data,rcrossrdot);
    magR2 = inner_product(values->data, values->data);
    omega = sqrt(inner_product(rcrossrdot,rcrossrdot))/magR2;
    v = cbrt( omega );

    /** OPTV3: (5): hCoeffs **/

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3107-3111*/
    /* Update Hamiltonian coefficients as per |Skerr| */
    s1dotLN = inner_product(values->data+6, rcrossrdot)/(m1*m2);
    s2dotLN = inner_product(values->data+9, rcrossrdot)/(m1*m2);
    chiS = 0.5 * (s1dotLN + s2dotLN);
    chiA = 0.5 * (s1dotLN - s2dotLN);

    // OPTV3: Because SpinAlignedEOBVersion = 2
    tplspin = (1.-2.*seobParams->eobParams->eta) * chiS + (m1 - m2)/(m1 + m2) * chiA; /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3119*/

    /* Update hlm coefficients */
    // OPTV3: SpinAlignedEOBversion = 2
    XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( seobParams->eobParams->hCoeffs, m1, m2, seobParams->eobParams->eta, tplspin, chiS, chiA, 3 ); /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3133*/

    /** OPTV3: (6): polarDynamics **/
    /* Calculate the polar data */
    /* Calculate the orbital angular momentum */
    cross_product( values->data, values->data+3, rcrossp );

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3150-3153*/
    polarDynamics.data[0] = sqrt(magR2);
    polarDynamics.data[1] = values->data[13] + values->data[12]; // OPTV3: phiMod.data[i] + phiDMod.data[i];
    polarDynamics.data[2] = inner_product(values->data, values->data+3) / polarDynamics.data[0];
    polarDynamics.data[3] = sqrt(inner_product(rcrossp, rcrossp));

    /** OPTV3: (7) hLM **/

    XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform_exact( &hLM, &polarDynamics, values, v, ham, 2, 2, seobParams ); /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3155*/
    XLALSimIMRSpinEOBNonQCCorrection( &hNQC, values, omega, seobParams->nqcCoeffs ); /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3160*/

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3165-3167*/
    hLM *= hNQC;
    *hp2= hLM;
    *hm2= conjl(hLM);

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3178*/
    XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform_exact( &hLM, &polarDynamics, values, v, ham, 2, 1, seobParams );

    /*OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3183-3185*/
    *hp1  = hLM;
    *hm1 = conjl(hLM);
    *h0  = 0.0;
}

static void XLALEOBSpinPrecCalcHplusHcross(
                  SpinEOBParams * seobParams, /**<< EOB Parameters */
                  REAL8Vector * values, /**<< Dynamical variables */
                  REAL8 aI2P, /**<< alpha Euler angle from inertial to precessing frame */
                  REAL8 bI2P, /**<< beta Euler angle from inertial to precessing frame */
                  REAL8 gI2P, /**<< gamma Euler angle from inertial to precessing frame */
                  REAL8 alJtoI, /**<< alpha Euler angle from final-J frame to Inertial frame*/
                  REAL8 betJtoI, /**<< beta Euler angle from final-J frame to Inertial frame*/
                  REAL8 gamJtoI, /**<< gamma Euler angle from final-J frame to Inertial frame*/
                  REAL8 JframeEx[3], /**<< x-axis of the total-angular-momentum frame */
                  REAL8 JframeEy[3], /**<< y-axis of the total-angular-momentum frame */
                  REAL8 JframeEz[3], /**<< z-axis of the total-angular-momentum frame */
                  COMPLEX16 Y[5], /**<< Spherical Harmonics */
                  REAL8 * hplus, /**<< OUTPUT: h_+ */
                  REAL8 * hcross /**<< OUTPUT: h_x */
                  ){
    COMPLEX16 hTSp2, hTSp1, hTS0, hTSm1, hTSm2;
    REAL8 aP2J, bP2J, gP2J;
    REAL8 LNhx, LNhy, LNhz;
    COMPLEX16 hxx[5], hresult[5], hresult2[5];
    REAL8 rcrossrdot[3], magrcrossrdot;
    REAL8 rdot[3] = {0.,0.,0.};

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3311*/
    XLALEOBSpinPrecCalcHModes(&hTSp2, &hTSp1, &hTS0, &hTSm1, &hTSm2, seobParams,values);
    /** OPTV3: rotate by alpha, beta, gamma */

    /* First prepare the rotation angles */
    hxx[0]= hTSm2; hxx[1]=hTSm1; hxx[2]=hTS0; hxx[3]=hTSp1, hxx[4]=hTSp2;

    /* Calculate dr/dt */
    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3006, 3311*/
    XLALSpinPrecHcapRvecDerivative_exact(0.0, values->data, rdot, seobParams);

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3030, 3323*/
    cross_product(values->data,rdot,rcrossrdot);

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3033*/
    magrcrossrdot = sqrt(inner_product(rcrossrdot,rcrossrdot));

    LNhx = rcrossrdot[0]/magrcrossrdot;
    LNhy = rcrossrdot[1]/magrcrossrdot;
    LNhz = rcrossrdot[2]/magrcrossrdot;

        /*Rounding LNhz creates the possible issue of causing the input for the acos() call in EulerAnglesP2J
        to be out of the [-1,1] interval required by the function, so it has been commented out.
        The other two rounding statements are being removed because they do not appear in the current version
        of the unoptimized code. -Tom Adams*/
//        if (fabs(LNhx) < 1.e-7)
//            LNhx = 0.0;
//        if (fabs(LNhy) < 1.e-7)
//            LNhy = 0.0;
//        if (fabs(LNhz-1.0) < 1.e-7)
//            LNhz = 1.0;

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3342-3358, 3375-3377*/
    EulerAnglesP2J(&aP2J,&bP2J,&gP2J,aI2P, bI2P, gI2P, LNhx, LNhy, LNhz, JframeEx, JframeEy, JframeEz );
    aP2J = -aP2J;
    gP2J = -gP2J;

    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 3570,3602*/
    /* OPTV3: From XLALSimInspiralPrecessionRotateModes; LALSimInspiralPrecess.c; line(s): 62-66*/
    COMPLEX16 W[5][5];
    int l = 2;
    for(int m=0; m<2*l+1; m++){
     hresult[m]=0.+0.*I;
        for(int mp=0; mp<2*l+1; mp++){
            W[m][mp]=XLALWignerDMatrix( l, mp-l, m-l, aP2J , bP2J , gP2J);
            hresult[m] += hxx[mp] * W[m][mp];
//            hresult[m] += hxx[mp] * XLALWignerDMatrix( l, mp-l, m-l, aP2J , bP2J , gP2J);
//            printf("[%.4e|%.4e]= hresult[%d] += hxx[%d] * W[%d][%d] = [%.4e|%.4e] * [%.4e|%.4e]\n",cabs(hresult[m]),carg(hresult[m]),m,mp,m,mp, cabs(hxx[mp]), carg(hxx[mp]), cabs(W[m][mp]),carg(W[m][mp]));
        }
    }
    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; 3966*/
    /* OPTV3: From XLALSimInspiralPrecessionRotateModes; LALSimInspiralPrecess.c; 62-66*/
    for(int m=0; m<2*l+1; m++){
        hresult2[m]=0.+0.*I;
        for(int mp=0; mp<2*l+1; mp++)
            hresult2[m] += hresult[mp] * XLALWignerDMatrix( l, mp-l, m-l, alJtoI , betJtoI , gamJtoI);
    }
    /* OPTV3: From XLALSimIMRSpinEOBWaveformAll; LALSimIMRSpinPrecEOB.c; line(s): 4022-4028*/
    COMPLEX16 x11 = 0.0+0.0*I;
    for (int i=0; i<2*l+1; i++) x11 += Y[i]*hresult2[i];
    /** OPTV3: Now, evaluate the real & imaginary parts to get hplus & hcross */
    *hplus  = creal(x11);
    *hcross = cimag(x11);
}

static void XLALEOBSpinPrecGenerateHTSModesFromEOMSoln(
                              COMPLEX16 * hTSm2, /**<< OUTPUT: Complex h array, mode m=2, length=retLen */
                              COMPLEX16 * hTSm1, /**<< OUTPUT: Complex h array, mode m=1, length=retLen */
                              COMPLEX16 * hTS0, /**<< OUTPUT: Complex h array, mode m=0, length=retLen */
                              COMPLEX16 * hTSp1, /**<< OUTPUT: Complex h array, mode m=-1, length=retLen */
                              COMPLEX16 * hTSp2, /**<< OUTPUT: Complex h array, mode m=-2, length=retLen */
                              int retLen, /**<< returned length from ODE solution */
                              REAL8Array * dynamics, /**<< Dynamical Variables */
                              SpinEOBParams * seobParams /**<< EOB Parameters */
                              ){

    REAL8 valuesdata[14] = {0.};
    REAL8Vector values;
    values.length=14;
    values.data= valuesdata;

    for(int t = 0; t<retLen; t++){
        for (int i=1; i<=14; i++)
            values.data[i-1]=dynamics->data[t+i*retLen];
        XLALEOBSpinPrecCalcHModes(hTSp2+t,hTSp1+t,hTS0+t,hTSm1+t,hTSm2+t,seobParams,&values);
    }
}

static void XLALEOBSpinPrecGenerateHplusHcrossTSFromEOMSoln(
                              REAL8Vector * h,              /**<< OUTPUT: Vector containing time, hplus, and hcross */
                              int retLen,                   /**<< Length of ODE solution */
                              REAL8Vector * AlphaI2PVec,    /**<< Vector of Euler angle alpha, from Inertial to precessing frame, across time */
                              REAL8Vector * BetaI2PVec,     /**<< Vector of Euler angle beta, from Inertial to precessing frame, across time */
                              REAL8Vector * GammaI2PVec,    /**<< Vector of Euler angle gamma, from Inertial to precessing frame, across time */
                              REAL8 Jx[3],                  /**<< x-axis of the total-angular-momentum frame */
                              REAL8 Jy[3],                  /**<< y-axis of the total-angular-momentum frame */
                              REAL8 Jz[3],                  /**<< z-axis of the total-angular-momentum frame */
                              COMPLEX16 Y[5],               /**<< Spherical Harmonics */
                              REAL8Array * dynamics,        /**<< Dynamical Variables */
                              SpinEOBParams * seobParams    /**<< EOB Parameters */
                              ){
    REAL8  alJtoI = -atan2(Jy[2], -Jx[2]);
    REAL8  betJtoI = acos(Jz[2]);
    REAL8  gamJtoI = -atan2(Jz[1], -Jz[0]);

    REAL8 valuesdata[14] = {0.};
    REAL8Vector values;
    values.length=14;
    values.data= valuesdata;
    for(int t = 0; t<retLen; t++){
        for (int i=1; i<=14; i++)
            values.data[i-1]=dynamics->data[t+i*retLen];
        h->data[t]=dynamics->data[t];
        XLALEOBSpinPrecCalcHplusHcross(seobParams,&values,AlphaI2PVec->data[t],BetaI2PVec->data[t],GammaI2PVec->data[t], alJtoI,betJtoI,gamJtoI,Jx,Jy,Jz,Y,h->data+t+retLen,h->data+t+2*retLen);
    }
}
#endif // _LALSIMIMRSPINPRECEOBWFGEN_C
