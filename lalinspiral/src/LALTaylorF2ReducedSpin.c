/*
*  Copyright (C) 2011 P. Ajith
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <lal/LALInspiral.h>

/**
 * Generate the "reduced-spin templates" proposed in http://arxiv.org/abs/1107.1267
 */
int XLALTaylorF2ReducedSpin(REAL4Vector *signalvec, 
		InspiralTemplate *params) {

    REAL8 df, shft, phi0, amp0, amp, f, m, eta, delta, chi_s, chi_a, chi, Psi;
    REAL8 psiNewt, psi2, psi3, psi4, psi5, psi6, psi6L, psi7, psi3S, psi4S, psi5S;
    REAL8 alpha2, alpha3, alpha4, alpha5, alpha6, alpha6L, alpha7, alpha3S, alpha4S, alpha5S; 
    REAL8 v, v2, v3, v4, v5, v6, v7, v0, mSevenBySix, piM, oneByThree; 
    INT4 i, j, n, nBy2;

    /* check inputs */
    if (!signalvec || !(signalvec->data) || !params) {
        XLALPrintError(LALINSPIRALH_MSGENULL);
        XLAL_ERROR(XLAL_EFAULT);
    }
    if ((signalvec->length < 2) || (params->fCutoff <= params->fLower)  
            || (params->mass1 <= 0) || (params->mass2 <= 0)
            || (params->spin1[0] != 0) || (params->spin1[1] != 0) 
            || (params->spin2[0] != 0) || (params->spin2[1] != 0)) {
        XLALPrintError(LALINSPIRALH_MSGECHOICE);
        XLAL_ERROR(XLAL_EDOM);
    }

	/* fill the waveform with zeros */
	memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));

    /* compute total mass (secs), mass ratio and the reduced spin parameter */
    m = (params->mass1+params->mass2)*LAL_MTSUN_SI;
    eta = params->mass1*params->mass2/pow(params->mass1+params->mass2,2.);
    delta = (params->mass1-params->mass2)/(params->mass1+params->mass2);
    chi_s = (params->spin1[2] + params->spin2[2])/2.;
    chi_a = (params->spin1[2] - params->spin2[2])/2.;
    chi = chi_s*(1. - 76.*eta/113.) + delta*chi_a;

    /* freq resolution and the low-freq bin */
    df = params->tSampling/signalvec->length;
    n = signalvec->length;

    /* extrinsic parameters */
    phi0  = params->startPhase;
    amp0 = pow(m,5./6.)*sqrt(5.*eta/24.)/(pow(LAL_PI,2./3.)*params->distance/LAL_C_SI);
    shft = -2.*LAL_PI * ((REAL8)signalvec->length/params->tSampling + 
            params->nStartPad/params->tSampling + params->startTime);
    // UNUSED!!: REAL8 t0 = params->startTime;
    phi0 = params->startPhase;

    /* spin terms in the amplitude and phase (in terms of the reduced
     * spin parameter */
    psi3S = 113.*chi/3.;
    psi4S = 63845.*(-81. + 4.*eta)*chi*chi/(8.*pow(-113. + 76.*eta, 2.)); 
    psi5S = -565.*(-146597. + 135856.*eta + 17136.*eta*eta)*chi/(2268.*(-113. + 76.*eta)); 

    alpha3S = (113.*chi)/24.; 
    alpha4S = (12769.*pow(chi,2)*(-81. + 4.*eta))/(32.*pow(-113. + 76.*eta,2)); 
    alpha5S = (-113.*chi*(502429. - 591368.*eta + 1680*eta*eta))/(16128.*(-113 + 76*eta)); 

    /* coefficients of the phase at PN orders from 0 to 3.5PN */
    psiNewt = 3./(128.*eta);
    psi2 = 3715./756. + 55.*eta/9.;
    psi3 = psi3S - 16.*LAL_PI;
    psi4 = 15293365./508032. + 27145.*eta/504. + 3085.*eta*eta/72. + psi4S;
    psi5 = (38645.*LAL_PI/756. - 65.*LAL_PI*eta/9. + psi5S);
    psi6 = 11583231236531./4694215680. - (640.*LAL_PI*LAL_PI)/3. - (6848.*LAL_GAMMA)/21. 
             + (-5162.983708047263 + 2255.*LAL_PI*LAL_PI/12.)*eta 
             + (76055.*eta*eta)/1728. - (127825.*eta*eta*eta)/1296.;
    psi6L = -6848./21.;
    psi7 = (77096675.*LAL_PI)/254016. + (378515.*LAL_PI*eta)/1512.  
             - (74045.*LAL_PI*eta*eta)/756.;

    /* amplitude coefficients */
    alpha2 = 1.1056547619047619 + (11*eta)/8.; 
    alpha3 = -2*LAL_PI + alpha3S; 
    alpha4 = 0.8939214212884228 + (18913*eta)/16128. + (1379*eta*eta)/1152. + alpha4S; 
    alpha5 = (-4757*LAL_PI)/1344. + (57*eta*LAL_PI)/16. + alpha5S; 
    alpha6 = -58.601030974347324 + (3526813753*eta)/2.7869184e7 - 
                (1041557*eta*eta)/258048. + (67999*eta*eta*eta)/82944. + 
                (10*pow(LAL_PI,2))/3. - (451*eta*pow(LAL_PI,2))/96.; 
    alpha6L = 856/105.; 
    alpha7 = (-5111593*LAL_PI)/2.709504e6 - (72221*eta*LAL_PI)/24192. - 
                (1349*eta*eta*LAL_PI)/24192.; 

    /* select the terms according to the PN order chosen */
    switch (params->order) {
        case LAL_PNORDER_THREE:
            psi7 = 0.;
            alpha7 = 0.;
            break;
        case LAL_PNORDER_TWO_POINT_FIVE:
            psi6 = 0.;
            psi6L = 0.;
            psi7 = 0.;
            alpha6 = 0.;
            alpha6L = 0.;
            alpha7 = 0.;
            break;
        case LAL_PNORDER_TWO:
            psi5 = 0.;
            psi6 = 0.;
            psi6L = 0.;
            psi7 = 0.;
            alpha5 = 0.;
            alpha6 = 0.;
            alpha6L = 0.;
            alpha7 = 0.;
            break;
        case LAL_PNORDER_ONE_POINT_FIVE:
            psi4 = 0.;
            psi5 = 0.;
            psi6 = 0.;
            psi6L = 0.;
            psi7 = 0.;
            alpha4 = 0.;
            alpha5 = 0.;
            alpha6 = 0.;
            alpha6L = 0.;
            alpha7 = 0.;
            break;
        case LAL_PNORDER_ONE:
            psi3 = 0.;
            psi4 = 0.;
            psi5 = 0.;
            psi6 = 0.;
            psi6L = 0.;
            psi7 = 0.;
            alpha3 = 0.;
            alpha4 = 0.;
            alpha5 = 0.;
            alpha6 = 0.;
            alpha6L = 0.;
            alpha7 = 0.;
            break;
        case LAL_PNORDER_NEWTONIAN:
            psi2 = 0.;
            psi3 = 0.;
            psi4 = 0.;
            psi5 = 0.;
            psi6 = 0.;
            psi6L = 0.;
            psi7 = 0.;
            alpha2 = 0.;
            alpha3 = 0.;
            alpha4 = 0.;
            alpha5 = 0.;
            alpha6 = 0.;
            alpha6L = 0.;
            alpha7 = 0.;
            break;
        default:
            break;
    }
    
    /* fill the zero and Nyquist */
    signalvec->data[0] = 0.;
    signalvec->data[n/2] = 0.;

    mSevenBySix = -7./6.;
    piM = LAL_PI*m;
    oneByThree = 1./3.;
    nBy2 = n/2;

    v0 = pow(LAL_PI*m*params->fLower, 1./3.);

    for (i=1; i<nBy2; i++) {

        /* this is the index of the imaginary part */
        j = n-i;

        /* fourier frequency corresponding to this bin */
      	f = i * df;
    
        /* PN expansion parameter */
        v = pow(piM*f, oneByThree);

        v2 = v*v;   v3 = v2*v;  v4 = v3*v;  v5 = v4*v;  v6 = v5*v;  v7 = v6*v;

        if ((f < params->fLower) || (f > params->fCutoff)) {
            amp = 0.;
            Psi = 0.;
        }
        else {

            /* compute the phase and amplitude */
            Psi = psiNewt*pow(v, -5.)*(1. 
                    + psi2*v2 + psi3*v3 + psi4*v4 
                    + psi5*v5*(1.+3.*log(v/v0)) 
                    + (psi6 + psi6L*log(4.*v))*v6 + psi7*v7); 

            amp = amp0*pow(f, mSevenBySix)*(1. 
                    + alpha2*v2 + alpha3*v3 + alpha4*v4 
                    + alpha5*v5 + (alpha6 + alpha6L*(LAL_GAMMA+log(4.*v)) )*v6 
                    + alpha7*v7); 

        }

        signalvec->data[i] = (REAL4) (amp * cos(Psi+shft*f+phi0+LAL_PI/4.));   /* real */
        signalvec->data[j] = (REAL4) -(amp * sin(Psi+shft*f+phi0+LAL_PI/4.));  /* imag */

    }    

    params->fFinal = params->fCutoff;

    return XLAL_SUCCESS;
}

/**
 * Generate two orthogonal "reduced-spin" templates
 */
int XLALTaylorF2ReducedSpinTemplates(REAL4Vector *signalvec1,
		REAL4Vector *signalvec2,
		InspiralTemplate *params) {

	/* check inputs */
	if (!signalvec1 || !signalvec2 || !(signalvec1->data) || !(signalvec2->data)) {
    	XLALPrintError(LALINSPIRALH_MSGENULL);
    	XLAL_ERROR(LALINSPIRALH_ENULL);
  	}

  	/* generate one waveform with startPhase specified by the user */
  	if (!XLALTaylorF2ReducedSpin(signalvec1, params))
    	XLAL_ERROR(XLAL_EFUNC);

  	/* generate another waveform orthogonal to it */
  	params->startPhase += LAL_PI_2;
  	if (!XLALTaylorF2ReducedSpin(signalvec2, params))
    	XLAL_ERROR(XLAL_EFUNC);

    return XLAL_SUCCESS;
}

/**
 * Compute the chirp time of the "reduced-spin" templates
 */
REAL8 XLALChirpTimeReducedSpin(REAL8 v, REAL8 m1, REAL8 m2, REAL8 spin1, 
        REAL8 spin2, UINT4 pnOrder) {

    REAL8 chis, chia, chis2, chia2;
    REAL8 tau, tk[8], eta2, eta3;
    UINT4 k;

    REAL8 m = m1+m2; 
    REAL8 eta = m1*m2/(m*m);
    REAL8 delta = (m1-m2)/m;

    chis  = (spin1+spin2)/2.;
    chia  = (spin1-spin2)/2.;

    eta2 = eta*eta;
    eta3 = eta2*eta;
    chis2 = chis*chis;
    chia2 = chia*chia;

    /* chirp time coefficients up to 3.5PN  */
    tk[0] = (5.*m*LAL_MTSUN_SI)/(256.*pow(v,8)*eta); 
    tk[1] = 0.;
    tk[2] = 2.9484126984126986 + (11*eta)/3.;
    tk[3] = (-32*LAL_PI)/5. + (226*chia*delta)/15. + chis*(15.066666666666666 - (152*eta)/15.);
    tk[4] = 6.020630590199042 + ((233*chis*chia)/24. - (719*chia*chis)/24.)*delta + 
       chia*chia*(4.854166666666667 - 20*eta) + chis2*(-14.979166666666666 - eta/12.) + 
       chis*chis*(4.854166666666667 + (7*eta)/12.) + (5429*eta)/504. + (617*eta2)/72. + 
       chia2*(-14.979166666666666 + 60*eta);
    tk[5] = (-7729*LAL_PI)/252. + (13*LAL_PI*eta)/3. + delta*((146597*chia)/756. + (28*chia*eta)/3.) + 
       chis*(193.91137566137567 - (4852*eta)/27. - (68*eta2)/3.);
    tk[6] = -428.291776175525 + (128*LAL_PI*LAL_PI)/3. + (6848*LAL_GAMMA)/105. + (3147553127*eta)/3.048192e6 - 
       (451*LAL_PI*LAL_PI*eta)/12. - (15211*eta2)/1728. + (25565*eta3)/1296. + (6848*log(4*v))/105.;  
    tk[7] = (-15419335*LAL_PI)/127008. - (75703*LAL_PI*eta)/756. + (14809*LAL_PI*eta2)/378.;

    /* compute chirp time. return it */
    tau = 1.;
    for (k = 2; k<=pnOrder; k++) {
        tau = tau + tk[k]*pow(v, k);
    }
    tau = tau*tk[0];

    return (tau);
}



