/*
 * Copyright (C) 2013 J. Lucaccioni, L. Wade
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */
/**
 * @addtogroup LALSimNeutronStarEOS_c
 * @{
 */
/**
 * @name Creation routines for dynamic polytrope equations of state
 * @{
 */

#include <lal/LALSimReadData.h>
#include <gsl/gsl_interp.h>
#include <lal/LALSimNeutronStar.h>

/** @cond */

/* Calculates the energy density given the pressure using an n-piece dynamic polytrope. */
static double non_causal_Nparams(double p, double esec[], double psec[], double gsec[], size_t nsec){
    double e;
    int k;
    // Find the polytrope section
    for (k = 0; k <= (int)nsec - 2; ++k)
        if (p < psec[k+1])
            break;
    // Eq. (5) of PRD 97, 123019 (2018)
    e = (esec[k] - psec[k] / (gsec[k] - 1)) * pow(p / psec[k], 1 / gsec[k]) + p / (gsec[k] - 1);
    return e;
}

/* Calculates the energy density given the pressure for an n-piece causal dynamic polytrope. */
static double causal_Nparams(double p, double esec[], double psec[], double usec[], double vsec[], size_t nsec){
    double e;
    int k;
    // Find the polytrope section
    for (k = 0; k <= (int)nsec - 2; ++k)
        if (p < psec[k+1])
            break;
    // Eq. (14) of PRD 97, 123019 (2018)
    e = esec[k] + p - psec[k] + (psec[k] * usec[k] / (1 + vsec[k+1])) * (pow(p / psec[k], 1 + vsec[k+1]) - 1);
    return e;
}

/** @endcond */

/**
 * @brief Reads dynamic analytic eos parameters to make an eos.
 * @details Reads an array of eos parameters to construct either the causal 
 * analytic eos or non-causal analytic polytrope generically outlined in PRD 97,
 * 123019 (2018). The models are dynamic because the stitching pressures are
 * also parameter choices.
 * @param[in] parameters[] Array of dynamic analytic eos parameters.
 * [param0, log10(p1), param1, log10(p2), ... , log10(pN), paramN], where
 * pressure is in SI units, and the params are either gammas (adiabatic
 * indexes) for the non-causal polytrope model or the vs (speed function
 * exponents) for the causal analytic model.
 * @param[in] nsec The number of sections (pressure sub-domains) to stitch to
 * crust eos.
 * @param[in] causal Option to use causal version (0=non-causal, 1=causal).
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSDynamicAnalytic(double parameters[], size_t nsec, int causal){
    // Check that causal flag appropriately assigned
    if (causal!=0 && causal!=1)
	XLAL_ERROR_NULL(XLAL_EINVAL,"Did not specify which approach to take, Causal or Non-Causal");

    // Size is the number of polytrope pieces to stitch to low-density EOS
    if (nsec<1)
        XLAL_ERROR_NULL(XLAL_EINVAL,"Number of polytrope pieces should be at least 1\n");

    // Declares the max pressure(pmax).  Units in geom.
    // Chosen to match spectral decomposition and correspond to e0 and p0 for
    // SLY EOS just below 0.5*rho_nuc. See Sec IID of PRD 98 063004 (2018) for
    // more details. Note: Make adjustable in the future.
    double pmax = 9.829054605e-8;
    double p0 = 4.43784199e-13;

    // Pull out and convert pressure parameters
    double psec[nsec];
    psec[0]=p0;
    for (int i=1; i < (int)nsec; i++)
	// Odd indices, SI -> geom (so c=1)
	psec[i] = pow(10,parameters[i*2-1])*LAL_G_C4_SI;

    if (psec[nsec-1]>pmax)
        XLAL_ERROR_NULL(XLAL_EINVAL,"Highest p is set larger than %e, the limit at which EOS is generated\n",pmax);

    // Fill in arrays with SLy data. Units in geom.
    // These values agree with LALSimNeutronStarEOS_SLY.dat
    double SLYdat_p[] =  {2.497300e-31,   1.592353e-30,
          1.015332e-29,   6.474064e-29,   4.128057e-28,
          2.632173e-27,   1.678353e-26,   1.070168e-25,
          6.823712e-25,   4.351004e-24,   1.915235e-23,
          8.060015e-23,   3.231442e-22,   4.345220e-22,
          1.185661e-21,   3.166995e-21,   8.312020e-21,
          2.151541e-20,   5.516008e-20,   7.219725e-20,
          1.345952e-19,   2.502695e-19,   3.411564e-19,
          4.160967e-19,   5.668037e-19,   1.050983e-18,
          1.946632e-18,   3.604079e-18,   4.678197e-18,
          6.363735e-18,   8.659043e-18,   1.177398e-17,
          1.601262e-17,   2.068090e-17,   2.812536e-17,
          3.823860e-17,   4.915329e-17,   6.683492e-17,
          9.088690e-17,   1.198055e-16,   1.235236e-16,
          1.679755e-16,   2.145757e-16,   2.389499e-16,
          2.718344e-16,   3.695792e-16,   4.805438e-16,
          6.228231e-16,   6.448839e-16,   6.519069e-16,
          6.900794e-16,   7.517173e-16,   7.841887e-16,
          8.122810e-16,   8.948228e-16,   1.780309e-15,
          2.831705e-15,   4.352574e-15,   6.442724e-15,
          9.210148e-15,   1.726355e-14,   2.961343e-14,
          4.760077e-14,   7.280619e-14,   1.069739e-13,
          1.786341e-13,   3.178976e-13,   4.169395e-13,
          4.437842e-13,   8.192049e-13,   2.464440e-12,
          5.749505e-12,   1.129311e-11,   1.962992e-11,
	  3.125683e-11,   4.664233e-11,   6.623341e-11,
          9.045725e-11,   1.197311e-10,   1.544499e-10,
          1.949854e-10,   2.417178e-10,   2.949858e-10,
          3.551282e-10,   4.224506e-10,   4.972587e-10,
          5.798254e-10,   6.704065e-10,   7.692419e-10,
          8.765628e-10,   9.925676e-10,   1.117413e-09,
          1.251265e-09,   1.394370e-09,   1.546730e-09,
          1.708591e-09,   1.880037e-09,   2.061232e-09,   2.252177e-09};

    // Ensure p0 is in the crust eos range
    size_t SLYdat_p_size = sizeof(SLYdat_p) / sizeof(*SLYdat_p);
    if (psec[0] > SLYdat_p[SLYdat_p_size - 1])
        XLAL_ERROR_NULL(XLAL_EINVAL,"p0 is set larger than %e, the limit of crust data",SLYdat_p[SLYdat_p_size - 1]);
    if (psec[0] < SLYdat_p[0])
        XLAL_ERROR_NULL(XLAL_EINVAL,"p0 is set smaller than %e, the lowest value of crust data",SLYdat_p[0]);

    double SLYdat_e[] =  {9.768004e-24,   3.088914e-23,
          9.768005e-23,   3.088915e-22,   9.768010e-22,
          3.088918e-21,   9.768031e-21,   3.088931e-20,
          9.768115e-20,   3.088985e-19,   7.759067e-19,
          1.949099e-18,   4.896221e-18,   6.163579e-18,
          1.229692e-17,   2.454210e-17,   4.898236e-17,
          9.773109e-17,   1.950244e-16,   2.455315e-16,
          3.892642e-16,   6.169263e-16,   7.768370e-16,
          9.006424e-16,   1.192376e-15,   1.890371e-15,
          3.094528e-15,   4.906359e-15,   5.964589e-15,
          7.510611e-15,   9.797770e-15,   1.233538e-14,
          1.553357e-14,   1.881882e-14,   2.462712e-14,
          3.100969e-14,   3.743924e-14,   4.918149e-14,
          6.194543e-14,   7.622245e-14,   8.111554e-14,
          1.052007e-13,   1.264362e-13,   1.370769e-13,
          1.557999e-13,   2.029004e-13,   2.471392e-13,
          3.112816e-13,   3.195299e-13,   3.221410e-13,
          3.362136e-13,   3.585373e-13,   3.701139e-13,
          3.800336e-13,   4.248215e-13,   1.571192e-12,
          2.360837e-12,   3.331525e-12,   4.497152e-12,
          5.875556e-12,   9.359056e-12,   1.398990e-11,
          2.000679e-11,   2.762451e-11,   3.691961e-11,
          5.352249e-11,   7.810201e-11,   9.194762e-11,
	  9.546290e-11,   1.258593e-10,   1.895055e-10,
          2.539315e-10,   3.194416e-10,   3.863255e-10,
          4.548506e-10,   5.252990e-10,   5.979231e-10,
          6.729756e-10,   7.507236e-10,   8.312713e-10,
          9.150863e-10,   1.002169e-09,   1.092816e-09,
          1.187250e-09,   1.285694e-09,   1.388296e-09,
          1.495280e-09,   1.606794e-09,   1.723060e-09,
          1.844228e-09,   1.970445e-09,   2.101934e-09,
          2.238770e-09,   2.381175e-09,   2.529225e-09,
          2.683140e-09,   2.842997e-09,   3.008942e-09,   3.181052e-09};

    // Insist that the size of SLYdat_p is the same as SLYdat_e
    ((void)sizeof(char[1 - 2*!!(sizeof(SLYdat_p) != sizeof(SLYdat_e))]));

    // Check the ammount of SLy data that will be included
    size_t limit;
    // Already checked that p0 > smallest p of crust eos, so start limit at 1
    for (limit = 1; limit < SLYdat_p_size; ++limit)
	if (psec[0] < SLYdat_p[limit])
            break;

    // 100 points is arbitrarily chosen to resolve high-density eos
    size_t ndat = 100 + limit;

    // Declaring array and allocating memory space
    double *edat;
    double *pdat;

    pdat = XLALCalloc(ndat, sizeof(*pdat));
    edat = XLALCalloc(ndat, sizeof(*edat));
    if (pdat == NULL || edat == NULL) {
        XLALFree(pdat);
	XLALFree(edat);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }

    double logpmax = log(pmax);
    double logp0 = log(psec[0]);
    double dlogp = (logpmax-logp0)/(ndat-limit);

    // Add SLy data for crust
    for (int i = 0; i < (int)limit; i++){
        pdat[i] = SLYdat_p[i];
        edat[i] = SLYdat_e[i];
    }

    // Find energy density at p0
    double b = SLYdat_e[limit-1];
    double dedp = (SLYdat_e[limit]-SLYdat_e[limit-1])/(SLYdat_p[limit]-SLYdat_p[limit-1]);
    double x = psec[0]-SLYdat_p[limit-1];
    double e0 = (dedp*x) + b;

    // Find starting energy's of each section
    double esec[nsec];
    esec[0] = e0;

    if (causal == 1){
        double vsec[nsec+1];
        double usec[nsec];
	// Solving Eq. (11) of PRD 97, 123019 (2018) for Y0
        usec[0] = dedp-1;
	// Solving Y(p0) = e^v0 for v0 = ln(Y0)
	vsec[0] = log(usec[0]);
	// Fill the rest of the v array
	vsec[1] = parameters[0];
        for (int i = 1; i < (int)nsec; i++)
            // Even indices
            vsec[i+1] = parameters[i*2];
        for (int k=1; k < (int)nsec; ++k){
	    // Eq. (15) of PRD 97, 123019 (2018)
            usec[k] = usec[k-1] * pow(psec[k] / psec[k-1],vsec[k]);
	    // Eq. (16) of PRD 97, 123019 (2018)
            esec[k] = esec[k-1] + psec[k] - psec[k-1] + ( psec[k-1] * usec[k-1] / (1 + vsec[k])) * (pow(psec[k] / psec[k-1], 1 + vsec[k]) - 1);
	}
        // Compute polytrope data beyond the crust
        for (int i = limit; i < (int)ndat; i++){
            pdat[i] = exp(logp0 + (dlogp*(i-limit)));
            edat[i] = causal_Nparams(pdat[i],esec,psec,usec,vsec,nsec);
        }
    }else{
        double gsec[nsec];
        gsec[0]=parameters[0];
        for (int i = 1; i < (int)nsec; i++)
	    // Even indices
            gsec[i] = parameters[i*2];
        for (int k = 0; k < (int)nsec-1; ++k)
            // Eq. (6) of PRD 97, 123019 (2018)
    	    esec[k+1] = (esec[k] - psec[k] / (gsec[k] - 1)) * pow(psec[k+1] / psec[k], 1 / gsec[k]) + psec[k+1] / (gsec[k] - 1);
        // Compute polytrope data beyond the crust
        for (int i = limit; i < (int)ndat; i++){
            pdat[i] = exp(logp0 + (dlogp*(i-limit)));
            edat[i] = non_causal_Nparams(pdat[i],esec,psec,gsec,nsec);
        }
    }

    // Updating eos structure
    LALSimNeutronStarEOS * eos;
    double *nbdat = NULL;
    double *mubdat = NULL;
    double *muedat = NULL;
    double *yedat = NULL;
    double *cs2dat = NULL;
    double *hdat = NULL;
    size_t ncol = 2;
    /* Updating eos structure */
    eos = eos_alloc_tabular(nbdat, edat, pdat, mubdat, muedat, hdat, yedat, cs2dat, ndat, ncol);


    XLALFree(edat);
    XLALFree(pdat);

    return eos;
}

/**
 * @brief Reads 5 dynamic polytrope eos parameters to make an eos.
 * @details Reads 5 dynamic polytrope eos parameters to construct a 3-piece
 * dynamic polytrope eos, generically outlined in PRD 97, 123019 (2018), and
 * stitched to a low-density SLy eos crust.
 * @param[in] g0 The adiabatic index of the first polytrope.
 * @param[in] log10p1_si The log10 of the first dividing pressure in SI units.
 * @param[in] g1 The adiabatic index of the second polytrope.
 * @param[in] log10p2_si The log10 of the second dividing pressure in SI units.
 * @param[in] g2 The adiabatic index of the third polytrope.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceDynamicPolytrope(double g0, double log10p1_si, double g1, double log10p2_si, double g2){
    double params[]={g0,log10p1_si,g1,log10p2_si,g2};
    LALSimNeutronStarEOS * eos;
    // 3-piece (non-causal) dynamic polytrope
    eos = XLALSimNeutronStarEOSDynamicAnalytic(params, 3, 0);
    return eos;
}

/**
 * @brief Reads 5 causal analytic eos parameters to make an eos.
 * @details Reads 5 causal analytic eos parameters to construct a 3-piece
 * causal eos, generically outlined in PRD 97, 123019 (2018), and stitched to
 * a low-density SLy eos crust.
 * @param[in] v1 The first sound function exponent.
 * @param[in] log10p1_si The log10 of the first dividing pressure in SI units.
 * @param[in] v2 The second sound function exponent.
 * @param[in] log10p2_si The log10 of the second dividing pressure in SI units.
 * @param[in] v3 The third sound function exponent.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceCausalAnalytic(double v1, double log10p1_si, double v2, double log10p2_si, double v3){
    double params[]={v1,log10p1_si,v2,log10p2_si,v3};
    LALSimNeutronStarEOS * eos;
    // 3-piece causal analytic eos
    eos = XLALSimNeutronStarEOSDynamicAnalytic(params, 3, 1);
    return eos;
}

/**
 * @brief Check that EOS has enough points (>4) in M-R space to interpolate.
 * As the TOV equations are integrated from pmin to pmax, masses are
 * calculated. Once the mass turns over (i.e. m(p_i) > m(p_{i+1})), the
 * integration is terminated. It is possible therefore to have less than 4
 * points in M-R space. We demand, then, that in the interval [pmin,pmax],
 * m(pmin) = m(p_0) < m(p_1) < m(p_2) < m(p_3).
 * @details Reads 3-piece dynamic polytrope (non-causal) or analytic (causal)
 * eos parameters and checks that the eos has at least 4 data points before
 * turning over in M-R space, which abruptly terminates the eos.
 * @param[in] p0 First model param, either g0 (non-causal) or v1 (causal).
 * @param[in] log10p1_si The log10 of the first dividing pressure in SI units.
 * @param[in] p1 Second model param, either g1 (non-causal) or v2 (causal).
 * @param[in] log10p2_si The log10 of the second dividing pressure in SI units.
 * @param[in] p2 Third model param, either g2 (non-causal) or v3 (causal).
 * @param[in] causal Option to use causal version (0=non-causal, 1=causal).
 * @return XLAL_SUCCESS or XLAL_FAILURE.
 */
int XLALSimNeutronStarEOS3PDViableFamilyCheck(double p0, double log10p1_si, double p1, double log10p2_si, double p2, int causal){
    // Check that causal flag appropriately assigned
    if (causal!=0 && causal!=1)
        XLAL_ERROR(XLAL_EINVAL,"Did not specify which approach to take, Causal or Non-Causal");
    double pdat;
    double mdat;
    double mdat_prev;
    double rdat;
    double kdat;

    LALSimNeutronStarEOS *eos=NULL;
    if (causal == 1)
        eos = XLALSimNeutronStarEOS3PieceCausalAnalytic(p0, log10p1_si, p1, log10p2_si, p2);
    else
        eos = XLALSimNeutronStarEOS3PieceDynamicPolytrope(p0, log10p1_si, p1, log10p2_si, p2);

    // Initialize previous value for mdat comparison
    mdat_prev = 0.0;

    // Ensure mass turnover does not happen too soon
    const double logpmin = 75.5;
    double logpmax = log(XLALSimNeutronStarEOSMaxPressure(eos));
    double dlogp = (logpmax - logpmin) / 100.;
    // Need at least four points
    for (int i = 0; i < 4; ++i) {
        pdat = exp(logpmin + i * dlogp);
        XLALSimNeutronStarTOVODEIntegrate(&rdat, &mdat, &kdat, pdat, eos);
        // Determine if maximum mass has been found
        if (mdat <= mdat_prev){
            // EOS has too few points to create family
            // Clean up
            XLALDestroySimNeutronStarEOS(eos);
            return XLAL_FAILURE;
        }
        mdat_prev = mdat;
    }

    XLALDestroySimNeutronStarEOS(eos);

    return XLAL_SUCCESS;
}

/** @} */
/** @} */
