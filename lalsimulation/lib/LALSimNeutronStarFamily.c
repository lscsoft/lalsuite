/*
 * Copyright (C) 2013 J. Creighton
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
 * @author Jolien Creighton
 * @addtogroup LALSimNeutronStarFamily_c
 * @brief Provides routines for one-parameter families of neutron stars of
 * a given equation of state.
 * @{
 */

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_min.h>
GSL_VAR const gsl_interp_type * lal_gsl_interp_steffen;

#include <lal/LALStdlib.h>
#include <lal/LALSimNeutronStar.h>

/** @cond */

/* Contents of the neutron star family structure. */
struct tagLALSimNeutronStarFamily {
    double *pdat;
    double *mdat;
    double *mbdat;
    double *rdat;
    double *k2dat;
    double *k3dat;
    double *k4dat;
    size_t ndat;
    gsl_interp *m_of_p_interp;
    gsl_interp *mb_of_m_interp;
    gsl_interp *p_of_m_interp;
    gsl_interp *r_of_m_interp;
    gsl_interp *k2_of_m_interp;
    gsl_interp *k3_of_m_interp;
    gsl_interp *k4_of_m_interp;
    gsl_interp_accel *m_of_p_acc;
    gsl_interp_accel *mb_of_m_acc;
    gsl_interp_accel *p_of_m_acc;
    gsl_interp_accel *r_of_m_acc;
    gsl_interp_accel *k2_of_m_acc;
    gsl_interp_accel *k3_of_m_acc;
    gsl_interp_accel *k4_of_m_acc;
};

struct tagFamMultiParts{
  int number_of_branches;
  double pmin;
  double pmax;
  LALSimNeutronStarFamily ** fam_branch;
};

/* gsl function for use in finding the maximum neutron star mass */
static double fminimizer_gslfunction(double x, void * params);
static double fminimizer_gslfunction(double x, void * params)
{
    LALSimNeutronStarEOS * eos = params;
    double r, m, k;
    XLALSimNeutronStarTOVODEIntegrate(&r, &m, &k, x, eos);
    return -m; /* maximum mass is minimum negative mass */
}

/* gsl function for use in finding the maximum neutron star mass */
static double fmaximizer_multi_gslfunction(double x, void * params);
static double fmaximizer_multi_gslfunction(double x, void * params)
{
    EOSMultiParts * eos = params;
    double r, m, mb, k2, k3, k4;
    XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
        &r, &m, &mb, &k2, &k3, &k4, x, eos, 1e-6, 0);
    return -m; /* maximum mass is minimum negative mass */
}

/* gsl function for use in finding the minimum neutron star mass */
static double fminimizer_multi_gslfunction(double x, void * params);
static double fminimizer_multi_gslfunction(double x, void * params)
{
    EOSMultiParts * eos = params;
    double r, m, mb, k2, k3, k4;
    XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
        &r, &m, &mb, &k2, &k3, &k4, x, eos, 1e-6, 0);
    return m;
}

/** @endcond */

/**
 * @brief Frees the memory associated with a pointer to a neutron star family.
 * @param fam Pointer to the neutron star family structure to be freed.
 */
void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam)
{
    if (fam) {
        gsl_interp_accel_free(fam->k2_of_m_acc);
	gsl_interp_accel_free(fam->k3_of_m_acc);
	gsl_interp_accel_free(fam->k4_of_m_acc);
        gsl_interp_accel_free(fam->r_of_m_acc);
        gsl_interp_accel_free(fam->p_of_m_acc);
	gsl_interp_accel_free(fam->mb_of_m_acc);
	gsl_interp_accel_free(fam->m_of_p_acc);
        gsl_interp_free(fam->k2_of_m_interp);
	gsl_interp_free(fam->k3_of_m_interp);
	gsl_interp_free(fam->k4_of_m_interp);
        gsl_interp_free(fam->r_of_m_interp);
        gsl_interp_free(fam->p_of_m_interp);
	gsl_interp_free(fam->mb_of_m_interp);
	gsl_interp_free(fam->m_of_p_interp);
        LALFree(fam->k2dat);
	LALFree(fam->k3dat);
	LALFree(fam->k4dat);
        LALFree(fam->rdat);
        LALFree(fam->mdat);
	LALFree(fam->mbdat);
        LALFree(fam->pdat);
        LALFree(fam);
    }
    return;
}

/**
 * @brief Frees the memory associated with a pointer to a neutron star multi-branch family.
 * @param fam Pointer to the neutron star multi-branch family structure to be freed.
 */
void XLALDestroySimNeutronStarMultiBranchFamily(FamMultiParts * fam)
{
    if (fam) {
	for(int i = 0; i < fam->number_of_branches; i++){
	    LALSimNeutronStarFamily * branch = fam->fam_branch[i];
	    XLALDestroySimNeutronStarFamily(branch);
        }
        LALFree(fam);
    }
    return;
}

/**
 * @brief Creates a neutron star family structure for a given equation of state.
 * @details
 * A neutron star family is a one-parameter family of neturon stars for a
 * fixed equation of state.  The one parameter is the neutron star central
 * pressure, or, equivalently, the mass of the neutron star.  The family
 * is terminated at the maximum neutron star mass for the specified equation
 * of state, so the mass can be used as the family parameter.
 * @param eos Pointer to the Equation of State structure.
 * @return A pointer to the neutron star family structure.
 */
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(
    LALSimNeutronStarEOS * eos)
{
    LALSimNeutronStarFamily * fam;
    const size_t ndatmax = 100;
    const double logpmin = 75.5;
    double logpmax;
    double dlogp;
    size_t ndat = ndatmax;
    size_t i;

    /* allocate memory */
    fam = LALMalloc(sizeof(*fam));
    if (!fam)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    fam->pdat = LALMalloc(ndat * sizeof(*fam->pdat));
    fam->mdat = LALMalloc(ndat * sizeof(*fam->mdat));
    fam->rdat = LALMalloc(ndat * sizeof(*fam->rdat));
    fam->k2dat = LALMalloc(ndat * sizeof(*fam->k2dat));
    if (!fam->mdat || !fam->rdat || !fam->k2dat)
        XLAL_ERROR_NULL(XLAL_ENOMEM);

    /* compute data tables */
    logpmax = log(XLALSimNeutronStarEOSMaxPressure(eos));
    dlogp = (logpmax - logpmin) / ndat;
    for (i = 0; i < ndat; ++i) {
        fam->pdat[i] = exp(logpmin + i * dlogp);
        XLALSimNeutronStarTOVODEIntegrate(&fam->rdat[i], &fam->mdat[i],
            &fam->k2dat[i], fam->pdat[i], eos);
        /* determine if maximum mass has been found */
        if (i > 0 && fam->mdat[i] <= fam->mdat[i-1])
            break;
    }

    if (i < ndat) {
        /* replace the ith point with the maximum mass */
        const double epsabs = 0.0, epsrel = 1e-6;
        double a = fam->pdat[i - 2];
        double x = fam->pdat[i - 1];
        double b = fam->pdat[i];
        double fa = -fam->mdat[i - 2];
        double fx = -fam->mdat[i - 1];
        double fb = -fam->mdat[i];
        int status;
        gsl_function F;
        gsl_min_fminimizer * s;
        F.function = &fminimizer_gslfunction;
        F.params = eos;
        s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
        gsl_min_fminimizer_set_with_values(s, &F, x, fx, a, fa, b, fb);
        do {
            status = gsl_min_fminimizer_iterate(s);
            x = gsl_min_fminimizer_x_minimum(s);
            a = gsl_min_fminimizer_x_lower(s);
            b = gsl_min_fminimizer_x_upper(s);
            status = gsl_min_test_interval(a, b, epsabs, epsrel);
        } while (status == GSL_CONTINUE);
        gsl_min_fminimizer_free(s);
        fam->pdat[i] = x;
        XLALSimNeutronStarTOVODEIntegrate(&fam->rdat[i], &fam->mdat[i],
            &fam->k2dat[i], fam->pdat[i], eos);

        /* resize arrays */
        if(fam->pdat[i] <= fam->pdat[i-1]){
            fam->pdat[i-1] = fam->pdat[i];
            ndat = i;
        }
        else{
            ndat = i + 1;
        }

        fam->pdat = LALRealloc(fam->pdat, ndat * sizeof(*fam->pdat));
        fam->mdat = LALRealloc(fam->mdat, ndat * sizeof(*fam->mdat));
        fam->rdat = LALRealloc(fam->rdat, ndat * sizeof(*fam->rdat));
        fam->k2dat = LALRealloc(fam->k2dat, ndat * sizeof(*fam->k2dat));
    }
    fam->ndat = ndat;

    /* setup interpolators */

    fam->p_of_m_acc = gsl_interp_accel_alloc();
    fam->r_of_m_acc = gsl_interp_accel_alloc();
    fam->k2_of_m_acc = gsl_interp_accel_alloc();

    fam->p_of_m_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    fam->r_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    fam->k2_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);

    gsl_interp_init(fam->p_of_m_interp, fam->mdat, fam->pdat, ndat);
    gsl_interp_init(fam->r_of_m_interp, fam->mdat, fam->rdat, ndat);
    gsl_interp_init(fam->k2_of_m_interp, fam->mdat, fam->k2dat, ndat);

    return fam;
}



FamMultiParts * XLALCreateSimNeutronStarFamilyPT(EOSMultiParts * eos){

    FamMultiParts *fam;
    fam = LALMalloc(sizeof(*fam));

    double logpmin = 75.5; // pmin = 6.16e32
    int nbranch_max = 4; //FIXME: Is 4 enough? Too many?
    int ndat = 100;
    int nbranch = 1;
    int min_interp_points = 5; // Minimum number of points required for all gls interpolators

    fam->pmin = exp(logpmin);
    fam->pmax = XLALSimNeutronStarEOSMultiPartsMaxPressureGeometerized(eos);

    //FIXME: Is 4 enough? Too many?
    fam->fam_branch = (LALSimNeutronStarFamily **) LALMalloc(sizeof(LALSimNeutronStarFamily *) * nbranch_max);
    fam->number_of_branches = nbranch;

    double pmax_pascal, pmax_geo, logpmax, dlogp;
    //double pmax_nuc, pmax_cgs;
    pmax_geo = fam->pmax;
    pmax_pascal = pmax_geo / LAL_G_C4_SI; // in Pascal
    //pmax_cgs = pmax_pascal * 10; // in CGS units
    //pmax_nuc = pmax_cgs / 1.6022e33; // in nuclear units
    logpmax = log(pmax_pascal);

    int turnover = 0;
    int cont = 0; // If not incremented, branch loop will stop after one iteration
    int b = 0;
    int j = 0;

    for (b = 0; b < nbranch; b++){

	LALSimNeutronStarFamily * fam_branch_i = LALMalloc(sizeof(LALSimNeutronStarFamily));

        fam_branch_i->pdat = LALMalloc(ndat * sizeof(*fam_branch_i->pdat));
        fam_branch_i->mdat = LALMalloc(ndat * sizeof(*fam_branch_i->mdat));
	fam_branch_i->mbdat = LALMalloc(ndat * sizeof(*fam_branch_i->mbdat));
        fam_branch_i->rdat = LALMalloc(ndat * sizeof(*fam_branch_i->rdat));
        fam_branch_i->k2dat = LALMalloc(ndat * sizeof(*fam_branch_i->k2dat));
	fam_branch_i->k3dat = LALMalloc(ndat * sizeof(*fam_branch_i->k3dat));
	fam_branch_i->k4dat = LALMalloc(ndat * sizeof(*fam_branch_i->k4dat));
	fam_branch_i->ndat = ndat;

        if (!fam_branch_i->mdat || !fam_branch_i->rdat || !fam_branch_i->k2dat
	    || !fam_branch_i->mbdat || !fam_branch_i->k3dat || !fam_branch_i->k4dat)
            XLAL_ERROR_NULL(XLAL_ENOMEM);

	dlogp = (logpmax - logpmin) / ndat; // in Pascal

        for (j = 0; j < ndat; ++j) {

            fam_branch_i->pdat[j] = exp(logpmin + j * dlogp); // pressure in Pascal

	    XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
	        &fam_branch_i->rdat[j], &fam_branch_i->mdat[j],
		&fam_branch_i->mbdat[j], &fam_branch_i->k2dat[j],
		&fam_branch_i->k3dat[j], &fam_branch_i->k4dat[j],
		fam_branch_i->pdat[j], eos, 1e-6, 0);

	    // End of stable branch; Coming down from max mass
	    if (j > 0 && turnover == 0 && fam_branch_i->mdat[j] <= fam_branch_i->mdat[j-1])
		turnover = j-1;

	    // End of unstable branch; Rising from min mass of unstable branch
            if (j > 2 && turnover > 0 && fam_branch_i->mdat[j] >= fam_branch_i->mdat[j-1]){
                // Find precise minimum mass of next branch
                const double epsabs = 0.0, epsrel = 1e-6;
                double a = fam_branch_i->pdat[j - 2];
                double x = fam_branch_i->pdat[j - 1];
                double c = fam_branch_i->pdat[j];
                double fa = fam_branch_i->mdat[j - 2];
                double fx = fam_branch_i->mdat[j - 1];
                double fc = fam_branch_i->mdat[j];
                int status;
                gsl_function F;
                gsl_min_fminimizer * s;
                F.function = &fminimizer_multi_gslfunction;
                F.params = eos;
                s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
                gsl_min_fminimizer_set_with_values(s, &F, x, fx, a, fa, c, fc);
                do {
                    status = gsl_min_fminimizer_iterate(s);
                    x = gsl_min_fminimizer_x_minimum(s);
                    a = gsl_min_fminimizer_x_lower(s);
                    c = gsl_min_fminimizer_x_upper(s);
                    status = gsl_min_test_interval(a, c, epsabs, epsrel);
                } while (status == GSL_CONTINUE);
                gsl_min_fminimizer_free(s);
                // Set new logpmin to be at minimum mass
	        logpmin = log(x);
		cont = 1;
                ndat -= j-1;
	        break;
	    }
	}

	if (turnover > min_interp_points){

            // Find precise maximum mass of branch
            int i = turnover+1;
            /* replace the ith point with the maximum mass */
            const double epsabs = 0.0, epsrel = 1e-6;
            double a = fam_branch_i->pdat[i - 2];
            double x = fam_branch_i->pdat[i - 1];
            double c = fam_branch_i->pdat[i];
            double fa = -fam_branch_i->mdat[i - 2];
            double fx = -fam_branch_i->mdat[i - 1];
            double fc = -fam_branch_i->mdat[i];
            int status;
            gsl_function F;
            gsl_min_fminimizer * s;
            F.function = &fmaximizer_multi_gslfunction;
            F.params = eos;
            s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
            gsl_min_fminimizer_set_with_values(s, &F, x, fx, a, fa, c, fc);
            do {
                status = gsl_min_fminimizer_iterate(s);
                x = gsl_min_fminimizer_x_minimum(s);
                a = gsl_min_fminimizer_x_lower(s);
                c = gsl_min_fminimizer_x_upper(s);
                status = gsl_min_test_interval(a, c, epsabs, epsrel);
            } while (status == GSL_CONTINUE);
            gsl_min_fminimizer_free(s);
            fam_branch_i->pdat[i] = x;
	    XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
                &fam_branch_i->rdat[i], &fam_branch_i->mdat[i],
                &fam_branch_i->mbdat[i], &fam_branch_i->k2dat[i],
		&fam_branch_i->k3dat[i], &fam_branch_i->k4dat[i],
                fam_branch_i->pdat[i], eos, 1e-6, 0);

            /* resize arrays */
            if(fam_branch_i->pdat[i] <= fam_branch_i->pdat[i-1]){
                fam_branch_i->pdat[i-1] = fam_branch_i->pdat[i];
                turnover = i;
            }
            else{
                turnover = i + 1;
            }

	    // Reset length
            fam_branch_i->ndat = turnover;

	    // Reset size of arrays
            fam_branch_i->pdat = LALRealloc(fam_branch_i->pdat, turnover * sizeof(*fam_branch_i->pdat));
            fam_branch_i->mdat = LALRealloc(fam_branch_i->mdat, turnover * sizeof(*fam_branch_i->mdat));
	    fam_branch_i->mbdat = LALRealloc(fam_branch_i->mbdat, turnover * sizeof(*fam_branch_i->mbdat));
            fam_branch_i->rdat = LALRealloc(fam_branch_i->rdat, turnover * sizeof(*fam_branch_i->rdat));
            fam_branch_i->k2dat = LALRealloc(fam_branch_i->k2dat, turnover * sizeof(*fam_branch_i->k2dat));
	    fam_branch_i->k3dat = LALRealloc(fam_branch_i->k3dat, turnover * sizeof(*fam_branch_i->k3dat));
	    fam_branch_i->k4dat = LALRealloc(fam_branch_i->k4dat, turnover * sizeof(*fam_branch_i->k4dat));

	    // Set up p(m), r(m), and k(m) interpolators
	    fam_branch_i->m_of_p_acc = gsl_interp_accel_alloc();
            fam_branch_i->p_of_m_acc = gsl_interp_accel_alloc();
	    fam_branch_i->mb_of_m_acc = gsl_interp_accel_alloc();
            fam_branch_i->r_of_m_acc = gsl_interp_accel_alloc();
            fam_branch_i->k2_of_m_acc = gsl_interp_accel_alloc();
	    fam_branch_i->k3_of_m_acc = gsl_interp_accel_alloc();
	    fam_branch_i->k4_of_m_acc = gsl_interp_accel_alloc();

	    fam_branch_i->m_of_p_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);
            fam_branch_i->p_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);
	    fam_branch_i->mb_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);
            fam_branch_i->r_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);
            fam_branch_i->k2_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);
	    fam_branch_i->k3_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);
	    fam_branch_i->k4_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, turnover);

	    gsl_interp_init(fam_branch_i->m_of_p_interp, fam_branch_i->pdat, fam_branch_i->mdat, turnover);
            gsl_interp_init(fam_branch_i->p_of_m_interp, fam_branch_i->mdat, fam_branch_i->pdat, turnover);
	    gsl_interp_init(fam_branch_i->mb_of_m_interp, fam_branch_i->mdat, fam_branch_i->mbdat, turnover);
            gsl_interp_init(fam_branch_i->r_of_m_interp, fam_branch_i->mdat, fam_branch_i->rdat, turnover);
            gsl_interp_init(fam_branch_i->k2_of_m_interp, fam_branch_i->mdat, fam_branch_i->k2dat, turnover);
	    gsl_interp_init(fam_branch_i->k3_of_m_interp, fam_branch_i->mdat, fam_branch_i->k3dat, turnover);
	    gsl_interp_init(fam_branch_i->k4_of_m_interp, fam_branch_i->mdat, fam_branch_i->k4dat, turnover);

            if (cont==1){
                turnover = 0; // Reset turnover location
                nbranch += 1; // Increment the number of branches to allow loop to continue
                cont = 0; // Reset continue variable
            }

	} else if (ndat > min_interp_points){

            // Set up p(m), r(m), and k(m) interpolators
	    fam_branch_i->m_of_p_acc = gsl_interp_accel_alloc();
            fam_branch_i->p_of_m_acc = gsl_interp_accel_alloc();
	    fam_branch_i->mb_of_m_acc = gsl_interp_accel_alloc();
            fam_branch_i->r_of_m_acc = gsl_interp_accel_alloc();
            fam_branch_i->k2_of_m_acc = gsl_interp_accel_alloc();
	    fam_branch_i->k3_of_m_acc = gsl_interp_accel_alloc();
	    fam_branch_i->k4_of_m_acc = gsl_interp_accel_alloc();

	    fam_branch_i->m_of_p_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
            fam_branch_i->p_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
	    fam_branch_i->mb_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
            fam_branch_i->r_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
            fam_branch_i->k2_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
	    fam_branch_i->k3_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
	    fam_branch_i->k4_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);

	    gsl_interp_init(fam_branch_i->m_of_p_interp, fam_branch_i->pdat, fam_branch_i->mdat, ndat);
            gsl_interp_init(fam_branch_i->p_of_m_interp, fam_branch_i->mdat, fam_branch_i->pdat, ndat);
	    gsl_interp_init(fam_branch_i->mb_of_m_interp, fam_branch_i->mdat, fam_branch_i->mbdat, ndat);
            gsl_interp_init(fam_branch_i->r_of_m_interp, fam_branch_i->mdat, fam_branch_i->rdat, ndat);
            gsl_interp_init(fam_branch_i->k2_of_m_interp, fam_branch_i->mdat, fam_branch_i->k2dat, ndat);
	    gsl_interp_init(fam_branch_i->k3_of_m_interp, fam_branch_i->mdat, fam_branch_i->k3dat, ndat);
	    gsl_interp_init(fam_branch_i->k4_of_m_interp, fam_branch_i->mdat, fam_branch_i->k4dat, ndat);

	}

        fam->fam_branch[b] = fam_branch_i;

	// Exit if too many branches
	if (nbranch > nbranch_max){
	    nbranch = nbranch_max;
	    break;
	}
    }

    // FIXME: Reset number of branches
    //if (nbranch < nbranch_max)
    //    fam->fam_branch = (LALSimNeutronStarFamily **) LALRealloc(sizeof(LALSimNeutronStarFamily *) * nbranch);
    fam->number_of_branches = nbranch;

    return fam;
}


int XLALSimNeutronStarFamNumberOfBranches(FamMultiParts *fam)
{
    return fam->number_of_branches;
}

double XLALSimNeutronStarFamBranchMinMass(int branch, FamMultiParts *fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->mdat[0];
}

double XLALSimNeutronStarFamBranchMinCentralPressure(int branch, FamMultiParts *fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->pdat[0];
}

double XLALSimNeutronStarFamBranchMaxMass(int branch, FamMultiParts *fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->mdat[b->ndat - 1];
}

double XLALSimNeutronStarFamBranchMaxCentralPressure(int branch, FamMultiParts *fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->pdat[b->ndat - 1];
}

double XLALSimNeutronStarFamBranchRadius(double m, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double r;
    r = gsl_interp_eval(b->r_of_m_interp, b->mdat, b->rdat, m,
        b->r_of_m_acc);
    return r;
}

double XLALSimNeutronStarFamBranchCentralPressure(double m, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double p;
    p = gsl_interp_eval(b->p_of_m_interp, b->mdat, b->pdat, m,
        b->p_of_m_acc);
    return p;
}

double XLALSimNeutronStarFamBranchMass(double p, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double m;
    m = gsl_interp_eval(b->m_of_p_interp, b->pdat, b->mdat, p,
        b->m_of_p_acc);
    return m;
}

double XLALSimNeutronStarFamBranchBaryonicMass(double m, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double mb;
    mb = gsl_interp_eval(b->mb_of_m_interp, b->mdat, b->mbdat, m,
        b->mb_of_m_acc);
    return mb;
}

double XLALSimNeutronStarFamBranchLoveNumberK2(double m, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double k2;
    k2 = gsl_interp_eval(b->k2_of_m_interp, b->mdat, b->k2dat, m,
        b->k2_of_m_acc);
    return k2;
}

double XLALSimNeutronStarFamBranchLoveNumberK3(double m, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double k3;
    k3 = gsl_interp_eval(b->k3_of_m_interp, b->mdat, b->k3dat, m,
        b->k3_of_m_acc);
    return k3;
}

double XLALSimNeutronStarFamBranchLoveNumberK4(double m, int branch, FamMultiParts * fam)
{
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double k4;
    k4 = gsl_interp_eval(b->k4_of_m_interp, b->mdat, b->k4dat, m,
        b->k4_of_m_acc);
    return k4;
}


/**
 * @brief Returns the minimum mass of a neutron star family.
 * @param fam Pointer to the neutron star family structure.
 * @return The maximum mass of the neutron star family (kg).
 */
double XLALSimNeutronStarFamMinimumMass(LALSimNeutronStarFamily *fam)
{
    return fam->mdat[0];
}

/**
 * @brief Returns the maximum mass of a neutron star family.
 * @param fam Pointer to the neutron star family structure.
 * @return The maximum mass of the neutron star family (kg).
 */
double XLALSimNeutronStarMaximumMass(LALSimNeutronStarFamily * fam)
{
    return fam->mdat[fam->ndat - 1];
}

/**
 * @brief Returns the central pressure of a neutron star of mass @a m.
 * @param m The mass of the neutron star (kg).
 * @param fam Pointer to the neutron star family structure.
 * @return The central pressure of the neutron star (Pa).
 */
double XLALSimNeutronStarCentralPressure(double m,
    LALSimNeutronStarFamily * fam)
{
    double p;
    p = gsl_interp_eval(fam->p_of_m_interp, fam->mdat, fam->pdat, m,
        fam->p_of_m_acc);
    return p;
}

/**
 * @brief Returns the radius of a neutron star of mass @a m.
 * @param m The mass of the neutron star (kg).
 * @param fam Pointer to the neutron star family structure.
 * @return The radius of the neutron star (m).
 */
double XLALSimNeutronStarRadius(double m, LALSimNeutronStarFamily * fam)
{
    double r;
    r = gsl_interp_eval(fam->r_of_m_interp, fam->mdat, fam->rdat, m,
        fam->r_of_m_acc);
    return r;
}

/**
 * @brief Returns the tidal Love number k2 of a neutron star of mass @a m.
 * @param m The mass of the neutron star (kg).
 * @param fam Pointer to the neutron star family structure.
 * @return The dimensionless tidal Love number k2.
 */
double XLALSimNeutronStarLoveNumberK2(double m, LALSimNeutronStarFamily * fam)
{
    double k2;
    k2 = gsl_interp_eval(fam->k2_of_m_interp, fam->mdat, fam->k2dat, m,
        fam->k2_of_m_acc);
    return k2;
}

/** @} */
