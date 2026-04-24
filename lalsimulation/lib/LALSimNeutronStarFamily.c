/*
 * Copyright (C) 2013 J. Creighton
 * Copyright (C) 2026 P. Davis, M. Oertel, L. Suleiman, L. Wade
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
 * @author Jolien Creighton, Philip Davis, Micaela Oertel, Lami Suleiman, Leslie Wade
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
  int mtov;
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

/* Function used by GSLminimizer to find the central pressure for the
 * exact maximum or minimum neutron star mass from EOSMultiParts structure.
 */
static double fextrimizer(double x, void * params)
{
    EOSMultiParts * eos = params;
    double r, m, k2;
    XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(
        &r, &m, &k2, x, eos, 1e-6);
    return m;
}

static double get_central_pressure_mext(EOSMultiParts *eos, int index, double *pdat, double *mdat, int factor){
    // Recalculate the maximum mass of a branch
    const double epsabs = 0.0, epsrel = 1e-6;
    double a = pdat[index-1];
    double x = pdat[index];
    double c = pdat[index+1];
    double fa = factor*mdat[index-1];
    double fx = factor*mdat[index];
    double fc = factor*mdat[index+1];
    int status;
    gsl_function F;
    gsl_min_fminimizer * s;
    F.function = &fextrimizer;
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
    return x;
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

    /* allocate memory (calloc zeroes all fields so unused pointers are NULL) */
    fam = LALCalloc(1, sizeof(*fam));
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
    dlogp = (logpmax - logpmin) / (ndat-1);
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

/**
 * @brief Creates a neutron star family structure (FamMultiParts) which can accomodate
 * multiple branches for a given EOSMultiParts equation of state structure.
 * @details
 * The neutron star family contains the astrophysical parameter solution for a
 * given equation of state. When the structure is created, it solves the TOV+Love number
 * ODEs with XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance if min_fam = 0, or with
 * XLALSimNeutronStarTOVODEMiniIntegrate if min_fam=1. The minimum central pressure
 * to solve the neutron star astrophysical parameters is min(Pc) = exp(75.5) Pa.
 * The maximum central pressure is that defined as the largest pressure of the equation of state.
 * The neutron star family can accomodate multiple branches and provides interpolating
 * function for astrophysical quantities only for stable branches (continuously increasing
 * mass in a sequence).
 * @param eos Pointer to the Equation of State structure.
 * @param min_fam Flag to choose to construct the neutron star family with only the minimum
 * of astrophysical parameters, M, R and k2 (using XLALSimNeutronStarTOVODEMiniIntegrate).
 * @return A pointer to the neutron star family structure.
 */
FamMultiParts * XLALCreateSimNeutronStarFamilyPT(EOSMultiParts * eos, int min_fam){

    if(min_fam!=0 && min_fam!=1) XLAL_ERROR_NULL(XLAL_EDOM);

    FamMultiParts *fam;
    fam = LALCalloc(1, sizeof(*fam));

    int ndat = 100;
    int min_interp_points = 5; // Minimum number of points required for all gls interpolators
    double logpmin = 75.5; // pmin = 6.16e32
    fam->pmin = exp(logpmin);
    fam->pmax = XLALSimNeutronStarEOSMultiPartsMaxPressure(eos);
    double logpmax = log(fam->pmax);
    double dlogp = (logpmax - logpmin) / (ndat-1); // in Pascal
    fam->mtov = 0;

    double *pdat, *rdat, *mdat, *mbdat, *k2dat, *k3dat, *k4dat;
    pdat = LALMalloc(ndat * sizeof(*pdat));
    mdat = LALMalloc(ndat * sizeof(*mdat));
    rdat = LALMalloc(ndat * sizeof(*rdat));
    mbdat = LALMalloc(ndat * sizeof(*mbdat));
    k2dat = LALMalloc(ndat * sizeof(*k2dat));
    k3dat = LALMalloc(ndat * sizeof(*k3dat));
    k4dat = LALMalloc(ndat * sizeof(*k4dat));

    int nb_stable_branches = 0;
    int flag_unstable = 0;
    int flag_stable = 0;
    int index_begin_stable_branch[100] = { [ 0 ... 99 ] = 0};
    int index_end_stable_branch[100] = { [ 0 ... 99 ] = 0};
    for (int j = 0; j < ndat; ++j) {
        // find pdat being careful with first and last sample
    	if(j==0) pdat[j] = fam->pmin; // exact min
	else if(j==ndat-1) pdat[j] = fam->pmax; // exact max
	else pdat[j] = exp(logpmin + j * dlogp); // pressure in Pascal
        // Solve TOV+Love number ODEs
        if(min_fam==1){
            XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(
                &rdat[j], &mdat[j],
                &k2dat[j], pdat[j], eos, 1e-6);
        } else {
            XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
                &rdat[j], &mdat[j],
                &mbdat[j], &k2dat[j],
                &k3dat[j], &k4dat[j],
                pdat[j], eos, 1e-6);
        }
        // Determine stable and unstable branches
        if (j > 0){
            if (mdat[j] > mdat[j-1]){ // Stable branch
                if (flag_stable == 0) {
                    index_begin_stable_branch[nb_stable_branches] = j-1;
                    nb_stable_branches += 1;
                }
                flag_stable = 1;
                flag_unstable = 0;
            } else {
                if (flag_unstable == 0){
                    if (nb_stable_branches != 0) index_end_stable_branch[nb_stable_branches-1] = j-1;
                }
                flag_stable = 0;
                flag_unstable = 1;
            }
        }

    }

    // If the sequence ends on a stable branch
    fam->number_of_branches = nb_stable_branches;
    fam->mtov = 0;
    if (nb_stable_branches != 0) {
  	// if the last branch is stable, append the index for the end of the stable branch as the last index of TOV data
        if (flag_stable == 1) index_end_stable_branch[nb_stable_branches-1] = ndat-1;
        fam->fam_branch = (LALSimNeutronStarFamily **) LALMalloc(sizeof(LALSimNeutronStarFamily *) * nb_stable_branches);
        int *ndat_branch;
        ndat_branch = LALMalloc(nb_stable_branches * sizeof(*ndat_branch));
        for (int b = 0; b < nb_stable_branches; b++){
            ndat_branch[b] = index_end_stable_branch[b] - index_begin_stable_branch[b] + 1;
            LALSimNeutronStarFamily * fam_branch_i = LALCalloc(1, sizeof(LALSimNeutronStarFamily));
            if(min_fam==1){
                fam_branch_i->mbdat = NULL;
                fam_branch_i->k3dat = NULL;
                fam_branch_i->k4dat = NULL;
                fam_branch_i->mb_of_m_interp = NULL;
                fam_branch_i->k3_of_m_interp = NULL;
                fam_branch_i->k4_of_m_interp = NULL;
                fam_branch_i->mb_of_m_acc = NULL;
                fam_branch_i->k3_of_m_acc = NULL;
                fam_branch_i->k4_of_m_acc = NULL;
            }


            // Redefine ndat as the number of points per branch
            fam_branch_i->ndat = ndat_branch[b];
            fam_branch_i->pdat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->pdat));
            fam_branch_i->mdat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->mdat));
            fam_branch_i->rdat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->rdat));
            fam_branch_i->k2dat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->k2dat));
            if(min_fam==0){
                fam_branch_i->mbdat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->mbdat));
                fam_branch_i->k3dat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->k3dat));
                fam_branch_i->k4dat = LALMalloc(ndat_branch[b] * sizeof(*fam_branch_i->k4dat));
            }

            if (!fam_branch_i->mdat || !fam_branch_i->rdat || !fam_branch_i->k2dat)
                    XLAL_ERROR_NULL(XLAL_ENOMEM);
            if(min_fam==0){
                if (!fam_branch_i->mbdat || !fam_branch_i->k3dat || !fam_branch_i->k4dat)
                    XLAL_ERROR_NULL(XLAL_ENOMEM);
            }

            // Append the data into the branch
            for (int j = 0 ; j < ndat_branch[b] ; j++){
                fam_branch_i->pdat[j] = pdat[index_begin_stable_branch[b] + j];
                fam_branch_i->mdat[j] = mdat[index_begin_stable_branch[b] + j];
                fam_branch_i->rdat[j] = rdat[index_begin_stable_branch[b] + j];
                fam_branch_i->k2dat[j] = k2dat[index_begin_stable_branch[b] + j];
                if (min_fam == 0) {
                    fam_branch_i->mbdat[j] = mbdat[index_begin_stable_branch[b] + j];
                    fam_branch_i->k3dat[j] = k3dat[index_begin_stable_branch[b] + j];
                    fam_branch_i->k4dat[j] = k4dat[index_begin_stable_branch[b] + j];
                }
            }

            // Recalculate the maximum mass of a branch
            if (index_end_stable_branch[b] < ndat-1) { // if an unstable branch follows the current stable branch
                double pc_max = get_central_pressure_mext(eos, index_end_stable_branch[b], pdat, mdat, -1.0);
                if(min_fam==1) XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(
                                    &fam_branch_i->rdat[ndat_branch[b]-1],
                                    &fam_branch_i->mdat[ndat_branch[b]-1],
                                    &fam_branch_i->k2dat[ndat_branch[b]-1],
                                    pc_max, eos, 1e-6);
                else XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
                        &fam_branch_i->rdat[ndat_branch[b]-1],
                        &fam_branch_i->mdat[ndat_branch[b]-1],
                        &fam_branch_i->mbdat[ndat_branch[b]-1],
                        &fam_branch_i->k2dat[ndat_branch[b]-1],
                        &fam_branch_i->k3dat[ndat_branch[b]-1],
                        &fam_branch_i->k4dat[ndat_branch[b]-1],
                        pc_max, eos, 1e-6);
                fam_branch_i->pdat[ndat_branch[b]-1] = pc_max;

                if (b == nb_stable_branches -1 ) fam->mtov = 1; // if the last stable branch is followed by an unstable branch
            }


            // Recalculate the minimum mass of a branch
            if (index_begin_stable_branch[b] != 0) {
                double pc_min = get_central_pressure_mext(eos, index_begin_stable_branch[b], pdat, mdat, 1.0);
                if(min_fam==1){
                    XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(
                        &fam_branch_i->rdat[0],
                        &fam_branch_i->mdat[0],
                        &fam_branch_i->k2dat[0],
                        pc_min, eos, 1e-6);

                } else {
                    XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(
                        &fam_branch_i->rdat[0],
                        &fam_branch_i->mdat[0],
                        &fam_branch_i->mbdat[0],
                        &fam_branch_i->k2dat[0],
                        &fam_branch_i->k3dat[0],
                        &fam_branch_i->k4dat[0],
                        pc_min, eos, 1e-6);
                }
                fam_branch_i->pdat[0] = pc_min;
            }

            if (ndat_branch[b] > min_interp_points) {
                // Set up p(m), r(m), and k(m) interpolators
                fam_branch_i->m_of_p_acc = gsl_interp_accel_alloc();
                fam_branch_i->p_of_m_acc = gsl_interp_accel_alloc();
                fam_branch_i->r_of_m_acc = gsl_interp_accel_alloc();
                fam_branch_i->k2_of_m_acc = gsl_interp_accel_alloc();

                fam_branch_i->m_of_p_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);
                fam_branch_i->p_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);
                fam_branch_i->r_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);
                fam_branch_i->k2_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);

                gsl_interp_init(fam_branch_i->m_of_p_interp, fam_branch_i->pdat, fam_branch_i->mdat, ndat_branch[b]);
                gsl_interp_init(fam_branch_i->p_of_m_interp, fam_branch_i->mdat, fam_branch_i->pdat, ndat_branch[b]);
                gsl_interp_init(fam_branch_i->r_of_m_interp, fam_branch_i->mdat, fam_branch_i->rdat, ndat_branch[b]);
                gsl_interp_init(fam_branch_i->k2_of_m_interp, fam_branch_i->mdat, fam_branch_i->k2dat, ndat_branch[b]);

                if(min_fam==0){
                    fam_branch_i->mbdat = LALRealloc(fam_branch_i->mbdat, ndat_branch[b] * sizeof(*fam_branch_i->mbdat));
                    fam_branch_i->k3dat = LALRealloc(fam_branch_i->k3dat, ndat_branch[b] * sizeof(*fam_branch_i->k3dat));
                    fam_branch_i->k4dat = LALRealloc(fam_branch_i->k4dat, ndat_branch[b] * sizeof(*fam_branch_i->k4dat));

                    fam_branch_i->mb_of_m_acc = gsl_interp_accel_alloc();
                    fam_branch_i->k3_of_m_acc = gsl_interp_accel_alloc();
                    fam_branch_i->k4_of_m_acc = gsl_interp_accel_alloc();

                    fam_branch_i->mb_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);
                    fam_branch_i->k3_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);
                    fam_branch_i->k4_of_m_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat_branch[b]);

                    gsl_interp_init(fam_branch_i->mb_of_m_interp, fam_branch_i->mdat, fam_branch_i->mbdat, ndat_branch[b]);
                    gsl_interp_init(fam_branch_i->k3_of_m_interp, fam_branch_i->mdat, fam_branch_i->k3dat, ndat_branch[b]);
                    gsl_interp_init(fam_branch_i->k4_of_m_interp, fam_branch_i->mdat, fam_branch_i->k4dat, ndat_branch[b]);
                }
            }
        fam->fam_branch[b] = fam_branch_i;
        }
        LALFree(ndat_branch);
    } else printf("No stable branch has been found in the TOV solution\n");

    LALFree(pdat);
    LALFree(mdat);
    LALFree(rdat);
    LALFree(mbdat);
    LALFree(k2dat);
    LALFree(k3dat);
    LALFree(k4dat);

    return fam;
}
/**
 * @brief Returns the number of branches of a neutron star family.
 * @param fam Pointer to the neutron star family structure.
 * @return The number of branches (stable hydrostatic configurations) for a given neutron star family.
 */
int XLALSimNeutronStarFamNumberOfBranches(FamMultiParts *fam)
{
    return fam->number_of_branches;
}


/**
 * @brief Returns the mass for the last stable configuration Mtov of a neutron star family.
 * @param fam Pointer to the neutron star family structure.
 * @return The TOV mass limit of a neutron star family, zero if none was found.
 */
double XLALSimNeutronStarFamMassTOVLimit(FamMultiParts *fam)
{
    double mtov;
    double mmax_branch = 0.;
    double mmax_branch_temp = 0.;
    if (fam->mtov == 1) {
        int nb_branch = XLALSimNeutronStarFamNumberOfBranches(fam);
        for (int b = 0 ; b < nb_branch; b++){
            mmax_branch_temp = XLALSimNeutronStarFamBranchMaxMass(b, fam)
            if (mmax_branch_temp > mmax_branch)
                mmax_branch = mmax_branch_temp;
        }
        mtov = mmax_branch;
    } else mtov = 0.;
    return mtov;
}

/**
 * @brief Returns the minimum mass of the ith branch of a neutron star family.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The minimum mass of the ith branch for a given neutron star family.
 */
double XLALSimNeutronStarFamBranchMinMass(int branch, FamMultiParts *fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->mdat[0];
}


/**
 * @brief Returns the maximum mass for a given branch of a neutron star family.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The maximum mass for a given branch of a neutron star family.
 */
double XLALSimNeutronStarFamBranchMaxMass(int branch, FamMultiParts *fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->mdat[b->ndat - 1];
}


/**
 * @brief Returns the minimum central central pressure for a given branch of a neutron star family.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The minimum central pressure for a given branch of a neutron star family.
 */
double XLALSimNeutronStarFamBranchMinCentralPressure(int branch, FamMultiParts *fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->pdat[0];
}

/**
 * @brief Returns the maximum central central pressure for a given branch of a neutron star family.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The maximum central pressure for a given branch of a neutron star family.
 */
double XLALSimNeutronStarFamBranchMaxCentralPressure(int branch, FamMultiParts *fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];

    return b->pdat[b->ndat - 1];
}


/**
 * @brief Returns the radius corresponding to a given mass in the neutron star family.
 * @param m Mass at which the radius is interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The radius interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchRadius(double m, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double r;
    r = gsl_interp_eval(b->r_of_m_interp, b->mdat, b->rdat, m,
        b->r_of_m_acc);
    return r;
}



/**
 * @brief Returns the central pressure corresponding to a given mass in the neutron star family.
 * @param m Double for the mass at which the central pressure is interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The central pressure interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchCentralPressure(double m, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double p;
    p = gsl_interp_eval(b->p_of_m_interp, b->mdat, b->pdat, m,
        b->p_of_m_acc);
    return p;
}



/**
 * @brief Returns the mass corresponding to a given central pressure in the neutron star family.
 * @param p Double for the central pressure at which the mass if interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The mass interpolated at central pressure p, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchMass(double p, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double m;
    m = gsl_interp_eval(b->m_of_p_interp, b->pdat, b->mdat, p,
        b->m_of_p_acc);
    return m;
}


/**
 * @brief Returns the baryonic mass corresponding to a given (gravitationnal) mass in the neutron star family.
 * @param m Double for the mass at which the baryonic mass is interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The baryonic mass interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchBaryonicMass(double m, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    /* In case user used mini solver */
    if (b == NULL || b->mb_of_m_interp == NULL)
        XLAL_ERROR_REAL8(XLAL_EINVAL,"branch does not exist or mb never evaluated");
    double mb;
    mb = gsl_interp_eval(b->mb_of_m_interp, b->mdat, b->mbdat, m,
        b->mb_of_m_acc);
    return mb;
}


/**
 * @brief Returns k2 Love number corresponding to a given mass in the neutron star family.
 * @param m Double for the mass at which the baryonic mass is interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The k2 Love number interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchLoveNumberK2(double m, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    double k2;
    k2 = gsl_interp_eval(b->k2_of_m_interp, b->mdat, b->k2dat, m,
        b->k2_of_m_acc);
    return k2;
}


/**
 * @brief Returns k3 Love number corresponding to a given mass in the neutron star family.
 * @param m Double for the mass at which the baryonic mass is interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The k3 Love number interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchLoveNumberK3(double m, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    /* In case user used mini solver */
    if (b == NULL || b->k3_of_m_interp == NULL)
	XLAL_ERROR_REAL8(XLAL_EINVAL,"branch does not exist or k3 never evaluated");
    double k3;
    k3 = gsl_interp_eval(b->k3_of_m_interp, b->mdat, b->k3dat, m,
        b->k3_of_m_acc);
    return k3;
}


/**
 * @brief Returns k4 Love number corresponding to a given mass in the neutron star family.
 * @param m Double for the mass at which the baryonic mass is interpolated.
 * @param branch Integer for the branch number.
 * @param fam Pointer to the neutron star family structure.
 * @return The k4 Love number interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBranchLoveNumberK4(double m, int branch, FamMultiParts * fam)
{
    if (branch < 0 || branch >= fam->number_of_branches) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarFamily * b = fam->fam_branch[branch];
    /* In case user used mini solver */
    if (b == NULL || b->k4_of_m_interp == NULL)
        XLAL_ERROR_REAL8(XLAL_EINVAL,"branch does not exist or k4 never evaluated");
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
