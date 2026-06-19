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
#include <lal/AVFactories.h>
#include <lal/LALSimNeutronStar.h>

/** @cond */

/* Contents of the neutron star family structure. */
typedef struct tagFamBranch FamBranch;

struct tagFamBranch {
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

struct tagLALSimNeutronStarFamily{
  int number_of_branches;
  int mtov;
  double pcmin;
  double pcmax;
  FamBranch ** fam_branch;
};

/* Function used by GSLminimizer to find the central pressure for the
 * exact maximum or minimum neutron star mass from LALSimNeutronStarEOS structure.
 */
static double fextrimizer(double x, void * params)
{
    LALSimNeutronStarEOS * eos = params;
    double r, m, k2;
    XLALSimNeutronStarTOVODEIntegrateWithTolerance(
        &r, &m, &k2, x, eos, 1e-6);
    return m;
}

static double get_central_pressure_mext(LALSimNeutronStarEOS *eos, int index, double *pdat, double *mdat, int factor){
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

/*
 * Frees the memory associated with a pointer to a neutron star branch.
 */
static void XLALDestroySimNeutronStarBranch(FamBranch * branch)
{
    if (branch) {
        if (branch->k2_of_m_acc) gsl_interp_accel_free(branch->k2_of_m_acc);
        if (branch->k3_of_m_acc) gsl_interp_accel_free(branch->k3_of_m_acc);
        if (branch->k4_of_m_acc) gsl_interp_accel_free(branch->k4_of_m_acc);
        if (branch->r_of_m_acc) gsl_interp_accel_free(branch->r_of_m_acc);
        if (branch->p_of_m_acc) gsl_interp_accel_free(branch->p_of_m_acc);
        if (branch->mb_of_m_acc) gsl_interp_accel_free(branch->mb_of_m_acc);
        if (branch->m_of_p_acc) gsl_interp_accel_free(branch->m_of_p_acc);

        if (branch->k2_of_m_interp) gsl_interp_free(branch->k2_of_m_interp);
        if (branch->k3_of_m_interp) gsl_interp_free(branch->k3_of_m_interp);
        if (branch->k4_of_m_interp) gsl_interp_free(branch->k4_of_m_interp);
        if (branch->r_of_m_interp) gsl_interp_free(branch->r_of_m_interp);
        if (branch->p_of_m_interp) gsl_interp_free(branch->p_of_m_interp);
        if (branch->mb_of_m_interp) gsl_interp_free(branch->mb_of_m_interp);
        if (branch->m_of_p_interp) gsl_interp_free(branch->m_of_p_interp);

        LALFree(branch->k2dat);
        LALFree(branch->k3dat);
        LALFree(branch->k4dat);
        LALFree(branch->rdat);
        LALFree(branch->mdat);
        LALFree(branch->mbdat);
        LALFree(branch->pdat);

        LALFree(branch);
    }
    return ;
}

/**
 * @brief Frees the memory associated with a pointer to a neutron star family.
 * @param fam Pointer to the neutron star family structure to be freed.
 */
void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam)
{
    if (fam) {
        if (fam->fam_branch) {
            for (int i = 0; i < fam->number_of_branches; ++i)
                XLALDestroySimNeutronStarBranch(fam->fam_branch[i]);
            LALFree(fam->fam_branch);
        }
        LALFree(fam);
    }
    return;

}


/**
 * @brief Creates a neutron star family structure (LALSimNeutronStarFamily) which can accomodate
 * multiple branches for a given LALSimNeutronStarEOS equation of state structure and minimum
 * central pressure in Pa.
 * @details
 * The neutron star family contains the astrophysical parameter solution for a
 * given equation of state. When the structure is created, it solves the TOV+Love number
 * ODEs with XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance if min_fam = 0, or with
 * XLALSimNeutronStarTOVODEMiniIntegrate if min_fam=1. The minimum central pressure
 * is an input from the user.
 * The maximum central pressure is that defined as the largest pressure of the equation of state.
 * The neutron star family can accomodate multiple branches and provides interpolating
 * function for astrophysical quantities only for stable branches (continuously increasing
 * mass in a sequence).
 * @param eos Pointer to the Equation of State structure.
 * @param min_fam Flag to choose to construct the neutron star family
 * with only the minimum number of neutron star astrophysical parameters
 * (0: all NS parameters, 1: only mass, radius and k2 Love number).
 * @param logPcmin Logarithm of the minimum central pressure in Pa
 * @return A pointer to the neutron star family structure.
 */
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamilyWithPcmin(LALSimNeutronStarEOS * eos, int min_fam, double logPcmin){

    if(min_fam!=0 && min_fam!=1) XLAL_ERROR_NULL(XLAL_EDOM);

    LALSimNeutronStarFamily *fam;
    fam = LALCalloc(1, sizeof(*fam));

    int ndat = 100;
    int min_interp_points = 5; // Minimum number of points required for all gls interpolators
    double logpcmin = logPcmin;
    fam->pcmin = exp(logpcmin);
    fam->pcmax = XLALSimNeutronStarEOSMaxPressure(eos);
    double logpcmax = log(fam->pcmax);
    double dlogp = (logpcmax - logpcmin) / (ndat-1); // in Pascal
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
    	if(j==0) pdat[j] = fam->pcmin; // exact min
	else if(j==ndat-1) pdat[j] = fam->pcmax; // exact max
	else pdat[j] = exp(logpcmin + j * dlogp); // pressure in Pascal
        // Solve TOV+Love number ODEs
        if(min_fam==1){
            XLALSimNeutronStarTOVODEIntegrateWithTolerance(
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
        fam->fam_branch = (FamBranch **) LALMalloc(sizeof(FamBranch *) * nb_stable_branches);
        int *ndat_branch;
        ndat_branch = LALMalloc(nb_stable_branches * sizeof(*ndat_branch));
        for (int b = 0; b < nb_stable_branches; b++){
            ndat_branch[b] = index_end_stable_branch[b] - index_begin_stable_branch[b] + 1;
            FamBranch * fam_branch_i = LALCalloc(1, sizeof(FamBranch));
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

            /* On allocation failure, free everything before returning. */
            if (!fam_branch_i->pdat || !fam_branch_i->mdat || !fam_branch_i->rdat || !fam_branch_i->k2dat ||
                (min_fam==0 && (!fam_branch_i->mbdat || !fam_branch_i->k3dat || !fam_branch_i->k4dat))) {
                XLALDestroySimNeutronStarBranch(fam_branch_i);
                LALFree(ndat_branch);
                LALFree(fam->fam_branch);
                LALFree(fam);
                LALFree(pdat);
                LALFree(mdat);
                LALFree(rdat);
                LALFree(mbdat);
                LALFree(k2dat);
                LALFree(k3dat);
                LALFree(k4dat);
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
                double pc_max = get_central_pressure_mext(eos, index_end_stable_branch[b], pdat, mdat, -1);
                if(min_fam==1) XLALSimNeutronStarTOVODEIntegrateWithTolerance(
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
                double pc_min = get_central_pressure_mext(eos, index_begin_stable_branch[b], pdat, mdat, 1);
                if(min_fam==1){
                    XLALSimNeutronStarTOVODEIntegrateWithTolerance(
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
 * @brief Creates a neutron star family structure (LALSimNeutronStarFamily) which can accomodate
 * multiple branches for a given LALSimNeutronStarEOS equation of state structure.
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
 * @param min_fam Flag to choose to construct the neutron star family
 * with only the minimum number of neutron star astrophysical parameters
 * (0: all NS parameters, 1: only mass, radius and k2 Love number).
 * @return A pointer to the neutron star family structure.
 */
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(LALSimNeutronStarEOS * eos, int min_fam){

    double logpcmin = 75.5; // pcmin = 6.16e32
    if (XLALSimNeutronStarEOSMaxPressure(eos) < exp(logpcmin)) XLAL_ERROR_NULL(XLAL_EDOM);
    LALSimNeutronStarFamily * fam = XLALCreateSimNeutronStarFamilyWithPcmin(eos, min_fam, logpcmin);

    return fam;
}
/**
 * @brief Returns the number of branches of a neutron star family.
 * @param fam Pointer to the neutron star family structure.
 * @return The number of branches (stable hydrostatic configurations) for a given neutron star family.
 */
int XLALSimNeutronStarFamNumberOfBranches(LALSimNeutronStarFamily *fam)
{
    return fam->number_of_branches;
}

/*
 * Returns the ith branch of a neutron star family.
*/
static FamBranch * XLALSimNeutronStarFamSelectBranch(LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_NULL(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    return fam->fam_branch[branch_id];
}


static double find_max_mass_for_multiple_branches(LALSimNeutronStarFamily *fam){
    double mmax_branch = 0.;
    double mmax_branch_temp = 0.;
    int nb_branch = XLALSimNeutronStarFamNumberOfBranches(fam);
    for (int b = 0 ; b < nb_branch; b++){
        mmax_branch_temp = XLALSimNeutronStarFamMaxMassPerBranch(fam, b);
        if (mmax_branch_temp > mmax_branch)
            mmax_branch = mmax_branch_temp;
    }
    return mmax_branch;
}

static double find_min_mass_for_multiple_branches(LALSimNeutronStarFamily *fam){
    double mmin_branch = 1e100;
    double mmin_branch_temp = 1e100;
    int nb_branch = XLALSimNeutronStarFamNumberOfBranches(fam);
    for (int b = 0; b < nb_branch; b++){
        mmin_branch_temp = XLALSimNeutronStarFamMinMassPerBranch(fam, b);
        if (mmin_branch_temp < mmin_branch)
            mmin_branch = mmin_branch_temp;
    }
    return mmin_branch;
}

/**
 * @brief Returns the mass for the last stable configuration Mtov of a neutron star family in kg.
 * @param fam Pointer to the neutron star family structure.
 * @return The TOV mass limit in kg of a neutron star family, zero if none was found.
 */
double XLALSimNeutronStarFamMassTOVLimit(LALSimNeutronStarFamily *fam)
{
    double mtov;
    if (fam->mtov == 1)
        mtov = find_max_mass_for_multiple_branches(fam);
    else mtov = 0.;
    return mtov;
}

/*
 * Returns the minimum mass of a neutron star branch in kg.
 */
static double XLALSimNeutronStarBranchMinMass(FamBranch * branch)
{
    return branch->mdat[0];
}

/**
 * @brief Returns the minimum mass of the ith branch of a neutron star family in kg.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The minimum mass in kg of the ith branch for a given neutron star family .
 */
double XLALSimNeutronStarFamMinMassPerBranch(LALSimNeutronStarFamily *fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchMinMass(b);
}


/**
 * @brief Returns the minimum mass of a neutron star family in kg.
 * @param fam Pointer to the neutron star family structure.
 * @return The minimum mass in kg of a neutron star family.
 */
double XLALSimNeutronStarFamMinMass(LALSimNeutronStarFamily *fam)
{
    return find_min_mass_for_multiple_branches(fam);
}


/*
 * Returns the maximum mass of a neutron star branch in kg.
 */
static double XLALSimNeutronStarBranchMaxMass(FamBranch * branch)
{
    return branch->mdat[branch->ndat - 1];
}

/**
 * @brief Returns the maximum mass for the ith branch of a neutron star family in kg.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The maximum mass in kg for the ith branch of a neutron star family.
 */
double XLALSimNeutronStarFamMaxMassPerBranch(LALSimNeutronStarFamily *fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchMaxMass(b);
}


/**
 * @brief Returns the maximum mass of a neutron star family in kg.
 * @param fam Pointer to the neutron star family structure.
 * @return The maximum mass in kg of a neutron star family.
 */
double XLALSimNeutronStarFamMaxMass(LALSimNeutronStarFamily *fam)
{
    return find_max_mass_for_multiple_branches(fam);
}

/*
 * Returns the minimum central central pressure of a neutron star branch in Pa.
 */
static double XLALSimNeutronStarBranchMinCentralPressure(FamBranch * branch)
{
    return branch->pdat[0];
}


/**
 * @brief Returns the minimum central central pressure for the ith branch of a neutron star family in Pa.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The minimum central pressure in Pa for the ith branch of a neutron star family.
 */
double XLALSimNeutronStarFamMinCentralPressurePerBranch(LALSimNeutronStarFamily *fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchMinCentralPressure(b);
}

/**
 * @brief Returns the minimum central central pressure of a neutron star family in Pa.
 * @param fam Pointer to the neutron star family structure.
 * @return The minimum central pressure in Pa of a neutron star family.
 */
double XLALSimNeutronStarFamMinCentralPressure(LALSimNeutronStarFamily *fam)
{
    return XLALSimNeutronStarFamMinCentralPressurePerBranch(fam, 0);
}


/*
 * Returns the maximum central central pressure of a neutron star branch in Pa.
 */
static double XLALSimNeutronStarBranchMaxCentralPressure(FamBranch * branch)
{
    return branch->pdat[branch->ndat - 1];
}

/**
 * @brief Returns the maximum central central pressure for a given branch of a neutron star family in Pa.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The maximum central pressure in Pa for a given branch of a neutron star family.
 */
double XLALSimNeutronStarFamMaxCentralPressurePerBranch(LALSimNeutronStarFamily *fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchMaxCentralPressure(b);
}

/**
 * @brief Returns the maximum central central pressure of a neutron star family in Pa.
 * @param fam Pointer to the neutron star family structure.
 * @return The maximum central pressure in Pa of a neutron star family.
 */
double XLALSimNeutronStarFamMaxCentralPressure(LALSimNeutronStarFamily *fam)
{
    int branch_number = XLALSimNeutronStarFamNumberOfBranches(fam) - 1;
    return XLALSimNeutronStarFamMaxCentralPressurePerBranch(fam, branch_number);
}


/* FUNCTIONS RELATED TO INTERPOLATION IN NS FAMILY */

/*
 * Returns the radius in m corresponding to a given mass in kg in the neutron star branch.
 */
static double XLALSimNeutronStarBranchRadiusOfMass(double m, FamBranch * branch)
{
    if (!branch || !branch->r_of_m_interp || !branch->r_of_m_acc || !branch->mdat || !branch->rdat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double r;
    r = gsl_interp_eval(branch->r_of_m_interp, branch->mdat, branch->rdat, m,
        branch->r_of_m_acc);
    return r;
}

/**
 * @brief Returns the radius in m corresponding to a given mass in kg for the ith branch of the neutron star family.
 * @param m Mass in kg at which the radius is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The radius in m interpolated at mass m in kg, for the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamRadiusOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double mmin = XLALSimNeutronStarFamMinMassPerBranch(fam, branch_id);
    double mmax = XLALSimNeutronStarFamMaxMassPerBranch(fam, branch_id);
    if (m < mmin || m > mmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input mass m = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", m, mmin, mmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchRadiusOfMass(m, b);
}

/* Function to find the number of twins given a mass value */
static int find_number_of_twins_at_mass(double m, LALSimNeutronStarFamily * fam){
    int twins = 0;
    for (int i = 0; i < XLALSimNeutronStarFamNumberOfBranches(fam); i++){
        if (m >= XLALSimNeutronStarFamMinMassPerBranch(fam, i) && m <= XLALSimNeutronStarFamMaxMassPerBranch(fam, i))
            twins +=1;
    }
    return twins;
}

static int * find_branch_twins_at_mass(double m, LALSimNeutronStarFamily * fam, int twins){
    int * branch_indices = XLALCalloc(twins, sizeof(*branch_indices));
    int k = 0;
    for (int i = 0; i < XLALSimNeutronStarFamNumberOfBranches(fam); i++)
        if (m >= XLALSimNeutronStarFamMinMassPerBranch(fam, i) && m <= XLALSimNeutronStarFamMaxMassPerBranch(fam, i)){
            branch_indices[k++] = i;
    }
    return branch_indices;
}

/**
 * @brief Returns the radius in m corresponding to a given mass in kg.
 * @details This function can accomodate twins stars, it returns a REAL8Vector
 * whose length corresponds to the number of twins in the sequence. Requires
 * deallocation of memory with XLALDestroyREAL8Vector.
 * @param m Mass in kg at which the radius is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The radius in m interpolated at mass m in kg.
 */
REAL8Vector * XLALSimNeutronStarFamRadiusOfMass(double m, LALSimNeutronStarFamily * fam){
    int twins = find_number_of_twins_at_mass(m, fam);
    REAL8Vector * r = XLALCreateREAL8Vector(twins);
    int * twin_indices = find_branch_twins_at_mass(m, fam, twins);
    for (int i = 0; i < twins; i++)
        r->data[i] = XLALSimNeutronStarFamRadiusOfMassPerBranch(m, fam, twin_indices[i]);
    LALFree(twin_indices);
    return r;
}

/*
 * Returns the central pressure in Pa corresponding to a given mass in kg in the neutron star branch.
 */
static double XLALSimNeutronStarBranchCentralPressureOfMass(double m, FamBranch * branch)
{
    if (!branch || !branch->p_of_m_interp || !branch->p_of_m_acc || !branch->mdat || !branch->pdat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double p;
    p = gsl_interp_eval(branch->p_of_m_interp, branch->mdat, branch->pdat, m,
        branch->p_of_m_acc);
    return p;
}


/**
 * @brief Returns the central pressure in Pa corresponding to a given mass in kg
 * for the ith branch of the neutron star family.
 * @param m Mass in kg at which the central pressure is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The central pressure in Pa interpolated at mass m, for the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamCentralPressureOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double mmin = XLALSimNeutronStarFamMinMassPerBranch(fam, branch_id);
    double mmax = XLALSimNeutronStarFamMaxMassPerBranch(fam, branch_id);
    if (m < mmin || m > mmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input mass m = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", m, mmin, mmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchCentralPressureOfMass(m, b);
}

/**
 * @brief Returns the central pressure in Pa corresponding to a given mass in kg.
 * @details This function can accomodate twins stars, it returns a REAL8Vector
 * whose length corresponds to the number of twins in the sequence. Requires deallocation
 * of memory with XLALDestroyREAL8Vector.
 * @param m Mass in kg at which the central pressure is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The central pressure in Pa interpolated at mass m in kg.
 */
REAL8Vector * XLALSimNeutronStarFamCentralPressureOfMass(double m, LALSimNeutronStarFamily * fam){
    int twins = find_number_of_twins_at_mass(m, fam);
    REAL8Vector * pc = XLALCreateREAL8Vector(twins);
    int * twin_indices = find_branch_twins_at_mass(m, fam, twins);
    for (int i = 0; i < twins; i++)
        pc->data[i] = XLALSimNeutronStarFamCentralPressureOfMassPerBranch(m, fam, twin_indices[i]);
    LALFree(twin_indices);
    return pc;
}

/*
 * Returns the mass in kg corresponding to a given central pressure in Pa in the neutron star branch.
 */
static double XLALSimNeutronStarBranchMassOfCentralPressure(double p, FamBranch * branch)
{
    if (!branch || !branch->m_of_p_interp || !branch->m_of_p_acc || !branch->mdat || !branch->pdat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double m;
    m = gsl_interp_eval(branch->m_of_p_interp, branch->pdat, branch->mdat, p,
        branch->m_of_p_acc);
    return m;
}

/**
 * @brief Returns the mass in kg corresponding to a given central pressure in Pa
 * for the ith branch of the neutron star family.
 * @param p Central pressure in Pa at which the mass is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The mass in kg interpolated at central pressure p, for the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamMassOfCentralPressurePerBranch(double p, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double pcmin = XLALSimNeutronStarFamMinCentralPressurePerBranch(fam, branch_id);
    double pcmax = XLALSimNeutronStarFamMaxCentralPressurePerBranch(fam, branch_id);
    if (p < pcmin || p > pcmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input central pressure pc = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", p, pcmin, pcmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchMassOfCentralPressure(p, b);
}

/* This function looks for the branch index where a central pressure p exists */
static int find_branch_of_pressure(double p, LALSimNeutronStarFamily * fam){
    for (int i = 0; i < XLALSimNeutronStarFamNumberOfBranches(fam); i++){
        if (p >= XLALSimNeutronStarFamMinCentralPressurePerBranch(fam, i) &&
            p <= XLALSimNeutronStarFamMaxCentralPressurePerBranch(fam, i))
            return i;
    }
    return -1;
}
/**
 * @brief Returns the mass in kg corresponding to a given central pressure in Pa.
 * @param p Central pressure in Pa at which the mass is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The mass in kg interpolated at central pressure in Pa.
 */
double XLALSimNeutronStarFamMassOfCentralPressure(double p, LALSimNeutronStarFamily * fam){
    int index_branch = find_branch_of_pressure(p, fam);
    if (index_branch < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "Central pressure p = %.16e is not in any branch.", p);
    return XLALSimNeutronStarFamMassOfCentralPressurePerBranch(p, fam, index_branch);
}

/*
 * Returns the baryonic mass in kg corresponding to a given
 * (gravitationnal) mass in kg in the neutron star branch.
 */
static double XLALSimNeutronStarBranchBaryonicMassOfMass(double m, FamBranch * branch)
{
    /* In case user used mini solver */
    if (!branch || !branch->mb_of_m_interp || !branch->mb_of_m_acc || !branch->mdat || !branch->mbdat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double mb;
    mb = gsl_interp_eval(branch->mb_of_m_interp, branch->mdat, branch->mbdat, m,
        branch->mb_of_m_acc);
    return mb;
}

/**
 * @brief Returns the baryonic mass in kg corresponding to a given
 * (gravitationnal) mass in kg for the ith branch of the neutron star family.
 * @param m (Gravitational) Mass in kg at which the baryonic mass is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The baryonic mass in kg interpolated at mass m, in the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamBaryonicMassOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double mmin = XLALSimNeutronStarFamMinMassPerBranch(fam, branch_id);
    double mmax = XLALSimNeutronStarFamMaxMassPerBranch(fam, branch_id);
    if (m < mmin || m > mmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input mass m = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", m, mmin, mmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchBaryonicMassOfMass(m, b);
}

/**
 * @brief Returns the baryonic mass in kg corresponding to a given
 * (gravitationnal) mass in kg.
 * @details This function can accomodate twins stars, it returns a REAL8Vector
 * whose length corresponds to the number of twins in the sequence. Requires deallocation
 * of memory with XLALDestroyREAL8Vector.
 * @param m Mass in kg at which the baryonic mass is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The baryonic mass in kg interpolated at mass m in kg.
 */
REAL8Vector * XLALSimNeutronStarFamBaryonicMassOfMass(double m, LALSimNeutronStarFamily * fam){
    int twins = find_number_of_twins_at_mass(m, fam);
    REAL8Vector * mb = XLALCreateREAL8Vector(twins);
    int * twin_indices = find_branch_twins_at_mass(m, fam, twins);
    for (int i = 0; i < twins; i++)
        mb->data[i] = XLALSimNeutronStarFamBaryonicMassOfMassPerBranch(m, fam, twin_indices[i]);
    LALFree(twin_indices);
    return mb;
}

/*
 * Returns k2 Love number (dimensionless) corresponding to a
 * given mass in kg in the neutron star branch.
*/
static double XLALSimNeutronStarBranchLoveNumberK2OfMass(double m, FamBranch * branch)
{
    if (!branch || !branch->k2_of_m_interp || !branch->k2_of_m_acc || !branch->mdat || !branch->k2dat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double k2;
    k2 = gsl_interp_eval(branch->k2_of_m_interp, branch->mdat, branch->k2dat, m,
        branch->k2_of_m_acc);
    return k2;
}

/**
 * @brief Returns k2 Love number (dimensionless) corresponding to a given mass in kg
 * for the ith branch of the neutron star family.
 * @param m Mass in kg at which the k2 Love number is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The k2 Love number (dimensionless) interpolated at mass m for the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamLoveNumberK2OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double mmin = XLALSimNeutronStarFamMinMassPerBranch(fam, branch_id);
    double mmax = XLALSimNeutronStarFamMaxMassPerBranch(fam, branch_id);
    if (m < mmin || m > mmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input mass m = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", m, mmin, mmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchLoveNumberK2OfMass(m, b);
}

/**
 * @brief Returns k2 Love number (dimensionless) corresponding to a given mass in kg.
 * @details This function can accomodate twins stars, it returns a REAL8Vector
 * whose length corresponds to the number of twins in the sequence. Requires deallocation
 * of memory with XLALDestroyREAL8Vector.
 * @param m Mass in kg at which the k2 Love number is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The k2 Love number (dimensionless) interpolated at mass m.
 */
REAL8Vector * XLALSimNeutronStarFamLoveNumberK2OfMass(double m, LALSimNeutronStarFamily * fam){
    int twins = find_number_of_twins_at_mass(m, fam);
    REAL8Vector * k2 = XLALCreateREAL8Vector(twins);
    int * twin_indices = find_branch_twins_at_mass(m, fam, twins);
    for (int i = 0; i < twins; i++)
        k2->data[i] = XLALSimNeutronStarFamLoveNumberK2OfMassPerBranch(m, fam, twin_indices[i]);
    LALFree(twin_indices);
    return k2;
}

/*
 * Returns k3 Love number (dimensionless) corresponding to a given mass
 * in kg in the neutron star branch.
 */
static double XLALSimNeutronStarBranchLoveNumberK3OfMass(double m, FamBranch * branch)
{
    /* In case user used mini solver */
    if (!branch || !branch->k3_of_m_interp || !branch->k3_of_m_acc || !branch->mdat || !branch->k3dat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double k3;
    k3 = gsl_interp_eval(branch->k3_of_m_interp, branch->mdat, branch->k3dat, m,
        branch->k3_of_m_acc);
    return k3;
}


/**
 * @brief Returns k3 Love number (dimensionless) corresponding to a given mass in kg
 * for the ith branch of the neutron star family.
 * @param m Mass in kg at which the k3 Love number is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The k3 Love number (dimensionless) interpolated at mass m for the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamLoveNumberK3OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double mmin = XLALSimNeutronStarFamMinMassPerBranch(fam, branch_id);
    double mmax = XLALSimNeutronStarFamMaxMassPerBranch(fam, branch_id);
    if (m < mmin || m > mmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input mass m = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", m, mmin, mmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchLoveNumberK3OfMass(m, b);
}

/**
 * @brief Returns k3 Love number (dimensionless) corresponding to a given mass in kg.
 * @details This function can accomodate twins stars, it returns a REAL8Vector
 * whose length corresponds to the number of twins in the sequence. Requires deallocation
 * of memory with XLALDestroyREAL8Vector.
 * @param m Mass in kg at which the k3 Love number is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The k3 Love number (dimensionless) interpolated at mass m.
 */
REAL8Vector * XLALSimNeutronStarFamLoveNumberK3OfMass(double m, LALSimNeutronStarFamily * fam){
    int twins = find_number_of_twins_at_mass(m, fam);
    REAL8Vector * k3 = XLALCreateREAL8Vector(twins);
    int * twin_indices = find_branch_twins_at_mass(m, fam, twins);
    for (int i = 0; i < twins; i++)
        k3->data[i] = XLALSimNeutronStarFamLoveNumberK3OfMassPerBranch(m, fam, twin_indices[i]);
    LALFree(twin_indices);
    return k3;
}

/*
 * Returns k4 Love number (dimensionless) corresponding to a given mass
 * in kg in the neutron star branch.
 */
static double XLALSimNeutronStarBranchLoveNumberK4OfMass(double m, FamBranch * branch)
{
    /* In case user used mini solver */
    if (!branch || !branch->k4_of_m_interp || !branch->k4_of_m_acc || !branch->mdat || !branch->k4dat)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid family branch (uninitialized interpolation)");
    double k4;
    k4 = gsl_interp_eval(branch->k4_of_m_interp, branch->mdat, branch->k4dat, m,
        branch->k4_of_m_acc);
    return k4;
}

/**
 * @brief Returns k4 Love number (dimensionless) corresponding to a given mass
 * in kg for the ith branch of the neutron star family.
 * @param m Mass in kg at which the k4 Love number is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @param branch_id Integer for the branch number.
 * @return The k4 Love number (dimensionless) interpolated at mass m for the ith branch of the neutron star family.
 */
double XLALSimNeutronStarFamLoveNumberK4OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id)
{
    if (branch_id < 0 || branch_id >= fam->number_of_branches)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID branch number of LALSimNeutronStarFamily structure is incorrect.");
    double mmin = XLALSimNeutronStarFamMinMassPerBranch(fam, branch_id);
    double mmax = XLALSimNeutronStarFamMaxMassPerBranch(fam, branch_id);
    if (m < mmin || m > mmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input mass m = %.16e is beyond the Family branch interpolation range [%.16e:%.16e].", m, mmin, mmax);
    FamBranch * b = XLALSimNeutronStarFamSelectBranch(fam, branch_id);
    return XLALSimNeutronStarBranchLoveNumberK4OfMass(m, b);
}

/**
 * @brief Returns k4 Love number (dimensionless) corresponding to a given mass in kg.
 * @details This function can accomodate twins stars, it returns a REAL8Vector
 * whose length corresponds to the number of twins in the sequence. Requires deallocation
 * of memory with XLALDestroyREAL8Vector.
 * @param m Mass in kg at which the k4 Love number is interpolated.
 * @param fam Pointer to the neutron star family structure.
 * @return The k4 Love number (dimensionless) interpolated at mass m.
 */
REAL8Vector * XLALSimNeutronStarFamLoveNumberK4OfMass(double m, LALSimNeutronStarFamily * fam){
    int twins = find_number_of_twins_at_mass(m, fam);
    REAL8Vector * k4 = XLALCreateREAL8Vector(twins);
    int * twin_indices = find_branch_twins_at_mass(m, fam, twins);
    for (int i = 0; i < twins; i++)
        k4->data[i] = XLALSimNeutronStarFamLoveNumberK4OfMassPerBranch(m, fam, twin_indices[i]);
    LALFree(twin_indices);
    return k4;
}

/** @} */
