/*
 * Copyright (C) 2013 J. Creighton, B. Lackey
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
 * @author Jolien Creighton, Benjamin Lackey, Philip Davis, Micaela Oertel, Lami Suleiman, Leslie Wade
 * @addtogroup LALSimNeutronStar_h      Header LALSimNeutronStar.h
 * @ingroup lalsimulation_general
 * @brief Provides routines for neutron star physical parameters.
 * @{
 * @defgroup LALSimNeutronStarEOS_c     Module LALSimNeutronStarEOS.c
 * @defgroup LALSimNeutronStarTOV_c     Module LALSimNeutronStarTOV.c
 * @defgroup LALSimNeutronStarFamily_c  Module LALSimNeutronStarFamily.c
 * @}
 * @{
 */

#ifndef _LALSIMNEUTRONSTAR_H
#define _LALSIMNEUTRONSTAR_H

#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>


/* CONSTANTS */

/** Factor to convert density in kg/m^3 to geometrized units of m^-2. */
#define LAL_G_C2_SI ((LAL_G_SI) / ((double)(LAL_C_SI) * (double)(LAL_C_SI)))

/** Factor to convert pressure in Pa to geometrized units of m^-2. */
#define LAL_G_C4_SI ((LAL_G_C2_SI) / ((double)(LAL_C_SI) * (double)(LAL_C_SI)))

/** Nuclear density in kg m^-3. */
#define LAL_NUCLEAR_DENSITY_SI 2.8e17

/** Nuclear density in geometrized units of m^-2. */
#define LAL_NUCLEAR_DENSITY_GEOM_SI ((LAL_NUCLEAR_DENSITY_SI) * (LAL_G_C2_SI))



/* EOS ROUTINES */

/** Incomplete type for the neutron star Equation of State (EOS). */
typedef struct tagLALSimNeutronStarEOS LALSimNeutronStarEOS;
/** Recognised names of equations of state */
extern const char * const lalSimNeutronStarEOSNames[111];

/** Incomplete type for a neutron star family having a particular EOS. */
typedef struct tagLALSimNeutronStarFamily LALSimNeutronStarFamily;



/* FUNCTIONS FOR STRUCTURE CREATION AND DESTRUCTION */
LALSimNeutronStarEOS * XLALCreateSimNeutronStarEOS(int nb_pieces);
void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos);


/* FUNCTIONS RELATED TO AN EQUATION OF STATE PIECE */
int XLALSimNeutronStarEOSNumberPieces(LALSimNeutronStarEOS * eos);



/* FUNCTIONS FOR KEY EQUATION OF STATE VALUES */
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinPressurePerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinPressureGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinPressure(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxPressurePerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxPressureGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinPseudoEnthalpy(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMinEnergyDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinEnergyDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinEnergyDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinEnergyDensity(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMaxEnergyDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxEnergyDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxEnergyDensity(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMinRestMassDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinRestMassDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinRestMassDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinRestMassDensity(LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMaxRestMassDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxRestMassDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxRestMassDensity(LALSimNeutronStarEOS * eos);

/* FUNCTIONS FOR INTERPOLATED EQUATION OF STATE VALUES */

double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrizedPerPiece(double p, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressurePerPiece(double p, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressure(double p, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrizedPerPiece(double p, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSEnergyDensityOfPressurePerPiece(double p, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPressure(double p, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrizedPerPiece(double p, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressurePerPiece(double p, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressure(double p, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrizedPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpy(double h, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrizedPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(double h, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrizedPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(double h, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrizedPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyPerPiece(double h, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpy(double h, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrizedPerPiece(double e,  LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPressureOfEnergyDensityPerPiece(double e, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrized(double e, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfEnergyDensity(double e, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrizedPerPiece(double rho, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPressureOfRestMassDensityPerPiece(double rho, LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrized(double rho, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfRestMassDensity(double rho, LALSimNeutronStarEOS * eos);

// /* DEPRECATED: old misspelling of Geometrized and missing "OfPseudoEnthalpy" for coherence*/
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSound(double h, LALSimNeutronStarEOS * eos);



/* FUNCTIONS TO FILL THE EQUATION OF STATE STRUCTURE */

LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabulatedDataChoiceDirtyPT(
    const REAL8Vector *nbdat, const REAL8Vector *edat, const REAL8Vector *pdat,
    const REAL8Vector *mubdat, const REAL8Vector *muedat, const REAL8Vector *hdat,
    const REAL8Vector *yedat, const REAL8Vector *cs2dat, int dirty);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabulatedData(
    const REAL8Vector *nbdat, const REAL8Vector *edat, const REAL8Vector *pdat,
    const REAL8Vector *mubdat, const REAL8Vector *muedat, const REAL8Vector *hdat,
    const REAL8Vector *yedat, const REAL8Vector *cs2dat);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name);



/* FUNCTIONS FOR EOS PARAMETRISATIONS */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSPolytrope(double Gamma,
    double reference_pressure_si, double reference_density_si);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterPiecewisePolytrope(double
    logp1_si, double gamma1, double gamma2, double gamma3);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSSpectralDecomposition(double gamma[], int size);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterSpectralDecomposition(double SDgamma0, double SDgamma1, double SDgamma2, double SDgamma3);
int XLALSimNeutronStarEOS4ParamSDGammaCheck(double g0, double g1, double g2, double g3);
int XLALSimNeutronStarEOS4ParamSDViableFamilyCheck(double g0, double g1, double g2, double g3);
int XLALSimNeutronStarEOS3PDViableFamilyCheck(double p0, double log10p1_si, double p1, double log10p2_si, double p2, int causal);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSDynamicAnalytic(double parameters[],
    size_t nsec, int causal);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceDynamicPolytrope(double g0,
    double log10p1_si, double g1, double log10p2_si, double g2);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceCausalAnalytic(double v1,
    double log10p1_si, double v2, double log10p2_si, double v3);


/* FUNCTIONS FOR NS ASTROPHYSICAL SOLVER */
/* TOV and Love Number solver ROUTINES */
int XLALSimNeutronStarVirialODEIntegrate(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarVirialODEIntegrateWithTolerance(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

int XLALSimNeutronStarTOVODEExtendedIntegrate(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, LALSimNeutronStarEOS *eos);

int XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, LALSimNeutronStarEOS *eos);

int XLALSimNeutronStarTOVODEIntegrateWithTolerance(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

/* NEUTRON STAR ASTROPHYSICAL PARAMETER'S RELATED ROUTINES */

void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam);

LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamilyWithPcmin(LALSimNeutronStarEOS * eos, int min_fam, double logPcmin);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(LALSimNeutronStarEOS * eos, int min_fam);
int XLALSimNeutronStarFamNumberOfBranches(LALSimNeutronStarFamily *fam);
double XLALSimNeutronStarFamMassTOVLimit(LALSimNeutronStarFamily *fam);

double XLALSimNeutronStarFamMinMassPerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMinMass(LALSimNeutronStarFamily *fam);
double XLALSimNeutronStarFamMaxMassPerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMaxMass(LALSimNeutronStarFamily *fam);
double XLALSimNeutronStarFamMinCentralPressurePerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMinCentralPressure(LALSimNeutronStarFamily *fam);
double XLALSimNeutronStarFamMaxCentralPressurePerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMaxCentralPressure(LALSimNeutronStarFamily *fam);

double XLALSimNeutronStarFamRadiusOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);
REAL8Vector * XLALSimNeutronStarFamRadiusOfMass(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarFamCentralPressureOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);
REAL8Vector * XLALSimNeutronStarFamCentralPressureOfMass(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarFamMassOfCentralPressurePerBranch(double p, LALSimNeutronStarFamily * fam, int branch_id);
double XLALSimNeutronStarFamMassOfCentralPressure(double p, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarFamBaryonicMassOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);
REAL8Vector * XLALSimNeutronStarFamBaryonicMassOfMass(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarFamLoveNumberK2OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);
REAL8Vector * XLALSimNeutronStarFamLoveNumberK2OfMass(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarFamLoveNumberK3OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);
REAL8Vector * XLALSimNeutronStarFamLoveNumberK3OfMass(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarFamLoveNumberK4OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);
REAL8Vector * XLALSimNeutronStarFamLoveNumberK4OfMass(double m, LALSimNeutronStarFamily * fam);

#endif /* _LALSIMNEUTRONSTAR_H */

/** @} */
