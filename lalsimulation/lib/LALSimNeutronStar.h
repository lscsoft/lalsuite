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
typedef struct tagLALSimEOSMultiParts LALSimEOSMultiParts;
/** Recognised names of equations of state */
extern const char * const lalSimNeutronStarEOSNames[111];

/** Incomplete type for a neutron star family having a particular EOS. */
typedef struct tagLALSimNeutronStarFamily LALSimNeutronStarFamily;
typedef struct tagFamMultiParts FamMultiParts;

/* FUNCTIONS FOR KEY EOS VALUES */

void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos);
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPressureGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinPressureGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinPressure(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinEnergyDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinRestMassDensityGeometrized(LALSimNeutronStarEOS * eos);

void XLALDestroySimNeutronStarEOSMultiParts(LALSimEOSMultiParts * eos);
char *XLALSimNeutronStarEOSMultiPartsName(LALSimEOSMultiParts * eos);
int XLALSimNeutronStarEOSMultiPartsNumber(LALSimEOSMultiParts * eos);
LALSimNeutronStarEOS * XLALSimNeutronStarEOSPart(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsMaxPressureGeometrized(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinPressureGeometrized(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPressure(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinPressure(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinPseudoEnthalpy(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPseudoEnthalpy(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinAcausalPseudoEnthalpy(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxEnergyDensityGeometrized(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinEnergyDensityGeometrized(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxRestMassDensityGeometrized(LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinRestMassDensityGeometrized(LALSimEOSMultiParts * eos);

double XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometrized(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxPressure(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinPressure(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxPseudoEnthalpy(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinAcausalPseudoEnthalpy(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxEnergyDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinEnergyDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxRestMassDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinRestMassDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id);

/* FUNCTIONS FOR INTERPOLATED EOS VALUES WITH GEOMETRIZED UNITS */

double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSoundGeometrized(double h,
    LALSimNeutronStarEOS * eos);

/* DEPRECATED: old misspelling of Geometrized */
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h, LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSoundGeometrized(double h,
    LALSimEOSMultiParts * eos, int piece_id);

double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundGeometrizedOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos);

/* FUNCTIONS FOR INTERPOLATED EOS VALUES WITH SI UNITS */

double XLALSimNeutronStarEOSEnergyDensityOfPressure(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressure(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressure(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSound(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfEnergyDensity(double e,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfRestMassDensity(double rho,
    LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressure(double p,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressure(double p,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressure(double p,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSound(double h,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfEnergyDensity(double e,
    LALSimEOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfRestMassDensity(double rho,
    LALSimEOSMultiParts * eos, int piece_id);

double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressure(double p,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressure(double p,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressure(double p,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundOfPseudoEnthalpy(double h, LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfEnergyDensity(double e,
    LALSimEOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfRestMassDensity(double rho,
    LALSimEOSMultiParts * eos);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname);

LALSimEOSMultiParts *XLALSimNeutronStarEOSMultiPartsByName(const char *name);
LALSimEOSMultiParts *XLALSimNeutronStarEOSFromFilePhaseTransition(const char *fname);

#ifndef SWIG /* exclude from SWIG interface: double* array params are not handled correctly by SWIG */
/* Python users should use XLALSimNeutronStarEOSFromArrays instead, which takes REAL8Vector inputs */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabData(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat);

/* Python users should use XLALSimNeutronStarEOSFromArraysPhaseTransition instead, which takes REAL8Vector inputs */
LALSimEOSMultiParts *XLALSimNeutronStarEOSFromTabDataPhaseTransition(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat);
LALSimEOSMultiParts *XLALSimNeutronStarEOSFromTabDataPhaseTransitionChoiceDirtyPT(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat, int dirty);
#endif /* SWIG */

LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromArrays(
    const REAL8Vector *energy_density, const REAL8Vector *pressure);

LALSimEOSMultiParts *XLALSimNeutronStarEOSFromArraysPhaseTransition(
    const REAL8Vector *energy_density, const REAL8Vector *pressure);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSPolytrope(double Gamma,
    double reference_pressure_si, double reference_density_si);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterPiecewisePolytrope(double
    logp1_si, double gamma1, double gamma2, double gamma3);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSSpectralDecomposition(double
    gamma[], int size);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterSpectralDecomposition(
    double SDgamma0, double SDgamma1, double SDgamma2, double SDgamma3);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSDynamicAnalytic(double parameters[],
    size_t nsec, int causal);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceDynamicPolytrope(double g0,
    double log10p1_si, double g1, double log10p2_si, double g2);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceCausalAnalytic(double v1,
    double log10p1_si, double v2, double log10p2_si, double v3);
int XLALSimNeutronStarEOS4ParamSDGammaCheck(double g0, double g1, double g2, double g3);
int XLALSimNeutronStarEOS4ParamSDViableFamilyCheck(double g0, double g1, double g2, double g3);
int XLALSimNeutronStarEOS3PDViableFamilyCheck(double p0, double log10p1_si, double p1, double log10p2_si, double p2, int causal);


/* TOV and Love Number solver ROUTINES */

int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarTOVODEIntegrateWithTolerance(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

int XLALSimNeutronStarVirialODEIntegrate(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarVirialODEIntegrateWithTolerance(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

void XLALSimNeutronStarTOVODEExtendedIntegrate(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, LALSimEOSMultiParts *eos);

void XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, LALSimEOSMultiParts * eos, double epsrel);

void XLALSimNeutronStarTOVODEMiniIntegrate(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, LALSimEOSMultiParts *eos);

void XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, LALSimEOSMultiParts * eos, double epsrel);

/* NEUTRON STAR ASTROPHYSICAL PARAMETER'S RELATED ROUTINES */

void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam);
void XLALDestroySimNeutronStarMultiBranchFamily(FamMultiParts * fam);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(LALSimNeutronStarEOS * eos);
FamMultiParts * XLALCreateSimNeutronStarFamilyPT(LALSimEOSMultiParts * eos, int min_fam);
FamMultiParts * XLALCreateSimNeutronStarFamilyPTWithPcmin(LALSimEOSMultiParts * eos, int min_fam, double logPcmin);

int XLALSimNeutronStarFamNumberOfBranches(FamMultiParts *fam);
double XLALSimNeutronStarFamMassTOVLimit(FamMultiParts *fam);
double XLALSimNeutronStarFamBranchMinMass(int branch, FamMultiParts *fam);
double XLALSimNeutronStarFamBranchMinCentralPressure(int branch, FamMultiParts *fam);
double XLALSimNeutronStarFamBranchMaxMass(int branch, FamMultiParts *fam);
double XLALSimNeutronStarFamBranchMaxCentralPressure(int branch, FamMultiParts *fam);
double XLALSimNeutronStarFamBranchRadius(double m, int branch, FamMultiParts * fam);
double XLALSimNeutronStarFamBranchCentralPressure(double m, int branch, FamMultiParts * fam);
double XLALSimNeutronStarFamBranchMass(double p, int branch, FamMultiParts * fam);
double XLALSimNeutronStarFamBranchBaryonicMass(double m, int branch, FamMultiParts * fam);
double XLALSimNeutronStarFamBranchLoveNumberK2(double m, int branch, FamMultiParts * fam);
double XLALSimNeutronStarFamBranchLoveNumberK3(double m, int branch, FamMultiParts * fam);
double XLALSimNeutronStarFamBranchLoveNumberK4(double m, int branch, FamMultiParts * fam);

double XLALSimNeutronStarFamMinimumMass(LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarMaximumMass(LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarCentralPressure(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarRadius(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarLoveNumberK2(double m, LALSimNeutronStarFamily * fam);

#endif /* _LALSIMNEUTRONSTAR_H */

/** @} */
