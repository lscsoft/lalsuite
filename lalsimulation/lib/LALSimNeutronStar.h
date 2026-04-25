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
typedef struct tagEOSMultiParts EOSMultiParts;
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

void XLALDestroySimNeutronStarEOSMultiParts(EOSMultiParts * eos);
char *XLALSimNeutronStarEOSMultiPartsName(EOSMultiParts * eos);
int XLALSimNeutronStarEOSMultiPartsNumber(EOSMultiParts * eos);
LALSimNeutronStarEOS * XLALSimNeutronStarEOSPart(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsMaxPressureGeometrized(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinPressureGeometrized(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPressure(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinPressure(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinPseudoEnthalpy(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPseudoEnthalpy(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinAcausalPseudoEnthalpy(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxEnergyDensityGeometrized(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinEnergyDensityGeometrized(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMaxRestMassDensityGeometrized(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsMinRestMassDensityGeometrized(EOSMultiParts * eos);

double XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometrized(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxPressure(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinPressure(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxPseudoEnthalpy(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinAcausalPseudoEnthalpy(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxEnergyDensityGeometrized(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinEnergyDensityGeometrized(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMaxRestMassDensityGeometrized(EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceMinRestMassDensityGeometrized(EOSMultiParts * eos, int piece_id);

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

double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressureGeometrized(double p,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressureGeometrized(double p,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometrized(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometrized(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometrized(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressureGeometrized(double p,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSoundGeometrized(double h,
    EOSMultiParts * eos, int piece_id);

double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressureGeometrized(double p,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressureGeometrized(double p,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpyGeometrized(double h,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpyGeometrized(double h,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpyGeometrized(double h,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressureGeometrized(double p,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundGeometrizedOfPseudoEnthalpy(double h,
    EOSMultiParts * eos);

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
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressure(double p,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressure(double p,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSound(double h,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfEnergyDensity(double e,
    EOSMultiParts * eos, int piece_id);
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfRestMassDensity(double rho,
    EOSMultiParts * eos, int piece_id);

double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressure(double p,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressure(double p,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpy(double h,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressure(double p,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundOfPseudoEnthalpy(double h, EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfEnergyDensity(double e,
    EOSMultiParts * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfRestMassDensity(double rho,
    EOSMultiParts * eos);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname);

EOSMultiParts *XLALSimNeutronStarEOSMultiPartsByName(const char *name);
EOSMultiParts *XLALSimNeutronStarEOSFromFilePhaseTransition(const char *fname);

#ifndef SWIG /* exclude from SWIG interface: double* array params are not handled correctly by SWIG */
/* Python users should use XLALSimNeutronStarEOSFromArrays instead, which takes REAL8Vector inputs */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabData(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat);

/* Python users should use XLALSimNeutronStarEOSFromArraysPhaseTransition instead, which takes REAL8Vector inputs */
EOSMultiParts *XLALSimNeutronStarEOSFromTabDataPhaseTransition(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat);
#endif /* SWIG */

LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromArrays(
    const REAL8Vector *energy_density, const REAL8Vector *pressure);

EOSMultiParts *XLALSimNeutronStarEOSFromArraysPhaseTransition(
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
    double central_pressure_si, EOSMultiParts *eos);

void XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, EOSMultiParts * eos, double epsrel);

void XLALSimNeutronStarTOVODEMiniIntegrate(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, EOSMultiParts *eos);

void XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, EOSMultiParts * eos, double epsrel);

/* NEUTRON STAR ASTROPHYSICAL PARAMETER'S RELATED ROUTINES */

void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam);
void XLALDestroySimNeutronStarMultiBranchFamily(FamMultiParts * fam);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(LALSimNeutronStarEOS * eos);
FamMultiParts * XLALCreateSimNeutronStarFamilyPT(EOSMultiParts * eos, int min_fam);
FamMultiParts * XLALCreateSimNeutronStarFamilyPTWithPcmin(EOSMultiParts * eos, int min_fam, double logPcmin);

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
