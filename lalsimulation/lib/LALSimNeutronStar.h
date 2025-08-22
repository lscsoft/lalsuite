/*
 * Copyright (C) 2013 J. Creighton, B. Lackey
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
 * @author Jolien Creighton, Benjamin Lackey
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


/* CONSTANTS */

/** Factor to convert density in kg/m^3 to geometerized units of m^-2. */
#define LAL_G_C2_SI ((LAL_G_SI) / ((double)(LAL_C_SI) * (double)(LAL_C_SI)))

/** Factor to convert pressure in Pa to geometerized units of m^-2. */
#define LAL_G_C4_SI ((LAL_G_C2_SI) / ((double)(LAL_C_SI) * (double)(LAL_C_SI)))

/** Nuclear density in kg m^-3. */
#define LAL_NUCLEAR_DENSITY_SI 2.8e17

/** Nuclear density in geometrized units of m^-2. */
#define LAL_NUCLEAR_DENSITY_GEOM_SI ((LAL_NUCLEAR_DENSITY_SI) * (LAL_G_C2_SI))



/* EOS ROUTINES */

/** Incomplete type for the neutron star Equation of State (EOS). */
typedef struct tagLALSimNeutronStarEOS LALSimNeutronStarEOS;
//CUTER-dev
typedef struct tagEOSMultiParts EOSMultiParts;
/** Recognised names of equations of state */
extern const char * const lalSimNeutronStarEOSNames[111];

/** Incomplete type for a neutron star family having a particular EOS. */
typedef struct tagLALSimNeutronStarFamily LALSimNeutronStarFamily;
typedef struct tagFamMultiParts FamMultiParts;

void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos);
void XLALDestroySimNeutronStarEOSMultiParts(EOSMultiParts * eos);
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname);
//CUTER-dev
EOSMultiParts *XLALSimNeutronStarEOSFromFilePhaseTransition(const char *fname);
EOSMultiParts *XLALSimNeutronStarEOSFromFilePT(const char *fname);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabData(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat);

EOSMultiParts *XLALSimNeutronStarEOSFromTabDataPhaseTransition( double *nbdat, double *edat, double *pdat,
                                                                    double *mubdat, double *muedat, double *hdat,
                                                                    double *yedat, double *cs2dat, size_t ndat);

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

double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPressure(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPressureGeometerized(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMinEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsMinEnthalpy(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsMaxPseudoEnthalpy(EOSMultiParts * eos);
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsMinAcausalPseudoEnthalpy(EOSMultiParts * eos);
int XLALSimNeutronStarEOSMultiPartsNumber(EOSMultiParts * eos);
LALSimNeutronStarEOS * XLALSimNeutronStarEOSPart(EOSMultiParts * eos, int part_number);

double XLALSimNeutronStarEOSEnergyDensityOfPressure(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressure(double p,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressure(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressure(double p,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpyGeometerized(double
    h, EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressure(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressure(double p,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSSpeedOfSound(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsSpeedOfSound(double h,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSPressureOfEnergyDensity(double e,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMutliPartsPressureOfEnergyDensity(double e,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSPressureOfRestMassDensity(double rho,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsPressureOfRestMassDensity(double rho,
    EOSMultiParts * eos, int part_number);

double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressureGeometerized(double p,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressureGeometerized(double p,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMutliPartsPressureOfPseudoEnthalpyGeometerized(double h,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double
    h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpyGeometerized(double
    h, EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double
    h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressureGeometerized(double p,
    EOSMultiParts * eos, int part_number);
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundGeometerized(double h,
    EOSMultiParts * eos, int part_number);

/* TOV ROUTINES */

int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarTOVODEIntegrateWithTolerance(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos, double epsrel);

int XLALSimNeutronStarVirialODEIntegrate(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarVirialODEIntegrateWithTolerance(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos, double epsrel);

//CUTER-dev
int XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
             double *love_number_k2, double *love_number_k3, double *love_number_k4,
             double central_pressure_si,
             EOSMultiParts * eos,
             double epsrel, int flag_mini);
//CUTER-dev
int XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(double *radius, double *mass,
             double *love_number_k2, double central_pressure_si,
             EOSMultiParts * eos,
             double epsrel);

/* MASS-RADIUS TYPE RELATIONSHIP ROUTINES */

void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam);
void XLALDestroySimNeutronStarMultiBranchFamily(FamMultiParts * fam);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(
    LALSimNeutronStarEOS * eos);
FamMultiParts * XLALCreateSimNeutronStarFamilyPT(
    EOSMultiParts * eos, int min_fam);

int XLALSimNeutronStarFamNumberOfBranches(FamMultiParts *fam);
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
double XLALSimNeutronStarCentralPressure(double m,
    LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarRadius(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarLoveNumberK2(double m, LALSimNeutronStarFamily * fam);

#endif /* _LALSIMNEUTRONSTAR_H */

/** @} */
