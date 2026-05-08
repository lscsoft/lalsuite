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
typedef struct tagLALSimNeutronStarEOSPiece LALSimNeutronStarEOSPiece;
typedef struct tagLALSimNeutronStarEOS LALSimNeutronStarEOS;
/** Recognised names of equations of state */
extern const char * const lalSimNeutronStarEOSNames[111];

/** Incomplete type for a neutron star family having a particular EOS. */
typedef struct tagLALSimNeutronStarBranch LALSimNeutronStarBranch;
typedef struct tagLALSimNeutronStarFamily LALSimNeutronStarFamily;
//TODO remove the struct for static ones ???
void XLALDestroySimNeutronStarEOSPiece(LALSimNeutronStarEOSPiece * eos);
void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos);

/* FUNCTIONS FOR KEY EOS-PIECE VALUES */
char *XLALSimNeutronStarEOSPieceName(LALSimNeutronStarEOSPiece * eos);

double XLALSimNeutronStarEOSPieceMinPressureGeometrized(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMinPressure(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMaxPressureGeometrized(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMaxPressure(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMinPseudoEnthalpy(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMinAcausalPseudoEnthalpy(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMaxPseudoEnthalpy(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMinEnergyDensityGeometrized(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMinEnergyDensity(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMaxEnergyDensityGeometrized(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMaxEnergyDensity(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMinRestMassDensityGeometrized(LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceMaxRestMassDensityGeometrized(LALSimNeutronStarEOSPiece * eos); //TODO should there be the unit version ?

/* FUNCTIONS FOR KEY EOS VALUES */

char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos);
int XLALSimNeutronStarEOSNumberPieces(LALSimNeutronStarEOS * eos);
LALSimNeutronStarEOSPiece * XLALSimNeutronStarEOSSelectPiece(LALSimNeutronStarEOS * eos, int piece_id);

double XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinPressurePerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxPressurePerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinEnergyDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinEnergyDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxEnergyDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMinRestMassDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id);

double XLALSimNeutronStarEOSMinPressureGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinPressure(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPressureGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinEnergyDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinEnergyDensity(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxEnergyDensity(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinRestMassDensityGeometrized(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(LALSimNeutronStarEOS * eos);


/* FUNCTIONS FOR INTERPOLATED EOS-PIECE VALUES */

double XLALSimNeutronStarEOSPiecePseudoEnthalpyOfPressureGeometrized(double p, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPiecePseudoEnthalpyOfPressure(double p, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceEnergyDensityOfPressureGeometrized(double p, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceEnergyDensityOfPressure(double p, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceEnergyDensityDerivOfPressureGeometrized(double p, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceEnergyDensityDerivOfPressure(double p, LALSimNeutronStarEOSPiece * eos); // TODO see function in EOS.c

double XLALSimNeutronStarEOSPiecePressureOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPiecePressureOfPseudoEnthalpy(double h, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceEnergyDensityOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceEnergyDensityOfPseudoEnthalpy(double h, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceRestMassDensityOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceRestMassDensityOfPseudoEnthalpy(double h, LALSimNeutronStarEOSPiece * eos);

double XLALSimNeutronStarEOSPieceSpeedOfSoundOfPseudoEnthalpyGeometrized(double h, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPieceSpeedOfSoundOfPseudoEnthalpy(double h, LALSimNeutronStarEOSPiece * eos);

double XLALSimNeutronStarEOSPiecePressureOfEnergyDensityGeometrized(double e, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPiecePressureOfEnergyDensity(double e, LALSimNeutronStarEOSPiece * eos);

double XLALSimNeutronStarEOSPiecePressureOfRestMassDensity(double rho, LALSimNeutronStarEOSPiece * eos);
double XLALSimNeutronStarEOSPiecePressureOfRestMassDensityGeometrized(double rho, LALSimNeutronStarEOSPiece * eos);


/* FUNCTIONS FOR INTERPOLATED EOS VALUES */

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



// /* DEPRECATED: old misspelling of Geometrized */
// double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos);
// double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h, LALSimNeutronStarEOS * eos);


LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name);

#ifndef SWIG /* exclude from SWIG interface: double* array params are not handled correctly by SWIG */
/* Python users should use XLALSimNeutronStarEOSFromArrays instead, which takes REAL8Vector inputs */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabDataChoiceDirtyPT(double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat, int dirty);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabData( double *nbdat, double *edat, double *pdat,
    double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat);
#endif /* SWIG */

LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromArrays(
    const REAL8Vector *energy_density, const REAL8Vector *pressure);







/* FUNCTIONS FOR EOS PARAMETRISATIONS */

//TODO make the wrapper for those

LALSimNeutronStarEOS *XLALSimNeutronStarEOSPolytrope(double Gamma,
    double reference_pressure_si, double reference_density_si);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterPiecewisePolytrope(double
    logp1_si, double gamma1, double gamma2, double gamma3);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSSpectralDecomposition(double gamma[], int size);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterSpectralDecomposition(double SDgamma0, double SDgamma1, double SDgamma2, double SDgamma3);
int XLALSimNeutronStarEOS4ParamSDGammaCheck(double g0, double g1, double g2, double g3);
int XLALSimNeutronStarEOS4ParamSDViableFamilyCheck(double g0, double g1, double g2, double g3);
int XLALSimNeutronStarEOS3PDViableFamilyCheck(double p0, double log10p1_si, double p1, double log10p2_si, double p2, int causal);

//TODO finish this
LALSimNeutronStarEOS *XLALSimNeutronStarEOSDynamicAnalytic(double parameters[],
    size_t nsec, int causal);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceDynamicPolytrope(double g0,
    double log10p1_si, double g1, double log10p2_si, double g2);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS3PieceCausalAnalytic(double v1,
    double log10p1_si, double v2, double log10p2_si, double v3);


/* FUNCTIONS FOR NS ASTROPHYSICAL SOLVER */
/* TOV and Love Number solver ROUTINES */
int XLALSimNeutronStarTOVODEIntegrateWithToleranceEOSPiece(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOSPiece * eos, double epsrel);
int XLALSimNeutronStarTOVODEIntegrateWithTolerance(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos, double epsrel);
int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarVirialODEIntegrate(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos);

int XLALSimNeutronStarVirialODEIntegrateWithTolerance(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

void XLALSimNeutronStarTOVODEExtendedIntegrate(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, LALSimNeutronStarEOS *eos);

void XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
    double *love_number_k2, double *love_number_k3, double *love_number_k4,
    double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

void XLALSimNeutronStarTOVODEMiniIntegrate(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, LALSimNeutronStarEOS *eos);

void XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(double *radius, double *mass, double *love_number_k2,
    double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel);

/* NEUTRON STAR ASTROPHYSICAL PARAMETER'S RELATED ROUTINES */

void XLALDestroySimNeutronStarBranch(LALSimNeutronStarBranch * branch);
void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam);

LALSimNeutronStarBranch * XLALCreateSimNeutronStarBranch(LALSimNeutronStarEOSPiece * eos);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamilyWithPcmin(LALSimNeutronStarEOS * eos, int min_fam, double logPcmin);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(LALSimNeutronStarEOS * eos, int min_fam);
int XLALSimNeutronStarFamNumberOfBranches(LALSimNeutronStarFamily *fam);
LALSimNeutronStarBranch * XLALSimNeutronStarFamSelectBranch(LALSimNeutronStarFamily * fam, int branch_id);
double XLALSimNeutronStarFamMassTOVLimit(LALSimNeutronStarFamily *fam);


double XLALSimNeutronStarBranchMinMass(LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamMinMassPerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMinMass(LALSimNeutronStarFamily *fam);


double XLALSimNeutronStarBranchMaxMass(LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamMaxMassPerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMaxMass(LALSimNeutronStarFamily *fam);

double XLALSimNeutronStarBranchMinCentralPressure(LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamMinCentralPressurePerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMinCentralPressure(LALSimNeutronStarFamily *fam);

double XLALSimNeutronStarBranchMaxCentralPressure(LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamMaxCentralPressurePerBranch(LALSimNeutronStarFamily *fam, int branch_id);
double XLALSimNeutronStarFamMaxCentralPressure(LALSimNeutronStarFamily *fam);


double XLALSimNeutronStarBranchRadiusOfMass(double m, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamRadiusOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);

double XLALSimNeutronStarBranchCentralPressureOfMass(double m, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamCentralPressureOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);

double XLALSimNeutronStarBranchMassOfCentralPressure(double p, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamMassOfCentralPressurePerBranch(double p, LALSimNeutronStarFamily * fam, int branch_id);

double XLALSimNeutronStarBranchBaryonicMassOfMass(double m, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamBaryonicMassOfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);

double XLALSimNeutronStarBranchLoveNumberK2OfMass(double m, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamLoveNumberK2OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);

double XLALSimNeutronStarBranchLoveNumberK3OfMass(double m, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamLoveNumberK3OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);

double XLALSimNeutronStarBranchLoveNumberK4OfMass(double m, LALSimNeutronStarBranch * branch);
double XLALSimNeutronStarFamLoveNumberK4OfMassPerBranch(double m, LALSimNeutronStarFamily * fam, int branch_id);

//TODO finish this
// double XLALSimNeutronStarFamMinimumMass(LALSimNeutronStarFamily * fam);
// double XLALSimNeutronStarMaximumMass(LALSimNeutronStarFamily * fam);
// double XLALSimNeutronStarCentralPressure(double m, LALSimNeutronStarFamily * fam);
// double XLALSimNeutronStarRadius(double m, LALSimNeutronStarFamily * fam);
// double XLALSimNeutronStarLoveNumberK2(double m, LALSimNeutronStarFamily * fam);

#endif /* _LALSIMNEUTRONSTAR_H */

/** @} */
