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
 * @addtogroup LALSimNeutronStarEOS_c
 * @brief Provides routines for handling neutron star equations of state.
 * @{
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/FileIO.h>
#include <lal/LALSimNeutronStar.h>
#include <lal/LALSimReadData.h>

/** @cond */

/* Enumeration for type of equation of state data storage. */
enum {
    LALSIM_NEUTRON_STAR_EOS_DATA_TYPE_TABULAR,
    LALSIM_NEUTRON_STAR_EOS_DATA_TYPE_PIECEWISE_POLYTROPE,
    LALSIM_NEUTRON_STAR_EOS_DATA_TYPE_NUMBER
};

/* Incomplete types for the equation of state data storage. */
typedef struct tagLALSimNeutronStarEOSDataTabular
    LALSimNeutronStarEOSDataTabular;
typedef struct tagLALSimNeutronStarEOSDataPolytrope
    LALSimNeutronStarEOSDataPolytrope;
typedef struct tagLALSimNeutronStarEOSDataPiecewisePolytrope
    LALSimNeutronStarEOSDataPiecewisePolytrope;

typedef union tagLALSimNeutronStarEOSData {
    LALSimNeutronStarEOSDataTabular *tabular;
    LALSimNeutronStarEOSDataPolytrope *polytrope;
    LALSimNeutronStarEOSDataPiecewisePolytrope *piecewisePolytrope;
} LALSimNeutronStarEOSData;

/* Contents of the equation of state structure. */
struct tagLALSimNeutronStarEOS {
    char name[LALNameLength];
    double pmax;
    double hmax;
    double hMinAcausal; /* Minimum pseudo-enthalpy at which EOS becomes acausal (speed of sound > 1) */
    double (*e_of_p) (double p, LALSimNeutronStarEOS * myself);
    double (*h_of_p) (double p, LALSimNeutronStarEOS * myself);
    double (*p_of_h) (double h, LALSimNeutronStarEOS * myself);
    double (*e_of_h) (double h, LALSimNeutronStarEOS * myself);
    double (*rho_of_h) (double h, LALSimNeutronStarEOS * myself);
    double (*p_of_e) (double e, LALSimNeutronStarEOS * myself);
    double (*p_of_rho) (double rho, LALSimNeutronStarEOS * myself);
    double (*dedp_of_p) (double p, LALSimNeutronStarEOS * myself);
    double (*v_of_h) (double h, LALSimNeutronStarEOS * myself);
    void (*free) (LALSimNeutronStarEOS * myself);
    int datatype;
    LALSimNeutronStarEOSData data;
};

/** @endcond */

/** Recognised equations of state names */
const char * const lalSimNeutronStarEOSNames[111] = {
        "ALF1", "ALF2", "ALF3", "ALF4",
        "AP1", "AP2", "AP3", "AP4", "APR4_EPP",
        "BBB2", "BGN1H1", "BPAL12",
        "BSK19", "BSK20", "BSK21",
        "ENG", "FPS", "GNH3",
        "GS1", "GS2",
        "H1", "H2", "H3", "H4", "H5", "H6", "H7",
        "MPA1", "MS1B", "MS1B_PP", "MS1_PP", "MS1", "MS2",
        "PAL6", "PCL2", "PS",
        "QMC700",
        "SLY4", "SLY",
        "SQM1", "SQM2", "SQM3",
        "WFF1", "WFF2", "WFF3",
        /* From here, EOSs are coming from CompOSE via Ian */
        "APR", "BHF_BBB2",
        "KDE0V", "KDE0V1", "RS", "SK255", "SK272",
        "SKA", "SKB", "SKI2", "SKI3", "SKI4", "SKI5", "SKI6",
        "SKMP", "SKOP",
        "SLY2", "SLY230A", "SLY9",
        /* This EOS was added by request from
         * http://user.numazu-ct.ac.jp/~sumi/eos/HQC18_submit
         */
        "HQC18",
        /* These EOS files are part of the new EOS framework via Micaela Oertel's group and the LUTH group */
        "GMSR_BSK14_BSK24", "GMSR_DHSL59_BSK24", "GMSR_DHSL69_BSK24", "GMSR_F0_BSK24", "GMSR_H1_BSK24", "GMSR_H2_BSK24",
        "GMSR_H3_BSK24", "GMSR_H4_BSK24", "GMSR_H5_BSK24", "GMSR_LN55_BSK24", "GMSR_SLY5_BSK24", "GPPVA_DD2_BSK24",
        "GPPVA_DDME2_BSK24", "GPPVA_FSU2_BSK24", "GPPVA_FSU2H_BSK24", "GPPVA_NL3WRL55_BSK24",
        "KDE0V_BSK24", "KDE0V1_BSK24", "PCP_BSK24_BSK24", "RG_SLY4_BSK24", "RS_BSK24",
        "SK255_BSK24", "SKA_BSK24", "SKB_BSK24", "SKI2_BSK24", "SKI3_BSK24", "SKI4_BSK24", "SKI6_BSK24",
        "SKOP_BSK24", "SLY2_BSK24", "SLY9_BSK24", "SLY230A_BSK24",
        "XMLSLZ_DDLZ1_BSK24", "XMLSLZ_DDME2_BSK24", "XMLSLZ_DDMEX_BSK24", "XMLSLZ_GM1_BSK24", "XMLSLZ_MTVTC_BSK24",
        "XMLSLZ_NL3_BSK24", "XMLSLZ_PKDD_BSK24", "XMLSLZ_TM1_BSK24", "XMLSLZ_TW99_BSK24", "ABHT_QMC_RMF1_META",
        "ABHT_QMC_RMF2_META", "ABHT_QMC_RMF3_META", "ABHT_QMC_RMF4_META", "BL_CHIRAL_META"
    };

/**
 * @name Destruction routine
 * @{
 */

/**
 * @brief Frees the memory associated with a pointer to an EOS structure.
 * @param eos Pointer to the EOS structure to be freed.
 */
void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos)
{
    eos->free(eos);
    return;
}

/** @} */

/* Tabular Equation of State Code. */
#include "LALSimNeutronStarEOSTabular.c"

/* Piecewise-Polytrope Equation of State Code */
#include "LALSimNeutronStarEOSPiecewisePolytrope.c"

/* Spectral Decomposition Equation of State Code */
#include "LALSimNeutronStarEOSSpectralDecomposition.c"

/* Dynamic Piecewise-Polytrope Equation of State Code */
#include "LALSimNeutronStarEOSDynamicPolytrope.c"

/**
 * @name Routines to access equation of state variables
 * @{
 */

/**
 * @brief The name of the equation of state.
 * @param[in] eos Pointer to the EOS structure.
 * @return Pointer to a string containing the name of the equation of state.
 * @warning The pointer returned might be shallow and might be left
 * dangling if the @a eos structure is freed.
 */
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos)
{
    return eos->name;
}

/**
 * @brief Returns the maximum pressure of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The maximum pressure of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS *
    eos)
{
    return eos->pmax;
}

/**
 * @brief Returns the maximum pressure of the EOS in Pa.
 * @param eos Pointer to the EOS structure.
 * @return The maximum pressure of the EOS in Pa.
 */
double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos)
{
    double pmax;
    pmax = XLALSimNeutronStarEOSMaxPressureGeometerized(eos);
    pmax /= LAL_G_C4_SI;
    return pmax;
}

/**
 * @brief Returns the maximum pseudo enthalpy of the EOS (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The maximum pseudo enthalpy of the EOS (dimensionless).
 */
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos)
{
    return eos->hmax;
}

/**
 * @brief Returns the minimum pseudo-enthalpy at which EOS becomes acausal
 * (speed of sound > speed of light) (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The minimum pseudo-enthalpy at which EOS becomes acausal
 * (speed of sound > speed of light).
 */
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS *
    eos)
{
    return eos->hMinAcausal;
}

/* FUNCTIONS WITH GEOMETERIZED UNITS */

/**
 * @brief Returns the energy density in geometerized units (m^-2) at a given
 * pressure in geometerized units (m^-2).
 * @param p Pressure in geometerized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The energy density in geometerized units (m^-2).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos)
{
    double e;
    e = eos->e_of_p(p, eos);
    return e;
}

/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure in
 * geometerized units (m^-2).
 * @param p Pressure in geometerized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos)
{
    double h;
    h = eos->h_of_p(p, eos);
    return h;
}

/**
 * @brief Returns the pressure in geometerized units (m^-2) at a given value of
 * the dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The pressure in geometerized units (m^-2).
 */
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h,
    LALSimNeutronStarEOS * eos)
{
    double p;
    p = eos->p_of_h(h, eos);
    return p;
}

/**
 * @brief Returns the energy density in geometerized units (m^-2) at a given
 * value of the dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The energy density in geometerized units (m^-2).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double
    h, LALSimNeutronStarEOS * eos)
{
    double e;
    e = eos->e_of_h(h, eos);
    return e;
}

/**
 * @brief Returns the rest mass density in geometerized units (m^-2) at a given
 * value of the dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The rest mass density in geometerized units (m^-2).
 */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double
    h, LALSimNeutronStarEOS * eos)
{
    double rho;
    rho = eos->rho_of_h(h, eos);
    return rho;
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometerized
 * units (m^-2).
 * @param p Pressure in geometerized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos)
{
    double dedp;
    dedp = eos->dedp_of_p(p, eos);
    return dedp;
}

/**
 * @brief Returns the speed of sound in geometerized units (dimensionless)
 * at a given value of the pseudo-enthalpy (dimensionless).
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The speed of sound in geometerized units (dimensionless).
 */
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h,
    LALSimNeutronStarEOS * eos)
{
    double v;
    v = eos->v_of_h(h, eos);
    return v;
}

/* FUNCTIONS WITH SI UNITS */

/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa).
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPressure(double p,
    LALSimNeutronStarEOS * eos)
{
    double e;
    p *= LAL_G_C4_SI;
    e = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(p, eos);
    e /= LAL_G_C4_SI;
    return e;
}

/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure (Pa).
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressure(double p,
    LALSimNeutronStarEOS * eos)
{
    double h;
    p *= LAL_G_C4_SI;
    h = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(p, eos);
    return h;
}

/**
 * @brief Returns the pressure (Pa) at a given value of the dimensionless
 * pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The pressure (Pa).
 */
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos)
{
    double p;
    p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos);
    p /= LAL_G_C4_SI;
    return p;
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos)
{
    double e;
    e = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h,
        eos);
    e /= LAL_G_C4_SI;
    return e;
}

/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos)
{
    double rho;
    rho =
        XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(h,
        eos);
    rho /= LAL_G_C2_SI;
    return rho;
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure (Pa).
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressure(double p,
    LALSimNeutronStarEOS * eos)
{
    double dedp;
    p *= LAL_G_C4_SI;
    dedp =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(p, eos);
    return dedp;
}

/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless).
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSSpeedOfSound(double h, LALSimNeutronStarEOS * eos)
{
    double v;
    v = XLALSimNeutronStarEOSSpeedOfSoundGeometerized(h, eos);
    v *= LAL_C_SI;
    return v;
}

/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3.
 * @param e energy density in J/m^3
 * @param eos Pointer to the EOS structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfEnergyDensity(double e,
    LALSimNeutronStarEOS * eos)
{
    double p;
    e *= LAL_G_C4_SI;
    p = eos->p_of_e(e, eos);
    p /= LAL_G_C4_SI;
    return p;
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the EOS structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfRestMassDensity(double rho,
    LALSimNeutronStarEOS * eos)
{
    double p;
    rho *= LAL_G_C2_SI;
    p = eos->p_of_rho(rho, eos);
    p /= LAL_G_C4_SI;
    return p;
}


/** @} */
/** @} */
