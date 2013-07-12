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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
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
    double (*dedp_of_p) (double p, LALSimNeutronStarEOS * myself);
    double (*v_of_h) (double h, LALSimNeutronStarEOS * myself);
    void (*free) (LALSimNeutronStarEOS * myself);
    int datatype;
    LALSimNeutronStarEOSData data;
};

/** @endcond */

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


/** @} */
/** @} */
