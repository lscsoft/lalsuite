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
 * @addtogroup LALSimNeutronStar_h      Header LALSimNeutronStar.h
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

/** Incomplete type for a neutron star family having a particular EOS. */
typedef struct tagLALSimNeutronStarFamily LALSimNeutronStarFamily;

void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos);
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos);

LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname);
LALSimNeutronStarEOS *XLALSimNeutronStarEOSPolytrope(double Gamma,
    double reference_pressure_si, double reference_density_si);
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterPiecewisePolytrope(double
    logp1_si, double gamma1, double gamma2, double gamma3);

double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS *
    eos);
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS *
    eos);

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

double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double
    h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double
    h, LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p,
    LALSimNeutronStarEOS * eos);
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h,
    LALSimNeutronStarEOS * eos);

/* TOV ROUTINES */

int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos);

/* MASS-RADIUS TYPE RELATIONSHIP ROUTINES */

void XLALDestroySimNeutronStarFamily(LALSimNeutronStarFamily * fam);
LALSimNeutronStarFamily * XLALCreateSimNeutronStarFamily(
    LALSimNeutronStarEOS * eos);

double XLALSimNeutronStarMaximumMass(LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarCentralPressure(double m,
    LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarRadius(double m, LALSimNeutronStarFamily * fam);
double XLALSimNeutronStarLoveNumberK2(double m, LALSimNeutronStarFamily * fam);

#endif /* _LALSIMNEUTRONSTAR_H */

/** @} */
