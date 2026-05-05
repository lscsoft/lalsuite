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
    double hmin;
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

struct tagLALSimEOSMultiParts{
  char name[LALNameLength];
  int number_of_parts;
  double pmax;
  double pmin;
  double hmax;
  double hmin;
  double hMinAcausal;
  LALSimNeutronStarEOS ** eos_piece;
  void (*free) (LALSimEOSMultiParts * myself);
};


/* This function finds the id number of the piece EoS containing
 * the pressure value p, in a LALSimEOSMultiParts structure.
 */
static int find_eos_piece_pressure(double p, LALSimEOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    for (int j = 0 ; j < number_pieces; j++){
        double pmin =  XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, j);
        double pmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometrized(eos, j);
        if (p >= pmin && p <= pmax) return j;
    }
    XLALPrintError("Pressure not found in EoS piece.\n");
    return XLAL_FAILURE;
}


/* This function finds the id number of the piece EoS containing
 * the pseudo-enthalpy value h, in a LALSimEOSMultiParts structure.
 */
static int find_eos_piece_enthalpy(double h, LALSimEOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    for (int j = 0 ; j < number_pieces; j++){
        double hmin =  XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, j);
        double hmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxPseudoEnthalpy(eos, j);
        if (h >= hmin && h <= hmax) return j;
    }
    XLALPrintError("Enthalpy not found in EoS piece.\n");
    return XLAL_FAILURE;
}


/* This function finds the id number of the piece EoS containing
 * the rest-mass density value rho, in a LALSimEOSMultiParts structure.
 */
static int find_eos_piece_rest_mass_density(double rho, LALSimEOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    for (int j = 0 ; j < number_pieces; j++){
        double hmin =  XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, j);
        double hmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxPseudoEnthalpy(eos, j);
        double rhomin = XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometrized(hmin, eos, j);
        double rhomax = XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometrized(hmax, eos, j);
        if (rho >= rhomin && rho <= rhomax) return j;
    }
    XLALPrintError("Rest mass density not found in EoS piece.\n");
    return XLAL_FAILURE;
}

/* This function finds the id number of the piece EoS containing
 * the energy density value e, in a LALSimEOSMultiParts structure.
 */
static int find_eos_piece_energy_density(double e, LALSimEOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    for (int j = 0 ; j < number_pieces; j++){
        double hmin =  XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, j);
        double hmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxPseudoEnthalpy(eos, j);
        double emin = XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometrized(hmin, eos, j);
        double emax = XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometrized(hmax, eos, j);
        if (e >= emin && e <= emax) return j;
    }
    XLALPrintError("Energy density not found in EoS piece.\n");
    return XLAL_FAILURE;
}

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
    if (eos == NULL)
        return;

    if (eos->free != NULL)
        eos->free(eos);
    return;
}


/**
 * @brief Frees the memory associated with a pointer to an LALSimEOSMultiParts structure.
 * @param eos Pointer to the LALSimEOSMultiParts structure to be freed.
 */
void XLALDestroySimNeutronStarEOSMultiParts(LALSimEOSMultiParts * eos)
{
    if (!eos) return;

    if (eos->eos_piece) {
        for (int i = 0; i < eos->number_of_parts; i++) {
            if (eos->eos_piece[i] && eos->eos_piece[i]->free) eos->eos_piece[i]->free(eos->eos_piece[i]);
        }
        LALFree(eos->eos_piece);
    }

    LALFree(eos);
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

/* FUNCTIONS FOR KEY EOS VALUES */

/**
 * @brief The name of the equation of state.
 * @param[in] eos Pointer to the EOS structure.
 * @return Pointer to a string containing the name of the equation of state.
 * @warning The pointer returned might be shallow and might be left
 * dangling if the @a eos structure is freed.
 */
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos)
{
    return XLALStringDuplicate(eos->name);
}

/**
 * @brief Returns the maximum pressure of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The maximum pressure of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMaxPressureGeometrized(LALSimNeutronStarEOS *
    eos)
{
    return eos->pmax;
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSMaxPressureGeometrized */
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSMaxPressureGeometrized");
    return XLALSimNeutronStarEOSMaxPressureGeometrized(eos);
}

/**
 * @brief Returns the minimum pressure of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The minimum pressure of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMinPressureGeometrized(LALSimNeutronStarEOS * eos)
{
    return eos->p_of_h(XLALSimNeutronStarEOSMinPseudoEnthalpy(eos), eos);
}

/**
 * @brief Returns the maximum pressure of the EOS in Pa.
 * @param eos Pointer to the EOS structure.
 * @return The maximum pressure of the EOS in Pa.
 */
double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMaxPressureGeometrized(eos)/LAL_G_C4_SI;
}


/**
 * @brief Returns the minimum pressure of the EOS in Pa.
 * @param eos Pointer to the EOS structure.
 * @return The minimum pressure of the EOS in Pa.
 */
double XLALSimNeutronStarEOSMinPressure(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinPressureGeometrized(eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the EOS (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The maximum pseudo-enthalpy of the EOS (dimensionless).
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

/**
 * @brief Returns the minimum pseudo-enthalpy of the EOS.
 * @param eos Pointer to the EOS structure.
 * @return The minimum pseudo-enthalpy of the EOS.
 */
double XLALSimNeutronStarEOSMinPseudoEnthalpy(LALSimNeutronStarEOS *
    eos)
{
    return eos->hmin;
}


/**
 * @brief Returns the maximum energy density of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The maximum energy density of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(LALSimNeutronStarEOS *
    eos)
{
    return eos->e_of_p(XLALSimNeutronStarEOSMaxPressureGeometrized(eos), eos);
}


/**
 * @brief Returns the minimum energy density of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The minimum energy density of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMinEnergyDensityGeometrized(LALSimNeutronStarEOS *
    eos)
{
    return eos->e_of_p(XLALSimNeutronStarEOSMinPressureGeometrized(eos), eos);
}

/**
 * @brief Returns the maximum rest-mass density of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The maximum rest-mass density of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(LALSimNeutronStarEOS *
    eos)
{
    return eos->rho_of_h(XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos), eos);
}


/**
 * @brief Returns the minimum rest-mass density of the EOS in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The minimum rest-mass density of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMinRestMassDensityGeometrized(LALSimNeutronStarEOS *
    eos)
{
    return eos->rho_of_h(XLALSimNeutronStarEOSMinPseudoEnthalpy(eos), eos);
}

/**
 * @brief The name of the equation of state.
 * @param[in] eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return Pointer to a string containing the name of the equation of state.
 * @warning The pointer returned might be shallow and might be left
 * dangling if the @a eos structure is freed.
 */
char *XLALSimNeutronStarEOSMultiPartsName(LALSimEOSMultiParts * eos)
{
    return XLALStringDuplicate(eos->name);
}

/**
 * @brief Returns the number of pieces in the equation of state
 * divided by phase transitions in the LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns N+1.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The number of pieces in the LALSimEOSMultiParts structure.
 */
int XLALSimNeutronStarEOSMultiPartsNumber(LALSimEOSMultiParts * eos)
{
    return eos->number_of_parts;
}

/**
 * @brief Returns the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function returns a LALSimNeutronStarEOS
 * structure for the equation of piece with ID number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id ID number for the equation of state piece.
 * @return Pointer of type LALSimNeutronStarEOS, for the ith
 * piece of the LALSimEOSMultiParts structure.
 */
LALSimNeutronStarEOS * XLALSimNeutronStarEOSPart(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id < 0 || piece_id >= eos->number_of_parts) {
        XLAL_ERROR_NULL(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = eos->eos_piece[piece_id];
    return eos_piece;
}


/**
 * @brief Returns the maximum pressure of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The maximum pressure of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMaxPressureGeometrized(LALSimEOSMultiParts * eos)
{
    return eos->pmax;
}

/**
 * @brief Returns the minimum pressure of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The minimum pressure of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMinPressureGeometrized(LALSimEOSMultiParts * eos)
{
    return eos->pmin;
}

/**
 * @brief Returns the maximum pressure of the equation of state in Pa.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The maximum pressure of the equation of state in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsMaxPressure(LALSimEOSMultiParts * eos)
{
    return XLALSimNeutronStarEOSMultiPartsMaxPressureGeometrized(eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the minimum pressure of the equation of state in Pa.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The minimum pressure of the equation of state in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsMinPressure(LALSimEOSMultiParts * eos)
{
    return XLALSimNeutronStarEOSMultiPartsMinPressureGeometrized(eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the minimum pseudo-enthalpy of the equation of state.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The minimum pseudo-enthalpy of the equation of state.
 */
double XLALSimNeutronStarEOSMultiPartsMinPseudoEnthalpy(LALSimEOSMultiParts * eos)
{
    return eos->hmin;
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the equation of state (dimensionless).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The maximum pseudo-enthalpy of the equation of state (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsMaxPseudoEnthalpy(LALSimEOSMultiParts * eos)
{
    return eos->hmax;
}

/**
 * @brief Returns the minimum pseudo-enthalpy at which equation of state becomes acausal
 * (speed of sound > speed of light) (dimensionless).
 * @param eos Pointer to the equation of state structure.
 * @return The minimum pseudo-enthalpy at which equation of state becomes acausal
 * (speed of sound > speed of light).
 */
double XLALSimNeutronStarEOSMultiPartsMinAcausalPseudoEnthalpy(LALSimEOSMultiParts * eos)
{
    return eos->hMinAcausal;
}

/**
 * @brief Returns the maximum energy density of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The maximum energy density of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMaxEnergyDensityGeometrized(LALSimEOSMultiParts * eos)
{
    int max_id_piece = XLALSimNeutronStarEOSMultiPartsNumber(eos) - 1 ;
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, max_id_piece);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrized(eos_piece);
    return eos_piece->e_of_p(pmax, eos_piece);
}

/**
 * @brief Returns the minimum energy density of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The minimum energy density of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMinEnergyDensityGeometrized(LALSimEOSMultiParts * eos)
{
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, 0);
    double pmin =XLALSimNeutronStarEOSMinPressureGeometrized(eos_piece);
    return eos_piece->e_of_p(pmin, eos_piece);
}

/**
 * @brief Returns the maximum rest-mass density of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The maximum rest-mass of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMaxRestMassDensityGeometrized(LALSimEOSMultiParts * eos)
{
    int max_id_piece = XLALSimNeutronStarEOSMultiPartsNumber(eos) - 1 ;
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, max_id_piece);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos_piece);
    return eos_piece->rho_of_h(hmax, eos_piece);
}

/**
 * @brief Returns the minimum rest-mass of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The minimum rest-mass of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMinRestMassDensityGeometrized(LALSimEOSMultiParts * eos)
{
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, 0);
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos_piece);
    return eos_piece->rho_of_h(hmin, eos_piece);
}


/**
 * @brief Returns the maximum pressure of the ith equation of state piece
 * of LALSimEOSMultiParts structure in geometrized units m^-2.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pressure (geometrized units)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum pressure of the equation of state piece
 * in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometrized(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMaxPressureGeometrized(eos_piece);
}

/**
 * @brief Returns the minimum pressure of the ith equation of state piece
 * of LALSimEOSMultiParts structure in geometrized units m^-2.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the minimum pressure (geometrized units)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum pressure of the equation of state piece
 * in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    double hmin = XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id);
    return XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometrized(hmin, eos, piece_id);
}

/**
 * @brief Returns the maximum pressure of the ith equation of sstate piece
 * of LALSimEOSMultiParts structure in Pa.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum pressure of the equation of state piece in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxPressure(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    double pmax;
    pmax = XLALSimNeutronStarEOSMaxPressureGeometrized(eos_piece);
    pmax /= LAL_G_C4_SI;
    return pmax;
}

/**
 * @brief Returns the minimum pressure of the ith equation of state piece
 * of LALSimEOSMultiParts structure in Pa.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the minimum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum pressure of the equation of state piece in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinPressure(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    double hmin = XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id);
    double pmin;
    pmin = XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometrized(hmin, eos, piece_id);
    pmin /= LAL_G_C4_SI;
    return pmin;
}

/**
 * @brief Returns the minimum pseudo-enthalpy of the ith equation of state piece
 * of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the minimum pseudo-enthalpy
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum pseudo-enthalpy of the equation of state piece.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(LALSimEOSMultiParts * eos, int piece_id){
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMinPseudoEnthalpy(eos_piece);
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the ith equation of state piece
 * of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pseudo-enthalpy
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum pseudo-enthalpy of the equation of state piece.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxPseudoEnthalpy(LALSimEOSMultiParts * eos, int piece_id){
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos_piece);
}

/**
 * @brief Returns the minimum acausal pseudo-enthalpy of the ith equation
 * of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the minimum acausal pseudo-enthalpy
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum acausal pseudo-enthalpy of the equation of state piece.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinAcausalPseudoEnthalpy(LALSimEOSMultiParts * eos, int piece_id){
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(eos_piece);
}

/**
 * @brief Returns the maximum energy density of the ith equation of sstate piece
 * of LALSimEOSMultiParts structure in geometrized units.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum energy density of the equation of state piece in geometrized units.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxEnergyDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the minimum energy density of the ith equation of sstate piece
 * of LALSimEOSMultiParts structure in geometrized units.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum energy density of the equation of state piece in geometrized units.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinEnergyDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMinEnergyDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the maximum rest-mass density of the ith equation of sstate piece
 * of LALSimEOSMultiParts structure in geometrized units.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum rest-mass density of the equation of state piece in geometrized units.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxRestMassDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the minimum rest-mass density of the ith equation of sstate piece
 * of LALSimEOSMultiParts structure in geometrized units.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum rest-mass density of the equation of state piece in geometrized units.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinRestMassDensityGeometrized(LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMinRestMassDensityGeometrized(eos_piece);
}
/* FUNCTIONS FOR INTERPOLATED EOS VALUES WITH GEOMETRIZED UNITS */

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * pressure in geometrized units (m^-2).
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The energy density in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos)
{
    return eos->e_of_p(p, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized */
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized");
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(p, eos);
}


/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure in
 * geometrized units (m^-2).
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos)
{
    return eos->h_of_p(p, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized");
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(p, eos);
}


/**
 * @brief Returns the pressure in geometrized units (m^-2) at a given value of
 * the dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The pressure in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(double h,
    LALSimNeutronStarEOS * eos)
{
    return eos->p_of_h(h, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized");
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos);
}


/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * value of the dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The energy density in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimNeutronStarEOS * eos)
{
    return eos->e_of_h(h, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized");
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos);
}


/**
 * @brief Returns the rest mass density in geometrized units (m^-2) at a given
 * value of the dimensionless pseudo-enthalpy.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The rest mass density in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimNeutronStarEOS * eos)
{
    return eos->rho_of_h(h, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized");
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h, eos);
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometrized
 * units (m^-2).
 * @param p Pressure in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos)
{
    return eos->dedp_of_p(p, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized */
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized");
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
}


/**
 * @brief Returns the speed of sound in geometrized units (dimensionless)
 * at a given value of the pseudo-enthalpy (dimensionless).
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The speed of sound in geometrized units (dimensionless).
 */
double XLALSimNeutronStarEOSSpeedOfSoundGeometrized(double h,
    LALSimNeutronStarEOS * eos)
{
    return eos->v_of_h(h, eos);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSSpeedOfSoundGeometrized */
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSSpeedOfSoundGeometrized");
    return XLALSimNeutronStarEOSSpeedOfSoundGeometrized(h, eos);
}


/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * pressure in geometrized units (m^-2), in the ith equation of state piece
 * of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the energy density interpolated at
 * a given pressure p, for the equation of state piece number piece_id.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && p < XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pressure is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->e_of_p(p, eos_piece);
}

/**
 * @brief Returns the pseudo-enthalpy at a given pressure in geometrized units (m^-2),
 * in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the pseudo-enthalpy interpolated at
 * a given pressure p, for the equation of state piece number piece_id.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pseudo-enthalpy (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && p < XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pressure is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->h_of_p(p, eos_piece);
}


/**
 * @brief Returns the pressure in geometrized units (m^-2) at a
 * given pseudo-enthalpy, in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the geometrized pressure
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->p_of_h(h, eos_piece);
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a
 * given pseudo-enthalpy, in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the geometrized energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->e_of_h(h, eos_piece);
}


/**
 * @brief Returns the rest-mass density in geometrized units (m^-2) at a
 * given pseudo-enthalpy, in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the geometrized rest-mass density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The rest mass density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->rho_of_h(h, eos_piece);
}


/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometrized
 * units (m^-2), in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the gradient of the energy density
 * interpolated at a given pressure, for the equation of state piece number piece_id.
 * @param p Pressure in geometrized units (m^-2).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && p < XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pressure is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->dedp_of_p(p, eos_piece);
}


/**
 * @brief Returns the speed of sound in geometrized units (dimensionless)
 * at a given value of the pressure in geometrized units (m^-2),
 * in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the gradient of the energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The speed of sound in geometrized units (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSoundGeometrized(double h,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->v_of_h(h, eos_piece);
}


/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * pressure in geometrized units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the energy density
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->e_of_p(p, eos_piece);
}

/**
 * @brief Returns the pseudo-enthalpy at a given pressure in geometrized units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the pseudo-enthalpy
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The pseudo-enthalpy (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->h_of_p(p, eos_piece);
}


/**
 * @brief Returns the pressure in geometrized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the pressure
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The pressure in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpyGeometrized(double h,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->p_of_h(h, eos_piece);
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the energy density
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->e_of_h(h, eos_piece);
}


/**
 * @brief Returns the rest-mass density in geometrized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * rest-mass density at the pseudo-enthalpy value h, for the correct
 * equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The rest mass density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->rho_of_h(h, eos_piece);
}


/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometrized
 * units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the
 * energy density at the pressure value p, for the correct
 * equation of state piece.
 * @param p Pressure in geometrized units (m^-2).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->dedp_of_p(p, eos_piece);
}


/**
 * @brief Returns the speed of sound in geometrized units (dimensionless)
 * at a given value of the pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * speed of sound at the pseudo-enthalpy value h, for the correct
 * equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The speed of sound in geometrized units (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundGeometrizedOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->v_of_h(h, eos_piece);
}


/* FUNCTIONS FOR INTERPOLATED EOS VALUES WITH SI UNITS */

/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa).
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPressure(double p,
    LALSimNeutronStarEOS * eos)
{
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(p, eos)/LAL_G_C4_SI;
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
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(p, eos);
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
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos)/LAL_G_C4_SI;
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
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h,eos)/LAL_G_C4_SI;
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
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h,eos)/LAL_G_C2_SI;
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
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
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
    return XLALSimNeutronStarEOSSpeedOfSoundGeometrized(h, eos)*LAL_C_SI;
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
    e *= LAL_G_C4_SI;
    return eos->p_of_e(e, eos)/LAL_G_C4_SI;
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
    rho *= LAL_G_C2_SI;
    return eos->p_of_rho(rho, eos)/LAL_G_C4_SI;
}


/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa), in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the energy density
 * interpolated at a given pressure p, for the equation of state
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressure(double p,
    LALSimEOSMultiParts * eos, int piece_id)
{
    p *= LAL_G_C4_SI;
    if (piece_id >= eos->number_of_parts || piece_id < 0){
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && p < XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pressure is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(p, eos_piece)/LAL_G_C4_SI;
}


/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure (Pa), in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the pseudo-enthalpy
 * interpolated at a given pressure p, for the equation of state
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressure(double p,
    LALSimEOSMultiParts * eos, int piece_id)
{
    p *= LAL_G_C4_SI;
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && p < XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pressure is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(p, eos_piece);
}

/**
 * @brief Returns the pressure (Pa) at a given value of the dimensionless
 * pseudo-enthalpy, in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the pressure
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure (Pa).
 */
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos_piece)/LAL_G_C4_SI;
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * dimensionless pseudo-enthalpy, in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos_piece)/LAL_G_C4_SI;
}

/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy, in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the rest-mass density
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0){
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h, eos_piece)/LAL_G_C2_SI;
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure (Pa),
 * in the ith equation of state piece of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the gradient of the energy density
 * interpolated at a given pressure p, for the equation of state
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressure(double p,
    LALSimEOSMultiParts * eos, int piece_id)
{
    p *= LAL_G_C4_SI;
    if (piece_id >= eos->number_of_parts || piece_id < 0){
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && p < XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pressure is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos_piece);
}

/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless), in the ith equation of state piece
 * of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the speed of sound
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSound(double h,
    LALSimEOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0){
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && h < XLALSimNeutronStarEOSMultiPartsPieceMinPseudoEnthalpy(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input pseudo-enthalpy is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSSpeedOfSoundGeometrized(h, eos_piece)*LAL_C_SI;
}


/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3, in the ith equation of state piece
 * of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the pressure
 * interpolated at a given energy density e, for the equation of state
 * piece number piece_id.
 * @param e energy density in J/m^3
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfEnergyDensity(double e,
    LALSimEOSMultiParts * eos, int piece_id)
{
    e *= LAL_G_C4_SI;
    if (piece_id >= eos->number_of_parts || piece_id < 0){
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && e < XLALSimNeutronStarEOSMultiPartsPieceMinEnergyDensityGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input energy density is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->p_of_e(e, eos_piece)/LAL_G_C4_SI;
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3, in the ith equation of state piece
 * of LALSimEOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the pressure
 * interpolated at a given rest-mass density rho, for the equation of state
 * piece number piece_id.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfRestMassDensity(double rho,
    LALSimEOSMultiParts * eos, int piece_id)
{
    rho *= LAL_G_C2_SI;
    if (piece_id >= eos->number_of_parts || piece_id < 0){
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimEOSMultiParts structure is incorrect.");
    }
    if (piece_id > 0 && rho < XLALSimNeutronStarEOSMultiPartsPieceMinRestMassDensityGeometrized(eos, piece_id)) {
        XLAL_ERROR_REAL8(XLAL_EDOM, "Input rest-mass density is below the EOS piece interpolation range.");
    }
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, piece_id);
    return eos_piece->p_of_rho(rho, eos_piece)/LAL_G_C4_SI;
}


/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the energy density
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressure(double p,
    LALSimEOSMultiParts * eos)
{
    p *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(p, eos_piece)/LAL_G_C4_SI;
}

/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure (Pa).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the pseudo-enthalpy
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressure(double p,
    LALSimEOSMultiParts * eos)
{
    p *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(p, eos_piece);
}


/**
 * @brief Returns the pressure (Pa) at a given value of the dimensionless
 * pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the pressure
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The pressure (Pa).
 */
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos_piece)/LAL_G_C4_SI;
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the energy density
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h,eos_piece)/LAL_G_C4_SI;
}


/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * rest-mass density at the pseudo-enthalpy value h,
 * for the correct equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpy(double h,
    LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h, eos_piece)/LAL_G_C2_SI;
}


/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure (Pa).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the
 * gradient of the energy density with respect to the pressure
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressure(double p,
    LALSimEOSMultiParts * eos)
{
    p *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos_piece);
}


/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless).
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * speed of sound at the pseudo-enthalpy value h, for the correct
 * equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundOfPseudoEnthalpy(double h, LALSimEOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return XLALSimNeutronStarEOSSpeedOfSoundGeometrized(h, eos_piece)*LAL_C_SI;
}


/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the energy density value belongs to, then it interpolates the
 * pressure at the energy density value e, for the correct
 * equation of state piece.
 * @param e energy density in J/m^3
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPressureOfEnergyDensity(double e,
    LALSimEOSMultiParts * eos)
{
    e *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_energy_density(e, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->p_of_e(e, eos_piece)/LAL_G_C4_SI;
}


/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3.
 * @details If the equation of state contains N phase transitions,
 * the LALSimEOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N. This function finds automatically which equation of state
 * piece the rest-mass density value belongs to, then it interpolates the
 * pressure at the rest-mass density value rho, for the correct
 * equation of state piece.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the equation of state LALSimEOSMultiParts structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPressureOfRestMassDensity(double rho,
    LALSimEOSMultiParts * eos)
{
    rho *= LAL_G_C2_SI;
    int item_piece = find_eos_piece_rest_mass_density(rho, eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, item_piece);
    return eos_piece->p_of_rho(rho, eos_piece)/LAL_G_C4_SI;
}


/** @} */
/** @} */
