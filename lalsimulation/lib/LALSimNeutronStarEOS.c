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

struct tagEOSMultiParts{
  char name[LALNameLength];
  int number_of_parts;
  double pmax;
  double pmin;
  double hmax;
  double hmin;
  double hMinAcausal;
  LALSimNeutronStarEOS ** eos_part;
  void (*free) (EOSMultiParts * myself);
};


/* This function finds the id number of the piece EoS containing
 * the pressure value p, in a EOSMultiParts structure.
 */
static int find_eos_piece_pressure(double p, EOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    int item_piece = 0;
    for (int j = 0 ; j < number_pieces; j++){
        double pmin =  XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometerized(eos, j);
        double pmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometerized(eos, j);
        if (p >= pmin && p <= pmax) {
            item_piece = j;
            break;
        }
    }
    return item_piece;
}


/* This function finds the id number of the piece EoS containing
 * the pseudo-enthalpy value h, in a EOSMultiParts structure.
 */
static int find_eos_piece_enthalpy(double h, EOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    int item_piece = 0;
    for (int j = 0 ; j < number_pieces; j++){
        double hmin =  XLALSimNeutronStarEOSMultiPartsPieceMinEnthalpy(eos, j);
        double hmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxEnthalpy(eos, j);
        if (h >= hmin && h <= hmax) {
            item_piece = j;
            break;
        }
    }
    return item_piece;
}


/* This function finds the id number of the piece EoS containing
 * the rest-mass density value rho, in a EOSMultiParts structure.
 */
static int find_eos_piece_rest_mass_density(double rho, EOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    int item_piece = 0;
    for (int j = 0 ; j < number_pieces; j++){
        double hmin =  XLALSimNeutronStarEOSMultiPartsPieceMinEnthalpy(eos, j);
        double hmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxEnthalpy(eos, j);
        double rhomin = XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometerized(hmin, eos, j);
        double rhomax = XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometerized(hmax, eos, j);
        if (rho >= rhomin && rho <= rhomax) {
            item_piece = j;
            break;
        }
    }
    return item_piece;
}

/* This function finds the id number of the piece EoS containing
 * the energy density value e, in a EOSMultiParts structure.
 */
static int find_eos_piece_energy_density(double e, EOSMultiParts *eos)
{
    int number_pieces = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    int item_piece = 0;
    for (int j = 0 ; j < number_pieces; j++){
        double hmin =  XLALSimNeutronStarEOSMultiPartsPieceMinEnthalpy(eos, j);
        double hmax =  XLALSimNeutronStarEOSMultiPartsPieceMaxEnthalpy(eos, j);
        double emin = XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometerized(hmin, eos, j);
        double emax = XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometerized(hmax, eos, j);
        if (e >= emin && e <= emax) {
            item_piece = j;
            break;
        }
    }
    return item_piece;
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
    eos->free(eos);
    return;
}

/**
 * @brief Frees the memory associated with a pointer to an EOSMultiParts structure.
 * @param eos Pointer to the EOSMultiParts structure to be freed.
 */
void XLALDestroySimNeutronStarEOSMultiParts(EOSMultiParts * eos)
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
double XLALSimNeutronStarEOSMinEnthalpy(LALSimNeutronStarEOS *
    eos)
{
    return eos->hmin;
}


/**
 * @brief The name of the equation of state.
 * @param[in] eos Pointer to the equation of state EOSMultiParts structure.
 * @return Pointer to a string containing the name of the equation of state.
 * @warning The pointer returned might be shallow and might be left
 * dangling if the @a eos structure is freed.
 */
char *XLALSimNeutronStarEOSMultiPartsName(EOSMultiParts * eos)
{
    return eos->name;
}

/**
 * @brief Returns the number of pieces in the equation of state
 * divided by phase transitions in the EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns N+1.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The number of pieces in the EOSMultiParts structure.
 */
int XLALSimNeutronStarEOSMultiPartsNumber(EOSMultiParts * eos)
{
    return eos->number_of_parts;
}

/**
 * @brief Returns the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function returns a LALSimNeutronStarEOS
 * structure for the equation of piece with ID number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id ID number for the equation of state piece.
 * @return Pointer of type LALSimNeutronStarEOS, for the ith
 * piece of the EOSMultiParts structure.
 */
LALSimNeutronStarEOS * XLALSimNeutronStarEOSPart(EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts) XLAL_ERROR_NULL(XLAL_EDOM);
    LALSimNeutronStarEOS * eos_part = eos->eos_part[piece_id];
    return eos_part;
}


/**
 * @brief Returns the maximum pressure of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The maximum pressure of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMaxPressureGeometerized(EOSMultiParts * eos)
{
    return eos->pmax;
}

/**
 * @brief Returns the minimum pressure of the equation of state in geometrized units m^-2.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The minimum pressure of the equation of state in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsMinPressureGeometerized(EOSMultiParts * eos)
{
    return eos->pmin;
}

/**
 * @brief Returns the maximum pressure of the equation of state in Pa.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The maximum pressure of the equation of state in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsMaxPressure(EOSMultiParts * eos)
{
    double pmax;
    pmax = XLALSimNeutronStarEOSMultiPartsMaxPressureGeometerized(eos);
    pmax /= LAL_G_C4_SI;
    return pmax;
}

/**
 * @brief Returns the minimum pressure of the equation of state in Pa.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The minimum pressure of the equation of state in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsMinPressure(EOSMultiParts * eos)
{
    double pmin;
    pmin = XLALSimNeutronStarEOSMultiPartsMinPressureGeometerized(eos);
    pmin /= LAL_G_C4_SI;
    return pmin;
}

/**
 * @brief Returns the minimum pseudo-enthalpy of the equation of state.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The minimum pseudo-enthalpy of the equation of state.
 */
double XLALSimNeutronStarEOSMultiPartsMinEnthalpy(EOSMultiParts * eos)
{
    return eos->hmin;
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the equation of state (dimensionless).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The maximum pseudo-enthalpy of the equation of state (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsMaxEnthalpy(EOSMultiParts * eos)
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
double XLALSimNeutronStarEOSMultiPartsMinAcausalPseudoEnthalpy(EOSMultiParts * eos)
{
    return eos->hMinAcausal;
}

/**
 * @brief Returns the maximum pressure of the ith equation of state piece
 * of EOSMultiParts structure in geometrized units m^-2.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the maximum pressure (geometrized units)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum pressure of the equation of state piece
 * in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometerized(EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMaxPressureGeometerized(eos_part);
}

/**
 * @brief Returns the minimum pressure of the ith equation of state piece
 * of EOSMultiParts structure in geometrized units m^-2.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the minimum pressure (geometrized units)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum pressure of the equation of state piece
 * in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometerized(EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double hmin = XLALSimNeutronStarEOSMultiPartsPieceMinEnthalpy(eos, piece_id);
    return XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometerized(hmin, eos, piece_id);
}

/**
 * @brief Returns the maximum pressure of the ith equation of sstate piece
 * of EOSMultiParts structure in Pa.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the maximum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum pressure of the equation of state piece in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxPressure(EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    double pmax;
    pmax = XLALSimNeutronStarEOSMaxPressureGeometerized(eos_part);
    pmax /= LAL_G_C4_SI;
    return pmax;
}

/**
 * @brief Returns the minimum pressure of the ith equation of state piece
 * of EOSMultiParts structure in Pa.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the minimum pressure (Pa)
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum pressure of the equation of state piece in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinPressure(EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double hmin = XLALSimNeutronStarEOSMultiPartsPieceMinEnthalpy(eos, piece_id);
    double pmin;
    pmin = XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometerized(hmin, eos, piece_id);
    pmin /= LAL_G_C4_SI;
    return pmin;
}

/**
 * @brief Returns the minimum pseudo-enthalpy of the ith equation of state piece
 * of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the minimum pseudo-enthalpy
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum pseudo-enthalpy of the equation of state piece.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinEnthalpy(EOSMultiParts * eos, int piece_id){
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMinEnthalpy(eos_part);
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the ith equation of state piece
 * of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the maximum pseudo-enthalpy
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The maximum pseudo-enthalpy of the equation of state piece.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMaxEnthalpy(EOSMultiParts * eos, int piece_id){
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos_part);
}

/**
 * @brief Returns the minimum acausal pseudo-enthalpy of the ith equation
 * of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the minimum acausal pseudo-enthalpy
 * for the equation of state piece number piece_id.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The minimum acausal pseudo-enthalpy of the equation of state piece.
 */
double XLALSimNeutronStarEOSMultiPartsPieceMinAcausalPseudoEnthalpy(EOSMultiParts * eos, int piece_id){
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    return XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(eos_part);
}


/* FUNCTIONS FOR INTERPOLATED EOS VALUES WITH GEOMETERIZED UNITS */

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


/**
 * @brief Returns the energy density in geometerized units (m^-2) at a given
 * pressure in geometerized units (m^-2), in the ith equation of state piece
 * of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the energy density interpolated at
 * a given pressure p, for the equation of state piece number piece_id.
 * @param p Pressure in geometerized units (m^-2)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressureGeometerized(double p,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double e;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    e = eos_part->e_of_p(p, eos_part);
    return e;
}

/**
 * @brief Returns the pseudo-enthalpy at a given pressure in geometerized units (m^-2),
 * in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the pseudo-enthalpy interpolated at
 * a given pressure p, for the equation of state piece number piece_id.
 * @param p Pressure in geometerized units (m^-2)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pseudo-enthalpy (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressureGeometerized(double p,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double h;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    h = eos_part->h_of_p(p, eos_part);
    return h;
}


/**
 * @brief Returns the pressure in geometerized units (m^-2) at a
 * given pseudo-enthalpy, in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the geometrized pressure
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpyGeometerized(double h,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double p;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    p = eos_part->p_of_h(h, eos_part); // TODO put some safeguards if interpolatino cannot be made ??
    return p;
}

/**
 * @brief Returns the energy density in geometerized units (m^-2) at a
 * given pseudo-enthalpy, in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the geometrized energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpyGeometerized(double
    h, EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double e;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    e = eos_part->e_of_h(h, eos_part);
    return e;
}


/**
 * @brief Returns the rest-mass density in geometerized units (m^-2) at a
 * given pseudo-enthalpy, in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the geometrized rest-mass density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The rest mass density in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpyGeometerized(double
    h, EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double rho;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    rho = eos_part->rho_of_h(h, eos_part);
    return rho;
}


/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometerized
 * units (m^-2), in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the gradient of the energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param p Pressure in geometerized units (m^-2).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressureGeometerized(double p,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double dedp;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    dedp = eos_part->dedp_of_p(p, eos_part);
    return dedp;
}


/**
 * @brief Returns the speed of sound in geometerized units (dimensionless)
 * at a given value of the pressure in geometerized units (m^-2),
 * in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the gradient of the energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The speed of sound in geometerized units (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSoundGeometerized(double h,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double v;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    v = eos_part->v_of_h(h, eos_part);
    return v;
}


/**
 * @brief Returns the energy density in geometerized units (m^-2) at a given
 * pressure in geometerized units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the energy density
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure in geometerized units (m^-2)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The energy density in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressureGeometerized(double p,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_pressure(p, eos);
    double e;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    e = eos_part->e_of_p(p, eos_part);
    return e;
}

/**
 * @brief Returns the pseudo-enthalpy at a given pressure in geometerized units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the pseudo-enthalpy
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure in geometerized units (m^-2)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The pseudo-enthalpy (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressureGeometerized(double p,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_pressure(p, eos);
    double h;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    h = eos_part->h_of_p(p, eos_part);
    return h;
}


/**
 * @brief Returns the pressure in geometerized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the pressure
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The pressure in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpyGeometerized(double h,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double p;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    p = eos_part->p_of_h(h, eos_part); // TODO put some safeguards if interpolatino cannot be made ??
    return p;
}

/**
 * @brief Returns the energy density in geometerized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the energy density
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The energy density in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpyGeometerized(double
    h, EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double e;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    e = eos_part->e_of_h(h, eos_part);
    return e;
}


/**
 * @brief Returns the rest-mass density in geometerized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * rest-mass density at the pseudo-enthalpy value h, for the correct
 * equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The rest mass density in geometerized units (m^-2).
*/
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpyGeometerized(double
    h, EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double rho;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    rho = eos_part->rho_of_h(h, eos_part);
    return rho;
}


/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometerized
 * units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the
 * energy density at the pressure value p, for the correct
 * equation of state piece.
 * @param p Pressure in geometerized units (m^-2).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressureGeometerized(double p,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_pressure(p, eos);
    double dedp;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    dedp = eos_part->dedp_of_p(p, eos_part);
    return dedp;
}


/**
 * @brief Returns the speed of sound in geometerized units (dimensionless)
 * at a given value of the pressure in geometerized units (m^-2).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * speed of sound at the pseudo-enthalpy value h, for the correct
 * equation of state piece.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The speed of sound in geometerized units (dimensionless).
*/
double XLALSimNeutronStarEOSMultiPartsSpeedOfSoundGeometerized(double h,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double v;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    v = eos_part->v_of_h(h, eos_part);
    return v;
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


/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa), in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the energy density
 * interpolated at a given pressure p, for the equation of state
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPressure(double p,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double e;
    p *= LAL_G_C4_SI;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    e = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(p, eos_part);
    e /= LAL_G_C4_SI;
    return e;
}


/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure (Pa), in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the pseudo-enthalpy
 * interpolated at a given pressure p, for the equation of state
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsPiecePseudoEnthalpyOfPressure(double p,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double h;
    p *= LAL_G_C4_SI;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    h = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(p, eos_part);
    return h;
}

/**
 * @brief Returns the pressure (Pa) at a given value of the dimensionless
 * pseudo-enthalpy, in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the pressure
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure (Pa).
 */
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double p;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos_part);
    p /= LAL_G_C4_SI;
    return p;
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * dimensionless pseudo-enthalpy, in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the energy density
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double e;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    e = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h,
        eos_part);
    e /= LAL_G_C4_SI;
    return e;
}

/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy, in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the rest-mass density
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSMultiPartsPieceRestMassDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts || piece_id < 0) XLAL_ERROR_REAL8(XLAL_EDOM);
    double rho;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    rho =
        XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(h,
        eos_part);
    rho /= LAL_G_C2_SI;
    return rho;
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure (Pa),
 * in the ith equation of state piece of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the gradient of the energy density
 * interpolated at a given pressure p, for the equation of state
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsPieceEnergyDensityDerivOfPressure(double p,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts) XLAL_ERROR_REAL8(XLAL_EDOM);
    double dedp;
    p *= LAL_G_C4_SI;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    dedp =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(p, eos_part);
    return dedp;
}

/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless), in the ith equation of state piece
 * of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the speed of sound
 * interpolated at a given pseudo-enthalpy h, for the equation of state
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSMultiPartsPieceSpeedOfSound(double h, EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts) XLAL_ERROR_REAL8(XLAL_EDOM);
    double v;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    v = XLALSimNeutronStarEOSSpeedOfSoundGeometerized(h, eos_part);
    v *= LAL_C_SI;
    return v;
}


/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3, in the ith equation of state piece
 * of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the pressure
 * interpolated at a given energy density e, for the equation of state
 * piece number piece_id.
 * @param e energy density in J/m^3
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfEnergyDensity(double e,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts) XLAL_ERROR_REAL8(XLAL_EDOM);
    double p;
    e *= LAL_G_C4_SI;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    p = eos_part->p_of_e(e, eos_part);
    p /= LAL_G_C4_SI;
    return p;
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3, in the ith equation of state piece
 * of EOSMultiParts structure.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1, and the function returns the pressure
 * interpolated at a given rest-mass density rho, for the equation of state
 * piece number piece_id.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @param piece_id Integer to the equation of state piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPiecePressureOfRestMassDensity(double rho,
    EOSMultiParts * eos, int piece_id)
{
    if (piece_id >= eos->number_of_parts) XLAL_ERROR_REAL8(XLAL_EDOM);
    double p;
    rho *= LAL_G_C2_SI;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, piece_id);
    p = eos_part->p_of_rho(rho, eos_part);
    p /= LAL_G_C4_SI;
    return p;
}


/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the energy density
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPressure(double p,
    EOSMultiParts * eos)
{
    p *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_pressure(p, eos);
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    double e;
    e = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(p, eos_part);
    e /= LAL_G_C4_SI;
    return e;
}

/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure (Pa).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the pseudo-enthalpy
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsPseudoEnthalpyOfPressure(double p,
    EOSMultiParts * eos)
{
    p *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_pressure(p, eos);
    double h;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    h = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(p, eos_part);
    return h;
}


/**
 * @brief Returns the pressure (Pa) at a given value of the dimensionless
 * pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the pressure
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The pressure (Pa).
 */
double XLALSimNeutronStarEOSMultiPartsPressureOfPseudoEnthalpy(double h,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double p;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos_part);
    p /= LAL_G_C4_SI;
    return p;
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the energy density
 * at the pseudo-enthalpy value h, for the correct equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSMultiPartsEnergyDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double e;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    e = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h,
        eos_part);
    e /= LAL_G_C4_SI;
    return e;
}


/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * rest-mass density at the pseudo-enthalpy value h,
 * for the correct equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSMultiPartsRestMassDensityOfPseudoEnthalpy(double h,
    EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double rho;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    rho =
        XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(h,
        eos_part);
    rho /= LAL_G_C2_SI;
    return rho;
}


/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure (Pa).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pressure value belongs to, then it interpolates the
 * gradient of the energy density with respect to the pressure
 * at the pressure value p, for the correct equation of state piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSMultiPartsEnergyDensityDerivOfPressure(double p,
    EOSMultiParts * eos)
{
    p *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_pressure(p, eos);
    double dedp;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    dedp =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(p, eos_part);
    return dedp;
}


/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless).
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * speed of sound at the pseudo-enthalpy value h, for the correct
 * equation of state piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSMultiPartsSpeedOfSound(double h, EOSMultiParts * eos)
{
    int item_piece = find_eos_piece_enthalpy(h, eos);
    double v;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    v = XLALSimNeutronStarEOSSpeedOfSoundGeometerized(h, eos_part);
    v *= LAL_C_SI;
    return v;
}


/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the energy density value belongs to, then it interpolates the
 * pressure at the energy density value e, for the correct
 * equation of state piece.
 * @param e energy density in J/m^3
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPressureOfEnergyDensity(double e,
    EOSMultiParts * eos)
{
    e *= LAL_G_C4_SI;
    int item_piece = find_eos_piece_energy_density(e, eos);
    double p;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    p = eos_part->p_of_e(e, eos_part);
    p /= LAL_G_C4_SI;
    return p;
}


/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3.
 * @details If the equation of state contains N phase transitions,
 * the EOSMultiParts structure contains N+1 pieces with ID number
 * from 0 to N-1. This function finds automatically which equation of state
 * piece the rest-mass density value belongs to, then it interpolates the
 * pressure at the rest-mass density value rho, for the correct
 * equation of state piece.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the equation of state EOSMultiParts structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSMultiPartsPressureOfRestMassDensity(double rho,
    EOSMultiParts * eos)
{
    rho *= LAL_G_C2_SI;
    int item_piece = find_eos_piece_rest_mass_density(rho, eos);
    double p;
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, item_piece);
    p = eos_part->p_of_rho(rho, eos_part);
    p /= LAL_G_C4_SI;
    return p;
}


/** @} */
/** @} */
