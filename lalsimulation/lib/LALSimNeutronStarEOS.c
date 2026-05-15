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
struct tagEOSPiece {
    char name[LALNameLength];
    double pmin;
    double pmax;
    double hmax;
    double hmin;
    double hMinAcausal; /* Minimum pseudo-enthalpy at which EOS becomes acausal (speed of sound > 1) */
    double (*e_of_p) (double p, struct tagEOSPiece * myself);
    double (*h_of_p) (double p, struct tagEOSPiece * myself);
    double (*p_of_h) (double h, struct tagEOSPiece * myself);
    double (*e_of_h) (double h, struct tagEOSPiece * myself);
    double (*rho_of_h) (double h, struct tagEOSPiece * myself);
    double (*p_of_e) (double e, struct tagEOSPiece * myself);
    double (*p_of_rho) (double rho, struct tagEOSPiece * myself);
    double (*dedp_of_p) (double p, struct tagEOSPiece * myself);
    double (*v_of_h) (double h, struct tagEOSPiece * myself);
    void (*free) (struct tagEOSPiece * myself);
    int datatype;
    LALSimNeutronStarEOSData data;
};

struct tagLALSimNeutronStarEOS{
  char name[LALNameLength];
  int number_of_pieces;
  struct tagEOSPiece ** eos_piece;
  void (*free) (LALSimNeutronStarEOS * myself);
};


/* This function finds the id number of the piece EoS containing
 * the pressure value p, in a LALSimNeutronStarEOS structure.
 */
static int find_eos_piece_pressure(double p, LALSimNeutronStarEOS *eos)
{
    double pmin_eos = XLALSimNeutronStarEOSMinPressureGeometrized(eos);
    double pmax_eos = XLALSimNeutronStarEOSMaxPressureGeometrized(eos);
    if (p < pmin_eos || p > pmax_eos) {
        XLALPrintError("Pressure p = %.16e outside the interpolation range of EoS [%.16e:%.16e].\n", p, pmin_eos, pmax_eos);
        return XLAL_FAILURE;
    }
    for (int j = 0 ; j < XLALSimNeutronStarEOSNumberPieces(eos); j++){
        double pmin = XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(eos, j);
        double pmax = XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(eos, j);
        if (p >= pmin && p <= pmax) return j;
    }
    XLALPrintError("Pressure not found in EoS piece.\n");
    return XLAL_FAILURE;
}


/* This function finds the id number of the piece EoS containing
 * the pseudo-enthalpy value h, in a LALSimNeutronStarEOS structure.
 */
static int find_eos_piece_enthalpy(double h, LALSimNeutronStarEOS *eos)
{
    double hmin_eos = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos);
    double hmax_eos = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos);
    if (h < hmin_eos || h > hmax_eos) {
        XLALPrintError("Pseudo-enthalpy h = %.16e outside the interpolation range of EoS [%.16e:%.16e].\n", h, hmin_eos, hmax_eos);
        return XLAL_FAILURE;
    }
    for (int j = 0 ; j < XLALSimNeutronStarEOSNumberPieces(eos); j++){
        double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, j);
        double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, j);
        if (h >= hmin && h <= hmax) return j;
    }
    XLALPrintError("Enthalpy not found in EoS piece.\n");
    return XLAL_FAILURE;
}


/* This function finds the id number of the piece EoS containing
 * the rest-mass density value rho, in a LALSimNeutronStarEOS structure.
 */
static int find_eos_piece_rest_mass_density(double rho, LALSimNeutronStarEOS *eos)
{
    double rhomin_eos = XLALSimNeutronStarEOSMinRestMassDensityGeometrized(eos);
    double rhomax_eos = XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(eos);
    if (rho < rhomin_eos || rho > rhomax_eos) {
        XLALPrintError("Rest-mass density rho = %.16e outside the interpolation range of EoS [%.16e:%.16e].\n", rho, rhomin_eos, rhomax_eos);
        return XLAL_FAILURE;
    }
    for (int j = 0 ; j < XLALSimNeutronStarEOSNumberPieces(eos); j++){
        double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, j);
        double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, j);
        double rhomin = XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrizedPerPiece(hmin, eos, j);
        double rhomax = XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrizedPerPiece(hmax, eos, j);
        if (rho >= rhomin && rho <= rhomax) return j;
    }
    XLALPrintError("Rest mass density not found in EoS piece.\n");
    return XLAL_FAILURE;
}

/* This function finds the id number of the piece EoS containing
 * the energy density value e, in a LALSimNeutronStarEOS structure.
 */
static int find_eos_piece_energy_density(double e, LALSimNeutronStarEOS *eos)
{
    double emin_eos = XLALSimNeutronStarEOSMinEnergyDensityGeometrized(eos);
    double emax_eos = XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(eos);
    if (e < emin_eos || e > emax_eos) {
        XLALPrintError("Energy density e = %.16e outside the interpolation range of EoS [%.16e:%.16e].\n", e, emin_eos, emax_eos);
        return XLAL_FAILURE;
    }
    for (int j = 0 ; j < XLALSimNeutronStarEOSNumberPieces(eos); j++){
        double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, j);
        double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, j);
        double emin = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrizedPerPiece(hmin, eos, j);
        double emax = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrizedPerPiece(hmax, eos, j);
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
 * @name Creation and destruction routines
 * @{
 */

/**
 * @brief Allocated the memory associated with a pointer to an LALSimNeutronStarEOS structure.
 * This structure can accomodate phase transition. If the EOS contains N phase transitions
 * the EOS structure contains N+1 pieces with ID number from 0 to N.
 * @param nb_pieces Number of pieces of the EOS.
 * @return Pointer to the EOS structure to be freed.
 */
LALSimNeutronStarEOS * XLALCreateSimNeutronStarEOS(int nb_pieces)
{
    LALSimNeutronStarEOS *eos = NULL;
    eos = LALCalloc(1, sizeof(*eos));
    if (!eos) return NULL;
    eos->number_of_pieces = nb_pieces;
    eos->eos_piece = XLALCalloc(nb_pieces, sizeof(struct tagEOSPiece *));
    if (!eos->eos_piece) {
        LALFree(eos);
        return NULL;
    }
    return eos;
}


/**
 * @brief Frees the memory associated with a pointer to an LALSimNeutronStarEOS structure.
 * @param eos Pointer to the EOS structure to be freed.
 */
void XLALDestroySimNeutronStarEOS(LALSimNeutronStarEOS * eos)
{
    if(eos){
        if (eos->eos_piece) {
            for (int i = 0; i < eos->number_of_pieces; i++)
                if (eos->eos_piece[i]) eos->eos_piece[i]->free(eos->eos_piece[i]);
            LALFree(eos->eos_piece);
        }
        LALFree(eos);
    }

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



/* FUNCTIONS RELATED TO AN EQUATION OF STATE PIECE */

/**
 * @brief Returns the number of pieces in the Equation Of State (EOS)
 * structure LALSimNeutronStarEOS.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N, and the function returns N+1.
 * @param eos Pointer to the EOS structure.
 * @return The number of pieces in the EOS structure.
 */
int XLALSimNeutronStarEOSNumberPieces(LALSimNeutronStarEOS * eos)
{
    return eos->number_of_pieces;
}

/**
 * @brief Returns the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. This function returns a structure tagEOSPiece (only
 * internally defined) structure for the equation of piece with ID number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id ID number for the EOS piece.
 * @return Pointer of type LALSimNeutronStarEOSPiece, for the ith
 * piece of the LALSimNeutronStarEOS structure.
 */
static struct tagEOSPiece * XLALSimNeutronStarEOSSelectPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id < 0 || piece_id >= eos->number_of_pieces)
        XLAL_ERROR_NULL(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    return eos->eos_piece[piece_id];
}


/* FUNCTIONS FOR KEY EQUATION OF STATE VALUES */


/**
 * @brief The name of the Equation Of State (EOS).
 * @param[in] eos Pointer to the LALSimNeutronStarEOS structure.
 * @return Pointer to a string containing the name of the EOS.
 * @warning The pointer returned might be shallow and might be left
 * dangling if the @a eos structure is freed.
 */
char *XLALSimNeutronStarEOSName(LALSimNeutronStarEOS * eos)
{
    return XLALStringDuplicate(eos->name);
}

/**
 * @brief Returns the minimum pressure of the Equation Of State (EOS) piece in geometrized units m^-2.
 * @param eos Pointer to the EOS-piece structure.
 * @return The minimum pressure of the EOS piece in geometrized units m^-2.
 */
static double XLALSimNeutronStarEOSPieceMinPressureGeometrized(struct tagEOSPiece * eos)
{
    return eos->pmin;
}

/**
 * @brief Returns the minimum pressure of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in geometrized units m^-2.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum pressure (geometrized units)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum pressure of the EOS piece in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMinPressureGeometrized(eos_piece);
}

/**
 * @brief Returns the minimum pressure of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in Pa.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum pressure (Pa)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum pressure of the EOS piece in Pa.
 */
double XLALSimNeutronStarEOSMinPressurePerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the minimum pressure of the Equation Of State (EOS) in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The minimum pressure of the EOS in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSMinPressureGeometrized(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(eos, 0);
}

/**
 * @brief Returns the minimum pressure of the Equation Of State (EOS) in Pa.
 * @param eos Pointer to the EOS structure.
 * @return The minimum pressure of the EOS in Pa.
 */
double XLALSimNeutronStarEOSMinPressure(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinPressureGeometrized(eos)/LAL_G_C4_SI;
}


/**
 * @brief Returns the maximum pressure of the Equation Of State (EOS) piece in geometrized units m^-2.
 * @param eos Pointer to the EOS-piece structure.
 * @return The maximum pressure of the EOS piece in geometrized units m^-2.
 */
static double XLALSimNeutronStarEOSPieceMaxPressureGeometrized(struct tagEOSPiece * eos)
{
    return eos->pmax;
}

/**
 * @brief Returns the maximum pressure of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in geometrized units m^-2.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the maximum pressure (geometrized units)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum pressure of the EOS piece in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMaxPressureGeometrized(eos_piece);
}

/**
 * @brief Returns the maximum pressure of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in Pa.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the maximum pressure (Pa)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum pressure of the EOS piece in Pa.
 */
double XLALSimNeutronStarEOSMaxPressurePerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the maximum pressure of the Equation Of State (EOS) in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The maximum pressure of the EOS in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSMaxPressureGeometrized(LALSimNeutronStarEOS * eos)
{
    int number_of_pieces = XLALSimNeutronStarEOSNumberPieces(eos);
    return XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(eos, number_of_pieces-1);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSMaxPressureGeometrized */
double XLALSimNeutronStarEOSMaxPressureGeometerized(LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSPieceMaxPressureGeometrized");
    return XLALSimNeutronStarEOSMaxPressureGeometrized(eos);
}

/**
 * @brief Returns the maximum pressure of the Equation Of State (EOS) in Pa.
 * @param eos Pointer to the EOS structure.
 * @return The maximum pressure of the EOS in Pa.
 */
double XLALSimNeutronStarEOSMaxPressure(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMaxPressureGeometrized(eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the minimum pseudo-enthalpy of the Equation Of State (EOS) piece (dimensionless).
 * @param eos Pointer to the EOS-piece structure.
 * @return The minimum pseudo-enthalpy of the EOS piece (dimensionless).
 */
static double XLALSimNeutronStarEOSPieceMinPseudoEnthalpy(struct tagEOSPiece * eos)
{
    return eos->hmin;
}

/**
 * @brief Returns the minimum pseudo-enthalpy (dimensionless) of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum pseudo-enthalpy
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum pseudo-enthalpy of the EOS piece (dimensionless).
 */
double XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id){
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMinPseudoEnthalpy(eos_piece);
}

/**
 * @brief Returns the minimum pseudo-enthalpy of the Equation Of State (EOS) (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The minimum pseudo-enthalpy of the EOS (dimensionless).
 */
double XLALSimNeutronStarEOSMinPseudoEnthalpy(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, 0);
}

/**
 * @brief Returns the minimum pseudo-enthalpy at which Equation Of State (EOS) piece becomes acausal
 * (speed of sound > speed of light) (dimensionless).
 * @param eos Pointer to the EOS-piece structure.
 * @return The minimum pseudo-enthalpy at which EOS piece becomes acausal
 * (speed of sound > speed of light) (dimensionless).
 */
static double XLALSimNeutronStarEOSPieceMinAcausalPseudoEnthalpy(struct tagEOSPiece *
    eos)
{
    return eos->hMinAcausal;
}

/**
 * @brief Returns the minimum acausal pseudo-enthalpy (dimensionless) of the ith Equation Of State (EOS)
 * piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum acausal pseudo-enthalpy
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum acausal pseudo-enthalpy of the EOS piece (dimensionless).
 */
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id){
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMinAcausalPseudoEnthalpy(eos_piece);
}

/* Minimum pseudo-enthalpy at which EOS becomes acausal (speed of sound > 1).
 * If the EOS is always causal, return some large value hmax instead. */
static double eos_min_acausal_pseudo_enthalpy_tabular(double hmax,
    LALSimNeutronStarEOS * eos)
{
    size_t i;
    double h_im1, h_i;
    double v_im1, v_i;
    double m;   /* slope for linear interpolation */
    double hMinAcausal = hmax;  /* default large number for EOS that is always causal */
    int number_of_pieces = XLALSimNeutronStarEOSNumberPieces(eos);
    for (int j = 0; j < number_of_pieces; j++) {
        struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, j);
        h_im1 = exp(eos_piece->data.tabular->log_hdat[0]);
        v_im1 = eos_piece_v_of_h_tabular(h_im1, eos_piece);
        for (i = 1; i < eos_piece->data.tabular->ndat; i++) {
            h_i = exp(eos_piece->data.tabular->log_hdat[i]);
            v_i = eos_piece_v_of_h_tabular(h_i, eos_piece);
            if (v_i > 1.0) {
                /* solve vsound(h) = 1 */
                m = (v_i - v_im1) / (h_i - h_im1);
                hMinAcausal = h_im1 + (1.0 - v_im1) / m;
                return hMinAcausal;
            }
            h_im1 = h_i;
            v_im1 = v_i;
        }
    }
    return hMinAcausal;
}

/**
 * @brief Returns the minimum pseudo-enthalpy at which Equation Of State (EOS) becomes acausal
 * (speed of sound > speed of light) (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The minimum pseudo-enthalpy at which EOS becomes acausal
 * (speed of sound > speed of light) (dimensionless).
 */
double XLALSimNeutronStarEOSMinAcausalPseudoEnthalpy(LALSimNeutronStarEOS * eos)
{
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos);
    return eos_min_acausal_pseudo_enthalpy_tabular(hmax, eos);
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the Equation Of State (EOS) piece (dimensionless).
 * @param eos Pointer to the EOS-piece structure.
 * @return The maximum pseudo-enthalpy of the EOS piece (dimensionless).
 */
static double XLALSimNeutronStarEOSPieceMaxPseudoEnthalpy(struct tagEOSPiece * eos)
{
    return eos->hmax;
}

/**
 * @brief Returns the maximum pseudo-enthalpy (dimensionless) of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the maximum pseudo-enthalpy
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum pseudo-enthalpy of the EOS piece (dimensionless).
 */
double XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(LALSimNeutronStarEOS * eos, int piece_id){
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMaxPseudoEnthalpy(eos_piece);
}

/**
 * @brief Returns the maximum pseudo-enthalpy of the Equation Of State (EOS) (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The maximum pseudo-enthalpy of the EOS (dimensionless).
 */
double XLALSimNeutronStarEOSMaxPseudoEnthalpy(LALSimNeutronStarEOS * eos)
{
    int number_of_pieces = XLALSimNeutronStarEOSNumberPieces(eos);
    return XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, number_of_pieces-1);
}


/**
 * @brief Returns the minimum energy density of the Equation Of State (EOS)
 * piece in geometrized units m^-2.
 * @param eos Pointer to the EOS-piece structure.
 * @return The minimum energy density of the EOS piece in geometrized units m^-2.
 */
static double XLALSimNeutronStarEOSPieceMinEnergyDensityGeometrized(struct tagEOSPiece * eos)
{
    return eos->e_of_p(XLALSimNeutronStarEOSPieceMinPressureGeometrized(eos), eos);
}

/**
 * @brief Returns the minimum energy density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in geometrized units (m^-2).
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum energy density (m^-2)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum energy density of the EOS piece in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSMinEnergyDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMinEnergyDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the minimum energy density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in J/m^3.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum energy density (J/m^3)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum energy density of the EOS piece in J/m^3.
 */
double XLALSimNeutronStarEOSMinEnergyDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSMinEnergyDensityGeometrizedPerPiece(eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the minimum energy density of the Equation Of State (EOS) in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The minimum energy density of the EOS in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSMinEnergyDensityGeometrized(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinEnergyDensityGeometrizedPerPiece(eos, 0);
}

/**
 * @brief Returns the minimum energy density of the Equation Of State (EOS) in J m^-3.
 * @param eos Pointer to the EOS structure.
 * @return The minimum energy density of the EOS in J m^-3.
 */
double XLALSimNeutronStarEOSMinEnergyDensity(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinEnergyDensityGeometrized(eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the maximum energy density of the Equation Of State (EOS)
 * piece in geometrized units m^-2.
 * @param eos Pointer to the EOS-piece structure.
 * @return The maximum energy density of the EOS piece in geometrized units m^-2.
 */
static double XLALSimNeutronStarEOSPieceMaxEnergyDensityGeometrized(struct tagEOSPiece *  eos)
{
    return eos->e_of_p(XLALSimNeutronStarEOSPieceMaxPressureGeometrized(eos), eos);
}

/**
 * @brief Returns the maximum energy density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in geometrized units (m^-2).
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the maximum energy density (m^-2)
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum energy density of the EOS piece in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMaxEnergyDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the maximum energy density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in J/m^3.
 * @details If the EOS contains N phase transitions,
 * the LALSimNeutronStarEOS structure contains N+1 pieces with ID number
 * from 0 to N, and the function returns the maximum energy density J/m^3
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum energy density of the EOS piece in J/m^3.
 */
double XLALSimNeutronStarEOSMaxEnergyDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSMaxEnergyDensityGeometrizedPerPiece(eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the maximum energy density of the Equation Of State (EOS) in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The maximum energy density of the EOS in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(LALSimNeutronStarEOS * eos)
{
    int max_id_piece = XLALSimNeutronStarEOSNumberPieces(eos) - 1 ;
    return XLALSimNeutronStarEOSMaxEnergyDensityGeometrizedPerPiece(eos, max_id_piece);
}

/**
 * @brief Returns the maximum energy density of the Equation Of State (EOS) in J m^-3.
 * @param eos Pointer to the EOS structure.
 * @return The maximum energy density of the EOS in J m^-3.
 */
double XLALSimNeutronStarEOSMaxEnergyDensity(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the minimum rest-mass density of the Equation Of State (EOS)
 * piece in geometrized units m^-2.
 * @param eos Pointer to the EOS-piece structure.
 * @return The minimum rest-mass density of the EOS piece in geometrized units m^-2.
 */
static double XLALSimNeutronStarEOSPieceMinRestMassDensityGeometrized(struct tagEOSPiece *
    eos)
{
    return eos->rho_of_h(XLALSimNeutronStarEOSPieceMinPseudoEnthalpy(eos), eos);
}

/**
 * @brief Returns the minimum rest-mass density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in geometrized units.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum rest-mass density
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum rest-mass density of the EOS piece in geometrized units.
 */
double XLALSimNeutronStarEOSMinRestMassDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMinRestMassDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the minimum rest-mass density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in kg/m^3.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the minimum rest-mass density
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The minimum rest-mass density of the EOS piece in kg/m^3.
 */
double XLALSimNeutronStarEOSMinRestMassDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSMinRestMassDensityGeometrizedPerPiece(eos, piece_id)/LAL_G_C2_SI;
}

/**
 * @brief Returns the minimum rest-mass of the Equation Of State (EOS) in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The minimum rest-mass of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMinRestMassDensityGeometrized(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinRestMassDensityGeometrizedPerPiece(eos, 0);
}


/**
 * @brief Returns the minimum rest-mass of the Equation Of State (EOS) in kg/m^3.
 * @param eos Pointer to the EOS structure.
 * @return The minimum rest-mass of the EOS in kg/m^3.
 */
double XLALSimNeutronStarEOSMinRestMassDensity(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMinRestMassDensityGeometrized(eos)/LAL_G_C2_SI;
}

/**
 * @brief Returns the maximum rest-mass density of the Equation Of State (EOS)
 * piece in geometrized units m^-2.
 * @param eos Pointer to the EOS-piece structure.
 * @return The maximum rest-mass density of the EOS piece in geometrized units m^-2.
 */
static double XLALSimNeutronStarEOSPieceMaxRestMassDensityGeometrized(struct tagEOSPiece * eos)
{
    return eos->rho_of_h(XLALSimNeutronStarEOSPieceMaxPseudoEnthalpy(eos), eos);
}

/**
 * @brief Returns the maximum rest-mass density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in geometrized units.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the maximum rest-mass density
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum rest-mass density of the EOS piece in geometrized units.
 */
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrizedPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceMaxRestMassDensityGeometrized(eos_piece);
}

/**
 * @brief Returns the maximum rest-mass density of the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure in kg/m^3.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the maximum rest-mass density
 * for the EOS piece number piece_id.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The maximum rest-mass density of the EOS piece in kg/m^3.
 */
double XLALSimNeutronStarEOSMaxRestMassDensityPerPiece(LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSMaxRestMassDensityGeometrizedPerPiece(eos, piece_id)/LAL_G_C2_SI;
}

/**
 * @brief Returns the maximum rest-mass density of the Equation Of State (EOS)
 * in geometrized units m^-2.
 * @param eos Pointer to the EOS structure.
 * @return The maximum rest-mass of the EOS in geometrized units m^-2.
 */
double XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(LALSimNeutronStarEOS * eos)
{
    int max_id_piece = XLALSimNeutronStarEOSNumberPieces(eos) - 1 ;
    return XLALSimNeutronStarEOSMaxRestMassDensityGeometrizedPerPiece(eos, max_id_piece);
}

/**
 * @brief Returns the maximum rest-mass density of the Equation Of State (EOS)
 * in kg/m^3.
 * @param eos Pointer to the EOS structure.
 * @return The maximum rest-mass of the EOS in kg/m^3.
 */
double XLALSimNeutronStarEOSMaxRestMassDensity(LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(eos)/LAL_G_C2_SI;
}

/* FUNCTIONS FOR INTERPOLATED EOS VALUES */

/**
 * @brief Returns the dimensionless pseudo-enthalpy at a given pressure in
 * geometrized units (m^-2) for an Equation Of State (EOS) piece.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS-piece structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
static double XLALSimNeutronStarEOSPiecePseudoEnthalpyOfPressureGeometrized(double p,
    struct tagEOSPiece * eos)
{
    return eos->h_of_p(p, eos);
}

/**
 * @brief Returns the pseudo-enthalpy at a given pressure in geometrized units (m^-2),
 * in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pseudo-enthalpy interpolated at
 * a given pressure p, for the EOS piece number piece_id.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pseudo-enthalpy (dimensionless).
*/
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrizedPerPiece(double p,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double pmin = XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(eos, piece_id);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(eos, piece_id);
    if (p < pmin || p > pmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pressure p = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", p, pmin, pmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPiecePseudoEnthalpyOfPressureGeometrized(p, eos_piece);
}

/**
 * @brief Returns the pseudo-enthalpy (dimensionless) at a given pressure in Pa,
 * in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pseudo-enthalpy
 * interpolated at a given pressure p, for the EOS piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressurePerPiece(double p,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrizedPerPiece(p, eos, piece_id);
}

/**
 * @brief Returns the pseudo-enthalpy (dimensionless) at a given pressure in geometrized units (m^-2).
 * @details The function finds automatically which EOS
 * piece the pressure value belongs to, then it interpolates the pseudo-enthalpy
 * at the pressure value p, for the correct EOS piece.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The pseudo-enthalpy (dimensionless).
*/
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos)
{
    double pmin = XLALSimNeutronStarEOSMinPressureGeometrized(eos);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrized(eos);
    if (p < pmin || p > pmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pressure p = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", p, pmin, pmax);
    int item_piece = find_eos_piece_pressure(p, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrizedPerPiece(p, eos, item_piece);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized");
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(p, eos);
}

/**
 * @brief Returns the pseudo-enthalpy (dimensionless) at a given pressure (Pa).
 * @details The function finds automatically which EOS
 * piece the pressure value belongs to, then it interpolates the pseudo-enthalpy
 * at the pressure value p, for the correct EOS piece.
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @return The pseudo-enthalpy (dimensionless).
 */
double XLALSimNeutronStarEOSPseudoEnthalpyOfPressure(double p, LALSimNeutronStarEOS * eos)
{
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(p, eos);
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * pressure in geometrized units (m^-2) for an Equation Of State (EOS) piece.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS-piece structure.
 * @return The energy density in geometrized units (m^-2).
 */
static double XLALSimNeutronStarEOSPieceEnergyDensityOfPressureGeometrized(double p,
    struct tagEOSPiece * eos)
{
    return eos->e_of_p(p, eos);
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * pressure in geometrized units (m^-2), in the ith Equation Of State (EOS) piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the energy density interpolated at
 * a given pressure p, for the EOS piece number piece_id.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrizedPerPiece(double p,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double pmin = XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(eos, piece_id);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(eos, piece_id);
    if (p < pmin || p > pmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pressure p = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", p, pmin, pmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceEnergyDensityOfPressureGeometrized(p, eos_piece);
}

/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa),
 * in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the energy density
 * interpolated at a given pressure p, for the EOS
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The energy density (J m^3).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPressurePerPiece(double p,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrizedPerPiece(p, eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * pressure in geometrized units (m^-2).
 * @details The function finds automatically which EOS
 * piece the pressure value belongs to, then it interpolates the energy density
 * at the pressure value p, for the correct EOS piece.
 * @param p Pressure in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos)
{
    double pmin = XLALSimNeutronStarEOSMinPressureGeometrized(eos);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrized(eos);
    if (p < pmin || p > pmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pressure p = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", p, pmin, pmax);
    int item_piece = find_eos_piece_pressure(p, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrizedPerPiece(p, eos, item_piece);
}


/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized */
double XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized");
    return XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(p, eos);
}
/**
 * @brief Returns the energy density (J m^-3) at a given pressure (Pa).
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function finds automatically which EOS
 * piece the pressure value belongs to, then it interpolates the energy density
 * at the pressure value p, for the correct EOS piece.
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
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometrized
 * units (m^-2) for an Equation Of State (EOS) piece.
 * @param p Pressure in geometrized units (m^-2).
 * @param eos Pointer to the EOS-piece structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
static double XLALSimNeutronStarEOSPieceEnergyDensityDerivOfPressureGeometrized(double p,
    struct tagEOSPiece * eos)
{
    return eos->dedp_of_p(p, eos);
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometrized
 * units (m^-2), in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the gradient of the energy density
 * interpolated at a given pressure, for the EOS piece number piece_id.
 * @param p Pressure in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
*/
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrizedPerPiece(double p,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double pmin = XLALSimNeutronStarEOSMinPressureGeometrizedPerPiece(eos, piece_id);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrizedPerPiece(eos, piece_id);
    if (p < pmin || p > pmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pressure p = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", p, pmin, pmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceEnergyDensityDerivOfPressureGeometrized(p, eos_piece);
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in Pa,
 * in the ith EOS piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the gradient of the energy density
 * interpolated at a given pressure p, for the EOS
 * piece number piece_id.
 * @param p Pressure (Pa).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
 */
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressurePerPiece(double p,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    p *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrizedPerPiece(p,eos, piece_id);
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in geometrized
 * units (m^-2).
 * @details The function finds automatically which EOS
 * piece the pressure value belongs to, then it interpolates the
 * gradient of the energy density at the pressure value p, for the correct
 * EOS piece.
 * @param p Pressure in geometrized units (m^-2).
 * @param eos Pointer to the EOS structure.
 * @return The gradient of the energy density with respect to the pressure
 * (dimensionless).
*/
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(double p,
    LALSimNeutronStarEOS * eos)
{
    double pmin = XLALSimNeutronStarEOSMinPressureGeometrized(eos);
    double pmax = XLALSimNeutronStarEOSMaxPressureGeometrized(eos);
    if (p < pmin || p > pmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pressure p = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", p, pmin, pmax);
    int item_piece = find_eos_piece_pressure(p, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrizedPerPiece(p, eos, item_piece);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized */
double XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(double p, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized");
    return XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
}

/**
 * @brief Returns the gradient of the energy density with respect to the
 * pressure (dimensionless) at a given value of the pressure in Pa.
 * @details The function finds automatically which EOS
 * piece the pressure value belongs to, then it interpolates the
 * gradient of the energy density with respect to the pressure
 * at the pressure value p, for the correct EOS piece.
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
 * @brief Returns the pressure in geometrized units (m^-2) at a given value of
 * the dimensionless pseudo-enthalpy for an Equation Of State (EOS) piece.
 * @param h The value of the pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS-piece structure.
 * @return The pressure in geometrized units (m^-2).
 */
static double XLALSimNeutronStarEOSPiecePressureOfPseudoEnthalpyGeometrized(double h,
    struct tagEOSPiece * eos)
{
    return eos->p_of_h(h, eos);
}

/**
 * @brief Returns the pressure in geometrized units (m^-2) at a
 * given pseudo-enthalpy (dimensionless), in the ith Equation Of State (EOS)
 * piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the geometrized pressure
 * interpolated at a given pseudo-enthalpy h, for the EOS piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pressure in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrizedPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, piece_id);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, piece_id);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", h, hmin, hmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPiecePressureOfPseudoEnthalpyGeometrized(h, eos_piece);
}

/**
 * @brief Returns the pressure in Pa at a given pseudo-enthalpy,
 * in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pressure in Pa
 * interpolated at a given pseudo-enthalpy h, for the EOS piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pressure in Pa
*/
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrizedPerPiece(h, eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the pressure in geometrized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the pressure
 * at the pseudo-enthalpy value h, for the correct EOS piece.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the EOS structure.
 * @return The pressure in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(double h,
    LALSimNeutronStarEOS * eos)
{
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", h, hmin, hmax);
    int item_piece = find_eos_piece_enthalpy(h, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrizedPerPiece(h, eos, item_piece);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized");
    return XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos);
}

/**
 * @brief Returns the pressure (Pa) at a given value of the dimensionless
 * pseudo-enthalpy.
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the pressure
 * at the pseudo-enthalpy value h, for the correct EOS piece.
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
 * @brief Returns the energy density in geometrized units (m^-2) at a given
 * value of the dimensionless pseudo-enthalpy for an Equation Of State (EOS) piece.
 * @param h The value of the pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS-piece structure.
 * @return The energy density in geometrized units (m^-2).
 */
static double XLALSimNeutronStarEOSPieceEnergyDensityOfPseudoEnthalpyGeometrized(double
    h, struct tagEOSPiece * eos)
{
    return eos->e_of_h(h, eos);
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a given value of the
 * dimensionless pseudo-enthalpy, in the ith EOS piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the energy density
 * interpolated at a given pseudo-enthalpy h, for the EOS
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The energy density (m^-2).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrizedPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, piece_id);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, piece_id);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", h, hmin, hmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceEnergyDensityOfPseudoEnthalpyGeometrized(h, eos_piece);
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * pseudo-enthalpy (dimensionless), in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the energy density
 * interpolated at a given pseudo-enthalpy h, for the EOS
 * piece number piece_id.
 * @param h The value of the pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrizedPerPiece(h, eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the energy density in geometrized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the energy density
 * at the pseudo-enthalpy value h, for the correct EOS piece.
 * @param h Pseudo-enthalpy (dimensionless)
 * @param eos Pointer to the EOS structure.
 * @return The energy density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimNeutronStarEOS * eos)
{
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", h, hmin, hmax);
    int item_piece = find_eos_piece_enthalpy(h, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrizedPerPiece(h, eos, item_piece);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized");
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos);
}

/**
 * @brief Returns the energy density (J m^-3) at a given value of the
 * pseudo-enthalpy (dimensionless).
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the energy density
 * at the pseudo-enthalpy value h, for the correct EOS piece.
 * @param h The value of the pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The energy density (J m^-3).
 */
double XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h,eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the rest mass density in geometrized units (m^-2) at a given
 * value of the dimensionless pseudo-enthalpy for an Equation Of State (EOS) piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS-piece structure.
 * @return The rest mass density in geometrized units (m^-2).
 */
static double XLALSimNeutronStarEOSPieceRestMassDensityOfPseudoEnthalpyGeometrized(double
    h, struct tagEOSPiece * eos)
{
    return eos->rho_of_h(h, eos);
}

/**
 * @brief Returns the rest-mass density in geometrized units (m^-2) at a
 * given pseudo-enthalpy, in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the geometrized rest-mass density
 * interpolated at a given pseudo-enthalpy h, for the EOS piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The rest mass density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrizedPerPiece(double
    h, LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, piece_id);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, piece_id);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", h, hmin, hmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceRestMassDensityOfPseudoEnthalpyGeometrized(h, eos_piece);
}

/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy, in the ith EOS piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the rest-mass density
 * interpolated at a given pseudo-enthalpy h, for the EOS
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrizedPerPiece(h, eos, piece_id)/LAL_G_C2_SI;
}

/**
 * @brief Returns the rest-mass density in geometrized units (m^-2) at a
 * given pseudo-enthalpy.
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * rest-mass density at the pseudo-enthalpy value h, for the correct
 * EOS piece.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The rest mass density in geometrized units (m^-2).
*/
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(double
    h, LALSimNeutronStarEOS * eos)
{
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", h, hmin, hmax);
    int item_piece = find_eos_piece_enthalpy(h, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrizedPerPiece(h, eos, item_piece);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized");
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h, eos);
}

/**
 * @brief Returns the rest mass density (kg m^-3) at a given value of the
 * dimensionless pseudo-enthalpy.
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * rest-mass density at the pseudo-enthalpy value h,
 * for the correct EOS piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The rest mass density (kg m^-3), which is the number density of
 * baryons times the baryon rest mass.
 */
double XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(double h,
    LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h, eos)/LAL_G_C2_SI;
}

/**
 * @brief Returns the speed of sound in geometrized units (dimensionless)
 * at a given value of the pseudo-enthalpy (dimensionless) for an Equation Of State (EOS) piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS-piece structure.
 * @return The speed of sound in geometrized units (dimensionless).
 */
static double XLALSimNeutronStarEOSPieceSpeedOfSoundOfPseudoEnthalpyGeometrized(double h,
    struct tagEOSPiece * eos)
{
    return eos->v_of_h(h, eos);
}

/**
 * @brief Returns the speed of sound in geometrized units (dimensionless)
 * at a given value of the pressure in geometrized units (m^-2),
 * in the ith Equation Of State (EOS) piece of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the gradient of the energy density
 * interpolated at a given pseudo-enthalpy h, for the EOS piece number piece_id.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The speed of sound in geometrized units (dimensionless).
*/
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrizedPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpyPerPiece(eos, piece_id);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpyPerPiece(eos, piece_id);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", h, hmin, hmax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPieceSpeedOfSoundOfPseudoEnthalpyGeometrized(h, eos_piece);
}

/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless), in the ith EOS piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the speed of sound
 * interpolated at a given pseudo-enthalpy h, for the EOS
 * piece number piece_id.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyPerPiece(double h,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    return XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrizedPerPiece(h, eos, piece_id)*LAL_C_SI;
}

/**
 * @brief Returns the speed of sound in geometrized units (dimensionless)
 * at a given value of the pseudo-enthalpy (dimensionless).
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * speed of sound at the pseudo-enthalpy value h, for the correct
 * EOS piece.
 * @param h Pseudo-enthalpy (dimensionless).
 * @param eos Pointer to the EOS structure.
 * @return The speed of sound in geometrized units (dimensionless).
*/
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrized(double h,
    LALSimNeutronStarEOS * eos)
{
    double hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos);
    double hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos);
    if (h < hmin || h > hmax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input pseudo-enthalpy h = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", h, hmin, hmax);
    int item_piece = find_eos_piece_enthalpy(h, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrizedPerPiece(h, eos, item_piece);
}

/* DEPRECATED: old misspelling of XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrized */
double XLALSimNeutronStarEOSSpeedOfSoundGeometerized(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSSpeedOfSoundGeometrized");
    return XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrized(h, eos);
}

/**
 * @brief Returns the speed of sound (m s^-1) at a given value of the
 * pseudo-enthalpy (dimensionless).
 * @details The function finds automatically which EOS
 * piece the pseudo-enthalpy value belongs to, then it interpolates the
 * speed of sound at the pseudo-enthalpy value h, for the correct
 * EOS piece.
 * @param h The value of the dimensionless pseudo-enthalpy.
 * @param eos Pointer to the EOS structure.
 * @return The speed of sound (m s^-1).
 */
double XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpy(double h, LALSimNeutronStarEOS * eos)
{
    return XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpyGeometrized(h, eos)*LAL_C_SI;
}
/* DEPRECATED: old version of XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpy */
double XLALSimNeutronStarEOSSpeedOfSound(double h, LALSimNeutronStarEOS * eos)
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimNeutronStarEOSSpeedOfSound");
    return XLALSimNeutronStarEOSSpeedOfSoundOfPseudoEnthalpy(h, eos);
}

/**
 * @brief Returns the pressure in geometrized units (m^-2) at a given
 * energy density in geometrized units (m^-2).
 * @param e energy density in geometrized units (m^-2)
 * @param eos Pointer to the EOS structure.
 * @return The pressure in geometrized units (m^-2).
 */
static double XLALSimNeutronStarEOSPiecePressureOfEnergyDensityGeometrized(double e,
    struct tagEOSPiece * eos)
{
    return eos->p_of_e(e,eos);
}

/**
 * @brief Returns the pressure in geometrized units (m^-2) at a given
 * energy density in geometrized units (m^-2), in the ith EOS piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pressure
 * interpolated at a given energy density e, for the EOS
 * piece number piece_id.
 * @param e energy density in m^-2
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pressure in m^-2.
 */
double XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrizedPerPiece(double e,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double emin = XLALSimNeutronStarEOSMinEnergyDensityGeometrizedPerPiece(eos, piece_id);
    double emax = XLALSimNeutronStarEOSMaxEnergyDensityGeometrizedPerPiece(eos, piece_id);
    if (e < emin || e > emax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input energy density e = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", e, emin, emax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPiecePressureOfEnergyDensityGeometrized(e, eos_piece);
}

/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3, in the ith EOS piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pressure
 * interpolated at a given energy density e, for the EOS
 * piece number piece_id.
 * @param e energy density in J/m^3
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfEnergyDensityPerPiece(double e,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    e *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrizedPerPiece(e, eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the pressure in geometrized units (m^-2) at a given
 * energy density in J/m^3.
 * @details The function finds automatically which EOS
 * piece the energy density value belongs to, then it interpolates the
 * pressure at the energy density value e, for the correct
 * EOS piece. Note that if the user provides an energy density value that
 * falls within the phase transition jump, the value of the energu density is
 * automatically set to the maximum value of the EOS piece preceeding the phase
 * transition.
 * @param e energy density in J/m^3
 * @param eos Pointer to the EOS structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrized(double e,
    LALSimNeutronStarEOS * eos)
{
    double emin = XLALSimNeutronStarEOSMinEnergyDensityGeometrized(eos);
    double emax = XLALSimNeutronStarEOSMaxEnergyDensityGeometrized(eos);
    if (e < emin || e > emax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input energy density e = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", e, emin, emax);
    int number_of_pieces = XLALSimNeutronStarEOSNumberPieces(eos);
    if (number_of_pieces > 1){
        for (int i = 0; i < number_of_pieces-1; i++){
            struct tagEOSPiece * eos_low = XLALSimNeutronStarEOSSelectPiece(eos, i);
            struct tagEOSPiece * eos_high = XLALSimNeutronStarEOSSelectPiece(eos, i+1);
            if (e > XLALSimNeutronStarEOSPieceMaxEnergyDensityGeometrized(eos_low)
                && e < XLALSimNeutronStarEOSPieceMinEnergyDensityGeometrized(eos_high)){
                    e = XLALSimNeutronStarEOSPieceMaxEnergyDensityGeometrized(eos_low);
                }
        }
    }
    int item_piece = find_eos_piece_energy_density(e, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrizedPerPiece(e, eos, item_piece);
}

/**
 * @brief Returns the pressure in Pa at a given
 * energy density in J/m^3.
 * @details The function finds automatically which EOS
 * piece the energy density value belongs to, then it interpolates the
 * pressure at the energy density value e, for the correct
 * EOS piece. Note that if the user provides an energy density value that
 * falls within the phase transition jump, the value of the energu density is
 * automatically set to the maximum value of the EOS piece preceeding the phase
 * transition.
 * @param e energy density in J/m^3
 * @param eos Pointer to the EOS structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfEnergyDensity(double e,
    LALSimNeutronStarEOS * eos)
{
    e *= LAL_G_C4_SI;
    return XLALSimNeutronStarEOSPressureOfEnergyDensityGeometrized(e, eos)/LAL_G_C4_SI;
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the EOS structure.
 * @return The pressure in Pa.
 */
static double XLALSimNeutronStarEOSPiecePressureOfRestMassDensityGeometrized(double rho,
    struct tagEOSPiece * eos)
{
    return eos->p_of_rho(rho, eos);
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3, in the ith EOS piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pressure
 * interpolated at a given rest-mass density rho, for the EOS
 * piece number piece_id.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrizedPerPiece(double rho,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    if (piece_id >= eos->number_of_pieces || piece_id < 0)
        XLAL_ERROR_REAL8(XLAL_EDOM, "The ID piece number of LALSimNeutronStarEOS structure is incorrect.");
    double rhomin = XLALSimNeutronStarEOSMinRestMassDensityGeometrizedPerPiece(eos, piece_id);
    double rhomax = XLALSimNeutronStarEOSMaxRestMassDensityGeometrizedPerPiece(eos, piece_id);
    if (rho < rhomin || rho > rhomax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input rest-mass density e = %.16e is beyond the EOS piece interpolation range [%.16e:%.16e].", rho, rhomin, rhomax);
    struct tagEOSPiece * eos_piece = XLALSimNeutronStarEOSSelectPiece(eos, piece_id);
    return XLALSimNeutronStarEOSPiecePressureOfRestMassDensityGeometrized(rho, eos_piece);
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3, in the ith EOS piece
 * of LALSimNeutronStarEOS structure.
 * @details If the EOS contains N phase transitions that divide the EOS
 * into pieces, the LALSimNeutronStarEOS structure contains N+1 pieces
 * with ID number from 0 to N. The function returns the pressure
 * interpolated at a given rest-mass density rho, for the EOS
 * piece number piece_id.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the EOS structure.
 * @param piece_id Integer to the EOS piece ID number.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfRestMassDensityPerPiece(double rho,
    LALSimNeutronStarEOS * eos, int piece_id)
{
    rho *= LAL_G_C2_SI;
    return XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrizedPerPiece(rho, eos, piece_id)/LAL_G_C4_SI;
}

/**
 * @brief Returns the pressure in geometrized units (m^-2) at a given
 * rest-mass density in geometrized units.
 * @details The function finds automatically which EOS
 * piece the rest-mass density value belongs to, then it interpolates the
 * pressure at the rest-mass density value rho, for the correct
 * EOS piece.
 * @param rho rest-mass density in geometrized units
 * @param eos Pointer to the EOS structure.
 * @return The pressure in geometrized units (m^-2).
 */
double XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrized(double rho,
    LALSimNeutronStarEOS * eos)
{
    double rhomin = XLALSimNeutronStarEOSMinRestMassDensityGeometrized(eos);
    double rhomax = XLALSimNeutronStarEOSMaxRestMassDensityGeometrized(eos);
    if (rho < rhomin || rho > rhomax)
        XLAL_ERROR_REAL8(XLAL_EDOM,
                         "Input rest-mass density e = %.16e is beyond the EOS interpolation range [%.16e:%.16e].", rho, rhomin, rhomax);
    int number_of_pieces = XLALSimNeutronStarEOSNumberPieces(eos);
    if (number_of_pieces > 1){
        for (int i = 0; i < number_of_pieces-1; i++){
            struct tagEOSPiece * eos_low = XLALSimNeutronStarEOSSelectPiece(eos, i);
            struct tagEOSPiece * eos_high = XLALSimNeutronStarEOSSelectPiece(eos, i+1);
            if (rho > XLALSimNeutronStarEOSPieceMaxRestMassDensityGeometrized(eos_low)
                && rho < XLALSimNeutronStarEOSPieceMinRestMassDensityGeometrized(eos_high)){
                    rho = XLALSimNeutronStarEOSPieceMaxRestMassDensityGeometrized(eos_low);
                }
        }
    }
    int item_piece = find_eos_piece_rest_mass_density(rho, eos);
    if (item_piece < 0) XLAL_ERROR_REAL8(XLAL_EFUNC);
    return XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrizedPerPiece(rho, eos, item_piece);
}

/**
 * @brief Returns the pressure in Pa at a given
 * rest-mass density in kg/m^3.
 * @details The function finds automatically which EOS
 * piece the rest-mass density value belongs to, then it interpolates the
 * pressure at the rest-mass density value rho, for the correct
 * EOS piece.
 * @param rho rest-mass density in kg/m^3
 * @param eos Pointer to the EOS structure.
 * @return The pressure in Pa.
 */
double XLALSimNeutronStarEOSPressureOfRestMassDensity(double rho,
    LALSimNeutronStarEOS * eos)
{
    rho *= LAL_G_C2_SI;
    return XLALSimNeutronStarEOSPressureOfRestMassDensityGeometrized(rho, eos)/LAL_G_C4_SI;
}

/** @} */
/** @} */
