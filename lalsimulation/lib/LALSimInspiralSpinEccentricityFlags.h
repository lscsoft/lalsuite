/*
 * Copyright (C) 2024 Khun Sang Phukon, Nathan Johnson-McDaniel, Anuradha Gupta
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


#ifndef _LALSIMINSPIRALSPINECCENTRICITYFLAGS_H
#define _LALSIMINSPIRALSPINECCENTRICITYFLAGS_H


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/* Enumeration for Eccentricity order */
typedef enum tagLALSimInspiralSpinEccPolyEccOrder{
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_0  = 0,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_2  = 2,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_4  = 4,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_6  = 6,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_8  = 8,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_10  = 10,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_12  = 12,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_14  = 14,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_16  = 16,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_18  = 18,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_20  = 20,
    LAL_SIM_INSPIRAL_SPINECC_ECC_ORDER_ALL  = -1
}LALSimInspiralSpinEccPolyEccOrder;

/* Enumeration for Eccentricity enhancement function */
typedef enum tagEnhancementFunction{
    EF_Arun_etal,
    EF_LoutrelYunes_SuperAsym,
    EF_LoutrelYunes_HyperAsym
}EnhancementFunction;

/*
 * Evolution type
 * FullSeries: Store all intermediate values in a timeseries
 * FinalValue: Store only initial value and final value in a timeseries
 * */
typedef enum tagEvolutionType{
    FullSeries,
    FinalValue
}EvolutionType;

/*
 * Eccentricity evolution 2PN spin flag
 * Exclude: Exclude the 2PN spin terms in the eccentricity evolution, so zero eccentricity orbits stay at zero eccentricity
 * Include: Include the 2PN spin terms in the eccentricity evolution, so the eccentricity time derivative doesn't vanish at zero eccentricity
 * */
typedef enum tagEccEvol2PNSpinFlag{
    Exclude,
    Include
}EccEvol2PNSpinFlag;

/*
 * Eccentricity evolution output flag
 * ErrorOnStop: Raise an error if a stopping condition besides the frequency bound is triggered, except for the energy and angular momentum flux, xbar, and
 *              second time derivative of omega conditions when evolving as far forward as possible
 * WarningOnStop: Only print a warning if a stopping condition besides the frequency bound is triggered, and return the output up to (FullSeries) or at
 *                (FinalValue) the final timestep before the stopping condition is triggered
 * */
typedef enum tagEccEvolAlwaysOutput{
    ErrorOnStop,
    WarningOnStop
}EccEvolAlwaysOutput;


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALSPINECCENTRICITYFLAGS_H */
