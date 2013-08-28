/*
*  Copyright (C) 2007 Jolien Creighton, Robert Adam Mercer, Tania Regimbau
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
 * \author Regimbau Tania
 * \addtogroup AstroOmega_h
 *
 * \brief Compute the energy density spectrum of stochastic backgrounds produced
 * by cosmological population of astrophysical sources.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/AstroOmega.h>
 * \endcode
 *
 * \heading{Error conditions}
 * the errors that may occur in this module are integration errors already defined in Integrate.h
 *
 * \heading{Structures}
 * These are function pointers corresponding to the spectral energy density of a single source.
 * \code
 * typedef void (REAL8LALSDensity) (REAL8 *output, REAL8 input);
 * \endcode
 * These are input structures corresponding to the model parameters (the cosmological model parameters and the source model parameters)
 *
 * cosmological model parameters:
 *
 * \code
 * typedef struct
 * tagAstroOmegaCosmoParams
 * {
 * REAL8   ho; Hubble parameter
 * REAL8   density_matter; density parameter of matter
 * REAL8   density_vacuum; density parameter of vacuum
 * REAL8   density_k; density parameter of curvature
 * }
 * AstroOmegaCosmoParams;
 * \endcode
 *
 * source parameters
 *
 * \code
 * typedef struct
 * tagAstroOmegaSourceParams
 * {
 * REAL8LALSDensity   *SDensitySource; single spectral energy density
 * REAL8              numax; frequency cutoff in the source frame
 * REAL8              lambda; mass fraction of source progenitors expressed in inverse solar masses.
 * }
 * AstroOmegaSourceParams;
 * \endcode
 *
 * model parameters (cosmological + source)
 *
 * \code
 * typedef struct
 * tagAstroOmegaParams
 * {
 * AstroOmegaCosmoParams          cosmoparams;
 * AstroOmegaSourceParams         sourceparams;
 * void                           *extraparams;
 * }
 * AstroOmegaParams;
 * \endcode
 *
 * @{
 */

#ifndef _ASTROOMEGA_H
#define _ASTROOMEGA_H
#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>

#ifdef __cplusplus
extern "C" {
#endif

/*type corresponding to the spectral energy density of a single source*/
typedef void (REAL8LALSDensity) (REAL8 *output, REAL8 input);

/*MODEL PARAMETERS*/

/*cosmological model*/
typedef struct
tagAstroOmegaCosmoParams
 {
   REAL8   ho;
   REAL8   density_matter;
   REAL8   density_vacuum;
   REAL8   density_k;
 }
AstroOmegaCosmoParams;
/*source model*/
/*in the general case, the user should define previously the single spectral energy density*/
typedef struct
tagAstroOmegaSourceParams
 {
   REAL8LALSDensity   *SDensitySource;
   REAL8              numax;
   REAL8              lambda;
 }
AstroOmegaSourceParams;


typedef struct
tagAstroOmegaParams
 {
   AstroOmegaCosmoParams          cosmoparams;
   AstroOmegaSourceParams         sourceparams;
   void                           *extraparams;
 }
AstroOmegaParams;



/*functions returning $\Omega _{gw}(\nu _{o})$*/

void
LALAstroOmega (
    LALStatus    *status,
    REAL8        *result,
    REAL8         nu,
    void         *params
    );


/** @} */

#ifdef __cplusplus
}
#endif

#endif /* _ASTROOMEGA_H */
