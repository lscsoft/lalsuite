/*
*  Copyright (C) 2007 Badri Krishnan, Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: RngMedBias.h
 *
 * Authors:  Krishnan, B.
 *
 *
 * History:   Created by Krishnan Feb 24, 2004
 *
 *
 *-----------------------------------------------------------------------
 */

#ifndef _RNGMEDBIAS_H
#define _RNGMEDBIAS_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
\addtogroup RngMedBias_h
\author Krishnan, B., Itoh, Y.

\brief Routine for finding bias in median for exponential distribution

\heading{Synopsis}

\code
#include <lal/RngMedBias.h>
\endcode

*/
/*@{*/

/**\name Error Codes */
/*@{*/
#define RNGMEDBIASH_ENULL 1		/**< Null pointer */
#define RNGMEDBIASH_EVAL 5		/**< Invalid value */
/*@}*/

/** \cond DONT_DOXYGEN */
#define RNGMEDBIASH_MSGENULL "Null pointer"
#define RNGMEDBIASH_MSGEVAL  "Invalid value"
/** \endcond */

/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */



/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */
REAL8 XLALRngMedBias ( INT4 blkSize );


// ------------------------------ obsolte and deprecated LAL-interface functions --------------------
void LALRngMedBias (LALStatus   *status,
                    REAL8       *biasFactor,
                    INT4        blkSize
                    );

/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _RNGMEDBIAS_H */








