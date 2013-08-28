/*
*  Copyright (C) 2007 Stas Babak, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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

#ifndef _LALNOISEMODELS_H
#define _LALNOISEMODELS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/RealFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup LALNoiseModels_h
 * \author Sathyaprakash, B.S.
 *
 * \brief Header file for model noise generation codes.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/LALNoiseModels.h>
 * \endcode
 *
 * This header file covers routines that are used in
 * synthetic background noise  expected in various
 * detectors and signals with random parameters in background noise.
 *
 */
/* @{ */
/**\name Error Codes */ /*@{*/
#define LALNOISEMODELSH_ENULL 	1	/**< Arguments contained an unexpected null pointer */
#define LALNOISEMODELSH_EMEM 	2	/**< Memory allocation error */
#define LALNOISEMODELSH_ECHOICE 4	/**< Invalid choice for an input parameter */
#define LALNOISEMODELSH_EDIV0 	8	/**< Division by zero */
#define LALNOISEMODELSH_ESIZE 	16	/**< Invalid input size */
/*@} */

/** \cond DONT_DOXYGEN */
#define LALNOISEMODELSH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALNOISEMODELSH_MSGEMEM "Memory allocation error"
#define LALNOISEMODELSH_MSGECHOICE "Invalid choice for an input parameter"
#define LALNOISEMODELSH_MSGEDIV0 "Division by zero"
#define LALNOISEMODELSH_MSGESIZE "Invalid input size"
/** \endcond */

    enum enumDetector
    {
        geo,
        ligo,
        tama,
        virgo
    };

    typedef struct
            tagAddVectorsIn
            {
                REAL4Vector *v1;
                REAL4Vector *v2;
                REAL8       a1;
                REAL8       a2;
            }
    AddVectorsIn;


    typedef struct
            tagStatsREAL4VectorOut
            {
                REAL8 mean;
                REAL8 var;
                REAL8 stddev;
                REAL8 min;
                REAL8 max;
            }
    StatsREAL4VectorOut;



    /* Function prototypes */

    void
            LALNoiseSpectralDensity
            (
             LALStatus   *status,
             REAL8Vector *psd,
             void        (*NoisePsd)(LALStatus *status, REAL8 *shf, REAL8 f),
             REAL8       f
            );

    void
            LALEGOPsd
            (
             LALStatus *status,
             REAL8     *shf,
             REAL8     x
            );


    void
            LALGEOPsd
            (
             LALStatus *status,
             REAL8     *shf,
             REAL8     x
            );

    void
            LALAdvLIGOPsd
            (
             LALStatus *status,
             REAL8     *shf,
             REAL8     x
            );

    REAL8
            XLALLIGOIPsd
            (
             REAL8     f
            );

    void
            LALLIGOIPsd
            (
             LALStatus *status,
             REAL8     *shf,
             REAL8     x
            );

    void
            LALTAMAPsd
            (
             LALStatus *status,
             REAL8     *shf,
             REAL8     x
            );

    void
            LALVIRGOPsd
            (
             LALStatus *status,
             REAL8     *shf,
             REAL8     x
            );


    void
            LALColoredNoise
            (
             LALStatus   *status,
             REAL4Vector *noisy,
             REAL8Vector psd
            );

    void
            LALAddVectors
            (
             LALStatus *status,
             REAL4Vector *vector,
             AddVectorsIn in);


    void
            LALStatsREAL4Vector
            (
             LALStatus *status,
             StatsREAL4VectorOut *out,
             REAL4Vector *vector
            );


    REAL8 XLALETBPsd( REAL8 f );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _LALNOISEMODELS_H */
