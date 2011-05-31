
/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpACTD.h
 *
 * Author: McKechan, D. J. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \defgroup FindChirpACTD_h FindChirpACTD_h
 * \ingroup CBC_findchirp
 */

/**
\author McKechan, D. J. A.
\file
\ingroup FindChirpACTD_h

\brief Provides structures and functions for amplitude corrected time domain
templates using AmpCorPPN.

\heading{Synopsis}
\code
#include <lal/FindChirpACTD.h>
\endcode

*/


#ifndef _FINDCHIRPACTDH_H
#define _FINDCHIRPACTDH_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


NRCSID (FINDCHIRPACTDH, "$Id$");

/**\name Error Codes */ /*@{*/
#define FINDCHIRPACTDH_EQMAS 1
#define FINDCHIRPACTDH_MSGEQMAS "AmpCorPPN template equal mass"
/*@}*/

/** Define number of vectors, 6 for 0.5PN. */
#define NACTDVECS (3)

void
LALFindChirpACTDTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *theTmplt,
    FindChirpTmpltParams       *params
    );

void
LALFindChirpACTDNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpTmpltParams       *tmpltParams,
    FindChirpDataParams        *params
    );

void
LALFindChirpACTDFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

REAL4  XLALFindChirpACTDInnerProduct(
    COMPLEX8Vector *a,
    COMPLEX8Vector *b,
    COMPLEX8       *wtilde,
    REAL4           lower,
    REAL4           deltaT,
    UINT4           numPoints
    );


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPACTDH_H */
