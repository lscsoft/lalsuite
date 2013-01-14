
/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpACTD.h
 *
 * Author: McKechan, D. J. A.
 *
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPACTDH_H
#define _FINDCHIRPACTDH_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup FindChirpACTD_h
 \author McKechan, D. J. A.

\brief Provides structures and functions for amplitude corrected time domain
templates using AmpCorPPN.

\heading{Synopsis}
\code
#include <lal/FindChirpACTD.h>
\endcode

*/
/*@{*/

/**\name Error Codes */
/*@{*/
#define FINDCHIRPACTDH_EQMAS 1	/**< AmpCorPPN template equal mass */
/*@}*/

/** \cond DONT_DOXYGEN */
#define FINDCHIRPACTDH_MSGEQMAS "AmpCorPPN template equal mass"
/** \endcond */


/** Define number of vectors, 6 for 0.5PN. */
#define NACTDVECS (3)

/*@}*/


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPACTDH_H */
