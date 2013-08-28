/*
*  Copyright (C) 2007 Bernd Machenschalk, Patrick Brady, Reinhard Prix
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

#ifndef _LALRUNNINGMEDIAN_H
#define _LALRUNNINGMEDIAN_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


/**
 * \addtogroup LALRunningMedian_h
 * \author Somya D. Mohanty, B. Machenschalk
 *
 * \brief Provides routines to efficiently calculate the running median
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALRunningMedian.h>
 * \endcode
 *
 * This header covers routines to efficiently calculate the
 * running median of REAL4 and REAL8 sequences
 *
 * The routine <tt>LALDRunningMedian()</tt> calculates the running medians of a
 * REAL8Sequence. The routine <tt>LALSRunningMedian()</tt> does the same for a REAL4Sequence.
 * \c input ist a REAL4/REAL8Sequence containing the input array, \c blocksize
 * is the length of the block the medians are calculated of.
 * With n being the lenght of the input array and b being the blocksize,
 * the medians array must be a REAL4/REAL8 sequence of length (n-b+1).
 * <tt>LALDRunningMedian2()</tt> and <tt>LALSRunningMedian2()</tt> are a
 * different implentation of the same algorithm. It should behave exactly like
 * <tt>LALDRunningMedian()</tt>, but has proven to be a
 * little faster and more stable. Check if it works for you.
 *
 * ### Algorithm ###
 *
 * For a detailed description of the algorithm see the
 * LIGO document T-030168-00-D, Somya D. Mohanty:
 * Efficient Algorithm for computing a Running Median
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define LALRUNNINGMEDIANH_EMALOC1 1		/**< Could not allocate indexblock */
#define LALRUNNINGMEDIANH_EMALOC2 2		/**< Could not allocate checks */
#define LALRUNNINGMEDIANH_EMALOC3 3		/**< Could not allocate checks4shift */
#define LALRUNNINGMEDIANH_EMALOC4 4		/**< Could not allocate nodeaddresses */
#define LALRUNNINGMEDIANH_EMALOC5 5		/**< Could not aloocate first node */
#define LALRUNNINGMEDIANH_EMALOC6 6		/**< Could not allocate node */
#define LALRUNNINGMEDIANH_ECV 7			/**< Could not create output vector (LALCreateVector() failed) */
#define LALRUNNINGMEDIANH_ENULL 8		/**< Invalid input: NULL pointer. */
#define LALRUNNINGMEDIANH_EZERO 9		/**< Invalid input: block length must be \>2 */
#define LALRUNNINGMEDIANH_ELARGE 10		/**< Invalid input: block length larger than imput length */
#define LALRUNNINGMEDIANH_EIMED 11		/**< Invalid input: wrong size of median array */
/*@}*/

/** \cond DONT_DOXYGEN */
#define LALRUNNINGMEDIANH_MSGEMALOC1 "Could not allocate indexblock"
#define LALRUNNINGMEDIANH_MSGEMALOC2 "Could not allocate checks"
#define LALRUNNINGMEDIANH_MSGEMALOC3 "Could not allocate checks4shift"
#define LALRUNNINGMEDIANH_MSGEMALOC4 "Could not allocate nodeaddresses"
#define LALRUNNINGMEDIANH_MSGEMALOC5 "Could not aloocate first node"
#define LALRUNNINGMEDIANH_MSGEMALOC6 "Could not allocate node"
#define LALRUNNINGMEDIANH_MSGECV     "Could not create output vector (LALCreateVector failed)"
#define LALRUNNINGMEDIANH_MSGENULL   "Invalid input: NULL pointer."
#define LALRUNNINGMEDIANH_MSGEZERO   "Invalid input: block length must be >2"
#define LALRUNNINGMEDIANH_MSGELARGE  "Invalid input: block length larger than imput length"
#define LALRUNNINGMEDIANH_MSGEIMED   "Invalid input: wrong size of median array"
/** \endcond */

/* Structures. */

/**
 * This is the parameter structure for the LALRunningMedian functions.
 * Currently the only parameter supported is the blocksize, the number
 * of elements a single median is calculated from. The current
 * implementation requires the blocksize to be \< 2.
 */
typedef struct tagLALRunningMedianPar
{
  UINT4 blocksize;	/**< the number of elements a single median is calculated from */
}
LALRunningMedianPar;


/* Function prototypes. */

/** See LALRunningMedian_h for documentation */
void
LALDRunningMedian( LALStatus *status,
		   REAL8Sequence *medians,
		   const REAL8Sequence *input,
		   LALRunningMedianPar param);

/** See LALRunningMedian_h for documentation */
void
LALSRunningMedian( LALStatus *status,
		   REAL4Sequence *medians,
		   const REAL4Sequence *input,
		   LALRunningMedianPar param);

/** See LALRunningMedian_h for documentation */
void
LALDRunningMedian2( LALStatus *status,
		    REAL8Sequence *medians,
		    const REAL8Sequence *input,
		    LALRunningMedianPar param);

/** See LALRunningMedian_h for documentation */
void
LALSRunningMedian2( LALStatus *status,
		    REAL4Sequence *medians,
		    const REAL4Sequence *input,
		    LALRunningMedianPar param);

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _LALRUNNINGMEDIAN_H */
