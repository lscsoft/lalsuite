/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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
#ifndef _PHMD_H
#define _PHMD_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup PHMD_h Header PHMD.h
 * \ingroup pkg_pulsarHough
 * \author Sintes, A. M.
 * \date 2001
 * \brief Conversion from peaks in a spectrum into a partial Hough map derivative
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/PHMD.h>
 * \endcode
 *
 * The Hough map is an histogram, thus additive. It can be seen as the sum of several
 * partial Hough maps constructed using just one periodogram, or equivalently, as
 * the sum of partial Hough map derivatives (\c phmd) and then integrating the
 * result.
 *
 * A  \c phmd can be represented by a set of borders, here called <em>left</em> and \e right.
 * They indicate the beginning and the end of the annuli.
 * The position of the so-called left borders should be marked with \f$+1\f$, and
 * the position of the right borders should be marked with \f$-1\f$ in the \c phmd.
 * To obtain a partial Hough map, one needs to integrate each row of the \c phmd
 * from left to right.
 *
 * The representation of a  \c phmd is simplified by considering
 * pointers to the borders in a pre-calculated look-up-table, plus some extra information about
 * their character and edge effects when clipping on a finite patch.
 *
 */
/*@{*/

/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

#include <lal/LUT.h>

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/**\name Error Codes */
/*@{*/
#define PHMDH_ENULL 1
#define PHMDH_ESIZE 2
#define PHMDH_ESZMM 4
#define PHMDH_EINT  6
#define PHMDH_ESAME 8
#define PHMDH_EFREQ 10
#define PHMDH_EVAL 12

#define PHMDH_MSGENULL "Null pointer"
#define PHMDH_MSGESIZE "Invalid input size"
#define PHMDH_MSGESZMM "Size mismatch"
#define PHMDH_MSGEINT  "Invalid interval"
#define PHMDH_MSGESAME "Input/Output data vectors are the same"
#define PHMDH_MSGEFREQ "Invalid frequency"
#define PHMDH_MSGEVAL  "Invalid value"
/*@}*/


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

/**  Hough Map derivative pixel type */
typedef REAL8 HoughDT; /* for weighted hough maps */
/* typedef CHAR  HoughDT; */
/*  typedef INT2  HoughDT; */


/**
 * \brief This structure stores the ``peak-gram''
 */
  typedef struct tagHOUGHPeakGram{
    INT2    timeIndex;  /**< The time index of the peak-gram */
    REAL8   deltaF;     /**< Frequency resolution: <tt>df=1/TCOH</tt> */
    UINT8   fBinIni;    /**< Frequency index of the first element of the spectrum covered by this peak-gram; it can be seen as an offset */
    UINT8   fBinFin;    /**< Frequency index of the last element of the spectrum covered by this peak-gram */
    UINT4   length;     /**< Number of peaks present in the peak-gram */
    INT4    *peak;      /**< The peak indices relative to \c fBinIni, i.e., the zero peak  corresponds to \c fBinIni */
  } HOUGHPeakGram;

/**
 * \brief This structure stores a partial Hough map derivative
 */
typedef struct tagHOUGHphmd{
  UINT8          fBin;  	/**< Frequency bin of this partial map derivative */
  UINT2          lengthLeft; 	/**< Exact number of \e Left borders */
  UINT2          lengthRight;	/**< Exact number of \e Right borders */
  UINT2          maxNBorders; 	/**< Maximun number of borders of each type (for memory allocation purposes),
                                 * i.e.\ length of <tt>*leftBorderP</tt> and <tt>*rightBorderP</tt>
                                 */
  HOUGHBorder    **leftBorderP; /**< Pointers to borders */
  HOUGHBorder    **rightBorderP;/**< Pointers to borders */
  UINT2          ySide;  	/**< number of elements of firstColumn */
  UCHAR          *firstColumn; 	/**< Number of elements of \c firstColumn */
  HoughDT        weight; 	/**< First column border, containing the edge effects  when clipping on a finite patch */
} HOUGHphmd;

/*
 * 11. Extern Global variables. (discouraged)
 */

/*
 * 12. Functions Declarations (i.e., prototypes).
 */

void LALHOUGHPeak2PHMD (LALStatus    *status,
			HOUGHphmd    *phmd,
			HOUGHptfLUT  *lut,
			HOUGHPeakGram *pg
			);

/*@}*/
#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _PHMD_H */
