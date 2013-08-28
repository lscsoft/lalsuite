/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#ifndef _ZPGFILTER_H
#define _ZPGFILTER_H

#include <lal/LALStdlib.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup ZPGFilter_h
 * \author Creighton, T. D.
 *
 * \brief Provides routines to manipulate ZPG filters.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/ZPGFilter.h>
 * \endcode
 *
 * This header covers routines that create, destroy, and
 * transform objects of type <tt>\<datatype\>ZPGFilter</tt>, where
 * <tt>\<datatype\></tt> is either \c COMPLEX8 or \c COMPLEX16.
 * Generically, these data types can be used to store any rational
 * complex function in a factored form.  Normally this function is a
 * filter response, or "transfer function" \f$T(z)\f$, expressed in terms
 * of a complex frequency parameter \f$z=\exp(2\pi if\Delta t)\f$, where
 * \f$\Delta t\f$ is the sampling interval.  The rational function is
 * factored as follows:
 * \f[
 * T(f) = g\times\frac{\prod_k (z-a_k)}{\prod_l (z-b_l)}
 * \f]
 * where \f$g\f$ is the gain, \f$a_k\f$ are the (finite) zeros, and \f$b_l\f$ are the
 * (finite) poles.  It should be noted that rational functions always
 * have the same number of zeros as poles if one includes the point
 * \f$z=\infty\f$; any excess in the number of finite zeros or poles in the
 * rational expression simply indicates that there is a corresponding
 * pole or zero of that order at infinity.  It is also worth pointing out
 * that the "gain" is just the overall prefactor of this rational
 * function, and is not necessarily equal to the actual gain of the
 * transfer function at any particular frequency.
 *
 * Another common complex frequency space is the \f$w\f$-space, obtained
 * from the \f$z\f$-space by the bilinear transformation:
 * \f[
 * w = i\left(\frac{1-z}{1+z}\right) = \tan(\pi f\Delta t) , \quad
 * z = \frac{1+iw}{1-iw} \; .
 * \f]
 * Other variables can also be used to represent the complex frequency
 * plane.  The <tt>\<datatype\>ZPGFilter</tt> structure can be used to
 * represent the transfer function in any of these spaces by transforming
 * the coordinates of the zeros and poles, and incorporating any residual
 * factors into the gain.  Care must be taken to include any zeros or
 * poles that are brought in from infinity by the transformation, and to
 * remove any zeros or poles which were sent to infinity.  Thus the
 * number of zeros and poles of the <tt>\<datatype\>ZPGFilter</tt> is not
 * necessarily constant under transformations!  Routines invoking the
 * <tt>\<datatype\>ZPGFilter</tt> data types should document which complex
 * variable is assumed.
 *
 */
/*@{*/

/**
 * @{
 * \defgroup CreateZPGFilter_c Module CreateZPGFilter.c
 * \defgroup DestroyZPGFilter_c Module DestroyZPGFilter.c
 * \defgroup BilinearTransform_c Module BilinearTransform.c
 * @}
 */

/** \name Error Codes */
/*@{*/
#define ZPGFILTERH_ENUL 1	/**< Unexpected null pointer in arguments */
#define ZPGFILTERH_EOUT 2	/**< Output handle points to a non-null pointer */
#define ZPGFILTERH_EMEM 3	/**< Memory allocation error */
#define ZPGFILTERH_EBAD 4	/**< Bad filter parameters */
/*@}*/
/*@}*/

#define ZPGFILTERH_MSGENUL "Unexpected null pointer in arguments"
#define ZPGFILTERH_MSGEOUT "Output handle points to a non-null pointer"
#define ZPGFILTERH_MSGEMEM "Memory allocation error"
#define ZPGFILTERH_MSGEBAD "Bad filter parameters"


/* ---------- Function prototypes. ---------- */
COMPLEX8ZPGFilter *XLALCreateCOMPLEX8ZPGFilter(INT4 numZeros, INT4 numPoles);
COMPLEX16ZPGFilter *XLALCreateCOMPLEX16ZPGFilter(INT4 numZeros, INT4 numPoles);
void XLALDestroyCOMPLEX8ZPGFilter( COMPLEX8ZPGFilter *filter );
void XLALDestroyCOMPLEX16ZPGFilter( COMPLEX16ZPGFilter *filter );
int XLALWToZCOMPLEX8ZPGFilter( COMPLEX8ZPGFilter *filter );
int XLALWToZCOMPLEX16ZPGFilter( COMPLEX16ZPGFilter *filter );

void
LALCreateCOMPLEX8ZPGFilter( LALStatus         *status,
			    COMPLEX8ZPGFilter **output,
			    INT4              numZeros,
			    INT4              numPoles );

void
LALCreateCOMPLEX16ZPGFilter( LALStatus          *status,
			     COMPLEX16ZPGFilter **output,
			     INT4               numZeros,
			     INT4               numPoles );

void
LALDestroyCOMPLEX8ZPGFilter( LALStatus         *status,
			     COMPLEX8ZPGFilter **input );

void
LALDestroyCOMPLEX16ZPGFilter( LALStatus          *status,
			      COMPLEX16ZPGFilter **input );

void
LALWToZCOMPLEX8ZPGFilter( LALStatus         *status,
			  COMPLEX8ZPGFilter *filter );

void
LALWToZCOMPLEX16ZPGFilter( LALStatus          *status,
			   COMPLEX16ZPGFilter *filter );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _ZPGFILTER_H */
