/*
*  Copyright (C) 2007 Jolien Creighton, Drew Keppel
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

#ifndef _INTEGRATE_H
#define _INTEGRATE_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * \addtogroup Integrate_h
 *
 * \brief Integrates a function.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Integrate.h>
 * \endcode
 *
 * This header covers the routines for integrating a function.
 *
 * ### Description ###
 *
 * The routine \c LALSRombergIntegrate() performs the integral specified by the
 * structure \c input and the result is returned as \c result.  Any
 * additional parameters (other than the integration variable \f$x\f$) can be passed
 * as \c params.  The routine \c LALSRombergIntegrate() does not use
 * \c params but just passes it to the integrand.  The routine
 * \c LALDRombergIntegrate() is the same but for double precision.
 *
 * ### Operating Instructions ###
 *
 * The following program performs the integral \f$\int_0^2F(x)dx\f$ where
 * \f$F(x)=x^4\log(x+\sqrt{x^2+1})\f$.
 *
 * \code
 * #include <math.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/Integrate.h>
 *
 * static void F( LALStatus *s, REAL4 *y, REAL4 x, void *p )
 * {
 * REAL4 x2 = x*x;
 * REAL4 x4 = x2*x2;
 * INITSTATUS(s);
 * ASSERT( !p, s, 1, "Non-null pointer" );
 * y = x4 * log( x + sqrt( x2 + 1 ) );
 * RETURN( s );
 * }
 *
 * int main ()
 * {
 * const REAL4       epsilon = 1e-6;
 * const long double expect  = 8.153364119811650205L;
 * static LALStatus  status;
 * SIntegrateIn      intinp;
 * REAL4             result;
 *
 * intinp.function = F;
 * intinp.xmin     = 0;
 * intinp.xmax     = 2;
 * intinp.type     = ClosedInterval;
 *
 * LALSRombergIntegrate( &status, &result, &intinp, NULL );
 * if ( fabs( result - expect ) > epsilon * fabs( expect ) )
 * {
 * // integration did not achieve desired accuracy --- exit failure
 * return 1;
 * }
 *
 * return 0;
 * }
 * \endcode
 *
 * ### Algorithm ###
 *
 * This is an implementation of the Romberg integrating function \c qromb in
 * Numerical Recipes [\ref ptvf1992].
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define INTEGRATEH_ENULL 1		/**< Null pointer */
#define INTEGRATEH_ETYPE 2		/**< Unknown integral type */
#define INTEGRATEH_EIDOM 4		/**< Invalid domain */
#define INTEGRATEH_EMXIT 8		/**< Maximum iterations exceeded */
/*@}*/

/** \cond DONT_DOXYGEN */
#define INTEGRATEH_MSGENULL "Null pointer"
#define INTEGRATEH_MSGETYPE "Unknown integral type"
#define INTEGRATEH_MSGEIDOM "Invalid domain"
#define INTEGRATEH_MSGEMXIT "Maximum iterations exceeded"
/** \endcond */

/**
 * Type of integral.
 * The types of integration are the following:
 *
 * I. \c #ClosedInterval
 * indicates that the integral should be computed on equal-spaced domain
 * intervals including the boundary.
 *
 * II. \c #OpenInterval indicates that
 * the integral should be computed on intervals of the domain not including the
 * boundary.
 *
 * III. \c #SingularLowerLimit indicates that the integral
 * should be evaluated on an open interval with a transformation so that a
 * inverse-square-root singularity at the lower limit can be integrated.
 *
 * IV. \c #SingularUpperLimit is the same as above but for a singularity at the
 * upper limit.
 *
 * V. \c #InfiniteDomainPow indicates that the integral
 * should be evaluated over an semi-infinite domain---appropriate when both
 * limits have the same sign (though one is very large) and when the integrand
 * vanishes faster than \f$x^{-1}\f$ at infinity.
 *
 * VI. \c #InfiniteDomainExp
 * indicates that the integral should be evaluated over an infinite domain
 * starting at \c xmin and going to infinity (\c xmax is ignored)---the
 * integrand should vanish exponentially for large \f$x\f$.
 */
typedef enum
{
  ClosedInterval,     /**< evaluate integral on a closed interval             */
  OpenInterval,       /**< evaluate integral on an open interval              */
  SingularLowerLimit, /**< integrate an inv sqrt singularity at lower limit   */
  SingularUpperLimit, /**< integrate an inv sqrt singularity at upper limit   */
  InfiniteDomainPow,  /**< integrate infinite domain with power-law falloff   */
  InfiniteDomainExp   /**< integrate infinite domain with exponential falloff */
}
IntegralType;

/** These are input structures to the integration routines. */
typedef struct
tagSIntegrateIn
{
  void (*function)(LALStatus *s, REAL4 *y, REAL4 x, void *p);	/**<  The function to integrate */
  REAL4         xmax;	/**<  The maximum value of the domain of integration */
  REAL4         xmin;	/**< The minimum value of the domain of integration */
  IntegralType  type;	/**<  The type of integration */
}
SIntegrateIn;

/** These are input structures to the integration routines. */
typedef struct
tagDIntegrateIn
{
  void (*function)(LALStatus *s, REAL8 *y, REAL8 x, void *p);	/**<  The function to integrate */
  REAL8         xmax;	/**<  The maximum value of the domain of integration */
  REAL8         xmin;	/**< The minimum value of the domain of integration */
  IntegralType  type;	/**<  The type of integration */
}
DIntegrateIn;

/* ----- Integrate.c ----- */
/** \see See \ref Integrate_h for documentation */
void
LALSRombergIntegrate (
    LALStatus       *status,
    REAL4        *result,
    SIntegrateIn *input,
    void         *params
    );

/** \see See \ref Integrate_h for documentation */
void
LALDRombergIntegrate (
    LALStatus       *status,
    REAL8        *result,
    DIntegrateIn *input,
    void         *params
    );

/** \see See \ref Integrate_h for documentation */
REAL8
XLALREAL8RombergIntegrate (
    REAL8 (*f)(REAL8 x, void *params),
    void *params,
    REAL8 xmin,
    REAL8 xmax,
    IntegralType type
    );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRATE_H */
