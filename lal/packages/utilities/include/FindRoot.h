/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
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

#ifndef _FINDROOT_H
#define _FINDROOT_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * \addtogroup FindRoot_h
 *
 * \brief Root finding routines.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/FindRoot.h>
 * \endcode
 *
 * This header covers the routines for root finding.
 *
 * ### Description ###
 *
 * The routine <tt>LALSBracketRoot()</tt> expands the specified domain until a root
 * is contained.  The routine <tt>LALDBracketRoot()</tt> is the same but for a
 * double-precision function.
 *
 * The routine <tt>LALSBisectionFindRoot()</tt> bisects the domain (which must contain one
 * root) until the root is found with the desired accuracy.  The routine
 * <tt>LALDBisectionFindRoot()</tt> is the same but for a double-precision function.
 *
 * ### Operating Instructions ###
 *
 * Suppose we want to find the root of the function \f$y = F(x;y_0) = y_0 + x^2\f$.
 * Define the function:
 * \code
 * static void F( LALStatus *status, REAL4 *y, REAL4 x, void *y0 )
 * {
 * INITSTATUS(status);
 * ASSERT( y0, status, 1, "Null pointer" );
 * y = *(REAL4 *)y0 + x*x;
 * RETURN( status );
 * }
 * \endcode
 *
 * Then use the following code to bracket and find the root \f$x_0=1\f$ where
 * \f$F(x_0;y_0=-1)=0\f$:
 * \code
 * static LALStatus status;
 * SFindRootIn      input;
 * REAL4            y0;
 * REAL4            x0;
 *
 * y0             = -1;
 * input.function = F;
 * input.xmin     = 0.1;
 * input.xmax     = 0.2;
 * input.xacc     = 1e-5;
 *
 * /\* expand domain until a root is bracketed *\/
 * LALSBracketRoot( &status, &input, &y0 );
 *
 * /\* bisect domain until root is found *\/
 * LALSBisectionFindRoot( &status, &x0, &input, &y0 );
 * \endcode
 *
 * ### Algorithm ###
 *
 * This is an implementation of the root bracketing and bisection finding
 * routines \c zbrac and \c rtbis in Numerical Recipes \cite ptvf1992.
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define FINDROOTH_ENULL 1		/**< Null pointer */
#define FINDROOTH_EIDOM 2		/**< Invalid initial domain */
#define FINDROOTH_EMXIT 4		/**< Maximum iterations exceeded */
#define FINDROOTH_EBRKT 8		/**< Root not bracketed */
/*@}*/

/** \cond DONT_DOXYGEN */
#define FINDROOTH_MSGENULL "Null pointer"
#define FINDROOTH_MSGEIDOM "Invalid initial domain"
#define FINDROOTH_MSGEMXIT "Maximum iterations exceeded"
#define FINDROOTH_MSGEBRKT "Root not bracketed"
/** \endcond */

/**
 * These are function pointers to functions that map REAL4 numbers to REAL4 numbers.
 */
typedef struct
tagSFindRootIn
{
  void (*function)(LALStatus *s, REAL4 *y, REAL4 x, void *p);	/**< The function to find the root of */
  REAL4 xmax;	/**< The maximum value of the domain interval to look for the root */
  REAL4 xmin;	/**< The minimum value of the domain interval to look for the root */
  REAL4 xacc;	/**< The accuracy desired for the root */
}
SFindRootIn;

/**
 * These are function pointers to functions that map REAL8 numbers to REAL8 numbers.
 */
typedef struct
tagDFindRootIn
{
  void (*function)(LALStatus *s, REAL8 *y, REAL8 x, void *p);	/**< The function to find the root of */
  REAL8 xmax;	/**< The maximum value of the domain interval to look for the root */
  REAL8 xmin;	/**< The minimum value of the domain interval to look for the root */
  REAL8 xacc;	/**< The accuracy desired for the root */
}
DFindRootIn;

/* ----- FindRoot.c ----- */
void
LALSBracketRoot (
    LALStatus      *status,
    SFindRootIn *inout,
    void        *params
    );

int
XLALDBracketRoot(
    REAL8 (*y)(REAL8, void *),
    REAL8 *xmin,
    REAL8 *xmax,
    void *params
);

void
LALDBracketRoot (
    LALStatus      *status,
    DFindRootIn *inout,
    void        *params
    );

void
LALSBisectionFindRoot (
    LALStatus      *status,
    REAL4       *root,
    SFindRootIn *input,
    void        *params
    );

REAL8
XLALDBisectionFindRoot (
    REAL8 (*y)(REAL8, void *),
    REAL8 xmin,
    REAL8 xmax,
    REAL8 xacc,
    void *params
);

void
LALDBisectionFindRoot (
    LALStatus      *status,
    REAL8       *root,
    DFindRootIn *input,
    void        *params
    );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _FINDROOT_H */
