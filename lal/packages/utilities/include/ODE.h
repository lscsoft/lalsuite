/*
*  Copyright (C) 2007 Jolien Creighton
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

#ifndef _ODE_H
#define _ODE_H

#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * \addtogroup ODE_h
 * \author J. D. E. Creighton
 *
 * \brief Routines for solving ordinary differential equations (ODEs).
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/ODE.h>
 * \endcode
 *
 * These routines solve ordinary differential equations (ODEs) of the form:
 * \f[
 * \dot{\mathbf{x}} = {\mathbf{f}}({\mathbf{x}},t,\ldots)
 * \f]
 * where \f$\mathbf{f}\f$ is a specified vector-valued function.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define ODEH_ENULL 001		/**< Null pointer. */
#define ODEH_ESAME 002		/**< Same data pointer. */
#define ODEH_ESIZE 004		/**< Invalid size. */
#define ODEH_ESZMM 010		/**< Size mismatch. */
#define ODEH_ENSTP 020		/**< Step number mismatch. */
/*@}*/

/** \cond DONT_DOXYGEN */
#define ODEH_MSGENULL "Null pointer."
#define ODEH_MSGESAME "Same data pointer."
#define ODEH_MSGESIZE "Invalid size."
#define ODEH_MSGESZMM "Size mismatch."
#define ODEH_MSGENSTP "Step number mismatch."
/** \endcond */

/** The independent variables of the ODE (parameters to the ODE function) */
typedef struct
tagREAL4ODEIndep
{
  REAL4  t;	/**< The independent parameter (e.g., time) that is evolved */
  void  *aux;	/**< Storage for auxiliary variables used internally in the ODE routine */
}
REAL4ODEIndep;

/** The parameters for the ODE step integrator */
typedef struct
tagREAL4ODEParams
{
  void ( *ode )( LALStatus *, REAL4Vector *, REAL4Vector *, REAL4ODEIndep * );	/**< Pointer to the function that computes the RHS of the ODE */
  REAL4ODEIndep       *indep;	/**< The independent variables of used by this function */
  REAL4                tstep;	/**< The suggested time step to use */
  REAL4Vector         *xdot;	/**< The value of the LHS of the ODE at the current time */
  REAL4Vector         *xerr;	/**< The estimated errors from the last ODE step */
  REAL4                eps;	/**< The allowed fractional error per step; if \c eps is \< \c LAL_REAL4_EPS, then the latter is used */
  REAL4VectorSequence *dx;	/**< Workspace storage for use in the the step integrator; the length of the sequence depends on which integrator is used */
}
REAL4ODEParams;

/** See ODE_h for documentation */
void LALSRungeKutta4(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    );

/** See ODE_h for documentation */
void LALSRungeKutta5(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    );

/** See ODE_h for documentation */
void LALSRungeKutta5Adapt(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _ODE_H */
