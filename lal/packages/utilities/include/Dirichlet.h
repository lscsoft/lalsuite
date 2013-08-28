/*
*  Copyright (C) 2007 John Whelan
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

#ifndef  _DIRICHLET_H
#define  _DIRICHLET_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup Dirichlet_h
 * \author UTB Relativity Group; contact whelan@phys.utb.edu
 *
 * \brief Provides prototype and error code information for <tt>LALDirichlet()</tt>,
 * a routine which calculates the values of the Dirichlet kernel
 * \f${\cal D}_N(x)\f$.
 *
 * \heading{Synopsis}
 * \code
 * #include "Dirichlet.h"
 * \endcode
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define DIRICHLETH_ENULLPIN    1		/**< Null pointer to input parameters */
#define DIRICHLETH_ENVALUE     2		/**< Dirichlet parameter N less than or equal to zero */
#define DIRICHLETH_ESIZE       3		/**< Length parameter less than or equal to zero */
#define DIRICHLETH_EDELTAX     4		/**< Spacing of x values less than or equal to zero */
#define DIRICHLETH_ENULLPOUT   5		/**< Null pointer to ouput vector */
#define DIRICHLETH_ESIZEMM     6		/**< Length of data member of output vector does not equal length specified in input parameters */
#define DIRICHLETH_ENULLPDOUT  7		/**< Null pointer to data member of output vector */
/*@}*/

/** \cond DONT_DOXYGEN */
#define DIRICHLETH_MSGENULLPIN   "Null pointer to input parameters"
#define DIRICHLETH_MSGENVALUE    "Dirichlet parameter N less than or equal to zero"
#define DIRICHLETH_MSGESIZE      "Length parameter less than or equal to zero"
#define DIRICHLETH_MSGEDELTAX    "Spacing of x values less than or equal to zero"
#define DIRICHLETH_MSGENULLPOUT  "Null pointer to ouput vector"
#define DIRICHLETH_MSGESIZEMM    "Length of data member of output vector does not equal length specified in input parameters"
#define DIRICHLETH_MSGENULLPDOUT "Null pointer to data member of output vector"
/** \endcond */


/** Contains parameters that specify the Dirichlet kernel \f$\mathcal{D}_N(x)\f$ */
typedef struct tagDirichletParameters{
  UINT4	 n;       /**< Dirichlet parameter \f$N\f$ */
  UINT4	 length;  /**< Specified length of output vector */
  REAL8	 deltaX;  /**< Spacing of \f$x\f$ values */
} DirichletParameters;

void
LALDirichlet(LALStatus*,
	     REAL4Vector*,
	     const DirichletParameters*);

/*@}*/

#ifdef  __cplusplus
}
#endif /* C++ protection */
#endif /* _DIRICHLET_H */
