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

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <math.h>
#include <lal/Dirichlet.h>


/**
 * \brief Calculates the values of the Dirichlet kernel [\ref pw93].
 * \f[
 * {\cal D}_N(x):=
 * \left\{
 * \begin{array}{cl}
 * (-1)^{x(N-1)} & \quad x=0, \pm 1, \pm 2, \cdots\\
 * \frac{\sin(N\pi x)}{N\sin(\pi x)} & \quad \textrm{otherwise}\ .
 * \end{array}
 * \right.
 * \f]
 *
 * The magnitude of the Dirichlet kernel is \f$1/N\f$ times the magnitude
 * of the discrete Fourier transform of the discrete \f$N\f$-point rectangular
 * window.
 *
 * ### Notes ###
 *
 * The Dirichlet kernel is needed for a rigorous (i.e., exact)
 * calculation of the standard cross-correlation statistic, evaluated in
 * discrete time.  However, <tt>LALStochasticOptimalFilter()</tt> and
 * other routines in the \c stochastic package, as currently
 * implemented, do not make use of the Dirichlet kernel; these routines
 * are expected to be used in the large observation time continuum limit
 * approximation, for which \f${\cal D}_N(x)\f$ can effectively be replaced
 * by a Dirac delta function.
 *
 */
void
LALDirichlet(LALStatus*                 status,
             REAL4Vector*               output,
             const DirichletParameters* parameters )

{
  UINT4  i;
  UINT4  n;         /* Dirichlet parameter N */
  UINT4  length;   /* specified length of output vector */
  REAL8 deltaX;    /* spacing of x values */
  REAL8 x;
  REAL8 top;
  REAL8 bot;

  /* initialize status structure */
  INITSTATUS(status);

  /* check that pointer to input parameters is not null */
  ASSERT(parameters != NULL,status,
         DIRICHLETH_ENULLPIN,DIRICHLETH_MSGENULLPIN);

  /* check that Dirichlet parameter N is > 0 */
  ASSERT(parameters->n > 0,status,
         DIRICHLETH_ENVALUE,DIRICHLETH_MSGENVALUE);

  /* assign valid data to Dirichlet parameter N */
  n=parameters->n;

  /* check that length parameter is greater than zero */
  ASSERT(parameters->length > 0,status,
         DIRICHLETH_ESIZE,DIRICHLETH_MSGESIZE);

  /* assign valid data to specified length */
  length=parameters->length;

  /* check that spacing of x values is greater than zero */
  ASSERT(parameters->deltaX > 0,status,DIRICHLETH_EDELTAX,
         DIRICHLETH_MSGEDELTAX);

  /* assign valid data to delta x */
  deltaX=parameters->deltaX;

  /* check that pointer to output vector is not null */
  ASSERT(output != NULL,status,
         DIRICHLETH_ENULLPOUT,DIRICHLETH_MSGENULLPOUT);

  /* check length of output vector agrees with length specified in
     input parameters */
  ASSERT(output->length==length,status,DIRICHLETH_ESIZEMM,
         DIRICHLETH_MSGESIZEMM);

  /* check that pointer to data member of output vector is not null */
  ASSERT(output->data != NULL,status,
         DIRICHLETH_ENULLPDOUT,DIRICHLETH_MSGENULLPDOUT);

  /* everything okay here ----------------------------------------------*/

  /* calculate the values of the Dirichlet kernel D_N(x) */

  output->data[0] = 1;  /* D_N(0)=1 */

  for ( i = 1; i < length; i++) {
    x = i * deltaX;

    if ( x-floor(x) == 0 ) {  /* D_N(x)=(-1)^(x * (N-1)) */
      output->data[i] = pow(-1.0, ( (INT4)(x * (n-1)) ) );
    }
    else {
      top = sin(LAL_PI * n * deltaX * i);
      bot = n * sin(LAL_PI * deltaX * i);
      output->data[i] = top/bot;
    }
  }

  return;
}
