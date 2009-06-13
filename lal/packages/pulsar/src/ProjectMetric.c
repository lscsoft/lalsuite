/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix
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

/**
 * \author Creighton, T. D.
 * \date 2000
 * \file
 * \ingroup PulsarMetric
 * \brief Project out the zeroth dimension of a metric.
 *
 * $Id$
 *

 \par Description

 This function takes a metric \f$g_{\alpha\beta}\f$, where
 \f$\alpha,\beta=0,1,\ldots,n\f$, and computes the projected metric
 \f$\gamma_{ij}\f$ on the subspace \f$i,j=1,\ldots,n\f$, as described in the
 header StackMetric.h.

The argument \a *metric stores the metric components in the manner
used by the functions CoherentMetric() and
StackMetric(), and \a errors indicates whether error
estimates are included in \a *metric.  Thus \a *metric is a
vector of length \f$(n+1)(n+2)/2\f$ if \a errors is zero, or of length
\f$(n+1)(n+2)\f$ if \a errors is nonzero; see CoherentMetric.c
for the indexing scheme.

Upon return, \a *metric stores the components of \f$\gamma_{ij}\f$ in
the same manner as above, with the physically meaningless components
\f$\gamma_{\alpha0} = \gamma_{0\alpha}\f$ (and their uncertainties) set
identically to zero.

\par Algorithm

The function simply implements \ref eq_gij_gab "Eq. (2)" in header
StackMetric.h.  The formula used to convert uncertainties
\f$s_{\alpha\beta}\f$ in the metric components \f$g_{\alpha\beta}\f$ into
uncertainties \f$\sigma_{ij}\f$ in \f$\gamma_{ij}\f$ is:
\f[
\sigma_{ij} = s_{ij}
	+ s_{0i}\left|\frac{g_{0j}}{g_{00}}\right|
	+ s_{0j}\left|\frac{g_{0i}}{g_{00}}\right|
	+ s_{00}\left|\frac{g_{0i}g_{0j}}{(g_{00})^2}\right| \; .
\f]
Note that if the metric is highly degenerate, one may find that one or
more projected components are comparable in magnitude to their
estimated numerical uncertainties.  This can occur when the
observation times are very short or very long compared to the
timescales over which the timing derivatives are varying.  In the
former case, one is advised to use analytic approximations or a
different parameter basis.  In the latter case, the degenerate
components are often not relevant for data analysis, and can be
effectively set to zero.

Technically, starting from a full metric
\f$g_{\alpha\beta}(\mathbf{\lambda})\f$, the projection
\f$\gamma_{ij}(\vec\lambda)\f$ is the metric of a subspace
\f$\{\vec\lambda\}\f$ passing through the point \f$\mathbf{\lambda}\f$ on a plane
orthogonal to the \f$\lambda^0\f$ axis.  In order for \f$\gamma_{ij}\f$ to
measure the \em maximum distance between points \f$\vec\lambda\f$, it
is important to evaluate \f$g_{\alpha\beta}\f$ at the value of \f$\lambda^0\f$
that gives the largest possible separations.  For the pulsar search
formalism discussed in the header StackMetric.h, this is always
achieved by choosing the largest value of \f$\lambda^0=f_\mathrm{max}\f$
that is to be covered in the search.

*/

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/StackMetric.h>

NRCSID(PROJECTMETRICC,"$Id$");

/* <lalVerbatim file="ProjectMetricCP"> */
void
LALProjectMetric( LALStatus *stat, REAL8Vector *metric, BOOLEAN errors )
{ /* </lalVerbatim> */
  UINT4 s;     /* The number of parameters before projection. */
  UINT4 i, j;  /* Indecies. */
  REAL8 *data; /* The metric data array. */

  INITSTATUS(stat,"ProjectMetric",PROJECTMETRICC);

  /* Check that data exist. */
  ASSERT(metric,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(metric->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  data=metric->data;

  /* Make sure the metric's length is compatible with some
     dimensionality s. */
  for(s=0,i=1;i<metric->length;s++){
    i=s*(s+1);
    if(!errors)
      i=i>>1;
  }
  s--; /* counteract the final s++ */
  ASSERT(i==metric->length,stat,STACKMETRICH_EBAD,
	 STACKMETRICH_MSGEBAD);

  /* Now project out the zeroth component of the metric. */
  for(i=1;i<s;i++)
    for(j=1;j<=i;j++){
      INT4 i0 = (i*(i+1))>>1;
      INT4 myj0 = (j*(j+1))>>1;
      INT4 ij = i0+j;
      if(errors){
	data[2*ij]-=data[2*i0]*data[2*myj0]/data[0];
	data[2*ij+1]+=data[2*i0+1]*fabs(data[2*myj0]/data[0])
	  + data[2*myj0+1]*fabs(data[2*i0]/data[0])
	  + data[1]*fabs(data[2*i0]/data[0])*fabs(data[2*myj0]/data[0]);
      } else
	data[ij]-=data[i0]*data[myj0]/data[0];
    }

  /* Set all irrelevant metric coefficients to zero. */
  for(i=0;i<s;i++){
    INT4 i0 = i*(i+1)>>1;
    if(errors)
      data[2*i0]=data[2*i0+1]=0.0;
    else
      data[i0]=0.0;
  }

  /* And that's it! */
  RETURN(stat);
}
