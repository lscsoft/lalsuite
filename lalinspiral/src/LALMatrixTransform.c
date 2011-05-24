/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton
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
\author Sathyaprakash, B. S.
\file
\ingroup LALInspiralBank_h

\brief A routine to transform a second rank tensor under a given transformation.

\heading{Prototypes}

<tt>LALMatrixTransform()</tt>
<ul>
   <li> <tt>n,</tt> Input, dimension of the matrix (currently, and possibly always, only 3)
   </li><li> <tt>data1,</tt> Input, transformation matrix
   </li><li> <tt>data2,</tt> Input, matrix whose transformation is required
   </li><li> <tt>data3,</tt> Output, transformed matrix</li>
</ul>

\heading{Description}
Given the matrix of transformation in \c data1 and a second rank tensor
\c data2, this routine computes the transformed tensor in \c data3.

\heading{Algorithm}
\f[ C_{ij} = A_{im} A_{jl}  B_{ml}.\f]

\heading{Uses}
None.

\heading{Notes}

*/



#include <lal/LALInspiralBank.h>

NRCSID(LALMATRIXTRANSFORMC, "$Id$");



void LALMatrixTransform (LALStatus *status,
                         INT4 n,
                         REAL8 **data1,
                         REAL8 **data2,
                         REAL8 **data3)
{

   INT4 i, j, l, m;

   INITSTATUS(status, "LALMatrixTransform", LALMATRIXTRANSFORMC);
   ATTATCHSTATUSPTR(status);
   ASSERT (data1,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (data2,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (data3,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
      data3[i][j] = 0.0;
      for (l=0; l<n; l++) {
      for (m=0; m<n; m++) {
	 data3[i][j] += data1[i][m]*data2[m][l]*data1[j][l];
      }}
   }}
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
