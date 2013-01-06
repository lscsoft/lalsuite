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


#include <lal/LALInspiralBank.h>

#define Dim 3


/** \ingroup LALInspiralBank_h
 * \brief Uses \f$g^{ij} = 1/(2g) e^{ikl} e^{jab} g_{ka} g_{lb}\f$ to compute the inverse.
 * \author Sathyaprakash, B. S.
 *
 * \heading{Description}
 *
 * Uses \f$g^{ij} = 1/(2g) e^{ikl} e^{jab} g_{ka} g_{lb}\f$ to compute
 * the inverse; though not efficient, it is good enough for the
 * 3-d matrix that we have. Prevents the need for having a new library.
 *
 * \heading{Notes}
 * Do not generalise to more than 3 dimensions.
 */
void LALInverse3(LALStatus *status,	/**< LAL status pointer */
                 REAL8     **inverse,	/**< [out] inverse of the \f$(3\times 3)\f$ input matrix */
                 REAL8     **matrix	/**< [in] matrix whose inverse is required */
                 )
{

   REAL8 epsilon[Dim][Dim][Dim] = {{
                              { 0, 0, 0},
                              { 0, 0, 1},
                              { 0,-1, 0}},
                             {{ 0, 0,-1},
                              { 0, 0, 0},
                              { 1, 0, 0}},
                             {{ 0, 1, 0},
                              {-1, 0, 0},
                              { 0, 0, 0}}};
   REAL8 det,x;
   int i,j,p,q,a,b;

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (inverse,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (matrix,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   LALDeterminant3(status->statusPtr, &det, matrix); CHECKSTATUSPTR(status);

   ASSERT (det!=0,  status, LALINSPIRALBANKH_EDIV0, LALINSPIRALBANKH_MSGEDIV0);

   for (i=0; i<Dim; i++) {
   for (j=0; j<Dim; j++) {
      x = 0;
      for (a=0; a<Dim; a++) { for (b=0; b<Dim; b++) {
      for (p=0; p<Dim; p++) { for (q=0; q<Dim; q++) {
         x+=epsilon[j][a][p] * epsilon[i][b][q] * matrix[a][b] * matrix[p][q];
      }}}}
      inverse[i][j] = x/(2.*det);
   }}
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

#undef Dim
