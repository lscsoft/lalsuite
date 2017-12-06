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

#include <lal/LALNoiseModels.h>

/**
 * \author Sathyaprakash, B. S.
 * \ingroup LALNoiseModels_h
 * \brief Function to add two vectors with weights.
 *
 * ### Description ###
 *
 * Given weights \c A1 and \c A2 as in \c AddVectorsIn
 * and vectors \c v1 and \c v2 this code returns vector \c v
 * given by
 *
 * <tt>
 * v[i] = A1 v1[i] + A2 v2[i];
 * </tt>
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * none
 * \endcode
 *
 */
void
LALAddVectors
   (
   LALStatus   *status,
   REAL4Vector *vector,
   AddVectorsIn in
   )
{
   INT4 n, i;

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (vector,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (vector->length >= 2, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (vector->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (vector->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (in.v1->length == in.v2->length, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (in.v1->length == vector->length, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   n=vector->length;
   for (i=0; i<n; i++)
   {
      vector->data[i] = in.a1 * in.v1->data[i] + in.a2 * in.v2->data[i];
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
