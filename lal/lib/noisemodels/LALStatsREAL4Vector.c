/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash
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
 * \brief Module to compute the mean, rms, minimum and maximum of a \c REAL4Vector.
 *
 */
void
LALStatsREAL4Vector
   (
   LALStatus           *status,
   StatsREAL4VectorOut *out,
   REAL4Vector         *vector
   )

{

   INT4 i, n;
   REAL8 x;

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (vector->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (vector->length > 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
   ASSERT (out,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   out->mean = 0.;
   out->var = 0.;
   n = vector->length;
   out->max = vector->data[0];
   out->min = vector->data[0];

   for (i=0; i<n; i++)
   {
      x = vector->data[i];
      if (out->max < x) out->max = x;
      if (out->min > x) out->min = x;
      out->mean+=x;
   }
   out->mean/=(REAL8) n;

   for (i=0; i<n; i++)
   {
      x = vector->data[i]-out->mean;
      out->var+=x*x;
   }
   out->var /=(REAL8) n;
   out->stddev = sqrt(out->var);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
