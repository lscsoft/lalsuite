/*
*  Copyright (C) 2007 Thomas Cokelaer
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

/**
 * \ingroup LALInspiralBank_h
 * \brief Module to check whether a point with coordinates (x,y) is inside
 * a polygon defined by the vectors (vx, vy), which size (n) must be
 * provided. The functions returns 1 if the point is inside or 0 otherwise.
 *
 * \author Cokelaer. T
 *
 * \heading{Notes}
 * Tested in matlab codes and some BCV tests within lal/lalapps.
 *
 */
void LALInsidePolygon(  LALStatus          *status,	/**< LAL status pointer */
                        REAL4              *inputx,	/**< [in] two arrays of floats defining the polygon */
                        REAL4              *inputy,	/**< [in] two arrays of floats defining the polygon */
                        INT4               n,		/**< [in] the size of the vectors */
                        REAL4              x,		/**< [in] the coordinate of the point */
                        REAL4              y,		/**< [in] the coordinate of the point */
                        INT4               *valid	/**< [out] 0 if outside and 1 if inside */
                        )
{


   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (n>=3,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   {
     int i, j, c = 0;
     for (i = 0, j = n-1; i < n; j = i++) {
       if ((((inputy[i] <= y) && (y < inputy[j])) ||
	    ((inputy[j] <= y) && (y < inputy[i]))) &&
	   (x < (inputx[j] - inputx[i]) * (y - inputy[i]) / (inputy[j] - inputy[i]) + inputx[i]))
	 c = !c;
     }
     *valid = c;
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);

}
