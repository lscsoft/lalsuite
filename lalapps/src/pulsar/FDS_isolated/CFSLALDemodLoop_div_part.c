/*
*  Copyright (C) 2007 Bernd Machenschalk
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

          {

            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
	    REAL4 accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */
	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
	    REAL4 XRes=0.0, XIms=0.0; /* sums of Xa.re and Xa.im */

            for(k=0; k < klim ; k++)
	      {
                XRes = tempFreq1 * XRes + (*Xalpha_k).re * accFreq;
                XIms = tempFreq1 * XIms + (*Xalpha_k).im * accFreq;

		accFreq *= (REAL4)tempFreq1;
                tempFreq1 --;
                Xalpha_k ++;
              } /* for k < klim */

	    accFreq = 1.0 / accFreq;

	    XRes *= accFreq;
	    XIms *= accFreq;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
