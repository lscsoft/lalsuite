/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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

/*----------------------------------------------------------------------- 
 * 
 * File Name: olapredfcn_bars.c
 *
 * Author: John T. Whelan
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "olapredfcn.h"

/* This information (except for the ALLEGRO latitude and longitude) is
   taken from http://igec.lnl.infn.it/cgi-bin/browser.pl?Level=0,3,1
   Note that our convention for azimuth is now clockwise from north,
   which agrees with the IGEC convention
*/

const LALFrDetector lalCachedBars[LALNumCachedBars]
= { { "AURIGA",
      ( 11.0 + ( 56.0 + (54.0/60.0) ) / 60.0 ) * LAL_PI_180,
      ( 45.0 + ( 21.0 + (12.0/60.0) ) / 60.0 ) * LAL_PI_180,
      0,
      0, 44.0*LAL_PI_180, 0, 0 } ,
    { "NAUTILUS",
      ( 12.0 + ( 40.0 + (21.0/60.0) ) / 60.0 ) * LAL_PI_180,
      ( 41.0 + ( 49.0 + (26.0/60.0) ) / 60.0 ) * LAL_PI_180,
      0,
      0, 44.0*LAL_PI_180, 0, 0 } ,
    { "EXPLORER",
      (  6.0 + 12.0 / 60.0 ) * LAL_PI_180,
      ( 46.0 + 27.0 / 60.0 ) * LAL_PI_180,
      0,
      0, 39.0*LAL_PI_180, 0, 0 } ,
    /* Taken from Finn and Lazzarini gr-qc/0104040 */
    { "ALLEGRO ",
     -( 91.0 + ( 10.0 + (43.766/60.0) ) / 60.0 ) * LAL_PI_180,
      ( 30.0 + ( 24.0 + (45.110/60.0) ) / 60.0 ) * LAL_PI_180,
      0,
      0, -40.0*LAL_PI_180, 0, 0 } ,
    { "NIOBE",
      (115.0 + 49.0 / 60.0 ) * LAL_PI_180,
     -( 31.0 + 56.0 / 60.0 ) * LAL_PI_180,
      0,
      0, 0.0*LAL_PI_180, 0, 0 }
};
