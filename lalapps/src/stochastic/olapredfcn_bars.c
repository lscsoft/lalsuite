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
      ( 11.0L + ( 56.0L + (54.0L/60.0L) ) / 60.0L ) * LAL_PI_180,
      ( 45.0L + ( 21.0L + (12.0L/60.0L) ) / 60.0L ) * LAL_PI_180,
      0,
      0, 44.0*LAL_PI_180, 0, 0 } ,
    { "NAUTILUS",
      ( 12.0L + ( 40.0L + (21.0L/60.0L) ) / 60.0L ) * LAL_PI_180,
      ( 41.0L + ( 49.0L + (26.0L/60.0L) ) / 60.0L ) * LAL_PI_180,
      0,
      0, 44.0*LAL_PI_180, 0, 0 } ,
    { "EXPLORER",
      (  6.0L + 12.0L / 60.0L ) * LAL_PI_180,
      ( 46.0L + 27.0L / 60.0L ) * LAL_PI_180,
      0,
      0, 39.0*LAL_PI_180, 0, 0 } ,
    /* Taken from Finn and Lazzarini gr-qc/0104040 */
    { "ALLEGRO ",
     -( 91.0L + ( 10.0L + (43.766L/60.0L) ) / 60.0L ) * LAL_PI_180,
      ( 30.0L + ( 24.0L + (45.110L/60.0L) ) / 60.0L ) * LAL_PI_180,
      0,
      0, -40.0*LAL_PI_180, 0, 0 } ,
    { "NIOBE",
      (115.0L + 49.0L / 60.0L ) * LAL_PI_180,
     -( 31.0L + 56.0L / 60.0L ) * LAL_PI_180,
      0,
      0, 0.0*LAL_PI_180, 0, 0 }
};
