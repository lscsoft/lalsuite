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

/* This information is taken from
   http://igec.lnl.infn.it/cgi-bin/browser.pl?Level=0,3,1
   Note that the IGEC convention for azimuth is clockwise from North 
   and ours is counterclockwise from East
*/

const LALFrDetector lalCachedBars[LALNumCachedBars]
= { { "AURIGA",
       11.0L + ( 56.0L + (54.0L/60.0L) ) / 60.0L,    
       45.0L + ( 21.0L + (12.0L/60.0L) ) / 60.0L,
      0,
      0, (90.0-44.0)*LAL_PI_180, 0, 0 } ,
    { "NAUTILUS",
       12.0L + ( 40.0L + (21.0L/60.0L) ) / 60.0L,
       41.0L + ( 49.0L + (26.0L/60.0L) ) / 60.0L,
      0,
      0, (90.0-44.0)*LAL_PI_180, 0, 0 } ,
    { "EXPLORER",
        6.0L + 12.0L / 60.0L,    
       46.0L + 27.0L / 60.0L,
      0,
      0, (90.0-39.0)*LAL_PI_180, 0, 0 } ,
    { "ALLEGRO (pre-2000)",
      268.0L + 50.0L / 60.0L - 360.0L,
       30.0L + 27.0L / 60.0L,    
      0,
      0, (90.0+40.0)*LAL_PI_180, 0, 0 } ,
    { "NIOBE",
    -( 31.0L + 56.0L / 60.0L),
      115.0L + 49.0L / 60.0L,
      0,
      0, (90.0+ 0.0)*LAL_PI_180, 0, 0 } ,
    /* Taken from Finn and Lazzarini gr-qc/0104040 */
    { "ALLEGRO (post-2000)",
    -( 91.0L + ( 10.0L + (43.766L/60.0L) ) / 60.0L),
       30.0L + ( 24.0L + (45.110L/60.0L) ) / 60.0L,
      0,
      0, (90.0+40.0)*LAL_PI_180, 0, 0 }
};
