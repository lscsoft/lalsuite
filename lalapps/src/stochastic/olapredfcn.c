/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
*                2010 Nickolas Fotopoulos
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
 * File Name: olapredfcn.c
 *
 * Author: Brown, D. A.
 * 
 * 
 *-----------------------------------------------------------------------
 */

/**
 * \file
 * \ingroup lalapps_inspiral
 *
 *
 * <dl>
 *
 * <dt>Name</dt><dd>
 * <tt>lalapps_olapredfcn</tt> --- computes overlap reduction function given
 * a pair of known detectors.</dd>
 *
 * <dt>Synopsis</dt><dd>
 * <tt>lalapps_olapredfcn</tt>
 * [<tt>-h</tt>]
 * [<tt>-q</tt>]
 * [<tt>-v</tt>]
 * [<tt>-d</tt> <i>debugLevel</i>]
 * <tt>-s</tt> <i>siteID1</i>
 * [<tt>-a</tt> <i>azimuth1</i>]
 * <tt>-t</tt> <i>siteID2</i>
 * [<tt>-b</tt> <i>azimuth2</i>]
 * [<tt>-f</tt> <i>fLow</i>]
 * <tt>-e</tt> <i>deltaF</i>
 * <tt>-n</tt> <i>numPoints</i>
 * <tt>-o</tt> <i>outfile</i></dd>
 *
 * <dt>Description</dt><dd>
 * <tt>lalapps_olapredfcn</tt> computes the overlap reduction function
 * \f$\gamma(f)\f$ for a pair of known gravitational wave detectors. It uses
 * the LAL function <tt>LALOverlapReductionFunction()</tt>, which is
 * documented in the LAL Software Documentation under the
 * \c stochastic package.</dd>
 *
 * <dt>Options</dt><dd>
 * <dl>
 * <dt><tt>-h</tt></dt><dd>
 * Print a help message.</dd>
 *
 * <dt><tt>-q</tt></dt><dd>
 * Run silently (redirect standard input and error to <tt>/dev/null</tt>).</dd>
 *
 * <dt><tt>-v</tt></dt><dd>
 * Run in verbose mode.</dd>
 *
 * <dt><tt>-d</tt> <i>debugLevel</i></dt><dd>
 * Set the LAL debug level to <i>debugLevel</i>.</dd>
 *
 * <dt><tt>-s</tt> <i>siteID1</i>, <tt>-t</tt> <i>siteID2</i></dt><dd>
 * Use detector sites identified by <i>siteID1</i> and <i>siteID2</i>; ID
 * numbers between \c LALNumCachedDetectors (defined in the
 * \c tools package of LAL) refer to detectors cached in the constant
 * array <tt>lalCachedDetectors[]</tt>.  (At this point, these are all
 * interferometers.)  Additionally, the five resonant bar detectors of the
 * IGEC collaboration can be specified. The bar geometry data (summarized
 * in \ref table_cachedBars "this table" is used by the function
 * <tt>LALCreateDetector()</tt> from the \c tools package of LAL to
 * generate the Cartesian position vector and response tensor which are
 * used to calculate the overlap reduction function.  The ID numbers for
 * the bars depend on the value of \c LALNumCachedDetectors; the
 * correct ID numbers can be obtained by with the command
 * \code
 * > lalapps_olapredfcn -h
 * \endcode</dd>
 *
 * <dt><tt>-a</tt> <i>azimuth1</i>, <tt>-b</tt> <i>azimuth2</i></dt><dd>
 * If <i>siteID1</i> (<i>siteID2</i>) is a bar detector, assume it has an
 * azimuth of <i>azimuth1</i> (<i>azimuth2</i>) degrees East of North
 * rather than the default IGEC orientation given in
 * \ref table_cachedBars "this table".  Note that this convention, measuring
 * azimuth in degrees clockwise from North is not the same as that used in
 * LAL (which comes from the frame spec).  Note also that any specified
 * azimuth angle is ignored if the corresponding detector is an
 * interferometer.</dd>
 *
 * <dt><tt>-f</tt> <i>fLow</i></dt><dd>
 * Begin the frequency series at a frequency of <i>fLow</i>\,Hz; if this
 * is omitted, the default value of 0\,Hz is used.</dd>
 *
 * <dt><tt>-e</tt> <i>deltaF</i></dt><dd>
 * Construct the frequency series with a frequency spacing of
 * <i>deltaF</i>\,Hz</dd>
 *
 * <dt><tt>-n</tt> <i>numPoints</i></dt><dd>
 * Construct a frequency series with <i>numPoints</i> points.</dd>
 *
 * <dt><tt>-o</tt> <i>outfile</i></dt><dd>
 * Write the output to file \e outfile. The format of this file is
 * that output by the routine <tt>LALPrintFrequencySeries()</tt> in the
 * \c support package of LAL, which consists of a header describing
 * metadata followed by two-column rows, each containing the doublet
 * \f$\{f,\gamma(f)\}\f$.</dd>
 * </dl>
 *
 * \anchor table_cachedBars
 * <table class="doxtable" width="100%" align="center">
 * <caption align="top" style="text-align: left; font-weight: normal;">
 * Location and orientation data for the five IGEC resonant bar
 * detectors, stored in the <tt>lalCachedBars[]</tt> array. The data are
 * taken from
 * <tt>http://igec.lnl.infn.it/cgi-bin/browser.pl?Level=0,3,1</tt>
 * except for the latitude and longitude of ALLEGRO, which were taken from
 * Finn \& Lazzarini, gr-qc/0104040.  Note that the elevation above the
 * WGS-84 reference ellipsoid and altitude angle for each bar is not given,
 * and therefore set to zero.
 * </caption>
 * <tr><th>Name</th><th>Longitude</th><th>Latitude</th><th>Azimuth</th></tr>
 * <tr><td>AURIGA</td><td>\f$11^\circ56'54"\f$E</td><td>\f$45^\circ21'12"\f$N</td><td>N\f$44^\circ\f$E</td></tr>
 * <tr><td>NAUTILUS</td><td>\f$12^\circ40'21"\f$E</td><td>\f$41^\circ49'26"\f$N</td><td>N\f$44^\circ\f$E</td></tr>
 * <tr><td>EXPLORER</td><td>\f$6^\circ12'\f$E</td><td>\f$46^\circ27'\f$N</td><td>N\f$39^\circ\f$E</td></tr>
 * <tr><td>ALLEGRO</td><td>\f$91^\circ10'43.\!\!"766\f$W</td><td>\f$30^\circ24'45.\!\!"110\f$N</td><td>N\f$40^\circ\f$W</td></tr>
 * <tr><td>NIOBE</td><td>\f$115^\circ49'\f$E</td><td>\f$31^\circ56'\f$S</td><td>N\f$0^\circ\f$E</td></tr>
 * </table>
 *
 * <dt>Example</dt><dd>
 * To compute the overlap reduction function for LIGO Hanford and LIGO
 * Livingston, with a resolution of 1\,Hz from 0\,Hz to 1024\,Hz:
 *
 * \code
 * > lalapps_olapredfcn -s 0 -t 1 -e 1 -n 1025 -o LHOLLO.dat
 * \endcode
 *
 * To compute the overlap reduction function for ALLEGRO in its optimal
 * orientation of \f$72.\!\!^\circ08\f$ West of South (see Finn \& Lazzarini,
 * gr-qc/0104040) and LIGO Livingston, with a resolution of 0.5\,Hz
 * from 782.5\,Hz to 1032\,Hz (assuming \c lalNumCachedBars is 6):
 *
 * \code
 * > lalapps_olapredfcn -s 9 -a 252.08 -t 1 -f 782.5 -e 0.5 -n 500 -o ALLEGROLHO.dat
 * \endcode</dd>
 *
 * <dt>Author</dt><dd>
 * John T. Whelan
 * </dd>
 * </dl>
 */

#include "olapredfcn.h"

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose  = OLAPREDFCNH_FALSE;
REAL8 optDeltaF     = -1;
UINT4 optLength     = 0;
REAL8 optF0         = 0.0;
UINT4 optDetector1  = LALNumCachedDetectors;
UINT4 optDetector2  = LALNumCachedDetectors;
REAL4 optAzimuth1   = OLAPREDFCNH_OOR;
REAL4 optAzimuth2   = OLAPREDFCNH_OOR;
CHAR optFile[LALNameLength] = "";

int main( int argc, char *argv[] )
{
  LALStatus                status = blank_status;
  
  OverlapReductionFunctionParameters   parameters;
  REAL4FrequencySeries     overlap;

  LALDetectorPair     detectors;

  lal_errhandler = LAL_ERR_EXIT;

  olapredfcn_parse_options( argc, argv );

  if ( !optFile[0] ) {
    fprintf( stderr, "ERROR: No output file specified\n" );
    return OLAPREDFCNH_EARG;
  }

  if (optLength == 0) {
    fprintf( stderr, "ERROR: Zero length specified\n");
    return OLAPREDFCNH_EARG;
  }

  if (optDeltaF <= 0) {
    fprintf( stderr, "ERROR: Non-positive frequency spacing specified\n");
    return OLAPREDFCNH_EARG;
  }

  if (optF0 < 0) {
    fprintf( stderr, "ERROR: Negative start frequency specified\n");
    return OLAPREDFCNH_EARG;
  }

  if (optDetector1 >= LALNumCachedDetectors) {
    fprintf( stderr, "ERROR: Detector 1 unknown\n");
    return OLAPREDFCNH_EARG;
  }
  if (optDetector2 >= LALNumCachedDetectors) {
    fprintf( stderr, "ERROR: Detector 2 unknown\n");
    return OLAPREDFCNH_EARG;
  }

  if ( (( optAzimuth1 != OLAPREDFCNH_OOR ) &&
        (( optAzimuth1 < -360.0 ) || ( optAzimuth1 > 360.0 ))) ||
       (( optAzimuth2 != OLAPREDFCNH_OOR ) &&
        (( optAzimuth2 < -360.0 ) && ( optAzimuth2 > 360.0 )))) {
    fprintf( stderr, "ERROR: azimuth must be between -360 and 360\n" );
    return OLAPREDFCNH_EARG;
  }

  parameters.length = optLength;
  parameters.f0 = optF0;
  parameters.deltaF = optDeltaF;

  /* set detector from known detectors */
  detectors.detectorOne = lalCachedDetectors[optDetector1];
  detectors.detectorTwo = lalCachedDetectors[optDetector2];
  if ( optVerbose ) {
    fprintf( stderr, "Using cached detector %s for detector 1\n",
      detectors.detectorOne.frDetector.name );
    fprintf( stderr, "Using cached detector %s for detector 2\n",
      detectors.detectorTwo.frDetector.name );
  }

  /* bars are uniquely rotatable */
  if ( (( detectors.detectorOne.type != LALDETECTORTYPE_CYLBAR ) &&
        ( optAzimuth1 != OLAPREDFCNH_OOR )) ||
       (( detectors.detectorTwo.type != LALDETECTORTYPE_CYLBAR ) &&
        ( optAzimuth2 != OLAPREDFCNH_OOR )) ) {
    fprintf( stderr, "ERROR: Can only set azimuth for bar detectors\n" );
    return OLAPREDFCNH_EARG;
  }

  if ( detectors.detectorOne.type == LALDETECTORTYPE_CYLBAR ) {
    if ( optAzimuth1 != OLAPREDFCNH_OOR ) {
      if ( optVerbose ) {
        fprintf( stderr, "Changing azimuth to %.3f degrees East of North\n",
          optAzimuth1 );
      }
      detectors.detectorOne.frDetector.xArmAzimuthRadians =
        optAzimuth1 * LAL_PI_180;
    }
    else if ( optVerbose ) {
      fprintf( stderr, "Using IGEC azimuth of %.3f degrees East of North\n",
        detectors.detectorOne.frDetector.xArmAzimuthRadians * LAL_180_PI );
    }
  }
  if ( detectors.detectorTwo.type == LALDETECTORTYPE_CYLBAR ) {
    if ( optAzimuth2 != OLAPREDFCNH_OOR ) {
      if ( optVerbose ) {
        fprintf( stderr, "Changing azimuth to %.3f degrees East of North\n",
          optAzimuth2 );
      }
      detectors.detectorTwo.frDetector.xArmAzimuthRadians =
        optAzimuth2 * LAL_PI_180;
    }
    else if ( optVerbose ) {
      fprintf( stderr, "Using IGEC azimuth of %.3f degrees East of North\n",
        detectors.detectorTwo.frDetector.xArmAzimuthRadians * LAL_180_PI );
    }
  }

  overlap.data = NULL;
  LAL_CALL(
     LALSCreateVector(&status, &(overlap.data), optLength),
     &status );
  LAL_CALL(
     LALOverlapReductionFunction(&status, &overlap, &detectors, &parameters),
     &status );
  LAL_CALL(
     LALPrintFrequencySeries( &overlap, optFile ),
     &status );

  if ( optVerbose ) {
   fprintf(stderr,
     "======== Overlap Reduction Function Written to File %s ========\n",
     optFile);
  }

  LAL_CALL(
     LALSDestroyVector(&status, &(overlap.data)),
     &status );

  LALCheckMemoryLeaks();
  return OLAPREDFCNH_ENOM;
}
