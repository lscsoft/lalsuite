/*-----------------------------------------------------------------------
 *
 * File Name:  SiteInfo.h
 *
 * Author: J.D. Romano
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * SiteInfo.h
 *
 * SYNOPSIS
 * #include <lal/SiteInfo.h>
 *
 * DESCRIPTION
 * Site location and orientation information for the 5 major interferometers
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 *
 * REVISION HISTORY
 *
 * $Log$
 * Revision 1.4  2000/11/09 23:37:20  jolien
 * Support for LAM; Cleaned up LALConfig.h; Put headers in lal subdir.
 *
 * Revision 1.3  2000/09/08 23:19:25  jolien
 * Changed LALDebugLevel to lalDebugLevel.
 *
 * Revision 1.2  2000/07/23 01:14:48  jolien
 * LAL-prefix in front of functions and Status; debuglevel is now lalDebugLevel.
 *
 * Revision 1.1  2000/06/19 22:33:35  jolien
 * Added LALCrossCorr and SiteInfo.
 *
 *
 *-----------------------------------------------------------------------
 */

#ifndef _SITEINFO_H
#define _SITEINFO_H

typedef enum { LHO, LLO, VIRGO, GEO600, TAMA300 } Detector;

/*
 *
 * LIGO Hanford Observatory (LHO)
 * Data taken from:
 *
 * William Althouse, Larry Jones, Albert Lazzarini (1999)
 * "Determination of Global and Local Coordinate Axes for the LIGO Sites"
 * LIGO-T980044-08-E
 *
 */

/* latitude N */
#define LHO_LAT_DEG  46
#define LHO_LAT_MIN  27
#define LHO_LAT_SEC  18.528

/* longitude E */
#define LHO_LON_DEG 240
#define LHO_LON_MIN  35
#define LHO_LON_SEC  32.4343

/* elevation (meters) relative to the WGS84 ellipsoid */
#define LHO_ELEV_SI 142.554

/* arm azimuth CCW from E (geodetic N = 90 degrees) */
#define LHO_ARM_X_AZI_RAD (125.9994*M_PI/180.0)
#define LHO_ARM_Y_AZI_RAD (215.9994*M_PI/180.0)

/* arm length (meters) */
#define LHO_ARM_LEN_SI 4000.0


/*
 *
 * LIGO Livingston Observatory (LLO)
 * Data taken from:
 *
 * William Althouse, Larry Jones, Albert Lazzarini (1999)
 * "Determination of Global and Local Coordinate Axes for the LIGO Sites"
 * LIGO-T980044-08-E
 *
 */

/* latitude N */
#define LLO_LAT_DEG  30
#define LLO_LAT_MIN  33
#define LLO_LAT_SEC  46.4196

/* longitude E */
#define LLO_LON_DEG 269
#define LLO_LON_MIN  13
#define LLO_LON_SEC  32.7346

/* elevation (meters) relative to the WGS84 ellipsoid */
#define LLO_ELEV_SI (-6.574)

/* arm azimuth CCW from E (geodetic N = 90 degrees) */
#define LLO_ARM_X_AZI_RAD (197.7165*M_PI/180)
#define LLO_ARM_Y_AZI_RAD (287.7165*M_PI/180)

/* arm length (meters) */
#define LLO_ARM_LEN_SI 4000.0


/*
 *
 * VIRGO Interferometer (VIRGO)
 * Data provided by: 
 *
 * Benoit Mours
 * E-mail:  mours@lapp.in2p3.fr
 *
 */

/* latitude N */
#define VIRGO_LAT_DEG  43
#define VIRGO_LAT_MIN  37
#define VIRGO_LAT_SEC  53.0921

/* longitude E */
#define VIRGO_LON_DEG  10
#define VIRGO_LON_MIN  30
#define VIRGO_LON_SEC  16.1878

/* elevation (meters) relative to the WGS84 ellipsoid */
#define VIRGO_ELEV_SI  51.884

/* arm azimuth CCW from E (geodetic N = 90 degrees) */
#define VIRGO_ARM_X_AZI_RAD  (70.5674*M_PI/180.0)
#define VIRGO_ARM_Y_AZI_RAD (160.5674*M_PI/180.0)

/* arm length (meters) */
#define VIRGO_ARM_LEN_SI 3000.0


/*
 *
 * GEO-600 Interferometer (GEO600)
 * Data taken from:
 *
 * http://www.geo600.uni-hannover.de/geo600/project/location.html
 *
 */

/* latitude N */
#define GEO600_LAT_DEG  52
#define GEO600_LAT_MIN  14
#define GEO600_LAT_SEC  42.528

/* longitude E */
#define GEO600_LON_DEG   9
#define GEO600_LON_MIN  48
#define GEO600_LON_SEC  25.894

/* elevation (meters) relative to the WGS84 ellipsoid */
#define GEO600_ELEV_SI (114.425)

/* arm azimuth CCW from E (geodetic N = 90 degrees) */
#define GEO600_ARM_X_AZI_RAD (115.9431*M_PI/180)
#define GEO600_ARM_Y_AZI_RAD  (21.6117*M_PI/180)

/* arm length (meters) */
#define GEO600_ARM_LEN_SI 600.0


/*
 *
 * TAMA Interferometer (TAMA300)
 * Data taken from: August, 1999
 *
 * Masa-Katsu Fujimoto (1995)
 * unpublished
 * E-mail:  fujimoto@gravity.mtk.nao.ac.jp 
 *
 */

/* latitude N */
#define TAMA300_LAT_DEG  35
#define TAMA300_LAT_MIN  40
#define TAMA300_LAT_SEC  35.6

/* longitude E */
#define TAMA300_LON_DEG 139
#define TAMA300_LON_MIN  32
#define TAMA300_LON_SEC   9.8

/* elevation (meters) relative to the WGS84 ellipsoid */
#define TAMA300_ELEV_SI  90

/* arm azimuth CCW from E (geodetic N = 90 degrees) */
#define TAMA300_ARM_X_AZI_RAD (180*M_PI/180.0)
#define TAMA300_ARM_Y_AZI_RAD (270*M_PI/180.0)

/* arm length (meters) */
#define TAMA300_ARM_LEN_SI 300.0


#endif /* _SITEINFO_H */

