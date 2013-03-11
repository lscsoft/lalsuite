/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Matt Pitkin
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

/**
 * \author Matt Pitkin, Bernd Machenschalk
 * \date 2007
 * \file
 * \ingroup pulsarTODO
 * \brief Functions to calculate binary system time delays and read TEMPO pulsar parameter files
 *
 *
 * The main function in this code - <tt>XLALBinaryPulsarDeltaT</tt> - is for
 *  calculating the time delay on a pulsar signal caused by its orbit within a
 * binary system. It tranforms a signal from the binary system barycentre to the
 * pulsar proper time. It is equivalent to the <tt>LALBarycenter()</tt>
 * functions for transforming the Earth reference frame to the solar system
 * barycentre. It copies the functions from the standard pulsar timing software
 * <tt>TEMPO</tt> and can use the various binary system models BT, BT1P, BT2P,
 * BTX, ELL1, DD and MSS.
 *
 * Included is a function to read in a set of pulsar parameters from a
 * standard <tt>TEMPO</tt> parameter (.par) file -
 * <tt>XLALReadTEMPOParFile</tt>.
 *
 * Also included are function to convert times given in the TT, TDB or TCB
 * frames (given as a modified Julian Date - MJD) into a GPS time.
 */

#ifndef _BINARYPULSARTIMING_H
#define _BINARYPULSARTIMING_H

#include <ctype.h>
#include <unistd.h>

#include <lal/LALStdlib.h>
#include <lal/StringVector.h>
#include <lal/LALBarycenter.h>
#include <lal/Date.h>

#ifdef __cplusplus
extern "C" {
#endif

/**\name Error Codes */ /*@{*/
#define BINARYPULSARTIMINGH_ENULLINPUT 1
#define BINARYPULSARTIMINGH_ENULLOUTPUT 2
#define BINARYPULSARTIMINGH_ENULLPARAMS 3
#define BINARYPULSARTIMINGH_ENULLBINARYMODEL 4
#define BINARYPULSARTIMINGH_ENULLPULSARANDPATH 5
#define BINARYPULSARTIMINGH_EPARFILEERROR 6
#define BINARYPULSARTIMINGH_EFAIL 7
#define BINARYPULSARTIMINGH_ENAN 8
#define BINARYPULSARTIMINGH_EPARFAIL 9

#define BINARYPULSARTIMINGH_MSGENULLINPUT "Input was Null"
#define BINARYPULSARTIMINGH_MSGENULLOUTPUT "Output was Null"
#define BINARYPULSARTIMINGH_MSGENULLPARAMS "Params was Null"
#define BINARYPULSARTIMINGH_MSGNULLBINARYMODEL "Binary model is Null or not specified - you should\
not be in the binary timing routine"
#define BINARYPULSARTIMINGH_MSGENULLPULSARANDPATH "Path to pulsar.par file not specified"
#define BINARYPULSARTIMINGH_MSGEPARFILEERROR ".par file path is wrong or file doesn't exist"
#define BINARYPULSARTIMINGH_MSGEFAIL "Time delay computation failed"
#define BINARYPULSARTIMINGH_MSGENAN "Output is NaN!"
#define BINARYPULSARTIMINGH_MSGEPARFAIL "Failed to read in .par file"

#define GPS0MJD 44244.0 /* start of GPS time (MJD 44244) */
/* the difference between TDT (or TT) and TAI (see e.g. Eqn 15 of T.-Y. Huang et
 * al, A&A, 220, p. 329, 1989) */
#define TDT_TAI 32.184
/* the difference between TDT/TT and the GPS epoch */
#define GPS_TDT (TDT_TAI + XLAL_EPOCH_GPS_TAI_UTC)
                        
/*@}*/

/** A structure to contain all pulsar parameters and associated errors.
    The structure does not have to be used for a binary pulsar, but can just contain the parameters for
   an isolated pulsar. All parameters are in the same units as given by TEMPO.
 */
typedef struct
tagBinaryPulsarParams
{
  CHAR *name;   /**< pulsar name */
  CHAR *jname;  /**< pulsar J name */
  CHAR *bname;  /**< pulsar B name */
  
  CHAR *model;  /**< TEMPO binary model e.g. BT, DD, ELL1 */

  REAL8 f0;     /**< spin frequency (Hz) */
  REAL8 f1;     /**< frequency first derivative (Hz/s) */
  REAL8 f2;     /**< frequency second derivative (Hz/s^2) */
  REAL8 f3;     /**< frequency third derivative (Hz/s^3) */
  REAL8 f4;     /**< frequency fourth derivative (Hz/s^4) */
  REAL8 f5;     /**< frequency fifth derivative (Hz/s^5) */
  REAL8 f6;     /**< frequency sixth derivative (Hz/s^6) */
  REAL8 f7;     /**< frequency seventh derivative (Hz/s^7) */
  REAL8 f8;     /**< frequency eighth derivative (Hz/s^8) */
  REAL8 f9;     /**< frequency ninth derivative (Hz/s^9) */
  
  REAL8 ra;     /**< right ascension (rads) */
  REAL8 dec;    /**< declination (rads) */
  REAL8 pmra;   /**< proper motion in RA (rads/s) */
  REAL8 pmdec;  /**< proper motion in dec (rads/s) */

  REAL8 posepoch; /**< position epoch */
  REAL8 pepoch;   /**< period/frequency epoch */
  
  /* all parameters will be in the same units as used in TEMPO */

  /* Keplerian parameters */
  REAL8 e;      /**< orbital eccentricity */
  REAL8 Pb;     /**< orbital period (days) */
  REAL8 w0;     /**< logitude of periastron (deg) */
  REAL8 x;      /**< projected semi-major axis/speed of light (light secs) */
  REAL8 T0;     /**< time of orbital perisastron as measured in TDB (MJD) */
  
  /* add extra parameters for the BT1P and BT2P models which contain two and
     three orbits respectively (only the first one can be relativistic, so the
     others only have the Keplerian parameters) */
  REAL8 e2, e3;
  REAL8 Pb2, Pb3;
  REAL8 w02, w03;
  REAL8 x2, x3;
  REAL8 T02, T03;

  REAL8 xpbdot;

  /* for low eccentricity orbits (ELL1 model) use Laplace parameters */
  /* (eps1 = e*sin(w), eps2 = e*cos(w)) instead of e, w.             */
  /* see Appendix A, Ch. Lange etal, MNRAS (2001)                    */
  REAL8 eps1;       /**< e*sin(w) */
  REAL8 eps2;       /**< e*cos(w) */
  REAL8 eps1dot;
  REAL8 eps2dot;
  REAL8 Tasc;  /**< time of the ascending node (used rather than T0) */
  UINT4 nEll;  /**< set to zero if have eps time derivs (default) set to 1 if have wdot */

  /* Post-Keplarian parameters */
  /* for Blandford-Teukolsky (BT) model */
  REAL8 wdot;   /**< precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
  REAL8 gamma;  /**< gravitational redshift and time dilation parameter (s)*/
  REAL8 Pbdot;  /**< rate of change of Pb (dimensionless) */
  REAL8 xdot;   /**< rate of change of x(=asini/c) - optional */
  REAL8 edot;   /**< rate of change of e */

  /* for Epstein-Haugan (EH) model */
  REAL8 s; /**< Shapiro 'shape' parameter sin i */
  CHAR *sstr;
  
  /* for Damour-Deruelle (DD) model */
  /*REAL8 r; Shapiro 'range' parameter - defined internally as Gm2/c^3 */
  REAL8 dr;
  REAL8 dth;    /* (10^-6) */
  REAL8 a0, b0; /**< abberation delay parameters */

  /* for DD (General Relativity) (DDGR) - assumes GR is correct model */
  /* we do not need wdot, gamma, Pbdot, s, r, xdot and edot */
  REAL8 M;      /**< M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
  REAL8 m2;     /**< companion mass */

  /* for the DDS model we need the shapmax parameter */
  REAL8 shapmax;
  
  /* orbital frequency coefficients for BTX model (only for one orbit at the
     moment i.e. a two body system) */
  REAL8 *fb; /**< orbital frequency coefficients for BTX model */
  INT4 nfb;     /**< the number of fb coefficients */

  REAL8 px;     /**< pulsar parallax (converted from milliarcsecs to rads) */
  REAL8 dist;   /**< pulsar distance (in kiloparsecs) */

  REAL8 DM;     /**< dispersion measure */
  REAL8 DM1;    /**< first derivative of dispersion measure */

  /* sinusoidal parameters for fitting quasi-periodic timing noise */
  REAL8 *waveCos;
  REAL8 *waveSin;
  REAL8 wave_om;  /**< fundamental frequency timing noise terms */
  REAL8 waveepoch;
  INT4 nwaves;
  
  /* gravitational wave parameters */
  REAL8 h0;     /**< gravitational wave amplitude */
  REAL8 cosiota;/**< cosine of the pulsars orientation angle */
  REAL8 psi;    /**< polarisation angle */
  REAL8 phi0;   /**< initial phase */
  REAL8 Aplus;  /**< 0.5*h0*(1+cos^2iota) */
  REAL8 Across; /**< h0*cosiota */
  
  /* pinned superfluid gw parameters*/
  REAL8 I21;    /**< parameter for pinsf model.**/
  REAL8 I31;    /**< parameter for pinsf model.**/
  REAL8 r;      /**< parameter for pinsf model.**/
  REAL8 lambda; /**< this is a longitude like angle between pinning axis and
                     line of sight */
  REAL8 costheta;  /**< angle between rotation axis and pinning axis */
  
  /* parameters for Kopeikin terms */
  REAL8 daop;   /**< parameter for the Kopeikin annual orbital parallax */
  INT4 daopset; /**< set if daop is set from the par file */
  REAL8 kin;    /**< binary inclination angle */
  INT4 kinset;
  REAL8 kom;
  INT4 komset;
  
  /******** errors read in from a .par file **********/
  REAL8 f0Err;
  REAL8 f1Err;
  REAL8 f2Err;
  REAL8 f3Err;
  REAL8 f4Err;
  REAL8 f5Err;
  REAL8 f6Err;
  REAL8 f7Err;
  REAL8 f8Err;
  REAL8 f9Err;
  
  REAL8 pepochErr;
  REAL8 posepochErr;

  REAL8 raErr;
  REAL8 decErr;
  REAL8 pmraErr;
  REAL8 pmdecErr;

  REAL8 eErr, e2Err, e3Err;
  REAL8 PbErr, Pb2Err, Pb3Err;
  REAL8 w0Err, w02Err, w03Err;
  REAL8 xErr, x2Err, x3Err;
  REAL8 T0Err, T02Err, T03Err;

  REAL8 xpbdotErr;

  REAL8 eps1Err;
  REAL8 eps2Err;
  REAL8 eps1dotErr;
  REAL8 eps2dotErr;
  REAL8 TascErr;

  REAL8 wdotErr;
  REAL8 gammaErr;
  REAL8 PbdotErr;
  REAL8 xdotErr;
  REAL8 edotErr;

  REAL8 sErr;
  REAL8 shapmaxErr;
  
  /*REAL8 rErr; Shapiro 'range' parameter - defined internally as Gm2/c^3 */
  REAL8 drErr;
  REAL8 dthErr;
  REAL8 a0Err, b0Err;

  REAL8 MErr;
  REAL8 m2Err;

  REAL8 *fbErr;

  REAL8 pxErr;
  REAL8 distErr;

  REAL8 DMErr;
  REAL8 DM1Err;

  /* gravitational wave parameters */
  REAL8 h0Err;
  REAL8 cosiotaErr;
  REAL8 psiErr;
  REAL8 phi0Err;
  REAL8 AplusErr;
  REAL8 AcrossErr;
  REAL8 I21Err;
  REAL8 I31Err;
  REAL8 rErr;
  REAL8 lambdaErr;
  REAL8 costhetaErr;
  
  /* timing noise fitting parameters */
  REAL8 wave_omErr;
  
  CHAR *units; /**< The time system used e.g. TDB */
  CHAR *ephem; /**< The JPL solar system ephemeris used e.g. DE405 */
}BinaryPulsarParams;

/** structure containing the Kopeikin terms */
typedef struct
tagKopeikinTerms
{
  REAL8 DK011;
  REAL8 DK012;
  REAL8 DK013;
  REAL8 DK014;
  
  REAL8 DK021;
  REAL8 DK022;
  REAL8 DK023;
  REAL8 DK024;
  
  REAL8 DK031;
  REAL8 DK032;
  REAL8 DK033;
  REAL8 DK034;
  
  REAL8 DK041;
  REAL8 DK042;
  REAL8 DK043;
  REAL8 DK044;
}KopeikinTerms;

/** structure containing the input parameters for the binary delay function */
typedef struct
tagBinaryPulsarInput
{
  REAL8 tb;    /**< Time of arrival (TOA) at the SSB */
  
  EarthState earth; /**< The current Earth state (for e.g. calculating 
                         Kopeikin terms) */
}BinaryPulsarInput;

/** structure containing the output parameters for the binary delay function */
typedef struct
tagBinaryPulsarOutput
{
  REAL8 deltaT; /**< deltaT to add to TDB in order to account for binary */
}BinaryPulsarOutput;

/**** DEFINE FUNCTIONS ****/
/** \brief This function will iteratively calculate the eccentric anomaly from
  * Kelper's equation
  * 
  * The equation is solved using a Newton-Raphson technique and the S9
  * starting value in Odell & Gooding  1986 CeMec 38 307. This is taken from the
  * TEMPO2 code T2model.C 
  */
void
XLALComputeEccentricAnomaly( REAL8 phase, REAL8 ecc, REAL8 *u);

/** \brief This function will compute the effect of binary parameters on the
  * pulsar parallax
  * 
  * This function is based on the terms given in Kopeikin, Ap. J. Lett, 439,
  * 1995. The computation is copied from the KopeikinTerms function in the
  * T2model.C file of TEMPO2.
  */
void XLALComputeKopeikinTerms( KopeikinTerms *kop,
                               BinaryPulsarParams *params,
                               BinaryPulsarInput *input );

/** function to calculate the binary system delay
 */
void
XLALBinaryPulsarDeltaT( BinaryPulsarOutput   *output,
                        BinaryPulsarInput    *input,
                        BinaryPulsarParams   *params );

void
LALBinaryPulsarDeltaT( LALStatus            *status,
                       BinaryPulsarOutput   *output,
                       BinaryPulsarInput    *input,
                       BinaryPulsarParams   *params );

/** function to read in a TEMPO parameter file
 */
void
XLALReadTEMPOParFile( BinaryPulsarParams    *output,
                      CHAR                  *pulsarAndPath );

void
LALReadTEMPOParFile( LALStatus              *status,
                     BinaryPulsarParams    *output,
                     CHAR                  *pulsarAndPath );

/** \brief This function will read in a TEMPO-style parameter correlation matrix
 *
 * This function will read in a TEMPO-style parameter correlation matrix file,
 * which contains the correlations between parameters as fit by TEMPO. An
 * example the format would be:
\verbatim
     RA     DEC    F0
RA   1.000
DEC  0.954  1.000
F0   -0.007 0.124  1.000
\endverbatim
 *
 * In the output all parameter names will sometimes be
 * converted to a more convenient naming convention. If non-diagonal parameter
 * correlation values are +/-1 then they will be converted to be +/-0.99999 to
 * avoid some problems of the matrix becoming singular. The output matrix will
 * have both the upper and lower triangle completed.
 *
 * \param cormat [out] A REAL8 array into which the correlation matrix will be output
 * \param corfile [in] A string containing the path and filename of the
 * 			TEMPO-style correlation matrix file
 */
LALStringVector *XLALReadTEMPOCorFile( REAL8Array *cormat, CHAR *corfile );


/** \brief Convert a string containing an angle in "hours:minutes:seconds" format into radians
 *
 * This function will covert a string containing an angle given in "hours:minutes:seconds"
 * format (e.g. a right ascension) into radians. It requires that the hours value is positive
 * and the minutes and seconds values are between 0 to 60. An example would be:
 *
 * rads = XLALhmsToRads( "12:05:07.765" );
 */
REAL8 XLALhmsToRads( const CHAR *hms );


/** \brief Convert a string containing an angle in "degrees:minutes:seconds" format into radians
 *
 * This function will covert a string containing an angle given in "degrees:minutes:seconds"
 * format (e.g. a declination) into radians. It requires that the minutes and seconds values
 * are between 0 to 60. An example would be:
 *
 * rads = XLALdmsToRads( "-06:52:16.875" );
 */
REAL8 XLALdmsToRads( const CHAR *dms );


/** \deprecated Use XLALdmsToRads() or XLALhmsToRads() instead.
 *
 * A function to convert RA and Dec in format dd:mm:ss.ss or ddmmss.ss into the
 * number of degrees as a float degs is the string containing the
 * dd/hh:mm:ss.sss coords is either ra/RA or dec/DEC.
 */
REAL8
LALDegsToRads(CHAR *degs, const CHAR *coords);

/** Functions for converting times given in Terrestrial time TT, TDB, or TCB in
 * MJD to times in GPS - this is important for epochs given in <tt>.par</tt>
 * files which are in TDB(MJD). TT and GPS are different by a factor of 51.184
 * secs, this is just the historical factor of 32.184 secs between TT and TAI
 * (International Atomic Time) and the other 19 seconds come from the leap
 * seonds added between the TAI and UTC up to the point of definition of GPS
 * time at UTC 01/01/1980.
 */
/*@{*/
REAL8
XLALTTMJDtoGPS(REAL8 MJD);

REAL8
XLALTDBMJDtoGPS(REAL8 MJD);

REAL8
XLALTCBMJDtoGPS(REAL8 MJD);
/*@}*/

/** function to print out all the pulsar parameters read in from a par file
 */
void
PrintPulsarParameters( BinaryPulsarParams params );

#ifdef __cplusplus
}
#endif

#endif /* _BINARYPULSARTIMING_H */
