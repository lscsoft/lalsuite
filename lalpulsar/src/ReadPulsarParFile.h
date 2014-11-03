/*
*  Copyright (C) 2013 Bernd Machenschalk, Jolien Creighton, Matt Pitkin
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
 * \author Matt Pitkin
 * \date 2013
 * \file
 * \ingroup pulsarTODO
 * \brief Functions to read TEMPO pulsar parameter files
 *
 * Here we define a function to read in pulsar parameters from a standard <tt>TEMPO(2)</tt> parameter
 * (.par) file - <tt>XLALReadTEMPOParFile</tt>.
 *
 * We also define a structure for holding pulsar parameters. This is a linked list of structures containing
 * the parameter name (which is always converted to an uppercase string), the parameter value (which can be any
 * data type) and the parameter uncertainty/error (a \c REAL8 value) if given.
 *
 */

#ifndef _READPULSARPARFILE_H
#define _READPULSARPARFILE_H

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
#define READPULSARPARFILEH_ENULLOUTPUT 1

#define READPULSARPARFILEH_MSGENULLOUTPUT "Output was Null"
/*@}*/

#define PULSAR_HASHTABLE_SIZE 512
#define PULSAR_PARNAME_MAX 128
#define DAYSTOSECS 86400.0 /* number of seconds in an SI day */

#define GPS0MJD 44244.0 /* start of GPS time (MJD 44244) */
/* the difference between TDT (or TT) and TAI (see e.g. Eqn 15 of T.-Y. Huang et
 * al, A&A, 220, p. 329, 1989) */
#define TDT_TAI 32.184
/* the difference between TDT/TT and the GPS epoch */
#define GPS_TDT (TDT_TAI + XLAL_EPOCH_GPS_TAI_UTC)

/** An enumerated type for denoting the type of a variable. Several LAL types are supported. */
typedef enum {
  PULSARTYPE_UINT4_t = 0,
  PULSARTYPE_REAL8_t,
  PULSARTYPE_REAL8Vector_t,
  PULSARTYPE_string_t,
  PULSARTYPE_void_ptr_t
} PulsarParamType;


extern size_t PulsarTypeSize[5];


/** \brief The PulsarParam list node structure
 *
 * This structure contains a pulsar parameter defined by the name.
 *
 * This should only be accessed using the accessor functions below.
*/
typedef struct tagPulsarParam {
  CHAR                        name[PULSAR_PARNAME_MAX]; /**< Parameter name */
  void                        *value;                   /**< Parameter value */
  void                        *err;                     /**< Parameter error/uncertainty */
  UINT4                       *fitFlag;                 /**< Set to 1 if the parameter has been fit in the par file */
  PulsarParamType             type;                     /**< Parameter type e.g. REAL8, CHAR, INT4 */
  struct tagPulsarParam       *next;
} PulsarParam;


/** \brief The PulsarParameters structure to contain a set of pulsar parameters
 *
 * This is implemented as a linked list of PulsarParam structures and should
 * only be accessed using the accessor functions below. This is intended as a replacement to the
 * \c BinaryPulsarParams structure.
 */
typedef struct tagPulsarParameters {
  PulsarParam       *head;                              /**< A linked list of \c PulsarParam structures */
  INT4              nparams;                            /**< The total number of parameters in the structure */
  PulsarParam       *hash_table[PULSAR_HASHTABLE_SIZE]; /**< Hash table of parameters */
} PulsarParameters;


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

  CHAR *units; /**< The time system used e.g. TDB */
  CHAR *ephem; /**< The JPL solar system ephemeris used e.g. DE405 */

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
  REAL8 Pb;     /**< orbital period (seconds) */
  REAL8 w0;     /**< logitude of periastron (radians) */
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
  REAL8 M;      /**< M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) (in kg) */
  REAL8 m2;     /**< companion mass (in kg) */

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
  REAL8 lambda; /**< this is a longitude like angle between pinning axis and line of sight */
  REAL8 costheta;  /**< angle between rotation axis and pinning axis */

  /* complex amplitude and phase parameters for l=2, m=1 and 2 harmonics */
  REAL8 C22;
  REAL8 C21;
  REAL8 phi22;
  REAL8 phi21;

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
  REAL8 lambdaErr;
  REAL8 costhetaErr;
  REAL8 C22Err;
  REAL8 C21Err;
  REAL8 phi22Err;
  REAL8 phi21Err;

  /* timing noise fitting parameters */
  REAL8 wave_omErr;

  REAL8 cgw; /**< The speed of gravitational waves as a fraction of the speed of light <tt>LAL_C_SI</tt> */
  REAL8 cgwErr;
}BinaryPulsarParams;


/* DEFINE FUNCTIONS */
/** \brief Get the required parameter value from the \c PulsarParameters structure
 *
 * This function will return a void pointer to the parameter value. This should be cast into the required
 * variable type once returned.
 */
void *PulsarGetParam( const PulsarParameters *pars, const CHAR *name );

/** \brief Get the required parameter error value from the \c PulsarParameters structure
 *
 * This function will return a void pointer to the parameter error value. The parameter error will be either
 * a \c REAL8 or a \c REAL8Vector and should be cast accordingly once returned.
 */
void *PulsarGetParamErr( const PulsarParameters *pars, const CHAR *name );

/** \brief Get the fit flag array for a given parameter from the \c PulsarParameters structure
 *
 * This function will return a \c UINT4 array to the parameter fit flag.
 */
UINT4 *PulsarGetParamFitFlag( const PulsarParameters *pars, const CHAR *name );

/** \brief Return a \c REAL8 parameter error value
 *
 * This function will call \c PulsarGetParamErr for a \c REAL8 parameter and properly cast it for returning.
 */
REAL8 PulsarGetREAL8ParamErr( const PulsarParameters *pars, const CHAR *name );

/** \brief Return a \c REAL8Vector parameter error value
 *
 * This function will call \c PulsarGetParamErr for a \c REAL8Vector parameter and properly cast it for returning.
 */
REAL8Vector *PulsarGetREAL8VectorParamErr( const PulsarParameters *pars, const CHAR *name );

/** \brief Return an individual \c REAL8 value from a \c REAL8Vector parameter
 *
 * This function will call \c PulsarGetParam for a \c REAL8Vector parameter's errors and then return a specific
 * value from within the vector. The input \c name must be of the form e.g. \c FB1, which will get the \c REAL8Vector
 * for the \c FB parameter and return the value from the \c 1 index.
 */
REAL8 PulsarGetREAL8VectorParamErrIndividual( const PulsarParameters *pars, const CHAR *name );

/** \brief Get the required parameter's type
 */
PulsarParamType PulsarGetParamType( const PulsarParameters *pars, const char *name );

/** \brief Return a \c REAL8 parameter
 *
 * This function will call \c PulsarGetParam for a \c REAL8 parameter and properly cast it for returning.
 */
REAL8 PulsarGetREAL8Param( const PulsarParameters *pars, const CHAR *name );

/** \brief Return a string parameter
 *
 * This function will call \c PulsarGetParam for a string parameter and properly cast it for returning.
 */
CHAR *PulsarGetStringParam( const PulsarParameters *pars, const CHAR *name );

/** \brief Return a \c REAL8Vector parameter
 *
 * This function will call \c PulsarGetParam for a \c REAL8Vector parameter and properly cast it for returning.
 */
REAL8Vector *PulsarGetREAL8VectorParam( const PulsarParameters *pars, const CHAR *name );

/** \brief Return an individual \c REAL8 value from a \c REAL8Vector parameter
 *
 * This function will call \c PulsarGetParam for a \c REAL8Vector parameter and then return a specific value
 * from within the vector. The input \c name must be of the form e.g. \c FB1, which will get the \c REAL8Vector
 * for the \c FB parameter and return the value from the \c 1 index.
 */
REAL8 PulsarGetREAL8VectorParamIndividual( const PulsarParameters *pars, const CHAR *name );

/** \brief Add a parameter and value to the \c PulsarParameters structure
 *
 * This function adds a new parameter, and associated value to the \c PulsarParameters structure. If the parameter
 * already exists then the old value is replaced with the new value.
 */
void PulsarAddParam( PulsarParameters *pars, const CHAR *name, void *value, PulsarParamType type );

/** \brief Free all the parameters from a \c PulsarParameters structure */
void PulsarClearParams( PulsarParameters *pars );

/** \brief Remove a given parameter from a \c PulsarParameters structure */
void PulsarRemoveParam( PulsarParameters *pars, const CHAR *name );

/** \brief Set the value of a parameter in the \c PulsarParameters structure
 *
 * Set the value of the parameter given by \c name in the \c PulsarParameters structure. The parameter must already
 * exist in the structure, otherwise it should be added using \c PulsarAddParam().
 */
void PulsarSetParam( PulsarParameters* pars, const CHAR *name, void *value );

/** \brief Set the value of the error of a parameter in the \c PulsarParameters structure
 *
 * Set the value of the error on the parameter given by \c name in the \c PulsarParameters structure.
 * The parameter must already exist in the structure, otherwise it should be added using \c PulsarAddParam().
 * If the .par file contains a 1 before the error value (because that parameter has been included in the
 * TEMPO(2) fitting procedure) then that must be input as the fit flag (this can be a vector for e.g. FB
 * values with multiple parameters, in which case \c nfits will be the number of values in that vector).
 */
void PulsarSetParamErr( PulsarParameters* pars, const CHAR *name, void *value, UINT4 fitFlag, UINT4 nfits, UINT4 len );

/** \brief Check for the existence of the parameter \c name in the \c PulsarParameters structure */
int PulsarCheckParam( PulsarParameters *pars, const CHAR *name );

/** \brief Function to free memory from pulsar parameters */
void PulsarFreeParams( PulsarParameters *par );

/** Conversion functions from units used in TEMPO parameter files into SI units */

/** Convert the input string into a double precision floating point number */
void ParConvToFloat( const CHAR *in, void *out );
/** Convert the input string into a unsigned integer number */
void ParConvToInt( const CHAR *in, void *out );
/** Copy the input string into the output pointer */
void ParConvToString( const CHAR *in, void *out );
/** Convert the input string from degrees to radians */
void ParConvDegsToRads( const CHAR *in, void *out );
/** Convert the input string from milliarcsecs per year to radians per second */
void ParConvMasPerYrToRadPerSec( const CHAR *in, void *out );
/** Convert the input string from seconds to radians */
void ParConvSecsToRads( const CHAR *in, void *out );
/** Convert the input string from arcseconds to radians */
void ParConvArcsecsToRads( const CHAR *in, void *out );
/** Convert the input string from milliarcsecs to radians */
void ParConvMasToRads( const CHAR *in, void *out );
/** Convert the input string from 1/acrseconds to 1/radians */
void ParConvInvArcsecsToInvRads( const CHAR *in, void *out );
/** Convert the input string from days to seconds */
void ParConvDaysToSecs( const CHAR *in, void *out );
/** Convert the binary system parameter from a string to a double, but  make the check (as performed by TEMPO2)
 * that this is > 1e-7 then it's in units of 1e-12, so needs converting by that factor. It also checks if the
 * number is too large (> 10000) and if so sets to to zero.
 */
void ParConvBinaryUnits( const CHAR *in, void *out );
/** Convert the input string from a TT MJD value into a GPS time */
void ParConvMJDToGPS( const CHAR *in, void *out );
/** Convert the input string from degrees per year to radians per second */
void ParConvDegPerYrToRadParSec( const CHAR *in, void *out );
/** Convert the input string from solar masses to kilograms */
void ParConvSolarMassToKg( const CHAR *in, void *out );
/** Convert a right ascension input string in the format "hh:mm:ss.s" into radians */
void ParConvRAToRads( const CHAR *in, void *out );
/** Convert a declination input string in the format "dd:mm:ss.s" into radians */
void ParConvDecToRads( const CHAR *in, void *out );
/** Convert an input string from microseconds into seconds */
void ParConvMicrosecToSec( const CHAR *in, void *out );


/** \brief Read in the parameters from a TEMPO(2) parameter file into a \c PulsarParameters structure
 *
 * This function will read in a TEMPO(2) parameter file into a \c PulsarParameters structure. The structure of this
 * function is similar to those in the TEMPO2 code \c readParfile.C and this is intended to supersede the
 * \c XLALReadTEMPOParFile function.
 *
 * \param pulsarAndPath [in] The path to the pulsar parameter file
 */
PulsarParameters *XLALReadTEMPOParFileNew( const CHAR *pulsarAndPath );

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
 *                      TEMPO-style correlation matrix file
 */
LALStringVector *XLALReadTEMPOCorFile( REAL8Array *cormat, CHAR *corfile );

/** function to print out all the pulsar parameters read in from a par file */
void PrintPulsarParameters( BinaryPulsarParams params );

/** \brief Convert a string containing an angle in "hours:minutes:seconds" format into radians
 *
 * This function will covert a string containing an angle given in "hours:minutes:seconds"
 * format (e.g. a right ascension) into radians. It requires that the hours value is positive
 * and the minutes and seconds values are between 0 to 60. Hours are limited to be between
 * 0 and 24 hours. An example would be:
 *
 * rads = XLALhmsToRads( "12:05:07.765" );
 */
REAL8 XLALhmsToRads( const CHAR *hms );


/** \brief Convert a string containing an angle in "degrees:minutes:seconds" format into radians
 *
 * This function will covert a string containing an angle given in "degrees:minutes:seconds"
 * format (e.g. a declination) into radians. It requires that the minutes and seconds values
 * are between 0 to 60. Degrees are allowed to be any positive of negative integer. An
 * example would be:
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

#ifdef __cplusplus
}
#endif

#endif /* _READPULSARPARFILE_H */
