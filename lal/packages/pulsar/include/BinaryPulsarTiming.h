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
 * \ingroup pulsar
 * \brief Functions to calculate binary system time delays and read TEMPO pulsar parameter files
 *
 * $Id$
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

/* LAL code to take into account binary pulsar motion */
/*  for definitions of different models see Taylor and Weisberg (1989)
    and TEMPO software documentation

    Also contains function to read TEMPO .par files to obtain parameters
    and errors on parameters (if available) */
/* <lalVerbatim file="BinaryPulsarTimingHV">
   Author: Pitkin, M. D.
   $Id$
*/
/* Matt Pitkin 29/04/04 */

/* <lalLaTeX>
   \section{Header \texttt{BinaryPulsarTiming.h}}
   label{ss:BinaryPulsarTiming.h}

   Calculates time delay to a signal from a pulsar in a binary system.

   \subsection*{Synopsis}
   \begin{verbatim}
   #include <lal/BinaryPulsarTiming.h>
   \end{verbatim}

   \noindent This header covers routines for calculating the time delay to a
   signal from a pulsar in a binary system. The are also routines for reading
   pulsar data from TEMPO .par file formats.
   </lalLaTeX> */

#ifndef _BINARYPULSARTIMING_H
#define _BINARYPULSARTIMING_H

#include <lal/LALStdlib.h>

NRCSID (BINARYPULSARTIMINGH,"$Id$");

#ifdef __cplusplus
extern "C" {
#endif

/* <lalLaTeX>
   \subsection*{Error conditions}
   </lalLaTeX> */

/* <lalErrTable> */
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

/* </lalErrTable> */

/* <lalLaTeX>
   \subsection*{Structures}
   \idx[Type]{BinaryPulsarParams}
   \idx[Type]{BinaryPulsarInput}
   \idx[Type]{BinaryPulsarOutput}

   \begin{verbatim}
   typedef struct tagBinaryPulsarParams BinaryPulsarParams;
   \end{verbatim}

   This structure contains all the pulsars parameters. The structure does not
   have to be used for a binary pulsar, but can just contain the parameters for
   an isolated pulsar. All parameters are in the same units as given by TEMPO.

   \begin{verbatim}
   typedef struct tagBinaryPulsarInput BinaryPulsarInput;
   \end{verbatim}

   This structure contains the input time at which the binary correction needs
   to be calculated.

   \begin{verbatim}
   typedef struct tagBinaryPulsarOutput BinaryPulsarOutput;
   \end{verbatim}

   This structure contains the binary time delay for the input time.

   </lalLaTeX> */

/**** DEFINE STRUCTURES ****/
/** A structure to contain all pulsar parameters and associated errors */
typedef struct
tagBinaryPulsarParams
{
  CHAR *name;   /**< pulsar name */

  CHAR *model;  /**< TEMPO binary model e.g. BT, DD, ELL1 */

  REAL8 f0;     /**< spin frequency (Hz) */
  REAL8 f1;     /**< frequency first derivative (Hz/s) */
  REAL8 f2;     /**< frequency second derivative (Hz/s^2) */
  REAL8 f3;     /**< frequency third derivative (Hz/s^3) */
  REAL8 f4;     /**< frequency fourth derivative (Hz/s^4) */
  REAL8 f5;     /**< frequency fifth derivative (Hz/s^5) */

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

  REAL8 xpbdot;	/**< (10^-12) */

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
  REAL8 Pbdot;  /**< rate of change of Pb (dimensionless 10^-12) */
  REAL8 xdot;   /**< rate of change of x(=asini/c) - optional (10^-12)*/
  REAL8 edot;   /**< rate of change of e (10^-12)*/

  /* for Epstein-Haugan (EH) model */
  REAL8 s; /**< Shapiro 'shape' parameter sin i */

  /* for Damour-Deruelle (DD) model */
  /*REAL8 r; Shapiro 'range' parameter - defined internally as Gm2/c^3 */
  REAL8 dr;
  REAL8 dth;    /* (10^-6) */
  REAL8 a0, b0; /**< abberation delay parameters */

  /* for DD (General Relativity) (DDGR) - assumes GR is correct model */
  /* we do not need wdot, gamma, Pbdot, s, r, xdot and edot */
  REAL8 M;      /**< M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
  REAL8 m2;     /**< companion mass */

  /* orbital frequency coefficients for BTX model (only for one orbit at the
     moment i.e. a two body system) */
  REAL8 fb[12]; /**> orbital frequency coefficients for BTX model */
  INT4 nfb;     /**> the number of fb coefficients */

  REAL8 px;     /**> pulsar parallax (in milliarcsecs) */
  REAL8 dist;   /**> pulsar distance (in kiloparsecs) */

  REAL8 DM;     /**> dispersion measure */
  REAL8 DM1;    /**> first derivative of dispersion measure */

  /******** errors read in from a .par file **********/
  REAL8 f0Err;
  REAL8 f1Err;
  REAL8 f2Err;
  REAL8 f3Err;
  REAL8 f4Err;
  REAL8 f5Err;

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

  /*REAL8 rErr; Shapiro 'range' parameter - defined internally as Gm2/c^3 */
  REAL8 drErr;
  REAL8 dthErr;
  REAL8 a0Err, b0Err;

  REAL8 MErr;
  REAL8 m2Err;

  REAL8 fbErr[12];

  REAL8 pxErr;
  REAL8 distErr;

  REAL8 DMErr;
  REAL8 DM1Err;
}BinaryPulsarParams;

/** structure containing the input parameters for the binary delay function */
typedef struct
tagBinaryPulsarInput
{
  REAL8 tb;    /**< Time of arrival (TOA) at the SSB */
}BinaryPulsarInput;

/** structure containing the output parameters for the binary delay function */
typedef struct
tagBinaryPulsarOutput
{
  REAL8 deltaT;	/**< deltaT to add to TDB in order to account for binary */
}BinaryPulsarOutput;

/* <lalLaTeX>
   \newpage\input{BinaryPulsarTimingC}
   </lalLaTeX> */

/**** DEFINE FUNCTIONS ****/
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

/** A function to convert RA and Dec in format dd:mm:ss.ss or ddmmss.ss into the
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
REAL8
LALTTMJDtoGPS(REAL8 MJD);

REAL8
LALTDBMJDtoGPS(REAL8 MJD);

REAL8
LALTCBMJDtoGPS(REAL8 MJD);

/** function to print out all the pulsar parameters read in from a par file
 */
void
PrintPulsarParameters( BinaryPulsarParams params );

#ifdef __cplusplus
}
#endif

#endif /* _BINARYPULSARTIMING_H */
