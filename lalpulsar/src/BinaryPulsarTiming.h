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
 * The main function in this code - <tt>XLALBinaryPulsarDeltaT</tt> - is for
 * calculating the time delay on a pulsar signal caused by its orbit within a
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

#include <lal/ReadPulsarParFile.h>

#ifdef __cplusplus
extern "C" {
#endif


/**\name Error Codes */ /*@{*/
#define BINARYPULSARTIMINGH_ENULLINPUT 1
#define BINARYPULSARTIMINGH_ENULLOUTPUT 2
#define BINARYPULSARTIMINGH_ENULLPARAMS 3
#define BINARYPULSARTIMINGH_ENULLBINARYMODEL 4
#define BINARYPULSARTIMINGH_EFAIL 5
#define BINARYPULSARTIMINGH_ENAN 6

#define BINARYPULSARTIMINGH_MSGENULLINPUT "Input was Null"
#define BINARYPULSARTIMINGH_MSGENULLOUTPUT "Output was Null"
#define BINARYPULSARTIMINGH_MSGENULLPARAMS "Params was Null"
#define BINARYPULSARTIMINGH_MSGNULLBINARYMODEL "Binary model is Null or not specified - you should\
not be in the binary timing routine"
#define BINARYPULSARTIMINGH_MSGEFAIL "Time delay computation failed"
#define BINARYPULSARTIMINGH_MSGENAN "Output is NaN!"

/*@}*/


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
/**
 * \brief This function will iteratively calculate the eccentric anomaly from
 * Kelper's equation
 *
 * The equation is solved using a Newton-Raphson technique and the S9
 * starting value in Odell & Gooding  1986 CeMec 38 307. This is taken from the
 * TEMPO2 code T2model.C
 */
void
XLALComputeEccentricAnomaly( REAL8 phase, REAL8 ecc, REAL8 *u);

/**
 * \brief This function will compute the effect of binary parameters on the
 * pulsar parallax
 *
 * This function is based on the terms given in Kopeikin, Ap. J. Lett, 439,
 * 1995. The computation is copied from the KopeikinTerms function in the
 * T2model.C file of TEMPO2.
 */
void XLALComputeKopeikinTerms( KopeikinTerms *kop,
                               BinaryPulsarParams *params,
                               BinaryPulsarInput *input );

/**
 * function to calculate the binary system delay
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


#ifdef __cplusplus
}
#endif

#endif /* _BINARYPULSARTIMING_H */
