/*
*  Copyright (C) 2007 Jolien Creighton
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
 * \author Dupuis, R.J.
 * \file
 * \ingroup pulsarTODO
 *
 * Provides routines for finding the best fit of the measured data to the strain expected from
 * non-precessing pulsar.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/FitToPulsar.h>
 * \endcode
 *
 * The model to be fitted to the data after \c LALFineHeterodyneToPulsar has been applied is
 * \f{equation}{
 * y(t;\mathrm{a}) = F_{+}(t;\psi)h_{0} (1 + \cos^{2}\iota)e^{i2\phi_{0}} - 2iF_{\times}(t;\psi) h_{0} \cos\iota e^{i2\phi_{0}}
 * \f}
 *
 * The reduced set of data points is fitted to this model by minimizing
 * \f$\chi^2\f$ over \f$h_{0}\f$, \f$\phi_{0}\f$, \f$\iota\f$, and \f$\psi\f$.
 *
 * \f{equation}{
 * \chi^2(\mathrm{a}) = \sum_{k}\left|\frac{B_{k} - y(t;\mathrm{a})}{\sigma_{k}^{2}}\right|^2
 * \f}
 *
 * The minimization of \f$\chi^2\f$ is done in two steps <tt>LALCoarseFitToPulsar()</tt> and <tt>LALFineFitToPulsar()</tt>.
 *
 * More documentation soon.
 *
 * ### Error conditions ###
 *
 *
 * ### Structures ###
 *
 *
 * ### Structure \c CoarseFitInput ###
 *
 * This structure stores locked data to be fitted by model.
 * <dl>
 * <dt><tt>COMPLEX16Vector *B</tt></dt><dd> heterodyned, averaged and resampled data</dd>
 * <dt><tt>COMPLEX16Vector *var</tt></dt><dd> variance of the rFactor points that were averaged</dd>
 * <dt><tt>LIGOTimeGPS *t</tt></dt><dd> time stamp for each data point (not necessarily with equal time steps)</dd>
 * </dl>
 *
 * ### Structure \c CoarseFitOutput ###
 *
 * This structure stores the results from the coarse fit of parameters.
 * <dl>
 * <dt><tt>REAL8 h0</tt></dt><dd> best fit h0</dd>
 * <dt><tt>REAL8 eh0[3]</tt></dt><dd> standard error for h0, min standard error, max standard error</dd>
 * <dt><tt>REAL8 cosIota</tt></dt><dd> best fit cosIota</dd>
 * <dt><tt>REAL8 phase</tt></dt><dd> best fit phase</dd>
 * <dt><tt>REAL8 psi</tt></dt><dd> best fit psi (polarization angle)</dd>
 * <dt><tt>REAL8 chiSquare</tt></dt><dd> min value of chi square</dd>
 * <dt><tt>REAL8Vector *mChiSquare</tt></dt><dd> matrix with chi square values</dd>
 * </dl>
 *
 * ### Structure \c CoarseFitParams ###
 *
 * This structure stores the parameters for the coarse fit.
 * <dl>
 * <dt><tt>REAL8 meshH0[3]</tt></dt><dd> min h0, delta h0, number of steps</dd>
 * <dt><tt>REAL8 meshCosIota[3]</tt></dt><dd> min cosIota, delta cosIota, number of steps</dd>
 * <dt><tt>REAL8 meshPhase[3]</tt></dt><dd> min phase, delta phase, number of steps</dd>
 * <dt><tt>REAL8 meshPsi[3]</tt></dt><dd> min psi, delta psi, number of steps</dd>
 * <dt><tt>LALSource pulsarSrc</tt></dt><dd> describes sky position of pulsar</dd>
 * <dt><tt>LALDetector detector</tt></dt><dd> detector</dd>
 * </dl>
 *
 */

#ifndef _FITTOPULSAR_H
#define _FITTOPULSAR_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/Units.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**\name Error Codes */ /*@{*/
#define FITTOPULSARH_ENULLINPUT 1
#define FITTOPULSARH_ENULLOUTPUT 2
#define FITTOPULSARH_ENULLPARAMS 3
#define FITTOPULSARH_EVECSIZE 4
#define FITTOPULSARH_EMESH 5
#define FITTOPULSARH_EVAR 6
#define FITTOPULSARH_EMAXCHI 7
#define FITTOPULSARH_EDIVZERO 8

#define FITTOPULSARH_MSGENULLINPUT "Input was Null"
#define FITTOPULSARH_MSGENULLOUTPUT "Output was Null"
#define FITTOPULSARH_MSGENULLPARAMS "Params was Null"
#define FITTOPULSARH_MSGEVECSIZE "Input vectors were not the same length"
#define FITTOPULSARH_MSGEMESH "Mesh paramaters supplied were invalid"
#define FITTOPULSARH_MSGEVAR "Variance vector in Input had invalid values"
#define FITTOPULSARH_MSGEMAXCHI "The minimum value of chiSquare was greater than INICHISQU"
#define FITTOPULSARH_MSGEDIVZERO "Attempted to divide by zero"

/*@}*/

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

/****** DEFINE NEW STRUCTURES AND TYPES ************/
typedef struct
tagCoarseFitInput
{
  COMPLEX16Vector *B;     /* heterodyned, averaged and resampled data */
  COMPLEX16Vector *var;   /* variance of the rFactor points that were averaged */
  LIGOTimeGPS *t;        /* time stamp for each data point (not necessarily with equal time steps)*/

} CoarseFitInput;

typedef struct
tagCoarseFitOutput
{
  REAL8 h0;              /* best fit h0 */
  REAL8	eh0[3];		 /* standard error for h0, min standard error, max standard error */
  REAL8 cosIota;         /* best fit cosIota */
  REAL8 phase;           /* best fit phase */
  REAL8 psi;             /* best fit psi (polarization angle) */
  REAL8 chiSquare;       /* value of min chi square */
  REAL8Vector *mChiSquare;  /* matrix with chi square values*/
} CoarseFitOutput;

typedef struct
tagCoarseFitParams
{
  REAL8 meshH0[3];	  /* min h0, delta h0, number of steps */
  REAL8 meshCosIota[3];   /* min cosIota, delta cosIota, number of steps */
  REAL8 meshPhase[3];     /* min phase, delta phase, number of steps */
  REAL8 meshPsi[3];       /* min psi, delta psi, number of steps */
  LALSource pulsarSrc;    /* describes sky position of pulsar */
  LALDetector detector;   /* detector */
} CoarseFitParams;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALCoarseFitToPulsar	(	LALStatus                       *status,
		   		CoarseFitOutput                 *output,
		   		CoarseFitInput 			*input,
		   		CoarseFitParams                 *params);

#ifdef  __cplusplus
}
#endif

#endif /* _FITTOPULSAR_H */
