/*
*  Copyright (C) 2007  Rejean Dupuis
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
 *
 * ### Header \ref FitToPulsar.h ###
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
 * Functions in this package calculate the posterior probability of the following model given the data ({Bk})
 * \f{equation}{
 * y(t;\mathrm{a}) = F_{+}(t;\psi)h_{0} (1 + \cos^{2}\iota)e^{i2\phi_{0}} - 2iF_{\times}(t;\psi) h_{0} \cos\iota e^{i2\phi_{0}}
 * \f}
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
#define FITTOPULSARH_ETIM 9


#define FITTOPULSARH_MSGENULLINPUT "Input was Null"
#define FITTOPULSARH_MSGENULLOUTPUT "Output was Null"
#define FITTOPULSARH_MSGENULLPARAMS "Params was Null"
#define FITTOPULSARH_MSGEVECSIZE "Input vectors were not the same length"
#define FITTOPULSARH_MSGEMESH "Mesh paramaters supplied were invalid"
#define FITTOPULSARH_MSGEVAR "Variance vector in Input had invalid values"
#define FITTOPULSARH_MSGEMAXCHI "The minimum value of chiSquare was greater than INICHISQU"
#define FITTOPULSARH_MSGEDIVZERO "Attempted to divide by zero"
#define FITTOPULSARH_MSGETIM "Problem with time stamps.. there are not M data chunks of length N"


/*@}*/

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

/****** DEFINE NEW STRUCTURES AND TYPES ************/
typedef struct
tagFitInputStudentT
{
  COMPLEX16Vector *B;     /* heterodyned, averaged and resampled data */
  LIGOTimeGPS *t;        /* time stamp for each data point (not necessarily with equal time steps)*/
  INT4	N; /* number of B_k data points that we want to combine together*/  
} FitInputStudentT;

typedef struct
tagCoarseFitInput
{
  COMPLEX16Vector *B;     /* heterodyned, averaged and resampled data */
  COMPLEX16Vector *var;   /* variance of the rFactor points that were averaged */
  LIGOTimeGPS *t;        /* time stamp for each data point (not necessarily with equal time steps)*/
  INT4	N; /* number of B_k data points that we want to combine together*/     
} CoarseFitInput;

typedef struct
tagCoarseFitOutput
{
  REAL8 h0;              /* best fit h0 */
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

typedef struct 
tagPulsarPdfs
{
  REAL4Vector *pdf;       /* pdf of h given data -- marginalizing over other parameters */
  REAL4Vector *cdf;	  /* cdf of h given data -- marginalizing over the other parameters */
  REAL4Vector *pdfPhase;
  REAL4Vector *cdfPhase;
  REAL4Vector *pdfPsi;
  REAL4Vector *cdfPsi;
  REAL4Vector *pdfCosIota;
  REAL4Vector *cdfCosIota;
  /* later also include the 3 contour plots here.. */
} PulsarPdfs;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALCoarseFitToPulsar	(	LALStatus                       *status,
		   		CoarseFitOutput                 *output,
		   		CoarseFitInput 			*input,
		   		CoarseFitParams                 *params); 	
				
void
LALPulsarMarginalize  ( 	LALStatus              *status,
		    		PulsarPdfs             *output,
		    		CoarseFitOutput        *input,
		    		CoarseFitParams        *params );  	
								
void
LALFitToPulsarStudentT	(	LALStatus                       *status,
		   		CoarseFitOutput               *output,
		   		FitInputStudentT		*input,
		   		CoarseFitParams                 *params);							
				
				

#ifdef  __cplusplus
}
#endif

#endif /* _FITTOPULSAR_H */
