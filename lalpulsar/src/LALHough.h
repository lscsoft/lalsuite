/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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
#ifndef _LALHOUGH_H
#define _LALHOUGH_H
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALHough_h Header LALHough.h
 * \ingroup pkg_pulsarHough
 * \author Alicia Sintes, Badri Krishnan
 *
 * \brief Routines for building and updating the space of partial
 * Hough map derivatives and related functions needed for the construction
 * of  total Hough maps at different frequencies and
 * possible residual spin down parameters.
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/LALHough.h>
 * \endcode
 *
 * As we mention before,
 * the issue is to build histograms, the Hough map (HM), in the
 * parameter space: for each intrinsic frequency  \f$ f_0 \f$, each residual
 * spin-down parameter, and each refined sky location inside the patch.
 * Notice, from the master equation, that the effect of the residual
 * spin-down parameter is just a change in \f$ F_0\f$ , and, at any given
 * time, \f$  F_0 \f$ can be considered constant.
 * Also, the Hough map is a histogram, thus additive. It can be seen as
 * the sum of several partial Hough maps constructed using just one periodogram.
 *
 * Therefore, we can construct  the HM for any  \f$ f_0\f$  and spin-down
 * value by adding together, at different times, partial Hough maps (PHM)
 * corresponding to different \f$ F_0\f$  values (or equivalently, adding their
 * derivatives and then integrating the result).
 *
 * In practice this means that in order to obtain the HM for a given
 * frequency and all possible residual spin-down parameters, we  have to construct
 * a CYLINDER of around the frequency \f$ f_0\f$ .   All of the \e phmd coming
 * from data demodulated with the same parameters.
 * The coordinates of the \e phmd locate the position of the source in
 * the sky, and by summing along different directions inside the cylinder we refine
 * the spin-down value.
 * To analyze another frequency, for all possible spin-down parameters,
 * we just need to add a new line to the cylinder (and remove another one, in a
 * circular buffer)
 * and then proceed making all the possible sums again.
 *
 * For the case of only 1 spin-down parameter we have to sum
 * following straight lines whose slope is related to the grid in the
 * residual spin-down parameter. We can distinguish (at most) as
 * many lines as the number of the different periodograms used.
 */
/*@{*/


/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/SFTfileIO.h>
#include <lal/LUT.h>
#include <lal/PHMD.h>
#include <lal/HoughMap.h>
#include <lal/PulsarDataTypes.h>
#include <lal/NormalizeSFTRngMed.h>

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/**\name Error Codes */
/*@{*/
#define LALHOUGHH_ENULL 1
#define LALHOUGHH_ESIZE 2
#define LALHOUGHH_ESZMM 4
#define LALHOUGHH_EINT  6
#define LALHOUGHH_ESAME 8
#define LALHOUGHH_EFREQ 10
#define LALHOUGHH_EVAL  12
#define LALHOUGHH_EMEM  32

#define LALHOUGHH_MSGENULL "Null pointer"
#define LALHOUGHH_MSGESIZE "Invalid input size"
#define LALHOUGHH_MSGESZMM "Size mismatch"
#define LALHOUGHH_MSGEINT  "Invalid interval"
#define LALHOUGHH_MSGESAME "Input/Output pointers are the same"
#define LALHOUGHH_MSGEFREQ "Invalid frequency"
#define LALHOUGHH_MSGEVAL  "Invalid value"
#define LALHOUGHH_MSGEMEM  "malloc() failed / Out of memory"
/*@}*/


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */

/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

/**
 * This structure stores the frequency indexes of the partial-Hough map
 * derivatives at different time stamps that have to be combined to form a Hough map
 * for a given (null or) residual spin-down parameters
 */
typedef struct tagUINT8FrequencyIndexVector{
  UINT4      length;  /**< number of elements */
  REAL8      deltaF;  /**< frequency resolution */
  UINT8      *data;   /**< the frequency indexes */
} UINT8FrequencyIndexVector;

/**
 * This structure stores a set of frequency-index vectors. Every set
 * corresponds to a different spin-down residual value; there will thus be as many
 * sets as many spin-down residuals one wants to search over with the hough stage.
 */
typedef struct tagUINT8FrequencyIndexVectorSequence{
  UINT4                          length;        /**< number of elements */
  UINT4                          vectorLength;  /**< frequency resolution */
  UINT8FrequencyIndexVector      *freqIndV;     /**< the frequency indexes */
} UINT8FrequencyIndexVectorSequence;

/** This structure contains a vector of peak-grams (for the different time stamps) */
typedef struct tagHOUGHPeakGramVector{
  UINT4             length; /**< number of elements */
  HOUGHPeakGram     *pg;    /**< the Peakgrams */
} HOUGHPeakGramVector;

/** This structure contains a vector of partial look up tables (for the different  time stamps) */
typedef struct tagHOUGHptfLUTVector{
  UINT4            length; /**< number of elements */
  HOUGHptfLUT     *lut;    /**< the partial Look Up Tables */
} HOUGHptfLUTVector;

/** This structure contains a vector of Hough maps */
typedef struct tagHOUGHMapTotalVector{
  UINT4            length; /**< number of elements */
  HOUGHMapTotal    *ht;    /**< the Hough maps */
} HOUGHMapTotalVector;

/**
 * This structure contains a vector sequence of partial-Hough maps
 * derivatives (for different time stamps and different frequencies) representing
 * a circular buffer for the frequency indexes
 */
typedef struct tagPHMDVectorSequence{
  UINT4       nfSize;    /**< number of different frequencies */
  UINT4       length;    /**< number of elements for each frequency */
  UINT8       fBinMin;   /**< frequency index of smallest intrinsic frequency in circular buffer */
  REAL8       deltaF;    /**< frequency resolution */
  UINT4       breakLine; /**< Mark \f$\in\f$[0, \c nfSize) (of the circular buffer) pointing to the starting of the \c fBinMin line */
  HOUGHphmd   *phmd;     /**< the partial Hough map derivatives */
} PHMDVectorSequence;

/** This structure stores the residual spin-down parameters at a given time */
typedef struct tagHOUGHResidualSpinPar{
  REAL8          deltaF;   	/**< Frequency resolution;  df=1/TCOH */
  REAL8          timeDiff; 	/**< \f$T_{\hat N}(t)-T_{\hat N}(\hat t_0)\f$: time difference */
  REAL8Vector    spinRes; 	/**< length: Maximum order of spdwn parameter *data: pointer to residual Spin parameter set fk */
} HOUGHResidualSpinPar;

/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */
void LALHOUGHConstructSpacePHMD (LALStatus            *status,
				 PHMDVectorSequence   *phmdVS,
				 HOUGHPeakGramVector  *pgV,
				 HOUGHptfLUTVector    *lutV
				 );

void LALHOUGHupdateSpacePHMDup (LALStatus            *status,
				PHMDVectorSequence   *phmdVS,
				HOUGHPeakGramVector  *pgV,
				HOUGHptfLUTVector    *lutV
				);

void LALHOUGHupdateSpacePHMDdn (LALStatus            *status,
				PHMDVectorSequence   *phmdVS,
				HOUGHPeakGramVector  *pgV,
				HOUGHptfLUTVector    *lutV
				);

void LALHOUGHConstructHMT  (LALStatus                  *status,
			    HOUGHMapTotal              *ht,
			    UINT8FrequencyIndexVector  *freqInd,
			    PHMDVectorSequence         *phmdVS
			    );

void LALHOUGHComputeFBinMap (LALStatus             *status,
			     UINT8                 *fBinMap,
			     UINT8                 *f0Bin,
			     HOUGHResidualSpinPar  *rs
			     );

void LALHOUGHConstructHMT_W  (LALStatus                  *status,
			      HOUGHMapTotal              *ht,
			      UINT8FrequencyIndexVector  *freqInd,
			      PHMDVectorSequence         *phmdVS
			      );

void LALHOUGHWeighSpacePHMD  (LALStatus            *status,
			      PHMDVectorSequence   *phmdVS,
			      REAL8Vector *weightV
			      );

void LALHOUGHInitializeWeights  (LALStatus            *status,
				 REAL8Vector *weightV
				 );

void LALHOUGHNormalizeWeights  (LALStatus            *status,
				REAL8Vector *weightV
				);

void LALHOUGHComputeAMWeights  (LALStatus         *status,
				REAL8Vector       *weightV,
				LIGOTimeGPSVector *timeV,
				LALDetector       *detector,
				EphemerisData     *edat,
				REAL8             alpha,
				REAL8             delta
				);

void LALHOUGHComputeMultiIFOAMWeights  (LALStatus          *status,
					REAL8Vector        *weightV,
					SFTCatalog         *catalog,
					EphemerisData      *edat,
					REAL8              alpha,
					REAL8              delta
					);


/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _LALHOUGH_H */
