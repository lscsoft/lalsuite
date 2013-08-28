/*
*  Copyright (C) 2007 Chad Hanna, David Churches, Duncan Brown, Jolien Creighton, Benjamin Owen, B.S. Sathyaprakash, Anand Sengupta, Craig Robinson , Thomas Cokelaer, Evan Ochsner
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

#ifndef _LALINSPIRALBANK_H
#define _LALINSPIRALBANK_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/LALNoiseModels.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup LALInspiralBank_h
 * \author Churches, D.K. and Sathyaprakash, B.S., Cokelaer, T.
 *
 * \brief %Header file for the template placement codes.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/LALInspiralBank.h>
 * \endcode
 *
 * This header file covers routines that are used in template placement.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define LALINSPIRALBANKH_ENULL      1	/**< Null pointer */
#define LALINSPIRALBANKH_EMEM       2	/**< Memory allocation failure */
#define LALINSPIRALBANKH_ECHOICE    3	/**< Invalid choice for an input parameter */
#define LALINSPIRALBANKH_EDIV0      4	/**< Division by zero */
#define LALINSPIRALBANKH_ESIZE      8	/**< Invalid input range */
#define LALINSPIRALBANKH_EFRANGE    16	/**< Limits outside range of frequency series */
#define LALINSPIRALBANKH_EORDER     32	/**< Inappropriate PN order */
#define LALINSPIRALBANKH_EGRIDSPACING 64	/**< Inappropriate grid spacing parameter [SquareNotOriented or Hexagonal] */
#define LALINSPIRALBANKH_EHEXAINIT 128	/**< Empty bank. abnormal behaviour in HexaBank generation. */
#define LALINSPIRALBANKH_EFCUT      5	/**< Inappropriate cutoff frequency [SchwarzISCO, BKLISCO, LightRing, ERD, FRD or LRD] */
#define LALINSPIRALBANKH_EFHIGH     6	/**< Final frequency is less than the low frequency cutoff. */
#define LALINSPIRALBANKH_ENUMFCUT   7	/**< Number of fcut must be greater or equal to 1 */
/*@}*/

/** \cond DONT_DOXYGEN */
#define LALINSPIRALBANKH_MSGENULL   "Null pointer"
#define LALINSPIRALBANKH_MSGEMEM    "Memory allocation failure"
#define LALINSPIRALBANKH_MSGECHOICE "Invalid choice for an input parameter"
#define LALINSPIRALBANKH_MSGEDIV0   "Division by zero"
#define LALINSPIRALBANKH_MSGESIZE   "Invalid input range"
#define LALINSPIRALBANKH_MSGEFRANGE "Limits outside range of frequency series"
#define LALINSPIRALBANKH_MSGEORDER  "Inappropriate PN order"
#define LALINSPIRALBANKH_MSGEGRIDSPACING "Inappropriate grid spacing parameter [SquareNotOriented or Hexagonal]"
#define LALINSPIRALBANKH_MSGEHEXAINIT "Empty bank. abnormal behaviour in HexaBank generation."
#define LALINSPIRALBANKH_MSGEFCUT   "Inappropriate cutoff frequency [SchwarzISCO, BKLISCO, LightRing, ERD, FRD or LRD]"
#define LALINSPIRALBANKH_MSGEFHIGH  "Final frequency is less than the low frequency cutoff."
#define LALINSPIRALBANKH_MSGENUMFCUT "Number of fcut must be greater or equal to 1"
/** \endcond */

/** UNDOCUMENTED */
typedef enum
{
  disable,
  enable
}
ComputeMoments;


/**
 * Choose templates either in the \f$(\tau_0,\tau_2)\f$ or \f$(\tau_0,\tau_3)\f$
 * space.  This is one of the members of the InspiralCoarseBankIn structure.
 *
 * This enum allows users to choose template bank either in the \f$(\tau_0, \tau_2)\f$
 * space of chirptimes (the choice made by #Tau0Tau2) or in the
 * \f$(\tau_0, \tau_3)\f$ space of chirptimes (the choice made by #Tau0Tau3).
 * This was implemented in releases before May 25, 2002. On May 25 we migrated to a
 * new, slightly faster, computation of the metric in which, at present, only the
 * choice \c Tau0Tau3 can be made. Since October 2003 a new choice \c Psi0Psi3
 * was added to handle BCV templates. In November 2007 two new choices were addded:
 * \c PTFIntrinctic is a PTF metric in only the intrinsic parameters (a \f$4
 * \times 4\f$ matrix), and \c PTFFull is the PTF metric in the full parameter
 * space (intrinsic and extrinsic parameters).
 */
typedef enum
{
  Tau0Tau2,	/**< \f$(\tau_0, \tau_2)\f$ space of chirptimes */
  Tau0Tau3,	/**< \f$(\tau_0, \tau_3)\f$ space of chirptimes */
  Psi0Psi3,	/**< for BCV templates */
  PTFIntrinsic,	/**< a PTF metric in only the intrinsic parameters (a \f$4\times 4\f$ matrix) */
  PTFFull	/**< PTF metric in the full parameter space (intrinsic and extrinsic parameters). */
}
CoordinateSpace;

/**
 * This enum is set by the user to specify the type of placement requested. It can be
 * <tt>Square, Hexagonal, SquareNotOriented, HexagonalNotOriented, S2BCV</tt>. The two first
 * align the ellipse along the eigen-vectors whereas the two next do not. The last is a
 * square placement which was being used during S2 and is therefore obsolete and should
 * not be used (feel free to remove it). Historically, we used the \c SquareNotOriented
 * placement until S4. Then, in S5, we switched to the \c Hexagonal placement,
 * which should be used for future searches.
 */
typedef enum
{
  SquareNotOriented,	/**< UNDOCUMENTED */
  Square,		/**< UNDOCUMENTED */
  HexagonalNotOriented,	/**< UNDOCUMENTED */
  Hexagonal,		/**< UNDOCUMENTED */
  HybridHexagonal,	/**< UNDOCUMENTED */
  S2BCV			/**< UNDOCUMENTED */
}
GridSpacing;

/**
 * This enum can take the following values <tt>In, Out, Below, Edge, Above</tt> and is used
 * \e only by the Hexagonal placement.  It simply specifies
 * the place of a point with respect to the parameter space. Edge, means that the ellipse
 * covers two boundaries(upper and lower).
 */
typedef enum
{
  In,		/**< UNDOCUMENTED */
  Above,	/**< UNDOCUMENTED */
  Below,	/**< UNDOCUMENTED */
  Out,		/**< UNDOCUMENTED */
  Edge		/**< UNDOCUMENTED */
}
Position;

/**
 * This enum is set to true or false, it is just a boolean variable for the
 * purpose of BCV placement but can be used in an other context.
 */
typedef enum
{
  False,
  True
}
InsidePolygonEnum;

/**
 * This enum is either <tt>fertile,sterile</tt>, and is a boolean expression used \e only
 * by the Hexagonal placement.
 */
typedef enum
{
  Sterile,
  Fertile
}
Generation;

/**
 * An enum that appears in the \c InspiralCoarseBankIn structure
 * which fixes the way templates are chosen: The choice
 * \c MinComponentMassMaxTotalMass means the minimum of the
 * component masses will be given by \c mMin and maximum total
 * mass is given by \c MMax of the \c InspiralBankCoarseIn structure.
 * The choice \c MinMaxComponentMass means the minimum of the
 * components masses will be again fixed by \c mMin and the
 * maximum of the component masses is fixed by \c mMax of the
 * \c InspiralCoarseIn structure below.
 */
typedef enum
{
  MinComponentMassMaxTotalMass,
  MinMaxComponentMass,
  MinMaxComponentTotalMass
}
InspiralBankMassRange;


/**
 * An enum that lists all the formulas that can be used to specify an upper
 * frequency cutoff.
 */
typedef enum
{
  FreqCut_SchwarzISCO,	/**<     the innermost stable circular orbit (ISCO) for a test particle orbiting a Schwarzschild black hole */
  FreqCut_BKLISCO,	/**<     a mass ratio dependent ISCO derived from estimates of the final spin of a merged black found in a
                         *       paper by Buonanno, Kidder, and Lehner (arXiv:0709.3839)
                         */
  FreqCut_LightRing,	/**<     the unstable circular orbit for photons orbiting a Schwarzschild black hole */
  FreqCut_ERD,		/**<     an effective ringdown frequency studied in Pan et al (arXiv:0704.1964) that was found to give good
                         *       fit between stationary-phase templates and  numerical relativity waveforms
                         */
  FreqCut_FRD,		/**<     the "Fundamental RingDown" frequency which is calculated from the Berti, Cardoso and Will (arXiv:gr-qc/0512160)
                         *       value for the \f$\omega_{220}\f$ QNM frequency using mass ratio dependent fits to the final BH mass and spin from Buonanno et al (arXiv:0706.3732)
                         */
  FreqCut_LRD		/**<     the "Lorentzian RingDown" frequency = 1.2*FRD which captures part of the Lorentzian tail from the decay of the QNMs */
}
FreqCut;


/**
 * Structure to store metric at various points the signal manifold.
 * We store the diagonalized metric together with the angle theta
 * between the \f$\tau_0\f$-axis and the semi-major axis of the ambiguity ellipse.
 */
typedef struct
tagInspiralMetric
{
  REAL8            G00;		/**< 00-component of the metric in \f$(\tau_0,\tau_{2(3)})\f$ coordinates */
  REAL8            G11;		/**< 11-component of the metric in \f$(\tau_0,\tau_{2(3)})\f$ coordinates */
  REAL8            G01;		/**< 01-component of the metric in \f$(\tau_0,\tau_{2(3)})\f$ coordinates */

  REAL8            g00;		/**< 00-component of the diagonalised metric */
  REAL8            g11;		/**< 11-component of the diagonalised metric */
  REAL8            theta;	/**< Angle from tau0 to semi-major axis of the ellipse */

  REAL4            Gamma[10];	/**< 3d metric co-efficients in \f$(t_C, \tau_0,\tau_{2(3)})\f$ coordinates;
                                 * Gamma[6] is a vector that stores the upper triangular part of the metric in
                                 * the space of parameters; For time domain searches, Gamma[0,...,5] stores
                                 * the following information :<br>
                                 *    Gamma[0] -> (tc,tc) metric component<br>
                                 *    Gamma[1] -> (tc,t0) metric component<br>
                                 *    Gamma[2] -> (tc,t3) metric component<br>
                                 *    Gamma[3] -> (t0,t0) metric component<br>
                                 *    Gamma[4] -> (t0,t3) metric component<br>
                                 *    Gamma[5] -> (t3,t3) metric component<br>
                                 * For spinBCV searches, (in 4 dimensions) Gamma[0,...,9] would be required
                                 */
  CoordinateSpace  space;	/**< The enum describing the coordinate space in which the metric is computed */
}
InspiralMetric;


/**
 * A grid of inspiral templates (ie a template list).
 * Structure returned by the coarse and fine bank generation routines.
 * Currently we generate an array of type \c InspiralTemplateList
 * which contains the coordinate markers (the parameter structure
 * \c InspiralTemplate defined in the \c inspiral package)
 * and the metric at each of those points. There is a desire to make this
 * a truly linked list at some time in the future.
 */
typedef struct
tagInspiralTemplateList
{
  INT4              ID;		/**< An unique integer ID of the template */
  InspiralTemplate  params;	/**< Value of the parameters at the lattice point */
  InspiralMetric    metric;	/**< metric at the lattice point */
  UINT4             nLayer;	/**< UNDOCUMENTED */
  struct tagInspiralTemplateList *next;	/**< pointer to next lattice point; but this is currently not filled by the bank code */
}
InspiralTemplateList;

/**
 * This is a structure needed in the inner workings of the \c LALInspiralHexagonalBank code.
 * It contains some part of CoarseBankIn and some other standard parameters.  It provides the
 * parameter space boundaries with the minimum and maximum values of mass parameters, the
 * minimal match, the space, massRange and gridSpacing parameter.
 */
typedef struct
tagHexaGridParam
{
  REAL4 x0Min;
  REAL4 x1Min;
  REAL4 x0Max;
  REAL4 x1Max;
  REAL4 mm;
  REAL4 mMin;
  REAL4 mMax;
  REAL4 etaMin;
  REAL4 MMin;
  REAL4 MMax;
  REAL4 fLower;
  GridSpacing 		gridSpacing;
  InspiralBankMassRange massRange;
  CoordinateSpace       space;
}
HexaGridParam;

/**
 * This is a structure needed in the inner workings of the \c LALInspiralHexagonalBank code.
 *
 * This structure checks the status of the placement. \c fertile tells if the
 * placement is still evolving or not. \c nTemplateMax is the number of maximum templates allowed,
 * which can be resized. And \c nTemplate is the number of template set. nTemplate can not
 * be higher than nTemplateMax.
 */
typedef struct
tagCellEvolution
{
  INT4 nTemplateMax;
  INT4 nTemplate;
  INT4 fertile;
}
CellEvolution;


/**
 * This is a structure needed in the inner workings of the \c LALInspiralHexagonalBank code.
 *
 * Similarly to the square placement, which uses InspiralList, we used a
 * linked list for the hexagonal placement. A different structure has been
 * implemented so as to simplify the complexity of the algorithm. It also set
 * an id to each cell which has been created. This id is unique to each
 * cell/template.
 */
typedef struct
tagCellList
{
  INT4 id;
  struct tagCellList *next;
}
CellList;

/**
 * This is a structure needed in the inner workings of the \c LALInspiralHexagonalBank code.
 *
 * Each cell is defined by this structure, which contains the position of
 * each cell in the tau0/tau3 parameter space, the metric at that point, and various
 * information such as the status of the cell. Is it still fertile ? what is its position
 * with respect to the parameter space and so on. child is a 6-length array with a link
 * to the 6 templates (hexagonal) around the current template that we are dealing with.
 */
typedef struct
tagInspiralCell
{
  INT4  ID;
  INT4  in;
  INT4  child[6];
  REAL4 t0;
  REAL4 t3;
  REAL4 dx0;
  REAL4 dx1;
  Generation status;
  Position position;
  Position RectPosition[5];
  InspiralMetric metric;
}
InspiralCell;


/**
 * This is a structure needed in the inner workings
 * of the \c LALInspiralCreateCoarseBank code.
 */
typedef struct
tagInspiralBankParams
{
  INT4           nparams;	/**< Number of parameters (currently fixed at 2, so this is as of now unused) */
  REAL8          minimalMatch;	/**< UNDOCUMENTED */
  REAL8          x0;		/**< the first coordinate, chosen to be always \f$\tau_0\f$ */
  REAL8          x1;		/**< the second coordinate, chosen to be either \f$\tau_2\f$ or \f$\tau_3\f$ */
  REAL8          dx0;		/**< increment in the x0-direction */
  REAL8          dx1;		/**< increment in the x1-direction */
  REAL8          x0Min;		/**< minimum value of the first coordinate as defined by the search region */
  REAL8          x0Max;		/**< maximum value of the first coordinate as defined by the search region */
  REAL8          x1Min;		/**< minimum value of the second coordinate as defined by the search region */
  REAL8          x1Max;		/**< maximum value of the second coordinate as defined by the search region */
  InspiralMetric *metric;	/**< pointer to the metric at the current location */
}
InspiralBankParams;


/**
 * Input for choosing a template bank. This is the structure that must
 * be filled by a routine calling the code \c InspiralCreateCoarseBank() or \c InspiralCreateBCVBank().
 * Unless BCV template bank is needed (that is, \c InspiralCreateBCVBank())  then one can ignore the
 * parameters <tt>psi0Min, psi0Max, psi3Min, psi3Max, alpha, numFcutTemplates.</tt>
 */
typedef struct
tagInspiralCoarseBankIn
{
  InspiralBankMassRange         massRange;	/**< enum that determines whether templates should be chosen using fixed ranges for component masses or
                                                 * to use minimum component mass and maximum totalmass
                                                 */
  CoordinateSpace               space;		/**< enum that decides whether to use \f$(\tau_0,\tau_2)\f$ or \f$(\tau_0,\tau_3)\f$ in constructing the template bank */
  REAL8                         mMin;		/**< minimum mass of components to search for */
  REAL8                         mMax;		/**< maximum mass of components to search for */
  REAL8                         MMax;		/**< alternatively, maximum total mass of binary to search for */
  REAL8                         MMin;		/**< UNDOCUMENTED */
  REAL8                         alpha;		/**< the BCV amplitude correction parameter */
  REAL8                         psi0Min;	/**< minimum value of the parameter \f$\psi_0\f$ */
  REAL8                         psi0Max;	/**< maximum value of the parameter \f$\psi_0\f$ */
  REAL8                         psi3Min;	/**< minimum value of the parameter \f$\psi_3\f$ */
  REAL8                         psi3Max;	/**< maximum value of the parameter \f$\psi_3\f$ */
  REAL8                         mmCoarse;	/**< Coarse grid minimal match */
  REAL8                         mmFine;		/**< Fine grid minimal match */
  REAL8                         fLower;		/**< Lower frequency cutoff */
  REAL8                         fUpper;		/**< Upper frequency cutoff */
  REAL8                         tSampling;	/**< Sampling rate */
  REAL8                         etamin;		/**< minimum value of eta in our search */
  REAL8				betaMin;	/**< UNDOCUMENTED */
  REAL8				betaMax;	/**< UNDOCUMENTED */
  REAL8                         chiMin;		/**< UNDOCUMENTED */
  REAL8                         chiMax;		/**< UNDOCUMENTED */
  REAL8                         kappaMin;	/**< UNDOCUMENTED */
  REAL8                         kappaMax;	/**< UNDOCUMENTED */
  INT4                          nPointsChi;	/**< UNDOCUMENTED */
  INT4                          nPointsKappa;	/**< UNDOCUMENTED */
  REAL8FrequencySeries          shf;		/**< Frequency series containing the PSD */
  UINT4				ShMaxSz;	/**< Maximum size of the power spectral density array for use in
                                                 * the computation of the metric in SBBH; typical values that
                                                 * assures that the code runs quickly are 1024-8192
                                                 */
  UINT4				iseed;		/**< See for random number generation in RandomBank algorithm */
  UINT4				nTIni;		/**< nTIni is an estimate for the number of templates that might
                                                 * be required; this is used in the random bank generation
                                                 * routine with a seed number of templates = nTIni*sqrt(nTIni)
                                                 */
  INT4                          iflso;		/**< iflso is an integer that tells whether to compute the moments
                                                 * using an upper limit defined by flso; this is not used anywhere
                                                 * at the moment
                                                 */

  INT4                          spinBank;	/**< spinBank=0:use Owen+Hanna bank,
                                                 * spinBank=1:use extended bank by AEI/Cardiff/Osaka,
                                                 * spinBank=2:use random bank algorithm
                                                 */
  UINT4                         numFcutTemplates;  /**< Number of templates required in the fCut (upper cutoff)
                                                    * dimension and the value of upper and lower cutoffs
                                                    */
  REAL4				HighGM;		/**< UNDOCUMENTED */
  REAL4				LowGM;		/**< UNDOCUMENTED */
  GridSpacing                   gridSpacing;	/**< Type of gridspacing required */

  /* post-Newtonian order( phase), approximant, and amplitude PN order */
  LALPNOrder                    order;		/**< Post-Newtonian order of the waveform */
  Approximant                   approximant;	/**< Approximant of the waveform */
  LALPNOrder                    ampOrder;

  /* parameters for different/multiple freq cutoffs */
  INT4                          numFreqCut;	/**< Number of different upper frequency cutoffs (spaced evenly between minFreqCut and maxFreqCut) to use when creating a template bank */
  FreqCut                       maxFreqCut;	/**< largest upper frequency cutoff to use */
  FreqCut                       minFreqCut;	/**< smallest upper frequency cutoff to use */

  InsidePolygonEnum             insidePolygon;	/**< UNDOCUMENTED */
  ComputeMoments                computeMoments;	/**< ComputeMoments tells whether to re-compute the moments
                                                 * using an upper limit defined by flso; This is done after
                                                 * the template bank is gnerated
                                                 */
}
InspiralCoarseBankIn;

/**
 * Inputs to the function that computes the moments of the PSD.
 * The moment is defined as:
 * \f[I(p) \equiv \int_{x_\textrm{min}}^{x_\textrm{max}}
 * \frac{x^{-p}}{S_h(x)} dx,\f]
 * where \f$x=f/f_0\f$ is a scaled frequency, \f$f_0\f$
 * being a fiducial frequency, taken in these routines
 * as the user supplied lower cutoff of the detector
 * response.
 */
typedef struct
tagInspiralMomentsIn
{
  REAL8                xmin;	/**< lower limit of the integral \f$x_\textrm{min}\f$ */
  REAL8                xmax;	/**< upper limit of the integral \f$x_\textrm{max}\f$ */
  REAL8                ndx;	/**< index \f$p\f$ (without the negative sign) in the moment integral as above */
  REAL8	               norm;	/**< norm to be used in computing the moment, the returned value is the above integral divided by the norm */
  REAL8FrequencySeries *shf;	/**< the frequency series containing the noise psd */
}
InspiralMomentsIn;


/**
 * Structure needed by the function \c LALInspiralCreateFineBank.
 * which computes a finer mesh around a given lattice point
 * using the value of the fine-mesh minimal match, coarse-mesh
 * minimal match and the metric at the current lattice point.
 */
typedef struct
tagInspiralFineBankIn
{
  InspiralTemplateList templateList;	/**< A list containing all the fine-mesh templates */
  InspiralCoarseBankIn coarseIn;	/**< input structure that contains useful necessary parameters to construct a fine-mesh  */
}
InspiralFineBankIn;


/**
 * Parameter structure that holds the moments of the PSD and other useful
 * constants required in the computation of the metric.
 */
typedef struct
tagInspiralMomentsEtc
{
/**
 * \name Expansion coefficients
 * Coefficients in the expansion of the phase of the Fourier transform of an inspiral waveform computed
 * in the stationary phase approximation. See documentation under the function \c LALInspiralComputeMetric() in this
 * Section for a description of these coefficients
 */
  /*@{*/
  REAL8 a01, a21, a22, a31, a41, a42, a43;
  /*@}*/
  REAL8 j[18];	/**< The required moments are all computed once and
                 * stored in this array; The required moments are from J(1) to J(17)
                 * (except J(2), J(3) and J(16) that are not required at 2PN order,
                 * however, they are computed since future extensions, planned in the
                 * near future, will require them); However, in C we need an array size
                 * 18 to use an array that has an index 18; To ease the notation we have
                 * therefore defined an over sized (by one element) array
                 */
}
InspiralMomentsEtc;

/** UNDOCUMENTED */
typedef struct
tagInspiralMomentsEtcBCV
{
  REAL8 n0, n15;
  REAL8 j[9];
  REAL8 i[23];
  REAL8 alpha;
  REAL8 fcut;

  REAL8 M1[2][2];
  REAL8 M2[2][2];
  REAL8 M3[2][2];
}
InspiralMomentsEtcBCV;


/**
 * Input structure to function LALRectangleVertices()
 */
typedef struct
tagRectangleIn
{
  REAL8 x0, y0, dx, dy, theta;
}
RectangleIn;

/**
 * Output structure to function LALRectangleVertices().
 */
typedef struct
tagRectangleOut
{
  REAL8 x1;
  REAL8 y1;
  REAL8 x2;
  REAL8 y2;
  REAL8 x3;
  REAL8 y3;
  REAL8 x4;
  REAL8 y4;
  REAL8 x5;
  REAL8 y5;
}
RectangleOut;

/** UNDOCUMENTED */
typedef struct
tagHexagonOut
{
  REAL8 x1;
  REAL8 y1;
  REAL8 x2;
  REAL8 y2;
  REAL8 x3;
  REAL8 y3;
  REAL8 x4;
  REAL8 y4;
  REAL8 x5;
  REAL8 y5;
  REAL8 x6;
  REAL8 y6;
  REAL8 x7;
  REAL8 y7;
}
HexagonOut;

/** UNDOCUMENTED */
typedef struct
tagPRIN
{
REAL4 ct;
REAL4 b;
}
PRIN;

/*@}*/ /* end:LALInspiralBank_h */

/* ---------- Function prototypes ---------- */

void
LALInspiralCreateCoarseBank (
    LALStatus              *status,
    InspiralTemplateList   **list,
    INT4                   *nlist,
    InspiralCoarseBankIn   bankIn
    );

void
LALInspiralCreatePNCoarseBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );




void
LALInspiralCreateBCVBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralCreateFlatBankS3S4 (
    LALStatus            *status,
    REAL4VectorSequence  *list,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn coarseIn
    );

void
LALExcludeTemplate(
    LALStatus            *status,
    INT4 *valid,
    InspiralBankParams   *bankParams,
    REAL4 x,
    REAL4 y
    );

void
LALInspiralBCVBankFcutS3S4 (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                *NList,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralBCVFcutBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                *NList,
    InspiralCoarseBankIn coarseIn
    );


void
LALEmpiricalPSItoMassesConversion(
    LALStatus            *status,
    InspiralTemplate    *params,
    UINT4               *valid,
    REAL4               lightring
);

void
LALPSItoMasses(
    LALStatus            *status,
    InspiralTemplate    *params,
    UINT4               *valid,
    REAL4               thisfreq
);

void
LALInspiralBCVRegularFcutBank (
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                *NList,
    InspiralCoarseBankIn coarseIn);

void
LALNudgeTemplatesToConstantTotalMassLine(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 nlist,
    InspiralCoarseBankIn coarseIn
    );




REAL4
XLALInspiralTau3FromTau0AndEqualMassLine(
    REAL4               tau0,
    REAL4               fL
    );

REAL4
XLALInspiralTau3FromNonEqualMass(
  REAL4               	m1,
  REAL4 		m2,
  REAL4			fL
);

REAL4
XLALInspiralTau0FromMEta(
  REAL4              	M,
  REAL4 		eta,
  REAL4			fL
);


REAL8
XLALInspiralBissectionLine (
  REAL8 x,
  REAL8 fL,
  REAL8 mMin,
  REAL8 mMax);

REAL8
XLALInspiralMFromTau0AndNonEqualMass(
  REAL8 tau0,
  REAL8 extremMass,
  REAL8 fL);

void
LALInspiralSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable    **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    );

#if 0
void
LALInspiralSpinBankBoundary(
    LALStatus            *status,
    NDTemplateBankInput  *input,
    NDTemplateBankOutput *output,
    INT2                 *flag
    );

void
LALInspiralSpinBankMetric(
    LALStatus           *status,
    NDTemplateBankInput *input,
    REAL4Array          *metric
    );
#endif

void
LALInspiralBankGeneration(
    LALStatus            *status,
    InspiralCoarseBankIn *in,
    SnglInspiralTable    **out,
    INT4                 *count
    );

void
LALInspiralCreateFlatBank (
    LALStatus            *status,
    REAL4VectorSequence  *list,
    InspiralBankParams   *bankParams
    );

void
LALInspiralCreateFineBank (
    LALStatus              *status,
    InspiralTemplateList   **outlist,
    INT4                   *nlist,
    InspiralFineBankIn     fineIn
    );

void
LALInspiralComputeMetric (
    LALStatus           *status,
    InspiralMetric      *metric,
    InspiralTemplate    *params,
    InspiralMomentsEtc  *moments
    );

int
XLALInspiralComputeMetric (
    InspiralMetric     *metric,
    InspiralMomentsEtc *moments,
    REAL8 fLower,
    LALPNOrder order,
    REAL8 t0,
    REAL8 t3
    );

void
LALInspiralComputeMetricBCV
(
 LALStatus             *status,
 InspiralMetric        *metric,
 REAL8FrequencySeries  *psd,
 InspiralTemplate      *params
);

void
LALInspiralLongestTemplateInBank (
    LALStatus            *status,
    UINT4                *templateLength,
    InspiralCoarseBankIn *coarseIn
    );

void
LALGetInspiralMoments (
    LALStatus            *status,
    InspiralMomentsEtc   *moments,
    REAL8FrequencySeries *psd,
    InspiralTemplate     *params
    );

int
XLALGetInspiralMoments (
    InspiralMomentsEtc   *moments,
    REAL8 fLower,
    REAL8 fCutoff,
    REAL8FrequencySeries *psd
    );

void
LALGetInspiralMomentsBCV (
    LALStatus               *status,
    InspiralMomentsEtcBCV   *moments,
    REAL8FrequencySeries    *psd,
    InspiralTemplate        *params
    );

void
LALInspiralMoments (
    LALStatus         *status,
    REAL8             *moment,
    InspiralMomentsIn pars
    );

REAL8
XLALInspiralMoments(
    REAL8 xmin,
    REAL8 xmax,
    REAL8 ndx,
    REAL8 norm,
    REAL8FrequencySeries *shf
    );

void
LALInspiralMomentsIntegrand
(
   LALStatus *status,
   REAL8  *integrand,
   REAL8  f,
   void   *pars
   );

void
LALInspiralSetSearchLimits (
    LALStatus            *status,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralNextTemplate (
    LALStatus          *status,
    InspiralBankParams *bankPars,
    InspiralMetric      metric
    );

void
LALInspiralComputeParams (
    LALStatus            *status,
    InspiralTemplate     *pars,
    InspiralBankParams   bankParams,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralValidParams (
    LALStatus            *status,
    INT4                 *valid,
    InspiralBankParams   bankParams,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralValidTemplate(
   LALStatus            *status,
   INT4                 *valid,
   InspiralBankParams   bankParams,
   InspiralCoarseBankIn coarseIn
   );

void
LALInspiralUpdateParams (
    LALStatus          *status,
    InspiralBankParams *bankParams,
    InspiralMetric     metric,
    REAL8              minimalMatch
    );

void
LALMatrixTransform (
    LALStatus *status,
    INT4      Dim,
    REAL8     **trans,
    REAL8     **buff1,
    REAL8     **mm3
    );

void
LALDeterminant3 (
    LALStatus *status,
    REAL8  *determinant,
    REAL8  **matrix
    );

void
LALInverse3(
        LALStatus *status,
        REAL8     **inverse,
        REAL8     **matrix
);

void
LALInspiralSetParams (
    LALStatus            *status,
    InspiralTemplate     *tempPars,
    InspiralCoarseBankIn coarseIn
    );

void
LALRectangleVertices(
   LALStatus *status,
   RectangleOut *out,
   RectangleIn *in
   );

void
LALHexagonVertices(
   LALStatus *status,
   HexagonOut *out,
   RectangleIn *in
   );





void
LALInsidePolygon(
   LALStatus            *status,
   REAL4                *inputx,
   REAL4               *inputy,
   INT4                 n,
   REAL4                x,
   REAL4                y,
   INT4                 *valid
   );

void
LALInspiralCreatePNCoarseBankHexa(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void
LALInspiralCreatePNCoarseBankHybridHexa(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    );

void
LALInitHexagonalBank(
    LALStatus         *status,
    InspiralCell      **cell,
    INT4              id,
    InspiralMomentsEtc        *moments,
    InspiralTemplate          *paramsIn,
    HexaGridParam             *gridParam,
    CellEvolution             *cellEvolution,
    CellList                  **cellList
    );


void
LALPopulateCell(
    LALStatus               *status,
    InspiralMomentsEtc      *moments,
    InspiralCell            **cell,
    INT4                    l,
    InspiralTemplate        *paramsIn,
    HexaGridParam           *gridParam,
    CellEvolution           *cellEvolution,
    CellList **cellList
    );

void
LALFindPosition(
    LALStatus               *status,
    REAL4                   dx0,
    REAL4                   dx1,
    Position                *position,
    InspiralTemplate        *paramsIn,
    HexaGridParam           *gridParam
    );


void
LALSPAValidPosition(
    LALStatus           *status,
    InspiralCell        **cell,
    INT4                id1,
    InspiralMomentsEtc  *moments,
    CellEvolution *cellEvolution,
    CellList **list
    );

void
GetPositionRectangle(
    LALStatus *status,
    InspiralCell **cell,
    INT4 id,
    InspiralTemplate *params,
    HexaGridParam           *gridParam,
    CellEvolution *cellEvolution,
    CellList **cellList,
    INT4 *valid
    );



void
LALListAppend(
    CellList ** headRef,
    INT4 id
    );

void
LALListDelete(
    CellList ** headRef,
    INT4 id
    );

UINT4
LALListLength(
    CellList *head
    );

void
LALSPAF(
	LALStatus 	*status,
	REAL4 		*result,
	REAL4 		x,
	void 		*t3
);

INT4 XLALInspiralComputePTFIntrinsicMetric (
    InspiralMetric             *metric,
	REAL8Vector				   *fullmetric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    );

INT4 XLALInspiralComputePTFFullMetric (
    InspiralMetric             *metric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    );

INT4 XLALInspiralComputePTFWaveform (
    REAL8Vector				   *ptfwave,
    InspiralTemplate           *params
    );

INT4 XLALInspiralComputePTFWDeriv (
    COMPLEX16Vector			   *Wderiv,
	REAL8FrequencySeries       *psd,
    InspiralTemplate           *params,
	INT4					   paramid,
	REAL8					   initdelta,
	REAL8					   tolerance
    );

INT4 XLALInspiralComputePTFQDeriv (
    REAL8VectorSequence		   *Qderiv,
    InspiralTemplate           *params
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRALBANK_H */
