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
#ifndef _LUT_H
#define _LUT_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup LUT_h Header LUT.h
 * \ingroup pkg_pulsarHough
 * \date 2001, 2003
 * \brief Provides structures and function prototypes required for the construction of look up tables
 * that are the core for building the Hough maps.
 * \author Sintes, A.M, Papa, M.A. and Krishnan, B.
 *
 * History:   Created by Sintes May 11, 2001
 * Modified by Badri Krishnan Feb 2003
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/LUT.h>
 * \endcode
 *
 * Our goal is the construction of Hough maps. In order to produce them
 * efficiently, the present implemetation makes use of <tt>lut</tt>s.
 * Here we provide the necessary routines for their construction and use.
 *
 * In principle, the subroutines provided are valid for
 * any Hough master equation of the form:
 * \f[ \nu-F_0  =\vec\xi (t) \cdot (\hat n -\hat N )\, ,\f]
 * where \f$\nu\f$ is the measured frequency of the signal at time \f$t\f$,
 * \f$F_0\f$ intrinsic frequency of the signal at that time, \f$\hat n\f$ location of the
 * souce in the sky, \f$\hat N\f$ the center of the sky patch used in the
 * demodulation procedure,
 * and
 * \f$\vec\xi (t)\f$ any vector.
 *
 * The  form of this vector \f$\vec\xi (t)\f$
 * depends on the demodulation procedure
 * used in the previous  step.  In our case this corresponds to
 * \f[\vec\xi (t) = \left( F_0+
 * \sum_k F_k \left[ \Delta T  \right]^k\right) \frac{\vec v(t)}{c}
 * + \left( \sum_k k F_k \left[ \Delta  T \right]^{k-1}\right)
 * \frac {\vec x(t)- \vec x(\hat t_0)}{c}\, ,\f]
 * and
 * \f[F_0 \equiv  f_0 + \sum_k \Delta f_k
 * \left[ \Delta T \right]^k \, , \f]
 * where
 * \f$\vec v(t)\f$ is the velocity of the detector, \f$\vec x(t)\f$ is the detector
 * position,
 * \f$ T_{\hat N}(t)\f$ is the  time  at
 * the solar system barycenter (for a given sky
 * location \f$\hat N\f$),
 * \f$\Delta T \equiv T_{\hat N}(t)-T_{\hat N}(\hat t_0)\f$,
 * \f$\Delta f_k = f_k-F_k\f$ the  residual spin-down parameter,
 * \f$F_k\f$ the  spin-down parameter used in the demodulation, and \f$f_0\f$, \f$f_k\f$ the
 * intrinsic  frequency and spin-down parameters of the source at time \f$\hat t_0\f$.
 *
 * Looking at the generic Hough master equation, one realizes that
 * for a fixed  time, a given value of \f$F_0\f$, and a measured frequency \f$\nu\f$
 * p(from a selected peak), the source could be located anywhere on a circle (whose
 * center points in the same direction of \f$\vec\xi (t)\f$ and is characterized by
 * \f$\phi\f$, the angle between \f$\hat n\f$ and \f$\vec\xi\f$).  Since the Hough transform is performed on a
 * set of spectra with discrete frequencies, a peak on the spectrum appearing at
 * \f$\nu\f$  could correspond to any source with a demodulated
 * frequency in a certain interval. As a consequence, the location of the sources
 * compatible with \f$F_0\f$ and \f$\nu\f$ is not a circle but an annulus with a certain
 * width.
 *
 * Our purpose is to map these annuli on a discrete space. An estimation of the
 * average thickness of the annuli  tells us that the vast majority of
 * annuli will be very thin, and therefore our algorithm should not be
 * optimized for drawing
 * thick annuli but for thin ones. Also, the mapping implementation should be one with a uniform
 * probability distribution in order to avoid discretization errors.
 * In order to remove border effects, we use a biunivocal mapping, which
 * requires
 * that a pixel in a partial Hough map can belong only to one annulus,
 * just touched by one peak of the spectrum. The criteria for the biunivocal mapping
 * is that if and only if the center of the pixel is inside the annulus, then the pixel
 * will be enhanced.
 *
 * In order to simplify  (reduce the computational cost of) this task we
 * construct look up tables (<tt>lut</tt>) where the  borders of these annuli are
 * marked for any possible \f$\nu -F_0\f$ value. Since we work on a discrete space these <tt>lut</tt> are
 * valid for many \f$F_0\f$ values.
 *
 * \image html LUTstereo.png "Fig. [fig_stereo]: Stereographic projection. [Credit to: D.Hilbert, S.Cohn-Vossen, P. Nemenyi, ``Geometry and Imagination'', Chelsea Publishing Company, New York 1952.]"
 * \image latex LUTstereo.eps "Stereographic projection. [Credit to: D.Hilbert, S.Cohn-Vossen, P. Nemenyi, ``Geometry and Imagination'', Chelsea Publishing Company, New York 1952.]"
 *
 * At this point we have already chosen a sky tiling to produce the Hough map  efficiently.
 * It consists of changing coordinates so that the center of the   patch is located at
 * \f$(0,-\pi/2)\f$ in (\f$\alpha-\delta\f$) (or in any other coordinate system), then we make use
 * of the stereographic projection (see Fig.\figref{fig_stereo}) and we take horizontal and vertical lines on
 * the projected plane at constant space separation.
 * This projection is advantageous because it avoids distortions, i.e. the
 * pixel size is almost constant independently of the sky location, and makes the
 * algorithm simpler. The stereographic projection has the property to map circles
 * on the celestial sphere to circles on the projected plane.
 *
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
# include <stdlib.h>
# include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/**\name Error Codes */
/*@{*/
#define LUTH_ENULL 1
#define LUTH_ESIZE 2
#define LUTH_ESZMM 4
#define LUTH_EINT  6
#define LUTH_ESAME 8
#define LUTH_EFREQ 10
#define LUTH_EVAL 12

#define LUTH_MSGENULL "Null pointer"
#define LUTH_MSGESIZE "Invalid input size"
#define LUTH_MSGESZMM "Size mismatch"
#define LUTH_MSGEINT  "Invalid interval"
#define LUTH_MSGESAME "Input/Output data vectors are the same"
#define LUTH_MSGEFREQ "Invalid frequency"
#define LUTH_MSGEVAL  "Invalid value"
/*@}*/


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


#define MAX(A, B)  (((A) < (B)) ? (B) : (A))
#define MIN(A, B)  (((A) < (B)) ? (A) : (B))
#define cot(A)  (1./tan(A))

/** \name Constant declarations */
/*@{*/

/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */

/**
 * Maximum ``error'' (as a fraction of the width of the thinnest annulus)
 * which allows to represent a circle by  a line.
 */
#define LINERR     0.001

/**
 * Maximum ``error'' (as a fraction of the width of the thinnest annulus)
 * which allows to consider two border equivalents!
 * It is relevant for determining the LUT frequency range validity
 */
#define PIXERR     0.5


/* Width of the thinnest annulus in terms of pixels */
/* #define PIXELFACTOR  2 */

/** Earth v_epicycle/c  TO BE CHANGED DEPENDING ON DETECTOR */
#define VEPI LAL_TWOPI * LAL_REARTH_SI / ( LAL_DAYSID_SI * LAL_C_SI )
/* #define VEPI 1.0e-06 */
/* vepi = 2*pi*R_earth /(day *c) */

/** Total detector velocity/c TO BE CHANGED DEPENDING ON DETECTOR */
#define VTOT LAL_TWOPI * LAL_AU_SI / ( LAL_YRSID_SI * LAL_C_SI )
/* #define VTOT 1.06e-04 */
/* vtot = 2*pi* 1AU / (year * c) */

/*@}*/

/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

/** To be changed to {<tt>INT2 COORType</tt>} if the number of pixels in the x-direction exceeds 255. */
typedef INT2 COORType;  /* typedef INT4 COORType; */ /* typedef  UCHAR COORType; */

/** This structure stores the border of a circle clipped on the projected plane. */
typedef struct tagHOUGHBorder{
  INT4  yUpper;		/**< upper y pixel affected by this border */
  INT4  yLower;    	/**< lower y pixel affected by this border and yUpper<yLower or yUpper<0 are possible */
  INT4  yCenter;   	/**< y pixel value of the center of the circle */
  UINT2     ySide; 	/**< length of xPixel */
  COORType *xPixel; 	/**< x pixel index to be marked */
} HOUGHBorder;


/**
 * This structure stores the border indexes corresponding to one frequency
 * bin plus the corrections to be added to the first column of the patch.
 */
typedef struct tagHOUGHBin2Border{
  INT2   leftB1;     	/**< Border index to be used (<em>start-border</em> `\f$+1\f$') */
  INT2   rightB1;    	/**< Border index to be used (<em>stop-border</em> `\f$-1\f$') */
  INT2   leftB2;     	/**< Border index  to be used (<em>start-border</em> `\f$+1\f$') */
  INT2   rightB2;	/**< Border index  to be used (<em>stop-border</em> `\f$-1\f$') */
  INT2   piece1max;  	/**< Interval limits of the (first piece) correction to the first column. */
  INT2   piece1min;	/**< If <tt>piece1min \f$>\f$ piece1max</tt> no corrections should be added */
  INT2   piece2max;	/**< Interval limits of the (second piece)  correction to the first column */
  INT2   piece2min;	/**< If <tt>piece2min \f$>\f$ piece2max</tt> no corrections should be added */
} HOUGHBin2Border;

/** This structure stores the patch-time-frequency <em>look up table</em>. */
typedef struct tagHOUGHptfLUT{
  INT2    timeIndex;  /**< time index of the \c lut */
  INT8    f0Bin;      /**< Frequency bin for which it has been constructed. */
  REAL8   deltaF;     /**< Frequency resolution <tt>df=1/TCOH</tt>, where <tt>1/TCOH</tt>
                       * is the coherent integration time used in teh demodulation procedure.
                       */
  INT8    nFreqValid; 	/**< Number of frequencies where the \c lut is valid */
  INT4    iniBin;     	/**< First bin affecting the patch with respect to \c f0 */
  INT4    nBin;       	/**< Exact number of bins affecting the patch */
  INT4    offset;	/**< Frequency bin corresponding to center of patch measured with respect to \c f0Bin (zero in modulated case) */
  UINT2   maxNBins;    /**< Maximum number of bins affecting the patch (for memory allocation purposes) */
  UINT2   maxNBorders; /**< Maximum number of borders affecting the patch (for memory allocation purposes) */
  HOUGHBorder      *border; /**< The annulus borders */
  HOUGHBin2Border  *bin;    /**< Bin to border correspondence. */
} HOUGHptfLUT;


/** This structure stores patch-frequency \e grid information. */
typedef struct tagHOUGHPatchGrid{
  REAL8   f0;         	/**< Frequency to construct grid */
  REAL8   deltaF;     	/**< Frequency resolution: <tt>df=1/TCOH</tt> */
  REAL8   deltaX;	/**< Longitudinal space resolution, x-direction */
  REAL8   xMin;     	/**< Patch limit, as the coordinate of the center of the first pixel */
  REAL8   xMax;		/**< Patch limit, as the coordinate of the center  of the last pixel. */
  UINT2   xSide;    	/**< Real number of pixels in the x direction (in the projected plane); it should be less than or equal to \c xSideMax. */
  REAL8   *xCoor;   	/**< Coordinates of the pixel centers */
  REAL8   deltaY;	/**< Longitudinal space resolution, y-direction */
  REAL8   yMin;     	/**< Patch limit, as center of the first pixel */
  REAL8   yMax;		/**< Patch limit, as center of the last pixel */
  UINT2   ySide;    	/**< Real number of pixels in the y-direction (in the projected
                         * plane). It should be less than or equal to \c ySideMax.
                         */
  REAL8   *yCoor;   	/**< Coordinates of the pixel centers. */
} HOUGHPatchGrid;

/** parameters needed for gridding the patch */
typedef struct tagHOUGHResolutionPar{
  INT8    f0Bin; 	/**< Frequency bin at which construct the grid */
  REAL8   deltaF;	/**< Frequency resolution: <tt>df=1/TCOH</tt> */
  REAL8   patchSkySizeX;/**< Patch size in radians along x-axis */
  REAL8   patchSkySizeY;/**< Patch size in radians along y-axis */
  REAL8   pixelFactor; /**< number of pixel that fit in the thinnest annulus*/
  REAL8   pixErr;   	/**< for validity of LUT as PIXERR */
  REAL8   linErr;   	/**< as LINERR circle ->line */
  REAL8   vTotC;    	/**< estimate value of v-total/C as VTOT */
} HOUGHResolutionPar;

/** required for constructing patch */
typedef struct tagHOUGHSizePar{
  INT8    f0Bin; 	/**< corresponding freq bin  */
  REAL8   deltaF;       /**< df=1/TCOH */
  REAL8   deltaX; 	/**< pixel size in the projected plane */
  REAL8   deltaY;
  UINT2   xSide;    	/**< number of pixels in the x direction (projected plane)*/
  UINT2   ySide;    	/**< number of pixels in the y direction */
  UINT2   maxNBins;    	/**< maximum number of bins affecting the patch; for memory allocation */
  UINT2   maxNBorders; 	/**< maximum number of borders affecting the patch; for memory allocation */
  INT8    nFreqValid; 	/**< number of frequencies where the LUT is valid */
  REAL8   epsilon; 	/**< max. angle (rad.) from the pole to consider a circle as a line in the projected plane */
} HOUGHSizePar;

/** Three dimensional Cartessian coordinates */
typedef struct tagREAL8Cart3Coor{
  REAL8  x;
  REAL8  y;
  REAL8  z;
} REAL8Cart3Coor;

/** Two dimensional Cartessian coordinates */
typedef struct tagREAL8Cart2Coor{
  REAL8  x;
  REAL8  y;
} REAL8Cart2Coor;

/** Two dimensional polar coordinates. */
typedef struct tagREAL8Polar2Coor{
  REAL8  alpha;
  REAL8  radius;
} REAL8Polar2Coor;

/** Polar coordinates of a unitary vector on the sphere */
typedef struct tagREAL8UnitPolarCoor{
  REAL8  alpha;  	/**< any value */
  REAL8  delta;  	/**< In the interval [\f$-\pi/2, \,  \pi/2\f$] */
} REAL8UnitPolarCoor;

/** Parameters needed to construct the partial look up table */
typedef struct tagHOUGHParamPLUT{
  INT8             f0Bin;   	/**< freq. bin for which it has been constructed */
  REAL8            deltaF;  	/**< Frequency resolution: <tt>df=1/TCOH</tt> */
  REAL8UnitPolarCoor xi;  	/**< Center of the circle on the celestial sphere,
                                 * \f$\xi\f$(alpha,delta) in the rotated coordinates
                                 */
  REAL8            cosDelta;    /**< \f$\Delta \cos(\phi)\f$ corresponding to one annulus. */
  INT4             offset;	/**< Frequency bin corresponding to center of patch; measured w.r.t. \c f0Bin */
  INT8             nFreqValid;	/**< Number of frequency bins for which the LUT is valid */
  REAL8            cosPhiMax0;	/**< \f$\max(\cos(\phi))\f$ of the \c f0Bin */
  REAL8            cosPhiMin0;	/**< \f$\min(\cos(\phi))\f$ of the \c f0Bin */
  REAL8            epsilon; 	/**< maximum angle (distance in radians) from the pole
                                 * to consider  a circle as a line in the projected plane
                                 */
} HOUGHParamPLUT;

/**
 * Demodulation parameters needed for the Hough transform; all
 * coordinates are assumed to be with respect to the same reference system.
 */
typedef struct tagHOUGHDemodPar{
  REAL8               deltaF;   /**< Frequency resolution: <tt>df=1/TCOH</tt> */
  REAL8UnitPolarCoor  skyPatch; /**< \f$N_{center}\f$ (alpha, delta): position of the center of the patch */
  REAL8Cart3Coor      veloC;    /**< \f$v(t)/c\f$ (x,y,z): Relative detector velocity */
  REAL8Cart3Coor      positC;   /**< \f$(x(t)-x(t0))/c\f$ (x,y,z): Position of the detector */
  REAL8               timeDiff; /**< \f$T_{\hat N} (t)-T_{\hat N} (\hat t_0)\f$: time difference */
  REAL8Vector         spin; 	/**< Spin down information. It includes the fields:
                                 * \c length: maximum order of spin-down parameter, and
                                 * <tt>*data</tt>: pointer to spin-down parameter set \f$F_k\f$
                                 */
} HOUGHDemodPar;

/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */
void LALHOUGHComputeSizePar (LALStatus  *status, /* demod case */
                   HOUGHSizePar        *out,
                   HOUGHResolutionPar  *in
		   );

void LALHOUGHComputeNDSizePar (LALStatus  *status, /* non -demod case*/
                   HOUGHSizePar          *out,
                   HOUGHResolutionPar    *in
		   );


void LALHOUGHFillPatchGrid (LALStatus   *status,
		   HOUGHPatchGrid      *out,
                   HOUGHSizePar        *par
		       );

void LALHOUGHCalcParamPLUT (LALStatus   *status, /* Demod. case */
                   HOUGHParamPLUT   *out,  /* parameters needed build LUT*/
		   HOUGHSizePar     *sizePar,
                   HOUGHDemodPar    *par  /* demodulation parameters */
			);

void LALNDHOUGHParamPLUT (LALStatus   *status, /* non-demod. case */
                   HOUGHParamPLUT   *out,  /* parameters needed build LUT*/
		   HOUGHSizePar     *sizePar,
                   HOUGHDemodPar    *par  /* demodulation parameters */
			);

void LALRotatePolarU(LALStatus            *status,
		     REAL8UnitPolarCoor   *out,
		     REAL8UnitPolarCoor   *in,
		     REAL8UnitPolarCoor   *par
		     );

void LALInvRotatePolarU(LALStatus            *status,
			REAL8UnitPolarCoor   *out,
			REAL8UnitPolarCoor   *in,
			REAL8UnitPolarCoor   *par
			);

void LALStereoProjectPolar(LALStatus           *status,
			   REAL8Polar2Coor     *out,
			   REAL8UnitPolarCoor  *in
			   );

void LALStereoProjectCart(LALStatus           *status,
			  REAL8Cart2Coor      *out,
			  REAL8UnitPolarCoor  *in
			  );

void LALStereoInvProjectPolar(LALStatus           *status,
			      REAL8UnitPolarCoor  *out,
			      REAL8Polar2Coor     *in
			      );

void LALStereoInvProjectCart(LALStatus           *status,
			     REAL8UnitPolarCoor  *out,
			     REAL8Cart2Coor      *in
			     );

void LALHOUGHConstructPLUT(LALStatus       *status,
			   HOUGHptfLUT     *lut,
			   HOUGHPatchGrid  *patch,
			   HOUGHParamPLUT  *par
			   );

/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _LUT_H */













