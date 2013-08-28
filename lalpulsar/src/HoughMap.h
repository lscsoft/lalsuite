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
#ifndef _HOUGHMAP_H
#define _HOUGHMAP_H

#ifdef  __cplusplus
extern "C" {
#endif

/*  * History:   Created by Sintes June 22, 2001
 *            Modified    August 6, 2001
 */

/**
 * \defgroup HoughMap_h Header HoughMap.h
 * \ingroup pkg_pulsarHough
 * \author Alicia M. Sintes and Badri Krishnan
 *
 * \brief Provides subroutines for initialization and construction of Hough-map derivatives and total Hough-maps.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/HoughMap.h>
 * \endcode
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

# include <lal/LUT.h>
# include <lal/PHMD.h>


/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/**\name Error Codes */
/*@{*/
#define HOUGHMAPH_ENULL 1
#define HOUGHMAPH_ESIZE 2
#define HOUGHMAPH_ESZMM 4
#define HOUGHMAPH_EINT  6
#define HOUGHMAPH_ESAME 8
#define HOUGHMAPH_EFREQ 10
#define HOUGHMAPH_EVAL 12

#define HOUGHMAPH_MSGENULL "Null pointer"
#define HOUGHMAPH_MSGESIZE "Invalid input size"
#define HOUGHMAPH_MSGESZMM "Size mismatch"
#define HOUGHMAPH_MSGEINT  "Invalid interval"
#define HOUGHMAPH_MSGESAME "Input/Output data vectors are the same"
#define HOUGHMAPH_MSGEFREQ "Invalid frequency"
#define HOUGHMAPH_MSGEVAL  "Invalid value"
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
 * Total Hough Map pixel type.
 * Depending of the number of maps to accumulate
 * change both types \c HoughDT and \c HoughTT to \c INT2 or \c UINT2 respectively.
 */
  typedef REAL8 HoughTT; /* for weighted hough maps only */
  /* Depending of the number of maps to accumulate, */
  /* if needed change both types  to INT2 or UINT2  */
/* typedef UCHAR HoughTT; */
/*typedef UINT2 HoughTT; */


/** This structure stores the Hough map derivative */
typedef struct tagHOUGHMapDeriv{
  UINT2     xSide;  /**< number of physical pixels in the x direction */
  UINT2     ySide;  /**< number of physical pixels in the y direction */
  HoughDT   *map ;  /**< the pixel count derivatives;
                     * the number of elements to allocate is ySide*(xSide+1)* */
} HOUGHMapDeriv;


/**  This structure stores the Hough map */
typedef struct tagHOUGHMapTotal{
  INT8               f0Bin;      /**< frequency bin for which it has been constructed */
  REAL8              deltaF;     /**< frequency resolution */
  UINT4              mObsCoh;    /**< ratio of observation time and coherent timescale */
  UINT4              nPG;        /**< number of peakgrams used  <tt><= mObsCoh</tt>; there could be gaps during the observation time */
  REAL8              patchSizeX; /**< x size of patch */
  REAL8              patchSizeY; /**< y size of patch */
  REAL8UnitPolarCoor skyPatch;   /**< Coordinates of the versor \f$\hat N_{center}\f$ (alpha, delta) pointing to the center of the sky patch */
  REAL8Vector        spinDem;    /**< Spin parameters used in the demodulation stage */
  REAL8Vector        spinRes;    /**< Refined spin parameters used in the Hough transform */
  REAL8Vector        dFdot;      /**< resolution in spindown parameters */
  UINT2              xSide;      /**< number of physical pixels in the x direction */
  UINT2              ySide;      /**< number of physical pixels in the y direction */
  HoughTT            *map;       /**< the pixel counts; the number of elements to allocate is ySide*xSide */
} HOUGHMapTotal;

/*
 * 11. Extern Global variables. (discouraged)
 */

/*
 * 12. Functions Declarations (i.e., prototypes).
 */

void LALHOUGHInitializeHD (LALStatus      *status,
			  HOUGHMapDeriv   *hd /* the Hough map derivative */
			  );

void LALHOUGHAddPHMD2HD (LALStatus      *status,
			 HOUGHMapDeriv  *hd,  /* the Hough map derivative */
			 HOUGHphmd      *phmd  /* info from a partial map */
			 );

void LALHOUGHAddPHMD2HD_W (LALStatus      *status,
			   HOUGHMapDeriv  *hd,  /* the Hough map derivative */
			   HOUGHphmd      *phmd  /* info from a partial map */
			   );

void LALHOUGHIntegrHD2HT (LALStatus       *status,
			  HOUGHMapTotal   *ht,     /* the total Hough map */
			  HOUGHMapDeriv   *hd /* the Hough map derivative */
			  );

void LALHOUGHInitializeHT (LALStatus      *status,
			  HOUGHMapTotal   *ht,     /* the total Hough map */
			  HOUGHPatchGrid  *patch      /* patch information */
			  );

void LALStereo2SkyLocation (LALStatus  *status,
			    REAL8UnitPolarCoor *sourceLocation, /* output*/
			    UINT2              xPos,
			    UINT2              yPos,
			    HOUGHPatchGrid    *patch,
			    HOUGHDemodPar     *parDem);

/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _HOUGHMAP_H */
