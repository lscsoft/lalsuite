/*
*  Copyright (C) 2013 Chris Messenger
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

#ifndef _SEMICOHERENT_H
#define _SEMICOHERENT_H

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDemod.h>
#include <lal/Date.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif
 
#define NPARAMS 9
#define NPARAMS_IL 3
  
 /** A single parameter dimensions boundaries
 */
typedef struct {
  REAL8 min;                        /**< the parameter space minimum */
  REAL8 max;                        /**< the parameter space maximium */
  REAL8 span;                       /**< the parameter space span */
  CHAR name[LALNameLength];         /**< string containing the name of the dimension */
} REAL8Dimension;

/** A vector of parameter space boundary information
 */
typedef struct {
  REAL8Dimension *data;             /**< the boundaries, span, etc for a single dimension */
  UINT4 ndim;                       /**< the number of dimensions */
} REAL8Space;

/** Stores the gridding parameters for a single dimension
 */
typedef struct {
  REAL8 min;                        /**< the starting points of the grid */
  REAL8 delta;                      /**< the grid spacings */
  REAL8 oneoverdelta;               /**< the inverse of the spacing */
  UINT4 length;                     /**< the number of templates in each dimension */
  CHAR name[LALNameLength];         /**< string containing the name of the dimension */
} Grid;

/** Stores the current location in a hyper-cubic parameter space
 */
typedef struct {
  REAL8 *x;                         /**< the location in parameter space */
  INT4 *idx;                        /**< the index of each dimension for this template */
  UINT4 ndim;                       /**< the dimension of the parameter space */
  UINT4 currentidx;                 /**< the current index value of the template */
} Template;

/** Stores the gridding parameters for a hypercubic grid of templates
 */
typedef struct {
  Grid *grid;                       /**< stores the parameters defining a single dimension */
  UINT4 ndim;                       /**< the number of dimensions */
  UINT4 *prod;                      /**< internal variable used to store the size of sub-dimensions */
  UINT4 max;                        /**< the maximum (total) number of templates */
  REAL8 mismatch;                   /**< the mismatch */
  INT4 Nr;
} GridParameters;

  /* LALDemod functions now put into CFSLALDemod.c */
  void GenerateSideBandTemplate(LALStatus *,BinarySourceParams *,SideBandTemplateParams *,SideBandTemplate **);
  void ComputeABCcoefficients (LALStatus *,ABCcoParams *,LALDetector *,ABCcoefficients **);
  void ReadTimeStamps(LALStatus *,CHAR *,INT4, SideBandTemplateParams **);
  void ComputeSideBandWindow(LALStatus *,ABCcoefficients *, CHAR *,SideBandTemplateParams **);
  void ComputeSideBandLikelihood(LALStatus *,SideBandMCMCVector *,SideBandDataSet *,SideBandTemplate **,SideBandTemplateParams *);
  void InitEphemeris(LALStatus *,EphemerisData *,const CHAR *,const CHAR *);
  void ReadSideBandPriors(LALStatus *status,CHAR *,SideBandMCMCRanges *,SideBandMCMCJumpProbs *);
  void SelectSideBandFrequencies (LALStatus *,SideBandDataSet **,SelectSideBandFrequencyParams *,SideBandDataSet **);
  void ReadSideBandData (LALStatus *,ReadSideBandDataParams *,SideBandDataSet **);
  void EstimateSideBandNoise (LALStatus *, SideBandDataSet **,EstimateSideBandNoiseParams *);
  REAL4 bessj(INT4, REAL4);
  REAL4 bessj1(REAL4);
  REAL4 bessj0(REAL4);
      

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
