/*
*  Copyright (C) 2007 Bernd Machenschalk, Chad Hanna, Jolien Creighton, Benjamin Owen
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

/*________________________________________________________________________
 *
 * File Name: TemplateBankGeneration.h
 *
 * Author: Hanna C. R.
 *
 *________________________________________________________________________
 */

#ifndef _TEMPLATEBANKGENERATION_H
#define _TEMPLATEBANKGENERATION_H

#include<lal/LALStdlib.h>
#include<lal/LALStatusMacros.h>
#if 0
#include<lal/LALInspiral.h>
#include<lal/LALInspiralBank.h>
#endif
#include<lal/LALDatatypes.h>
#include<lal/LIGOMetadataTables.h>

/**
 * \addtogroup TemplateBankGeneration_h
 * \author Hanna, C. R.
 *
 * \brief This header file includes all the necessary types and
 * function prototypes for LALNDTemplateBank() and LALMakeTemplateBank().
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/TemplateBankGeneration.h>
 * \endcode
 *
 * NDTemplateBank() provides a general way to tile a parameter space with
 * a constant metric ( currently only for less than 12 dimensions).
 * MakeTemplateBank() provides a general way for applications to
 * generate a template bank with suitable I/O.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define TEMPLATEBANKGENERATIONH_ENULL 1
#define TEMPLATEBANKGENERATIONH_MSGENULL "Unexpected NULL pointer to an input type"
/*@}*/

typedef enum {
  /* Binary Inspiral Searches 100-199 */
     BCVType,
     BCVSpinType,
     PrecessingType,
  /* Pulsar Searches 200-299 */
     Pulsar,
  /* Burst Searches 300-399 */
     Burst
  /* Other Searches 400-499 Etc... */
     } TemplateBankType;

typedef struct tagNDTemplateBankInput {
  REAL4                	minCoordinates[12]; 	/**< Psi0, Psi3, etc      */
  REAL4                	maxCoordinates[12]; 	/**< Psi0, Psi3, etc      */
  REAL4                	minParameters[12];  	/**< mass, eta, etc       */
  REAL4                	maxParameters[12];  	/**< mass, eta, etc       */
  INT2                 	dimension;          	/**< 3D?, 4D? -> ND!       */
  REAL4                	mm;                 	/**< mismatch              */
  TemplateBankType     	type;               	/**< whats the search?     */
  REAL8FrequencySeries *PSD;		   	/**< Power Spec Density    */
  REAL4			f0;			/**< Moment scaling freq   */
  } NDTemplateBankInput;



typedef struct tagNDTemplateBankOutput{
  REAL4				coordinateVals[12];
  REAL4 			parameterVals[12];
  struct tagNDTemplateBankOutput	*next;
  } NDTemplateBankOutput;



typedef void (*NDTemplateBankMetricPtr)( LALStatus *, NDTemplateBankInput *, REAL4Array *);
typedef void (*NDTemplateBankTestPtr)(LALStatus *, NDTemplateBankInput *, NDTemplateBankOutput *, INT2 *);

typedef struct tagNDTemplateBankFunctionPtrs {
  NDTemplateBankMetricPtr 	metric; 	/**< Ptr to metric function */
  NDTemplateBankTestPtr		test;		/**< Ptr boundary test fct  */
  } NDTemplateBankFunctionPtrs;


void
LALNDTemplateBank(
 	LALStatus *,
   	NDTemplateBankInput *,
        NDTemplateBankFunctionPtrs *,
       	NDTemplateBankOutput **);


/*@}*/ /* end:TemplateBankGeneration_h */

#if 0
void
LALMakeTemplateBank(
     	LALStatus *,
     	TemplateBankType *,
     	MakeTemplateBankInput *,
     	MetadataTable *);
     /* LALMakeTemplateBank(); */
#endif

#endif
