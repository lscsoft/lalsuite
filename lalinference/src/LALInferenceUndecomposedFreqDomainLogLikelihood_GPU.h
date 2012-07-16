/*
*  Copyright (C) 2011 Alex Ayerdi
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

/*-----------------------------------------------------------------------
 *
 * File Name: LALInferenceUndecomposedFreqDomainLogLikelihood_GPU.h
 *
 * Author: Ayerdi, A.
 *
 *
 *-----------------------------------------------------------------------
 */

#ifndef _LALInferenceUndecomposedFreqDomainLogLikelihood_GPU_H
#define _LALInferenceUndecomposedFreqDomainLogLikelihood_GPU_H

#define LAL_USE_OLD_COMPLEX_STRUCTS
#pragma GCC system_header
#define INT64_C
#define __STDC_CONSTANT_MACROS

#ifdef __cplusplus
extern "C" {
#include <lal/LALInferenceLikelihood.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
}
#endif

//REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood_GPU(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
  //                            LALInferenceTemplateFunction *_template);

#endif /* _LALInferenceUndecomposedFreqDomainLogLikelihood_GPU_H */
