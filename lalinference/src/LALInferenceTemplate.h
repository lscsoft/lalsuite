/*
 *
 *  LALInference:             Bayesian Followup        
 *  LALInferenceTemplate.h:   Template generation functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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
 * \file LALInference.h
 * \brief Main header file for LALInference common routines
 */

#ifndef LALInferenceTemplate_h
#define LALInferenceTemplate_h

#include <lal/LALInference.h>

void dumptemplateFreqDomain(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                            LALInferenceTemplateFunction *template, const char *filename);
void dumptemplateTimeDomain(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                            LALInferenceTemplateFunction *template, const char *filename);

void LALTemplateGeneratePPN(LALInferenceIFOData *IFOdata);
void templateStatPhase(LALInferenceIFOData *IFOdata);
void templateNullFreqdomain(LALInferenceIFOData *IFOdata);
void templateNullTimedomain(LALInferenceIFOData *IFOdata);
void templateLAL(LALInferenceIFOData *IFOdata);
void template3525TD(LALInferenceIFOData *IFOdata);
void templateSineGaussian(LALInferenceIFOData *IFOdata);
void templateDampedSinusoid(LALInferenceIFOData *IFOdata);
void templateSinc(LALInferenceIFOData *IFOdata);
void templateLALSTPN(LALInferenceIFOData *IFOdata);
void templateASinOmegaT(LALInferenceIFOData *IFOdata);
void templateLALGenerateInspiral(LALInferenceIFOData *IFOdata);

#endif
