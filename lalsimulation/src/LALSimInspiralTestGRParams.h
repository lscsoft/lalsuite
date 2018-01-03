/* Copyright (C) 2012 Walter Del pozzo, Evan Ochsner
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

#ifndef _LALSIMINSPIRALTESTGRPARAMS_H
#define _LALSIMINSPIRALTESTGRPARAMS_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h> 
#include <string.h>
#include <lal/LALMalloc.h>
#include <lal/LALError.h>

#if defined(__cplusplus)
extern "C" {
  #elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @addtogroup LALSimInspiral_h
 * @{
 */

/**
 * Linked list node for the testing GR parameters
 */
typedef struct tagLALSimInspiralTestGRParamData
{
    char name[32]; 	/**< Name of the variable */
    double value;  	/**< Value of the variable */
} LALSimInspiralTestGRParamData;

/**
 * Linked list of any number of parameters for testing GR
 */
typedef struct tagLALSimInspiralTestGRParam
{
    struct tagLALSimInspiralTestGRParamData *data; /**< Current variable */
    struct tagLALSimInspiralTestGRParam *next; /**< The next variable in linked list */
}  LALSimInspiralTestGRParam;

/** @} */

#ifdef SWIG
SWIGLAL(INOUT_STRUCTS(LALSimInspiralTestGRParam**, parameter));
#endif 

LALSimInspiralTestGRParam *XLALSimInspiralCreateTestGRParam(const char *name, double value);
int XLALSimInspiralAddTestGRParam(LALSimInspiralTestGRParam **parameter, const char *name, const double value);
int XLALSimInspiralSetTestGRParam(LALSimInspiralTestGRParam *parameter, const char *name, const double value);
double XLALSimInspiralGetTestGRParam(const LALSimInspiralTestGRParam *parameter, const char *name);
bool XLALSimInspiralTestGRParamExists(const LALSimInspiralTestGRParam *parameter, const char *name);
int XLALSimInspiralPrintTestGRParam(FILE *fp, LALSimInspiralTestGRParam *parameter);
void XLALSimInspiralDestroyTestGRParam(LALSimInspiralTestGRParam *parameter);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALTESTGRPARAMS_H */
