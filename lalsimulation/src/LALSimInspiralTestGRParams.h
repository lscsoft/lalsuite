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

/**
 * Function that creates the head node of the test GR parameters linked list.
 * It is initialized with a single parameter with given name and value
 */
 
#ifdef SWIG   // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(LALSimInspiralTestGRParam**, parameter));
#endif 

LALSimInspiralTestGRParam *XLALSimInspiralCreateTestGRParam(
        const char *name, /**< Name of first parameter in new linked list */
        double value 	 /**< Value of first parameter in new linked list */
        );

/**
 * Function that adds a prameter to the test GR parameters linked list. If the
 * parameter already exists, it throws an error.
 */
int XLALSimInspiralAddTestGRParam(
        LALSimInspiralTestGRParam **parameter, /**< Pointer to the head node of the linked list of parameters */
        const char *name, 			/**< Parameter name */
        const double value 			/**< Parameter value */
        );

/**
 * Function that sets the value of the desired parameter in the test GR
 * parameters linked list to 'value'.  Throws an error if the parameter
 * is not found
 */
int XLALSimInspiralSetTestGRParam(
        LALSimInspiralTestGRParam *parameter, /**< Linked list to be modified */
        const char *name, 		/**< Name of parameter to be modified */
        const double value 		/**< New value for parameter */
        );

/**
 * Function that returns the value of the desired parameters in the
 * test GR parameters linked list.  Aborts if the parameter is not found
 */
double XLALSimInspiralGetTestGRParam(
        const LALSimInspiralTestGRParam *parameter, /**< Linked list to retrieve from */
        const char *name 	/**< Name of parameter to be retrieved */
        );

/**
 * Function that checks whether the requested parameter exists within the
 * test GR parameters linked list.  Returns true (1) or false (0) accordingly
 */
bool XLALSimInspiralTestGRParamExists(
        const LALSimInspiralTestGRParam *parameter, /**< Linked list to check */
        const char *name 		/**< Parameter name to check for */
        );

/** Function that prints the whole test GR params linked list */
int XLALSimInspiralPrintTestGRParam(
        FILE *fp, 			/** FILE pointer to write to */
        LALSimInspiralTestGRParam *parameter /**< Linked list to print */
        );

/** Function that destroys the whole test GR params linked list */
void XLALSimInspiralDestroyTestGRParam(
        LALSimInspiralTestGRParam *parameter /**< Linked list to destroy */
        );

#endif /* _LALSIMINSPIRALTESTGRPARAMS_H */
