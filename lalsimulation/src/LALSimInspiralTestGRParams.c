/* Copyright (C) 2012 Walter Del Pozzo, Evan Ochsner and Salvatore Vitale
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
 
#include  <lal/LALSimInspiralTestGRParams.h>

/**
 * Function that creates the head node of the test GR parameters linked list.
 * It is initialized with a single parameter with given name and value
 */
LALSimInspiralTestGRParam *XLALSimInspiralCreateTestGRParam(
        const char *name, /**< Name of first parameter in new linked list */
        double value 	 /**< Value of first parameter in new linked list */
        )
{
        LALSimInspiralTestGRParam *parameter = (LALSimInspiralTestGRParam *)XLALMalloc(sizeof(LALSimInspiralTestGRParam));
        if (parameter) 
        {
            parameter->data =  (LALSimInspiralTestGRParamData *)XLALMalloc(sizeof(LALSimInspiralTestGRParamData));
            memcpy(parameter->data->name, name, 32);
            parameter->data->value = value;
        }
        parameter->next=NULL;
        return parameter;
}

/**
 * Function that adds a prameter to the test GR parameters linked list. If the
 * parameter already exists, it throws an error.
 */
int XLALSimInspiralAddTestGRParam(
        LALSimInspiralTestGRParam **parameter, /**< Pointer to the head node of the linked list of parameters */
        const char *name, 		/**< Parameter name */
        double value 			/**< Parameter value */
        )
{
    LALSimInspiralTestGRParam *temp;
    temp = *parameter;
    if (*parameter==NULL) 
    {
        temp = XLALSimInspiralCreateTestGRParam(name,value); 
        //temp->next=NULL;
        *parameter=temp;
    }
    else 
    {

        if (!XLALSimInspiralTestGRParamExists(*parameter, name))
        {
            temp = *parameter;
             while(temp->next!=NULL) {temp=temp->next;}
            LALSimInspiralTestGRParam *newParam = XLALSimInspiralCreateTestGRParam(name,value);        
            temp->next = newParam;
        }
        else 
        {
            XLALPrintError("XLAL Error - %s: parameter '%s' exists already! Not added to the structure\n",
                    __func__, name);
            XLAL_ERROR(XLAL_EINVAL);
        }
    }
    return XLAL_SUCCESS;
}

/**
 * Function that sets the value of the desired parameter in the test GR
 * parameters linked list to 'value'.  Throws an error if the parameter
 * is not found
 */
int XLALSimInspiralSetTestGRParam(
        LALSimInspiralTestGRParam *parameter, /**< Linked list to be modified */
        const char *name, 		/**< Name of parameter to be modified */
        const double value 		/**< New value for parameter */
        )
{
    if (XLALSimInspiralTestGRParamExists(parameter, name)) 
    {
        while(parameter)
        {
            if(!strcmp(parameter->data->name, name)) parameter->data->value = value;
            parameter=parameter->next;
        }
        return XLAL_SUCCESS;
    }
    else
    {
        XLALPrintError("XLAL Error - %s: parameter '%s' unknown!\n",
                __func__, name);
        XLAL_ERROR(XLAL_EINVAL);
    }
}

/**
 * Function that returns the value of the desired parameters in the
 * test GR parameters linked list.  Aborts if the parameter is not found
 */
double XLALSimInspiralGetTestGRParam(
        const LALSimInspiralTestGRParam *parameter, /**< Linked list to retrieve from */
        const char *name 	   /**< Name of parameter to be retrieved */
        )
{
    if (XLALSimInspiralTestGRParamExists(parameter, name)) 
        {
            while(parameter) 
            {
                if(!strcmp(parameter->data->name, name)) return parameter->data->value;
                parameter=parameter->next;
            }
        }
    else 
    {
        XLALPrintError("XLAL Error - %s: parameter '%s' unknown!\n",
                __func__, name);
        XLAL_ERROR(XLAL_EINVAL);
    }
    return 0.0; // Should not actually get here!
}

/**
 * Function that checks whether the requested parameter exists within the
 * test GR parameters linked list.  Returns true (1) or false (0) accordingly
 */
bool XLALSimInspiralTestGRParamExists(
        const LALSimInspiralTestGRParam *parameter, 	/**< Linked list to check */
        const char *name 		/**< Parameter name to check for */
        )
{
  if(!parameter) return false;
  while(parameter) {if(!strcmp(parameter->data->name, name)) return true; else parameter=parameter->next;}
  return false;
}

/** Function that prints the whole test GR params linked list */
int XLALSimInspiralPrintTestGRParam(
        FILE *fp, 			/** FILE pointer to write to */
        LALSimInspiralTestGRParam *parameter 	/**< Linked list to print */
        )
{
    if (parameter!=NULL)
    {
        while(parameter) 
        {
            fprintf(fp,"%s %10.5f\n",parameter->data->name,parameter->data->value);
            parameter=parameter->next;
        }
        return XLAL_SUCCESS;
    }
    else
    {
        XLALPrintError("XLAL Error - %s: parameter not allocated!\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }
}

/** Function that destroys the whole test GR params linked list */
void XLALSimInspiralDestroyTestGRParam(
        LALSimInspiralTestGRParam *parameter 	/**< Linked list to destroy */
        )
{
   LALSimInspiralTestGRParam *tmp;
   while(parameter){
	tmp=parameter->next;
	XLALFree(parameter->data);
	XLALFree(parameter);
	parameter=tmp;
	}
}
