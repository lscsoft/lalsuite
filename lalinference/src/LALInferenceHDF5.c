/*
 *  Copyright (C) 2016 John Veitch
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

#include <lal/LALInference.h>
#include <lal/H5FileIO.h>
#include "LALInferenceHDF5.h"

int XLALInferenceVariablesArray2H5Group(LALH5File *h5file, LALInferenceVariables *const *const varsArray, UINT4 N, const char *GroupName)
{
  LALH5File *groupPtr=NULL;
  LALInferenceVariableItem *varitem=NULL;
  UINT4 i=0,j=0;
  
  /* Sanity check input */
  if(!varsArray) {
    XLALPrintError("Received null varsArray pointer");
    XLAL_ERROR(XLAL_EFAULT);
  }
  if(!h5file)
  {
    XLALPrintError("Received null h5file pointer\n");
    XLAL_ERROR(XLAL_EFAULT);
  }
  if(N==0) return(0);

  /* H5GroupOpen */
  groupPtr = XLALH5GroupOpen(h5file, GroupName);
  
  char *column_names[varsArray[0]->dimension];
  UINT4 Nvary=0;
  char *fixed_names[varsArray[0]->dimension];
  UINT4 Nfixed=0;
  
  /* Build a list of PARAM and FIELD elements */
  
  for(varitem=varsArray[0]->head;varitem;varitem=varitem->next)
  {
    switch(varitem->vary){
      case LALINFERENCE_PARAM_LINEAR:
      case LALINFERENCE_PARAM_CIRCULAR:
      case LALINFERENCE_PARAM_OUTPUT:
      {
        column_names[Nvary++]=varitem->name;
        break;
      }
      case LALINFERENCE_PARAM_FIXED:
      {
        fixed_names[Nfixed++]=varitem->name;
        break;
      }
      default:
      {
        XLALPrintWarning("Unknown param vary type");
      }
    }
  }

  for(i=0;i<Nvary;i++)
  {
    LALInferenceVariableType type=LALInferenceGetVariableType(varsArray[0],column_names[i]);
    switch(type)
    {
      case LALINFERENCE_REAL8_t:
      {
        REAL8Vector *vector=XLALCreateREAL8Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        XLALH5FileWriteREAL8Vector(groupPtr, column_names[i], vector);
        XLALDestroyREAL8Vector(vector);
        break;
      }
      case LALINFERENCE_REAL4_t:
      {
        REAL4Vector *vector=XLALCreateREAL4Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        XLALH5FileWriteREAL4Vector(groupPtr, column_names[i], vector);
        XLALDestroyREAL4Vector(vector);
        break;
      }
      case LALINFERENCE_UINT4_t:
      {
        UINT4Vector *vector=XLALCreateUINT4Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        XLALH5FileWriteUINT4Vector(groupPtr, column_names[i], vector);
        XLALDestroyUINT4Vector(vector);
        break;
      }
      case LALINFERENCE_INT4_t:
      {
        INT4Vector *vector=XLALCreateINT4Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        XLALH5FileWriteINT4Vector(groupPtr, column_names[i], vector);
        XLALDestroyINT4Vector(vector);
        break;
      }
      case LALINFERENCE_COMPLEX8_t:
      {
        COMPLEX8Vector *vector=XLALCreateCOMPLEX8Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        XLALH5FileWriteCOMPLEX8Vector(groupPtr, column_names[i], vector);
        XLALDestroyCOMPLEX8Vector(vector);
        break;
      }
      case LALINFERENCE_COMPLEX16_t:
      {
        COMPLEX16Vector *vector=XLALCreateCOMPLEX16Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        XLALH5FileWriteCOMPLEX16Vector(groupPtr, column_names[i], vector);
        XLALDestroyCOMPLEX16Vector(vector);
        break;
      }
      default:
      {
        XLALPrintWarning("LALInferenceType %i for parameter %s not implemented for HDF5, ignoring\n",type,column_names[i]);
        break;
      }
    }
  }
  
  for(i=0;i<Nfixed;i++)
  {
    LALInferenceVariableType type=LALInferenceGetVariableType(varsArray[0],fixed_names[i]);
    LALTYPECODE laltype=LAL_D_TYPE_CODE;
    switch(type)
    {
      case LALINFERENCE_UINT4_t:
      {
        laltype = LAL_U4_TYPE_CODE;
        break;
      }
      case LALINFERENCE_REAL8_t:
      {
        laltype = LAL_D_TYPE_CODE;
        break;
      }
      case LALINFERENCE_REAL4_t:
      {
        laltype = LAL_S_TYPE_CODE;
        break;
      }
      case LALINFERENCE_INT4_t:
      {
        laltype = LAL_I4_TYPE_CODE;
        break;
      }
      default:
      {
        XLALPrintWarning("LALInferenceType %i for parameter %s not implemented for HDF5 attribute, ignoring\n",type,fixed_names[i]);
        break;
      }
    }
    XLALH5FileAddScalarAttribute(groupPtr, fixed_names[i], LALInferenceGetVariable(varsArray[0],fixed_names[i]), laltype);
  }

  return(XLAL_SUCCESS);
  
}
