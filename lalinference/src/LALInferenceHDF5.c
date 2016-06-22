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

const char LALInferenceHDF5PosteriorSamplesGroupName[]="posterior_samples";
const char LALInferenceHDF5NestedSamplesGroupName[]="nested_samples";

LALH5File *LALInferenceCreateHDF5GroupStructure(LALH5File *h5file, const char *codename, const char *runID)
{
  LALH5File *codeGroup = XLALH5GroupOpen(h5file, codename);
  LALH5File *runGroup = XLALH5GroupOpen(codeGroup, runID);
  return(runGroup);
}

int LALInferenceH5GroupToVariablesArray(LALH5File *group , LALInferenceVariables ***varsArray, UINT4 *N)
{
  char **dataset_names=NULL;
  char **fixed_names=NULL;
  UINT4 Nfixed=0;
  LALInferenceVariables **va=NULL;
  UINT4 i=0,j=0;
  UINT4 Nsamples=0;
  
  /* Read the group datasets in as arrays */
  for(i=0;dataset_names[i];i++)
  {
    char *pname=dataset_names[i];
    LALH5Dataset *dset = XLALH5DatasetRead(group, pname);
    LALTYPECODE LALtype = XLALH5DatasetQueryType(dset);
    LALInferenceParamVaryType varyType = XLALH5DatasetQueryINT4AttributeValue(dset,"vary_type");;
    switch(LALtype)
    {
      case(LAL_D_TYPE_CODE):
      {
        REAL8Vector *vector = XLALH5DatasetReadREAL8Vector(dset);
        for(j=0;j<vector->length;j++) LALInferenceAddVariable(va[j],pname,&(vector->data[j]),LALINFERENCE_REAL8_t,varyType);
        XLALDestroyREAL8Vector(vector);
        break;
      }
      case(LAL_S_TYPE_CODE):
      {
        REAL4Vector *vector = XLALH5DatasetReadREAL4Vector(dset);
        for(j=0;j<vector->length;j++) LALInferenceAddVariable(va[j],pname,&(vector->data[j]),LALINFERENCE_REAL4_t,varyType);
        XLALDestroyREAL4Vector(vector);
        break;
      }
      case(LAL_C_TYPE_CODE):
      {
        COMPLEX8Vector *vector = XLALH5DatasetReadCOMPLEX8Vector(dset);
        for(j=0;j<vector->length;j++) LALInferenceAddVariable(va[j],pname,&(vector->data[j]),LALINFERENCE_COMPLEX8_t,varyType);
        XLALDestroyCOMPLEX8Vector(vector);
        break;
      }
      case(LAL_Z_TYPE_CODE):
      {
        COMPLEX16Vector *vector = XLALH5DatasetReadCOMPLEX16Vector(dset);
        for(j=0;j<vector->length;j++) LALInferenceAddVariable(va[j],pname,&(vector->data[j]),LALINFERENCE_COMPLEX16_t,varyType);
        XLALDestroyCOMPLEX16Vector(vector);
        break;
      }
      case(LAL_I4_TYPE_CODE):
      {
        INT4Vector *vector = XLALH5DatasetReadINT4Vector(dset);
        for(j=0;j<vector->length;j++) LALInferenceAddVariable(va[j],pname,&(vector->data[j]),LALINFERENCE_INT4_t,varyType);
        XLALDestroyINT4Vector(vector);
        break;
      }
      case(LAL_U4_TYPE_CODE):
      {
        UINT4Vector *vector = XLALH5DatasetReadUINT4Vector(dset);
        for(j=0;j<vector->length;j++) LALInferenceAddVariable(va[j],pname,&(vector->data[j]),LALINFERENCE_UINT4_t,varyType);
        XLALDestroyUINT4Vector(vector);
        break;
      }
      default:
      {
        XLALPrintWarning("%s: Unknown type code %i\n",__func__,LALtype);
        break;
      }
    } /* End switch */
  }
  Nsamples = j;
  /* end loop over parameter names */
  
  /* Read fixed parameters from attributes */
  XLALH5FileGetAttributeNames(group, &fixed_names, &Nfixed);
  for(i=0;i<Nfixed;i++)
  {
    char *pname=fixed_names[i];
    LALTYPECODE LALtype = XLALH5FileQueryScalarAttributeType(group, pname);
    switch(LALtype)
    {
      case(LAL_D_TYPE_CODE):
      {
        REAL8 value=0.0;
        XLALH5FileQueryScalarAttributeValue(&value, group,pname);
        for(j=0;j<Nsamples;j++) LALInferenceAddREAL8Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_S_TYPE_CODE):
      {
        REAL4 value=0.0;
        XLALH5FileQueryScalarAttributeValue(&value, group,pname);
        for(j=0;j<Nsamples;j++) LALInferenceAddREAL4Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_C_TYPE_CODE):
      {
        COMPLEX8 value=0.0;
        XLALH5FileQueryScalarAttributeValue(&value, group,pname);
        for(j=0;j<Nsamples;j++) LALInferenceAddCOMPLEX8Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_Z_TYPE_CODE):
      {
        COMPLEX16 value=0.0;
        XLALH5FileQueryScalarAttributeValue(&value, group,pname);
        for(j=0;j<Nsamples;j++) LALInferenceAddCOMPLEX16Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_I4_TYPE_CODE):
      {
        INT4 value=0;
        XLALH5FileQueryScalarAttributeValue(&value, group,pname);
        for(j=0;j<Nsamples;j++) LALInferenceAddINT4Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_U4_TYPE_CODE):
      {
        UINT4 value=0;
        XLALH5FileQueryScalarAttributeValue(&value, group,pname);
        for(j=0;j<Nsamples;j++) LALInferenceAddUINT4Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      default:
      {
        XLALPrintWarning("%s: Unknown type code %i\n",__func__,LALtype);
        break;
      }
    } /* End switch */
  } /* End loop over fixed_params */

  
  /* Construct the array of LALInferenceVariables */
  *varsArray = va;
  *N = Nsamples;
  return(XLAL_SUCCESS);
}

int LALInferenceVariablesArray2H5Group(LALH5File *h5file, LALInferenceVariables *const *const varsArray, UINT4 N, const char *GroupName)
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
  if(!groupPtr)
  {
    XLAL_ERROR(XLAL_EFAILED,"Error adding group %s to HDF5 file\n",GroupName);
  }
  
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
    LALH5Dataset *dset=NULL;
    switch(type)
    {
      case LALINFERENCE_REAL8_t:
      {
        REAL8Vector *vector=XLALCreateREAL8Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        dset = XLALH5DatasetAllocREAL8Vector(groupPtr, column_names[i], vector);
        XLALDestroyREAL8Vector(vector);
        break;
      }
      case LALINFERENCE_REAL4_t:
      {
        REAL4Vector *vector=XLALCreateREAL4Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        dset = XLALH5DatasetAllocREAL4Vector(groupPtr, column_names[i], vector);
        XLALDestroyREAL4Vector(vector);
        break;
      }
      case LALINFERENCE_UINT4_t:
      {
        UINT4Vector *vector=XLALCreateUINT4Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        dset = XLALH5DatasetAllocUINT4Vector(groupPtr, column_names[i], vector);
        XLALDestroyUINT4Vector(vector);
        break;
      }
      case LALINFERENCE_INT4_t:
      {
        INT4Vector *vector=XLALCreateINT4Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        dset = XLALH5DatasetAllocINT4Vector(groupPtr, column_names[i], vector);
        XLALDestroyINT4Vector(vector);
        break;
      }
      case LALINFERENCE_COMPLEX8_t:
      {
        COMPLEX8Vector *vector=XLALCreateCOMPLEX8Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        dset = XLALH5DatasetAllocCOMPLEX8Vector(groupPtr, column_names[i], vector);
        XLALDestroyCOMPLEX8Vector(vector);
        break;
      }
      case LALINFERENCE_COMPLEX16_t:
      {
        COMPLEX16Vector *vector=XLALCreateCOMPLEX16Vector(N);
        for(j=0;j<N;j++) memcpy(&(vector->data[j]),LALInferenceGetVariable(varsArray[j],column_names[i]),LALInferenceTypeSize[type]);
        dset = XLALH5DatasetAllocCOMPLEX16Vector(groupPtr, column_names[i], vector);
        XLALDestroyCOMPLEX16Vector(vector);
        break;
      }
      default:
      {
        XLALPrintWarning("LALInferenceType %i for parameter %s not implemented for HDF5, ignoring\n",type,column_names[i]);
        break;
      }
    } /* End switch */
    if(!dset) XLAL_ERROR(XLAL_EFAILED,"Failed to write HDF5 Dataset for parameter %s\n",column_names[i]);
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(varsArray[0],column_names[i]);
    XLALH5DatasetAddScalarAttribute(dset, "vary_type", &vtype , LAL_I4_TYPE_CODE );
  }/* End loop over varying parameters */
  
  for(i=0;i<Nfixed;i++)
  {
    LALInferenceVariableToAttribute(groupPtr, varsArray[0], fixed_names[i]);
  } /* End loop over fixed parameters */

  return(XLAL_SUCCESS);
  
}

int LALInferenceVariableToAttribute(LALH5File *filePtr, LALInferenceVariables *vars, char *name)
{
  if(filePtr==NULL || vars==NULL)
  {
    XLAL_ERROR(XLAL_EFAULT, "%s: Received NULL pointer\n",__func__);
  }
  LALInferenceVariableType type=LALInferenceGetVariableType(vars,name);
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
      XLALPrintWarning("LALInferenceType %i for parameter %s not implemented for HDF5 attribute, ignoring\n",type,name);
      break;
    }
  } /* End switch */
  XLALH5FileAddScalarAttribute(filePtr, name, LALInferenceGetVariable(vars,name), laltype);
  return(XLAL_SUCCESS);
}