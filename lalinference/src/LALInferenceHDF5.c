/*
 *  Copyright (C) 2016 John Veitch and Leo Singer
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
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>
#include <stdlib.h>

const char LALInferenceHDF5PosteriorSamplesDatasetName[]="posterior_samples";
const char LALInferenceHDF5NestedSamplesDatasetName[]="nested_samples";


static void assert_not_reached(void) __attribute__ ((noreturn));
static void assert_not_reached(void)
{
#ifndef NDEBUG
  abort();
#endif
}



/* BEGIN COPIED FROM H5FILEIOLowLevel.c */

struct tagLALH5Object {
	hid_t object_id; /* this object's id must be first */
};

struct tagLALH5File {
	hid_t file_id; /* this object's id must be first */
	unsigned int mode;
	int is_a_group;
	char fname[FILENAME_MAX];
};

struct tagLALH5Dataset {
	hid_t dataset_id; /* this object's id must be first */
	hid_t parent_id;
	hid_t space_id;
	hid_t dtype_id; /* note: this is the in-memory type */
	char name[]; /* flexible array member must be last */
};

/* creates HDF5 float complex data type; use H5Tclose() to free */
typedef struct { float re; float im; } internal_float_complex_type;
static hid_t XLALH5TypeNativeFloatComplex(void)
{
	hid_t dtype_id;
	dtype_id = H5Tcreate(H5T_COMPOUND, sizeof(internal_float_complex_type));
	H5Tinsert(dtype_id, "r", HOFFSET(internal_float_complex_type, re), H5T_NATIVE_FLOAT);
	H5Tinsert(dtype_id, "i", HOFFSET(internal_float_complex_type, im), H5T_NATIVE_FLOAT);
	return dtype_id;
}

/* creates HDF5 double complex data type; use H5Tclose() to free */
typedef struct { double re; double im; } internal_double_complex_type;
static hid_t XLALH5TypeNativeDoubleComplex(void)
{
	hid_t dtype_id;
	dtype_id = H5Tcreate(H5T_COMPOUND, sizeof(internal_double_complex_type));
	H5Tinsert(dtype_id, "r", HOFFSET(internal_double_complex_type, re), H5T_NATIVE_DOUBLE);
	H5Tinsert(dtype_id, "i", HOFFSET(internal_double_complex_type, im), H5T_NATIVE_DOUBLE);
	return dtype_id;
}

/* END COPIED FROM H5FILEIOLowLevel.c */



static int LALInferenceH5VariableToAttribute(LALH5Dataset *dataset, LALInferenceVariables *vars, char *name);

LALH5File *LALInferenceH5CreateGroupStructure(LALH5File *h5file, const char *codename, const char *runID)
{
  LALH5File *codeGroup = XLALH5GroupOpen(h5file, codename);
  LALH5File *runGroup = XLALH5GroupOpen(codeGroup, runID);
  XLALH5FileClose(codeGroup);
  return(runGroup);
}

int LALInferenceH5DatasetToVariablesArray(LALH5Dataset *dataset, LALInferenceVariables ***varsArray, UINT4 *N)
{
  hid_t type_id = dataset->dtype_id;
  size_t type_size = H5Tget_size(type_id);
  UINT4 Nvary = H5Tget_nmembers(type_id);

  int vary[Nvary];
  char *column_names[Nvary];
  size_t column_offsets[Nvary];
  LALInferenceVariableType column_types[Nvary];
  int ret;

  for (UINT4 i = 0; i < Nvary; i ++)
  {
    column_names[i] = H5Tget_member_name(type_id, i);
    column_offsets[i] = H5Tget_member_offset(type_id, i);
    hid_t subtype_id = H5Tget_member_type(type_id, i);
    assert(subtype_id >= 0);
    H5T_class_t subtype_class = H5Tget_class(subtype_id);
    assert(subtype_class >= 0);
    hid_t subtype_native_type = H5Tget_native_type(subtype_id, H5T_DIR_ASCEND);
    assert(subtype_native_type >= 0);

    switch(subtype_class)
    {
      case H5T_INTEGER:
      case H5T_FLOAT:
        if (H5Tequal(subtype_native_type, H5T_NATIVE_DOUBLE))
          column_types[i] = LALINFERENCE_REAL8_t;
        else if (H5Tequal(subtype_native_type, H5T_NATIVE_FLOAT))
          column_types[i] = LALINFERENCE_REAL4_t;
        else if (H5Tequal(subtype_native_type, H5T_NATIVE_UINT))
          column_types[i] = LALINFERENCE_UINT4_t;
        else if (H5Tequal(subtype_native_type, H5T_NATIVE_INT))
          column_types[i] = LALINFERENCE_INT4_t;
        else
          assert_not_reached();
        break;
      case H5T_COMPOUND:
        assert(H5Tget_nmembers(subtype_id) == 2);
        if (H5Tequal(subtype_native_type, H5T_NATIVE_FLOAT))
          column_types[i] = LALINFERENCE_COMPLEX8_t;
        else if (H5Tequal(subtype_native_type, H5T_NATIVE_DOUBLE))
          column_types[i] = LALINFERENCE_COMPLEX16_t;
        else
          assert_not_reached();
        break;
      default:
        assert_not_reached();
    }

    H5Tclose(subtype_native_type);
    H5Tclose(subtype_id);
  }

  size_t nbytes = XLALH5DatasetQueryNBytes(dataset);
  char *data = XLALMalloc(nbytes);
  assert(data);
  ret = XLALH5DatasetQueryData(data, dataset);
  assert(ret == 0);

  char **fixed_names=NULL;
  UINT4 Nfixed=0;
  LALInferenceVariables **va=NULL;
  UINT4 Nsamples=XLALH5DatasetQueryNPoints(dataset);
  
  va=XLALCalloc(Nsamples,sizeof(LALInferenceVariables *));
  for(UINT4 i=0;i<Nsamples;i++) va[i]=XLALCalloc(1,sizeof(LALInferenceVariables));

  hsize_t vary_attr_type_dims[] = {-1};
  hid_t vary_attr_id = H5Aopen(dataset->dataset_id, "vary", H5P_DEFAULT);
  assert(vary_attr_id >= 0);
  hid_t vary_attr_type_id = H5Aget_type(vary_attr_id);
  assert(vary_attr_type_id >= 0);
  assert(H5Tget_class(vary_attr_type_id) == H5T_ARRAY);
  assert(H5Tget_array_ndims(vary_attr_type_id) == 1);
  H5Tget_array_dims2(vary_attr_type_id, &vary_attr_type_dims[0]);
  assert(vary_attr_type_dims[0] == Nvary);
  assert(H5Tget_size(vary_attr_type_id) == sizeof(vary));
  herr_t read_status = H5Aread(vary_attr_id, vary_attr_type_id, vary);
  assert(read_status >= 0);
  H5Tclose(vary_attr_type_id);
  H5Aclose(vary_attr_id);

  /* Read the group datasets in as arrays */
  for(UINT4 i=0;i<Nsamples;i++)
    for(UINT4 j=0;j<Nvary;j++)
      LALInferenceAddVariable(va[i], column_names[j], data + type_size * i + column_offsets[j], column_types[j], vary[j]);
  XLALFree(data);

  for (UINT4 i=0;i<Nvary;i++)
    free(column_names[i]);

  /* Read fixed parameters from attributes */
  {
    /* FIXME: no counterpart of XLALH5FileGetAttributeNames for datasets */
    LALH5File fakefile = {.file_id = dataset->dataset_id};
    XLALH5FileGetAttributeNames(&fakefile, &fixed_names, &Nfixed);
  }
  for(UINT4 i=0;i<Nfixed;i++)
  {
    char *pname=fixed_names[i];
    if (strcmp(pname, "vary") == 0)
      continue;
    LALTYPECODE LALtype = XLALH5DatasetQueryScalarAttributeType(dataset, pname);
    switch(LALtype)
    {
      case(LAL_D_TYPE_CODE):
      {
        REAL8 value=0.0;
        XLALH5DatasetQueryScalarAttributeValue(&value, dataset,pname);
        for(UINT4 j=0;j<Nsamples;j++) LALInferenceAddREAL8Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_S_TYPE_CODE):
      {
        REAL4 value=0.0;
        XLALH5DatasetQueryScalarAttributeValue(&value, dataset,pname);
        for(UINT4 j=0;j<Nsamples;j++) LALInferenceAddREAL4Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_C_TYPE_CODE):
      {
        COMPLEX8 value=0.0;
        XLALH5DatasetQueryScalarAttributeValue(&value, dataset,pname);
        for(UINT4 j=0;j<Nsamples;j++) LALInferenceAddCOMPLEX8Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_Z_TYPE_CODE):
      {
        COMPLEX16 value=0.0;
        XLALH5DatasetQueryScalarAttributeValue(&value, dataset,pname);
        for(UINT4 j=0;j<Nsamples;j++) LALInferenceAddCOMPLEX16Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_I4_TYPE_CODE):
      {
        INT4 value=0;
        XLALH5DatasetQueryScalarAttributeValue(&value, dataset,pname);
        for(UINT4 j=0;j<Nsamples;j++) LALInferenceAddINT4Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
        break;
      }
      case(LAL_U4_TYPE_CODE):
      {
        UINT4 value=0;
        XLALH5DatasetQueryScalarAttributeValue(&value, dataset,pname);
        for(UINT4 j=0;j<Nsamples;j++) LALInferenceAddUINT4Variable(va[j],pname,value,LALINFERENCE_PARAM_FIXED);
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
  for(UINT4 i=0;i<Nfixed;i++) XLALFree(fixed_names[i]);
  XLALFree(fixed_names);
  return(XLAL_SUCCESS);
}

static void enum_insert(hid_t type, const char *name, int value)
{
  herr_t ret = H5Tenum_insert(type, name, &value);
  assert(ret >= 0);
}

static hid_t create_vary_type(void)
{
  hid_t type = H5Tenum_create(H5T_NATIVE_INT);
  assert(type >= 0);

	enum_insert(type, "linear", LALINFERENCE_PARAM_LINEAR);
	enum_insert(type, "circular", LALINFERENCE_PARAM_CIRCULAR);
	enum_insert(type, "fixed", LALINFERENCE_PARAM_FIXED);
	enum_insert(type, "output", LALINFERENCE_PARAM_OUTPUT);

  return type;
}

int LALInferenceH5VariablesArrayToDataset(LALH5File *h5file, LALInferenceVariables *const *const varsArray, UINT4 N, const char *TableName)
{
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

  const char *column_names[varsArray[0]->dimension];
  UINT4 Nvary=0;
  hsize_t type_size=0;
  size_t column_offsets[varsArray[0]->dimension];
  size_t column_sizes[varsArray[0]->dimension];
  hid_t column_types[varsArray[0]->dimension];
  char *fixed_names[varsArray[0]->dimension];
  int vary[varsArray[0]->dimension];
  UINT4 Nfixed=0;

  hid_t native_floatcomplex = XLALH5TypeNativeFloatComplex();
  hid_t native_doublecomplex = XLALH5TypeNativeDoubleComplex();

  /* Build a list of PARAM and FIELD elements */
  for(LALInferenceVariableItem *varitem=varsArray[0]->head;varitem;varitem=varitem->next)
  {
    switch(varitem->vary){
      case LALINFERENCE_PARAM_LINEAR:
      case LALINFERENCE_PARAM_CIRCULAR:
      case LALINFERENCE_PARAM_OUTPUT:
      {
        hid_t tp;
        size_t sz;
        switch(varitem->type)
        {
          case LALINFERENCE_REAL8_t:
            tp = H5T_NATIVE_DOUBLE; sz = sizeof(REAL8); break;
          case LALINFERENCE_REAL4_t:
            tp = H5T_NATIVE_FLOAT; sz = sizeof(REAL4); break;
          case LALINFERENCE_UINT4_t:
            tp = H5T_NATIVE_UINT32; sz = sizeof(UINT4); break;
          case LALINFERENCE_INT4_t:
            tp = H5T_NATIVE_INT32; sz = sizeof(INT4); break;
          case LALINFERENCE_COMPLEX8_t:
            tp = native_floatcomplex; sz = sizeof(COMPLEX8); break;
          case LALINFERENCE_COMPLEX16_t:
            tp = native_doublecomplex; sz = sizeof(COMPLEX16); break;
          default:
            XLALPrintWarning("LALInferenceType %i for parameter %s not implemented for HDF5, ignoring\n",varitem->type,varitem->name);
            continue;
        } /* End switch */
        vary[Nvary] = varitem->vary;
        column_types[Nvary]=tp;
        column_sizes[Nvary]=sz;
        column_offsets[Nvary]=type_size;
        type_size+=sz;
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

  /* Create table */
  LALH5Dataset *dataset = LALCalloc(1, sizeof(LALH5Dataset));
  assert(dataset);
  dataset->dtype_id = H5Tcreate(H5T_COMPOUND, type_size);
  assert(dataset->dtype_id >= 0);
  for (UINT4 i = 0; i < Nvary; i ++)
  {
    herr_t ret = H5Tinsert(
      dataset->dtype_id, column_names[i], column_offsets[i], column_types[i]);
    assert(ret >= 0);
  }
  {
    hsize_t dims[] = {N};
    dataset->space_id = H5Screate_simple(1, dims, NULL);
    assert(dataset->space_id >= 0);
  }
  dataset->dataset_id = H5Dcreate2(
    h5file->file_id, TableName, dataset->dtype_id, dataset->space_id,
  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(dataset->dataset_id >= 0);

  /* Gather together data in one big array */
  char *data = XLALCalloc(N, type_size);
  assert(data);
  for(UINT4 i=0;i<N;i++)
  {
    for(UINT4 j=0;j<Nvary;j++)
    {
      void *var = LALInferenceGetVariable(varsArray[i], column_names[j]);
      memcpy(data + type_size * i + column_offsets[j], var, column_sizes[j]);
    }
  }

  /* Write data */
  {
    int ret = XLALH5DatasetWrite(dataset, data);
    assert(ret == XLAL_SUCCESS);
  }
  XLALFree(data);

  hid_t vary_type_id = create_vary_type();
  assert(vary_type_id >= 0);
  hsize_t vary_array_type_dims[] = {Nvary};
  hid_t vary_array_type_id = H5Tarray_create2(
    vary_type_id, 1, vary_array_type_dims);
  assert(vary_array_type_id >= 0);
  hid_t space_id = H5Screate(H5S_SCALAR);
  assert(space_id >= 0);
  hid_t attr_id = H5Acreate2(
    dataset->dataset_id, "vary", vary_array_type_id, space_id,
    H5P_DEFAULT, H5P_DEFAULT);
  assert(attr_id >= 0);
  herr_t attr_write_status = H5Awrite(attr_id, vary_array_type_id, &vary);
  assert(attr_write_status >= 0);
  H5Aclose(attr_id);
  H5Sclose(space_id);
  H5Tclose(vary_array_type_id);
  H5Tclose(vary_type_id);

  /* Write attributes, if any */
  for(UINT4 i=0;i<Nfixed;i++)
  {
    int ret = LALInferenceH5VariableToAttribute(dataset, varsArray[0], fixed_names[i]);
    assert(ret == XLAL_SUCCESS);
  }

  H5Tclose(native_floatcomplex);
  H5Tclose(native_doublecomplex);

  XLALH5DatasetFree(dataset);

  return(XLAL_SUCCESS);
}

static int LALInferenceH5VariableToAttribute(LALH5Dataset *dataset, LALInferenceVariables *vars, char *name)
{
  if(dataset==NULL || vars==NULL)
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
  XLALH5DatasetAddScalarAttribute(dataset, name, LALInferenceGetVariable(vars,name), laltype);
  return(XLAL_SUCCESS);
}

