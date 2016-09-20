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


const char LALInferenceHDF5PosteriorSamplesDatasetName[] = "posterior_samples";
const char LALInferenceHDF5NestedSamplesDatasetName[] = "nested_samples";


static void assert_not_reached(void) __attribute__ ((noreturn));
static void assert_not_reached(void)
{
#ifndef NDEBUG
    abort();
#endif
}


static void LALInferenceH5VariableToAttribute(
    LALH5Generic gdataset, LALInferenceVariables *vars, char *name);


LALH5File *LALInferenceH5CreateGroupStructure(
    LALH5File *h5file, const char *codename, const char *runID)
{
    LALH5File *codeGroup = XLALH5GroupOpen(h5file, codename);
    LALH5File *runGroup = XLALH5GroupOpen(codeGroup, runID);
    XLALH5FileClose(codeGroup);
    return(runGroup);
}


int LALInferenceH5DatasetToVariablesArray(
    LALH5Dataset *dataset, LALInferenceVariables ***varsArray, UINT4 *N)
{
    size_t type_size = XLALH5TableQueryRowSize(dataset);
    size_t Nvary = XLALH5TableQueryNColumns(dataset);
    LALH5Generic gdataset = {.dset = dataset};

    int vary[Nvary];
    char *column_names[Nvary];
    size_t column_offsets[Nvary];
    LALInferenceVariableType column_types[Nvary];
    int ret;

    for (size_t i = 0; i < Nvary; i ++)
    {
        size_t column_name_len = XLALH5TableQueryColumnName(
            NULL, 0, dataset, i);
        column_names[i] = malloc(column_name_len + 1);
        XLALH5TableQueryColumnName(
            column_names[i], column_name_len + 1, dataset, i);
        column_offsets[i] = XLALH5TableQueryColumnOffset(dataset, i);
        LALTYPECODE column_type = XLALH5TableQueryColumnType(dataset, i);

        switch (column_type)
        {
            case LAL_D_TYPE_CODE:
                column_types[i] = LALINFERENCE_REAL8_t; break;
            case LAL_S_TYPE_CODE:
                column_types[i] = LALINFERENCE_REAL4_t; break;
            case LAL_U4_TYPE_CODE:
                column_types[i] = LALINFERENCE_UINT4_t; break;
            case LAL_I4_TYPE_CODE:
                column_types[i] = LALINFERENCE_INT4_t; break;
            case LAL_Z_TYPE_CODE:
                column_types[i] = LALINFERENCE_COMPLEX16_t; break;
            case LAL_C_TYPE_CODE:
                column_types[i] = LALINFERENCE_COMPLEX8_t; break;
            default:
                assert_not_reached();
        }
    }

    size_t nbytes = XLALH5DatasetQueryNBytes(dataset);
    char *data = XLALMalloc(nbytes);
    assert(data);
    ret = XLALH5DatasetQueryData(data, dataset);
    assert(ret == 0);

    LALInferenceVariables **va = NULL;
    UINT4 Nsamples = XLALH5DatasetQueryNPoints(dataset);

    va = XLALCalloc(Nsamples, sizeof(LALInferenceVariables *));
    for (size_t i = 0; i < Nsamples; i++)
        va[i] = XLALCalloc(1, sizeof(LALInferenceVariables));

    for (UINT4 i = 0; i < Nvary; i ++)
    {
        char pname[] = "FIELD_NNN_VARY";
        snprintf(pname, sizeof(pname), "FIELD_%d_VARY", i);
        INT4 value;
        ret = XLALH5AttributeQueryScalarValue(&value, gdataset, pname);
        assert(ret == 0);
        vary[i] = value;
    }

    /* Read the group datasets in as arrays */
    for (UINT4 i = 0; i < Nsamples; i++)
        for (UINT4 j = 0; j < Nvary; j++)
            LALInferenceAddVariable(
                va[i], column_names[j],
                data + type_size * i + column_offsets[j],
                column_types[j], vary[j]);
    XLALFree(data);

    for (size_t i = 0; i < Nvary; i++)
        free(column_names[i]);

    size_t Nfixed = XLALH5AttributeQueryN(gdataset);

    for (size_t i = 0; i < Nfixed; i ++)
    {
        int len = XLALH5AttributeQueryName(NULL, 0, gdataset, i);
        char pname[len + 1];
        char value[16]; /* Big enough to hold largest supported LAL type */
        XLALH5AttributeQueryName(pname, sizeof(pname), gdataset, i);

        /* Skip the "vary" attribute as well as any attribute
         * associated with the H5TB interface
         * (https://www.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html). */
        if (strcmp(pname, "CLASS") == 0 ||
            strcmp(pname, "VERSION") == 0 ||
            strcmp(pname, "TITLE") == 0 ||
            strncmp(pname, "FIELD_", 6) == 0) continue;

        LALTYPECODE laltype = XLALH5AttributeQueryScalarType(gdataset, pname);
        LALInferenceVariableType lalinftype;

        switch (laltype)
        {
            case LAL_D_TYPE_CODE:
                lalinftype = LALINFERENCE_REAL8_t; break;
            case LAL_S_TYPE_CODE:
                lalinftype = LALINFERENCE_REAL4_t; break;
            case LAL_U4_TYPE_CODE:
                lalinftype = LALINFERENCE_UINT4_t; break;
            case LAL_I4_TYPE_CODE:
                lalinftype = LALINFERENCE_INT4_t; break;
            case LAL_Z_TYPE_CODE:
                lalinftype = LALINFERENCE_COMPLEX16_t; break;
            case LAL_C_TYPE_CODE:
                lalinftype = LALINFERENCE_COMPLEX8_t; break;
            default:
                XLALPrintWarning(
                    "%s: Unknown type code %i\n", __func__, laltype);
                continue;
        }

        XLALH5AttributeQueryScalarValue(&value, gdataset, pname);
        for (UINT4 j = 0; j < Nsamples; j++)
            LALInferenceAddVariable(
                va[j], pname, value, lalinftype, LALINFERENCE_PARAM_FIXED);
    } /* End loop over fixed_params */

    /* Construct the array of LALInferenceVariables */
    *varsArray = va;
    *N = Nsamples;
    return(XLAL_SUCCESS);
}


int LALInferenceH5VariablesArrayToDataset(
    LALH5File *h5file, LALInferenceVariables *const *const varsArray, UINT4 N,
    const char *TableName)
{
    /* Sanity check input */
    if (!varsArray)
        XLAL_ERROR(XLAL_EFAULT, "Received null varsArray pointer");
    if (!h5file)
        XLAL_ERROR(XLAL_EFAULT, "Received null h5file pointer");
    if (N == 0)
        return 0;

    const char *column_names[varsArray[0]->dimension];
    UINT4 Nvary = 0;
    hsize_t type_size = 0;
    size_t column_offsets[varsArray[0]->dimension];
    size_t column_sizes[varsArray[0]->dimension];
    LALTYPECODE column_types[varsArray[0]->dimension];
    char *fixed_names[varsArray[0]->dimension];
    int vary[varsArray[0]->dimension];
    UINT4 Nfixed = 0;

    /* Build a list of PARAM and FIELD elements */
    for (LALInferenceVariableItem *varitem = varsArray[0]->head; varitem;
         varitem = varitem->next)
    {
        switch(varitem->vary)
        {
            case LALINFERENCE_PARAM_LINEAR:
            case LALINFERENCE_PARAM_CIRCULAR:
            case LALINFERENCE_PARAM_OUTPUT:
            {
                LALTYPECODE tp;
                size_t sz;
                switch (varitem->type)
                {
                    case LALINFERENCE_REAL8_t:
                        tp = LAL_D_TYPE_CODE; sz = sizeof(REAL8); break;
                    case LALINFERENCE_REAL4_t:
                        tp = LAL_S_TYPE_CODE; sz = sizeof(REAL4); break;
                    case LALINFERENCE_UINT4_t:
                        tp = LAL_U4_TYPE_CODE; sz = sizeof(UINT4); break;
                    case LALINFERENCE_INT4_t:
                        tp = LAL_I4_TYPE_CODE; sz = sizeof(INT4); break;
                    case LALINFERENCE_COMPLEX8_t:
                        tp = LAL_C_TYPE_CODE; sz = sizeof(COMPLEX8); break;
                    case LALINFERENCE_COMPLEX16_t:
                        tp = LAL_Z_TYPE_CODE; sz = sizeof(COMPLEX16); break;
                    default:
                        XLALPrintWarning(
                            "LALInferenceType %i for parameter %s not "
                            "implemented for HDF5, ignoring\n",
                            varitem->type, varitem->name);
                        continue;
                } /* End switch */
                vary[Nvary] = varitem->vary;
                column_types[Nvary] = tp;
                column_sizes[Nvary] = sz;
                column_offsets[Nvary] = type_size;
                type_size += sz;
                column_names[Nvary++] = varitem->name;
                break;
            }
            case LALINFERENCE_PARAM_FIXED:
                fixed_names[Nfixed++] = varitem->name;
                break;
            default:
                XLALPrintWarning("Unknown param vary type");
        }
    }

    /* Gather together data in one big array */
    char *data = XLALCalloc(N, type_size);
    assert(data);
    for (UINT4 i = 0; i < N; i++)
    {
        for (UINT4 j = 0; j < Nvary; j++)
        {
            void *var = LALInferenceGetVariable(varsArray[i], column_names[j]);
            memcpy(
                data + type_size * i + column_offsets[j], var, column_sizes[j]);
        }
    }

    /* Create table */
    LALH5Dataset *dataset = XLALH5TableAlloc(h5file, TableName, Nvary,
        column_names, column_types, column_offsets, type_size);
    assert(dataset);
    int ret = XLALH5TableAppend(
        dataset, column_offsets, column_sizes, N, type_size, data);
    assert(ret == 0);
    XLALFree(data);

    LALH5Generic gdataset = {.dset = dataset};
    for (UINT4 i = 0; i < Nvary; i ++)
    {
        INT4 value = vary[i];
        char pname[] = "FIELD_NNN_VARY";
        snprintf(pname, sizeof(pname), "FIELD_%d_VARY", i);
        ret = XLALH5AttributeAddScalar(
            gdataset, pname, &value, LAL_I4_TYPE_CODE);
        assert(ret == 0);
    }

    /* Write attributes, if any */
    for (UINT4 i = 0; i < Nfixed; i++)
        LALInferenceH5VariableToAttribute(
            gdataset, varsArray[0], fixed_names[i]);

    XLALH5DatasetFree(dataset);
    return XLAL_SUCCESS;
}


static void LALInferenceH5VariableToAttribute(
    LALH5Generic gdataset, LALInferenceVariables *vars, char *name)
{
    LALInferenceVariableType type = LALInferenceGetVariableType(vars, name);
    LALTYPECODE laltype;
    switch(type)
    {
        case LALINFERENCE_UINT4_t:
            laltype = LAL_U4_TYPE_CODE; break;
        case LALINFERENCE_REAL8_t:
            laltype = LAL_D_TYPE_CODE; break;
        case LALINFERENCE_REAL4_t:
            laltype = LAL_S_TYPE_CODE; break;
        case LALINFERENCE_INT4_t:
            laltype = LAL_I4_TYPE_CODE; break;
        default:
            XLALPrintWarning(
                "LALInferenceType %i for parameter %s not "
                "implemented for HDF5 attribute, ignoring", type, name);
            return;
    } /* End switch */
    XLALH5AttributeAddScalar(
        gdataset, name, LALInferenceGetVariable(vars,name), laltype);
}
