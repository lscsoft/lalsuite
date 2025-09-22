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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#ifndef LALInferenceHDF5_h
#define LALInferenceHDF5_h

#include <lal/LALInference.h>
#include <lal/H5FileIO.h>

/**
 * Returns 1 if a non-empty file exists, 0 otherwise
 */
int LALInferenceCheckNonEmptyFile(char *filename);

/**
 * Prints the size of the file
 */
int LALInferencePrintCheckpointFileInfo(char *filename);

int LALInferenceH5VariablesArrayToDataset(
    LALH5File *h5file, LALInferenceVariables *const *const varsArray, UINT4 N,
    const char *TableName);

int LALInferenceH5DatasetToVariablesArray(
    LALH5Dataset *dataset, LALInferenceVariables ***varsArray, UINT4 *N);

/**
 * Create a HDF5 heirarchy in the given LALH5File reference
 * /codename/runID/
 * Returns a LALH5File pointer to the runID group.
 */
LALH5File *LALInferenceH5CreateGroupStructure(
    LALH5File *h5file, const char *codename, const char *runID);

extern const char *const LALInferenceHDF5PosteriorSamplesDatasetName;
extern const char *const LALInferenceHDF5NestedSamplesDatasetName;

#endif /* LALInferenceHDF5_h */
