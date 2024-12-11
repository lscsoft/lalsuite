/*
 *  Copyright (C) 2024  Lorenzo Pompili, Michael Puerrer
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

/**
 * \author Lorenzo Pompili, Michael Puerrer
 *
 * \file
 *
 * \brief Auxiliary functions related to HDF5 waveform data files.
 */

#include <string.h>

#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include "LALSimIMRDataUtilities.h"

#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>


int ROM_check_version_number(LALH5File *file, INT4 version_major_in, INT4 version_minor_in, INT4 version_micro_in) {
  INT4 version_major;
  INT4 version_minor;
  INT4 version_micro;

  LALH5Generic gfile = {.file = file};
  XLALH5AttributeQueryScalarValue(&version_major, gfile, "version_major");
  XLALH5AttributeQueryScalarValue(&version_minor, gfile, "version_minor");
  XLALH5AttributeQueryScalarValue(&version_micro, gfile, "version_micro");

  if ((version_major_in != version_major) || (version_minor_in != version_minor) || (version_micro_in != version_micro)) {
    XLAL_ERROR(XLAL_EIO, "Expected ROM data version %d.%d.%d, but got version %d.%d.%d.",
    version_major_in, version_minor_in, version_micro_in, version_major, version_minor, version_micro);
  }
  else {
    XLALPrintInfo("Reading ROM data version %d.%d.%d.\n", version_major, version_minor, version_micro);
    return XLAL_SUCCESS;
  }
}

int ROM_check_canonical_file_basename(LALH5File *file, const char file_name[], const char attribute[]) {

  LALH5Generic gfile = {.file = file};
  int len = XLALH5AttributeQueryStringValue(NULL, 0, gfile, attribute) + 1;
  char *canonical_file_basename = XLALMalloc(len);
  XLALH5FileQueryStringAttributeValue(canonical_file_basename, len, file, attribute); 
  
  if (strcmp(canonical_file_basename, file_name) != 0) {
    XLAL_ERROR(XLAL_EIO, "Expected CANONICAL_FILE_BASENAME %s, but got %s.",
    file_name, canonical_file_basename);
  }
  else {
    XLALPrintInfo("ROM canonical_file_basename %s\n", canonical_file_basename);
  }
  XLALFree(canonical_file_basename);
  return XLAL_SUCCESS;
}

#endif /* LAL_HDF5_ENABLED */
