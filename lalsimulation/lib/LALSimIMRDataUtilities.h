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

#ifndef _LALSIMIMRDATAUTILITIES_H
#define _LALSIMIMRDATAUTILITIES_H

#include <lal/XLALError.h>

#ifdef LAL_HDF5_ENABLED

#include <lal/H5FileIO.h>

int ROM_check_version_number(LALH5File *file, INT4 version_major_in, INT4 version_minor_in, INT4 version_micro_in);
int ROM_check_canonical_file_basename(LALH5File *file, const char file_name[], const char attribute[]);

#endif

#endif /* _LALSIMIMRDATAUTILITIES_H */
