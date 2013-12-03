/*
 * data_output.h - SGWB Standalone Analysis Pipeline
 *               - Data Output Function Prototypes
 *
 * Copyright (C) 2002-2006,2009-2010 Adam Mercer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 */

#ifndef DATA_OUTPUT_H
#define DATA_OUTPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <getopt.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/LALFrStream.h>
#include <lal/LALFrameL.h>
#include <lal/LALStdio.h>
#include <lal/FrequencySeries.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Units.h>

#include <lalapps.h>
#include <processtable.h>

/* save out ccSpectra as a frame file */
void write_ccspectra_frame(COMPLEX8FrequencySeries *series,
    CHAR *ifo_one,
    CHAR *ifo_two,
    LIGOTimeGPS epoch,
    INT4 duration);

/* save out xml tables */
void save_xml_file(LALStatus *status,
    const CHAR *program_name,
    CHAR *output_path,
    CHAR *base_name,
    StochasticTable *stoch_table,
    MetadataTable proc_table,
    MetadataTable proc_params,
    ProcessParamsTable *this_proc_param,
    CHAR comment[LIGOMETA_COMMENT_MAX]);

#endif /* DATA_OUTPUT_H */

/*
 * vim: et
 */
