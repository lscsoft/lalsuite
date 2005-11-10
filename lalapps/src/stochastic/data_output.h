/*
 * data_output.h - SGWB Standalone Analysis Pipeline
 *               - Data Output Function Prototypes
 * 
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <getopt.h>

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameStream.h>
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
    LALLeapSecAccuracy accuracy,
    CHAR *program_name,
    CHAR *output_path,
    CHAR *base_name,
    StochasticTable *stochtable,
    MetadataTable proctable,
    MetadataTable procparams,
    ProcessParamsTable *this_proc_param,
    CHAR comment[LIGOMETA_COMMENT_MAX]);

/*
 * vim: et
 */
