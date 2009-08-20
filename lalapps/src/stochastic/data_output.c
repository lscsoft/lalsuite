/*
 * data_output.c - SGWB Standalone Analysis Pipeline
 *               - Data Output Functions
 * 
 * Copyright (C) 2002-2006,2009 Adam Mercer
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
 * $Id$
 */

#include "data_output.h"

NRCSID(DATAOUTPUTC, "$Id$");
RCSID("$Id$");

/* externally declared variables */
extern int middle_segment_flag;
extern int apply_mask_flag;
extern int high_pass_flag;
extern int overlap_hann_flag;
extern int recentre_flag;
extern int cc_spectra_flag;
extern int vrbflg;

/* save out ccSpectra as a frame file */
void write_ccspectra_frame(COMPLEX8FrequencySeries *series,
    CHAR *ifo_one,
    CHAR *ifo_two,
    LIGOTimeGPS epoch,
    INT4 duration)
{
  /* variables */
  CHAR hertz[] = "Hz";
  CHAR frame_comment[] = "$Id$";
  CHAR frame_type[] = "CCSPECTRA";
  CHAR source[FILENAME_MAX];
  CHAR fname[FILENAME_MAX];
  CHAR units[LALUnitNameSize];
  struct FrFile *frfile;
  struct FrameH *frame;
  struct FrVect *vect;
  struct FrProcData *proc;

  /* set frame filename */
  snprintf(source, sizeof(source), "%s%s", ifo_one, ifo_two);
  snprintf(fname, sizeof(fname), "%s-%s-%d-%d.gwf", source, \
      frame_type, epoch.gpsSeconds, duration);

  /* setup frame file */
  frfile = FrFileONew(fname, 0);

  /* set frame properties */
  frame = FrameHNew(source);
  frame->run = 0;
  frame->frame = 0;
  frame->GTimeS = epoch.gpsSeconds;
  frame->GTimeN = epoch.gpsNanoSeconds;
  frame->dt = duration;

  /* allocate memory for frame */
  proc = LALCalloc(1, sizeof(*proc));
  vect = FrVectNew1D(series->name, FR_VECT_8C, series->data->length, \
      series->deltaF, hertz, units);

  /* check that memory has been allocated */
  if (!vect)
  {
    LALFree(proc);
    FrVectFree(vect);
    fprintf(stderr, "unable to allocate memory for frame.\n");
    exit(1);
  }

  /* set start frequency */
  vect->startX[0] = series->f0;

  /* set frame properties */
  FrStrCpy(&proc->name, frame_type);
  FrStrCpy(&proc->comment, frame_comment);
  proc->next = frame->procData;
  frame->procData = proc;
  proc->classe = FrProcDataDef();
  proc->type = 2;
  proc->data = vect;
  proc->subType = 0;
  proc->tRange = duration;
  proc->fRange = series->data->length * series->deltaF;

  /* copy data into frame structure */
  memcpy(vect->dataD, series->data->data, \
      series->data->length * sizeof(*series->data->data));

  /* write frame */
  FrameWrite(frame, frfile);

  /* free frame */
  FrVectFree(vect); 
  vect=NULL;

  /* end frame file */
  FrFileOEnd(frfile);
}

/* save out xml tables */
void save_xml_file(LALStatus *status,
    CHAR *program_name,
    CHAR *output_path,
    CHAR *base_name,
    StochasticTable *stoch_table,
    MetadataTable proc_table,
    MetadataTable proc_params,
    ProcessParamsTable *this_proc_param,
    CHAR comment[LIGOMETA_COMMENT_MAX])
{
  /* variables */
  MetadataTable output_table;
  CHAR xml_file_name[FILENAME_MAX];
  LIGOLwXMLStream xml_stream;

  /* save out any flags to the process params table */
  if (middle_segment_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--middle-segment");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (apply_mask_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--apply-mask");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (high_pass_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--high-pass-filter");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (overlap_hann_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--overlap-hann");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (recentre_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--recentre");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (cc_spectra_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--cc-spectra");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (vrbflg)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        program_name);
    snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--verbose");
    snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }

  /* add the xml comment, if specified */
  if (!*comment)
  {
    snprintf(proc_table.processTable->comment, \
        LIGOMETA_COMMENT_MAX, " ");
  }
  else
  {
    snprintf(proc_table.processTable->comment, \
        LIGOMETA_COMMENT_MAX, "%s", comment);
  }

  /* delete empty first entry in process params table */
  this_proc_param = proc_params.processParamsTable;
  proc_params.processParamsTable = proc_params.processParamsTable->next;
  free(this_proc_param);

  /* set xml output file */
  if (output_path)
  {
    snprintf(xml_file_name, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml",
        output_path, base_name);
  }
  else
  {
    snprintf(xml_file_name, FILENAME_MAX * sizeof(CHAR), "%s.xml",
        base_name);
  }

  /* write out xml */
  if (vrbflg)
    fprintf(stdout, "Writing output XML file...\n");

  /* opening xml file stream */
  memset(&xml_stream, 0, sizeof(LIGOLwXMLStream));
  LAL_CALL(LALOpenLIGOLwXMLFile(status, &xml_stream, xml_file_name), status);

  /* write out process and process params tables */
  XLALGPSTimeNow(&proc_table.processTable->end_time);
  LAL_CALL(LALBeginLIGOLwXMLTable(status, &xml_stream, process_table), \
      status);
  LAL_CALL(LALWriteLIGOLwXMLTable(status, &xml_stream, proc_table, \
        process_table), status);
  LAL_CALL(LALEndLIGOLwXMLTable(status, &xml_stream), status);
  free(proc_table.processTable);

  /* write the process params table */
  LAL_CALL(LALBeginLIGOLwXMLTable(status, &xml_stream, \
        process_params_table), status);
  LAL_CALL(LALWriteLIGOLwXMLTable(status, &xml_stream, proc_params, \
        process_params_table), status);
  LAL_CALL(LALEndLIGOLwXMLTable(status, &xml_stream), status);

  /* write stochastic table */
  if (stoch_table)
  {
    output_table.stochasticTable = stoch_table;
    LAL_CALL(LALBeginLIGOLwXMLTable(status, &xml_stream, stochastic_table), \
        status);
    LAL_CALL(LALWriteLIGOLwXMLTable(status, &xml_stream, output_table, \
          stochastic_table), status);
    LAL_CALL(LALEndLIGOLwXMLTable(status, &xml_stream), status);
  }

  /* close xml file */
  LAL_CALL(LALCloseLIGOLwXMLFile(status, &xml_stream), status);

  return;
}

/*
 * vim: et
 */
