/*
 * Copyright (C) 2007 Badri Krishnan, Lucia Santamaria Lara, Robert Adam Mercer, Stephen Fairhurst
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


/**
 * \file ninja.c
 * \author Badri Krishnan
 * \brief Code for parsing and selecting numerical relativity
 *        waves in frame files
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include <lalapps.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Inject.h>
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>
#include <lal/UserInput.h>
#include <lalappsfrutils.h>
#include <lal/FrameStream.h>
#include <lal/LogPrintf.h>
#include <lal/LALFrameL.h>

#include <LALAppsVCSInfo.h>

#include <processtable.h>
#include "inspiral.h"

/* program info */
#define PROGRAM_NAME "lalapps_ninja"

/* verbose flag */
extern int vrbflg;

/*
 * structure definitions
 */

typedef struct {
  REAL8 massRatioMin;
  REAL8 massRatioMax;
  REAL8 sx1Min;
  REAL8 sx1Max;
  REAL8 sx2Min;
  REAL8 sx2Max;
  REAL8 sy1Min;
  REAL8 sy1Max;
  REAL8 sy2Min;
  REAL8 sy2Max;
  REAL8 sz1Min;
  REAL8 sz1Max;
  REAL8 sz2Min;
  REAL8 sz2Max;
  REAL8 freqStart22min;
  INT4 numGroups;
  NumRelGroup *grouplist;
} NrParRange;

typedef struct
{
  REAL8 massRatio; /**< Mass ratio m1/m2 where we assume m1 >= m2*/
  REAL8 spin1[3];  /**< Spin of m1 */
  REAL8 spin2[3];  /**< Spin of m2 */
  INT4  mode[2];   /**< l and m values */
  REAL8 freqStart22; /**< start frequency of 22 mode */
  CHAR  filename[LALNameLength]; /**< filename where data is stored */
  NumRelGroup group;
} NinjaMetaData;

/*
 * local helper function prototypes
 */

static int get_nr_metadata_from_framehistory(NinjaMetaData *data, FrHistory *history, CHAR *metadata_format);
static int get_mode_index_from_channel_name(INT4  *mode_l, INT4  *mode_m, CHAR  *name);
static int get_minmax_modes(INT4 *min, INT4 *max, FrameH *frame);
static int get_metadata_from_string(NinjaMetaData *data, CHAR *comment, CHAR *metadata_format);
static int metadata_in_range(NinjaMetaData *data, NrParRange *range);
static int parse_group_list ( NrParRange *range, CHAR *list);

/* main program entry */
int main(INT4 argc, CHAR *argv[])
{
  LALStatus status = blank_status;

  /* frame file stuff */
  LALCache *frGlobCache = NULL;
  LALCache *frInCache = NULL;
  FrameH *frame = NULL;
  FrFile *frFile = NULL;

  /* inspiral table stuff */
  SimInspiralTable *this_inj = NULL;
  LIGOLwXMLStream xmlfp;
  MetadataTable injections;
  MetadataTable proctable;
  //MetadataTable procparams;
  //ProcessParamsTable *this_proc_param = NULL;

  /* nrwave stuff */
  NinjaMetaData metaData;
  NrParRange range;

  /* counter */
  UINT4 k;

  /* user input variables */
  BOOLEAN uvar_help = 0;
  CHAR *uvar_nrDir = NULL;
  CHAR *uvar_pattern = NULL;
  CHAR *uvar_nrGroup = NULL;
  CHAR *uvar_outFile = NULL;
  CHAR *uvar_format = NULL;
  REAL8 uvar_minMassRatio = 1, uvar_maxMassRatio = 0;
  REAL8 uvar_minSx1 = -1, uvar_minSx2 = -1, uvar_maxSx1 = 1, uvar_maxSx2 = 1;
  REAL8 uvar_minSy1 = -1, uvar_minSy2 = -1, uvar_maxSy1 = 1, uvar_maxSy2 = 1;
  REAL8 uvar_minSz1 = -1, uvar_minSz2 = -1, uvar_maxSz1 = 1, uvar_maxSz2 = 1;
  REAL8 uvar_freqLo = 40;
  INT4 uvar_minMode = 2, uvar_maxMode = 2;

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;

  /* set default output file */
  uvar_outFile = (CHAR *)LALCalloc(1, FILENAME_MAX * sizeof(CHAR));
  strcpy(uvar_outFile, "ninja_out.xml");

  /* set default metadata format */
  uvar_format = (CHAR *)LALCalloc(1, 256 * sizeof(CHAR));
  strcpy(uvar_format, "NINJA1");

  LAL_CALL(LALRegisterBOOLUserVar(&status, "help", 'h', UVAR_HELP, "Print this message", &uvar_help), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "datadir", 'D', UVAR_REQUIRED, "Directory with NR data", &uvar_nrDir), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "pattern", 0, UVAR_OPTIONAL, "Filename pattern", &uvar_pattern), &status);

  LAL_CALL(LALRegisterSTRINGUserVar(&status, "outfile", 'o', UVAR_OPTIONAL, "Output xml filename", &uvar_outFile), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "format", 0, UVAR_OPTIONAL, "Metadata format", &uvar_format), &status);

  LAL_CALL(LALRegisterREALUserVar(&status, "min-mass-ratio", 0, UVAR_OPTIONAL, "Min. mass ratio", &uvar_minMassRatio), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-mass-ratio", 0, UVAR_OPTIONAL, "Max. mass ratio", &uvar_maxMassRatio), &status);

  LAL_CALL(LALRegisterREALUserVar(&status, "min-sx1", 0, UVAR_OPTIONAL, "Min. x-spin of first BH", &uvar_minSx1), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "min-sx2", 0, UVAR_OPTIONAL, "Min. x-Spin of second BH", &uvar_minSx2), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-sx1", 0, UVAR_OPTIONAL, "Max. x-spin of first BH", &uvar_maxSx1), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-sx2", 0, UVAR_OPTIONAL, "Max. x-spin of second BH", &uvar_maxSx2), &status);

  LAL_CALL(LALRegisterREALUserVar(&status, "min-sy1", 0, UVAR_OPTIONAL, "Min. y-spin of first BH", &uvar_minSy1), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "min-sy2", 0, UVAR_OPTIONAL, "Min. y-Spin of second BH", &uvar_minSy2), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-sy1", 0, UVAR_OPTIONAL, "Max. y-spin of first BH", &uvar_maxSy1), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-sy2", 0, UVAR_OPTIONAL, "Max. y-spin of second BH", &uvar_maxSy2), &status);

  LAL_CALL(LALRegisterREALUserVar(&status, "min-sz1", 0, UVAR_OPTIONAL, "Min. z-spin of first BH", &uvar_minSz1), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "min-sz2", 0, UVAR_OPTIONAL, "Min. z-Spin of second BH", &uvar_minSz2), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-sz1", 0, UVAR_OPTIONAL, "Max. z-spin of first BH", &uvar_maxSz1), &status);
  LAL_CALL(LALRegisterREALUserVar(&status, "max-sz2", 0, UVAR_OPTIONAL, "Max. z-spin of second BH", &uvar_maxSz2), &status);

  LAL_CALL(LALRegisterREALUserVar(&status, "freq-lo", 0, UVAR_OPTIONAL, "Lower cutoff frequency", &uvar_freqLo), &status);
  LAL_CALL(LALRegisterINTUserVar(&status, "min-mode", 0, UVAR_OPTIONAL, "Min mode value to be injected", &uvar_minMode), &status);
  LAL_CALL(LALRegisterINTUserVar(&status, "max-mode", 0, UVAR_OPTIONAL, "Max mode value to be injected", &uvar_maxMode), &status);

  LAL_CALL(LALRegisterSTRINGUserVar(&status, "nr-group", 0, UVAR_OPTIONAL, "NR group list (default=all)", &uvar_nrGroup), &status);

  /* read all command line variables */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0);

  /* check for supported metadata format */
  if (strcmp(uvar_format, "NINJA1") == 0);
  else if (strcmp(uvar_format, "NINJA2") == 0);
  else
  {
    fprintf(stderr, "Supported metadata formats are NINJA1 and NINJA2 (%s specified)\n", uvar_format);
    exit(1);
  }

  range.massRatioMin = uvar_minMassRatio;
  range.massRatioMax = uvar_maxMassRatio;

  range.sx1Min = uvar_minSx1;
  range.sx1Max = uvar_maxSx1;

  range.sx2Min = uvar_minSx2;
  range.sx2Max = uvar_maxSx2;

  range.sy1Min = uvar_minSy1;
  range.sy1Max = uvar_maxSy1;

  range.sy2Min = uvar_minSy2;
  range.sy2Max = uvar_maxSy2;

  range.sz1Min = uvar_minSz1;
  range.sz1Max = uvar_maxSz1;

  range.sz2Min = uvar_minSz2;
  range.sz2Max = uvar_maxSz2;

  parse_group_list(&range, uvar_nrGroup);

  LogPrintf(LOG_NORMAL, "Globbing frame files...");

  /* create a frame cache by globbing *.gwf in specified dir */
  frGlobCache = XLALCacheGlob(uvar_nrDir, uvar_pattern);

  frInCache = XLALCacheDuplicate(frGlobCache);

  XLALDestroyCache(frGlobCache);

  /* check we globbed at least one frame file */
  if (!frInCache->length)
  {
    fprintf(stderr, "error: no numrel frame files found\n");
    exit(1);
  }
  LogPrintfVerbatim(LOG_NORMAL, "found %d\n",frInCache->length);

  /* initialize head of simInspiralTable linked list to null */
  injections.simInspiralTable = NULL;

  LogPrintf(LOG_NORMAL, "Selecting frame files with right numrel parameters...");

  /* loop over frame files and select the ones with nr-params in the right range */
  for (k = 0; k < frInCache->length; k++)
  {
    frFile = XLALFrOpenURL(frInCache->list[k].url);

    frame = FrameRead(frFile);

    memset(&metaData, 0, sizeof(NinjaMetaData));
    get_nr_metadata_from_framehistory(&metaData, frame->history, uvar_format);

    /* if we find parameters in range then write to the siminspiral table */
    if (metadata_in_range(&metaData, &range))
    {
      REAL8 tmp;
      INT4 minMode, maxMode;

      /* alloc next element of inspiral table linked list */
      if (injections.simInspiralTable)
        this_inj = this_inj->next = (SimInspiralTable *)LALCalloc(1, sizeof(SimInspiralTable));
      else
        injections.simInspiralTable = this_inj = (SimInspiralTable *)LALCalloc(1, sizeof(SimInspiralTable));

      get_minmax_modes(&minMode,&maxMode,frame);

      /* eta = 1/(sqrt(mu) + 1/sqrt(mu))^2 where mu = m1/m2 */
      tmp = sqrt(metaData.massRatio) + (1.0 / sqrt(metaData.massRatio));
      this_inj->eta = 1.0 / (tmp * tmp);

      this_inj->spin1x = metaData.spin1[0];
      this_inj->spin1y = metaData.spin1[1];
      this_inj->spin1z = metaData.spin1[2];

      this_inj->spin2x = metaData.spin2[0];
      this_inj->spin2y = metaData.spin2[1];
      this_inj->spin2z = metaData.spin2[2];
      this_inj->f_lower = metaData.freqStart22;

      strcpy(this_inj->numrel_data, frInCache->list[k].url);

      this_inj->numrel_mode_min = uvar_minMode;
      this_inj->numrel_mode_max = uvar_maxMode;

    } /* end if (metadata is in range) */

    FrFileIEnd(frFile);

  } /* end loop over framefiles */
  LogPrintfVerbatim(LOG_NORMAL, "done\n");

  /* now write the output xml file */
  LogPrintf(LOG_NORMAL, "Writing xml output...");

  /* first the process table */
  proctable.processTable = (ProcessTable *)LALCalloc(1, sizeof(ProcessTable));
  XLALGPSTimeNow(&(proctable.processTable->start_time));

  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  snprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " ");

  memset(&xmlfp, 0, sizeof(LIGOLwXMLStream));
  LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlfp, uvar_outFile), &status);

  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, process_table), &status);

  LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, proctable, process_table), &status);
  LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);

#if 0
  /* now the process params table */
  LAL_CALL(LALUserVarGetProcParamsTable (&status, &this_proc_param, PROGRAM_NAME), &status);
  procparams.processParamsTable = this_proc_param;

  if (procparams.processParamsTable)
  {
    LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, process_params_table),
              &status);
    LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, procparams,
                                      process_params_table), &status);
    LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);
  }
#endif


  /* and finally the simInspiralTable itself */
  if (injections.simInspiralTable)
  {
    LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, sim_inspiral_table),
              &status);
    LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, injections,
                                      sim_inspiral_table), &status);
    LAL_CALL(LALEndLIGOLwXMLTable (&status, &xmlfp), &status);
  }
  LogPrintfVerbatim (LOG_NORMAL, "done\n");


  /* we are now done with the xml writing stuff */
  /* free memory and exit */

  LogPrintf(LOG_NORMAL, "Free memory and exiting...");

  /* close the various xml tables */
  while (injections.simInspiralTable)
  {
    this_inj = injections.simInspiralTable;
    injections.simInspiralTable = injections.simInspiralTable->next;
    LALFree(this_inj);
  }

#if 0
  while (procparams.processParamsTable)
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    LALFree(this_proc_param);
  }

  LALFree(proctable.processTable);
#endif

  /* close cache */
  /* LAL_CALL(LALFrClose(&status, &frStream), &status); */
  XLALDestroyCache(frInCache);

  /* close the injection file */
  LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlfp), &status);

  /* destroy all user input variables */
  if (range.grouplist != NULL)
    LALFree(range.grouplist);
  LAL_CALL(LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();
  LogPrintfVerbatim(LOG_NORMAL, "bye\n");

  return 0;
} /* main */


/* metadata is stored in the history field comment
   -- this function parses the comment to fill the metadata struct */
static int get_nr_metadata_from_framehistory(NinjaMetaData *data,
                                      FrHistory *history,
                                      CHAR *metadata_format)
{
  UINT4 stringlen = 128;
  CHAR *comment = NULL; /* the comments string */
  FrHistory *localhist;

  comment = LALMalloc(stringlen * sizeof(CHAR));

  localhist = history;
  while (localhist)
  {
    /* get history comment string and parse it   */
    /* The author-emails list can be > 128 chars */
    if (strlen(localhist->comment) + 1 > stringlen)
    {
      stringlen = strlen(localhist->comment) + 1;
      comment   = LALRealloc(comment, stringlen * sizeof(CHAR));
    }

    strcpy(comment,localhist->comment);
    get_metadata_from_string(data, comment, metadata_format);
    localhist = localhist->next;
  }

  LALFree(comment);
  return 0;
}


static int get_metadata_from_string(NinjaMetaData *data,
                             CHAR *comment,
                             CHAR *metadata_format)
{
  CHAR *token;
  CHAR *thiscomment = NULL;

  thiscomment = LALCalloc(1, (strlen(comment) + 1) * sizeof(CHAR));
  strcpy(thiscomment, comment);

  token = strtok(thiscomment, ":");

  if (strstr(token, "spin1x"))
  {
    token = strtok(NULL, ":");
    data->spin1[0] = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strstr(token, "spin1y"))
  {
    token = strtok(NULL, ":");
    data->spin1[1] = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strstr(token, "spin1z"))
  {
    token = strtok(NULL, ":");
    data->spin1[2] = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strstr(token, "spin2x"))
  {
    token = strtok(NULL, ":");
    data->spin2[0] = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strstr(token, "spin2y"))
  {
    token = strtok(NULL, ":");
    data->spin2[1] = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strstr(token, "spin2z"))
  {
    token = strtok(NULL, ":");
    data->spin2[2] = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strstr(token, "mass-ratio"))
  {
    token = strtok(NULL, ":");
    data->massRatio = atof(token);
    LALFree(thiscomment);
    return 0;
  }

  if (strcmp(metadata_format, "NINJA1") == 0)
  {
    if (strstr(token, "freqStart22"))
    {
      token = strtok(NULL, ":");
      data->freqStart22 = atof(token);
      LALFree(thiscomment);
      return 0;
    }
  }
  else if (strcmp(metadata_format, "NINJA2") == 0)
  {
    if (strstr(token, "freq_start_22"))
    {
      token = strtok(NULL, ":");
      data->freqStart22 = atof(token);
      LALFree(thiscomment);
      return 0;
    }
  }
  else
  {
    fprintf(stderr, "Supported metadata formats are NINJA1 and NINJA2 (%s specified)\n", metadata_format);
    exit(1);
  }

  if (strstr(token, "nr-group"))
  {
    token = strtok(NULL, ":");
    data->group = XLALParseNumRelGroupName(token);
    LALFree(thiscomment);
    return 0;
  }

  /* did not match anything */
  LALFree(thiscomment);
  return -1;

}


static int metadata_in_range(NinjaMetaData *data, NrParRange *range)
{

  INT4 ret, k;
  BOOLEAN flag = 0;
  BOOLEAN groupflag = 0;

  flag = (data->massRatio >= range->massRatioMin) && (data->massRatio <= range->massRatioMax);
  flag = flag && (data->spin1[0] >= range->sx1Min) && (data->spin1[0] <= range->sx1Max);
  flag = flag && (data->spin2[0] >= range->sx2Min) && (data->spin2[0] <= range->sx2Max);
  flag = flag && (data->spin1[1] >= range->sy1Min) && (data->spin1[1] <= range->sy1Max);
  flag = flag && (data->spin2[1] >= range->sy2Min) && (data->spin2[1] <= range->sy2Max);
  flag = flag && (data->spin1[2] >= range->sz1Min) && (data->spin1[2] <= range->sz1Max);
  flag = flag && (data->spin2[2] >= range->sz2Min) && (data->spin2[2] <= range->sz2Max);

  for (k = 0; k < range->numGroups; k++)
  {
    if (range->grouplist[k] == data->group)
      groupflag = 1;
  }

  /* if numgroups == 0 then user did not enter any groups and
     so we must select all groups */
  if (range->numGroups == 0)
    groupflag = 1;

  flag = flag && groupflag;

  if (flag)
    ret = 1;
  else
    ret = 0;

  return(ret);

}


static int get_minmax_modes(INT4 *min,
                     INT4 *max,
                     FrameH *frame)
{
  int ret = 1;
  INT4 mode_l = -1, mode_m, locmin, locmax;
  FrSimData *sim;

  locmin = 10;
  locmax = 0;
  sim = frame->simData;
  while (sim)
  {
    if (!get_mode_index_from_channel_name(&mode_l, &mode_m, sim->name))
    {
      if (locmin > mode_l)
        locmin = mode_l;
      if (locmax < mode_l)
        locmax = mode_l;
    }
    sim = sim->next;
  }

  *min = locmin;
  *max = locmax;

  return ret;
}

/* very hackish -- need to make this better */
static int get_mode_index_from_channel_name(INT4 *mode_l,
                                     INT4 *mode_m,
                                     CHAR *name)
{
  int ret = 1;
  CHAR *tmp;
  INT4 sign = 0;

  tmp = strstr(name, "hcross_");
  if (tmp)
  {
    tmp += strlen("hcross_") + 1;
    *mode_l = atoi(tmp);
    tmp = strstr(tmp, "m");
    tmp++;

    if (!strncmp(tmp, "p", 1))
      sign = 1;

    if (!strncmp(tmp, "n", 1))
      sign = -1;

    tmp++;
    *mode_m = sign * atoi(tmp);
    ret = 0;
  }

  tmp = strstr(name, "hplus_");
  if (tmp)
  {
    tmp += strlen("hplus_") + 1;
    *mode_l = atoi(tmp);
    tmp = strstr(tmp, "m");
    tmp++;

    if (!strncmp(tmp, "p", 1))
      sign = 1;

    if (!strncmp(tmp, "n", 1))
      sign = -1;

    tmp++;
    *mode_m = sign * atoi(tmp);

    ret = 0;
  }

  return ret;
}

/** take a list of numrel group names separated by ";" and parse it to
    get a vector of NumRelGroup */
static int parse_group_list(NrParRange *range,
                       CHAR *list)
{
  UINT4 numGroups = 0;
  NumRelGroup thisGroup;
  NumRelGroup *grouplist = NULL;
  CHAR *token;
  CHAR *thislist = NULL;

  /* look for the ";" token */
  if (list)
  {
    thislist = LALCalloc(1, (strlen(list) + 1) * sizeof(CHAR));
    strcpy(thislist, list);

    token = strtok(thislist, ";");

    while (token)
    {
      thisGroup = XLALParseNumRelGroupName(token);

      /* if the parsing was successful, add to list */
      if (thisGroup != NINJA_GROUP_LAST)
      {
        numGroups++;
        grouplist = LALRealloc(grouplist, numGroups * sizeof(*grouplist));
        grouplist[numGroups - 1] = thisGroup;
      }
      token = strtok(NULL, ";");
    }
    LALFree(thislist);
  } /* if(list) */

  /* copy to output */
  range->numGroups = numGroups;
  range->grouplist = grouplist;

  return numGroups;
}
