/*
 * helper.c - SGWB Standalone Analysis Pipeline
 *          - helper functions
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
#include <errno.h>
#include <getopt.h>

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/LALStdio.h>
#include <lal/PrintFTSeries.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SimulateSB.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>

#include <lalapps.h>
#include <processtable.h>

#include "stochastic.h"

NRCSID(HELPERC, "$Id$");
RCSID("$Id$");

/* external variables */
extern INT4 sampleRate;
extern INT4 resampleRate;
extern REAL4 geoHighPassFreq;
extern INT4  geoHighPassOrder;
extern REAL4 geoHighPassAtten;
extern int vrbflg;

/* read a time series */
REAL4TimeSeries *get_time_series(LALStatus *status,
    CHAR *ifo,
    CHAR *cacheFile,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 buffer)
{
  /* variables */
  REAL4TimeSeries *series;
  FrStream *stream = NULL;
  FrCache *frameCache = NULL;
  ResampleTSParams params;
  size_t length;

  /* calculate number of data points */
  length = delta_gps_to_float(status, end, start) * resampleRate;

  /* apply resample buffer */
  if (buffer)
  {
    start.gpsSeconds -= buffer;
    end.gpsSeconds += buffer;
  }

  /* set resample parameters */
  params.deltaT = 1./resampleRate;
  params.filterType = defaultButterworth;

  if (vrbflg)
    fprintf(stdout, "Opening frame cache \"%s\"...\n", cacheFile);

  /* open frame stream */
  LAL_CALL(LALFrCacheImport(status, &frameCache, cacheFile), status);
  LAL_CALL(LALFrCacheOpen(status, &stream, frameCache), status);
  LAL_CALL(LALDestroyFrCache(status, &frameCache), status);

  /* turn on checking for missing data */
  stream->mode = LAL_FR_VERBOSE_MODE;

  /* get the data */
  if (strncmp(ifo, "G1", 2) == 0)
    series = get_geo_data(status, stream, channel, start, end);
  else
    series = get_ligo_data(status, stream, channel, start, end);

  /* check for missing data */
  if (stream->state & LAL_FR_GAP)
  {
    fprintf(stderr, "Gap in data detected between GPS times %d s and %d s\n", \
        start.gpsSeconds, end.gpsSeconds);
    LAL_CALL(LALDestroyREAL4TimeSeries(status, series), status);
    exit(1);
  }

  /* clean up */
  LAL_CALL(LALFrClose(status, &stream), status);

  /* resample if required */
  if (resampleRate != sampleRate)
  {
    if (vrbflg)
      fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);

    LAL_CALL(LALResampleREAL4TimeSeries(status, series, &params), status);
  }

  /* remove resample buffer */
  if (buffer)
  {
    LAL_CALL(LALShrinkREAL4TimeSeries(status, series, \
          buffer * resampleRate, length), status);
  }

  return(series);
}

/* read a LIGO time series */
REAL4TimeSeries *get_ligo_data(LALStatus *status,
    FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  REAL4TimeSeries *series;
  FrChanIn channelIn;
  size_t length;

  /* set channelIn */
  channelIn.name = channel;
  channelIn.type = ADCDataChannel;

  if (vrbflg)
    fprintf(stderr, "Allocating memory for \"%s\" series...\n", channel);

  /* create and initialise time series */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, channel, start, 0, 0, \
        lalADCCountUnit, 0), status);

  if (vrbflg)
    fprintf(stderr, "Reading \"%s\" series metadata...\n", channel);

  /* get the series meta data */
  LAL_CALL(LALFrGetREAL4TimeSeriesMetadata(status, series, &channelIn, \
        stream), status);

  if (vrbflg)
    fprintf(stderr, "Resizing \"%s\" series...\n", channel);

  /* resize series to the correct number of samples */
  length = delta_gps_to_float(status, end, start) / series->deltaT;
  LAL_CALL(LALResizeREAL4TimeSeries(status, series, 0, length), status);

  if (vrbflg)
    fprintf(stdout, "Reading channel \"%s\"...\n", channel);

  /* seek to and read data */
  LAL_CALL(LALFrSeek(status, &start, stream), status);
  LAL_CALL(LALFrGetREAL4TimeSeries(status, series, &channelIn, stream), status);

  return(series);
}

/* read and high pass filter a GEO time series */
REAL4TimeSeries *get_geo_data(LALStatus *status,
    FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  PassBandParamStruc highPassParam;
  REAL4TimeSeries *series;
  REAL8TimeSeries *geo;
  FrChanIn channelIn;
  size_t length;
  size_t i;

  /* set channelIn */
  channelIn.name = channel;
  channelIn.type = ADCDataChannel;

  if (vrbflg)
    fprintf(stderr, "Allocating memory for \"%s\" series...\n", channel);

  /* create and initialise time series */
  LAL_CALL(LALCreateREAL8TimeSeries(status, &geo, channel, start, 0, 0, \
        lalADCCountUnit, 0), status);

  if (vrbflg)
    fprintf(stderr, "Reading \"%s\" series metadata...\n", channel);

  /* get the series meta data */
  LAL_CALL(LALFrGetREAL8TimeSeriesMetadata(status, geo, &channelIn, \
        stream), status);

  if (vrbflg)
    fprintf(stderr, "Resizing \"%s\" series...\n", channel);

  /* resize series to the correct number of samples */
  length = delta_gps_to_float(status, end, start) / series->deltaT;
  LAL_CALL(LALResizeREAL8TimeSeries(status, geo, 0, length), status);

  if (vrbflg)
    fprintf(stdout, "Reading channel \"%s\"...\n", channel);

  /* seek to and read data */
  LAL_CALL(LALFrSeek(status, &start, stream), status);
  LAL_CALL(LALFrGetREAL8TimeSeries(status, geo, &channelIn, stream), status);

  if (vrbflg)
    fprintf(stdout, "High pass filtering \"%s\"...\n", channel);

  /* high pass filter before casting to a REAL4 */
  highPassParam.nMax = geoHighPassOrder;
  highPassParam.f1 = -1;
  highPassParam.f2 = geoHighPassFreq;
  highPassParam.a1 = -1;
  highPassParam.a2 = geoHighPassAtten;
  LAL_CALL(LALButterworthREAL8TimeSeries(status, geo, &highPassParam), status);

  if (vrbflg)
    fprintf(stdout, "Casting \"%s\" as a REAL4...\n", channel);

  /* cast as a REAL4 */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, geo->name, geo->epoch, \
        geo->f0, geo->deltaT, geo->sampleUnits, geo->data->length), status);
  for (i = 0; i < series->data->length; i++)
    series->data->data[i] = (REAL4)geo->data->data[i];

  /* destroy geo series */
  LAL_CALL(LALDestroyREAL8TimeSeries(status, geo), status);

  return(series);
}

/* return the difference between two GPS times as REAL8 */
REAL8 delta_gps_to_float(LALStatus *status,
    LIGOTimeGPS end,
    LIGOTimeGPS start)
{
  /* varaibles */
  LALTimeInterval i;
  REAL8 d;

  /* get the time difference as a REAL8 */
  LAL_CALL(LALDeltaGPS(status, &i, &end, &start), status);
  LAL_CALL(LALIntervalToFloat(status, &d, &i), status);

  return(d);
}

/* wrapper function to return the spectrum */
REAL4FrequencySeries *omega_gw(LALStatus *status,
    REAL4 alpha,
    REAL8 fRef,
    REAL4 omegaRef,
    UINT4 length,
    REAL8 f0,
    REAL8 deltaF,
    LIGOTimeGPS time)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOmegaGWParameters params;

  /* create and initialise frequency series */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "OmegaGW", time, f0, \
        deltaF, lalDimensionlessUnit, length), status);

  /* set parameters */
  params.alpha = alpha;
  params.fRef = fRef;
  params.omegaRef = omegaRef;
  params.length = length;
  params.f0 = f0;
  params.deltaF = deltaF;

  /* calculate spectrum */
  LAL_CALL(LALStochasticOmegaGW(status, series, &params), status);

  return(series);
}

/* wrapper function to return overlap reduction function */
REAL4FrequencySeries *overlap_reduction_function(LALStatus *status,
    UINT4 length,
    REAL8 f0,
    REAL8 deltaF,
    INT4 siteOne,
    INT4 siteTwo,
    LIGOTimeGPS time)
{
  /* variables */
  REAL4FrequencySeries *series;
  OverlapReductionFunctionParameters params;
  LALDetectorPair detectors;

  /* create and initialise frequency series */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "Overlap", time, f0, \
        deltaF, lalDimensionlessUnit, length), status);

  /* set parameters */
  params.length = length;
  params.f0 = f0;
  params.deltaF = deltaF;

  /* set detectors */
  detectors.detectorOne = lalCachedDetectors[siteOne];
  detectors.detectorTwo = lalCachedDetectors[siteTwo];

  /* calculate overlap reduction function */
  LAL_CALL(LALOverlapReductionFunction(status, series, &detectors, &params), \
      status);

  return(series);
}

/* helper function to increment the seconds component of a gps time with
 * an integer */
LIGOTimeGPS increment_gps(LALStatus *status,
    LIGOTimeGPS time,
    INT4 increment)
{
  /* variables */
  LIGOTimeGPS result;
  LALTimeInterval interval;

  /* set interval */
  interval.seconds = increment;
  interval.nanoSeconds = 0;

  /* increment GPS time */
  LAL_CALL(LALIncrementGPS(status, &result, &time, &interval), status);

  return(result);
}

/* function to cut a time series between given start and end times */
REAL4TimeSeries *cut_time_series(LALStatus *status,
    REAL4TimeSeries *input,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  REAL4TimeSeries *series;
  INT4 length;
  INT4 first;

  /* calculate length of segment to cut */
  length = (INT4)delta_gps_to_float(status, end, start) / input->deltaT;

  /* get first bin */
  first = (INT4)((start.gpsSeconds - input->epoch.gpsSeconds) / input->deltaT);
  
  /* allocate memory */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, input->name, start, \
        input->f0, input->deltaT, input->sampleUnits, length), status);

  /* cut time series */
  LAL_CALL(LALCutREAL4TimeSeries(status, &series, input, first, length), \
      status);

  return(series);
}

/* function to save out ccSpectra as a frame file */
void write_ccspectra_frame(COMPLEX8FrequencySeries *series,
    CHAR *ifoOne,
    CHAR *ifoTwo,
    LIGOTimeGPS time,
    INT4 duration)
{
  /* variables */
  CHAR hertz[] = "Hz";
  CHAR comment[] = "$Id$";
  CHAR source[FILENAME_MAX];
  CHAR fname[FILENAME_MAX];
  CHAR units[LALUnitNameSize];
  struct FrFile *frfile;
  struct FrameH *frame;
  struct FrVect *vect;
  struct FrProcData *proc;

  /* set frame filename */
  LALSnprintf(source, sizeof(source), "%s%s", ifoOne, ifoTwo);
  LALSnprintf(fname, sizeof(fname), "%s-ccSpectra-%d-%d.gwf", source, \
      time.gpsSeconds, duration);

  /* setup frame file */
  frfile = FrFileONew(fname, 0);

  /* set frame properties */
  frame = FrameHNew(source);
  frame->run = 0;
  frame->frame = 0;
  frame->GTimeS = time.gpsSeconds;
  frame->GTimeN = time.gpsNanoSeconds;
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
  FrStrCpy(&proc->name, "CCSPECTRA");
  FrStrCpy(&proc->comment, comment);
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

/*
 * vim: et
 */
