/*
 * stochastic_dev.c - SGWB Standalone Analysis Pipeline
 *                  - Development Branch
 *
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <math.h>
#include <getopt.h>

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/LALStdio.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/SimulateSB.h>

#include <lalapps.h>
#include <processtable.h>

NRCSID(STOCHASTICC, "$Id$");
RCSID("$Id$");

/* cvs info */
#define PROGRAM_NAME "stochastic_dev"
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define CVS_SOURCE "$Source$"

/* xml process param table helper */
#define ADD_PROCESS_PARAM(pptype, format, ppvalue) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
    calloc(1, sizeof(ProcessParamsTable)); \
  LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME); \
  LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
      long_options[option_index].name); \
  LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype); \
  LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue);

/* window duration for PSD estimation */
#define PSD_WINDOW_DURATION 4

/* flags for getopt_long */
static int middle_segment_flag;
static int apply_mask_flag;
static int high_pass_flag;
static int overlap_hann_flag;
static int recentre_flag;
static int cc_spectra_flag;
static int debug_flag;
extern int vrbflg;

/* xml comment/tags */
CHAR comment[LIGOMETA_COMMENT_MAX];
CHAR *userTag = NULL;

/* xml tables */
MetadataTable proctable;
MetadataTable procparams;
ProcessParamsTable *this_proc_param;

/* parameters for the stochastic search */

/* sampling parameters */
INT4 resampleRate;

/* data parameters */
INT4 startTime = 0;
INT4 endTime = 0;
INT4 totalDuration = -1;
INT4 intervalDuration = -1;
INT4 segmentDuration = -1;
INT4 calibOffset = -1;
CHAR *frameCacheOne = NULL;
CHAR *frameCacheTwo = NULL;
CHAR *calCacheOne = NULL;
CHAR *calCacheTwo = NULL;
CHAR *channelOne = NULL;
CHAR *channelTwo = NULL;
CHAR *ifoOne = NULL;
CHAR *ifoTwo = NULL;
INT4 siteOne;
INT4 siteTwo;

/* frequency band */
REAL8 fMin = -1;
REAL8 fMax = -1;

/* omegaGW parameters */
REAL4 alpha = 0;
REAL4 fRef = 100;
REAL4 omegaRef = 1;

/* window parameters */
INT4 hannDuration = -1;

/* high pass filtering parameters */
REAL4 highPassFreq = -1;
REAL4 highPassAtten = -1;
INT4  highPassOrder = -1;

/* GEO scale factor */
REAL4 geoScaleFactor = 1e18;

/* GEO high pass filter parameters */
REAL4 geoHighPassFreq = -1;
INT4  geoHighPassOrder = -1;
REAL4 geoHighPassAtten = -1;

/* number of bins for frequency masking */
INT4 maskBin = -1;

/* output file */
CHAR *outputPath = NULL;

/* helper functions */

/* can't use round() as its not C89 */
static double myround(double x)
{
  return(x < 0 ? ceil(x - 0.5): floor(x + 0.5));
}

/* read a LIGO time series */
static REAL4TimeSeries *get_ligo_data(LALStatus *status,
    FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  REAL4TimeSeries *series;
  FrChanIn channel_in;
  size_t length;

  /* set channelIn */
  channel_in.name = channel;
  channel_in.type = ADCDataChannel;

  if (vrbflg)
    fprintf(stderr, "Allocating memory for \"%s\" series...\n", channel);

  /* create and initialise time series */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, channel, start, 0, 0, \
        lalADCCountUnit, 0), status);

  if (vrbflg)
    fprintf(stderr, "Reading \"%s\" series metadata...\n", channel);

  /* get the series meta data */
  XLALFrGetREAL4TimeSeriesMetadata(series, stream);

  if (vrbflg)
    fprintf(stderr, "Resizing \"%s\" series...\n", channel);

  /* resize series to the correct number of samples */
  length = floor((XLALDeltaFloatGPS(&end, &start) / series->deltaT) + 0.5);
  LAL_CALL(LALResizeREAL4TimeSeries(status, series, 0, length), status);

  if (vrbflg)
    fprintf(stdout, "Reading channel \"%s\"...\n", channel);

  /* seek to and read data */
  XLALFrSeek(stream, &start);
  XLALFrGetREAL4TimeSeries(series, stream);

  return(series);
}

/* read and high pass filter a GEO time series */
static REAL4TimeSeries *get_geo_data(LALStatus *status,
    FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  PassBandParamStruc high_pass_params;
  REAL4TimeSeries *series;
  REAL8TimeSeries *geo;
  FrChanIn channel_in;
  size_t length;
  size_t i;

  /* set channelIn */
  channel_in.name = channel;
  channel_in.type = ADCDataChannel;

  if (vrbflg)
    fprintf(stderr, "Allocating memory for \"%s\" series...\n", channel);

  /* create and initialise time series */
  LAL_CALL(LALCreateREAL8TimeSeries(status, &geo, channel, start, 0, 0, \
        lalADCCountUnit, 0), status);

  if (vrbflg)
    fprintf(stderr, "Reading \"%s\" series metadata...\n", channel);

  /* get the series meta data */
  XLALFrGetREAL8TimeSeriesMetadata(geo, stream);

  if (vrbflg)
    fprintf(stderr, "Resizing \"%s\" series...\n", channel);

  /* resize series to the correct number of samples */
  length = floor((XLALDeltaFloatGPS(&end, &start) / series->deltaT) + 0.5);
  LAL_CALL(LALResizeREAL8TimeSeries(status, geo, 0, length), status);

  if (vrbflg)
    fprintf(stdout, "Reading channel \"%s\"...\n", channel);

  /* seek to and read data */
  XLALFrSeek(stream, &start);
  XLALFrGetREAL8TimeSeries(geo, stream);

  if (vrbflg)
    fprintf(stdout, "High pass filtering \"%s\"...\n", channel);

  /* high pass filter before casting to a REAL4 */
  high_pass_params.nMax = geoHighPassOrder;
  high_pass_params.f1 = -1;
  high_pass_params.f2 = geoHighPassFreq;
  high_pass_params.a1 = -1;
  high_pass_params.a2 = geoHighPassAtten;
  XLALButterworthREAL8TimeSeries(geo, &high_pass_params);

  if (vrbflg)
    fprintf(stdout, "Casting \"%s\" as a REAL4...\n", channel);

  /* cast as a REAL4 */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, geo->name, geo->epoch, \
        geo->f0, geo->deltaT, geo->sampleUnits, geo->data->length), status);
  for (i = 0; i < series->data->length; i++)
    series->data->data[i] = (REAL4)geo->data->data[i];

  /* destroy geo series */
  XLALDestroyREAL8TimeSeries(geo);

  return(series);
}

/* read a time series */
static REAL4TimeSeries *get_time_series(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 buffer)
{
  /* variables */
  REAL4TimeSeries *series;
  FrStream *stream = NULL;
  FrCache *cache = NULL;
  size_t length;
  PassBandParamStruc high_pass_params;
  int mode = LAL_FR_VERBOSE_MODE;

  /* apply resample buffer if required */
  if (buffer)
  {
    start.gpsSeconds -= buffer;
    end.gpsSeconds += buffer;
  }

  if (vrbflg)
    fprintf(stdout, "Opening frame cache \"%s\"...\n", cache_file);

  /* open frame stream */
  cache = XLALFrImportCache(cache_file);
  stream = XLALFrCacheOpen(cache);
  XLALFrDestroyCache(cache);

  /* turn on checking for missing data */
  XLALFrSetMode(stream, mode);

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
    XLALDestroyREAL4TimeSeries(series);
    exit(1);
  }

  /* clean up */
  XLALFrClose(stream);

  /* resample if required */
  if (resampleRate)
  {
    if (vrbflg)
      fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);

    /* resample */
    XLALResampleREAL4TimeSeries(series, 1./resampleRate);
  }

  /* high pass fitering */
  if (high_pass_flag)
  {
    if (vrbflg)
    {
      fprintf(stdout, "Applying high pass filter to \"%s\"...\n", \
          series->name);
    }

    /* set high pass filter parameters */
    high_pass_params.nMax = highPassOrder;
    high_pass_params.f1 = -1;
    high_pass_params.f2 = highPassFreq;
    high_pass_params.a1 = -1;
    high_pass_params.a2 = highPassAtten;

    /* high pass filter */
    XLALButterworthREAL4TimeSeries(series, &high_pass_params);
  }

  /* remove resample buffer */
  if (buffer)
  {
    /* recover original start and end times */
    start.gpsSeconds += buffer;
    end.gpsSeconds -= buffer;

    /* calculate required length */
    length = floor((XLALDeltaFloatGPS(&end, &start) / series->deltaT) + 0.5);

    /* remove resample buffer */
    LAL_CALL(LALShrinkREAL4TimeSeries(status, series, \
          buffer / series->deltaT, length), status);
  }

  return(series);
}

/* wrapper function to return the spectrum */
static REAL4FrequencySeries *omega_gw(LALStatus *status,
    REAL4 exponent,
    REAL8 f_ref,
    REAL4 omega_ref,
    UINT4 length,
    REAL8 f0,
    REAL8 deltaF)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOmegaGWParameters omega_params;
  LIGOTimeGPS epoch;

  /* set epoch */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  /* create and initialise frequency series */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "OmegaGW", \
        epoch, f0, deltaF, lalDimensionlessUnit, length), status);

  /* set parameters */
  omega_params.alpha = exponent;
  omega_params.fRef = f_ref;
  omega_params.omegaRef = omega_ref;
  omega_params.length = length;
  omega_params.f0 = f0;
  omega_params.deltaF = deltaF;

  /* calculate spectrum */
  LAL_CALL(LALStochasticOmegaGW(status, series, &omega_params), status);

  return(series);
}

/* wrapper function to return overlap reduction function */
static REAL4FrequencySeries *overlap_reduction_function(LALStatus *status,
    UINT4 length,
    REAL8 f0,
    REAL8 deltaF,
    INT4 site_one,
    INT4 site_two)
{
  /* variables */
  REAL4FrequencySeries *series;
  OverlapReductionFunctionParameters overlap_params;
  LALDetectorPair detectors;
  LIGOTimeGPS epoch;

  /* set epoch */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  /* create and initialise frequency series */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "Overlap", \
        epoch, f0, deltaF, lalDimensionlessUnit, length), status);

  /* set parameters */
  overlap_params.length = length;
  overlap_params.f0 = f0;
  overlap_params.deltaF = deltaF;

  /* set detectors */
  detectors.detectorOne = lalCachedDetectors[site_one];
  detectors.detectorTwo = lalCachedDetectors[site_two];

  /* calculate overlap reduction function */
  LAL_CALL(LALOverlapReductionFunction(status, series, &detectors, \
        &overlap_params), status);

  return(series);
}

/* function to cut a time series between given start and end times */
static REAL4TimeSeries *cut_time_series(LALStatus *status,
    REAL4TimeSeries *input,
    LIGOTimeGPS start,
    UINT4 duration)
{
  /* variables */
  REAL4TimeSeries *series;
  INT4 length;
  INT4 first;

  /* calculate length of segment to cut */
  length = floor((duration / input->deltaT) + 0.5);

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
static void write_ccspectra_frame(COMPLEX8FrequencySeries *series,
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
  LALSnprintf(source, sizeof(source), "%s%s", ifo_one, ifo_two);
  LALSnprintf(fname, sizeof(fname), "%s-%s-%d-%d.gwf", source, \
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

/* function to return the data window */
static REAL4Window *data_window(REAL8 deltaT,
    INT4 length,
    INT4 hann_duration)
{
  /* variables */
  REAL4Window *window = NULL;
  REAL4Window *hann = NULL;
  INT4 hann_length;
  INT4 i;

  /* get length of hann segment requested */
  hann_length = (INT4)(hann_duration / deltaT);

  if (hann_length == 0)
  {
    /* rectangular window requested */
    window = XLALCreateRectangularREAL4Window(length);
  }
  else if (hann_length == length)
  {
    /* pure hann window requested */
    window = XLALCreateHannREAL4Window(length);
  }
  else if ((hann_length > 0) && (hann_length < length))
  {
    window = XLALCreateRectangularREAL4Window(length);
    hann =  XLALCreateHannREAL4Window(hann_length);

    /* construct tukey window */
    for (i = 0; i < hann_length / 2; i++)
      window->data->data[i] = hann->data->data[i];
    for (i = hann_length / 2; i < hann_length; i++)
      window->data->data[length - hann_length + i] = hann->data->data[i];

    /* free memory for hann window */
    XLALDestroyREAL4Window(hann);
  }
  else
  {
    fprintf(stderr, "Invalid hann_length to data_window()...\n");
    exit(1);
  }

  return(window);
}

/* wrapper function for calculating the inverse noise */
static REAL4FrequencySeries *inverse_noise(LALStatus *status,
    REAL4FrequencySeries *psd,
    COMPLEX8FrequencySeries *response)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticInverseNoiseInput input;
  StochasticInverseNoiseCalOutput output;

  /* allocate memory */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "calPSD", \
        response->epoch, response->f0, response->deltaF, \
        lalDimensionlessUnit, response->data->length), status);

  /* set input */
  input.unCalibratedNoisePSD = psd;
  input.responseFunction = response;

  /* set output */
  output.calibratedInverseNoisePSD = series;

  /* generate inverse noise */
  LAL_CALL(LALStochasticInverseNoiseCal(status, &output, &input), status);

  return(series);
}

/* wrapper function for calculating optimal filter */
static REAL4FrequencySeries *optimal_filter(LALStatus *status,
    REAL4FrequencySeries *overlap,
    REAL4FrequencySeries *omega,
    REAL4FrequencySeries *psd_one,
    REAL4FrequencySeries *psd_two,
    REAL4Window *window,
    REAL8 *sigma)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOptimalFilterNormalizationInput norm_input;
  StochasticOptimalFilterNormalizationOutput norm_output;
  StochasticOptimalFilterNormalizationParameters norm_params;
  StochasticOptimalFilterCalInput input;
  REAL4WithUnits normalisation;
  REAL4WithUnits variance;

  /* set parameters for normalisation */
  norm_params.fRef = fRef;
  norm_params.heterodyned = 0;
  norm_params.window1 = window->data;
  norm_params.window2 = window->data;

  /* set inputs for normalisation */
  norm_input.overlapReductionFunction = overlap;
  norm_input.omegaGW = omega;
  norm_input.inverseNoisePSD1 = psd_one;
  norm_input.inverseNoisePSD2 = psd_two;

  /* set normalisation output */
  norm_output.normalization = &normalisation;
  norm_output.variance = &variance;

  /* calculate variance and normalisation for the optimal filter */
  LAL_CALL(LALStochasticOptimalFilterNormalization(status, \
        &norm_output, &norm_input, &norm_params), status);

  /* get theoretical sigma */
  *sigma = sqrt((REAL8)(segmentDuration * variance.value * \
        pow(10, variance.units.powerOfTen)));

  /* allocate memory */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "filter", \
        psd_one->epoch, psd_one->f0, psd_one->deltaF, lalDimensionlessUnit, \
        psd_one->data->length), status);

  /* set input */
  input.overlapReductionFunction = overlap;
  input.omegaGW = omega;
  input.calibratedInverseNoisePSD1 = psd_one;
  input.calibratedInverseNoisePSD2 = psd_two;

  /* generate optimal filter */
  LAL_CALL(LALStochasticOptimalFilterCal(status, series, &input, \
        &normalisation), status);

  return(series);
}

/* wrapper function for estimating the PSD */
static REAL4FrequencySeries *estimate_psd(LALStatus *status,
    REAL4TimeSeries *series,
    REAL8 f0,
    INT4 shrink_length)
{
  /* variables */
  REAL4FrequencySeries *psd;
  REAL4Window *window;
  RealFFTPlan *plan = NULL;
  UINT4 length;
  UINT4 overlap;
  REAL8 deltaF;
  UINT4 psd_length;

  /* lengths */
  length = PSD_WINDOW_DURATION / series->deltaT;
  overlap = length / 2;
  deltaF = 1./(REAL8)PSD_WINDOW_DURATION;
  psd_length = (length / 2) + 1;

  /* allocate memory */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &psd, "psd", \
        series->epoch, series->f0, deltaF, lalDimensionlessUnit, \
        psd_length), status);

  /* create window for PSD estimation */
  window = XLALCreateHannREAL4Window(length);

  /* create fft plan for PSD estimation */
  plan = XLALCreateForwardREAL4FFTPlan(length, 0);

  /* esimate PSD */
  XLALREAL4AverageSpectrumWelch(psd, series, length, overlap, window, plan);

  /* destroy fft plan */
  XLALDestroyREAL4FFTPlan(plan);

  /* free memory for window */
  XLALDestroyREAL4Window(window);

  /* reduce to relevant frequency range */
  LAL_CALL(LALShrinkREAL4FrequencySeries(status, psd, (INT4)(f0/deltaF), \
        shrink_length), status);

  return(psd);
}

/* function to return a unity response function, for use with GEO data */
static COMPLEX8FrequencySeries *unity_response(LALStatus *status,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 deltaF,
    LALUnit units,
    INT4 length)
{
  /* variables */
  COMPLEX8FrequencySeries *response;
  int i;

  /* allocate memory */
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(status, &response, "response", \
        epoch, f0, deltaF, units, length), status);

  /* get unity response function */
  for (i = 0; i < length; i++)
  {
    response->data->data[i].re = 1;
    response->data->data[i].im = 0;
  }

  return(response);
}

/* wrapper function for generating response function for LIGO data */
static COMPLEX8FrequencySeries *ligo_response(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 deltaF,
    LALUnit units,
    INT4 length,
    INT4 offset)
{
  /* variables */
  COMPLEX8FrequencySeries *response;
  FrCache *cache = NULL;
  CalibrationUpdateParams calib_params;

  /* apply offset to epoch */
  epoch.gpsSeconds += offset;

  /* allocate memory */
  memset(&calib_params, 0, sizeof(CalibrationUpdateParams));
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(status, &response, "response", \
        epoch, 0, deltaF, units, length + (INT4)(f0/deltaF)), status);

  /* set ifo */
  calib_params.ifo = ifo;

  /* open calibration frame cache */
  cache = XLALFrImportCache(cache_file);

  /* generate response function */
  LAL_CALL(LALExtractFrameResponse(status, response, cache, &calib_params), \
      status);

  /* destory calibration frame cache */
  XLALFrDestroyCache(cache);

  /* reduce to required band */
  LAL_CALL(LALShrinkCOMPLEX8FrequencySeries(status, response, f0/deltaF, \
        length), status);

  return(response);
}

/* wrapper to unity_response and ligo_response for generating the
 * appropriate response for the given detector */
static COMPLEX8FrequencySeries *generate_response(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 deltaF,
    LALUnit units,
    INT4 length,
    INT4 offset)
{
  /* variables */
  COMPLEX8FrequencySeries *response = NULL;

  if (strncmp(ifo, "G1", 2) == 0)
  {
    /* generate response for GEO */
    response = unity_response(status, epoch, f0, deltaF, units, length);
  }
  else
  {
    /* generate response function for LIGO */
    response = ligo_response(status, ifo, cache_file, epoch, f0, deltaF, \
        units, length, offset);
  }

  return(response);
}

/* helper function to return the frequency mask */
static REAL4FrequencySeries *frequency_mask(LALStatus *status,
    REAL8 f0,
    REAL8 deltaF,
    INT4 length,
    INT4 bins)
{
  /* counters */
  INT4 i, j;

  /* variables */
  REAL4FrequencySeries *mask;
  LIGOTimeGPS epoch;
  INT4 nBins;
  INT4 numFMin;
  INT4 full_length;

  /* initialise time */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  /* extra bins */
  nBins = (bins - 1) / 2;
  numFMin = (INT4)(f0 / deltaF);

  /* get length for full frequency band */
  full_length = length + numFMin;

  /* allocate memory for frequency mask */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &mask, \
        "mask", epoch, 0, deltaF, lalDimensionlessUnit, \
        full_length), status);

  /* set all values to 1 */
  for (i = 0; i < full_length; i++)
    mask->data->data[i] = 1;

  /* remove multiples of 16 Hz */
  for (i = 0; i < full_length; i += (INT4)(16 / deltaF))
  {
    mask->data->data[i]= 0;

    for (j = 0; j < nBins; j++)
    {
      if ((i + 1 + j) < full_length)
        mask->data->data[i + 1 + j]= 0;
      if ((i - 1 - j) > 0 )
        mask->data->data[i - 1 - j]= 0;
    }
  }

  /* remove multiples of 60 Hz */
  for (i = 0; i < full_length; i += (INT4)(60 / deltaF))
  {
    mask->data->data[i] = 0;

    for (j = 0; j < nBins; j++)
    {
      if ((i + 1 + j) < full_length)
        mask->data->data[i + 1 + j]= 0;
      if ((i - 1 - j) > 0 )
        mask->data->data[i - 1 - j]= 0;
    }
  }

  /* get appropriate band */
  LAL_CALL(LALShrinkREAL4FrequencySeries(status, mask, numFMin, \
        length), status);

  return(mask);
}

/* wrapper function for performing zero pad and fft */
static COMPLEX8FrequencySeries *zero_pad_and_fft(LALStatus *status,
    REAL4TimeSeries *series,
    REAL8 deltaF,
    INT4 length,
    REAL4Window *window)
{
  /* variables */
  COMPLEX8FrequencySeries *zero_pad;
  RealFFTPlan *plan = NULL;
  SZeroPadAndFFTParameters zero_pad_params;

  /* create fft plan */
  plan = XLALCreateForwardREAL4FFTPlan(2 * series->data->length, 0);

  /* allocate memory */
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(status, &zero_pad, "zero_pad", \
        series->epoch, 0, deltaF, lalDimensionlessUnit, length), status);

  /* set zeropad parameters */
  zero_pad_params.fftPlan = plan;
  zero_pad_params.window = window->data;
  zero_pad_params.length = 2 * series->data->length;

  /* zero pad and fft */
  LAL_CALL(LALSZeroPadAndFFT(status, zero_pad, series, &zero_pad_params), \
      status);

  /* destroy fft plan */
  XLALDestroyREAL4FFTPlan(plan);

  return(zero_pad);
}

/* wrapper function for generating the cross correlation spectra */
static COMPLEX8FrequencySeries *cc_spectrum(LALStatus *status,
    COMPLEX8FrequencySeries *zero_pad_one,
    COMPLEX8FrequencySeries *zero_pad_two,
    COMPLEX8FrequencySeries *response_one,
    COMPLEX8FrequencySeries *response_two,
    REAL4FrequencySeries *opt_filter)
{
  /* variables */
  COMPLEX8FrequencySeries *series;
  StochasticCrossCorrelationCalInput cc_input;

  /* allocate memory */
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(status, &series, "cc_spectra", \
        opt_filter->epoch, opt_filter->f0, opt_filter->deltaF, \
        lalDimensionlessUnit, opt_filter->data->length), status);

  /* set inputs */
  cc_input.hBarTildeOne = zero_pad_one;
  cc_input.hBarTildeTwo = zero_pad_two;
  cc_input.responseFunctionOne = response_one;
  cc_input.responseFunctionTwo = response_two;
  cc_input.optimalFilter = opt_filter;
  
  /* calculate spectrum */
  LAL_CALL(LALStochasticCrossCorrelationSpectrumCal(status, series, \
        &cc_input, 1), status);

  return(series);
}

/* helper function to generate cross correlation statistic from cross
 * correlation spectra */
static REAL8 cc_statistic(COMPLEX8FrequencySeries *cc_spectra)
{
  /* variables */
  REAL8 cc_stat = 0;
  UINT4 i;

  /* sum up frequencies */
  for (i = 0; i < cc_spectra->data->length; i++)
  {
    cc_stat += cc_spectra->data->data[i].re;
  }

  /* normalise */
  cc_stat *= 2 * cc_spectra->deltaF;

  return(cc_stat);
}

/* helper function to save out xml tables */
static void save_xml_file(LALStatus *status,
    LALLeapSecAccuracy accuracy,
    CHAR *output_path,
    CHAR *base_name,
    StochasticTable *stochtable)
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
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--middle-segment");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (apply_mask_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--apply-mask");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (high_pass_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--high-pass-filter");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (overlap_hann_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--overlap-hann");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (recentre_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--recentre");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (cc_spectra_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--cc-spectra");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (vrbflg)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--verbose");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }

  /* add the xml comment, if specified */
  if (!*comment)
  {
    LALSnprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " ");
  }
  else
  {
    LALSnprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", \
        comment);
  }

  /* delete empty first entry in process params table */
  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free(this_proc_param);

  /* set xml output file */
  if (output_path)
  {
    LALSnprintf(xml_file_name, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml",
        output_path, base_name);
  }
  else
  {
    LALSnprintf(xml_file_name, FILENAME_MAX * sizeof(CHAR), "%s.xml",
        base_name);
  }

  /* write out xml */
  if (vrbflg)
    fprintf(stdout, "Writing output XML file...\n");

  /* opening xml file stream */
  memset(&xml_stream, 0, sizeof(LIGOLwXMLStream));
  LAL_CALL(LALOpenLIGOLwXMLFile(status, &xml_stream, xml_file_name), status);

  /* write out process and process params tables */
  LAL_CALL(LALGPSTimeNow(status, &(proctable.processTable->end_time), \
        &accuracy), status);
  LAL_CALL(LALBeginLIGOLwXMLTable(status, &xml_stream, process_table), \
      status);
  LAL_CALL(LALWriteLIGOLwXMLTable(status, &xml_stream, proctable, \
        process_table), status);
  LAL_CALL(LALEndLIGOLwXMLTable(status, &xml_stream), status);
  free(proctable.processTable);

  /* write the process params table */
  LAL_CALL(LALBeginLIGOLwXMLTable(status, &xml_stream, \
        process_params_table), status);
  LAL_CALL(LALWriteLIGOLwXMLTable(status, &xml_stream, procparams, \
        process_params_table), status);
  LAL_CALL(LALEndLIGOLwXMLTable(status, &xml_stream), status);

  /* write stochastic table */
  if (stochtable)
  {
    output_table.stochasticTable = stochtable;
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

/* function to perform the stochastic search */
static StochasticTable *stochastic_search(LALStatus *status,
    REAL4TimeSeries *series_one,
    REAL4TimeSeries *series_two,
    REAL4FrequencySeries *overlap,
    REAL4FrequencySeries *omega,
    REAL4Window *window,
    INT4 num_segments,
    INT4 filter_length,
    INT4 segs_in_interval,
    INT4 segment_length)
{
  /* counters */
  INT4 i, j, k;

  /* data structures */
  REAL4Vector *cal_psd_one;
  REAL4Vector *cal_psd_two;
  REAL4TimeSeries *segment_one = NULL;
  REAL4TimeSeries *segment_two = NULL;
  COMPLEX8FrequencySeries *response_one = NULL;
  COMPLEX8FrequencySeries *response_two = NULL;
  REAL4FrequencySeries *psd_one = NULL;
  REAL4FrequencySeries *psd_two = NULL;
  REAL4FrequencySeries *inv_psd_one = NULL;
  REAL4FrequencySeries *inv_psd_two = NULL;
  REAL4FrequencySeries *opt_filter = NULL;
  COMPLEX8FrequencySeries *zero_pad_one = NULL;
  COMPLEX8FrequencySeries *zero_pad_two = NULL;
  COMPLEX8FrequencySeries *cc_spectra = NULL;
  LALUnit countPerAttoStrain = {18,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* variables */
  INT4 num_intervals;
  INT4 segment_shift;
  INT4 middle_segment;
  LIGOTimeGPS seg_epoch;
  LIGOTimeGPS analysis_epoch;
  REAL8 delta_f;

  /* results */
  REAL8 sigma;
  REAL8 y;
  StochasticTable *stochHead = NULL;
  StochasticTable *thisStoch = NULL;

  /* debug file */
  CHAR debug_filename[FILENAME_MAX];

  /* initialise analysis_epoch */
  analysis_epoch.gpsSeconds = 0;
  analysis_epoch.gpsNanoSeconds = 0;

  /* calculate number of intervals, and required shift to get to next
   * interval */
  if (overlap_hann_flag)
  {
    num_intervals = (2 * (num_segments - 2)) - 1;
    segment_shift = segmentDuration / 2;
  }
  else
  {
    num_intervals = num_segments - 2;
    segment_shift = segmentDuration;
  }

  /* get middle segment number */
  middle_segment = (segs_in_interval - 1) / 2;

  /* get delta_f for optimal filter */
  delta_f = 1./(REAL8)PSD_WINDOW_DURATION;

  /* allocate memory for calibrated PSDs */
  cal_psd_one = XLALCreateREAL4Vector(filter_length);
  cal_psd_two = XLALCreateREAL4Vector(filter_length);

  if (vrbflg)
    fprintf(stdout, "Starting analysis loop...\n");

  /* loop over intervals */
  for (j = 0; j < num_intervals; j++)
  {	
    /* initialize average PSDs */
    for (i = 0; i < filter_length; i++)
    {
      cal_psd_one->data[i] = 0;
      cal_psd_two->data[i] = 0;
    }

    /* loop over segments in the interval */
    for (k = 0; k < segs_in_interval; k++)
    {
      /* set segment start time */
      LAL_CALL(LALAddFloatToGPS(status, &seg_epoch, &series_one->epoch, \
            (REAL8)((j * segment_shift) + (k * segmentDuration))), status);

      /* is this the analysis segment? */
      if (k == middle_segment)
        analysis_epoch = seg_epoch;

      if (vrbflg)
      {
        fprintf(stdout, "Request data at GPS time %d\n", \
            seg_epoch.gpsSeconds);
      }

      /* cut segments from series */
      segment_one = cut_time_series(status, series_one, seg_epoch, \
          segmentDuration);
      segment_two = cut_time_series(status, series_two, seg_epoch, \
          segmentDuration);

      /* save intermediate products */
      if (debug_flag)
      {
        LALSnprintf(debug_filename, FILENAME_MAX, "segment_1-%d.dat", \
            seg_epoch.gpsSeconds);
        LALSPrintTimeSeries(segment_one, debug_filename);
        LALSnprintf(debug_filename, FILENAME_MAX, "segment_2-%d.dat", \
            seg_epoch.gpsSeconds);
        LALSPrintTimeSeries(segment_two, debug_filename);
      }

      /* compute response */
      response_one = generate_response(status, ifoOne, calCacheOne, \
          seg_epoch, fMin, delta_f, countPerAttoStrain, filter_length, \
          calibOffset);
      response_two = generate_response(status, ifoTwo, calCacheTwo, \
          seg_epoch, fMin, delta_f, countPerAttoStrain, filter_length, \
          calibOffset);

      /* save intermediate products */
      if (debug_flag)
      {
        LALSnprintf(debug_filename, FILENAME_MAX, "response_1-%d.dat", \
            seg_epoch.gpsSeconds + calibOffset);
        LALCPrintFrequencySeries(response_one, debug_filename);
        LALSnprintf(debug_filename, FILENAME_MAX, "response_2-%d.dat", \
            seg_epoch.gpsSeconds + calibOffset);
        LALCPrintFrequencySeries(response_two, debug_filename);
      }

      /* check if on middle segment and if we want to include this in
       * the analysis */
      if ((k == middle_segment) && (!middle_segment_flag))
      {
        if (vrbflg)
          fprintf(stdout, "Ignoring middle segment..\n");
      }
      else
      {
        if (vrbflg)
          fprintf(stdout, "Estimating PSDs...\n");

        /* compute uncalibrated PSDs */
        psd_one = estimate_psd(status, segment_one, fMin, filter_length);
        psd_two = estimate_psd(status, segment_two, fMin, filter_length);

        if (vrbflg)
          fprintf(stdout, "Generating inverse noise...\n");

        /* compute inverse calibrate PSDs */
        inv_psd_one = inverse_noise(status, psd_one, response_one);
        inv_psd_two = inverse_noise(status, psd_two, response_two);

        /* sum over calibrated PSDs for average */
        for (i = 0; i < filter_length; i++)
        {
          cal_psd_one->data[i] += 1. / inv_psd_one->data->data[i];
          cal_psd_two->data[i] += 1. / inv_psd_two->data->data[i];
        }
      }
    }

    /* average calibrated PSDs and take inverse */
    for (i = 0; i < filter_length; i++)
    {
      /* average */
      if (!middle_segment_flag)
      {
        cal_psd_one->data[i] /= (REAL4)(segs_in_interval - 1);
        cal_psd_two->data[i] /= (REAL4)(segs_in_interval - 1);
      }
      else
      {
        cal_psd_one->data[i] /= (REAL4)segs_in_interval;
        cal_psd_two->data[i] /= (REAL4)segs_in_interval;
      }
      /* take inverse */
      inv_psd_one->data->data[i] = 1. / cal_psd_one->data[i];
      inv_psd_two->data->data[i] = 1. / cal_psd_two->data[i];
    }

    /* save intermediate products */
    if (debug_flag)
    {
      LALSnprintf(debug_filename, FILENAME_MAX, "inv_psd_1-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintFrequencySeries(inv_psd_one, debug_filename);
      LALSnprintf(debug_filename, FILENAME_MAX, "inv_psd_2-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintFrequencySeries(inv_psd_two, debug_filename);
    }

    if (vrbflg)
      fprintf(stdout, "Generating optimal filter...\n");

    /* build optimal filter */
    opt_filter = optimal_filter(status, overlap, omega, inv_psd_one, \
        inv_psd_two, window, &sigma);

    /* save intermediate products */
    if (debug_flag)
    {
      LALSnprintf(debug_filename, FILENAME_MAX, "opt_filter-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintFrequencySeries(opt_filter, debug_filename);
    }          

    if (vrbflg)
    {
      fprintf(stdout, "Analysing segment at GPS %d\n", \
          analysis_epoch.gpsSeconds);
    }

    /* cut analysis segment from full series */
    segment_one = cut_time_series(status, series_one, analysis_epoch, \
        segmentDuration);
    segment_two = cut_time_series(status, series_two, analysis_epoch, \
        segmentDuration);

    /* save intermediate products */
    if (debug_flag)
    {
      LALSnprintf(debug_filename, FILENAME_MAX, "analysis_segment_1-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintTimeSeries(segment_one, debug_filename);
      LALSnprintf(debug_filename, FILENAME_MAX, "analysis_segment_2-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALSPrintTimeSeries(segment_two, debug_filename);
    }

    if (vrbflg)
      fprintf(stdout, "Performing zero pad and FFT...\n");

    /* zero pad and fft */
    zero_pad_one = zero_pad_and_fft(status, segment_one, delta_f, \
        segment_length + 1, window);
    zero_pad_two = zero_pad_and_fft(status, segment_two, delta_f, \
        segment_length + 1, window);

    /* save intermediate products */
    if (debug_flag)
    {
      LALSnprintf(debug_filename, FILENAME_MAX, "zero_pad_1-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALCPrintFrequencySeries(zero_pad_one, debug_filename);
      LALSnprintf(debug_filename, FILENAME_MAX, "zero_pad_2-%d.dat", \
          analysis_epoch.gpsSeconds);
      LALCPrintFrequencySeries(zero_pad_two, debug_filename);
    }    

    if (vrbflg)
      fprintf(stdout, "Calculating cross correlation spectrum...\n");

    /* get response functions for analysis segments */
    response_one = generate_response(status, ifoOne, calCacheOne, \
        analysis_epoch, fMin, delta_f, countPerAttoStrain, filter_length, \
        calibOffset);
    response_two = generate_response(status, ifoTwo, calCacheTwo, \
        analysis_epoch, fMin, delta_f, countPerAttoStrain, filter_length, \
        calibOffset);

    /* calculate cc spectrum */
    cc_spectra = cc_spectrum(status, zero_pad_one, zero_pad_two, \
        response_one, response_two, opt_filter);

    if (cc_spectra_flag)
    {
      /* save out cc spectra as frame */
      if (vrbflg)
        fprintf(stdout, "Saving ccSpectra to frame...\n");
      write_ccspectra_frame(cc_spectra, ifoOne, ifoTwo, \
          analysis_epoch, segmentDuration);
    }
      
    /* cc statistic */
    y = cc_statistic(cc_spectra);

    /* display results */
    if (vrbflg)
    {
      fprintf(stdout, "Interval %d\n", j + 1);
      fprintf(stdout, "  GPS time  = %d\n", analysis_epoch.gpsSeconds);
      fprintf(stdout, "  y         = %e\n", y);
      fprintf(stdout, "  sigma     = %e\n", sigma);
    }

    /* allocate memory for table */
    if (!stochHead)
    {
      stochHead = thisStoch = (StochasticTable *) \
                  LALCalloc(1, sizeof(StochasticTable));
    }
    else
    {
      thisStoch = thisStoch->next = (StochasticTable *) \
                  LALCalloc(1, sizeof(StochasticTable));
    }

    /* populate columns */
    LALSnprintf(thisStoch->ifo_one, LIGOMETA_IFO_MAX, ifoOne);
    LALSnprintf(thisStoch->ifo_two, LIGOMETA_IFO_MAX, ifoTwo);
    LALSnprintf(thisStoch->channel_one, LIGOMETA_CHANNEL_MAX, channelOne);
    LALSnprintf(thisStoch->channel_two, LIGOMETA_CHANNEL_MAX, channelTwo);
    thisStoch->start_time.gpsSeconds = analysis_epoch.gpsSeconds;
    thisStoch->start_time.gpsNanoSeconds = analysis_epoch.gpsNanoSeconds;
    thisStoch->duration.gpsSeconds = segmentDuration;
    thisStoch->duration.gpsNanoSeconds = 0;
    thisStoch->f_min = fMin;
    thisStoch->f_max = fMax;
    thisStoch->cc_stat = y;
    thisStoch->cc_sigma = sigma;
  }

  /* clean up */
  XLALDestroyREAL4Vector(cal_psd_one);
  XLALDestroyREAL4Vector(cal_psd_two);
  XLALDestroyREAL4TimeSeries(segment_one);
  XLALDestroyREAL4TimeSeries(segment_two);
  XLALDestroyCOMPLEX8FrequencySeries(response_one);
  XLALDestroyCOMPLEX8FrequencySeries(response_two);
  XLALDestroyREAL4FrequencySeries(psd_one);
  XLALDestroyREAL4FrequencySeries(psd_two);
  XLALDestroyREAL4FrequencySeries(inv_psd_one);
  XLALDestroyREAL4FrequencySeries(inv_psd_two);
  XLALDestroyREAL4FrequencySeries(opt_filter);
  XLALDestroyCOMPLEX8FrequencySeries(zero_pad_one);
  XLALDestroyCOMPLEX8FrequencySeries(zero_pad_two);
  XLALDestroyCOMPLEX8FrequencySeries(cc_spectra);

  return(stochHead);
}

/* function for generating random noise */
static REAL4TimeSeries *generate_random_noise(LALStatus *status,
    INT4 duration,
    INT4 sample_rate)
{
  /* variables */
  REAL4TimeSeries *series;
  LIGOTimeGPS start;
  REAL8 deltaT;
  RandomParams *noise_params = NULL;
  struct timeval tv;

  /* initialise time */
  start.gpsSeconds = 0;
  start.gpsNanoSeconds = 0;

  /* get deltaT from sample_rate */
  deltaT = 1.0/(REAL8)(sample_rate);

  if (vrbflg)
    fprintf(stderr, "Allocating memory for random noise...\n");

  /* create and initialise time series */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, "noise", start, 0, \
        deltaT, lalSecondUnit, duration * sample_rate), status);

  /* get current time, for random seed */
  gettimeofday(&tv, NULL);

  if (vrbflg)
    fprintf(stderr, "Generating random noise...\n");

  /* generate noise_params */
  LAL_CALL(LALCreateRandomParams(status, &noise_params, tv.tv_usec), status);

  /* generate random noise */
  LAL_CALL(LALNormalDeviates(status, series->data, noise_params), status);

  /* destroy noise_params */
  LAL_CALL(LALDestroyRandomParams(status, &noise_params), status);

  return(series);
}

/* function for generating fake detector output */
static SSSimStochBGOutput *generate_fake_detector_output(LALStatus *status,
    REAL4TimeSeries *noise_one,
    REAL4TimeSeries *noise_two,
    REAL4FrequencySeries *omega,
    REAL8 deltaF)
{
  /* variables */
  REAL4TimeSeries *series_one;
  REAL4TimeSeries *series_two;
  SSSimStochBGInput input;
  SSSimStochBGParams params;
  SSSimStochBGOutput *output;
  COMPLEX8FrequencySeries *response_one = NULL;
  COMPLEX8FrequencySeries *response_two = NULL;
  LIGOTimeGPS start;
  INT4 freq_length;
  struct timeval tv;
  UINT4 i;

  /* initialise epoch */
  start.gpsSeconds = 0;
  start.gpsNanoSeconds = 0;

  /* get frequency length */
  freq_length = (noise_one->data->length / 2) + 1;

  if (vrbflg)
    fprintf(stderr, "Allocating memory for fake detector output...\n");

  /* create and initialise time series' */
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series_one, "one", \
        noise_one->epoch, noise_one->f0, noise_one->deltaT, \
        noise_one->sampleUnits, noise_one->data->length), status);
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series_two, "two", \
        noise_two->epoch, noise_two->f0, noise_two->deltaT, \
        noise_two->sampleUnits, noise_two->data->length), status);

  /* generate unity response functions */
  response_one = unity_response(status, start, 0, deltaF, \
      lalDimensionlessUnit, freq_length);
  response_two = unity_response(status, start, 0, deltaF, \
      lalDimensionlessUnit, freq_length);

  /* get current time, for random seed */
  gettimeofday(&tv, NULL);

  /* setup params for fake signal generation */
  params.length = noise_one->data->length;
  params.deltaT = noise_one->deltaT;
  params.seed = tv.tv_usec;
  params.detectorOne = lalCachedDetectors[LALDetectorIndexLHODIFF];
  params.detectorTwo = lalCachedDetectors[LALDetectorIndexLHODIFF];
  params.SSimStochBGTimeSeries1Unit = lalStrainUnit;
  params.SSimStochBGTimeSeries2Unit = lalStrainUnit;

  /* setup intputs for fake signal generation */
  input.omegaGW = omega;
  input.whiteningFilter1 = response_one;
  input.whiteningFilter2 = response_two;

  /* setup output for fake signal generation */
  output->SSimStochBG1 = series_one;
  output->SSimStochBG2 = series_two;

  /* generate fake signal */
  LAL_CALL(LALSSSimStochBGTimeSeries(status, output, &input, &params), \
      status);

  /* inject signal into noise */
  for (i = 0; i < series_one->data->length; i++)
  {
    series_one->data->data[i] += noise_one->data->data[i];
    series_two->data->data[i] += noise_two->data->data[i];
  }

  return(output);
}

/* display usage information */
static void display_usage()
{
  fprintf(stdout, "Usage: " PROGRAM_NAME " [options]\n");
  fprintf(stdout, " --help                        print this message\n");
  fprintf(stdout, " --version                     display version\n");
  fprintf(stdout, " --verbose                     verbose mode\n");
  fprintf(stdout, " --debug                       debug mode\n");
  fprintf(stdout, " --debug-level N               set lalDebugLevel\n");
  fprintf(stdout, " --user-tag STRING             set the user tag\n"); 
  fprintf(stdout, " --comment STRING              set the comment\n");
  fprintf(stdout, " --output-dir DIR              directory for output files\n");
  fprintf(stdout, " --cc-spectra                  save out cross correlation spectra\n");
  fprintf(stdout, " --gps-start-time N            GPS start time\n");
  fprintf(stdout, " --gps-end-time N              GPS end time\n");
  fprintf(stdout, " --interval-duration N         interval duration\n");
  fprintf(stdout, " --segment-duration N          segment duration\n");
  fprintf(stdout, " --resample-rate N             resample rate\n");
  fprintf(stdout, " --f-min N                     minimal frequency\n");
  fprintf(stdout, " --f-max N                     maximal frequency\n");
  fprintf(stdout, " --ifo-one IFO                 ifo for first stream\n");
  fprintf(stdout, " --ifo-two IFO                 ifo for second stream\n");
  fprintf(stdout, " --channel-one CHANNEL         channel for first stream\n");
  fprintf(stdout, " --channel-two CHANNEL         channel for second stream\n");
  fprintf(stdout, " --frame-cache-one FILE        cache file for first stream\n");
  fprintf(stdout, " --frame-cache-two FILE        cache file for second stream\n");
  fprintf(stdout, " --calibration-cache-one FILE  first stream calibration cache\n");
  fprintf(stdout, " --calibration-cache-two FILE  second stream calibration cache\n");
  fprintf(stdout, " --calibration-offset N        calibration offset\n");
  fprintf(stdout, " --apply-mask                  apply frequency masking\n");
  fprintf(stdout, " --mask-bin N                  number of bins to mask\n");
  fprintf(stdout, " --overlap-hann                overlaping hann windows\n");
  fprintf(stdout, " --hann-duration N             hann duration\n");
  fprintf(stdout, " --high-pass-filter            apply high pass filtering\n");
  fprintf(stdout, " --hpf-frequency N             high pass filter knee frequency\n");
  fprintf(stdout, " --hpf-attenuation N           high pass filter attenuation\n");
  fprintf(stdout, " --hpf-order N                 high pass filter order\n");
  fprintf(stdout, " --recentre                    recentre jobs\n");
  fprintf(stdout, " --middle-segment              use middle segment in PSD estimation\n");
  fprintf(stdout, " --geo-hpf-frequency N         GEO high pass filter knee frequency\n");
  fprintf(stdout, " --geo-hpf-attenuation N       GEO high pass filter attenuation\n");
  fprintf(stdout, " --geo-hpf-order N             GEO high pass filter order\n");
  fprintf(stdout, " --alpha N                     exponent on filter spectrum\n");
  fprintf(stdout, " --f-ref N                     reference frequency for filter spectrum\n");
  fprintf(stdout, " --omega0 N                    reference omega_0 for filter spectrum\n");
}

/* parse command line options */
static void parse_options(INT4 argc, CHAR *argv[])
{
  int c = -1;
  struct stat fileStatus;

  /* tempory variables */
  CHAR *channelOneTemp = NULL;
  CHAR *channelTwoTemp = NULL;

  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"middle-segment", no_argument, &middle_segment_flag, 1},
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"verbose", no_argument, &vrbflg, 1},
      {"recentre", no_argument, &recentre_flag, 1},
      {"cc-spectra", no_argument, &cc_spectra_flag, 1},
      {"debug", no_argument, &debug_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'a'},
      {"version", no_argument, 0, 'b'},
      {"debug-level", required_argument, 0, 'c'},
      {"user-tag", required_argument, 0, 'd'},
      {"comment", required_argument, 0, 'e'},
      {"output-dir", required_argument, 0, 'f'},
      {"gps-start-time", required_argument, 0, 'g'},
      {"gps-end-time", required_argument, 0, 'h'},
      {"interval-duration", required_argument, 0, 'i'},
      {"segment-duration", required_argument, 0, 'j'},
      {"resample-rate", required_argument, 0, 'k'},
      {"f-min", required_argument, 0, 'l'},
      {"f-max", required_argument, 0, 'm'},
      {"ifo-one", required_argument, 0, 'n'},
      {"ifo-two", required_argument, 0, 'o'},
      {"channel-one", required_argument, 0, 'p'},
      {"channel-two", required_argument, 0, 'q'},
      {"frame-cache-one", required_argument, 0, 'r'},
      {"frame-cache-two", required_argument, 0, 's'},
      {"calibration-cache-one", required_argument, 0, 't'},
      {"calibration-cache-two", required_argument, 0, 'u'},
      {"calibration-offset", required_argument, 0, 'v'},
      {"mask-bin", required_argument, 0, 'w'},
      {"hann-duration", required_argument, 0, 'x'},
      {"hpf-frequency", required_argument, 0, 'y'},
      {"hpf-attenuation", required_argument, 0, 'z'},
      {"hpf-order", required_argument, 0, 'A'},
      {"geo-hpf-frequency", required_argument, 0, 'E'},
      {"geo-hpf-attenuation", required_argument, 0, 'F'},
      {"geo-hpf-order", required_argument, 0, 'G'},
      {"alpha", required_argument, 0, 'H'},
      {"f-ref", required_argument, 0, 'I'},
      {"omega0", required_argument, 0, 'J'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, \
        "abc:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:" \
        "A:E:F:G:H:I:J:", long_options, &option_index);

    if (c == -1)
    {
      /* end of options, break loop */
      break;
    }

    switch(c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
        {
          break;
        }
        else
        {
          fprintf(stderr, "error parseing option %s with argument %s\n", \
              long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'a':
        /* help */
        display_usage();
        exit(0);
        break;

      case 'b':
        /* display version info and exit */
        fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
        exit(0);
        break;

      case 'c':
        /* debug level */
        set_debug_level( optarg );
        ADD_PROCESS_PARAM("string", "%s", optarg);
        break;

      case 'd':
        /* user tag */
        optarg_len = strlen(optarg) + 1;
        userTag = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(userTag, optarg, optarg_len);

        /* add to process_params table */
        this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                          calloc(1, sizeof(ProcessParamsTable));
        LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
            PROGRAM_NAME);
        LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--user-tag");
        LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
        LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, "%s", optarg);
        break;

      case 'e':
        /* xml comment */
        if (strlen(optarg) > LIGOMETA_COMMENT_MAX - 1)
        {
          fprintf(stderr, "invalid argument to --%s:\n" \
              "comment must be less than %d characters\n", \
              long_options[option_index].name, LIGOMETA_COMMENT_MAX);
          exit(1);
        }
        else
        {
          LALSnprintf(comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'f':
        /* directory for output files */
        optarg_len = strlen(optarg) + 1;
        outputPath = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(outputPath, optarg, optarg_len);
        if (stat(outputPath, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Directory does not exist: (%s specified)\n", \
              long_options[option_index].name, outputPath);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", outputPath);
        break;

      case 'g':
        /* start time */
        startTime = atoi(optarg);
        if (startTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is prior to 1 January 1994 00:00:00 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }
        if (startTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is after 14 September 2011 01:46:26 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%ld", startTime);
        break;

      case 'h':
        /* end time */
        endTime = atoi(optarg);
        if (endTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is prior to 1 January 1994 00:00:00 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }
        if (endTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is after 14 September 2011 01:46:26 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%ld", endTime);
        break;

      case 'i':
        /* interval duration */
        intervalDuration = atoi(optarg);
        if (intervalDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Interval duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, intervalDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", intervalDuration);
        break;

      case 'j':
        /* segment duration */
        segmentDuration = atoi(optarg);
        if (segmentDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Segment duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, segmentDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", segmentDuration);
        break;

      case 'k':
        /* resample rate */
        resampleRate = atoi(optarg);
        if (resampleRate < 2 || resampleRate > 16384 || resampleRate % 2)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Resample rate must be a power of 2 between 2 and 16384: " \
              "inclusive: (%d specified)\n", long_options[option_index].name, \
              resampleRate);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", resampleRate);

        break;

      case 'l':
        /* minimal frequency */
        fMin = atof(optarg);
        if (fMin < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Minimum frequency is less than 0 Hz (%f specified)\n", \
              long_options[option_index].name, fMin);
          exit(1);
        }
        /* check that min frequency can be represented by the
         * sampling rate of the data and round accordingly */
        if (fMin != myround(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMin = myround(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
          fprintf(stderr, "warning: fMin has been rounded to %f\n", fMin);
        }
        ADD_PROCESS_PARAM("float", "%e", fMin);
        break;

      case 'm':
        /* maximal frequency */
        fMax = atof(optarg);
        if (fMax < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Maximum frequency is less than 0 Hz (%f specified)\n", \
              long_options[option_index].name, fMax);
          exit(1);
        }
        /* check that the max frequency can be represented by the
         * sampling rate of the data and round accordingly */
        if (fMax != myround(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMax = myround(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
          fprintf(stderr, "warning: fMax has been rounded to %f\n", fMax);
        }
        ADD_PROCESS_PARAM("float", "%e", fMax);
        break;

      case 'n':
        /* ifo for first stream */
        optarg_len = strlen(optarg) + 1;
        ifoOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifoOne, optarg, optarg_len);

        /* set site id for ifo one */
        if (strncmp(ifoOne, "H1", 2) == 0)
          siteOne = 0;
        else if (strncmp(ifoOne, "H2", 2) == 0)
          siteOne = 0;
        else if (strncmp(ifoOne, "L1", 2) == 0)
          siteOne = 1;
        else if (strncmp(ifoOne, "G1", 2) == 0)
          siteOne = 3;
        else
        {
          fprintf(stderr, "First IFO not recognised...\n");
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", ifoOne);
        break;

      case 'o':
        /* ifo for second stream */
        optarg_len = strlen(optarg) + 1;
        ifoTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifoTwo, optarg, optarg_len);

        /* set site id for ifo two */
        if (strncmp(ifoTwo, "H1", 2) == 0)
          siteTwo = 0;
        else if (strncmp(ifoTwo, "H2", 2) == 0)
          siteTwo = 0;
        else if (strncmp(ifoTwo, "L1", 2) == 0)
          siteTwo = 1;
        else if (strncmp(ifoTwo, "G1", 2) == 0)
          siteOne = 3;
        else
        {
          fprintf(stderr, "Second IFO not recognised...\n");
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", ifoTwo);
        break;

      case 'p':
        /* channel one */
        optarg_len = strlen(optarg) + 4;
        channelOneTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(channelOneTemp, optarg, optarg_len);
        ADD_PROCESS_PARAM("string", "%s", channelOneTemp);
        break;

      case 'q':
        /* channel two */
        optarg_len = strlen(optarg) + 4;
        channelTwoTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(channelTwoTemp, optarg, optarg_len);
        ADD_PROCESS_PARAM("string", "%s", channelTwoTemp);
        break;

      case 'r':
        /* frame cache one */
        optarg_len = strlen(optarg) + 1;
        frameCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(frameCacheOne, optarg, optarg_len);
        if (stat(frameCacheOne, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheOne);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", frameCacheOne);
        break;

      case 's':
        /* frame cache two */
        optarg_len = strlen(optarg) + 1;
        frameCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(frameCacheTwo, optarg, optarg_len);
        if (stat(frameCacheTwo, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheTwo);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", frameCacheTwo);
        break;

      case 't':
        /* calibration cache one */
        optarg_len = strlen(optarg) + 1;
        calCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(calCacheOne, optarg, optarg_len);
        if (stat(calCacheOne, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheOne);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", calCacheOne);
        break;

      case 'u':
        /* calibration cache two */
        optarg_len = strlen(optarg) + 1;
        calCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(calCacheTwo, optarg, optarg_len);
        if (stat(calCacheTwo, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheTwo);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", calCacheTwo);
        break;

      case 'v':
        /* calibration time offset */
        calibOffset = atoi(optarg);
        ADD_PROCESS_PARAM("int", "%d", calibOffset);
        break;

      case 'w':
        /* number of bins to mask for frequency mask */
        maskBin = atoi(optarg);
        if (maskBin <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Number of bins to mask must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name, maskBin);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", maskBin);
        break;

      case 'x':
        /* hann window duration */
        hannDuration = atoi(optarg);
        if (hannDuration < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Hann duartion is less than 0: (%d specified)\n", \
              long_options[option_index].name, hannDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", hannDuration);
        break;

      case 'y':
        /* high pass knee filter frequency  */
        highPassFreq = atof(optarg);
        if (highPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "High pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              highPassFreq);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", highPassFreq);
        break;

      case 'z':
        /* high pass filter attenuation  */
        highPassAtten = atof(optarg);
        if ((highPassAtten < 0.0) || (highPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              highPassAtten);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", highPassAtten);
        break;

      case 'A':
        /* high pass filter order  */
        highPassOrder = atoi(optarg);
        if (highPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              highPassOrder);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", highPassOrder);
        break;

      case 'E':
        /* GEO high pass knee filter frequency */
        geoHighPassFreq = atof(optarg);
        if (geoHighPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "GEO high pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              geoHighPassFreq);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", geoHighPassFreq);
        break;

      case 'F':
        /* GEO high pass filter attenuation */
        geoHighPassAtten = atof(optarg);
        if ((geoHighPassAtten < 0.0) || (geoHighPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GEO high pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              geoHighPassAtten);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", geoHighPassAtten);
        break;

      case 'G':
        /* GEO high pass filter order */
        geoHighPassOrder = atoi(optarg);
        if (geoHighPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GEO high pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              geoHighPassOrder);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", geoHighPassOrder);
        break;

      case 'H':
        /* filter spectrum exponent */
        alpha = atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", alpha);
        break;

      case 'I':
        /* filter reference frequency */
        fRef = atof(optarg);
        if (fRef < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference frequency must be greater than 0: " \
              "(%f specified)\n", long_options[option_index].name, fRef);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", fRef);
        break;

      case 'J':
        /* filter reference omega */
        omegaRef = atof(optarg);
        if (omegaRef <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference omega_0 must be positive: (%f specified)\n", \
              long_options[option_index].name, omegaRef);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", omegaRef);
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing options\n");
        exit(1);
    }
  }

  if (optind < argc)
  {
    fprintf(stderr, "Extraneous command line arguments:\n");
    while(optind < argc)
    {
      fprintf(stderr, "%s\n", argv[optind++]);
    }
    exit(1);
  }

  /* check for required arguments */

  /* start/end time */
  if (startTime == 0)
  {
    fprintf(stderr, "--gps-start-time must be specified\n");
    exit(1);
  }
  if (endTime == 0)
  {
    fprintf(stderr, "--gps-end-time must be specified\n");
    exit(1);
  }

  /* interval duration */
  if (intervalDuration == -1)
  {
    fprintf(stderr, "--interval-duration must be specified\n");
    exit(1);
  }

  /* segment duration */
  if (segmentDuration == -1)
  {
    fprintf(stderr, "--segment-duration must be specified\n");
    exit(1);
  }

  /* min/max frequency */
  if (fMin == -1)
  {
    fprintf(stderr, "--f-min must be specified\n");
    exit(1);
  }
  if (fMax == -1)
  {
    fprintf(stderr, "--f-max must be specified\n");
    exit(1);
  }

  /* ifos */
  if (ifoOne == NULL)
  {
    fprintf(stderr, "--ifo-one must be specified\n");
    exit(1);
  }
  if (ifoTwo == NULL)
  {
    fprintf(stderr, "--ifo-two must be specified\n");
    exit(1);
  }

  /* channels */
  if (channelOne == NULL)
  {
    fprintf(stderr, "--channel-one must be specified\n");
    exit(1);
  }
  if (channelTwo == NULL)
  {
    fprintf(stderr, "--channel-two must be specified\n");
    exit(1);
  }

  /* frame cache */
  if (frameCacheOne == NULL)
  {
    fprintf(stderr, "--frame-cache-one must be specified\n");
    exit(1);
  }
  if (siteOne != siteTwo)
  {
    /* only need second frame cache if ifos differ */
    if (frameCacheTwo == NULL)
    {
      fprintf(stderr, "--frame-cache-two must be specified\n");
      exit(1);
    }
  }
  else
  {
    /* if site ids are the same then the frames for the different
     * detectors are in the same frame cache */
    frameCacheTwo = frameCacheOne;
  }

  /* calibration cache */
  if (strncmp(ifoOne, "G1", 2) != 0)
  {
    if (calCacheOne == NULL)
    {
      fprintf(stderr, "--calibration-cache-one must be specified\n");
      exit(1);
    }
  }
  if (strncmp(ifoTwo, "G1", 2) != 0)
  {
    if (calCacheTwo == NULL)
    {
      fprintf(stderr, "--calibration-cache-two must be specified\n");
      exit(1);
    }
  }

  /* calibration offset */
  if (calibOffset == -1)
  {
    fprintf(stderr, "--calibration-offset must be specified\n");
    exit(1);
  }

  /* mask */
  if (apply_mask_flag)
  {
    if (maskBin == -1)
    {
      fprintf(stderr, "--mask-bin must be specified\n");
      exit(1);
    }
  }

  /* hann duration */
  if (overlap_hann_flag)
  {
    if (hannDuration != -1)
    {
      fprintf(stderr, "Overlapping Hann windows specified, --hann-duration " \
          "will be ignored\n");
    }
  }
  else
  {
    if (hannDuration == -1)
    {
      fprintf(stderr, "--hann-duration must be specified\n");
      exit(1);
    }
  }

  /* high pass filter */
  if (high_pass_flag)
  {
    if (highPassFreq == -1)
    {
      fprintf(stderr, "--hpf-frequency must be specified\n");
      exit(1);
    }
    if (highPassAtten == -1)
    {
      fprintf(stderr, "--hpf-attenuation must be specified\n");
      exit(1);
    }
    if (highPassOrder == -1)
    {
      fprintf(stderr, "--hpf-order must be specified\n");
      exit(1);
    }
  }

  /* GEO high pass filter */
  if ((strncmp(ifoOne, "G1", 2) == 0) || (strncmp(ifoTwo, "G1", 2) == 0))
  {  
    if (geoHighPassFreq == -1)
    {
      fprintf(stderr, "--geo-hpf-frequency must be specified\n");
      exit(1);
    }
    if (geoHighPassAtten == -1)
    {
      fprintf(stderr, "--geo-hpf-attenuation must be specified\n");
      exit(1);
    }
    if (geoHighPassOrder == -1)
    {
      fprintf(stderr, "--geo-hpf-order must be specified\n");
      exit(1);
    }
  }

  /* check for sensible arguments */

  /* start time same as stop time */
  if (startTime == endTime)
  {
    fprintf(stderr, "Start time same as end time; no analysis to perform\n");
    exit(1);
  }

  /* stop time before start time */
  if (startTime > endTime)
  {
    fprintf(stderr, "Invalid start/end time; end time (%d) is before " \
        "start time (%d)\n", endTime, startTime);
    exit(1);
  }

  /* interval duration must be a least 3 times the segment duration */
  if ((intervalDuration / segmentDuration) < 3)
  {
    fprintf(stderr, "Invalid interval duration (%d): must be a least 3 times " \
        "the segment\nduration (%d)\n", intervalDuration, segmentDuration);
    exit(1);
  }

  /* interval duration must be an odd mutliple of segment duration */
  if (((intervalDuration / segmentDuration) % 2) != 1)
  {
    fprintf(stderr, "Invalid interval duration (%d): must be an odd " \
        "multiple of the segment\nduration (%d)\n", intervalDuration, \
        segmentDuration);
    exit(1);
  }

  /* min frequency same as max */
  if (fMin == fMax)
  {
    fprintf(stderr, "Minimum frequency same as maximum; no analysis to " \
        "perform\n");
    exit(1);
  }

  /* max frequency less than min */
  if (fMin > fMax)
  {
    fprintf(stderr, "Invalid frequency band; maximum frequency (%f Hz) is " \
        "before minimum\nfrequency (%f Hz)\n", fMax, fMin);
    exit(1);
  }

  /* filter reference frequency less than min */
  if (fRef < fMin)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is less than minimum " \
        "frequency (%f Hz)\n", fRef, fMin);
    exit(1);
  }

  /* filter reference frequency greater than max */
  if (fRef > fMax)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is greater than maximum " \
        "frequency (%f Hz)\n", fRef, fMax);
    exit(1);
  }

  /* set channels */
  strcpy(channelOne, ifoOne);
  strcpy(channelTwo, ifoTwo);
  strcat(channelOne, ":");
  strcat(channelTwo, ":");
  strcat(channelOne, channelOneTemp);
  strcat(channelTwo, channelTwoTemp);
  free(channelOneTemp);
  free(channelTwoTemp);

  return;
}

/* display usage information for fake data*/
static void display_usage_fake()
{
  fprintf(stdout, "Usage: " PROGRAM_NAME " [options]\n");
  fprintf(stdout, " --help                        print this message\n");
  fprintf(stdout, " --version                     display version\n");
  fprintf(stdout, " --verbose                     verbose mode\n");
  fprintf(stdout, " --debug                       debug mode\n");
  fprintf(stdout, " --debug-level N               set lalDebugLevel\n");
  fprintf(stdout, " --user-tag STRING             set the user tag\n"); 
  fprintf(stdout, " --comment STRING              set the comment\n");
  fprintf(stdout, " --output-dir DIR              directory for output files\n");
  fprintf(stdout, " --cc-spectra                  save out cross correlation spectra\n");
  fprintf(stdout, " --duration N                  duration of fake data to generate\n");
  fprintf(stdout, " --interval-duration N         interval duration\n");
  fprintf(stdout, " --segment-duration N          segment duration\n");
  fprintf(stdout, " --resample-rate N             resample rate\n");
  fprintf(stdout, " --f-min N                     minimal frequency\n");
  fprintf(stdout, " --f-max N                     maximal frequency\n");
  fprintf(stdout, " --apply-mask                  apply frequency masking\n");
  fprintf(stdout, " --mask-bin N                  number of bins to mask\n");
  fprintf(stdout, " --overlap-hann                overlaping hann windows\n");
  fprintf(stdout, " --hann-duration N             hann duration\n");
  fprintf(stdout, " --high-pass-filter            apply high pass filtering\n");
  fprintf(stdout, " --hpf-frequency N             high pass filter knee frequency\n");
  fprintf(stdout, " --hpf-attenuation N           high pass filter attenuation\n");
  fprintf(stdout, " --hpf-order N                 high pass filter order\n");
  fprintf(stdout, " --recentre                    recentre jobs\n");
  fprintf(stdout, " --middle-segment              use middle segment in PSD estimation\n");
  fprintf(stdout, " --alpha N                     exponent on filter spectrum\n");
  fprintf(stdout, " --f-ref N                     reference frequency for filter spectrum\n");
  fprintf(stdout, " --omega0 N                    reference omega_0 for filter spectrum\n");
}

/* parse command line options for fake data */
static void parse_options_fake(INT4 argc, CHAR *argv[])
{
  int c = -1;
  struct stat fileStatus;

  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"middle-segment", no_argument, &middle_segment_flag, 1},
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"verbose", no_argument, &vrbflg, 1},
      {"cc-spectra", no_argument, &cc_spectra_flag, 1},
      {"debug", no_argument, &debug_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'a'},
      {"version", no_argument, 0, 'b'},
      {"debug-level", required_argument, 0, 'c'},
      {"user-tag", required_argument, 0, 'd'},
      {"comment", required_argument, 0, 'e'},
      {"output-dir", required_argument, 0, 'f'},
      {"duration", required_argument, 0, 'g'},
      {"interval-duration", required_argument, 0, 'i'},
      {"segment-duration", required_argument, 0, 'j'},
      {"resample-rate", required_argument, 0, 'k'},
      {"f-min", required_argument, 0, 'l'},
      {"f-max", required_argument, 0, 'm'},
      {"mask-bin", required_argument, 0, 'w'},
      {"hann-duration", required_argument, 0, 'x'},
      {"hpf-frequency", required_argument, 0, 'y'},
      {"hpf-attenuation", required_argument, 0, 'z'},
      {"hpf-order", required_argument, 0, 'A'},
      {"alpha", required_argument, 0, 'H'},
      {"f-ref", required_argument, 0, 'I'},
      {"omega0", required_argument, 0, 'J'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, \
        "abc:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:" \
        "A:E:F:G:H:I:J:", long_options, &option_index);

    if (c == -1)
    {
      /* end of options, break loop */
      break;
    }

    switch(c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
        {
          break;
        }
        else
        {
          fprintf(stderr, "error parseing option %s with argument %s\n", \
              long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'a':
        /* help */
        display_usage();
        exit(0);
        break;

      case 'b':
        /* display version info and exit */
        fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
        exit(0);
        break;

      case 'c':
        /* debug level */
        set_debug_level( optarg );
        ADD_PROCESS_PARAM("string", "%s", optarg);
        break;

      case 'd':
        /* user tag */
        optarg_len = strlen(optarg) + 1;
        userTag = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(userTag, optarg, optarg_len);

        /* add to process_params table */
        this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                          calloc(1, sizeof(ProcessParamsTable));
        LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
            PROGRAM_NAME);
        LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--user-tag");
        LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
        LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, "%s", optarg);
        break;

      case 'e':
        /* xml comment */
        if (strlen(optarg) > LIGOMETA_COMMENT_MAX - 1)
        {
          fprintf(stderr, "invalid argument to --%s:\n" \
              "comment must be less than %d characters\n", \
              long_options[option_index].name, LIGOMETA_COMMENT_MAX);
          exit(1);
        }
        else
        {
          LALSnprintf(comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'f':
        /* directory for output files */
        optarg_len = strlen(optarg) + 1;
        outputPath = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(outputPath, optarg, optarg_len);
        if (stat(outputPath, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Directory does not exist: (%s specified)\n", \
              long_options[option_index].name, outputPath);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", outputPath);
        break;

      case 'g':
        /* duration */
        totalDuration = atoi(optarg);
        if (totalDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Total duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, totalDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", totalDuration);
        break;

      case 'i':
        /* interval duration */
        intervalDuration = atoi(optarg);
        if (intervalDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Interval duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, intervalDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", intervalDuration);
        break;

      case 'j':
        /* segment duration */
        segmentDuration = atoi(optarg);
        if (segmentDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Segment duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, segmentDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", segmentDuration);
        break;

      case 'k':
        /* resample rate */
        resampleRate = atoi(optarg);
        if (resampleRate < 2 || resampleRate > 16384 || resampleRate % 2)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Resample rate must be a power of 2 between 2 and 16384: " \
              "inclusive: (%d specified)\n", long_options[option_index].name, \
              resampleRate);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", resampleRate);

        break;

      case 'l':
        /* minimal frequency */
        fMin = atof(optarg);
        if (fMin < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Minimum frequency is less than 0 Hz (%f specified)\n", \
              long_options[option_index].name, fMin);
          exit(1);
        }
        /* check that min frequency can be represented by the
         * sampling rate of the data and round accordingly */
        if (fMin != myround(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMin = myround(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
          fprintf(stderr, "warning: fMin has been rounded to %f\n", fMin);
        }
        ADD_PROCESS_PARAM("float", "%e", fMin);
        break;

      case 'm':
        /* maximal frequency */
        fMax = atof(optarg);
        if (fMax < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Maximum frequency is less than 0 Hz (%f specified)\n", \
              long_options[option_index].name, fMax);
          exit(1);
        }
        /* check that the max frequency can be represented by the
         * sampling rate of the data and round accordingly */
        if (fMax != myround(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMax = myround(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
          fprintf(stderr, "warning: fMax has been rounded to %f\n", fMax);
        }
        ADD_PROCESS_PARAM("float", "%e", fMax);
        break;

      case 'x':
        /* hann window duration */
        hannDuration = atoi(optarg);
        if (hannDuration < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Hann duartion is less than 0: (%d specified)\n", \
              long_options[option_index].name, hannDuration);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", hannDuration);
        break;

      case 'y':
        /* high pass knee filter frequency  */
        highPassFreq = atof(optarg);
        if (highPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "High pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              highPassFreq);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", highPassFreq);
        break;

      case 'z':
        /* high pass filter attenuation  */
        highPassAtten = atof(optarg);
        if ((highPassAtten < 0.0) || (highPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              highPassAtten);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", highPassAtten);
        break;

      case 'A':
        /* high pass filter order  */
        highPassOrder = atoi(optarg);
        if (highPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              highPassOrder);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", highPassOrder);
        break;

      case 'H':
        /* filter spectrum exponent */
        alpha = atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", alpha);
        break;

      case 'I':
        /* filter reference frequency */
        fRef = atof(optarg);
        if (fRef < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference frequency must be greater than 0: " \
              "(%f specified)\n", long_options[option_index].name, fRef);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", fRef);
        break;

      case 'J':
        /* filter reference omega */
        omegaRef = atof(optarg);
        if (omegaRef <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference omega_0 must be positive: (%f specified)\n", \
              long_options[option_index].name, omegaRef);
          exit(1);
        }
        ADD_PROCESS_PARAM("float", "%e", omegaRef);
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing options\n");
        exit(1);
    }
  }

  if (optind < argc)
  {
    fprintf(stderr, "Extraneous command line arguments:\n");
    while(optind < argc)
    {
      fprintf(stderr, "%s\n", argv[optind++]);
    }
    exit(1);
  }

  /* check for required arguments */

  /* duration */
  if (totalDuration == -1)
  {
    fprintf(stderr, "--duration must be specified\n");
    exit(1);
  }
  
  /* interval duration */
  if (intervalDuration == -1)
  {
    fprintf(stderr, "--interval-duration must be specified\n");
    exit(1);
  }

  /* segment duration */
  if (segmentDuration == -1)
  {
    fprintf(stderr, "--segment-duration must be specified\n");
    exit(1);
  }

  /* min/max frequency */
  if (fMin == -1)
  {
    fprintf(stderr, "--f-min must be specified\n");
    exit(1);
  }
  if (fMax == -1)
  {
    fprintf(stderr, "--f-max must be specified\n");
    exit(1);
  }

  /* mask */
  if (apply_mask_flag)
  {
    if (maskBin == -1)
    {
      fprintf(stderr, "--mask-bin must be specified\n");
      exit(1);
    }
  }

  /* hann duration */
  if (overlap_hann_flag)
  {
    if (hannDuration != -1)
    {
      fprintf(stderr, "Overlapping Hann windows specified, --hann-duration " \
          "will be ignored\n");
    }
  }
  else
  {
    if (hannDuration == -1)
    {
      fprintf(stderr, "--hann-duration must be specified\n");
      exit(1);
    }
  }

  /* high pass filter */
  if (high_pass_flag)
  {
    if (highPassFreq == -1)
    {
      fprintf(stderr, "--hpf-frequency must be specified\n");
      exit(1);
    }
    if (highPassAtten == -1)
    {
      fprintf(stderr, "--hpf-attenuation must be specified\n");
      exit(1);
    }
    if (highPassOrder == -1)
    {
      fprintf(stderr, "--hpf-order must be specified\n");
      exit(1);
    }
  }

  /* check for sensible arguments */

  /* interval duration must be a least 3 times the segment duration */
  if ((intervalDuration / segmentDuration) < 3)
  {
    fprintf(stderr, "Invalid interval duration (%d): must be a least 3 times " \
        "the segment\nduration (%d)\n", intervalDuration, segmentDuration);
    exit(1);
  }

  /* interval duration must be an odd mutliple of segment duration */
  if (((intervalDuration / segmentDuration) % 2) != 1)
  {
    fprintf(stderr, "Invalid interval duration (%d): must be an odd " \
        "multiple of the segment\nduration (%d)\n", intervalDuration, \
        segmentDuration);
    exit(1);
  }

  /* min frequency same as max */
  if (fMin == fMax)
  {
    fprintf(stderr, "Minimum frequency same as maximum; no analysis to " \
        "perform\n");
    exit(1);
  }

  /* max frequency less than min */
  if (fMin > fMax)
  {
    fprintf(stderr, "Invalid frequency band; maximum frequency (%f Hz) is " \
        "before minimum\nfrequency (%f Hz)\n", fMax, fMin);
    exit(1);
  }

  /* filter reference frequency less than min */
  if (fRef < fMin)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is less than minimum " \
        "frequency (%f Hz)\n", fRef, fMin);
    exit(1);
  }

  /* filter reference frequency greater than max */
  if (fRef > fMax)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is greater than maximum " \
        "frequency (%f Hz)\n", fRef, fMax);
    exit(1);
  }

  return;
}

/* program entry point */
INT4 main(INT4 argc, CHAR *argv[])
{
  /* lal initialisation variables */
  LALStatus status = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  
  /* xml */
  CHAR baseName[FILENAME_MAX];
  StochasticTable *stochHead = NULL;
  StochasticTable *thisStoch = NULL;

  /* counter */
  INT4 i;

  /* data structures */
  REAL4TimeSeries *seriesOne;
  REAL4TimeSeries *seriesTwo;
  REAL4FrequencySeries *overlap;
  REAL4FrequencySeries *mask;
  REAL4FrequencySeries *omegaGW;
  REAL4Window *dataWindow;
  LIGOTimeGPS gpsStartTime;
  LIGOTimeGPS gpsEndTime;

  /* variables */
  INT4 numSegments;
  INT4 duration, durationEff, extrasec;
  INT4 segsInInt;
  INT4 segmentLength;
  INT4 padData;
  REAL8 deltaF;
  INT4 filterLength;
  INT4 numFMin;
  INT4 numFMax;

  /* error handler */
  status.statusPtr = NULL;
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level("3");

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc(1, sizeof(ProcessTable));
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->start_time), \
        &accuracy), &status);
  LAL_CALL(populate_process_table(&status, proctable.processTable, \
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &status);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) \
                    calloc(1, sizeof(ProcessParamsTable));
  memset(comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR));

  /* parse command line options */
  parse_options(argc, argv);

  /* get xml file basename */
  if (userTag)
  {
    LALSnprintf(baseName, FILENAME_MAX, "%s%s-STOCHASTIC_%s_%d-%d", \
        ifoOne, ifoTwo, userTag, startTime, (endTime - startTime));
  }
  else
  {
    LALSnprintf(baseName, FILENAME_MAX, "%s%s-STOCHASTIC-%d-%d", \
        ifoOne, ifoTwo, startTime, (endTime - startTime));
  }

  /* get number of segments */
  duration = endTime - startTime;
  numSegments = duration / segmentDuration;
  segsInInt = intervalDuration / segmentDuration;
  
  /* recentre */
  if (recentre_flag)
  {
    if (vrbflg)
      fprintf(stdout, "Recentring within data stream...\n");

    durationEff = numSegments * segmentDuration;
    extrasec = duration - durationEff;
    startTime += extrasec / 2;
    endTime = startTime + durationEff;
  }

  /* add a resample buffer, if required */
  if ((resampleRate) || (high_pass_flag))
    padData = 1;
  else
    padData = 0;

  /* initialise gps time structures */
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0;
  gpsEndTime.gpsSeconds = endTime;
  gpsEndTime.gpsNanoSeconds = 0;

  /* read data */
  seriesOne = get_time_series(&status, ifoOne, frameCacheOne, channelOne, \
      gpsStartTime, gpsEndTime, padData);
  seriesTwo = get_time_series(&status, ifoTwo, frameCacheTwo, channelTwo, \
      gpsStartTime, gpsEndTime, padData);

  /* check that the two series have the same sample rate */
  if (seriesOne->deltaT != seriesTwo->deltaT)
  {
    fprintf(stderr, "Series have different sample rates...\n");
    exit(1);
  }
  else
  {
    /* get resample rate, if required */
    if (!resampleRate)
      resampleRate = (INT4)(1./seriesOne->deltaT);
  }

  /* get deltaF for optimal filter */
  deltaF = 1./(REAL8)PSD_WINDOW_DURATION;

  /* get bins for min and max frequencies */
  numFMin = (INT4)(fMin / deltaF);
  numFMax = (INT4)(fMax / deltaF);

  /* get lengths */
  filterLength = numFMax - numFMin + 1;
  segmentLength = segmentDuration * resampleRate;

  if (vrbflg)
    fprintf(stdout, "Generating data segment window...\n");

  /* for overlapping hann windows, the hann window length is the segment
   * length */
  if (overlap_hann_flag)
    hannDuration = segmentDuration;

  /* create window for data */
  dataWindow = data_window(seriesOne->deltaT, segmentLength, hannDuration);

  /* generate overlap reduction function */
  if (vrbflg)
    fprintf(stdout, "Generating the overlap reduction function...\n");
  overlap = overlap_reduction_function(&status, filterLength, fMin, deltaF, \
      siteOne, siteTwo);

  /* generage omegaGW */
  if (vrbflg)
    fprintf(stdout, "Generating spectrum for optimal filter...\n");
  omegaGW = omega_gw(&status, alpha, fRef, omegaRef, filterLength, \
      fMin, deltaF);

  /* frequency mask */
  if (apply_mask_flag)
  {
    if (vrbflg)
      fprintf(stdout, "Applying frequency mask to spectrum..\n");

    /* generate frequency mask */
    mask = frequency_mask(&status, fMin, deltaF, filterLength, maskBin);

    /* apply mask to omegaGW */
    for (i = 0; i < filterLength; i++)
      omegaGW->data->data[i] *= mask->data->data[i];

    /* destroy frequency mask */
    XLALDestroyREAL4FrequencySeries(mask);
  }

  /* perform search */
  stochHead = stochastic_search(&status, seriesOne, seriesTwo, overlap, \
      omegaGW, dataWindow, numSegments, filterLength, segsInInt, \
      segmentLength);

  /* save out xml table */
  save_xml_file(&status, accuracy, outputPath, baseName, stochHead);

  /* cleanup */
  XLALDestroyREAL4FrequencySeries(overlap);
  XLALDestroyREAL4FrequencySeries(omegaGW);
  XLALDestroyREAL4Window(dataWindow);
  XLALDestroyREAL4TimeSeries(seriesOne);
  XLALDestroyREAL4TimeSeries(seriesTwo);

  /* free memory used in the stochastic xml table */
  while(stochHead)
  {
    thisStoch = stochHead;
		stochHead = stochHead->next;
		LALFree(thisStoch);
	}

  /* free calloc'd memory */
  if (strcmp(frameCacheOne, frameCacheTwo))
  {
    free(frameCacheOne);
    free(frameCacheTwo);
  }
  else
  {
    free(frameCacheOne);
  }
  free(calCacheOne);
  free(calCacheTwo);
  free(channelOne);
  free(channelTwo);
  free(ifoOne);
  free(ifoTwo);
  free(userTag);
  free(outputPath);

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  exit(0);
}

/*
 * vim: et
 */
