/*
 * stochastic.c - SGWB Standalone Analysis Pipeline
 *
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 * Tania Regimbau <Tania.Regimbau@astro.cf.ac.uk>
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

/* C99 prototypes */
double round(double x);

NRCSID(STOCHASTICC, "$Id$");
RCSID("$Id$");

/* cvs info */
#define PROGRAM_NAME "stochastic"
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

/* window duration for psd estimation */
#define PSD_WINDOW_DURATION 4

/* system error checking */
extern int errno;

/* variables for getopt options parsing */
extern char *optarg;
extern int optind;

/* flags for getopt_long */
static int middle_segment_flag;
static int inject_flag;
static int apply_mask_flag;
static int high_pass_flag;
static int overlap_hann_flag;
static int recentre_flag;
static int cc_spectra_flag;
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

/* monte carlo parameters */
REAL4 scaleFactor = -1;
INT4 seed = -1;
INT4 NLoop = 1;

/* window parameters */
INT4 hannDuration = -1;

/* high pass filtering parameters */
REAL4 highPassFreq = -1;
REAL4 highPassAtten = -1;
INT4  highPassOrder = -1;

/* GEO scale factor */
REAL4 geoScaleFactor = 1e18;

/* GEO high pass filter parameters */
REAL4 geoHighPassFreq = -1; /* 70; */
INT4  geoHighPassOrder = -1; /* 8; */
REAL4 geoHighPassAtten = -1; /* 0.9; */

/* number of bins for frequency masking */
INT4 maskBin = -1;

/* output file */
CHAR *outputFilePath = NULL;

/* helper functions */

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
  LAL_CALL(LALButterworthREAL8TimeSeries(status, geo, &high_pass_params), \
      status);

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
  ResampleTSParams resample_params;
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

    /* set resample parameters */
    resample_params.deltaT = 1./resampleRate;
    resample_params.filterType = defaultButterworth;

    /* resample */
    LAL_CALL(LALResampleREAL4TimeSeries(status, series, &resample_params), \
        status);
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
    LAL_CALL(LALButterworthREAL4TimeSeries(status, series, \
          &high_pass_params), status);
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
    REAL8 deltaF,
    LIGOTimeGPS gps_time)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOmegaGWParameters omega_params;

  /* create and initialise frequency series */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "OmegaGW", \
        gps_time, f0, deltaF, lalDimensionlessUnit, length), status);

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
    INT4 site_two,
    LIGOTimeGPS gps_time)
{
  /* variables */
  REAL4FrequencySeries *series;
  OverlapReductionFunctionParameters overlap_params;
  LALDetectorPair detectors;

  /* create and initialise frequency series */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "Overlap", \
        gps_time, f0, deltaF, lalDimensionlessUnit, length), status);

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

/* helper function to increment the seconds component of a gps time with
 * an integer */
static LIGOTimeGPS increment_gps(LALStatus *status,
    LIGOTimeGPS gps_time,
    INT4 increment)
{
  /* variables */
  LIGOTimeGPS result;
  LALTimeInterval interval;

  /* set interval */
  interval.seconds = increment;
  interval.nanoSeconds = 0;

  /* increment GPS time */
  LAL_CALL(LALIncrementGPS(status, &result, &gps_time, &interval), status);

  return(result);
}

/* function to cut a time series between given start and end times */
static REAL4TimeSeries *cut_time_series(LALStatus *status,
    REAL4TimeSeries *input,
    LIGOTimeGPS start,
    LIGOTimeGPS end)
{
  /* variables */
  REAL4TimeSeries *series;
  INT4 length;
  INT4 first;

  /* calculate length of segment to cut */
  length = floor((XLALDeltaFloatGPS(&end, &start) / input->deltaT) + 0.5);

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
    LIGOTimeGPS gps_time,
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
  LALSnprintf(fname, sizeof(fname), "%s-ccSpectra-%d-%d.gwf", source, \
      gps_time.gpsSeconds, duration);

  /* setup frame file */
  frfile = FrFileONew(fname, 0);

  /* set frame properties */
  frame = FrameHNew(source);
  frame->run = 0;
  frame->frame = 0;
  frame->GTimeS = gps_time.gpsSeconds;
  frame->GTimeN = gps_time.gpsNanoSeconds;
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
    REAL4FrequencySeries *psdOne,
    REAL4FrequencySeries *psdTwo,
    REAL4WithUnits normalisation)
{
  /* variables */
  REAL4FrequencySeries *series;
  StochasticOptimalFilterCalInput input;

  /* allocate memory */
  LAL_CALL(LALCreateREAL4FrequencySeries(status, &series, "filter", \
        psdOne->epoch, psdOne->f0, psdOne->deltaF, lalDimensionlessUnit, \
        psdOne->data->length), status);

  /* set input */
  input.overlapReductionFunction = overlap;
  input.omegaGW = omega;
  input.calibratedInverseNoisePSD1 = psdOne;
  input.calibratedInverseNoisePSD2 = psdTwo;

  /* generate optimal filter */
  LAL_CALL(LALStochasticOptimalFilterCal(status, series, &input, \
        &normalisation), status);

  return(series);
}

/* wrapper function for estimating the psd */
static REAL4FrequencySeries *estimate_psd(LALStatus *status,
    REAL4TimeSeries *series)
{
  /* variables */
  REAL4FrequencySeries *psd;
  REAL4Window *window;
  RealFFTPlan *plan = NULL;
  AverageSpectrumParams psd_params;
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

  /* create window for psd estimation */
  window = XLALCreateHannREAL4Window(length);

  /* create fft plan for psd estimation */
  plan = XLALCreateForwardREAL4FFTPlan(length, 0);

  /* set parameters */
  psd_params.window = window;
  psd_params.overlap = overlap;
  psd_params.method = useMean;
  psd_params.plan = plan;

  /* esimate psd */
  LAL_CALL(LALREAL4AverageSpectrum(status, psd, series, &psd_params), status);

  /* destroy fft plan */
  XLALDestroyREAL4FFTPlan(plan);  

  /* free memory for window */
  XLALDestroyREAL4Window(window);

  return(psd);
}  

/* display usage information */
static void display_usage()
{
  fprintf(stdout, "Usage: " PROGRAM_NAME " [options]\n");
  fprintf(stdout, " --help                        print this message\n");
  fprintf(stdout, " --version                     display version\n");
  fprintf(stdout, " --verbose                     verbose mode\n");
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
  fprintf(stdout, " --inject                      inject a signal into the data\n");
  fprintf(stdout, " --scale-factor N              scale factor for injection\n");
  fprintf(stdout, " --seed N                      random seed\n");
  fprintf(stdout, " --trials N                    number of trials for Monte Carlo\n");
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
      {"inject", no_argument, &inject_flag, 1},
      {"middle-segment", no_argument, &middle_segment_flag, 1},
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"verbose", no_argument, &vrbflg, 1},
      {"recentre", no_argument, &recentre_flag, 1},
      {"cc-spectra", no_argument, &cc_spectra_flag, 1},
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
      {"scale-factor", required_argument, 0, 'B'},
      {"seed", required_argument, 0, 'C'},
      {"trials", required_argument, 0, 'D'},
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
        "A:B:C:D:E:F:G:H:I:J:", long_options, &option_index);

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
        outputFilePath = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(outputFilePath, optarg, optarg_len);
        if ((stat(outputFilePath, &fileStatus) == -1) && (errno = ENOENT))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Directory does not exist: (%s specified)\n", \
              long_options[option_index].name, outputFilePath);
          exit(1);
        }
        ADD_PROCESS_PARAM("string", "%s", outputFilePath);
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
        if (fMin != round(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMin = round(fMin * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
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
        if (fMax != round(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION)
        {
          fMax = round(fMax * PSD_WINDOW_DURATION) / PSD_WINDOW_DURATION;
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
        if ((stat(frameCacheOne, &fileStatus) == -1) && (errno = ENOENT))
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
        if ((stat(frameCacheTwo, &fileStatus) == -1) && (errno = ENOENT))
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
        if ((stat(calCacheOne, &fileStatus) == -1) && (errno = ENOENT))
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
        if ((stat(calCacheTwo, &fileStatus) == -1) && (errno = ENOENT))
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

      case 'B':
        /* scale factor */
        scaleFactor = atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", scaleFactor);
        break;

      case 'C':
        /* seed */
        seed = atoi(optarg);
        ADD_PROCESS_PARAM("float", "%e", seed);
        break;

      case 'D':
        /* number of trials */
        NLoop = atoi(optarg);
        if (NLoop <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Number of trials to must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name, NLoop);
          exit(1);
        }
        ADD_PROCESS_PARAM("int", "%d", NLoop);
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

  /* injections */
  if (inject_flag)
  {
    if (scaleFactor == -1)
    {
      fprintf(stderr, "--scale-factor must be specified\n");
      exit(1);
    }
    if (seed == -1)
    {
      fprintf(stderr, "--seed must be specified\n");
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

/* program entry point */
INT4 main(INT4 argc, CHAR *argv[])
{
  /* lal initialisation variables */
  LALStatus status = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;

  /* xml */
  CHAR xmlFileName[FILENAME_MAX];
  LIGOLwXMLStream xmlStream;
  StochasticTable *stochHead = NULL;
  StochasticTable *thisStoch = NULL;
  MetadataTable outputTable;

  /* counters */
  INT4 i, j, segLoop, interLoop;

  /* results parameters */
  REAL8 y;
  REAL8 varTheo;

  /* input data */
  LIGOTimeGPS gpsStartTime;
  LIGOTimeGPS gpsEndTime;
  LIGOTimeGPS gpsSegStartTime;
  LIGOTimeGPS gpsSegEndTime;
  LIGOTimeGPS gpsAnalysisTime;
  REAL4TimeSeries *seriesOne;
  REAL4TimeSeries *seriesTwo;

  /* input data segment */
  INT4 numSegments;
  INT4 duration, durationEff, extrasec;
  INT4 segsInInt, numIntervals, segMiddle;
  INT4 segmentLength, intervalLength;
  INT4 segmentShift;
  INT4 padData;
  LIGOTimeGPS gpsCalibTime;
  REAL4TimeSeries *segmentOne;
  REAL4TimeSeries *segmentTwo;
  REAL4Vector *segOne[100], *segTwo[100];

  /* simulated signal structures */
  SSSimStochBGParams SBParams;
  SSSimStochBGInput SBInput;
  SSSimStochBGOutput SBOutput;
  REAL4TimeSeries *SimStochBGOne;
  REAL4TimeSeries *SimStochBGTwo;
  REAL4FrequencySeries *MComegaGW = NULL;
  COMPLEX8FrequencySeries *MCresponseOne;
  COMPLEX8FrequencySeries *MCresponseTwo;
  COMPLEX8Vector *MCrespOne[100], *MCrespTwo[100];
  INT4 MCLoop;
  INT4 MCfreqLength = 0;
  REAL8 MCdeltaF, MCdeltaT;
  LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* window for segment data streams */
  REAL4Window *dataWindow;

  /* response functions */
  COMPLEX8FrequencySeries *responseTempOne;
  COMPLEX8FrequencySeries *responseTempTwo;
  COMPLEX8FrequencySeries *responseOne;
  COMPLEX8FrequencySeries *responseTwo;
  COMPLEX8Vector *respOne[100], *respTwo[100];
  INT4 respLength;
  CalibrationUpdateParams calfacts;
  LALUnit countPerAttoStrain = {18,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};
  FrCache *calibCache = NULL;

  /* data structures for PSDs */
  REAL8 deltaF;
  INT4 filterLength;
  INT4 numFMin;
  INT4 numFMax;
  REAL4FrequencySeries *psdTempOne = NULL;
  REAL4FrequencySeries *psdTempTwo = NULL;
  REAL4FrequencySeries *psdOne;
  REAL4FrequencySeries *psdTwo;
  REAL4Vector *calPsdOne, *calPsdTwo;
  LALUnit psdUnits = {0,{0,0,1,0,0,0,2},{0,0,0,0,0,0,0}};

  /* calibrated inverse noise data structures */
  REAL4FrequencySeries *calInvPsdOne = NULL;
  REAL4FrequencySeries *calInvPsdTwo = NULL;

  /* zeropad and fft structures */
  SZeroPadAndFFTParameters zeroPadParams;
  RealFFTPlan *fftDataPlan = NULL;
  COMPLEX8FrequencySeries *hBarTildeOne;
  COMPLEX8FrequencySeries *hBarTildeTwo;
  UINT4 zeroPadLength;
  UINT4 fftDataLength;

  /* overlap reduction function */
  REAL4FrequencySeries *overlap;

  /* frequency mask structures */
  REAL4FrequencySeries *mask;
  INT4 Nbin;

  /* structures for optimal filter normalisation */
  StochasticOptimalFilterNormalizationInput normInput;
  StochasticOptimalFilterNormalizationOutput normOutput;
  StochasticOptimalFilterNormalizationParameters normParams;
  REAL4WithUnits normLambda;
  REAL4WithUnits normSigma;
  REAL8 lambda;

  /* structures for optimal filter */
  REAL4FrequencySeries *optFilter = NULL;

  /* spectrum structures */
  REAL4FrequencySeries *omegaGW;

  /* structures for CC spectrum and CC statistics */
  StochasticCrossCorrelationCalInput ccIn;
  BOOLEAN epochsMatch = 1;
  REAL4WithUnits ccStat;
  COMPLEX8FrequencySeries *ccSpectrum;

  /* error handler */
  status.statusPtr = NULL;
  lal_errhandler = LAL_ERR_EXIT;

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

  /* get xml file name */
  if (userTag)
  {
    LALSnprintf(xmlFileName, FILENAME_MAX, "%s%s-stochastic_%s_%d-%d.xml", \
        ifoOne, ifoTwo, userTag, startTime, endTime);
  }
  else
  {
    LALSnprintf(xmlFileName, FILENAME_MAX, "%s%s-stochastic-%d-%d.xml", \
        ifoOne, ifoTwo, startTime, endTime);
  }

  /* get number of segments */
  duration = endTime - startTime;
  numSegments = duration / segmentDuration;
  segsInInt = intervalDuration / segmentDuration;
  segMiddle = (segsInInt - 1) / 2;
  
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

  /* calculate number of intervals, and required shift to get to next
   * interval */
  if (overlap_hann_flag)
  {
    numIntervals = (2 * (numSegments - 2)) - 1;
    segmentShift = segmentDuration / 2;
  }
  else
  {
    numIntervals = numSegments - 2;
    segmentShift = segmentDuration;
  }

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
    fprintf(stderr, "series have different sample rates...\n");
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

  /* initialize calibration gps time structure */
  gpsCalibTime.gpsSeconds = startTime + calibOffset;
  gpsCalibTime.gpsNanoSeconds = 0;

  /* set length for data segments */
  intervalLength = intervalDuration * resampleRate;
  segmentLength = segmentDuration * resampleRate;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for data segments...\n");

  /* allocate memory for data segments */
  LAL_CALL(LALCreateREAL4TimeSeries(&status, &segmentOne, "segmentOne", \
        gpsStartTime, 0, 1./resampleRate, lalDimensionlessUnit, \
        segmentLength), &status);
  LAL_CALL(LALCreateREAL4TimeSeries(&status, &segmentTwo, "segmentTwo", \
        gpsStartTime, 0, 1./resampleRate, lalDimensionlessUnit, \
        segmentLength), &status);

  /* allocate memory for temporary segments - to be re-worked */
  for (i = 0; i < segsInInt; i++)
  {
    segOne[i]= segTwo[i] = NULL;
    LAL_CALL(LALCreateVector(&status, &(segOne[i]), segmentLength), &status);
    LAL_CALL(LALCreateVector(&status, &(segTwo[i]), segmentLength), &status);
    memset(segOne[i]->data, 0, segOne[i]->length * sizeof(*segOne[i]->data));
    memset(segTwo[i]->data, 0, segTwo[i]->length * sizeof(*segTwo[i]->data));
  }

  /* allocate memory for injections - to be re-worked */
  if (inject_flag)
  {
    if (vrbflg)
      fprintf(stdout, "Allocating memory for MC...\n");

    /* set  monte carlo parameters */
    MCdeltaT = 1.0 / resampleRate;
    MCdeltaF = (REAL8)resampleRate / (REAL8)segmentLength;
    MCfreqLength = (segmentLength / 2) + 1;

    /* create vectors to store the simulated signal */
    LAL_CALL(LALCreateREAL4TimeSeries(&status, &SimStochBGOne, \
          "Whitened-SimulatedSB1", gpsStartTime, 0, 1./resampleRate, \
          lalDimensionlessUnit, segmentLength), &status);
    LAL_CALL(LALCreateREAL4TimeSeries(&status, &SimStochBGTwo, \
          "Whitened-SimulatedSB2", gpsStartTime, 0, 1./resampleRate, \
          lalDimensionlessUnit, segmentLength), &status);

    /* generate omegaGW */
    MComegaGW = omega_gw(&status, alpha, fRef, omegaRef, MCfreqLength, 0, \
        MCdeltaF, gpsStartTime);

    /* response functions */
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &MCresponseOne, \
          "MCresponseOne", gpsCalibTime, 0, MCdeltaF, countPerStrain, \
          MCfreqLength), &status);
    LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &MCresponseTwo, \
          "MCresponseTwo", gpsCalibTime, 0, MCdeltaF, countPerStrain, \
          MCfreqLength), &status);

    /* allocate memory for temporary monte carlo segments */
    for (i = 0; i < segsInInt; i++)
    {
      MCrespOne[i]= MCrespTwo[i] = NULL;
      LAL_CALL(LALCCreateVector(&status, &(MCrespOne[i]), MCfreqLength), \
          &status);
      LAL_CALL(LALCCreateVector(&status, &(MCrespTwo[i]), MCfreqLength), \
          &status);
      memset(MCrespOne[i]->data, 0, \
          MCrespOne[i]->length * sizeof(*MCrespOne[i]->data));
      memset(MCrespTwo[i]->data, 0, \
          MCrespTwo[i]->length * sizeof(*MCrespTwo[i]->data));
    }
  }

  /* get bins for min and max frequencies */
  numFMin = (INT4)(fMin / deltaF);
  numFMax = (INT4)(fMax / deltaF);

  /* get lengths */
  filterLength = numFMax - numFMin + 1;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for reduced frequency band PSDs...\n");

  /* allocate memory for reduced frequency band PSDs */
  LAL_CALL(LALCreateREAL4FrequencySeries(&status, &psdOne, "psdOne", \
        gpsStartTime, fMin, deltaF, psdUnits, filterLength), &status);
  LAL_CALL(LALCreateREAL4FrequencySeries(&status, &psdTwo, "psdTwo", \
        gpsStartTime, fMin, deltaF, psdUnits, filterLength), &status);

  /* allocate memory for calibrated PSDs */
  calPsdOne = calPsdTwo = NULL;
  LAL_CALL(LALCreateVector(&status, &calPsdOne, filterLength), &status);
  LAL_CALL(LALCreateVector(&status, &calPsdTwo, filterLength), &status);
  memset(calPsdOne->data, 0, calPsdOne->length * sizeof(*calPsdOne->data));
  memset(calPsdTwo->data, 0, calPsdTwo->length * sizeof(*calPsdTwo->data));

  /* set parameters for response functions */
  respLength = (UINT4)(fMax / deltaF) + 1;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for response functions...\n");

  /* allocate memory for response functions */
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &responseTempOne, \
        "responseTempOne", gpsCalibTime, 0, deltaF, countPerAttoStrain, \
        respLength), &status);
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &responseTempTwo, \
        "responseTempTwo", gpsCalibTime, 0, deltaF, countPerAttoStrain, \
        respLength), &status);

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for reduced frequency band " \
        "response functions...\n");
  }

  /* allocate memory for reduced frequency band response functions */
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &responseOne, \
        "responseOne", gpsCalibTime, fMin, deltaF, countPerAttoStrain, \
        filterLength), &status);
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &responseTwo, \
        "responseTwo", gpsCalibTime, fMin, deltaF, countPerAttoStrain, \
        filterLength), &status);

  /* allocate memory for temporary response functions */
  for (i = 0; i < segsInInt; i++)
  {
    respOne[i]= respTwo[i] = NULL;
    LAL_CALL(LALCCreateVector(&status, &(respOne[i]), filterLength), &status);
    LAL_CALL(LALCCreateVector(&status, &(respTwo[i]), filterLength), &status);
    memset(respOne[i]->data, 0, respOne[i]->length * \
        sizeof(*respOne[i]->data));
    memset(respTwo[i]->data, 0, respTwo[i]->length * \
        sizeof(*respTwo[i]->data));
  }

  if (vrbflg)
    fprintf(stdout, "Generating data segment window...\n");

  /* for overlapping hann windows, the hann window length is the segment
   * length */
  if (overlap_hann_flag)
    hannDuration = segmentDuration;

  /* create window for data */
  dataWindow = data_window(seriesOne->deltaT, segmentLength, hannDuration);

  /* zeropad lengths */
  zeroPadLength = 2 * segmentLength;
  fftDataLength = (zeroPadLength / 2) + 1;

  /* create fft plan */
  fftDataPlan = XLALCreateForwardREAL4FFTPlan(zeroPadLength, 0);

  if (vrbflg)
    fprintf(stdout, "Allocating memory for zeropad...\n");

  /* allocate memory for zeropad */
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &hBarTildeOne, \
        "hBarTildeOne", gpsStartTime, 0, deltaF, lalDimensionlessUnit, \
        fftDataLength), &status);
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &hBarTildeTwo, \
        "hBarTildeTwo", gpsStartTime, 0, deltaF, lalDimensionlessUnit, \
        fftDataLength), &status);

  /* quantities needed to build the optimal filter */

  /* generate overlap reduction function */
  if (vrbflg)
    fprintf(stdout, "Generating the overlap reduction function...\n");
  overlap = overlap_reduction_function(&status, filterLength, fMin, deltaF, \
      siteOne, siteTwo, gpsStartTime);

  /* generage omegaGW */
  if (vrbflg)
    fprintf(stdout, "Generating spectrum for optimal filter...\n");
  omegaGW = omega_gw(&status, alpha, fRef, omegaRef, filterLength, \
      fMin, deltaF, gpsStartTime);

  /* frequency mask */
  if (apply_mask_flag)
  {
    /* extra bins */
    Nbin = (maskBin - 1) / 2;

    if (vrbflg)
      fprintf(stdout, "Allocating memory for frequency mask...\n");

    /* allocate memory for frequency mask */
    LAL_CALL(LALCreateREAL4FrequencySeries(&status, &mask, \
          "mask", gpsStartTime, fMin, deltaF, lalDimensionlessUnit, \
          respLength), &status);

    if (vrbflg)
      fprintf(stdout, "Generating frequency mask...\n");

    /* set all values to 1 */
    for (i = 0; i < respLength; i++)
      mask->data->data[i] = 1;

    if (vrbflg)
      fprintf(stdout, "Masking multiples of 16 Hz...\n");

    /* remove multiples of 16 Hz */
    for (i = 0; i < respLength; i += (UINT4)(16 / deltaF))
    {
      mask->data->data[i]= 0;

      for (j = 0; j < Nbin; j++)
      {
        if ((i + 1 + j) < respLength)
          mask->data->data[i + 1 + j]= 0;
        if ((i - 1 - j) > 0 )
          mask->data->data[i - 1 - j]= 0;
      }
    }

    if (vrbflg)
      fprintf(stdout, "Masking multiples of 60 Hz...\n");

    /* remove multiples of 60 Hz */
    for (i = 0; i < respLength; i += (UINT4)(60 / deltaF))
    {
      mask->data->data[i] = 0;

      for (j = 0; j < Nbin; j++)
      {
        if ((i + 1 + j) < respLength)
          mask->data->data[i + 1 + j]= 0;
        if ((i - 1 - j) > 0 )
          mask->data->data[i - 1 - j]= 0;
      }
    }

    if (vrbflg)
      fprintf(stdout, "Getting appropriate frequency band for mask...\n");

    /* get appropriate band */
    LAL_CALL(LALShrinkREAL4FrequencySeries(&status, mask, numFMin, \
          filterLength), &status);

    if (vrbflg)
      fprintf(stdout, "Applying frequency mask to spectrum..\n");

    /* apply mask to omegaGW */
    for (i = 0; i < filterLength; i++)
      omegaGW->data->data[i] *= mask->data->data[i];

    /* destroy frequency mask */
    XLALDestroyREAL4FrequencySeries(mask);
  }

  if (vrbflg)
    fprintf(stdout, "Allocating memory for CC Spectrum...\n");

  /* allocate memory for CC spectrum*/
  LAL_CALL(LALCreateCOMPLEX8FrequencySeries(&status, &ccSpectrum, \
        "ccSpectrum", gpsStartTime, fMin, deltaF, lalDimensionlessUnit, \
        filterLength), &status);

  if (vrbflg)
    fprintf(stdout, "Done with memory allocation...\n");

  /* main analysis loop */
  for (MCLoop = 0; MCLoop < NLoop; MCLoop++)
  {
    /* loop over intervals */
    for (interLoop = 0; interLoop < numIntervals; interLoop++)
    {	
      /* loop over segments */
      for (segLoop = 0; segLoop < segsInInt; segLoop++)
      {
        /* get segment start/end time */
        gpsSegStartTime = increment_gps(&status, gpsStartTime, \
            (interLoop * segmentShift) + (segLoop * segmentDuration));
        gpsSegEndTime = increment_gps(&status, gpsSegStartTime, \
            segmentDuration);
        segmentOne->epoch = gpsSegStartTime;
        segmentTwo->epoch = gpsSegStartTime;

        /* is this the analysis segment */
        if (segLoop == segMiddle)
        {
          gpsAnalysisTime = gpsSegStartTime;
        }

        if (vrbflg)
        {
          fprintf(stdout, "request data at GPS time %d\n", \
              gpsSegStartTime.gpsSeconds);
        }

        /* cut segments from series */
        segmentOne = cut_time_series(&status, seriesOne, gpsSegStartTime, \
            gpsSegEndTime);
        segmentTwo = cut_time_series(&status, seriesTwo, gpsSegStartTime, \
            gpsSegEndTime);

        /* store in memory */
        for (i = 0; i < segmentLength; i++)
        {
          segOne[segLoop]->data[i] = segmentOne->data->data[i];
          segTwo[segLoop]->data[i] = segmentTwo->data->data[i];
        }

        /* compute response functions */
        gpsCalibTime.gpsSeconds = gpsSegStartTime.gpsSeconds + calibOffset;
        responseOne->epoch = gpsCalibTime;
        responseTwo->epoch = gpsCalibTime;

        /* set response function to unity for GEO */
        if (strncmp(ifoOne, "G1", 2) == 0)
        {
          for (i = 0; i < filterLength; i++)
          {
            responseOne->data->data[i].re = 1;
            responseOne->data->data[i].im = 0;
          }
        }
        else
        {
          /* generate response function */
          responseTempOne->epoch = gpsCalibTime;
          if (vrbflg)
          {
            fprintf(stdout, "Getting response function for %s at %d\n", \
                ifoOne, gpsCalibTime.gpsSeconds);
          }
          memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
          calfacts.ifo = ifoOne;
          calibCache = XLALFrImportCache(calCacheOne);
          LAL_CALL(LALExtractFrameResponse(&status, responseTempOne, \
                calibCache, &calfacts), &status);
          XLALFrDestroyCache(calibCache);

          /* reduce to the optimal filter frequency range */
          for (i = 0; i < filterLength; i++)
          {
            responseOne->data->data[i] = responseTempOne->data->data[i + \
                                         numFMin];
          }
        }

        /* set response function to unity for GEO */
        if (strncmp(ifoTwo, "G1", 2) == 0)
        {
          for (i = 0; i < filterLength; i++)
          {
            responseTwo->data->data[i].re = 1;
            responseTwo->data->data[i].im = 0;
          }
        }
        else
        {
          /* generate response function */
          responseTempTwo->epoch = gpsCalibTime;
          if (vrbflg)
          {
            fprintf(stdout, "Getting response function for %s at %d\n", \
                ifoOne, gpsCalibTime.gpsSeconds);
          }
          memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
          calfacts.ifo = ifoTwo;
          calibCache = XLALFrImportCache(calCacheTwo);
          LAL_CALL(LALExtractFrameResponse(&status, responseTempTwo, \
                calibCache, &calfacts), &status);
          XLALFrDestroyCache(calibCache);

          /* reduce to the optimal filter frequency range */
          for (i = 0; i < filterLength; i++)
          {
            responseTwo->data->data[i] = responseTempTwo->data->data[i + \
                                         numFMin];
          }
        }

        /* store in memory */
        for (i = 0; i < filterLength; i++)
        {
          respOne[segLoop]->data[i] = responseOne->data->data[i];
          respTwo[segLoop]->data[i] = responseTwo->data->data[i];
        }

        if (inject_flag)
        {
          /* convert response function for use in the MC routine */
          MCresponseOne->epoch = gpsCalibTime;
          MCresponseTwo->epoch = gpsCalibTime;
          LAL_CALL(LALResponseConvert(&status, MCresponseOne, \
                responseTempOne), &status);
          LAL_CALL(LALResponseConvert(&status, MCresponseTwo, \
                responseTempTwo), &status);

          /* force DC to be 0 and nyquist to be real */
          MCresponseOne->data->data[0].re = 0;
          MCresponseTwo->data->data[0].re = 0;
          MCresponseOne->data->data[0].im = 0;
          MCresponseTwo->data->data[0].im = 0;
          MCresponseOne->data->data[MCfreqLength-1].im = 0;
          MCresponseTwo->data->data[MCfreqLength-1].im = 0;

          /* store in memory */
          for (i = 0; i < MCfreqLength; i++)
          {
            MCrespOne[segLoop]->data[i] = MCresponseOne->data->data[i];
            MCrespTwo[segLoop]->data[i] = MCresponseTwo->data->data[i];
          }
        }
      }

      /* initialize average PSDs */
      for (i = 0; i < filterLength; i++)
      {
        calPsdOne->data[i] = 0;
        calPsdTwo->data[i] = 0;
      }

      for (segLoop = 0; segLoop < segsInInt; segLoop++)
      {
        /* set segment start and end time */
        gpsSegStartTime.gpsSeconds = gpsStartTime.gpsSeconds + \
                                     (interLoop * intervalDuration) + \
                                     (segLoop * segmentDuration);
        gpsSegStartTime.gpsNanoSeconds = 0;
        gpsSegEndTime.gpsSeconds = gpsSegStartTime.gpsSeconds + \
                                   segmentDuration;
        segmentOne->epoch = gpsSegStartTime;
        segmentTwo->epoch = gpsSegStartTime;
        gpsCalibTime.gpsSeconds = gpsSegStartTime.gpsSeconds + calibOffset;
        responseOne->epoch = gpsCalibTime;
        responseTwo->epoch = gpsCalibTime;

        /* copy to temporary storage */
        for (i = 0; i < filterLength; i++)
        {
          responseOne->data->data[i] = respOne[segLoop]->data[i];
          responseTwo->data->data[i] = respTwo[segLoop]->data[i];
        }

        /* simulate signal */
        if (inject_flag)
        {
          /* copy to temporay storage */
          for (i = 0; i < MCfreqLength; i++)
          {
            MCresponseOne->data->data[i] = MCrespOne[segLoop]->data[i];
            MCresponseTwo->data->data[i] = MCrespTwo[segLoop]->data[i];
          }

          /* set parameters for monte carlo */
          SBParams.length = segmentLength;
          SBParams.deltaT = 1. / resampleRate;
          SBParams.detectorOne = lalCachedDetectors[siteOne];
          SBParams.detectorTwo = lalCachedDetectors[siteTwo];
          SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
          SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;
          SBParams.seed = seed;

          /* define input structure for SimulateSB */
          SBInput.omegaGW = MComegaGW;
          SBInput.whiteningFilter1 = MCresponseOne;
          SBInput.whiteningFilter2 = MCresponseTwo;

          /* define output structure for SimulateSB */
          SBOutput.SSimStochBG1 = SimStochBGOne;
          SBOutput.SSimStochBG2 = SimStochBGTwo;

          /* set epochs for monte carlo outputs */
          SimStochBGOne->epoch = gpsSegStartTime;
          SimStochBGTwo->epoch = gpsSegStartTime;

          /* perform monte carlo */
          LAL_CALL(LALSSSimStochBGTimeSeries(&status, &SBOutput, \
                &SBInput, &SBParams), &status);

          /* multiply by scale factor and inject into real data */
          for (i = 0; i < segmentLength ; i++)
          {
            segmentOne->data->data[i] = segOne[segLoop]->data[i] + \
                                        (scaleFactor * \
                                         SimStochBGOne->data->data[i]);
            segmentTwo->data->data[i] = segTwo[segLoop]->data[i] + \
                                        (scaleFactor * \
                                         SimStochBGTwo->data->data[i]);
          }

          /* increase seed */
          seed = seed + 2;
        }
        else
        {
          for (i = 0; i < segmentLength; i++)
          {
            segmentOne->data->data[i] = segOne[segLoop]->data[i];
            segmentTwo->data->data[i] = segTwo[segLoop]->data[i];
          }
        }

        /* store in memory */
        for (i = 0; i < segmentLength; i++)
        {
          segOne[segLoop]->data[i] = segmentOne->data->data[i];
          segTwo[segLoop]->data[i] = segmentTwo->data->data[i];
        }

        /* check if on middle segment and if we want to include this in
         * the analysis */
        if ((segLoop == segMiddle) && (middle_segment_flag == 0))
        {
          if (vrbflg)
            fprintf(stdout, "Ignoring middle segment..\n");
        }
        else
        {
          if (vrbflg)
            fprintf(stdout, "Estimating PSDs...\n");

          /* compute uncalibrated PSDs */
          psdTempOne = estimate_psd(&status, segmentOne);
          psdTempTwo = estimate_psd(&status, segmentTwo);

          if (vrbflg)
          {
            fprintf(stdout, "Getting appropriate frequency band for " \
                "PSDs..\n");
          }

          /* reduce to the optimal filter frequency range */
          for (i = 0; i < filterLength; i++)
          {
            psdOne->data->data[i] =  psdTempOne->data->data[i + numFMin];
            psdTwo->data->data[i] =  psdTempTwo->data->data[i + numFMin];
          }

          if (vrbflg)
            fprintf(stdout, "Generating inverse noise...\n");

          /* compute inverse calibrate PSDs */
          calInvPsdOne = inverse_noise(&status, psdOne, responseOne);
          calInvPsdTwo = inverse_noise(&status, psdTwo, responseTwo);

          /* sum over calibrated PSDs for average */
          for (i = 0; i < filterLength; i++)
          {
            calPsdOne->data[i] = calPsdOne->data[i] + \
                                 1. / calInvPsdOne->data->data[i];
            calPsdTwo->data[i] = calPsdTwo->data[i] + \
                                 1. / calInvPsdTwo->data->data[i];
          }
        }
      }

      /* average calibrated PSDs and take inverse */
      for (i = 0; i < filterLength; i++)
      {
        /* average */
        if (middle_segment_flag == 0)
        {
          calPsdOne->data[i] = calPsdOne->data[i] / (REAL4)(segsInInt - 1);
          calPsdTwo->data[i] = calPsdTwo->data[i] / (REAL4)(segsInInt - 1);
        }
        else
        {
          calPsdOne->data[i] = calPsdOne->data[i] / (REAL4)segsInInt;
          calPsdTwo->data[i] = calPsdTwo->data[i] / (REAL4)segsInInt;
        }
        /* take inverse */
        calInvPsdOne->data->data[i] = 1. / calPsdOne->data[i];
        calInvPsdTwo->data->data[i] = 1. / calPsdTwo->data[i];
      }

      /* set normalisation parameters */
      normParams.fRef = fRef;
      normParams.heterodyned = 0;
      normParams.window1 = dataWindow->data;
      normParams.window2 = dataWindow->data;

      /* set normalisation input */
      normInput.overlapReductionFunction = overlap;
      normInput.omegaGW = omegaGW;
      normInput.inverseNoisePSD1 = calInvPsdOne;
      normInput.inverseNoisePSD2 = calInvPsdTwo;

      /* set normalisation output */
      normOutput.normalization = &normLambda;
      normOutput.variance = &normSigma;

      if (vrbflg)
        fprintf(stdout, "Normalising optimal filter...\n");

      /* compute variance and normalisation for optimal filter */
      LAL_CALL(LALStochasticOptimalFilterNormalization(&status, \
            &normOutput, &normInput, &normParams), &status);
      lambda = (REAL8)(normLambda.value * \
          pow(10.,normLambda.units.powerOfTen));
      varTheo = (REAL8)(segmentDuration * normSigma.value * \
          pow(10.,normSigma.units.powerOfTen));

      if (vrbflg)
        fprintf(stdout, "Generating optimal filter...\n");

      /* build optimal filter */
      optFilter = optimal_filter(&status, overlap, omegaGW, calInvPsdOne, \
          calInvPsdTwo, normLambda);

      if (vrbflg)
      {
        fprintf(stdout, "Analysing segment at GPS %d\n", \
            gpsAnalysisTime.gpsSeconds);
      }

      /* copy to temporary storage */
      for (i = 0; i < segmentLength; i++)
      {
        segmentOne->data->data[i] = segOne[segMiddle]->data[i];
        segmentTwo->data->data[i] = segTwo[segMiddle]->data[i];
      }

      /* set zeropad parameters */
      zeroPadParams.fftPlan = fftDataPlan;
      zeroPadParams.window = dataWindow->data;
      zeroPadParams.length = zeroPadLength;

      /* zero pad and fft */
      LAL_CALL(LALSZeroPadAndFFT(&status, hBarTildeOne, segmentOne, \
            &zeroPadParams), &status);
      LAL_CALL(LALSZeroPadAndFFT(&status, hBarTildeTwo, segmentTwo, \
            &zeroPadParams), &status);

      /* set CC inputs */
      ccIn.hBarTildeOne = hBarTildeOne;
      ccIn.hBarTildeTwo = hBarTildeTwo;
      ccIn.responseFunctionOne = responseOne;
      ccIn.responseFunctionTwo = responseTwo;
      ccIn.optimalFilter = optFilter;

      if (cc_spectra_flag)
      {
        if (vrbflg)
          fprintf(stdout, "Generating cross correlation spectrum...\n");

        /* calculate cc spectrum */
        LAL_CALL(LALStochasticCrossCorrelationSpectrumCal(&status, \
              ccSpectrum, &ccIn, epochsMatch), &status);

        /* save out cc spectra as frame */
        if (vrbflg)
        {
          fprintf(stdout, "Saving ccSpectra to frame...\n");
          write_ccspectra_frame(ccSpectrum, ifoOne, ifoTwo, \
              gpsAnalysisTime, segmentDuration);
        }
      }
      
      /* cc statistic */
      LAL_CALL(LALStochasticCrossCorrelationStatisticCal(&status, &ccStat, \
            &ccIn, epochsMatch), &status);
      y = (REAL8)(ccStat.value * pow(10., ccStat.units.powerOfTen));

      /* save */
      if (vrbflg)
      {
        fprintf(stdout, "Interval %d\n", interLoop + 1);
        fprintf(stdout, "  GPS time  = %d\n", gpsAnalysisTime.gpsSeconds);
        fprintf(stdout, "  y         = %e\n", y);
        fprintf(stdout, "  sigmaTheo = %e\n", sqrt(varTheo));
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
      thisStoch->start_time.gpsSeconds = gpsAnalysisTime.gpsSeconds;
      thisStoch->start_time.gpsNanoSeconds = gpsAnalysisTime.gpsNanoSeconds;
      thisStoch->duration.gpsSeconds = segmentDuration;
      thisStoch->duration.gpsNanoSeconds = 0;
      thisStoch->f_min = fMin;
      thisStoch->f_max = fMax;
      thisStoch->cc_stat = y;
      thisStoch->cc_sigma = sqrt(varTheo);
    }
  }

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
  if (inject_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--inject");
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

  /* write out xml */
  if (vrbflg)
    fprintf(stdout, "Writing output XML file...\n");

  /* opening xml file stream */
  memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
  LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlStream, xmlFileName), &status);

  /* write out process and process params tables */
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->end_time), \
        &accuracy), &status);
  LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, process_table), \
      &status);
  LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, proctable, \
        process_table), &status);
  LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);
  free(proctable.processTable);

  /* write the process params table */
  LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, \
        process_params_table), &status);
  LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, procparams, \
        process_params_table), &status);
  LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);

  /* write stochastic table */
  if (stochHead)
  {
    outputTable.stochasticTable = stochHead;
    LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, stochastic_table), \
        &status);;
    LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, outputTable, \
          stochastic_table), &status);
    LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);
  }

  /* close xml file */
  LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);

  /* cleanup */
  XLALDestroyREAL4TimeSeries(segmentOne);
  XLALDestroyREAL4TimeSeries(segmentTwo);
  XLALDestroyREAL4FFTPlan(fftDataPlan);
  XLALDestroyREAL4FrequencySeries(psdTempOne);
  XLALDestroyREAL4FrequencySeries(psdTempTwo);
  XLALDestroyREAL4FrequencySeries(psdOne);
  XLALDestroyREAL4FrequencySeries(psdTwo);
  XLALDestroyREAL4Vector(calPsdOne);
  XLALDestroyREAL4Vector(calPsdTwo);
  XLALDestroyCOMPLEX8FrequencySeries(responseTempOne);
  XLALDestroyCOMPLEX8FrequencySeries(responseTempTwo);
  XLALDestroyCOMPLEX8FrequencySeries(responseOne);
  XLALDestroyCOMPLEX8FrequencySeries(responseTwo);
  XLALDestroyREAL4FrequencySeries(optFilter);
  XLALDestroyREAL4FrequencySeries(calInvPsdOne);
  XLALDestroyREAL4FrequencySeries(calInvPsdTwo);
  XLALDestroyREAL4FrequencySeries(overlap);
  XLALDestroyREAL4FrequencySeries(omegaGW);
  XLALDestroyREAL4Window(dataWindow);
  XLALDestroyCOMPLEX8FrequencySeries(hBarTildeOne);
  XLALDestroyCOMPLEX8FrequencySeries(hBarTildeTwo);
  if (inject_flag)
  {
    XLALDestroyREAL4TimeSeries(SimStochBGOne);
    XLALDestroyREAL4TimeSeries(SimStochBGTwo);
    XLALDestroyCOMPLEX8FrequencySeries(MCresponseOne);
    XLALDestroyCOMPLEX8FrequencySeries(MCresponseTwo);
    XLALDestroyREAL4FrequencySeries(MComegaGW);
  }
  for (i = 0; i <segsInInt; i++)
  {
    XLALDestroyCOMPLEX8Vector(respOne[i]);
    XLALDestroyCOMPLEX8Vector(respTwo[i]);
    XLALDestroyREAL4Vector(segOne[i]);
    XLALDestroyREAL4Vector(segTwo[i]);
    if (inject_flag)
    {
      XLALDestroyCOMPLEX8Vector(MCrespOne[i]);
      XLALDestroyCOMPLEX8Vector(MCrespTwo[i]);
    }
  }

  /* free memory used in the stochastic xml table */
  while (stochHead)
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

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  exit(0);
}

/*
 * vim: et
 */
