/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include "config.h"
#include "detresponse.h"

/*
 * TODO: make file names depend on name of source, and of detector
 */

static const INT4 grid_lim = NUM_RA * NUM_DEC;
static const INT4 dec_lim  = (NUM_DEC - 1)/2;

static const double rad2deg = 180./LAL_PI;
static const double eps = 23.5;
static const double vorb = 2.*LAL_PI*1.5e11/(365.25*24.*3600.);
static const double clight = 2.998e8;

static double vorbrel;

static double relval(double ra, double dec, int i, int nrelvals)
{
    int j;
    double latdet,vrotrel,radearth;
    double phipinit,dphip,dphiptot,phiprad;
    double phirotinit,dphirot,dphirottot,phirotrad;
    double numdays,coseps,sineps;
    double unitvect[3];
    double vearth[3];
    double vdet[3];
    double vtotal[3];
    double retval = 0.;
    
    
    latdet = 46.5;
    radearth = 6.38e6;
    vrotrel = 2.*LAL_PI*radearth/(24.*3600.)/.997/clight*cos(latdet/rad2deg);
        
    numdays = 59.;
    phipinit = 180. - 34./365.25*360.;
    dphip = numdays/365.25*360./(nrelvals-1);
    phirotinit = -phipinit - 4./24.*360.;
    dphirot = numdays/.997*360./(nrelvals-1);
        
    coseps = cos(eps/rad2deg);
    sineps = sin(eps/rad2deg);
    
    unitvect[0] = cos(dec/rad2deg) * cos(ra/rad2deg);
    unitvect[1] = cos(dec/rad2deg) * sin(ra/rad2deg);
    unitvect[2] = sin(dec/rad2deg);
 
    dphiptot = i*dphip;
    phiprad = (phipinit + dphiptot)/rad2deg;
    vearth[0] = -sin(phiprad)*vorbrel;
    vearth[1] = cos(phiprad)*coseps*vorbrel;
    vearth[2] = cos(phiprad)*sineps*vorbrel;
        
    dphirottot = i*dphirot;
    phirotrad = (phirotinit + dphirottot)/rad2deg;
    vdet[0] = -sin(phirotrad)*vrotrel;
    vdet[1] = cos(phirotrad)*vrotrel;
    vdet[2] = 0.;
    
    for (j = 0; j < 3; ++j)
    {
        vtotal[j] = vearth[j] + vdet[j];
        retval += vtotal[j] * unitvect[j];
    }
    
    return retval;
}

void compute_skygrid(LALStatus * status, char * format_arg)
{
  LALDetector             detector;
  LALFrDetector           fr_detector;
  LALSource               source;
  LALDetAndSource         det_and_source = {NULL, NULL};
  LALDetAMResponse        response;
  LALGPSandAcc            gps_and_acc;
  LALTimeIntervalAndNSample time_info;
  LALTimeInterval         time_interval;
  skygrid_t               grid_cros_sq;  
  skygrid_t               grid_plus_sq;
  skygrid_t               grid_sum_sq;
  skygrid_t               grid_relfreq;
  CHAR                    cross_file_name[LALNameLength];
  CHAR                    plus_file_name[LALNameLength];
  CHAR                    sum_file_name[LALNameLength];
  CHAR                    relfreq_file_name[LALNameLength];
  CHAR                    ser_cross_file_name[LALNameLength];
  CHAR                    ser_plus_file_name[LALNameLength];
  CHAR                    ser_sum_file_name[LALNameLength];
  CHAR                    ser_relfreq_file_name[LALNameLength];
  CHAR                    dottimestamp[LALNameLength];
  CHAR                    outfile_suffix[5];
  INT4                    i, j, k, cnt;
  FILE                  * timesfile;

  /* null out strings */
  cross_file_name[0] = '\0';
  plus_file_name[0] = '\0';
  sum_file_name[0] = '\0';
  relfreq_file_name[0] = '\0';

  if (args_info.detector_given)
    {
      if (mystrncasecmp(args_info.detector_arg,"lho",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
      else if (mystrncasecmp(args_info.detector_arg,"llo",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
      else if (mystrncasecmp(args_info.detector_arg,"virgo",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      else if (mystrncasecmp(args_info.detector_arg,"geo",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      else if (mystrncasecmp(args_info.detector_arg,"tama",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      else if (mystrncasecmp(args_info.detector_arg,"cit",LALNameLength) == 0)
        detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
      else if (mystrncasecmp(args_info.detector_arg,"test",LALNameLength) == 0)
        {
          set_detector_params(status, &fr_detector, &detector,
                              "TEST - North Pole",
                              0., LAL_PI_2, 0.,
                              0., LAL_PI_2,
                              0., 0.);
        }
      else
        {
          /* unknown detector -- exit with error */
          fprintf(stderr, "lalapps_detresponse: Unknown detector '%s'\n",
                  args_info.detector_arg);
          exit(2);
        }

    }
  
  if (verbosity_level & 4)
    PrintLALDetector(&detector);

  if (!(args_info.orientation_given))
    {
      fprintf(stderr, "lalapps_detresponse: ERROR: must specify orientation of source for whole sky mode\n");
      exit(4);
    }
  else
    {
      set_source_params(&source, args_info.source_name_arg,
                        0., 0., args_info.orientation_arg);
    }

  if (verbosity_level & 4)
    print_source(&source);

  /*
   * set up time info
   */
  if (args_info.start_time_sec_given &&
      args_info.nsample_given &&
      args_info.sampling_interval_given)
    {
      /* for time-averaged */
      time_info.accuracy             = LALLEAPSEC_STRICT;
      time_info.epoch.gpsSeconds     = args_info.start_time_sec_arg;
      time_info.epoch.gpsNanoSeconds = args_info.start_time_nanosec_arg;
      time_info.deltaT               = args_info.sampling_interval_arg;
      time_info.nSample              = args_info.nsample_arg;

      LALFloatToInterval(status, &time_interval, &(time_info.deltaT));

      /* for snapshot */
      gps_and_acc.accuracy = time_info.accuracy;
      gps_and_acc.gps.gpsSeconds = time_info.epoch.gpsSeconds;
      gps_and_acc.gps.gpsNanoSeconds = time_info.epoch.gpsNanoSeconds;
    }
  else
    {
      fprintf(stderr, "ERROR: time info incomplete\n");
      exit(5);
    }


  /*
   * set up output file names
   */
      
  /* first, the suffixes */
  if ((strncasecmp(args_info.format_arg, "am", LALNameLength) == 0) ||
      (strncasecmp(args_info.format_arg, "af", LALNameLength) == 0) ||
      (strncasecmp(args_info.format_arg, "al", LALNameLength) == 0))
    {
      (void)mystrlcpy(outfile_suffix, ".txt", LALNameLength);
    }
  else
    {
      (void)mystrlcpy(outfile_suffix, ".dat", LALNameLength);
    }
   
  (void)mystrlcpy(cross_file_name, "whole_sky_cross", LALNameLength);
  (void)mystrlcpy(plus_file_name, "whole_sky_plus", LALNameLength);
  (void)mystrlcpy(sum_file_name, "whole_sky_sum", LALNameLength);

  (void)strncat(cross_file_name, outfile_suffix, LALNameLength);
  (void)strncat(plus_file_name, outfile_suffix, LALNameLength);
  (void)strncat(sum_file_name, outfile_suffix, LALNameLength);  

  if (lalDebugLevel && (verbosity_level & 4))
    {
      printf("cross_file_name = %s\n", cross_file_name);
      printf("plus_file_name  = %s\n", plus_file_name);
      printf("sum_file_name   = %s\n", sum_file_name);
    }

  /*
   * compute response over whole sky
   */

  /* SNAPSHOT */
  det_and_source.pDetector = &detector;
  det_and_source.pSource   = &source;

  for (j = 0; j < NUM_RA; ++j)
    {
      source.equatorialCoords.longitude =
        (REAL8)j/(REAL8)NUM_RA * ((REAL8)LAL_TWOPI);
      
      for (i = -dec_lim; i <= dec_lim; ++i)
        {
          cnt = j*NUM_DEC + i + dec_lim;
          source.equatorialCoords.latitude = asin((REAL8)i/(REAL8)dec_lim);

          LALComputeDetAMResponse(status, &response, &det_and_source,
                                  &gps_and_acc);

          grid_cros_sq[cnt] = response.cross * response.cross;
          grid_plus_sq[cnt] = response.plus  * response.plus;
          grid_sum_sq[cnt]  = grid_cros_sq[cnt] + grid_plus_sq[cnt];
        }
    }

  skygrid_print(format_arg, &(gps_and_acc.gps), grid_cros_sq, cross_file_name, "w");
  skygrid_print(format_arg, &(gps_and_acc.gps), grid_plus_sq, plus_file_name, "w");
  skygrid_print(format_arg, &(gps_and_acc.gps), grid_sum_sq, sum_file_name, "w");


  /* TIME SERIES */
  (void)mystrlcpy(cross_file_name, "ser_whole_sky_cross", LALNameLength);
  (void)mystrlcpy(plus_file_name, "ser_whole_sky_plus", LALNameLength);
  (void)mystrlcpy(sum_file_name, "ser_whole_sky_sum", LALNameLength);
  (void)mystrlcpy(relfreq_file_name, "ser_whole_sky_relfreq", LALNameLength);
  
  /* zero out arrays */
  skygrid_zero(grid_cros_sq);
  skygrid_zero(grid_plus_sq);
  skygrid_zero(grid_sum_sq);
  skygrid_zero(grid_relfreq);
  
  timesfile = xfopen("times.txt", "w");
  
  for (k = 0; k < (int)time_info.nSample; ++k)
    {
      fprintf(timesfile, "%.9d\n", gps_and_acc.gps.gpsSeconds);
      snprintf(dottimestamp, LALNameLength, ".%.9d", gps_and_acc.gps.gpsSeconds);
      
      for (j = 0; j < NUM_RA; ++j)
        {
          source.equatorialCoords.longitude =
            (REAL8)j/(REAL8)NUM_RA * ((REAL8)LAL_TWOPI);
      
          for (i = -dec_lim; i <= dec_lim; ++i)
            {
              cnt = j*NUM_DEC + i + dec_lim;
              source.equatorialCoords.latitude = asin((REAL8)i/(REAL8)dec_lim);

              LALComputeDetAMResponse(status, &response, &det_and_source,
                                      &gps_and_acc);

              grid_cros_sq[cnt] = response.cross * response.cross;
              grid_plus_sq[cnt] = response.plus  * response.plus;
              grid_sum_sq[cnt]  = (grid_cros_sq[cnt] + grid_plus_sq[cnt]);
                
              grid_relfreq[cnt] += relval(source.equatorialCoords.latitude,
                                          source.equatorialCoords.longitude,
                                          k, (int)time_info.nSample);
            }
        }
        
        LALIncrementGPS(status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                      &time_interval);
       
        /* set up filenames */
        (void)mystrlcpy(ser_cross_file_name, cross_file_name, LALNameLength);
        (void)mystrlcpy(ser_plus_file_name, plus_file_name, LALNameLength);
        (void)mystrlcpy(ser_sum_file_name, sum_file_name, LALNameLength);
        (void)mystrlcpy(ser_relfreq_file_name, relfreq_file_name, LALNameLength);

        /* al files have no format suffix */
        if (!(strncasecmp(args_info.format_arg, "al", LALNameLength) == 0))
	        {
            (void)strncat(ser_cross_file_name, dottimestamp, LALNameLength);
            (void)strncat(ser_plus_file_name, dottimestamp, LALNameLength);
            (void)strncat(ser_sum_file_name, dottimestamp, LALNameLength);  
            (void)strncat(ser_relfreq_file_name, dottimestamp, LALNameLength);  
	        }

        (void)strncat(ser_cross_file_name, outfile_suffix, LALNameLength);
        (void)strncat(ser_plus_file_name, outfile_suffix, LALNameLength);
        (void)strncat(ser_sum_file_name, outfile_suffix, LALNameLength);  
        (void)strncat(ser_relfreq_file_name, outfile_suffix, LALNameLength);  
        
        skygrid_print(format_arg, &(gps_and_acc.gps), grid_cros_sq, ser_cross_file_name, "w");
        skygrid_print(format_arg, &(gps_and_acc.gps), grid_plus_sq, ser_plus_file_name, "w");
        skygrid_print(format_arg, &(gps_and_acc.gps), grid_sum_sq, ser_sum_file_name, "w");
        skygrid_print(format_arg, &(gps_and_acc.gps), grid_relfreq, ser_relfreq_file_name, "w");
        
    } 
  fclose(timesfile);    
  
  /* TIME-AVERAGED */
  (void)mystrlcpy(cross_file_name, "int_whole_sky_cross", LALNameLength);
  (void)mystrlcpy(plus_file_name, "int_whole_sky_plus", LALNameLength);
  (void)mystrlcpy(sum_file_name, "int_whole_sky_sum", LALNameLength);
  
  (void)strncat(cross_file_name, outfile_suffix, LALNameLength);
  (void)strncat(plus_file_name, outfile_suffix, LALNameLength);
  (void)strncat(sum_file_name, outfile_suffix, LALNameLength);  

  /* zero out arrays */
  skygrid_zero(grid_cros_sq);
  skygrid_zero(grid_plus_sq);
  skygrid_zero(grid_sum_sq);

  for (k = 0; k < (int)time_info.nSample; ++k)
    {
      for (j = 0; j < NUM_RA; ++j)
        {
          source.equatorialCoords.longitude =
            (REAL8)j/(REAL8)NUM_RA * ((REAL8)LAL_TWOPI);
      
          for (i = -dec_lim; i <= dec_lim; ++i)
            {
              cnt = j*NUM_DEC + i + dec_lim;
              source.equatorialCoords.latitude = asin((REAL8)i/(REAL8)dec_lim);

              LALComputeDetAMResponse(status, &response, &det_and_source,
                                      &gps_and_acc);

              grid_cros_sq[cnt] += response.cross * response.cross /
                time_info.nSample;
              grid_plus_sq[cnt] += response.plus  * response.plus /
                time_info.nSample;
              grid_sum_sq[cnt]  += (grid_cros_sq[cnt] + grid_plus_sq[cnt]) /
                time_info.nSample;
            }
        }

      LALIncrementGPS(status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                      &time_interval);
    }

  skygrid_print(format_arg, &(gps_and_acc.gps), grid_cros_sq, cross_file_name, "w");
  skygrid_print(format_arg, &(gps_and_acc.gps), grid_plus_sq, plus_file_name, "w");
  skygrid_print(format_arg, &(gps_and_acc.gps), grid_sum_sq, sum_file_name, "w");
  
  LALCheckMemoryLeaks();
  
  return;
}

REAL4 skygrid_avg(const skygrid_t response)
{
  INT4 i, j, cnt;
  REAL4 retval = 0.;

  for (j = 0; j < NUM_RA; ++j)
    {
      for (i = -dec_lim+1; i <= dec_lim-1; ++i)
        {
          cnt = j*NUM_DEC + i + dec_lim;
          retval += response[cnt];
        }
    }

  /*
  for (i = 0; i < lim; ++i)
    retval += response[i];
  */

  retval /= grid_lim;

  return retval;
}




void skygrid_square(skygrid_t square, const skygrid_t input)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    square[i] = (input[i]) * (input[i]);
  
}




void skygrid_sqrt(skygrid_t result, const skygrid_t input)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    result[i] = (REAL4)sqrt((double)(input[i]));
}



REAL4 skygrid_rms(const skygrid_t input)
{
  skygrid_t tmpgrid;

  skygrid_square(tmpgrid, input);
  return (REAL4)(sqrt(skygrid_avg(tmpgrid)));
}



INT4 skygrid_copy(skygrid_t dest, const skygrid_t src)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    dest[i] = src[i];

  return i;
}



void skygrid_print(const char * format,
                   const LIGOTimeGPS * gps,
                   const skygrid_t input, 
                   const char * filename,
                   const char * mode)
{
  INT4 i, j, fmtcode;
  FILE * outfile = NULL;
  
  if (strncasecmp(format, "am", LALNameLength) == 0)
    fmtcode = 0;
  else if (strncasecmp(format, "af", LALNameLength) == 0)
    fmtcode = 1;
  else if (strncasecmp(format, "b", LALNameLength) == 0)
    fmtcode = 3;
  else if (strncasecmp(format, "al", LALNameLength) == 0)
    fmtcode = 4;
  else 
    fmtcode = 0;

  /* al files are in one big contiguous chunk, rather
   * than separate files */
  if (fmtcode == 4)
    mode = "a";

  outfile = xfopen(filename, mode);
  
  if (fmtcode == 0) /* am */
    {
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%c %.9d\n", '%', gps->gpsSeconds);

      for (i = 0; i < NUM_RA; ++i)
        {
          for (j = 0; j < NUM_DEC; ++j)
            fprintf(outfile, "% 14.8e\t", input[i*NUM_DEC + j]);
          fprintf(outfile, "\n");
        }
    }
  else if (fmtcode == 1) /* af */
    {
      const int tot = NUM_RA * NUM_DEC;
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%c %.9d\n", '%', gps->gpsSeconds);
        
      for (i = 0; i < tot; ++i)
        fprintf(outfile, "% 14.8e\t", input[i]); 
      fprintf(outfile, "\n");
    }
  else if (fmtcode == 3) /* b */
    {
      /* no comments for binary output */
      const size_t num_el = NUM_DEC*NUM_RA;
      if (gps == (LIGOTimeGPS *)NULL)
        {
          fprintf(stderr, "must specify GPS timestamp for binary format\n");
          exit(22);
        }
      fwrite(&(gps->gpsSeconds), sizeof(INT4), 1, outfile);
      fwrite(&num_el, sizeof(INT4), 1, outfile);
      fwrite(input, sizeof(REAL4), num_el, outfile);
    }
  else if (fmtcode == 4) /* al */
    {
      const int tot = NUM_RA * NUM_DEC;
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%.9d", gps->gpsSeconds);
        
        for (i = 0; i < tot; ++i)
            fprintf(outfile, " %14.8e", input[i]); 
        fprintf(outfile, "\n");
    }

  xfclose(outfile);
}



void skygrid_fabs(skygrid_t absgrid, const skygrid_t input)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    absgrid[i] = fabs(input[i]);
}



void skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    sum[i] = a[i] + b[i];
}

void skygrid_subtract(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    sum[i] = a[i] - b[i];
}

void skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    result[i] = b * a[i];
}

void skygrid_zero(skygrid_t a)
{
  INT4 i;

  for (i = 0; i < grid_lim; ++i)
    a[i] = 0.;
}
