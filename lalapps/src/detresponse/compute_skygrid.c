/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include <stdio.h>
#include <ctype.h>
#include "config.h"
#include "detresponse.h"
#include "skygrid.h"


/*
 * TODO: 
 *   1. make file names depend on name of source, and of detector
 *       
 */

static const double rad2deg = 180./(double)LAL_PI;
static const double     eps = 23.5;
static const double  clight = 2.998e8;

static double vorbrel;

/*
 * compute instantaneous doppler factor
 * NOTE: for some unfathomable reason, LALDetectorVel() takes
 *   AvgVelPar as input, while LALAvgDetectorVel() takes
 *   VelocityPar as input. -> the structure AvgVelPar is now killed - Badri Krishnan
 */
/* gps = time when vel. is needed */   
/* inputs = detector, ephemerides */
static double doppler(LALStatus * status, 
                      REAL8         velocity[3],  
                      SkyPosition * source_loc)
{
  /* the reference frame in which everything is done is the barycentric,
  * with J2000 equatorial defining the x-y plane */
  double          e_source[3];
  double          doppler_factor = 0.;
  double          sin_theta;
  int             i;
  
  /* source_loc needs to be in Equatorial */
  if (source_loc->system != COORDINATESYSTEM_EQUATORIAL)
  {
    /* puke and go */
    exit(81);
  }
  
  /* NOTE: from the LALBarycenter() file, here's
   * how to compute the unit vector pointing to a 
   * source of RA=alpha, Dec=delta [the vector is
   * in the J2000 frame]
   *
   *     sinTheta=sin(LAL_PI/2.0-delta);
   *     // s is vector that points towards source 
   *     s[2]=cos(LAL_PI/2.0-delta); 
   *     // in Cartesian coords based on J2000
   *     s[1]=sinTheta*sin(alpha);
   *     // 0=x,1=y,2=z 
   *     s[0]=sinTheta*cos(alpha);   */
  
  sin_theta = sin(LAL_PI_2 - source_loc->latitude);   
  e_source[2] = cos(LAL_PI_2 - source_loc->latitude);
  e_source[1] = sin_theta * sin(source_loc->longitude);
  e_source[0] = sin_theta * cos(source_loc->longitude);
  
 for (i = 0; i < 3; ++i)
    doppler_factor += velocity[i] * e_source[i];

  /* printf("doppler_factor = %20.14e\n", doppler_factor); */
    
  return doppler_factor;
}


#if 0
static double approx_doppler(LIGOTimeGPS * gps,    
                             double alpha, double dalpha,
                             double delta, double ddelta)
{
  /* there are 31,556,925.9747 ephemeris seconds per year */
  const unsigned int secs_per_year = 31556926;
  /* autumnal equinox 2003 was at GPS time 748349233 */ 
  const unsigned int aut_equinox_2003 = 748349233;
  
  /* phi = angle along Earth's orbit */
  double theta; /* theta := pi/2 - delta, delta = Dec. of source */
  /* alpha = RA of source */
  
}
#endif

/*
 * relval() is deprecated
 */
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


/*
 * NOTES on gridding the sky:
 *   - want to not close the circle in RA
 *   - want to not reach the poles in Dec
 *
 */
void compute_skygrid(LALStatus * status, EphemerisData *p_ephemeris_data,
                     skygrid_t grid_cros_sq, skygrid_t grid_plus_sq,
                     skygrid_t grid_sum_sq, skygrid_t grid_relfreq)
{
  LALDetector             detector;
  LALFrDetector           fr_detector;
  LALSource               source;
  LALDetAndSource         det_and_source = {NULL, NULL};
  LALDetAMResponse        response;
  LALGPSandAcc            gps_and_acc;
  LIGOTimeGPS             start_time;
  LALTimeInterval         time_interval;
  LALLeapSecFormatAndAcc  leapsec_info = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4                    tmp_leapsecs;
  REAL8                   det_velocity[3];


  CHAR                    cross_file_name[LALNameLength];
  CHAR                    plus_file_name[LALNameLength];
  CHAR                    sum_file_name[LALNameLength];
  CHAR                    relfreq_file_name[LALNameLength];
  CHAR                    dottimestamp[LALNameLength];
  CHAR                    outfile_suffix[5];
  UINT4                   i, j, k, cnt, cnt_offset;
  UINT4                   num_ra, count_ra, start_ra, end_ra;
  UINT4                   num_dec, count_dec, start_dec, end_dec;
  REAL8                   Pi_num_ra;
  FILE                  * timesfile = NULL;
  
  /* null out strings */
  cross_file_name[0] = '\0';
  plus_file_name[0] = '\0';
  sum_file_name[0] = '\0';
  relfreq_file_name[0] = '\0';
  outfile_suffix[0] = '\0';
  
  /* useful numbers */
  num_ra = args_info.n_ra_arg;
  count_ra = args_info.count_ra_arg;
  start_ra = args_info.start_ra_arg;
  end_ra = args_info.count_ra_arg + start_ra;
  
  num_dec = args_info.n_dec_arg;
  count_dec = args_info.count_dec_arg;
  start_dec = args_info.start_dec_arg;
  end_dec = args_info.count_dec_arg + start_dec;
  
  cnt_offset = start_ra*count_dec + start_dec;
    
  Pi_num_ra = (REAL8)LAL_PI / (REAL8)num_ra;
  
  /* 
   * set up detector 
   */
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
  
  /* 
   * set up source 
   */
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
   * set up time info -- snapshot, or timeseries
   */
  if (args_info.timeseries_given)
  { 
    if (args_info.start_time_sec_given &&
        args_info.nsample_given &&
        args_info.sampling_interval_given)
    {
      gps_and_acc.accuracy = leapsec_info.accuracy;
      gps_and_acc.gps.gpsSeconds = args_info.start_time_sec_arg;
      gps_and_acc.gps.gpsNanoSeconds = args_info.start_time_nanosec_arg;
    
      start_time.gpsSeconds = gps_and_acc.gps.gpsSeconds;
      start_time.gpsNanoSeconds = gps_and_acc.gps.gpsNanoSeconds;
            
      LALFloatToInterval(status, &time_interval, 
                           &(args_info.sampling_interval_arg));    
    }
    else
    {
      fprintf(stderr, "ERROR: time info incomplete\n");
      exit(5);
    }
  } 
  else if (args_info.snapshot_given)
  {
    if (args_info.start_time_sec_given)
    {
      gps_and_acc.accuracy = leapsec_info.accuracy;
      gps_and_acc.gps.gpsSeconds = args_info.start_time_sec_arg;
      gps_and_acc.gps.gpsNanoSeconds = args_info.start_time_nanosec_arg;
    }
    else
    {
      fprintf(stderr, "ERROR: time info incomplete\n");
      exit(6);
    }
  }
  

  /*
   * set up output file suffixes
   */
  /* 'al' format has no suffix */
  if (tolower(args_info.format_arg[1]) == 'a')
  {
    if (tolower(args_info.format_arg[2]) != 'l')
      (void *)mystrlcpy(outfile_suffix, ".txt", LALNameLength);
  }
  else
  {
    (void *)mystrlcpy(outfile_suffix, ".dat", LALNameLength);
  }
  
  
  /*
   * compute response over whole sky
   */
  det_and_source.pDetector = &detector;
  det_and_source.pSource   = &source;
  
  if (args_info.snapshot_given)
  {
    /* 
     * SNAPSHOT 
     */
    snprintf(cross_file_name,LALNameLength,"%s/whole_sky_cross", args_info.output_dir_arg);
    snprintf(plus_file_name,LALNameLength,"%s/whole_sky_plus", args_info.output_dir_arg);
    snprintf(sum_file_name,LALNameLength,"%s/whole_sky_sum", args_info.output_dir_arg);
    snprintf(relfreq_file_name,LALNameLength,"%s/whole_sky_relfreq", args_info.output_dir_arg);
    
    (void)mystrlcat(cross_file_name, outfile_suffix, LALNameLength);
    (void)mystrlcat(plus_file_name, outfile_suffix, LALNameLength);
    (void)mystrlcat(sum_file_name, outfile_suffix, LALNameLength);  
    (void)mystrlcat(relfreq_file_name, outfile_suffix, LALNameLength);  
    
    if (lalDebugLevel && (verbosity_level & 4))
    {
      printf("cross_file_name   = %s\n", cross_file_name);
      printf("plus_file_name    = %s\n", plus_file_name);
      printf("sum_file_name     = %s\n", sum_file_name);
      printf("relfreq_file_name = %s\n", relfreq_file_name);
    }    
    
    /*
     * set leapsec info in ephemerides
     */
    LALLeapSecs(status, &tmp_leapsecs, &(gps_and_acc.gps), &leapsec_info);
    p_ephemeris_data->leap = (INT2)tmp_leapsecs; 
    

    LALDetectorVel(status, det_velocity, &(gps_and_acc.gps), detector, p_ephemeris_data);
    
    Pi_num_ra = (REAL8)LAL_PI/(REAL8)num_ra;
    for (j = start_ra; j < end_ra; ++j)
    {
      source.equatorialCoords.longitude = Pi_num_ra * (1. + 2.*(REAL8)j);

      if (lalDebugLevel)
      {
        printf("j = %d\t\t", j);
        printf("RA = % f\n", source.equatorialCoords.longitude);
      }
      
      for (i = start_dec; i < end_dec; ++i)
      {
        cnt = j*count_dec + i - cnt_offset;
        
        if (lalDebugLevel)
        {
          printf("\ti = %d\n", i);
          printf("\t\tcnt = %d\n", cnt);
        }
        
        source.equatorialCoords.latitude = 
          asin(-1. + (1. + 2.*i)/(REAL8)num_dec);

        /*
        printf("\ti = %d\n", i);
        printf("\tcnt = %d\n", cnt);
        printf("\tDec = % f\n", source.equatorialCoords.latitude / LAL_PI * 180.);
        */

        LALComputeDetAMResponse(status, &response, &det_and_source,
                                &gps_and_acc);
        
        grid_cros_sq->data[cnt] = response.cross * response.cross;
        grid_plus_sq->data[cnt] = response.plus  * response.plus;
        grid_sum_sq->data[cnt]  = grid_cros_sq->data[cnt] + grid_plus_sq->data[cnt];
        grid_relfreq->data[cnt] = doppler(status, det_velocity,
                                    &(source.equatorialCoords));
      }
    }

    skygrid_print(status, &(gps_and_acc.gps), grid_cros_sq, cross_file_name);
    skygrid_print(status, &(gps_and_acc.gps), grid_plus_sq, plus_file_name);
    skygrid_print(status, &(gps_and_acc.gps), grid_sum_sq, sum_file_name);
    skygrid_print(status, &(gps_and_acc.gps), grid_relfreq, relfreq_file_name);
  }
  else if (args_info.timeseries_given)
  {
    /* 
     * TIME SERIES 
     */
    CHAR                    ser_cross_file_name[LALNameLength];
    CHAR                    ser_plus_file_name[LALNameLength];
    CHAR                    ser_sum_file_name[LALNameLength];
    CHAR                    ser_relfreq_file_name[LALNameLength];
    CHAR		    times_file_name[LALNameLength];

    snprintf(cross_file_name,LALNameLength,"%s/ser_whole_sky_cross", args_info.output_dir_arg);
    snprintf(plus_file_name,LALNameLength,"%s/ser_whole_sky_plus", args_info.output_dir_arg);
    snprintf(sum_file_name,LALNameLength,"%s/ser_whole_sky_sum", args_info.output_dir_arg);
    snprintf(relfreq_file_name,LALNameLength,"%s/ser_whole_sky_relfreq", args_info.output_dir_arg);
    
    /* zero out arrays */
    skygrid_zero(status, grid_cros_sq);
    skygrid_zero(status, grid_plus_sq);
    skygrid_zero(status, grid_sum_sq);
    skygrid_zero(status, grid_relfreq);
    
    snprintf(times_file_name,LALNameLength,"%s/times.txt", args_info.output_dir_arg);
    timesfile = xfopen(times_file_name, "w");
    
    for (k = 0; k < (UINT4)args_info.nsample_arg; ++k)
    {
      fprintf(timesfile, "%.9d\n", gps_and_acc.gps.gpsSeconds);
      snprintf(dottimestamp, LALNameLength, ".%.9d", gps_and_acc.gps.gpsSeconds);
      
      /*
       * set leapsec info in ephemerides
       */
      LALLeapSecs(status, &tmp_leapsecs, &(gps_and_acc.gps), &leapsec_info);
      p_ephemeris_data->leap = (INT2)tmp_leapsecs; 
      


      LALDetectorVel(status, det_velocity, &(gps_and_acc.gps), detector, p_ephemeris_data);
      
      for (j = start_ra; j < end_ra; ++j)
      {
        source.equatorialCoords.longitude = Pi_num_ra * (1. + 2.*j);
  
        for (i = start_dec; i < end_dec; ++i)
        {
          cnt = j*count_dec + i - cnt_offset;
          source.equatorialCoords.latitude = 
            asin(-1. + (1. + 2.*i)/(REAL8)num_dec);
          
          LALComputeDetAMResponse(status, &response, &det_and_source,
                                  &gps_and_acc);
          
          grid_cros_sq->data[cnt] = response.cross * response.cross;
          grid_plus_sq->data[cnt] = response.plus  * response.plus;
          grid_sum_sq->data[cnt]  = (grid_cros_sq->data[cnt] + grid_plus_sq->data[cnt]);
          
          grid_relfreq->data[cnt] = doppler(status, det_velocity, 
                                      &(source.equatorialCoords));
        }
      }
      
      /* DEBUG */
      /*
       source.equatorialCoords.latitude = 0.1;
       source.equatorialCoords.longitude = 0.1;
       printf("%.8e\n", doppler(status, &(gps_and_acc.gps), 
                                &detectorvel_inputs, 
                                &(source.equatorialCoords)));
       */
      
      LALIncrementGPS(status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                      &time_interval);
      
      /* set up filenames */
      (void)mystrlcpy(ser_cross_file_name, cross_file_name, LALNameLength);
      (void)mystrlcpy(ser_plus_file_name, plus_file_name, LALNameLength);
      (void)mystrlcpy(ser_sum_file_name, sum_file_name, LALNameLength);
      (void)mystrlcpy(ser_relfreq_file_name, relfreq_file_name, LALNameLength);
      
      if ((args_info.format_arg[0]=='m') ||
          (args_info.format_arg[0]=='M'))
      {
        (void)mystrlcat(ser_cross_file_name, dottimestamp, LALNameLength);
        (void)mystrlcat(ser_plus_file_name, dottimestamp, LALNameLength);
        (void)mystrlcat(ser_sum_file_name, dottimestamp, LALNameLength);  
        (void)mystrlcat(ser_relfreq_file_name, dottimestamp, LALNameLength);  

        (void)mystrlcat(ser_cross_file_name, outfile_suffix, LALNameLength);
        (void)mystrlcat(ser_plus_file_name, outfile_suffix, LALNameLength);
        (void)mystrlcat(ser_sum_file_name, outfile_suffix, LALNameLength);  
        (void)mystrlcat(ser_relfreq_file_name, outfile_suffix, LALNameLength);  
      }
      
      skygrid_print(status, &(gps_and_acc.gps), grid_cros_sq, ser_cross_file_name);
      skygrid_print(status, &(gps_and_acc.gps), grid_plus_sq, ser_plus_file_name);
      skygrid_print(status, &(gps_and_acc.gps), grid_sum_sq, ser_sum_file_name);
      skygrid_print(status, &(gps_and_acc.gps), grid_relfreq, ser_relfreq_file_name);
    } 
    fclose(timesfile);    
  }
  else if (args_info.average_given)
  {
    /* TIME-AVERAGED */
    (void)mystrlcpy(cross_file_name, "avg_whole_sky_cross", LALNameLength);
    (void)mystrlcpy(plus_file_name, "avg_whole_sky_plus", LALNameLength);
    (void)mystrlcpy(sum_file_name, "avg_whole_sky_sum", LALNameLength);
    (void)mystrlcpy(relfreq_file_name, "avg_whole_sky_relfreq", LALNameLength);
    
    snprintf(dottimestamp, LALNameLength, ".%.9d", start_time.gpsSeconds);
    
    /* al files have no format suffix */
    if (!(strncasecmp(args_info.format_arg, "al", LALNameLength) == 0))
    {
      (void)mystrlcat(cross_file_name, dottimestamp, LALNameLength);
      (void)mystrlcat(plus_file_name, dottimestamp, LALNameLength);
      (void)mystrlcat(sum_file_name, dottimestamp, LALNameLength);  
      (void)mystrlcat(relfreq_file_name, dottimestamp, LALNameLength);  
    
      (void)mystrlcat(cross_file_name, outfile_suffix, LALNameLength);
      (void)mystrlcat(plus_file_name, outfile_suffix, LALNameLength);
      (void)mystrlcat(sum_file_name, outfile_suffix, LALNameLength); 
      (void)mystrlcat(relfreq_file_name, outfile_suffix, LALNameLength); 
    }
    
    /* zero out arrays */
    skygrid_zero(status, grid_cros_sq);
    skygrid_zero(status, grid_plus_sq);
    skygrid_zero(status, grid_sum_sq);
    skygrid_zero(status, grid_relfreq);
    
    for (k = 0; k < (UINT4)args_info.nsample_arg; ++k)
    {
      /*
       * set leapsec info in ephemerides
       */
      LALLeapSecs(status, &tmp_leapsecs, &(gps_and_acc.gps), &leapsec_info);
      p_ephemeris_data->leap = (INT2)tmp_leapsecs; 


      LALDetectorVel(status, det_velocity, &(gps_and_acc.gps), detector, p_ephemeris_data);
          
      for (j = start_ra; j < end_ra; ++j)
      {
        source.equatorialCoords.longitude = Pi_num_ra * (1. + 2.*j);

        for (i = start_dec; i < end_dec; ++i)
        {
          cnt = j*count_dec + i - cnt_offset;
          source.equatorialCoords.latitude = 
            asin(-1. + (1. + 2.*i)/(REAL8)num_dec);
          
          LALComputeDetAMResponse(status, &response, &det_and_source,
                                  &gps_and_acc);
          
          grid_cros_sq->data[cnt] += response.cross * response.cross /
            args_info.nsample_arg;
          grid_plus_sq->data[cnt] += response.plus  * response.plus /
            args_info.nsample_arg;
          grid_sum_sq->data[cnt]  += (grid_cros_sq->data[cnt] + grid_plus_sq->data[cnt]) /
            args_info.nsample_arg;
          grid_relfreq->data[cnt] += doppler(status, det_velocity, 
                                       &(source.equatorialCoords));
        }
      }
      
      LALIncrementGPS(status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                      &time_interval);
    }
    
    skygrid_print(status, &start_time, grid_cros_sq, cross_file_name);
    skygrid_print(status, &start_time, grid_plus_sq, plus_file_name);
    skygrid_print(status, &start_time, grid_sum_sq, sum_file_name);
    skygrid_print(status, &start_time, grid_relfreq, relfreq_file_name);
  }
  
  if (lalDebugLevel)
  {
    printf("grid_cros_sq->data = 0x%x\n", grid_cros_sq->data);
    printf("grid_plus_sq->data = 0x%x\n", grid_plus_sq->data);
    printf("grid_sum_sq->data = 0x%x\n", grid_sum_sq->data);
    printf("grid_relfreq->data = 0x%x\n", grid_relfreq->data);
  }
  
  return;
} /* END: compute_skygrid() */


/*
static double
ab1(double y_front, double y_back, double delta)
{
  return ((y_front - y_back) / (2. * delta));
}

static double
ab2(double y_center, double y_front, double y_back, double delta)
{
  return ((y_front + y_back - 2.*y_center)/(delta * delta));
}


approx_doppler_coeffs()
{
  double delta_ra = LAL_PI / (REAL8)num_ra;
  double delta_sin_dec = 1. / (REAL8)num_dec;
  double delta_dec = asin(delta_sin_dec);

  double a1, a2, b1, b2;

  for (j = start_ra; j < end_ra; ++j)
  {
    source.equatorialCoords.longitude = Pi_num_ra * (1. + 2.*j);

    for (i = start_dec; i < end_dec; ++i)
    {
      cnt = j*count_dec + i - cnt_offset;
      source.equatorialCoords.latitude = asin(-1. + (1. + 2.*i)/(REAL8)num_dec);
    }
  }
}
*/
