/* 
  Author: Pitkin, M. D.
  $Id$
*/

/* Matt Pitkin - 07/02/06 ------------- heterodyne_pulsar.c */

/* lalapps code to perform a coarse/fine heterodyne on LIGO or GEO data given a set of frequency
parameters. The heterodyne can either be a 2 stage process with the code run for a coarse
heterodyne (not taking into account the SSB and BSB time delays, but using the other frequency
params) and then rerun for the fine heterodyne (taking into account SSB and BSB time delays). The
code can also be used to update heterodyned data using new parameters (old and new parameter files
are required) i.e. it will take the difference of the original heterodyne phase and the new phase
and reheterodyne with this phase difference. This code does not perform any calibration of the data.
*/

#include "heterodyne_pulsar.h"

int main(int argc, char *argv[]){
  static LALStatus status;

  InputParams inputParams;
  HeterodyneParams hetParams;

  Filters iirFilters;

  FILE *fpin=NULL, *fpout=NULL;
  FrameCache cache;
  INT4 frcount=0, count=0;

  CHAR outputfile[256]="";
  CHAR channel[20]="";

  INT4Vector *starts=NULL, *stops=NULL; /* science segment start and stop times */
  INT4 numSegs=0;
  CHAR *smalllist=NULL; /* list of frame files for a science segment */

  REAL8Vector *times=NULL; /* times of data read from coarse heterodyne file */

  /* get input options */
  get_input_args(&inputParams, argc, argv);

  hetParams.heterodyneflag = inputParams.heterodyneflag; /* set type of heterodyne */

  /* read in pulsar data */
  LALReadTEMPOParFile(&status, &hetParams.het, inputParams.paramfile);

  if(inputParams.heterodyneflag == 1) /*if performing fine heterdoyne using same params as coarse*/
    hetParams.hetUpdate = hetParams.het;

  hetParams.samplerate = inputParams.samplerate;
  sprintf(hetParams.earthfile, "%s", inputParams.earthfile);
  sprintf(hetParams.sunfile,"%s" ,inputParams.sunfile);

  /* set detector */
  if (strcmp(inputParams.ifo,"H1") == 0 || strcmp(inputParams.ifo,"H2") == 0)
    hetParams.detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if (strcmp(inputParams.ifo,"L1") == 0)
    hetParams.detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (strcmp(inputParams.ifo,"GEO") == 0 || strcmp(inputParams.ifo,"G1") == 0)
    hetParams.detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else{
  fprintf(stderr,"Error DETECTOR must be either H1, H1, L1, or GEO\n");
    return 2;
  }

  if(inputParams.heterodyneflag == 2) /* if updating parameters read in updated par file */
    LALReadTEMPOParFile(&status, &hetParams.hetUpdate, inputParams.paramfileupdate);

  /* get science segment lists - allocate initial memory for starts and stops */
  starts = XLALCreateINT4Vector(MAXNUMFRAMES);
  stops = XLALCreateINT4Vector(MAXNUMFRAMES);
  numSegs = get_segment_list(starts, stops, inputParams.segfile);
  starts = XLALResizeINT4Vector(starts, numSegs);
  stops = XLALResizeINT4Vector(stops, numSegs);

  /* open input file */
  if((fpin = fopen(inputParams.datafile, "r")) == NULL){
    fprintf(stderr, "Error... Can't open input data file!\n");
    return 1;
  }

  if(inputParams.heterodyneflag == 0){ /* input comes from frame files so read in frame filenames */
    CHAR det[3]; /* detector from cache file */
    CHAR type[10]; /* frame type e.g. RDS_R_L3 - from cache file */
    
    frcount=0;
    while(fscanf(fpin, "%s%s%d%d file://localhost%s", det, type, &cache.starttime[frcount],
&cache.duration[frcount], cache.framelist[frcount]) != EOF){
      /* fscanf(fpin, "%s", framelist[frcount]);*/
      frcount++;
    }
    fclose(fpin);
  }

  /************************BIT THAT DOES EVERYTHING*********************************************/

  /* set filters - values held for the whole data set so we don't get lots of glitches from the
  filter ringing */
  if(inputParams.filterknee > 0.0)
    set_filters(&iirFilters, inputParams.filterknee, inputParams.samplerate);

  if(inputParams.heterodyneflag == 0){
    sprintf(outputfile, "%s/coarsehet_%s_%s_%d-%d", inputParams.outputdir, inputParams.pulsar,
    inputParams.ifo, starts->data[0], stops->data[numSegs-1]);
  }
  else{
    sprintf(outputfile, "%s/finehet_%s_%s", inputParams.outputdir, inputParams.pulsar,
    inputParams.ifo);
  }

  remove(outputfile); /* if output file already exists remove it */
  sprintf(channel, "%s:%s", inputParams.ifo, inputParams.channel);

  /* loop through file and read in data */
  /* if the file is a list of frame files read science segment at a time and perform the heterodyne
  on that set of data - this will only contain a real component - plus remember read in double
  precision for GEO and LIGO h(t) and, single precision for uncalibrated LIGO data. If the file is
  already heterodyned it will be complex and we can read the whole file in in one go as it should be
  significantly downsampled */
  do{
    COMPLEX16TimeSeries *data; /* data for heterodyning */
    COMPLEX16TimeSeries *resampData; /* resampled data */
    INT4 i;

    data = LALCalloc(1, sizeof(*data));

    if(inputParams.heterodyneflag == 0){/* i.e. reading from frame files */
      REAL8 gpstime;
      INT4 duration;
      REAL8TimeSeries *datareal;

      if((duration = stops->data[count] - starts->data[count]) > MAXDATALENGTH)
        duration = MAXDATALENGTH; /* if duration of science segment is large just get part of it */

      fprintf(stderr, "Getting data between %d and %d.\n", starts->data[count],
      starts->data[count]+duration);

      hetParams.timestamp = (REAL8)starts->data[count];
      hetParams.length = inputParams.samplerate * duration;
      gpstime = (REAL8)starts->data[count];

      smalllist = set_frame_files(&starts->data[count], &stops->data[count], cache, frcount,
      &count);

      /* make vector (make sure imaginary parts are set to zero) */
      data->data = NULL;
      data->data = XLALCreateCOMPLEX16Vector( inputParams.samplerate * duration );

      /* read in frame data */
      if((datareal = get_frame_data(smalllist, channel, gpstime,
        inputParams.samplerate * duration, duration)) == NULL){
        fprintf(stderr, "Error... could not open frame files between %d and %d.\n", (INT4)gpstime,
        (INT4)gpstime + duration);

        fprintf(stderr, "I've read in the data.\n");

        if(count < numSegs){
          count++; /* if not finished reading in all data try next set of frames */
          continue;
        }
        else
          break; /* if at the end of data anyway then break */
      }

      /* put data into COMPLEX16 vector and set imaginary parts to zero */
      for(i=0;i<inputParams.samplerate * duration;i++){
        data->data->data[i].re = (REAL8)datareal->data->data[i];
        data->data->data[i].im = 0.;
      }

      XLALDestroyREAL8Vector(datareal->data);
      LALFree(datareal);

      free(smalllist);
      count++;
    }
    else if(inputParams.heterodyneflag == 1 || inputParams.heterodyneflag == 2){/* i.e. reading
      from a heterodyned file */
      data->data = NULL;
      data->data = XLALCreateCOMPLEX16Vector( MAXLENGTH );

      times = XLALCreateREAL8Vector( MAXLENGTH );
      i=0;

      fprintf(stderr, "Reading heterodyned data from %s.\n", inputParams.datafile);

      while(fscanf(fpin, "%lf%lf%lf",&times->data[i],&data->data->data[i].re,
            &data->data->data[i].im) != EOF){
        if(i==0)
          hetParams.timestamp = times->data[i]; /* set initial time stamp */

        i++;
      };

      fclose(fpin);

      /* resize vector to actual size */
      data->data = XLALResizeCOMPLEX16Vector(data->data, i);
      hetParams.length = i;

      times = XLALResizeREAL8Vector(times, i);
    }
    else{
      fprintf(stderr, "Error... Heterodyne flag = %d, should be 0, 1, or 2.\n",
      inputParams.heterodyneflag);
      return 0;
    }

    data->deltaT = 1./inputParams.samplerate;
    data->epoch.gpsSeconds = (UINT4)floor(hetParams.timestamp);
    data->epoch.gpsNanoSeconds = (UINT4)(1.e9*(hetParams.timestamp - floor(hetParams.timestamp)));

    /* heterodyne data */
    heterodyne_data(data, times, hetParams);

    /* filter data */
    if(inputParams.filterknee > 0.) /* filter if knee frequency is not zero */
      filter_data(data, &iirFilters);

    if(inputParams.heterodyneflag==0){
      times = NULL;
      times = XLALCreateREAL8Vector( MAXLENGTH );
    }

    /* resample data and data times */
    resampData = resample_data(data, times, starts, stops, inputParams.samplerate,
inputParams.resamplerate, inputParams.heterodyneflag);

    XLALDestroyCOMPLEX16Vector( data->data );
    LALFree(data);

    /* calibrate */
    if(inputParams.calibrate)
      calibrate(resampData, times, inputParams.calibfiles, 2.0*hetParams.het.f0);

    /* remove outliers above our threshold */
    if(inputParams.stddevthresh != 0.)
      remove_outliers(resampData, times, inputParams.stddevthresh);

    /* output data */
    if((fpout = fopen(outputfile, "a"))==NULL){
      fprintf(stderr, "Error... can't open output file %s!\n", outputfile);
      return 0;
    }

    for(i=0;i<resampData->data->length;i++){
      fprintf(fpout, "%lf\t%le\t%le\n", times->data[i], resampData->data->data[i].re,
      resampData->data->data[i].im);
    }

    fclose(fpout);
    XLALDestroyCOMPLEX16Vector(resampData->data);
    LALFree(resampData);

    XLALDestroyREAL8Vector(times);
  }while(count < numSegs && inputParams.heterodyneflag==0);

  fprintf(stderr, "Heterodyning complete.\n");

  XLALDestroyINT4Vector(stops);
  XLALDestroyINT4Vector(starts);
  
  LALDestroyREAL8IIRFilter( &status, &iirFilters.filter1Re );
  LALDestroyREAL8IIRFilter( &status, &iirFilters.filter1Im ); 
  LALDestroyREAL8IIRFilter( &status, &iirFilters.filter2Re );
  LALDestroyREAL8IIRFilter( &status, &iirFilters.filter2Im );   
  LALDestroyREAL8IIRFilter( &status, &iirFilters.filter3Re );
  LALDestroyREAL8IIRFilter( &status, &iirFilters.filter3Im );
  
  return 0;
}

/* function to parse the input arguments */
void get_input_args(InputParams *inputParams, int argc, char *argv[]){
  struct option long_options[] =
  {
    { "help",                     no_argument,        0, 'h' },
    { "ifo",                      required_argument,  0, 'i' },
    { "pulsar",                   required_argument,  0, 'p' },
    { "heterodyne-flag",          required_argument,  0, 'z' },
    { "param-file",               required_argument,  0, 'f' },
    { "param-file-update",        required_argument,  0, 'g' },
    { "filter-knee",              required_argument,  0, 'k' },
    { "sample-rate",              required_argument,  0, 's' },
    { "resample-rate",            required_argument,  0, 'r' },
    { "data-file",                required_argument,  0, 'd' },
    { "channel",                  required_argument,  0, 'c' },
    { "output-dir",               required_argument,  0, 'o' },
    { "ephem-earth-file",         required_argument,  0, 'e' },
    { "ephem-sun-file",           required_argument,  0, 'S' },
    { "seg-file",                 required_argument,  0, 'l' },
    { "calibrate",                no_argument, &inputParams->calibrate, 1 },
    { "response-file",            required_argument,  0, 'R' },
    { "coefficient-file",         required_argument,  0, 'C' },
    { "sensing-function",         required_argument,  0, 'F' },
    { "open-loop-gain",           required_argument,  0, 'O' },
    { "stddev-thresh",            required_argument,  0, 'T' },
    { 0, 0, 0, 0 }
  };

  char args[] = "h:i:p:z:f:g:k:s:r:d:c:o:e:S:l:R:C:F:O:T";
  char *program = argv[0];

  /* set defaults */
  inputParams->filterknee = 1.0; /* filter knee frequency of 1 Hz */
  inputParams->resamplerate = 0.; /* resample to 1 Hz */
  inputParams->samplerate = 0.;
  inputParams->calibrate = 0;
  inputParams->stddevthresh = 0.;
  /*sprintf(inputParams->calibfiles.calibcoefficientfile, "%s", "");
  sprintf(inputParams->calibfiles.sensingfunctionfile, "%s", "");
  sprintf(inputParams->calibfiles.openloopgainfile, "%s", "");*/
  inputParams->calibfiles.calibcoefficientfile = NULL;
  inputParams->calibfiles.sensingfunctionfile = NULL;
  inputParams->calibfiles.openloopgainfile = NULL;

  /* get input arguments */
  while(1){
    int option_index = 0;
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch(c){
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "error parsing option %s with argument %s\n",
            long_options[option_index].name, optarg );
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'i': /* interferometer */
        sprintf(inputParams->ifo, "%s", optarg);
        break;
      case 'z': /* heterodyne flag - 0 for coarse, 1 for fine, 2 for update to params */
        inputParams->heterodyneflag = atoi(optarg);
        break;
      case 'p': /* pulsar name */
        sprintf(inputParams->pulsar, "%s", optarg);
        break;
      case 'f': /* initial heterodyne parameter file */
        sprintf(inputParams->paramfile, "%s", optarg);
        break;
      case 'g': /* secondary heterodyne parameter file - for updated parameters */
        sprintf(inputParams->paramfileupdate, "%s", optarg);
        break;
      case 'k': /* low-pass filter knee frequency */
        {/* find if the string contains a / and get its position */
          CHAR *loc=NULL;
          CHAR numerator[10]="", *denominator=NULL;
          INT4 n;

          if((loc = strchr(optarg, '/'))!=NULL){
            n = loc-optarg; /* length of numerator i.e. bit before / */
            strncpy(numerator, optarg, n);

            denominator = strdup(loc+1); /* set the denominator i.e. bit after / */

            inputParams->filterknee = atof(numerator)/atof(denominator);
          }
          else
            inputParams->filterknee = atof(optarg);
        }
        break;
      case 's': /* sample rate of input data */
        {/* find if the string contains a / and get its position */
          CHAR *loc=NULL;
          CHAR numerator[10]="", *denominator=NULL;
          INT4 n;

          if((loc = strchr(optarg, '/'))!=NULL){
            n = loc-optarg; /* length of numerator i.e. bit before / */
            strncpy(numerator, optarg, n);

            denominator = strdup(loc+1); /* set the denominator i.e. bit after / */

            inputParams->samplerate = atof(numerator)/atof(denominator);
          }
          else
            inputParams->samplerate = atof(optarg);
        }
        break;
      case 'r': /* resample rate - allow fractions e.g. 1/60 Hz*/
        {/* find if the string contains a / and get its position */
          CHAR *loc=NULL;
          CHAR numerator[10]="", *denominator=NULL;
          INT4 n;

          if((loc = strchr(optarg, '/'))!=NULL){
            n = loc-optarg; /* length of numerator i.e. bit before / */

            strncpy(numerator, optarg, n);

            denominator = strdup(loc+1); /* set the denominator i.e. bit after / */

            inputParams->resamplerate = atof(numerator)/atof(denominator);
          }
          else
            inputParams->resamplerate = atof(optarg);
        }
        break;
      case 'd': /* file containing list of frame files, or file with previously heterodyned data */
        sprintf(inputParams->datafile, "%s", optarg);
        break;
      case 'c': /* frame channel */
        sprintf(inputParams->channel, "%s", optarg);
        break;
      case 'o': /* output data directory */
        sprintf(inputParams->outputdir, "%s", optarg);
        break;
      case 'e': /* earth ephemeris file */
        sprintf(inputParams->earthfile, optarg);
        break;
      case 'S': /* sun ephemeris file */
        sprintf(inputParams->sunfile, optarg);
        break;
      case 'l':
        sprintf(inputParams->segfile, optarg);
        break;
      case 'R':
        inputParams->calibfiles.responsefunctionfile = optarg;
        break;
      case 'C':
        inputParams->calibfiles.calibcoefficientfile = optarg;
        break;
      case 'F':
        inputParams->calibfiles.sensingfunctionfile = optarg;
        break;
      case 'O':
        inputParams->calibfiles.openloopgainfile = optarg;
        break;
      case 'T':
        inputParams->stddevthresh = atof(optarg);
        break;
      case '?':
        fprintf(stderr, "unknown error while parsing options\n" );
      default:
        fprintf(stderr, "unknown error while parsing options\n" );
    }
  }

  /* set more defaults */
  /*if the resample rate hasn't been defined set equal to samplerate*/
  if(inputParams->resamplerate==0)
    inputParams->resamplerate = inputParams->samplerate;

  /* check that sample rate / resample rate is an integer */
  if(fmod(inputParams->samplerate/inputParams->resamplerate, 1.) != 0.){
    fprintf(stderr, "Error... invalid sample rates.\n");
    exit(0);
  }
}

/* heterodyne data function */
void heterodyne_data(COMPLEX16TimeSeries *data, REAL8Vector *times, HeterodyneParams hetParams){
  static LALStatus status;

  REAL8 phaseCoarse=0., phaseUpdate=0., deltaphase=0.;
  REAL8 t=0., tdt=0., T0=0., T0Update=0.;
  REAL8 dtpos=0.; /* time between position epoch and data timestamp */
  INT4 i=0;

  EphemerisData *edat=NULL;
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime  emit;

  BinaryPulsarInput binInput;
  BinaryPulsarOutput binOutput;
  
  COMPLEX16 dataTemp;
  
  /* set the position and frequency epochs if not already set */
  if(hetParams.het.pepoch == 0. && hetParams.het.posepoch != 0.)
    hetParams.het.pepoch = hetParams.het.posepoch;
  else if(hetParams.het.posepoch == 0. && hetParams.het.pepoch != 0.)
    hetParams.het.posepoch = hetParams.het.pepoch;

  if(hetParams.heterodyneflag > 0){
    if(hetParams.hetUpdate.pepoch == 0. && hetParams.hetUpdate.posepoch != 0.)
      hetParams.hetUpdate.pepoch = hetParams.hetUpdate.posepoch;
    else if(hetParams.het.posepoch == 0. && hetParams.hetUpdate.pepoch != 0.)
      hetParams.hetUpdate.posepoch = hetParams.hetUpdate.pepoch;
      
    T0Update = hetParams.hetUpdate.pepoch;
  }

  T0 = hetParams.het.pepoch;

  /* set up ephemeris files */
  if( hetParams.heterodyneflag > 0){
    edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));

    (*edat).ephiles.earthEphemeris = hetParams.earthfile;
    (*edat).ephiles.sunEphemeris = hetParams.sunfile;
    LALInitBarycenter(&status, edat);

    /* set up location of detector */
    baryinput.site.location[0] = hetParams.detector.location[0]/LAL_C_SI;
    baryinput.site.location[1] = hetParams.detector.location[1]/LAL_C_SI;
    baryinput.site.location[2] = hetParams.detector.location[2]/LAL_C_SI;
  }
  
  for(i=0;i<hetParams.length;i++){
    /*************************************************************************************/

    /* produce initial heterodyne phase for coarse heterodyne with no time delays */
    if(hetParams.heterodyneflag == 0)
      tdt = hetParams.timestamp + (REAL8)i/hetParams.samplerate - T0;
    else
      tdt = times->data[i] - T0;
    
    /* multiply by 2 to get gw phase */
    phaseCoarse = 2.*(hetParams.het.f0*tdt + 0.5*hetParams.het.f1*tdt*tdt +
    (1./6.)*hetParams.het.f2*tdt*tdt*tdt + (1./24.)*hetParams.het.f3*tdt*tdt*tdt*tdt);

    /**************************************************************************************/
    /* produce second phase for fine heterodyne */
    if(hetParams.heterodyneflag > 0){
      /* set up LALBarycenter */
      dtpos = hetParams.hetUpdate.posepoch - hetParams.timestamp;

      /* set up RA, DEC, and distance variables for LALBarycenter*/
      baryinput.alpha = hetParams.hetUpdate.ra + dtpos*hetParams.hetUpdate.pmra;
      baryinput.delta = hetParams.hetUpdate.dec + dtpos*hetParams.hetUpdate.pmdec;;
      baryinput.dInv = 0.0;  /* no parallax */

      t = times->data[i]; /* get data time */

      /* set leap seconds noting that for all runs prior to S5 that the number of leap seconds was
      13 and that 1 more leap seconds was added on 31st Dec 2005 24:00:00 i.e. GPS 820108813 */
      if(t <= 820108813)
        (*edat).leap = 13;
      else
        (*edat).leap = 14;

      baryinput.tgps.gpsSeconds = (UINT8)floor(t);
      baryinput.tgps.gpsNanoSeconds = (UINT8)floor((fmod(t,1.0)*1.e9));	 
      LALBarycenterEarth(&status, &earth, &baryinput.tgps, edat); 
      LALBarycenter(&status, &emit, &baryinput, &earth);

      /* if binary pulsar add extra time delay */
      if(hetParams.hetUpdate.model!=NULL){
        /* input SSB time into binary timing function */
        binInput.tb = t + emit.deltaT;

        /* calculate binary time delay */
        LALBinaryPulsarDeltaT(&status, &binOutput, &binInput, &hetParams.hetUpdate);

        /* add binary time delay */
        tdt = t - T0Update + emit.deltaT +  binOutput.deltaT;
      }
      else{
        tdt = t - T0Update + emit.deltaT;
      }
      /* multiply by 2 to get gw phase */
      phaseUpdate = 2.0*(hetParams.hetUpdate.f0*tdt + 0.5*hetParams.hetUpdate.f1*tdt*tdt +
      (1./6.)*hetParams.hetUpdate.f2*tdt*tdt*tdt + (1./24.)*hetParams.hetUpdate.f3*tdt*tdt*tdt*tdt);
    }

    /******************************************************************************************/
    dataTemp = data->data->data[i];

    /* perform heterodyne */
    if(hetParams.heterodyneflag == 0){
      deltaphase = 2.*LAL_PI*fmod(phaseCoarse, 1.);
      dataTemp.im = 0.; /* make sure imaginary part is zero */
    }
    if(hetParams.heterodyneflag > 0)
      deltaphase = 2.*LAL_PI*fmod(phaseUpdate - phaseCoarse, 1.);

    data->data->data[i].re = dataTemp.re*cos(-deltaphase) - dataTemp.im*sin(-deltaphase);
    data->data->data[i].im = dataTemp.re*sin(-deltaphase) + dataTemp.im*cos(-deltaphase);
  }

  if(hetParams.heterodyneflag > 0){
    LALFree(edat->ephemE);
    LALFree(edat->ephemS);
    LALFree(edat);
  }
}

/* function to extract the frame time and duration from the file name */
void get_frame_times(CHAR *framefile, REAL8 *gpstime, INT4 *duration){
  INT4 j=0;

  /* check framefile is not NULL */
  if(framefile == NULL){
    fprintf(stderr, "Error... Frame filename is not specified!\n");
    exit(1);
  }

  /* file names should end with GPS-length.gwf - *********-***.gwf, so we can
    locate the first '-' before the end of the file name */
  /* locate the start time and length of the data */
  for (j=strlen(framefile)-1; j>=0; j--)
    if (strstr((framefile+j), "-")!=NULL) break;

  /* get start times and durations (assumes the GPS time is 9 digits long) */
  *gpstime = atof(framefile+j-9);
  *duration=atoi(framefile+j+1);
}

/* function to read in frame data given a framefile and data channel */
REAL8TimeSeries *get_frame_data(CHAR *framefile, CHAR *channel, REAL8 time, REAL8 length, INT4
duration){
  REAL8TimeSeries *dblseries;

  FrFile *frfile=NULL;
  FrVect *frvect=NULL;

  LIGOTimeGPS epoch;
  INT4 i;

  dblseries = LALCalloc(1, sizeof(*dblseries));

  epoch.gpsSeconds = (UINT4)floor(time);
  epoch.gpsNanoSeconds = (UINT4)(1.e9*(time-floor(time)));
  dblseries->epoch = epoch;
  /* open frame file for reading */
  if((frfile = FrFileINew(framefile)) == NULL)
    return NULL; /* couldn't open frame file */

  strncpy( dblseries->name, channel, sizeof( dblseries->name ) - 1 );

  /* create data memory */
  dblseries->data = XLALCreateREAL8Vector( length );

  /* get data */
  /* read in frame data */
  if((frvect = FrFileIGetV(frfile, dblseries->name, time, (REAL8)duration)) == NULL)
    return NULL; /* couldn't read frame data */
  FrFileIEnd(frfile);

  /* fill into REAL8 vector */
  for(i=0;i<(INT4)length;i++){
    if(strstr(channel, "STRAIN") == NULL || strstr(channel, "DER_DATA") == NULL){ /* data is
uncalibrated single precision */
      dblseries->data->data[i] = (REAL8)frvect->dataF[i];
    }
    else /* data is calibrated h(t) */
      dblseries->data->data[i] = frvect->dataD[i];
  }

  FrVectFree(frvect);

  return dblseries;
}

/* function to set up three low-pass third order Butterworth IIR filters */
void set_filters(Filters *iirFilters, REAL8 filterKnee, REAL8 samplerate){
  static LALStatus status;
  COMPLEX16ZPGFilter *zpg=NULL;
  REAL4 wc;

  /* set zero pole gain values */
  wc = tan(LAL_PI * filterKnee/samplerate);
  LALCreateCOMPLEX16ZPGFilter(&status, &zpg, 0, 3);
  zpg->poles->data[0].re = wc*sqrt(3.)/2.;
  zpg->poles->data[0].im = wc*0.5;
  zpg->poles->data[1].re = 0.0;
  zpg->poles->data[1].im = wc;
  zpg->poles->data[2].re = -wc*sqrt(3.)/2.;
  zpg->poles->data[2].im = wc*0.5;
  zpg->gain.re = 0.0;
  zpg->gain.im = wc*wc*wc;
  LALWToZCOMPLEX16ZPGFilter( &status, zpg );

  /* create IIR filters */
  iirFilters->filter1Re = NULL;
  iirFilters->filter1Im = NULL;
  LALCreateREAL8IIRFilter( &status, &iirFilters->filter1Re, zpg );
  LALCreateREAL8IIRFilter( &status, &iirFilters->filter1Im, zpg );

  iirFilters->filter2Re = NULL;
  iirFilters->filter2Im = NULL;
  LALCreateREAL8IIRFilter( &status, &iirFilters->filter2Re, zpg );
  LALCreateREAL8IIRFilter( &status, &iirFilters->filter2Im, zpg );

  iirFilters->filter3Re = NULL;
  iirFilters->filter3Im = NULL;
  LALCreateREAL8IIRFilter( &status, &iirFilters->filter3Re, zpg );
  LALCreateREAL8IIRFilter( &status, &iirFilters->filter3Im, zpg );

  /* destroy zpg filter */
  LALDestroyCOMPLEX16ZPGFilter( &status, &zpg );
  LALFree(zpg);
}

/* function to low-pass filter the data using three third order Butterworth IIR filters */
void filter_data(COMPLEX16TimeSeries *data, Filters *iirFilters){
  COMPLEX16 tempData;
  INT4 i=0;

  /* filter the data */
  for(i=0;i<data->data->length;i++){
    tempData.re = 0.;
    tempData.im = 0.;

    tempData.re = LALDIIRFilter(data->data->data[i].re, iirFilters->filter1Re);
    tempData.im = LALDIIRFilter(data->data->data[i].im, iirFilters->filter1Im);

    data->data->data[i].re = LALDIIRFilter(tempData.re, iirFilters->filter2Re);
    data->data->data[i].im = LALDIIRFilter(tempData.im, iirFilters->filter2Im);

    tempData.re = LALDIIRFilter(data->data->data[i].re, iirFilters->filter3Re);
    tempData.im = LALDIIRFilter(data->data->data[i].im, iirFilters->filter3Im);

    data->data->data[i].re = tempData.re;
    data->data->data[i].im = tempData.im;
  }
}

/* function to average the data at one sample rate down to a new sample rate */
COMPLEX16TimeSeries *resample_data(COMPLEX16TimeSeries *data, REAL8Vector *times, INT4Vector
*starts, INT4Vector *stops, REAL8 sampleRate, REAL8 resampleRate, INT4 hetflag){
  COMPLEX16TimeSeries *series;
  INT4 i=0, j=0, k=0, length = floor((round(resampleRate*data->data->length))/sampleRate);
  COMPLEX16 tempData;
  INT4 size = (INT4)round(sampleRate/resampleRate);

  INT4 count=0;

  series = LALCalloc(1, sizeof(*series));
  series->data = XLALCreateCOMPLEX16Vector( length );

  if(hetflag == 0){ /* coarse heterodyne */
    for(i=0;i<data->data->length-size+1;i+=size){
      tempData.re = 0.;
      tempData.im = 0.;

      /* sum data */
      for(j=0;j<size;j++){
        tempData.re += data->data->data[i+j].re;
        tempData.im += data->data->data[i+j].im;
      }

      /* get average */
      series->data->data[count].re = tempData.re/(REAL8)size;
      series->data->data[count].im = tempData.im/(REAL8)size;

      times->data[count] = (REAL8)data->epoch.gpsSeconds + (REAL8)data->epoch.gpsNanoSeconds/1.e9 +
      (1./resampleRate)/2. + ((REAL8)count/resampleRate);

      count++;
    }
  }
  else if(hetflag > 0){ /* need to calculate how many chunks of data at the new sample rate will fit
into each science segment (starts and stops) */
    INT4 duration; /* duration of a science segment */
    INT4 remainder; /* number of data points lost */
    INT4 prevdur=0; /* duration of previous segment */

    count=0;

    for(i=0;i<starts->length;i++){
      if((duration = stops->data[i] - starts->data[i]) < size){
        prevdur += duration;
        continue; /* if segment is smaller than the number of samples needed then skip to next */
      }

      remainder = duration%size;

      for(j=prevdur+floor(remainder/2);j<prevdur+duration-ceil(remainder/2);j+=size){
        tempData.re = 0.;
        tempData.im = 0.;

        for(k=0;k<size;k++){
          tempData.re += data->data->data[j+k].re;
          tempData.im += data->data->data[j+k].im;
        }

        series->data->data[count].re = tempData.re/(REAL8)size;
        series->data->data[count].im = tempData.im/(REAL8)size;

        times->data[count] = times->data[j+size/2] - (1./sampleRate)/2.;
        count++;
      }

      prevdur += duration;

    }
  }

  times = XLALResizeREAL8Vector(times, count);
  series->data = XLALResizeCOMPLEX16Vector(series->data, count);

  /* create time stamps */
  series->deltaT = 1./resampleRate;
  series->epoch = data->epoch;
  series->data->length = count;

  return series;
}

/* read in science segment list file - returns the number of segments */
INT4 get_segment_list(INT4Vector *starts, INT4Vector *stops, CHAR *seglistfile){
  FILE *fp=NULL;
  INT4 i=0;
  long offset;
  CHAR jnkstr[256]; /* junk string to contain comment lines */
  INT4 num, dur; /* variable to contain the segment number and duration */

  if((fp=fopen(seglistfile, "r"))==NULL){
    fprintf(stderr, "Error... can't open science segment list file.\n");
    exit(0);
  }

  /* segment list files have comment lines starting with a # so want to ignore those lines */
  while(!feof(fp)){
    offset = ftell(fp); /* get current position of file stream */

    if(fscanf(fp, "%s", jnkstr) == EOF) /* scan in value and check if == to # */
      break; /* break if there is an extra line at the end of file containing nothing */

    if(strstr(jnkstr, "#")){
      fscanf(fp, "%*[^\n]");   /* if == # then skip to the end of the line */
      continue;
    }
    else{
      fseek(fp, offset, SEEK_SET); /* if line doesn't start with a # then it is data */
      fscanf(fp, "%d%d%d%d", &num, &starts->data[i], &stops->data[i], &dur); /*format is segwizard
type: num starts stops dur */
      i++;
    }
  }

  fclose(fp);

  return i;
}

/* get frame data for partcular science segment */
CHAR *set_frame_files(INT4 *starts, INT4 *stops, FrameCache cache, INT4 numFrames, INT4 *position){
  INT4 i=0;
  INT4 durlock; /* duration of locked segment */
  INT4 tempstart, tempstop;
  INT4 check=0;
  CHAR *smalllist=NULL;

  smalllist = (CHAR *)calloc(MAXLISTLENGTH, sizeof(CHAR));

  durlock = *stops - *starts;
  tempstart = *starts;
  tempstop = *stops;

  /* if lock stretch is long can't read whole thing in at once so only do a bit */
  if(durlock > MAXDATALENGTH)
    tempstop = tempstart + MAXDATALENGTH;

  for(i=0;i<numFrames;i++){
    if(tempstart >= cache.starttime[i] && tempstart < cache.starttime[i]+cache.duration[i] &&
       cache.starttime[i] < tempstop){
      if(check == 0){
        sprintf(smalllist, "%s %d %d ", cache.framelist[i], cache.starttime[i], cache.duration[i]);
        check++;
      }
      else{
        sprintf(smalllist, "%s %s %d %d ", smalllist, cache.framelist[i], cache.starttime[i],
                cache.duration[i]);
      }
      tempstart += cache.duration[i];
    }
    /* break out the loop when we have all the files for the segment */
    else if(cache.starttime[i] > tempstop)
      break;
  }

  if(durlock > MAXDATALENGTH){ /* set starts to its value plus MAXDATALENGTH */
    (*position)--; /* move position counter back one */
    *starts = tempstop;
  }

  if(strlen(smalllist) > MAXLISTLENGTH){
    fprintf(stderr, "Error... small list of frames files is too long.\n");
    exit(0);
  }

  return smalllist;
}

/* function to read in the calibration files and calibrate the data at that (gw) frequency */
void calibrate(COMPLEX16TimeSeries *series, REAL8Vector *datatimes, CalibrationFiles
calfiles, REAL8 frequency){
  FILE *fpcoeff=NULL;
  REAL8 Rfunc, Rphase;
  REAL8 C, Cphase; /* sensing function */
  REAL8 G, Gphase; /* open loop gain */
  INT4 i=0, j=0, k=0, ktemp=0, counter=0;
  
  COMPLEX16 tempData;

  long offset;
  CHAR jnkstr[256]; /* junk string to contain comment lines */
  
  if((fpcoeff = fopen(calfiles.calibcoefficientfile, "r"))==NULL){
    fprintf(stderr, "Error... can't open calibration coefficient file %s.\n\
Assume calibration coefficients are 1 and use the response funtcion.\n",
calfiles.responsefunctionfile);
    /* get response function values */
    get_calibration_values(&Rfunc, &Rphase, calfiles.responsefunctionfile, frequency);

    /* calibrate using response function */
    for(i=0;i<series->data->length;i++){
      tempData = series->data->data[i];
      series->data->data[i].re = Rfunc*(tempData.re*cos(Rphase)-tempData.im*sin(Rphase));
      series->data->data[i].im = Rfunc*(tempData.re*sin(Rphase)+tempData.im*cos(Rphase));
    }
  }
  else{ /* open the open loop gain and sensing function files to use for calibration */
    REAL8 times[MAXCALIBLENGTH];
    REAL8 alpha[MAXCALIBLENGTH], gamma[MAXCALIBLENGTH]; /* gamma = alpha*beta */
    COMPLEX16 Resp;

    /* open sensing function file for reading */
    get_calibration_values(&C, &Cphase, calfiles.sensingfunctionfile, frequency);

    /* get open loop gain values */
    get_calibration_values(&G, &Gphase, calfiles.openloopgainfile, frequency);

    /* read in calibration coefficients */
    while(!feof(fpcoeff)){
      offset = ftell(fpcoeff); /* get current position of file stream */

      if(fscanf(fpcoeff, "%s", jnkstr) == EOF) /* scan in value and check if == to % */
        break;

      if(strstr(jnkstr, "%")){
        fscanf(fpcoeff, "%*[^\n]");   /* if == % then skip to the end of the line */
        continue;
      }
      else{
        fseek(fpcoeff, offset, SEEK_SET); /* if line doesn't start with a % then it is data */
        fscanf(fpcoeff, "%lf%lf%lf", &times[i], &alpha[i], &gamma[i]);
        i++;

        if(i >= MAXCALIBLENGTH){
          fprintf(stderr, "Error... number of lines in calibration coefficient file is greater than %d.\n",
            MAXCALIBLENGTH);
          exit(0);
        }
      }
    }

    fclose(fpcoeff);
    
    /* check times aren't outside range of calib coefficients */
    if(times[0] > datatimes->data[series->data->length-1] || times[i-1] < datatimes->data[0]){
      fprintf(stderr, "Error... calibration coefficients outside range of data.\n");
      exit(0);
    }

    /* calibrate */
    for(j=0;j<series->data->length;j++){
      for(k=ktemp;k<i;k++){
        /* get alpha and gamma values closest to the data, assuming a value a minute */
        if(fabs(times[k] - datatimes->data[j]) <= 30.){
          /* if the coefficients were outside a range then they are bad so don't use */
          if((alpha[k] < ALPHAMIN || alpha[k] > ALPHAMAX) && j < series->data->length-1){
            ktemp = k;
            break;
          }
          else if((alpha[k] < ALPHAMIN || alpha[k] > ALPHAMAX) && j==series->data->length-1)
            break;

          Resp.re = (cos(Cphase) + gamma[k]*G*cos(Gphase - Cphase))/(alpha[k]*C);
          Resp.im = (-sin(Cphase) + gamma[k]*G*sin(Gphase - Cphase))/(alpha[k]*C);

          tempData = series->data->data[j];
          series->data->data[counter].re = tempData.re*Resp.re - tempData.im*Resp.im;
          series->data->data[counter].im = tempData.re*Resp.im + tempData.im*Resp.re;
          datatimes->data[counter] = datatimes->data[j];
          
          counter++;
          ktemp = k;
          break;
        }
      }
    }

    /* resize vectors incase any points have been vetoed by the alpha value cuts */
    series->data = XLALResizeCOMPLEX16Vector(series->data, counter);
    datatimes = XLALResizeREAL8Vector(datatimes, counter);
  }
}

void get_calibration_values(REAL8 *magnitude, REAL8 *phase, CHAR *calibfilename, REAL8
frequency){
  FILE *fp=NULL;
  INT4 i=0;
  long offset;
  CHAR jnkstr[256]; /* junk string to contain comment lines */

  REAL8 freq;

  /* open calibration file for reading */
  if((fp = fopen(calibfilename, "r"))==NULL){
    fprintf(stderr, "Error... can't open file %s.\n", calibfilename);
    exit(0);
  }

  /* calibration files can have lines starting with % at the top so ignore them */
  do{
    offset = ftell(fp); /* get current position of file stream */

    fscanf(fp, "%s", jnkstr); /* scan in value and check if == to % */
    if(strstr(jnkstr, "%")){
      fscanf(fp, "%*[^\n]");   /* if == % then skip to the end of the line */
      continue;
    }
    else{
      fseek(fp, offset, SEEK_SET); /* if line doesn't start with a % then it is data */
      fscanf(fp, "%lf%lf%lf", &freq, magnitude, phase);
      i++;

      if(i >= MAXCALIBLENGTH){
        fprintf(stderr, "Error... number of lines in calibration file is greater than %d.\n",
          MAXCALIBLENGTH);
        exit(0);
      }
    }
  }while(!feof(fp) && freq < frequency); /* stop when we've read in response for our
frequency*/

  fclose(fp);
}

/* function to remove outliers above a certain standard deviation threshold */
void remove_outliers(COMPLEX16TimeSeries *data, REAL8Vector *times, REAL8 stddevthresh){
  COMPLEX16 mean;
  COMPLEX16 stddev;
  INT4 i=0, j=0;

  /* calculate mean - could in reality just assume to be zero */
  mean.re = 0.;
  mean.im = 0.;
  /*for(i=0;i<data->data->length;i++){
    mean.re += data->data->data[i].re;
    mean.im += data->data->data[i].im;
  }
  
  mean.re /= (REAL8)data->data->length;
  mean.im /= (REAL8)data->data->length; */
  
  /* for now assume mean = zero as really large outliers could upset this */
  
  /* calculate standard deviation */
  stddev.re = 0.;
  stddev.im = 0.;
  for(i=0;i<data->data->length;i++){
    stddev.re += (data->data->data[i].re - mean.re)*(data->data->data[i].re - mean.re);
    stddev.im += (data->data->data[i].im - mean.im)*(data->data->data[i].im - mean.im);
  }
  
  stddev.re = sqrt(stddev.re/(REAL8)(data->data->length - 1));
  stddev.im = sqrt(stddev.im/(REAL8)(data->data->length - 1));
  
  /* exclude those points who's absolute value is greater than our stddevthreshold */
  for(i=0;i<data->data->length;i++){
    if(fabs(data->data->data[i].re) < stddev.re*stddevthresh && fabs(data->data->data[i].im) <
stddev.im*stddevthresh){
      data->data->data[j].re = data->data->data[i].re;
      data->data->data[j].im = data->data->data[i].im;
      times->data[j] = times->data[i];
      j++;
    }
  }
  
  /* resize data and times */
  data->data = XLALResizeCOMPLEX16Vector(data->data, j);
  times = XLALResizeREAL8Vector(times, j);
}
    
/* function to show memory useage borrowed from Greg's MakeSFTs code */
/*void printmemuse() {
  pid_t mypid=getpid();
  char commandline[256];
  fflush(NULL);
  sprintf(commandline,"cat /proc/%d/status | /bin/grep Vm | /usr/bin/fmt -140 -u", (int)mypid);
  system(commandline);
  fflush(NULL);
}*/
