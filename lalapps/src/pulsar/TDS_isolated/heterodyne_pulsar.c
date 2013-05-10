/*
*  Copyright (C) 2007 Matt Pitkin
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author M. D. Pitkin
 * \brief
 * lalapps code to perform a coarse/fine heterodyne on LIGO or GEO data given a set of pulsar
 * parameters. The heterodyne is a 2 stage process with the code run for a coarse
 * heterodyne (not taking into account the SSB and BSB time delays, but using the other frequency
 * params) and then rerun for the fine heterodyne (taking into account SSB and BSB time delays). The
 * code can also be used to update heterodyned data using new parameters (old and new parameter files
 * are required) i.e. it will take the difference of the original heterodyne phase and the new phase
 * and reheterodyne with this phase difference.
 */

/* Matt Pitkin - 07/02/06 ------------- heterodyne_pulsar.c */

/* lalapps code to perform a coarse/fine heterodyne on LIGO or GEO data given a set of pulsar
parameters. The heterodyne is a 2 stage process with the code run for a coarse
heterodyne (not taking into account the SSB and BSB time delays, but using the other frequency
params) and then rerun for the fine heterodyne (taking into account SSB and BSB time delays). The
code can also be used to update heterodyned data using new parameters (old and new parameter files
are required) i.e. it will take the difference of the original heterodyne phase and the new phase
and reheterodyne with this phase difference. */

#include "heterodyne_pulsar.h"

/* define a macro to round a number without having to use the C round function */
#define ROUND(a) (floor(a+0.5))

/* track memory usage under linux */
#define TRACKMEMUSE 0

/* routine to track memory usage in different categories */
#if TRACKMEMUSE
void printmemuse() {
   pid_t mypid=getpid();
   char commandline[256];
   fflush(NULL);
   snprintf(commandline, sizeof(commandline), "cat /proc/%d/status | /bin/grep\
 Vm | /usr/bin/fmt -140 -u", (int)mypid);
   system(commandline);
   fflush(NULL);
 }
#endif

/* verbose global variable */
INT4 verbose=0;

int main(int argc, char *argv[]){
  static LALStatus status;

  InputParams inputParams;
  HeterodyneParams hetParams;

  Filters iirFilters;

  FILE *fpin=NULL, *fpout=NULL;
  static FrameCache cache;
  INT4 count=0, frcount=0;

  CHAR outputfile[256]="";
  CHAR channel[128]="";
  CHAR *psrname = NULL;

  INT4Vector *starts=NULL, *stops=NULL; /* science segment start and stop times */
  INT4 numSegs=0;

  FilterResponse *filtresp=NULL; /* variable for the filter response function */

  CHAR *pos=NULL;

  /* set error handler */
  lalDebugLevel = 7;
  XLALSetErrorHandler(XLALAbortErrorHandler);

  #if TRACKMEMUSE
    fprintf(stderr, "Memory use at start of the code:\n"); printmemuse();
  #endif

  /* get input options */
  get_input_args(&inputParams, argc, argv);

  if( inputParams.verbose ) verbose=1;

  hetParams.heterodyneflag = inputParams.heterodyneflag; /* set type of heterodyne */

  /* read in pulsar data */
  XLALReadTEMPOParFile( &hetParams.het, inputParams.paramfile );

  /* set pulsar name - take from par file if available, or if not get from
     command line args */
  if( hetParams.het.jname )
    psrname = XLALStringDuplicate( hetParams.het.jname );
  else if ( hetParams.het.bname )
    psrname = XLALStringDuplicate( hetParams.het.bname );
  else if ( hetParams.het.name )
    psrname = XLALStringDuplicate( hetParams.het.name );
  else if ( inputParams.pulsar )
    psrname = XLALStringDuplicate( inputParams.pulsar );
  else{
    fprintf(stderr, "No pulsar name specified!\n");
    exit(0);
  }

  /* if there is an epoch given manually (i.e. not from the pulsar parameter
     file) then set it here and overwrite any other value - this is used, for
     example, with the pulsar hardware injections in which this should be set
     at 751680013.0 */
  if(inputParams.manualEpoch != 0.){
    hetParams.het.pepoch = inputParams.manualEpoch;
    hetParams.het.posepoch = inputParams.manualEpoch;
  }

  if(verbose){
    fprintf(stderr, "I've read in the pulsar parameters for %s.\n", psrname);
    fprintf(stderr, "alpha = %lf rads, delta = %lf rads.\n", hetParams.het.ra,
      hetParams.het.dec);
    fprintf(stderr, "f0 = %.1lf Hz, f1 = %.1e Hz/s, epoch = %.1lf.\n",
      hetParams.het.f0, hetParams.het.f1, hetParams.het.pepoch);

    fprintf(stderr, "I'm looking for gravitational waves at %.2lf times the \
pulsars spin frequency.\n", inputParams.freqfactor);
  }

  /*if performing fine heterdoyne using same params as coarse */
  if(inputParams.heterodyneflag == 1 || inputParams.heterodyneflag == 3)
    hetParams.hetUpdate = hetParams.het;

  hetParams.samplerate = inputParams.samplerate;

  /* set detector */
  hetParams.detector = *XLALGetSiteInfo( inputParams.ifo );

  if(verbose){  fprintf(stderr, "I've set the detector location for\
 %s.\n", inputParams.ifo); }

  if(inputParams.heterodyneflag == 2 || inputParams.heterodyneflag == 4){ /* if
    updating parameters read in updated par file */
    XLALReadTEMPOParFile( &hetParams.hetUpdate, inputParams.paramfileupdate );

    /* if there is an epoch given manually (i.e. not from the pulsar parameter
       file) then set it here and overwrite any other value */
    if(inputParams.manualEpoch != 0.){
      hetParams.hetUpdate.pepoch = inputParams.manualEpoch;
      hetParams.hetUpdate.posepoch = inputParams.manualEpoch;
    }

    if(verbose){
      fprintf(stderr, "I've read the updated parameters for %s.\n", psrname);
      fprintf(stderr, "alpha = %lf rads, delta = %lf rads.\n",
        hetParams.hetUpdate.ra, hetParams.hetUpdate.dec);
      fprintf(stderr, "f0 = %.1lf Hz, f1 = %.1e Hz/s, epoch = %.1lf.\n",
        hetParams.hetUpdate.f0, hetParams.hetUpdate.f1,
        hetParams.hetUpdate.pepoch);
    }
  }

  if( inputParams.heterodyneflag > 0 ){
    snprintf(hetParams.earthfile, sizeof(hetParams.earthfile), "%s",
      inputParams.earthfile);
    snprintf(hetParams.sunfile, sizeof(hetParams.sunfile), "%s",
      inputParams.sunfile);

    if( inputParams.timeCorrFile != NULL ){
      hetParams.timeCorrFile = XLALStringDuplicate( inputParams.timeCorrFile );

      if( hetParams.hetUpdate.units != NULL ){
        if ( !strcmp(hetParams.hetUpdate.units, "TDB") )
          hetParams.ttype = TIMECORRECTION_TDB; /* use TDB units i.e. TEMPO standard */
        else
          hetParams.ttype = TIMECORRECTION_TCB; /* default to TCB i.e. TEMPO2 standard */
      }
      else /* don't recognise units type, so default to the original code */
        hetParams.ttype = TIMECORRECTION_ORIGINAL;
    }
    else{
      hetParams.timeCorrFile = NULL;
      hetParams.ttype = TIMECORRECTION_ORIGINAL;
    }
  }

  /* get science segment lists - allocate initial memory for starts and stops */
  if( (starts = XLALCreateINT4Vector(1)) == NULL ||
      (stops = XLALCreateINT4Vector(1)) == NULL )
    {  XLALPrintError("Error, allocating segment list memory.\n");  }
  numSegs = get_segment_list(starts, stops, inputParams.segfile,
    inputParams.heterodyneflag);
  if( (starts = XLALResizeINT4Vector(starts, numSegs)) == NULL ||
      (stops = XLALResizeINT4Vector(stops, numSegs)) == NULL )
    {  XLALPrintError("Error, re-allocating segment list memory.\n");  }

  if(verbose){ fprintf(stderr, "I've read in the segment list.\n"); }

  /* open input file */
  if((fpin = fopen(inputParams.datafile, "r")) == NULL){
    fprintf(stderr, "Error... Can't open input data file!\n");
    return 1;
  }

  cache.starttime = NULL;
  cache.duration = NULL;
  cache.framelist = NULL;

  if(inputParams.heterodyneflag == 0 || inputParams.heterodyneflag == 3){
    /* input comes from frame files so read in frame filenames */
    CHAR det[10]; /* detector from cache file */
    CHAR type[256]; /* frame type e.g. RDS_R_L3 - from cache file */
    INT4 cachecount=0, ch=0;

    /* count the number of frame files in the cache file */
    while( (ch = fgetc(fpin) ) != EOF ){
      if( ch == '\n' ) cachecount++;
    }

    /* rewind file pointer */
    rewind(fpin);

    /* allocate memory for frame cache information */
    {
      INT4 ii=0;

      if( (cache.starttime = XLALCalloc(cachecount, sizeof(INT4))) == NULL ||
          (cache.duration = XLALCalloc(cachecount, sizeof(INT4))) == NULL ||
          (cache.framelist = XLALCalloc(cachecount, sizeof(CHAR *))) == NULL )
        {  XLALPrintError("Error allocating frame cache memory.\n");  }

      for( ii=0; ii<cachecount; ii++ ){
        if( (cache.framelist[ii] = XLALCalloc(MAXSTRLENGTH, sizeof(CHAR)))
          == NULL )
          {  XLALPrintError("Error allocating frame list memory.\n");  }
      }
    }

    frcount=0;
    while(fscanf(fpin, "%s%s%d%d file://localhost%s", det, type,
      &cache.starttime[frcount], &cache.duration[frcount],
      cache.framelist[frcount]) != EOF)
      {  frcount++;  }
    fclose(fpin);

    if( frcount != cachecount ){
      fprintf(stderr, "Error... There's been a problem reading in the frame \
data!\n");
      return 1;
    }

    cache.length = cachecount;

    if(verbose){  fprintf(stderr, "I've read in the frame list.\n");  }
  }

  if(inputParams.heterodyneflag == 1 || inputParams.heterodyneflag == 2 ||
    inputParams.heterodyneflag == 4 ){
    if(inputParams.filterknee == 0.){
      fprintf(stderr, "REMINDER: You aren't giving a filter knee frequency from\
 the coarse heterodyne stage! You could be reheterodyning data that has\
 already had the filter response removed, but are you sure this is what you're\
 doing?\n");
    }

    /* calculate the frequency and phase response of the filter used in the
       coarse heterodyne */
    filtresp = create_filter_response( inputParams.filterknee );

    /* reset the filter knee to zero so the filtering is not performed on the
       fine heterodyned data*/
    inputParams.filterknee = 0.;
  }

  /************************BIT THAT DOES EVERYTHING****************************/

  /* set filters - values held for the whole data set so we don't get lots of
     glitches from the filter ringing */
  if(inputParams.filterknee > 0.0){
    set_filters(&iirFilters, inputParams.filterknee, inputParams.samplerate);
    if(verbose){  fprintf(stderr, "I've set up the filters.\n");  }
  }

  /* the outputdir string contains the GPS start and end time of the analysis
     so use this in the filename */
  if((pos = strrchr(inputParams.outputdir, '/'))==NULL){
    fprintf(stderr, "Error... output directory path must contain a final \
directory of the form /GPS_START_TIME-GPS_END_TIME!\n");
    return 1;
  }

  if(inputParams.heterodyneflag == 0){
    snprintf(outputfile, sizeof(outputfile), "%s/coarsehet_%s_%s_%s",
      inputParams.outputdir, psrname, inputParams.ifo, pos+1);
    if(verbose){  fprintf(stderr, "I'm performing a coarse \
heterodyne.\n");  }
  }
  else{
    snprintf(outputfile, sizeof(outputfile), "%s/finehet_%s_%s",
      inputParams.outputdir, psrname, inputParams.ifo);
    if(verbose){  fprintf(stderr, "I'm performing a fine \
heterodyne.\n");  }
  }

  remove(outputfile); /* if output file already exists remove it */
  snprintf(channel, sizeof(channel), "%s", inputParams.channel);

  #if TRACKMEMUSE
    fprintf(stderr, "Memory use before entering main loop:\n"); printmemuse();
  #endif

  /* loop through file and read in data */
  /* if the file is a list of frame files read science segment at a time and
     perform the heterodyne on that set of data - this will only contain a real
     component - plus remember read in double precision for GEO and LIGO h(t)
     and, single precision for uncalibrated LIGO data. If the file is already
     heterodyned it will be complex and we can read the whole file in in one go
     as it should be significantly downsampled */
  do{
    COMPLEX16TimeSeries *data=NULL; /* data for heterodyning */
    COMPLEX16TimeSeries *resampData=NULL; /* resampled data */
    REAL8Vector *times=NULL; /*times of data read from coarse heterodyne file*/
    INT4 i;

    if(inputParams.heterodyneflag == 0 || inputParams.heterodyneflag == 3){
      /* i.e. reading from frame files */
      REAL8 gpstime;
      INT4 duration;
      REAL8TimeSeries *datareal=NULL;
      CHAR *smalllist=NULL; /* list of frame files for a science segment */
      LIGOTimeGPS epochdummy;

      epochdummy.gpsSeconds = 0;
      epochdummy.gpsNanoSeconds = 0;

      /* if the seg list has segment before the start time of the available
         data frame then increment the segment and continue */
      if( stops->data[count] <= cache.starttime[0] ){
        count++;
        continue;
      }
      /* if there are segments after the last available data from then break */
      if( cache.starttime[frcount-1] + cache.duration[frcount-1] <=
          starts->data[count] )
        break;

      if((duration = stops->data[count] - starts->data[count]) > MAXDATALENGTH)
        duration = MAXDATALENGTH; /* if duration of science segment is large
                                     just get part of it */

      fprintf(stderr, "Getting data between %d and %d.\n", starts->data[count],
        starts->data[count]+duration);

      hetParams.timestamp = (REAL8)starts->data[count];
      hetParams.length = inputParams.samplerate * duration;
      gpstime = (REAL8)starts->data[count];

      /* if there was no frame file for that segment move on */
      if((smalllist = set_frame_files(&starts->data[count], &stops->data[count],
        cache, frcount, &count))==NULL){
        /* if there was no frame file for that segment move on */
        fprintf(stderr, "Error... no frame files listed between %d and %d.\n",
          (INT4)gpstime, (INT4)gpstime + duration);

        if(count < numSegs){
          count++;/*if not finished reading in all data try next set of frames*/

          continue;
        }
        else
          break;
      }

      /* make vector (make sure imaginary parts are set to zero) */
      if( (data = XLALCreateCOMPLEX16TimeSeries( "", &epochdummy,
        hetParams.het.f0, 1./inputParams.samplerate, &lalSecondUnit,
        (INT4)inputParams.samplerate * duration )) == NULL )
        {  XLALPrintError("Error allocating data memory.\n");  }

      /* read in frame data */
      if( (datareal = get_frame_data(smalllist, channel, gpstime,
        inputParams.samplerate * duration, duration, inputParams.samplerate,
        inputParams.scaleFac, inputParams.highPass)) == NULL ){
        fprintf(stderr, "Error... could not open frame files between %d and \
%d.\n", (INT4)gpstime, (INT4)gpstime + duration);

        if( count < numSegs ){
          count++;/*if not finished reading in all data try next set of frames*/

          XLALDestroyCOMPLEX16TimeSeries( data );

          XLALFree( smalllist );

          continue;
        }
        else{
          break; /* if at the end of data anyway then break */
          XLALDestroyCOMPLEX16TimeSeries( data );
        }
      }

      /* put data into COMPLEX16 vector and set imaginary parts to zero */
      for( i=0;i<inputParams.samplerate * duration;i++ ){
        data->data->data[i] = (REAL8)datareal->data->data[i];
      }

      XLALDestroyREAL8TimeSeries( datareal );

      XLALFree( smalllist );

      count++;
    }
    else if( inputParams.heterodyneflag == 1 ||
      inputParams.heterodyneflag == 2 ||inputParams.heterodyneflag == 4 ){
      /* i.e. reading from a heterodyned file */
      REAL8 temptime=0.; /* temporary time storage variable */
      LIGOTimeGPS epochdummy;

      epochdummy.gpsSeconds = 0;
      epochdummy.gpsNanoSeconds = 0;

      if( (data = XLALCreateCOMPLEX16TimeSeries( "", &epochdummy,
          hetParams.het.f0, 1./inputParams.samplerate, &lalSecondUnit, 1 ))
          == NULL || (times = XLALCreateREAL8Vector( 1 )) == NULL )
        {  XLALPrintError("Error allocating memory for data.\n");  }
      i=0;

      fprintf(stderr, "Reading heterodyned data from %s.\n",
        inputParams.datafile);

      /* read in file - depends on if file is binary or not */
      if(inputParams.binaryinput){
        INT4 memcount=1;

        do{
          size_t rc;
          REAL8 reVal, imVal;
          rc = fread((void*)&times->data[i], sizeof(REAL8), 1, fpin);
          rc = fread((void*)&reVal, sizeof(REAL8), 1, fpin);
          rc = fread((void*)&imVal, sizeof(REAL8), 1, fpin);

          if( feof(fpin) || rc == 0 ) break;

          if(inputParams.scaleFac > 1.0){
            reVal *= inputParams.scaleFac;
            imVal *= inputParams.scaleFac;
          }

          data->data->data[i] = reVal + I * imVal;

          /* make sure data doesn't overlap previous data */
          if( times->data[i] > temptime ){
            temptime = times->data[i];
            i++;
          }
          else { continue; }

          /* if there is an error during read in then exit */
          if( ferror(fpin) ){
            fprintf(stderr, "Error... problem reading in binary data file!\n");
            exit(1);
          }

          /* dynamically allocate memory 100 lines at a time */
          if( ( i == 1 ) || ( i % 100 == 0 ) ){
            if( (times = XLALResizeREAL8Vector( times, 100*memcount )) == NULL
              || (data = XLALResizeCOMPLEX16TimeSeries( data, 0,
              100*memcount)) == NULL )
              {  XLALPrintError("Error resizing data memory.\n");  }
            memcount++;
          }
        }while( !feof(fpin) );
      }
      else{
        INT4 memcount=1;
        REAL8 reVal, imVal;

        while( fscanf(fpin, "%lf%lf%lf", &times->data[i], &reVal, &imVal) != EOF ){
          if( inputParams.scaleFac > 1.0 ){
            reVal *= inputParams.scaleFac;
            imVal *= inputParams.scaleFac;
          }

          data->data->data[i] = reVal + I * imVal;

          /* make sure data doesn't overlap previous data */
          if( times->data[i] > temptime ){
            temptime = times->data[i];
            i++;
          }
          else continue;

          /* dynamically allocate memory 100 lines at a time */
          if( ( i == 1 ) || ( i % 100 == 0 ) ){
            if( (times = XLALResizeREAL8Vector( times, 100*memcount )) == NULL
              || (data = XLALResizeCOMPLEX16TimeSeries( data, 0,
              100*memcount)) == NULL )
              {  XLALPrintError("Error resizing data memory.\n");  }
            memcount++;
          }
        }
      }

      fclose(fpin);

      hetParams.timestamp = times->data[0]; /* set initial time stamp */

      /* resize vector to actual size */
      if( (data = XLALResizeCOMPLEX16TimeSeries( data, 0, i )) == NULL ||
          (times = XLALResizeREAL8Vector(times, i)) == NULL )
        {  XLALPrintError("Error resizing data memory.\n");  }
      hetParams.length = i;

      if( verbose ) fprintf(stderr, "I've read in the fine heterodyne data.\n");
    }
    else{
      fprintf(stderr, "Error... Heterodyne flag = %d, should be 0, 1, 2, 3 or \
4.\n", inputParams.heterodyneflag);
      return 0;
    }

    XLALGPSSetREAL8(&data->epoch, hetParams.timestamp);

    /* heterodyne data */
    heterodyne_data(data, times, hetParams, inputParams.freqfactor, filtresp);
    if( verbose ){ fprintf(stderr, "I've heterodyned the data.\n"); }

    /* filter data */
    if( inputParams.filterknee > 0. ){/* filter if knee frequency is not zero */
      filter_data(data, &iirFilters);

      if( verbose ){  fprintf(stderr, "I've low pass filtered the \
data at %.2lf Hz\n", inputParams.filterknee);  }
    }

    if( inputParams.heterodyneflag==0 || inputParams.heterodyneflag==3 )
      if( (times = XLALCreateREAL8Vector( data->data->length )) == NULL )
        XLALPrintError("Error creating vector of data times.\n");

    /* resample data and data times */
    resampData = resample_data(data, times, starts, stops,
      inputParams.samplerate, inputParams.resamplerate,
      inputParams.heterodyneflag);
    if( verbose ){  fprintf(stderr, "I've resampled the data from \
%.2lf to %.4lf Hz\n", inputParams.samplerate, inputParams.resamplerate);  }

    XLALDestroyCOMPLEX16TimeSeries( data );

    /*perform outlier removal twice incase very large outliers skew the stddev*/
    if( inputParams.stddevthresh != 0. ){
      INT4 numOutliers=0;
      numOutliers = remove_outliers(resampData, times,
        inputParams.stddevthresh);
      if( verbose ){
        fprintf(stderr, "I've removed %lf%% of data above the threshold %.1lf \
sigma for 1st time.\n",
          100.*(double)numOutliers/(double)resampData->data->length,
          inputParams.stddevthresh);
      }
    }

    /* calibrate */
    if( inputParams.calibrate ){
      calibrate(resampData, times, inputParams.calibfiles,
        inputParams.freqfactor*hetParams.het.f0, inputParams.channel);
      if( verbose ){  fprintf(stderr, "I've calibrated the data at \
%.1lf Hz\n", inputParams.freqfactor*hetParams.het.f0);  }
    }

    /* remove outliers above our threshold */
    if( inputParams.stddevthresh != 0. ){
      INT4 numOutliers = 0;
      numOutliers = remove_outliers(resampData, times,
        inputParams.stddevthresh);
      if( verbose ){
        fprintf(stderr, "I've removed %lf%% of data above the threshold %.1lf \
sigma for 2nd time.\n",
          100.*(double)numOutliers/(double)resampData->data->length,
          inputParams.stddevthresh);
      }
    }

    /* output data */
    if( inputParams.binaryoutput ){
      if((fpout = fopen(outputfile, "ab"))==NULL){
        fprintf(stderr, "Error... can't open output file %s!\n", outputfile);
        return 0;
      }
    }
    else{
      if( (fpout = fopen(outputfile, "a")) == NULL ){
        fprintf(stderr, "Error... can't open output file %s!\n", outputfile);
        return 0;
      }
    }

    /* buffer the output, so that file system is not thrashed when outputing */
    /* buffer will be 1Mb */
    if( setvbuf(fpout, NULL, _IOFBF, 0x100000) )
      fprintf(stderr, "Warning: Unable to set output file buffer!");

    for( i=0;i<(INT4)resampData->data->length;i++ ){
      /* if data has been scaled then undo scaling for output */

      if( inputParams.binaryoutput ){
        size_t rc = 0;
        REAL8 tempreal, tempimag;

        tempreal = creal(resampData->data->data[i]);
        tempimag = cimag(resampData->data->data[i]);

        /*FIXME: maybe add header info to binary output - or output to frames!*/
        /* binary output will be same as ASCII text - time real imag */
        if( inputParams.scaleFac > 1.0 ){
          tempreal /= inputParams.scaleFac;
          tempimag /= inputParams.scaleFac;
        }

        rc = fwrite(&times->data[i], sizeof(REAL8), 1, fpout);
        rc = fwrite(&tempreal, sizeof(REAL8), 1, fpout);
        rc = fwrite(&tempimag, sizeof(REAL8), 1, fpout);

        if( ferror(fpout) || !rc ){
          fprintf(stderr, "Error... problem writing out data to binary \
file!\n");
          exit(1);
        }
      }
      else{
        if( inputParams.scaleFac > 1.0 ){
          fprintf(fpout, "%lf\t%le\t%le\n", times->data[i],
                  creal(resampData->data->data[i])/inputParams.scaleFac,
                  cimag(resampData->data->data[i])/inputParams.scaleFac);
        }
        else{
          fprintf(fpout, "%lf\t%le\t%le\n", times->data[i],
            creal(resampData->data->data[i]), cimag(resampData->data->data[i]));
        }
      }

    }
    if( verbose ){ fprintf(stderr, "I've output the data.\n"); }

    fclose(fpout);
    XLALDestroyCOMPLEX16TimeSeries( resampData );

    XLALDestroyREAL8Vector( times );
  }while( count < numSegs && (inputParams.heterodyneflag==0 ||
    inputParams.heterodyneflag==3) );

  #if TRACKMEMUSE
    fprintf(stderr, "Memory usage after completion of main loop:\n"); printmemuse();
  #endif

  fprintf(stderr, "Heterodyning complete.\n");

  XLALDestroyINT4Vector( stops );
  XLALDestroyINT4Vector( starts );

  if( inputParams.heterodyneflag == 0 || inputParams.heterodyneflag == 3){
    UINT4 ii=0, cachecount=cache.length;

    XLALFree( cache.starttime );
    XLALFree( cache.duration );

    for( ii=0; ii<cachecount; ii++ ) XLALFree(cache.framelist[ii]);

    XLALFree( cache.framelist );
  }

  if( inputParams.filterknee > 0. ){
    LALDestroyREAL8IIRFilter( &status, &iirFilters.filter1Re );
    LALDestroyREAL8IIRFilter( &status, &iirFilters.filter1Im );
    LALDestroyREAL8IIRFilter( &status, &iirFilters.filter2Re );
    LALDestroyREAL8IIRFilter( &status, &iirFilters.filter2Im );
    LALDestroyREAL8IIRFilter( &status, &iirFilters.filter3Re );
    LALDestroyREAL8IIRFilter( &status, &iirFilters.filter3Im );

    if( verbose ){ fprintf(stderr, "I've destroyed all filters.\n"); }
  }

  #if TRACKMEMUSE
    fprintf(stderr, "Memory use at the end of the code:\n"); printmemuse();
  #endif

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
    { "ephem-time-file",          required_argument,  0, 't' },
    { "seg-file",                 required_argument,  0, 'l' },
    { "calibrate",                no_argument,     NULL, 'A' },
    { "response-file",            required_argument,  0, 'R' },
    { "coefficient-file",         required_argument,  0, 'C' },
    { "sensing-function",         required_argument,  0, 'F' },
    { "open-loop-gain",           required_argument,  0, 'O' },
    { "stddev-thresh",            required_argument,  0, 'T' },
    { "freq-factor",              required_argument,  0, 'm' },
    { "scale-factor",             required_argument,  0, 'G' },
    { "high-pass-freq",           required_argument,  0, 'H' },
    { "manual-epoch",             required_argument,  0, 'M' },
    { "binary-input",             no_argument,     NULL, 'B' },
    { "binary-output",            no_argument,     NULL, 'b' },
    { "verbose",                  no_argument,     NULL, 'v' },
    { 0, 0, 0, 0 }
  };

  char args[] = "hi:p:z:f:g:k:s:r:d:c:o:e:S:t:l:R:C:F:O:T:m:G:H:M:ABbv";
  char *program = argv[0];

  /* set defaults */
  inputParams->pulsar = NULL;
  inputParams->filterknee = 0.; /* default is not to filter */
  inputParams->resamplerate = 0.; /* resample to 1 Hz */
  inputParams->samplerate = 0.;
  inputParams->calibrate = 0; /* default is not to calibrate */
  inputParams->verbose = 0; /* default is not to do verbose */
  inputParams->binaryinput = 0; /* default to NOT read in data from a binary
file */
  inputParams->binaryoutput = 0; /* default is to output data as ASCII text */
  inputParams->stddevthresh = 0.; /* default is not to threshold */
  inputParams->calibfiles.calibcoefficientfile = NULL;
  inputParams->calibfiles.sensingfunctionfile = NULL;
  inputParams->calibfiles.openloopgainfile = NULL;
  inputParams->calibfiles.responsefunctionfile = NULL;

  inputParams->freqfactor = 2.0; /* default is to look for gws at twice the
pulsar spin frequency */

  inputParams->scaleFac = 1.0; /* default scaling for calibrated GEO data */
  inputParams->highPass = 0.; /* default to not high-pass GEO data */

  inputParams->manualEpoch = 0.; /* default to zero i.e. it takes the epoch from
the pulsar parameter file */
  /* channel defaults to DARM_ERR */
  snprintf(inputParams->channel, sizeof(inputParams->channel), "DARM_ERR");

  inputParams->timeCorrFile = NULL;

  /* get input arguments */
  while(1){
    int option_index = 0;
    int c;

    c = getopt_long( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch(c){
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error parsing option %s with argument %s\n",
            long_options[option_index].name, optarg );
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'v': /* verbose output */
        inputParams->verbose = 1;
        break;
      case 'i': /* interferometer */
        snprintf(inputParams->ifo, sizeof(inputParams->ifo), "%s", optarg);
        break;
      case 'z': /* heterodyne flag - 0 for coarse, 1 for fine, 2 for update to
                   params, 3 for a one step fine heteroydne (like old code)*/
        inputParams->heterodyneflag = atoi(optarg);
        break;
      case 'p': /* pulsar name */
        inputParams->pulsar = XLALStringDuplicate( optarg );
        break;
      case 'A': /* calibration flag */
        inputParams->calibrate = 1;
        break;
      case 'B': /* input file is binary format */
        inputParams->binaryinput = 1;
        break;
      case 'b': /* output file is to be in binary format */
        inputParams->binaryoutput = 1;
        break;
      case 'f': /* initial heterodyne parameter file */
        snprintf(inputParams->paramfile, sizeof(inputParams->paramfile), "%s",
          optarg);
        break;
      case 'g': /*secondary heterodyne parameter file - for updated parameters*/
        snprintf(inputParams->paramfileupdate,
          sizeof(inputParams->paramfileupdate), "%s", optarg);
        break;
      case 'k': /* low-pass filter knee frequency */
        {/* find if the string contains a / and get its position */
          CHAR *loc=NULL;
          CHAR numerator[10]="", *denominator=NULL;
          INT4 n;

          if((loc = strchr(optarg, '/'))!=NULL){
            n = loc-optarg; /* length of numerator i.e. bit before / */
            XLALStringCopy(numerator, optarg, n+1);

            /*set the denominator i.e. the point after / */
            denominator = XLALStringDuplicate(loc+1);

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
            XLALStringCopy(numerator, optarg, n+1);

            denominator = XLALStringDuplicate(loc+1);

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
            XLALStringCopy(numerator, optarg, n+1);

            denominator = XLALStringDuplicate(loc+1);

            inputParams->resamplerate = atof(numerator)/atof(denominator);
          }
          else
            inputParams->resamplerate = atof(optarg);
        }
        break;
      case 'd': /* file containing list of frame files, or file with previously
                   heterodyned data */
        snprintf(inputParams->datafile, sizeof(inputParams->datafile), "%s",
          optarg);
        break;
      case 'c': /* frame channel */
        snprintf(inputParams->channel, sizeof(inputParams->channel), "%s",
          optarg);
        break;
      case 'o': /* output data directory */
        snprintf(inputParams->outputdir, sizeof(inputParams->outputdir), "%s",
          optarg);
        break;
      case 'e': /* earth ephemeris file */
        snprintf(inputParams->earthfile, sizeof(inputParams->earthfile), "%s",
          optarg);
        break;
      case 'S': /* sun ephemeris file */
        snprintf(inputParams->sunfile, sizeof(inputParams->sunfile), "%s",
          optarg);
        break;
      case 't': /* Einstein delay time correction file */
        inputParams->timeCorrFile = XLALStringDuplicate(optarg);
        break;
      case 'l':
        snprintf(inputParams->segfile, sizeof(inputParams->segfile), "%s",
          optarg);
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
      case 'm': /* this can be defined as a real number or fraction e.g. 2.0 or
                   4/3 */
        {/* find if the string contains a / and get its position */
          CHAR *loc=NULL;
          CHAR numerator[10]="", *denominator=NULL;
          INT4 n;

          if((loc = strchr(optarg, '/'))!=NULL){
            n = loc-optarg; /* length of numerator i.e. bit before / */

            XLALStringCopy(numerator, optarg, n+1);

            denominator = XLALStringDuplicate(loc+1);

            inputParams->freqfactor = atof(numerator)/atof(denominator);
          }
          else
            inputParams->freqfactor = atof(optarg);
        }
        break;
      case 'G':
        inputParams->scaleFac = atof(optarg);
        break;
      case 'H':
        inputParams->highPass = atof(optarg);
        break;
      case 'M':
        inputParams->manualEpoch = atof(optarg);
        break;
      case '?':
        fprintf(stderr, "unknown error while parsing options\n" );
      default:
        fprintf(stderr, "unknown error while parsing options\n" );
    }
  }

  /* set more defaults */
  /*if the sample rate or resample rate hasn't been defined give error */
  if(inputParams->samplerate==0. || inputParams->resamplerate==0.){
    fprintf(stderr, "Error... sample rate and/or resample rate has not been \
defined.\n");
    exit(1);
  }
  else if(inputParams->samplerate < inputParams->resamplerate){
    /* if sr is less than rsr give error */
    fprintf(stderr, "Error... sample rate is less than the resample rate.\n");
    exit(1);
  }

  /* FIXME: This check (fmod) doesn't work if any optimisation flags are set
    when compiling!! */
  /* check that sample rate / resample rate is an integer */
  /*if(fmod(inputParams->samplerate/inputParams->resamplerate, 1.) > 0.){
    fprintf(stderr, "Error... invalid sample rates.\n");
    exit(0);
  }*/

  if(inputParams->freqfactor < 0.){
    fprintf(stderr, "Error... frequency factor must be greater than zero.\n");
    exit(1);
  }

  /* check that we're not trying to set a binary file input for a coarse
     heterodyne */
  if(inputParams->binaryinput){
    if(inputParams->heterodyneflag == 0){
      fprintf(stderr, "Error... binary input should not be set for coarse \
heterodyne!\n");
      exit(1);
    }
  }
}

/* heterodyne data function */
void heterodyne_data(COMPLEX16TimeSeries *data, REAL8Vector *times,
  HeterodyneParams hetParams, REAL8 freqfactor, FilterResponse *filtresp){
  REAL8 phaseCoarse=0., phaseUpdate=0., deltaphase=0.;
  REAL8 t=0., t2=0., tdt=0., tdt2=0., T0=0., T0Update=0., tdt_2=0., tdt2_2=0.;
  REAL8 dtpos=0.; /* time between position epoch and data timestamp */
  INT4 i=0;

  EphemerisData *edat=NULL;
  TimeCorrectionData *tdat=NULL;
  BarycenterInput baryinput, baryinput2;
  EarthState earth, earth2;
  EmissionTime  emit, emit2;

  BinaryPulsarInput binInput, binInput2;
  BinaryPulsarOutput binOutput, binOutput2;

  COMPLEX16 dataTemp, dataTemp2;

  REAL8 df=0., fcoarse=0., ffine=0., resp=0., srate=0.;
  REAL8 filtphase=0.;
  UINT4 position=0, middle=0;

  REAL8 lyr_pc = LAL_PC_SI/LAL_LYR_SI; /* light years per parsec */

  REAL8 om = hetParams.het.wave_om, omu = hetParams.hetUpdate.wave_om;

  /* set the position, frequency and whitening epochs if not already set */
  if(hetParams.het.pepoch == 0. && hetParams.het.posepoch != 0.)
    hetParams.het.pepoch = hetParams.het.posepoch;
  else if(hetParams.het.posepoch == 0. && hetParams.het.pepoch != 0.)
    hetParams.het.posepoch = hetParams.het.pepoch;

  if(hetParams.het.waveepoch == 0. && hetParams.het.nwaves != 0)
    hetParams.het.waveepoch = hetParams.het.pepoch;

  if(hetParams.heterodyneflag == 1 || hetParams.heterodyneflag == 2 ||
    hetParams.heterodyneflag == 4 ){
    if(hetParams.hetUpdate.pepoch == 0. && hetParams.hetUpdate.posepoch != 0.)
      hetParams.hetUpdate.pepoch = hetParams.hetUpdate.posepoch;
    else if(hetParams.hetUpdate.posepoch == 0. &&
      hetParams.hetUpdate.pepoch != 0.)
      hetParams.hetUpdate.posepoch = hetParams.hetUpdate.pepoch;

    if(hetParams.hetUpdate.waveepoch == 0. && hetParams.hetUpdate.nwaves != 0)
      hetParams.hetUpdate.waveepoch = hetParams.hetUpdate.pepoch;

    T0Update = hetParams.hetUpdate.pepoch;
  }

  T0 = hetParams.het.pepoch;

  /* set up ephemeris files */
  if( hetParams.heterodyneflag > 0){
    XLAL_CHECK_VOID( (edat = XLALInitBarycenter( hetParams.earthfile,
                hetParams.sunfile )) != NULL, XLAL_EFUNC );

    /* get files containing Einstein delay correction look-up table */
    if ( hetParams.ttype != TIMECORRECTION_ORIGINAL ){
      XLAL_CHECK_VOID( (tdat = XLALInitTimeCorrections(
        hetParams.timeCorrFile ) ) != NULL, XLAL_EFUNC );
    }

    /* set up location of detector */
    baryinput.site.location[0] = hetParams.detector.location[0]/LAL_C_SI;
    baryinput.site.location[1] = hetParams.detector.location[1]/LAL_C_SI;
    baryinput.site.location[2] = hetParams.detector.location[2]/LAL_C_SI;

    /* set 1/distance using parallax or distance value is given (try parallax
       first) - in 1/secs */
    if( hetParams.heterodyneflag != 2 && hetParams.heterodyneflag != 4 ){ /* not
      using updated params */
      if( hetParams.het.px != 0. ){
        baryinput.dInv = ( 3600. / LAL_PI_180 )*hetParams.het.px /
          (LAL_C_SI*lyr_pc);
      }
      else if( hetParams.het.dist != 0. ){
        baryinput.dInv = 1./(hetParams.het.dist*1e3*LAL_C_SI*lyr_pc);
      }
      else
        baryinput.dInv = 0.; /* no parallax */
    }
    else{                                /* using updated params */
      if( hetParams.hetUpdate.px == 2 ){
        baryinput.dInv = ( 3600. / LAL_PI_180 ) * hetParams.hetUpdate.px /
          (LAL_C_SI*lyr_pc);
      }
      else if( hetParams.hetUpdate.dist != 0. ){
        baryinput.dInv = 1./(hetParams.hetUpdate.dist*1e3*LAL_C_SI*lyr_pc);
      }
      else
        baryinput.dInv = 0.; /* no parallax */
    }
  }

  for(i=0;i<hetParams.length;i++){

/******************************************************************************/
    REAL8 phaseWave = 0.; /* phase of any timing noise whitening parameters */

    /* produce initial heterodyne phase for coarse heterodyne with no time
       delays */
    if(hetParams.heterodyneflag == 0)
      tdt = hetParams.timestamp + (REAL8)i/hetParams.samplerate - T0;
    else if( hetParams.heterodyneflag == 1 || hetParams.heterodyneflag == 2){
      tdt = times->data[i] - T0;
      tdt_2 = times->data[i] - T0Update;
    }
    /* if doing one single heterodyne i.e. het flag = 3 then just calc
       phaseCoarse at all times */
    else if(hetParams.heterodyneflag == 3 || hetParams.heterodyneflag == 4 ){
      /* set up LALBarycenter */
      dtpos = hetParams.timestamp - hetParams.het.posepoch;

      /* set up RA, DEC, and distance variables for LALBarycenter*/
      baryinput.delta = hetParams.het.dec + dtpos*hetParams.het.pmdec;
      baryinput.alpha = hetParams.het.ra +
        dtpos*hetParams.het.pmra/cos(baryinput.delta);

      if( hetParams.heterodyneflag == 3 ) /*get data time */
        t = hetParams.timestamp + (REAL8)i/hetParams.samplerate;
      else
        t = times->data[i];

      XLALGPSSetREAL8(&baryinput.tgps, t);

      XLAL_CHECK_VOID( XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat,
        tdat, hetParams.ttype ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_VOID( XLALBarycenter( &emit, &baryinput, &earth ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );

      /* if binary pulsar add extra time delay */
      if(hetParams.het.model!=NULL){
        /* input SSB time into binary timing function */
        binInput.tb = t + emit.deltaT;
        binInput.earth = earth;

        /* calculate binary time delay */
        XLALBinaryPulsarDeltaT( &binOutput, &binInput, &hetParams.het );

        /* add binary time delay */
        tdt = (t - T0) + emit.deltaT +  binOutput.deltaT;
        tdt_2 = (t - T0Update) + emit.deltaT +  binOutput.deltaT;
      }
      else{
        tdt = t - T0 + emit.deltaT;
        tdt_2 = t - T0Update + emit.deltaT;
      }

      /* check if any timing noise whitening is used */
      if( hetParams.het.nwaves != 0 ){
        /* dtWave only doesn't include binary corrections as they would
         * also be sinusoidal terms */
        REAL8 dtWave = (XLALGPSGetREAL8(&emit.te) -
          hetParams.het.waveepoch)/86400.; /* in days */
        REAL8 tWave = 0.;

        for( INT4 k = 0; k < hetParams.het.nwaves; k++ ){
          tWave += hetParams.het.waveSin[k]*sin(om*(REAL8)(k+1.)*dtWave) +
            hetParams.het.waveCos[k]*cos(om*(REAL8)(k+1.)*dtWave);
        }
        phaseWave = hetParams.het.f0*freqfactor*tWave;
      }
    }

    tdt2 = tdt*tdt; /* tdt^2 for /slightly/ faster computation */
    tdt2_2 = tdt_2*tdt_2;

    /* multiply by 2 to get gw phase */
    phaseCoarse = freqfactor*(hetParams.het.f0*tdt +
      0.5*hetParams.het.f1*tdt2 + (1./6.)*hetParams.het.f2*tdt2*tdt +
      (1./24.)*hetParams.het.f3*tdt2*tdt2 +
      (1./120.)*hetParams.het.f4*tdt2*tdt2*tdt +
      (1./720.)*hetParams.het.f5*tdt2*tdt2*tdt2 +
      (1./5040.)*hetParams.het.f6*tdt2*tdt2*tdt2*tdt +
      (1./40320.)*hetParams.het.f7*tdt2*tdt2*tdt2*tdt2 +
      (1./362880.)*hetParams.het.f8*tdt2*tdt2*tdt2*tdt2*tdt +
      (1./3628800.)*hetParams.het.f9*tdt2*tdt2*tdt2*tdt2*tdt2);
    phaseCoarse += phaseWave;

    fcoarse = freqfactor*(hetParams.het.f0 + hetParams.het.f1*tdt +
      0.5*hetParams.het.f2*tdt2 + (1./6.)*hetParams.het.f3*tdt2*tdt +
      (1./24.)*hetParams.het.f4*tdt2*tdt2 +
      (1./120.)*hetParams.het.f5*tdt2*tdt2*tdt +
      (1./720.)*hetParams.het.f6*tdt2*tdt2*tdt2 +
      (1./5040.)*hetParams.het.f7*tdt2*tdt2*tdt2*tdt +
      (1./40320.)*hetParams.het.f8*tdt2*tdt2*tdt2*tdt2 +
      (1./362880.)*hetParams.het.f9*tdt2*tdt2*tdt2*tdt2*tdt);

    ffine = freqfactor*(hetParams.hetUpdate.f0 + hetParams.hetUpdate.f1*tdt_2 +
      0.5*hetParams.hetUpdate.f2*tdt2_2 + (1./6.)*hetParams.hetUpdate.f3*tdt2_2*tdt_2
      + (1./24.)*hetParams.hetUpdate.f4*tdt2_2*tdt2_2 +
      (1./120.)*hetParams.hetUpdate.f5*tdt2_2*tdt2_2*tdt_2 +
      (1./720.)*hetParams.hetUpdate.f6*tdt2_2*tdt2_2*tdt2_2 +
      (1./5040.)*hetParams.hetUpdate.f7*tdt2_2*tdt2_2*tdt2_2*tdt_2 +
      (1./40320.)*hetParams.hetUpdate.f8*tdt2_2*tdt2_2*tdt2_2*tdt2_2 +
      (1./362880.)*hetParams.hetUpdate.f9*tdt2_2*tdt2_2*tdt2_2*tdt2_2*tdt_2);

/******************************************************************************/
    /* produce second phase for fine heterodyne */
    if( hetParams.heterodyneflag == 1 || hetParams.heterodyneflag == 2 ||
      hetParams.heterodyneflag == 4 ){
      REAL8 tWave1 = 0., tWave2 = 0.;
      phaseWave = 0.;

      /* set up LALBarycenter */
      dtpos = hetParams.timestamp - hetParams.hetUpdate.posepoch;

      /* set up RA, DEC, and distance variables for LALBarycenter*/
      baryinput.delta = hetParams.hetUpdate.dec +
        dtpos*hetParams.hetUpdate.pmdec;
      baryinput.alpha = hetParams.hetUpdate.ra +
        dtpos*hetParams.hetUpdate.pmra/cos(baryinput.delta);

      t = times->data[i]; /* get data time */

      t2 = times->data[i] + 1.; /* just add a second to get the gradient */

      baryinput2 = baryinput;

      XLALGPSSetREAL8(&baryinput.tgps, t);

      XLALGPSSetREAL8(&baryinput2.tgps, t2);

      XLAL_CHECK_VOID( XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat,
        tdat, hetParams.ttype ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_VOID( XLALBarycenter( &emit, &baryinput, &earth ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );

      XLAL_CHECK_VOID( XLALBarycenterEarthNew( &earth2, &baryinput2.tgps, edat,
        tdat, hetParams.ttype ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_VOID( XLALBarycenter( &emit2, &baryinput2, &earth2 ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );

      /* if binary pulsar add extra time delay */
      if(hetParams.hetUpdate.model!=NULL){
        /* input SSB time into binary timing function */
        binInput.tb = t + emit.deltaT;
        binInput2.tb = t2 + emit2.deltaT;
        binInput.earth = binInput2.earth = earth;

        /* calculate binary time delay */
        XLALBinaryPulsarDeltaT( &binOutput, &binInput, &hetParams.hetUpdate );
        XLALBinaryPulsarDeltaT( &binOutput2, &binInput2, &hetParams.hetUpdate );

        /* add binary time delay */
        tdt = (t - T0Update) + emit.deltaT + binOutput.deltaT;
      }
      else{
        tdt = (t - T0Update) + emit.deltaT;
        binOutput.deltaT = 0.;
        binOutput2.deltaT = 0.;
      }

      if( hetParams.hetUpdate.nwaves != 0 ){
        REAL8 dtWave = (XLALGPSGetREAL8(&emit.te) -
          hetParams.hetUpdate.waveepoch)/86400.;

        for( INT4 k = 0; k < hetParams.hetUpdate.nwaves; k++ ){
          tWave1 += hetParams.hetUpdate.waveSin[k]*sin(omu*(REAL8)(k+1.)*dtWave) +
            hetParams.hetUpdate.waveCos[k]*cos(omu*(REAL8)(k+1.)*dtWave);

          tWave2 += hetParams.hetUpdate.waveSin[k]*sin(omu*(REAL8)(k+1.)*(dtWave+1.)) +
            hetParams.hetUpdate.waveCos[k]*cos(omu*(REAL8)(k+1.)*(dtWave+1.));
        }
        phaseWave = hetParams.hetUpdate.f0*freqfactor*tWave1;
      }

      if( filtresp != NULL ){
        /* calculate df  = f*(dt(t2) - dt(t))/(t2 - t) here (t2 - t) is 1 sec */
        df = fcoarse*(emit2.deltaT - emit.deltaT + binOutput2.deltaT -
          binOutput.deltaT + tWave2 - tWave1);
        df += (ffine - fcoarse);

        /*sample rate*/
        srate = (REAL8)filtresp->freqResp->length*filtresp->deltaf;
        middle = (UINT4)srate*FILTERFFTTIME/2;

        /* find filter frequency response closest to df */
        position = middle + (INT4)floor(df/filtresp->deltaf);
      }

      tdt2 = tdt*tdt;
      /* multiply by freqfactor to get gw phase */
      phaseUpdate = freqfactor*(hetParams.hetUpdate.f0*tdt +
        0.5*hetParams.hetUpdate.f1*tdt2 +
        (1./6.)*hetParams.hetUpdate.f2*tdt2*tdt +
        (1./24.)*hetParams.hetUpdate.f3*tdt2*tdt2 +
        (1./120.)*hetParams.hetUpdate.f4*tdt2*tdt2*tdt +
        (1./720.)*hetParams.hetUpdate.f5*tdt2*tdt2*tdt2 +
        (1./5040.)*hetParams.hetUpdate.f6*tdt2*tdt2*tdt2*tdt +
        (1./40320.)*hetParams.hetUpdate.f7*tdt2*tdt2*tdt2*tdt2 +
        (1./362880.)*hetParams.hetUpdate.f8*tdt2*tdt2*tdt2*tdt2*tdt +
        (1./3628800.)*hetParams.hetUpdate.f9*tdt2*tdt2*tdt2*tdt2*tdt2);
      phaseUpdate += phaseWave;
    }

/******************************************************************************/
    dataTemp = data->data->data[i];

    /* perform heterodyne */
    if(hetParams.heterodyneflag == 0  || hetParams.heterodyneflag == 3 ){
      deltaphase = 2.*LAL_PI*fmod(phaseCoarse, 1.);
      dataTemp = creal(dataTemp); /* make sure imaginary part is zero */
    }
    if(hetParams.heterodyneflag == 1 || hetParams.heterodyneflag == 2 ||
      hetParams.heterodyneflag == 4){
      deltaphase = 2.*LAL_PI*fmod(phaseUpdate - phaseCoarse, 1.);

      /* mutliply the data by the filters complex response function to remove
         its effect */
      /* do a linear interpolation between filter response function points */
      if( filtresp != NULL ){
        filtphase = filtresp->phaseResp->data[position] +
          ((filtresp->phaseResp->data[position+1] -
          filtresp->phaseResp->data[position]) / (filtresp->deltaf)) *
          ((REAL8)srate/2. + df - (REAL8)position*filtresp->deltaf);
        filtphase = fmod(filtphase - filtresp->phaseResp->data[middle],
          2.*LAL_PI);

        resp = filtresp->freqResp->data[position] +
          ((filtresp->freqResp->data[position+1] -
          filtresp->freqResp->data[position]) / (filtresp->deltaf)) *
          ((REAL8)srate/2. + df - (REAL8)position*filtresp->deltaf);

        dataTemp2 = dataTemp;

        dataTemp = (1./resp)*(creal(dataTemp2)*cos(filtphase) - cimag(dataTemp2)*sin(filtphase)) +
          I * (1./resp)*(creal(dataTemp2)*sin(filtphase) + cimag(dataTemp2)*cos(filtphase));
      }
    }

    data->data->data[i] = (creal(dataTemp)*cos(-deltaphase) - cimag(dataTemp)*sin(-deltaphase)) +
      I * (creal(dataTemp)*sin(-deltaphase) + cimag(dataTemp)*cos(-deltaphase));
  }

  if(hetParams.heterodyneflag > 0){
    XLALDestroyEphemerisData( edat );

    if ( hetParams.ttype != TIMECORRECTION_ORIGINAL )
      XLALDestroyTimeCorrectionData( tdat );
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
REAL8TimeSeries *get_frame_data(CHAR *framefile, CHAR *channel, REAL8 ttime,
  REAL8 length, INT4 duration, REAL8 samplerate, REAL8 scalefac,
  REAL8 highpass){
  REAL8TimeSeries *dblseries=NULL;

  FrFile *frfile=NULL;
  FrVect *frvect=NULL;

  LIGOTimeGPS epoch;
  INT4 i;

  /* open frame file for reading */
  if((frfile = FrFileINew(framefile)) == NULL)
    return NULL; /* couldn't open frame file */

  XLALGPSSetREAL8(&epoch, ttime);

  /* create data memory */
  if( (dblseries = XLALCreateREAL8TimeSeries( channel, &epoch, 0.,
       1./samplerate, &lalSecondUnit, (INT4)length )) == NULL )
    {  XLALPrintError("Error allocating data times series.\n");  }

  /* read in frame data */
  if((frvect = FrFileIGetV(frfile, channel, ttime, (REAL8)duration))==NULL){
    FrFileIEnd(frfile);
    XLALDestroyREAL8Vector(dblseries->data);
    XLALFree(dblseries);
    return NULL; /* couldn't read frame data */
  }

  /* check if there was missing data within the frame(s) */
  if(frvect->next){
    FrFileIEnd(frfile);
    XLALDestroyREAL8Vector(dblseries->data);
    XLALFree(dblseries);
    FrVectFree(frvect);
    return NULL; /* couldn't read frame data */
  }

  FrFileIEnd(frfile);

  /* fill in vector - checking for frame data type */
  if( frvect->type == FR_VECT_4R ){ /* data is float */
    /* check that data doesn't contain NaNs */
    if(isnan(frvect->dataF[0]) != 0){
      XLALDestroyREAL8Vector(dblseries->data);
      XLALFree(dblseries);
      FrVectFree(frvect);
      return NULL; /* couldn't read frame data */
    }

    for(i=0;i<(INT4)length;i++)
      dblseries->data->data[i] = scalefac*(REAL8)frvect->dataF[i];
  }
  else if( frvect->type == FR_VECT_8R ){ /* data is double */
    /* check that data doesn't contain NaNs */
    if(isnan(frvect->dataD[0]) != 0){
      XLALDestroyREAL8Vector(dblseries->data);
      XLALFree(dblseries);
      FrVectFree(frvect);
      return NULL; /* couldn't read frame data */
    }

    for(i=0;i<(INT4)length;i++)
      dblseries->data->data[i] = scalefac*(REAL8)frvect->dataD[i];
  }
  else{ /* channel name is not recognised */
    fprintf(stderr, "Error... Channel name %s is not recognised as a proper \
channel.\n", channel);
    exit(1); /* abort code */
  }

  /* if a high-pass filter is specified (>0) then filter data */
  if(highpass > 0.){
    PassBandParamStruc highpasspar;

    /* uses 8th order Butterworth, with 10% attenuation */
    highpasspar.nMax = 8;
    highpasspar.f1   = -1;
    highpasspar.a1   = -1;
    highpasspar.f2   = highpass;
    highpasspar.a2   = 0.9; /* this means 10% attenuation at f2 */

    XLALButterworthREAL8TimeSeries( dblseries, &highpasspar );
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
  LAL_CALL( LALCreateCOMPLEX16ZPGFilter(&status, &zpg, 0, 3), &status );
  zpg->poles->data[0] = (wc*sqrt(3.)/2.) + I*(wc*0.5);
  zpg->poles->data[1] = I * wc;
  zpg->poles->data[2] = -(wc*sqrt(3.)/2.) + I*(wc*0.5);
  zpg->gain = I * wc * wc * wc;
  LAL_CALL( LALWToZCOMPLEX16ZPGFilter( &status, zpg ), &status );

  /* create IIR filters */
  iirFilters->filter1Re = NULL;
  iirFilters->filter1Im = NULL;
  LAL_CALL( LALCreateREAL8IIRFilter( &status, &iirFilters->filter1Re, zpg ), &status);
  LAL_CALL( LALCreateREAL8IIRFilter( &status, &iirFilters->filter1Im, zpg ), &status );

  iirFilters->filter2Re = NULL;
  iirFilters->filter2Im = NULL;
  LAL_CALL( LALCreateREAL8IIRFilter( &status, &iirFilters->filter2Re, zpg ), &status );
  LAL_CALL( LALCreateREAL8IIRFilter( &status, &iirFilters->filter2Im, zpg ), &status) ;

  iirFilters->filter3Re = NULL;
  iirFilters->filter3Im = NULL;
  LAL_CALL( LALCreateREAL8IIRFilter( &status, &iirFilters->filter3Re, zpg ), &status );
  LAL_CALL( LALCreateREAL8IIRFilter( &status, &iirFilters->filter3Im, zpg ), &status );

  /* destroy zpg filter */
  LAL_CALL( LALDestroyCOMPLEX16ZPGFilter( &status, &zpg ), &status );
}

/* function to low-pass filter the data using three third order Butterworth IIR
   filters */
void filter_data(COMPLEX16TimeSeries *data, Filters *iirFilters){
  COMPLEX16 tempData;
  INT4 i=0;

  for(i=0;i<(INT4)data->data->length;i++){
    tempData = LALDIIRFilter(creal(data->data->data[i]), iirFilters->filter1Re) +
      I * LALDIIRFilter(cimag(data->data->data[i]), iirFilters->filter1Im);

    data->data->data[i] = LALDIIRFilter(creal(tempData), iirFilters->filter2Re) +
      I * LALDIIRFilter(cimag(tempData), iirFilters->filter2Im);

    tempData = LALDIIRFilter(creal(data->data->data[i]), iirFilters->filter3Re) +
      I * LALDIIRFilter(cimag(data->data->data[i]), iirFilters->filter3Im);

    data->data->data[i] = creal(tempData) + I * cimag(tempData);
  }

}

/* function to average the data at one sample rate down to a new sample rate */
COMPLEX16TimeSeries *resample_data(COMPLEX16TimeSeries *data,
  REAL8Vector *times, INT4Vector *starts, INT4Vector *stops, REAL8 sampleRate,
  REAL8 resampleRate, INT4 hetflag){
  COMPLEX16TimeSeries *series=NULL;
  INT4 i=0, j=0, k=0, length;
  COMPLEX16 tempData;
  INT4 size;

  INT4 count=0;
  LIGOTimeGPS epoch;

  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  if( sampleRate != resampleRate )
    length = (INT4)floor( ROUND( resampleRate*data->data->length )/sampleRate );
  else
    length = data->data->length;

  size = (INT4)ROUND( sampleRate/resampleRate );

  if( (series = XLALCreateCOMPLEX16TimeSeries( "", &epoch, 0., 1./resampleRate,
    &lalSecondUnit, length )) == NULL )
    {  XLALPrintError("Error creating time series for resampled data.\n");  }

  if( hetflag == 0 || hetflag == 3 ){ /* coarse heterodyne */
    for(i=0;i<(INT4)data->data->length-size+1;i+=size){
      tempData = 0.;

      /* sum data */
      for(j=0;j<size;j++){
        tempData += data->data->data[i+j];
      }

      /* get average */
      series->data->data[count] = tempData/(REAL8)size;

      /* take the middle point - if just resampling rather than averaging */
      /*series->data->data[count].re = data->data->data[i + size/2].re;
      series->data->data[count].im = data->data->data[i + size/2].im;*/
      if( sampleRate != resampleRate ){
        times->data[count] = (REAL8)data->epoch.gpsSeconds +
          (REAL8)data->epoch.gpsNanoSeconds/1.e9 + (1./resampleRate)/2. +
          ((REAL8)count/resampleRate);
      }
      else{
        times->data[count] = (REAL8)data->epoch.gpsSeconds +
          (REAL8)data->epoch.gpsNanoSeconds/1.e9 + ((REAL8)count/resampleRate);
      }

      count++;
    }
  }
  else if( hetflag == 1 || hetflag == 2 || hetflag == 4 ){ /* need to calculate
    how many chunks of data at the new sample rate will fit into each science
    segment (starts and stops) */
    INT4 duration=0;  /* duration of a science segment */
    INT4 rremainder=0; /* number of data points lost */
    INT4 prevdur=0;   /* duration of previous segment */
    INT4 frombeg=0, fromend=0;

    count=0;
    j=0;

    for( i=0;i<(INT4)starts->length;i++ ){
      /* find first bit of data within a segment */
      if( starts->data[i] <times->data[j] && stops->data[i] <= times->data[j] ){
        /* if the segmemt is before the jth data point then continue */
        continue;
      }
      if( starts->data[i] < times->data[j] && stops->data[i] > times->data[j] ){
        /* if the segment overlaps the jth data point then use as much as
           possible */
        starts->data[i] = times->data[j] - ((1./sampleRate)/2.);
      }
      if( starts->data[i] < times->data[times->length-1] && stops->data[i] >
        times->data[times->length-1] ){
        /* if starts is before the end of the data but stops is after the end of
           the data then shorten the segment */
        stops->data[i] = times->data[times->length-1] + ((1./sampleRate)/2.);
      }
      if( starts->data[i] >= times->data[times->length-1] ){
        /* segment is outside of data time, so exit */
        fprintf(stderr, "Segment %d to %d is outside the times of the data!\n",
          starts->data[i], stops->data[i]);
        fprintf(stderr, "End resampling\n");
        break;
      }

      while( times->data[j] < starts->data[i] )
        j++;

      starts->data[i] = times->data[j] - ((1./sampleRate)/2.);

      duration = stops->data[i] - starts->data[i];

      /* check duration is not negative */
      if( duration <= 0 )
        continue;

      /* check that data is contiguous within a segment - if not split in two */
      for( k=0;k<(INT4)(duration*sampleRate)-1;k++ ){
        /* check for break in the data */
        if( (times->data[j+k+1] - times->data[j+k]) > 1./sampleRate ){
          /* get duration of segment up until time of split */
          duration = times->data[j+k] - starts->data[i] + ((1./sampleRate)/2.);

          /* restart from this point as if it's a new segment */
          starts->data[i] = times->data[j+k+1] - ((1./sampleRate)/2.);

          /* check that the new point is still in the same segment or not - if
             we are then redo new segment, if not then move on to next segment*/
          if( starts->data[i] < stops->data[i] && duration > 0 )
            i--; /* reset i so it redoes the new segment */

          break;
        }
      }

      if( duration < size ){
        j += (INT4)(duration*sampleRate);
        continue; /* if segment is smaller than the number of samples needed
                     then skip to next */
      }

      rremainder = (INT4)(duration*sampleRate)%size;
      if( sampleRate != resampleRate ){
        frombeg = floor(rremainder/2);
        fromend = ceil(rremainder/2);
      }
      else{
        frombeg = 0;
        fromend = -1;
      }

      prevdur = j;

      for( j=prevdur+frombeg;j<prevdur + (INT4)ROUND(duration*sampleRate) -fromend-1 ; j+=size ){
        tempData = 0.;

        for( k=0;k<size;k++ ){
          tempData += data->data->data[j+k];
        }

        series->data->data[count] = tempData/(REAL8)size;

        if( sampleRate != resampleRate )
          times->data[count] = times->data[j+size/2] - (1./sampleRate)/2.;
        else
          times->data[count] = times->data[j];

        count++;
      }
    }
  }

  if( (INT4)times->length > count )
    if( (times = XLALResizeREAL8Vector( times, count )) == NULL )
      { XLALPrintError("Error resizing resampled times.\n");  }

  if( (INT4)series->data->length > count )
    if( (series = XLALResizeCOMPLEX16TimeSeries( series, 0, count )) == NULL )
      { XLALPrintError("Error resizing resampled data.\n");  }

  /* create time stamps */
  series->deltaT = 1./resampleRate;
  series->epoch = data->epoch;

  return series;
}

/* read in science segment list file - returns the number of segments */
INT4 get_segment_list(INT4Vector *starts, INT4Vector *stops, CHAR *seglistfile,
INT4 heterodyneflag){
  FILE *fp=NULL;
  INT4 i=0;
  long offset;
  CHAR jnkstr[256]; /* junk string to contain comment lines */
  INT4 num, dur; /* variable to contain the segment number and duration */
  INT4 linecount=0; /* number of lines in the segment file */
  INT4 ch=0;

  if((fp=fopen(seglistfile, "r"))==NULL){
    fprintf(stderr, "Error... can't open science segment list file.\n");
    exit(1);
  }

  /* count number of lines in the file */
  while ( (ch = fgetc(fp)) != EOF ){
    if ( ch == '\n' ) /* check for return at end of line */
      linecount++;
  }

  /* allocate memory for vectors */
  if( (starts = XLALResizeINT4Vector( starts, linecount )) == NULL ||
      (stops = XLALResizeINT4Vector( stops, linecount )) == NULL )
    {  XLALPrintError("Error resizing segment lists.\n");  }

  /* rewind file pointer */
  rewind(fp);

  /* segment list files have comment lines starting with a # so want to ignore
     those lines */
  while(!feof(fp)){
    offset = ftell(fp); /* get current position of file stream */

    if(fscanf(fp, "%s", jnkstr) == EOF) /* scan in value and check if == to # */
      break; /* break if there is an extra line at the end of file containing
                nothing */

    if(strstr(jnkstr, "#")){
       /* if == # then skip to the end of the line */
      if ( fscanf(fp, "%*[^\n]") == EOF ) break;
      continue;
    }
    else{
      fseek(fp, offset, SEEK_SET); /* if line doesn't start with a # then it is
                                      data */
      if( fscanf(fp, "%d%d%d%d", &num, &starts->data[i],
                 &stops->data[i], &dur) == EOF ) break;
      /*format is segwizard type: num starts stops dur */

      /* if performing a fine heterodyne remove the first 60 secs at the start
         of a segment to remove the filter response - if the segment is less
         than 60 secs then ignore it */
      if( heterodyneflag == 1 || heterodyneflag == 2 ){
        if(dur > 60){
          starts->data[i] += 60;
          dur -= 60;
        }
        else
          continue;
      }

      i++;
    }
  }

  fclose(fp);

  return i;
}

/* get frame data for partcular science segment */
CHAR *set_frame_files(INT4 *starts, INT4 *stops, FrameCache cache,
  INT4 numFrames, INT4 *position){
  INT4 i=0;
  INT4 durlock; /* duration of locked segment */
  INT4 tempstart, tempstop;
  INT4 check=0;
  CHAR *smalllist=NULL;

  if( (smalllist = XLALCalloc(MAXLISTLENGTH, sizeof(CHAR))) == NULL )
    {  XLALPrintError("Error allocating memory for small frame list.\n");  }

  durlock = *stops - *starts;
  tempstart = *starts;
  tempstop = *stops;

  /*if lock stretch is long can't read whole thing in at once so only do a bit*/
  if( durlock > MAXDATALENGTH )
    tempstop = tempstart + MAXDATALENGTH;

  for( i=0;i<numFrames;i++ ){
    if(tempstart >= cache.starttime[i]
        && tempstart < cache.starttime[i]+cache.duration[i]
        && cache.starttime[i] < tempstop){
      if ( (smalllist = XLALStringAppend(smalllist, cache.framelist[i]))
        == NULL ){
        fprintf(stderr, "Error... something wrong creating frame list with \
XLALStringAppend!\n");
        exit(1);
      }

      /* add a space between frame filenames */
      smalllist = XLALStringAppend(smalllist, " ");

      tempstart += cache.duration[i];
      check++;
    }
    /* break out the loop when we have all the files for the segment */
    else if( cache.starttime[i] > tempstop )
      break;
  }

  if( durlock > MAXDATALENGTH ){/* set starts to its value plus MAXDATALENGTH*/
    (*position)--; /* move position counter back one */
    *starts = tempstop;
  }

  /* if no data was found at all set small list to NULL */
  if( check == 0 ){
    XLALFree(smalllist);
    return NULL;
  }

  if( strlen(smalllist) > MAXLISTLENGTH ){
    fprintf(stderr, "Error... small list of frames files is too long.\n");
    exit(1);
  }

  return smalllist;
}

/* function to read in the calibration files and calibrate the data at that (gw)
frequency */
void calibrate(COMPLEX16TimeSeries *series, REAL8Vector *datatimes,
  CalibrationFiles calfiles, REAL8 frequency, CHAR *channel){
  FILE *fpcoeff=NULL;
  REAL8 Rfunc, Rphase;
  REAL8 C, Cphase; /* sensing function */
  REAL8 G, Gphase; /* open loop gain */
  INT4 i=0, j=0, k=0, ktemp=0, counter=0;

  COMPLEX16 tempData;

  long offset;
  CHAR jnkstr[256]; /* junk string to contain comment lines */
  int rc;

  if(calfiles.calibcoefficientfile == NULL){
    fprintf(stderr, "No calibration coefficient file.\n\
Assume calibration coefficients are 1 and use the response function.\n");
    /* get response function values */
    get_calibration_values(&Rfunc, &Rphase, calfiles.responsefunctionfile,
      frequency);

    /* calibrate using response function */
    for(i=0;i<(INT4)series->data->length;i++){
      tempData = series->data->data[i];
      series->data->data[i] = Rfunc*(creal(tempData)*cos(Rphase)-cimag(tempData)*sin(Rphase)) +
        I * Rfunc*(creal(tempData)*sin(Rphase)+cimag(tempData)*cos(Rphase));
    }
  }
  else{ /* open the open loop gain and sensing function files to use for
           calibration */
    REAL8 *times=NULL;
    REAL8 *alpha=NULL, *ggamma=NULL; /* gamma = alpha*beta */
    COMPLEX16 Resp;

    if((fpcoeff = fopen(calfiles.calibcoefficientfile, "r"))==NULL){
      fprintf(stderr, "Error... can't open calibration coefficient file %s.\n\
Assume calibration coefficients are 1 and use the response function.\n",
        calfiles.calibcoefficientfile);
      exit(1);
    }

    /* open sensing function file for reading */
    get_calibration_values(&C, &Cphase, calfiles.sensingfunctionfile,
      frequency);

    /* get open loop gain values */
    get_calibration_values(&G, &Gphase, calfiles.openloopgainfile, frequency);

    if( ( alpha = XLALCalloc( 100, sizeof(REAL8) ) ) == NULL ||
        ( ggamma = XLALCalloc( 100, sizeof(REAL8) ) ) == NULL ||
        ( times = XLALCalloc( 100, sizeof(REAL8) ) ) == NULL )
      {  XLALPrintError("Error allocating calibration data memory.\n");  }

    /* read in calibration coefficients */
    while(!feof(fpcoeff)){
      offset = ftell(fpcoeff); /* get current position of file stream */

      if(fscanf(fpcoeff, "%s", jnkstr) == EOF){ /* scan in value and check if ==
                                                  to % */
        break;
      }
      if(strstr(jnkstr, "%")){
        rc = fscanf(fpcoeff, "%*[^\n]");   /* if == % then skip to the end of the
                                              line */
        if( rc == EOF ) continue;

        continue;
      }
      else{
        fseek(fpcoeff, offset, SEEK_SET); /* if line doesn't start with a % then
                                             it is data */

        /* dynamically allocate memory (in increments of 100) */
        if( i != 0 && i%100 == 0 ){
          if( ( times = XLALRealloc( times, sizeof(REAL8)*(i+100))) == NULL ){
            fprintf(stderr, "Error... times memory allocation failed\n");
            exit(1);
          }
          if( ( alpha = XLALRealloc( alpha, sizeof(REAL8)*(i+100))) == NULL ){
            fprintf(stderr, "Error... alpha memory allocation failed\n");
            exit(1);
          }
          if( ( ggamma = XLALRealloc( ggamma, sizeof(REAL8)*(i+100))) == NULL){
            fprintf(stderr, "Error... gamma memory allocation failed\n");
            exit(1);
          }
        }

        if(fscanf(fpcoeff, "%lf%lf%lf", &times[i], &alpha[i], &ggamma[i])!= 3){
         fprintf(stderr, "Error... problem reading in values from calibration coefficient file!\n");
         exit(1);
        }
        i++;
      }
    }

    fclose(fpcoeff);

    /* check times aren't outside range of calib coefficients */
    if(times[0] > datatimes->data[series->data->length-1] || times[i-1] <
        datatimes->data[0]){
      fprintf(stderr, "Error... calibration coefficients outside range of data.\n");
      exit(1);
    }

    /* calibrate */
    for(j=0;j<(INT4)series->data->length;j++){
      for(k=ktemp;k<i;k++){
        /* get alpha and gamma values closest to the data, assuming a value a
           minute */
        if(fabs(times[k] - datatimes->data[j]) <= 30.){
          /* if the coefficients were outside a range then they are bad so don't
             use */
          if((alpha[k] < ALPHAMIN || alpha[k] > ALPHAMAX) && j <
              (INT4)series->data->length-1){
            ktemp = k;
            break;
          }
          else if((alpha[k] < ALPHAMIN || alpha[k] > ALPHAMAX) &&
            j==(INT4)series->data->length-1)
            break;

          /* response function for DARM_ERR is
              R(f) = (1 + \gamma*G)/\gamma*C */
          if(strcmp(channel, "LSC-DARM_ERR") == 0){
            Resp = ((cos(Cphase) + ggamma[k]*G*cos(Gphase - Cphase))/(ggamma[k]*C)) +
              I * ((-sin(Cphase) + ggamma[k]*G*sin(Gphase - Cphase))/(ggamma[k]*C));
          }
          /* response function for AS_Q is
              R(f) = (1 + \gamma*G)/\alpha*C */
          else if(strcmp(channel, "LSC-AS_Q") == 0){
            Resp = ((cos(Cphase) + ggamma[k]*G*cos(Gphase - Cphase))/(alpha[k]*C)) +
              I * ((-sin(Cphase) + ggamma[k]*G*sin(Gphase - Cphase))/(alpha[k]*C));
          }
          else{
            fprintf(stderr, "Error... data channel is not set. Give either AS_Q\
 or DARM_ERR!\n");
            exit(1);
          }

          tempData = series->data->data[j];
          series->data->data[counter] = tempData * Resp;
          datatimes->data[counter] = datatimes->data[j];

          counter++;
          ktemp = k;
          break;
        }
      }
    }

    /* free memory */
    XLALFree( times );
    XLALFree( alpha );
    XLALFree( ggamma );

    /*resize vectors incase any points have been vetoed by the alpha valuecuts*/
    if( (series->data = XLALResizeCOMPLEX16Vector(series->data, counter))
        == NULL ||
        (datatimes = XLALResizeREAL8Vector(datatimes, counter)) == NULL )
      {  XLALPrintError("Error resizing calibrated data.\n");  }
  }
}

void get_calibration_values(REAL8 *magnitude, REAL8 *phase, CHAR *calibfilename,
  REAL8 frequency){
  FILE *fp=NULL;
  long offset;
  CHAR jnkstr[256]; /* junk string to contain comment lines */
  REAL8 freq=0.;
  int rc;

  /* open calibration file for reading */
  if(calibfilename == NULL){
    fprintf(stderr, "Error... calibration filename has a null pointer\n");
    exit(1);
  }
  else if((fp = fopen(calibfilename, "r"))==NULL){
    fprintf(stderr, "Error... can't open calibration file %s.\n",
calibfilename);
    exit(1);
  }

  /*calibration files can have lines starting with % at the top so ignore them*/
  do{
    offset = ftell(fp); /* get current position of file stream */

    rc = fscanf(fp, "%s", jnkstr); /* scan in value and check if == to % */
    if(strstr(jnkstr, "%")){
      rc = fscanf(fp, "%*[^\n]");   /* if == % then skip to the end of the line */

      if ( rc == EOF ) continue;
      continue;
    }
    else{
      fseek(fp, offset, SEEK_SET); /* if line doesn't start with a % then it is
                                      data */
      if( fscanf(fp, "%lf%lf%lf", &freq, magnitude, phase) != 3 ){
        fprintf(stderr, "Error... problem reading data from calibration \
file!\n");
        exit(1);
      }
    }
  }while(!feof(fp) && freq < frequency); /* stop when we've read in response for
                                            our frequency*/
  fclose(fp);
}

/* function to remove outliers above a certain standard deviation threshold - it
   returns the number of outliers removed */
INT4 remove_outliers(COMPLEX16TimeSeries *data, REAL8Vector *times,
  REAL8 stddevthresh){
  COMPLEX16 mean;
  COMPLEX16 stddev;
  INT4 i=0, j=0, startlen=(INT4)data->data->length;

  /* calculate mean - could in reality just assume to be zero */
  mean = 0.;
  /*for(i=0;i<data->data->length;i++){
    mean.re += data->data->data[i].re;
    mean.im += data->data->data[i].im;
  }

  mean.re /= (REAL8)data->data->length;
  mean.im /= (REAL8)data->data->length; */

  /* for now assume mean = zero as really large outliers could upset this */

  /* calculate standard deviation */
  stddev = 0.;

  for(i=0;i<(INT4)data->data->length;i++){
    stddev += (creal(data->data->data[i]) - creal(mean))*(creal(data->data->data[i]) - creal(mean)) +
      I * ((cimag(data->data->data[i]) - cimag(mean))*(cimag(data->data->data[i]) - cimag(mean)));
  }

  stddev = sqrt(creal(stddev)/(REAL8)(data->data->length)) +
    I * sqrt(cimag(stddev)/(REAL8)(data->data->length));

  /* exclude those points who's absolute value is greater than our
     stddevthreshold */
  for(i=0;i<(INT4)data->data->length;i++){
    if(fabs(creal(data->data->data[i])) < creal(stddev)*stddevthresh &&
      fabs(cimag(data->data->data[i])) < cimag(stddev)*stddevthresh){
      data->data->data[j] = data->data->data[i];
      times->data[j] = times->data[i];
      j++;
    }
  }

  /* resize data and times */
  if( (data = XLALResizeCOMPLEX16TimeSeries(data, 0, j)) == NULL ||
      (times = XLALResizeREAL8Vector(times, j)) == NULL )
    {  XLALPrintError("Error resizing thresholded data.\n");  }

  return startlen - j;
}

FilterResponse *create_filter_response( REAL8 filterKnee ){
  int i = 0;
  int srate, ttime;

  FilterResponse *filtresp=NULL;
  Filters testFilters;

  COMPLEX16Vector *data=NULL;
  COMPLEX16Vector *fftdata=NULL;

  COMPLEX16FFTPlan *fftplan=NULL;

  REAL8 phase=0., tempphase=0.;
  INT4 count=0;

  /* check filter knee is > 0 otherwise return NULL */
  if( filterKnee <= 0. )
    return NULL;

  srate = 16;
  /* srate = 16384; */ /* sample at 16384 Hz */
  ttime = FILTERFFTTIME; /* have 200 second long data stretch - might need
    longer to increase resolution */

  /**** CREATE SET OF IIR FILTERS ****/
  set_filters(&testFilters, filterKnee, srate);

  /* allocate memory for filtresp */
  if( (filtresp = XLALMalloc(sizeof(FilterResponse))) == NULL )
    {  XLALPrintError("Error allocating memory for filter response.\n");  }

  /* create some data */
  if( (data = XLALCreateCOMPLEX16Vector(ttime*srate)) == NULL )
    { XLALPrintError("Error allocating data for filter response.\n");  }

  /* create impulse and perform filtering */
  for (i = 0;i<srate*ttime; i++){
    REAL8 dataTmpRe = 0., dataTmpIm = 0.;

    if(i==0){
      dataTmpRe = 1.;
      dataTmpIm = 1.;
    }

    dataTmpRe = LALDIIRFilter( dataTmpRe, testFilters.filter1Re );
    dataTmpRe = LALDIIRFilter( dataTmpRe, testFilters.filter2Re );
    dataTmpRe = LALDIIRFilter( dataTmpRe, testFilters.filter3Re );

    dataTmpIm = LALDIIRFilter( dataTmpIm, testFilters.filter1Im );
    dataTmpIm = LALDIIRFilter( dataTmpIm, testFilters.filter2Im );
    dataTmpIm = LALDIIRFilter( dataTmpIm, testFilters.filter3Im );

    data->data[i] = dataTmpRe + I * dataTmpIm;
  }

  /* FFT the data */
  if( (fftplan = XLALCreateForwardCOMPLEX16FFTPlan(srate*ttime, 1)) == NULL ||
      (fftdata = XLALCreateCOMPLEX16Vector(srate*ttime)) == NULL )
    {  XLALPrintError("Error creating FFT plan and data.\n");  }

  XLALCOMPLEX16VectorFFT(fftdata, data, fftplan);

  /* flip vector so that it's in ascending frequency */
  for(i=0;i<srate*ttime/2;i++){
    COMPLEX16 tempdata;

    tempdata = fftdata->data[i];

    fftdata->data[i] = fftdata->data[i+srate*ttime/2];
    fftdata->data[i+srate*ttime/2] = tempdata;
  }

  filtresp->srate = (REAL8)srate;

  if( (filtresp->freqResp = XLALCreateREAL8Vector(srate*ttime)) == NULL ||
      (filtresp->phaseResp = XLALCreateREAL8Vector(srate*ttime)) == NULL )
    {  XLALPrintError("Error allocating filter response vectors.\n");  }

  /* output the frequency and phase response */
  for(i=0;i<srate*ttime;i++){
    filtresp->freqResp->data[i] = sqrt(creal(fftdata->data[i])*creal(fftdata->data[i]) +
      cimag(fftdata->data[i])*cimag(fftdata->data[i]))/sqrt(2.);

    phase = atan2(creal(fftdata->data[i]), cimag(fftdata->data[i]));

    if(i==0){
      filtresp->phaseResp->data[i] = phase;
      continue;
    }
    else if(phase - tempphase < 0.)
      count++;

    filtresp->phaseResp->data[i] = phase + 2.*LAL_PI*(REAL8)count;

    tempphase = phase;
  }

  /* set frequency step */
  filtresp->deltaf = (REAL8)srate/(REAL8)filtresp->freqResp->length;

  XLALDestroyCOMPLEX16Vector(data);
  XLALDestroyCOMPLEX16Vector(fftdata);
  XLALDestroyCOMPLEX16FFTPlan(fftplan);

  return filtresp;
}
