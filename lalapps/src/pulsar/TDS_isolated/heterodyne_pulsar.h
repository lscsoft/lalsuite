/*
  Author: Pitkin, M. D.
  $Id$
*/ 

/* Matt Pitkin 09/02/06 -------------- heterodyne_pulsar.h */

/* header file for heterodyne_pulsar.c */

#ifndef _HETERODYNE_PULSAR_H
#define _HETERODYNE_PULSAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

#include <lal/BinaryPulsarTiming.h>

/* LAL headers */
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h> 

/* frame headers */
#include <FrIO.h>
#include <FrameL.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Usage format string. */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --ifo               name of ifo e.g. L1, H1, H2, G1\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
" --heterodyne-flag   (int) coarse heterodyne 0, fine heterodyne 1,\n\
                      parameter update 2\n"\
" --param-file        name of file containing initial pulsar parameters\n\
                      (.par file)\n"\
" --param-file-update name of file containing updated pulsar parameters\n"\
" --ephem-earth-file  name of file containing the earth ephemeris\n"\
" --ephem-sun-file    name of file containing the sun ephemeris\n"\
" --filter-knee       knee frequency of low-pass filter (don't filter if = 0)\n"\
" --sample-rate       sample rate of input data\n"\
" --resample-rate     sample rate for output data (i.e. 1, 1/60, 8192 etc) \n"\
" --data-file         file containing list of frame files (in frame cache format)\n\
                      or previously heterodyned data file\n"\
" --channel           frame data channel (i.e. LSC-DARM_ERR)\n"\
" --output-dir        directory for output data files\n"\
" --seg-file          name of file containing science segment list\n"\
" --calibrate         if specified calibrate data (no argument)\n"\
" --response-file     name of file containing the response function\n"\
" --coefficient-file  name of file containing the calibration coefficients\n"\
" --sensing-function  name of sensing function file (CAV_GAIN)\n"\
" --open-loop-gain    name of open loop gain file(OLOOP_GAIN)\n"\
" --stddev-thresh     standard deviation threshold veto for outliers\n"\
"\n"

#define MAXLENGTH 10000000 /* max number of lines in heterodyned data file */
#define MAXDATALENGTH 1000 /* maximum length of data to be read from frames */
#define MAXSTRLENGTH 100 /* maximum number of characters in a frame filename */
#define MAXNUMFRAMES 50000 /* maximum number of frame files in input data file */
#define MAXLISTLENGTH 20000 /* maximum length of a list of frames files */
#define MAXCALIBLENGTH 100000 /* maximum length of a calibration file */

#define ALPHAMIN 0.25 /* minimum acceptable value of alpha calib coefficient */
#define ALPHAMAX 2.0 /* maximum acceptable value of alpha calib coefficient */

/* define structures */

typedef struct tagCalibrationFiles{
  CHAR *responsefunctionfile;
  CHAR *calibcoefficientfile;
  CHAR *sensingfunctionfile;
  CHAR *openloopgainfile;
}CalibrationFiles;

/* structure to store data from a lal frame cache as output from LSCdataFind */
typedef struct tagFrameCache{
  CHAR framelist[MAXNUMFRAMES][MAXSTRLENGTH]; /* list of file names in frame cache file */
  INT4 duration[MAXNUMFRAMES]; /* duration of each frame file */
  INT4 starttime[MAXNUMFRAMES]; /* start time of each frame file */
}FrameCache;

typedef struct tagInputParams{
  CHAR ifo[3];
  CHAR pulsar[12];

  INT4 heterodyneflag;
  CHAR paramfile[256];
  CHAR paramfileupdate[256];

  CHAR earthfile[256];
  CHAR sunfile[256];

  REAL8 filterknee;
  REAL8 samplerate;
  REAL8 resamplerate;

  CHAR datafile[256];
  CHAR channel[20];

  CHAR outputdir[256];
  CHAR segfile[256];
  
  INT4 calibrate;
  CalibrationFiles calibfiles;
  
  REAL8 stddevthresh;
  
  INT4 verbose;
}InputParams;

typedef struct tagHeterodyneParams{
  BinaryPulsarParams het;
  BinaryPulsarParams hetUpdate;
  INT4 heterodyneflag;

  LALDetector detector;

  REAL8 samplerate;
  REAL8 timestamp;
  INT4 length;

  CHAR earthfile[256];
  CHAR sunfile[256];
}HeterodyneParams;

typedef struct tagFilters{
  REAL8IIRFilter *filter1Re; /* filters for real and imaginary parts of heterodyed data */
  REAL8IIRFilter *filter1Im;
  REAL8IIRFilter *filter2Re;
  REAL8IIRFilter *filter2Im;
  REAL8IIRFilter *filter3Re;
  REAL8IIRFilter *filter3Im;
}Filters;

/* define functions */
void get_input_args(InputParams *inputParams, int argc, char *argv[]);

void heterodyne_data(COMPLEX16TimeSeries *data, REAL8Vector *times, HeterodyneParams hetParams);

void set_filters(Filters *iirFilters, REAL8 filterKnee, REAL8 samplerate);

void filter_data(COMPLEX16TimeSeries *data, Filters *iirFilters);

COMPLEX16TimeSeries *resample_data(COMPLEX16TimeSeries *data, REAL8Vector *times, INT4Vector
*starts, INT4Vector *stops, REAL8 sampleRate, REAL8 resampleRate, INT4 heterodyneflag);

/* function to extract the frame time and duration from the file name */
void get_frame_times(CHAR *framefile, REAL8 *gpstime, INT4 *duration);

/* reads in a time series from frames */
REAL8TimeSeries *get_frame_data(CHAR *framefile, CHAR *channel, REAL8 time, REAL8 length, INT4
duration);

/* read in science segment list file - returns the number of segments */
INT4 get_segment_list(INT4Vector *starts, INT4Vector *stops, CHAR *seglistfile);

/* get frame data for partcular science segment */
CHAR *set_frame_files(INT4 *starts, INT4 *stops, FrameCache cache, INT4 numFrames, INT4 *position);

/* calibrate data */
void calibrate(COMPLEX16TimeSeries *series, REAL8Vector *datatimes, CalibrationFiles calfiles,
REAL8 frequency);

/* function to extract the correct calibration values from the files */
void get_calibration_values(REAL8 *magnitude, REAL8 *phase, CHAR *calibfilename, REAL8
frequency);

/* function to remove outliers above a certain standard deviation threshold */
void remove_outliers(COMPLEX16TimeSeries *data, REAL8Vector *times, REAL8 stddevthresh);

/*void printmemuse();*/

#ifdef  __cplusplus
}
#endif

#endif /* _HETERODYNE_PULSAR_H */
