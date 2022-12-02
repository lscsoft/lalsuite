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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/*
  Author: Pitkin, M. D.
*/

/* Matt Pitkin 09/02/06 -------------- heterodyne_pulsar.h */

/* header file for heterodyne_pulsar.c */

#ifndef _HETERODYNE_PULSAR_H
#define _HETERODYNE_PULSAR_H

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>

// definition taken from lal/lib/support/LogPrintf.c
#ifndef HAVE_LOCALTIME_R
#define localtime_r(timep, result) memcpy((result), localtime(timep), sizeof(struct tm))
#endif

#include <lal/BinaryPulsarTiming.h>
#include <lal/LALgetopt.h>
#include <lal/LogPrintf.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/SFTfileIO.h>
#include <lal/LALString.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/LALPulsarVCSInfo.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Usage format string. */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help (-h)              display this message\n"\
" --verbose (-v)           display all error messages\n"\
" --ifo (-i)               name of ifo e.g. L1, H1, H2, G1\n"\
" --pulsar (-p)            name of pulsar e.g. J0534+2200\n"\
" --heterodyne-flag (-z)   (int) coarse heterodyne 0, fine heterodyne 1,\n\
                          parameter update 2, full heterodyne in one go 3,\n\
                          re-heterodyne an already fine heterodyned file 4\n"\
" --param-file (-f)        name of file containing initial pulsar parameters\n\
                          (.par file)\n"\
" --param-file-update (-g) name of file containing updated pulsar parameters\n"\
" --manual-epoch (-M)      a hardwired epoch for the pulsar frequency and\n\
                          position (for use when dealing with hardware\n\
                          injections when this should be set to 751680013.0)\n"\
" --ephem-earth-file (-e)  name of file containing the earth ephemeris\n"\
" --ephem-sun-file (-S)    name of file containing the sun ephemeris\n"\
" --ephem-time-file (-t)   name of file containing the Einstein delay\n\
                          correction values\n"\
" --filter-knee (-k)       knee frequency of low-pass filter (don't filter\n\
                          if = 0)\n"\
" --sample-rate (-s)       sample rate of input data\n"\
" --resample-rate (-r)     sample rate for output data (i.e. 1, 1/60,\n\
                          8192 etc) \n"\
" --data-file (-d)         file containing list of frame files (in frame\n\
                          cache format) or previously heterodyned data\n\
                          file\n"\
" --data-chunk-length (-D)  maximum length of data (in seconds) to be read in\n\
                          at one time from a frame file, overriding\n\
                          MAXDATALENGTH\n"\
" --channel (-c)           frame data channel (i.e. LSC-DARM_ERR)\n"\
" --output-file (-o)       full path and filename for the output data\n"\
" --seg-file (-l)          name of file containing science segment list\n"\
" --calibrate (-A)         if specified calibrate data (no argument)\n"\
" --response-file (-R)     name of file containing the response function\n"\
" --coefficient-file (-C)  name of file containing the calibration\n\
                          coefficients\n"\
" --sensing-function (-F)  name of sensing function file (CAV_GAIN)\n"\
" --open-loop-gain (-O)    name of open loop gain file(OLOOP_GAIN)\n"\
" --stddev-thresh (-T)     standard deviation threshold veto for outliers\n"\
" --freq-factor (-m)       factor which mulitplies the pulsar spin frequency\n\
                          (default 2.0 i.e. a signal from a triaxial pulsar)\n"\
" --scale-factor (-G)      factor to scale the calibrated h(t) data by\n"\
" --high-pass-freq (-H)    high-pass frequency for calibrated h(t) data\n"\
" --binary-input (-B)      read in input data from binary file (for fine and\n\
                          update heterodynes only)\n"\
" --binary-output (-b)     output data to a binary file\n"\
" --legacy-input (-L)      if input heterodyned file is a legacy files produced\n\
                          before the introduction of headers, then set this\n\
                          flag\n"\
" --gzip-output (-Z)       if not outputting in binary format then gzip the\n\
                          output file. This will be done by default if the\n\
                          --output-file specified has the \".gz\" suffix, but\n\
                          if not this suffix will be appended\n"\
" --output-phase (-P)      if set, output the phase evolution to a text file\n\
                          (for debugging purposes)\n"\
"\n"

#define MAXDATALENGTH 256   /* maximum length of data to be read from frames */
#define MAXSTRLENGTH 1024   /* maximum number of characters in a frame filename */
#define MAXLISTLENGTH 25000 /* maximum length of a list of frames files */

#define ALPHAMIN 0.0 /* minimum acceptable value of alpha calib coefficient */
#define ALPHAMAX 2.0 /* maximum acceptable value of alpha calib coefficient */

#define FILTERFFTTIME 200

#define HEADERSIZE 2048 /* number of bytes in header for output files */

/* define structures */

typedef struct tagCalibrationFiles{
  CHAR *responsefunctionfile;
  CHAR *calibcoefficientfile;
  CHAR *sensingfunctionfile;
  CHAR *openloopgainfile;
}CalibrationFiles;

/* structure to store data from a lal frame cache as output from LSCdataFind */
typedef struct tagFrameCache{
  CHAR **framelist; /* list of file names in frame cache file */
  INT4 *duration; /* duration of each frame file */
  INT4 *starttime; /* start time of each frame file */
  UINT4 length;    /* length of cache */
}FrameCache;

typedef struct tagInputParams{
  CHAR ifo[3];
  CHAR *pulsar;

  INT4 heterodyneflag;
  CHAR paramfile[256];
  CHAR paramfileupdate[256];
  REAL8 manualEpoch;

  REAL8 freqfactor;

  CHAR earthfile[256];
  CHAR sunfile[256];
  CHAR *timeCorrFile;

  REAL8 filterknee;
  REAL8 samplerate;
  REAL8 resamplerate;

  CHAR datafile[256];
  CHAR channel[128];
  INT4 datachunklength;

  CHAR outputfile[256];
  CHAR segfile[256];

  INT4 calibrate;
  CalibrationFiles calibfiles;

  REAL8 stddevthresh;

  REAL8 scaleFac;
  REAL8 highPass;

  INT4 verbose;

  INT4 binaryinput;
  INT4 binaryoutput;
  INT4 gzipoutput;
  INT4 legacyinput;
  INT4 outputPhase;
}InputParams;

typedef struct tagHeterodyneParams{
  PulsarParameters *het;
  PulsarParameters *hetUpdate;

  INT4 heterodyneflag;

  LALDetector detector;

  REAL8 samplerate;
  REAL8 timestamp;
  INT4 length;

  CHAR earthfile[256];
  CHAR sunfile[256];
  CHAR *timeCorrFile;
  TimeCorrectionType ttype;
  INT4 outputPhase;
}HeterodyneParams;

typedef struct tagFilters{
  REAL8IIRFilter *filter1Re; /* filters for real and imaginary parts of heterodyed data */
  REAL8IIRFilter *filter1Im;
  REAL8IIRFilter *filter2Re;
  REAL8IIRFilter *filter2Im;
  REAL8IIRFilter *filter3Re;
  REAL8IIRFilter *filter3Im;
}Filters;

typedef struct tagFilterResponse{
  REAL8Vector *freqResp;
  REAL8Vector *phaseResp;

  REAL8 srate; /* sample rate */
  REAL8 deltaf; /* frequency step between successive points */
}FilterResponse;

/* define functions */
void get_input_args(InputParams *inputParams, int argc, char *argv[]);

void heterodyne_data(COMPLEX16TimeSeries *data, REAL8Vector *times, HeterodyneParams hetParams,
REAL8 freqfactor, FilterResponse *filtResp);

void set_filters(Filters *iirFilters, REAL8 filterKnee, REAL8 samplerate);

void filter_data(COMPLEX16TimeSeries *data, Filters *iirFilters);

COMPLEX16TimeSeries *resample_data(COMPLEX16TimeSeries *data, REAL8Vector *times, INT4Vector
*starts, INT4Vector *stops, REAL8 sampleRate, REAL8 resampleRate, INT4 heterodyneflag);

/* function to extract the frame time and duration from the file name */
void get_frame_times(CHAR *framefile, REAL8 *gpstime, INT4 *duration);

/* reads in a time series from frames */
REAL8TimeSeries *get_frame_data(LALCache *framecache, CHAR *channel, REAL8 gpstime,
  REAL8 length, INT4 duration, REAL8 samplerate, REAL8 scalefac,
  REAL8 highpass);

/* read in science segment list file - returns the number of segments */
INT4 get_segment_list(INT4Vector *starts, INT4Vector *stops, CHAR *seglistfile, INT4 heterodyneflag);

/* get frame data for partcular science segment */
LALCache *set_frame_files(INT4 *starts, INT4 *stops, LALCache *cache, INT4 *position, INT4 maxchunklength);

/* calibrate data */
void calibrate(COMPLEX16TimeSeries *series, REAL8Vector *datatimes, CalibrationFiles calfiles,
REAL8 frequency, CHAR *channel);

/* function to extract the correct calibration values from the files */
void get_calibration_values(REAL8 *magnitude, REAL8 *phase, CHAR *calibfilename, REAL8
frequency);

/* function to remove outliers above a certain standard deviation threshold - returns the number
of outliers removed */
INT4 remove_outliers(COMPLEX16TimeSeries *data, REAL8Vector *times, REAL8 stddevthresh);

FilterResponse *create_filter_response( REAL8 filterKnee );

/* free memory for filter response structure */
void destroy_filter_response( FilterResponse *filtresp );

#ifdef  __cplusplus
}
#endif

#endif /* _HETERODYNE_PULSAR_H */
