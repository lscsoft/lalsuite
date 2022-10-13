/*
 * Copyright (C) 2007 Gregory Mendell
 * Copyright (C) 2010, 2011, 2016 Bernd Machenschalk
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
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

#include "config.h"

#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 HPf;              /* High pass filtering frequency */
  INT4 T;                 /* SFT duration */
  char *stringT;          /* 12/27/05 gam; string with SFT duration */
  INT4 GPSStart;
  INT4 GPSEnd;
  INT4 makeGPSDirs;        /* 12/27/05 gam; add option to make directories based on gps time */
  INT4 sftVersion;         /* 12/28/05 gam; output SFT version */
  char *commentField;      /* 12/28/05 gam; string comment for version 2 SFTs */
  BOOLEAN htdata;          /* flag that indicates we're doing h(t) data */
  BOOLEAN makeTmpFile;     /* 01/09/06 gam */
  char *FrCacheFile;       /* Frame cache file */
  char *ChannelName;
  char *IFO;               /* 01/14/07 gam */
  char *SFTpath;           /* path to SFT file location */
  char *miscDesc;          /* 12/28/05 gam; string giving misc. part of the SFT description field in the filename */
  INT4 PSSCleaning;	   /* 1=YES and 0=NO*/
  REAL8 PSSCleanHPf;       /* Cut frequency for the bilateral highpass filter. It has to be used only if PSSCleaning is YES.*/
  INT4 PSSCleanExt;        /* Extend the timeseries at the beginning before calculating the autoregressive mean */
  INT4 windowOption;       /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
  REAL8 windowR;
  REAL8 overlapFraction;   /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 1.0). */
  BOOLEAN useSingle;       /* 11/19/05 gam; use single rather than double precision */
  char *frameStructType;   /* 01/10/07 gam */
};

/* GLOBAL VARIABLES */
extern REAL8 winFncRMS;
extern REAL8TimeSeries dataDouble;
extern REAL4TimeSeries dataSingle;

/* FUNCTION PROTOTYPES */
int WindowData( struct CommandLineArgsTag CLA );
int WindowDataTukey2( struct CommandLineArgsTag CLA );
int WindowDataHann( struct CommandLineArgsTag CLA );

int WindowData( struct CommandLineArgsTag CLA )
{
  REAL8 r = CLA.windowR;
  INT4 k, N, kl, kh;
  /* 10/05/12 gam */
  REAL8 win;
  /* This implementation of a Tukey window is describes
     in the Matlab tukey window documentation */

  /* initialize the RMS of the window function */
  winFncRMS = 0.0;

  /* 11/19/05 gam */
  if ( CLA.useSingle ) {
    N = dataSingle.data->length;
    kl = r / 2 * ( N - 1 ) + 1;
    kh = N - r / 2 * ( N - 1 ) + 1;
    for ( k = 1; k < kl; k++ ) {
      /* dataSingle.data->data[k-1]=dataSingle.data->data[k-1]*0.5*( (REAL4)( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ) ); */
      win = 0.5 * ( 1.0 + cos( LAL_TWOPI / r * ( k - 1 ) / ( N - 1 ) - LAL_PI ) );
      dataSingle.data->data[k - 1] *= ( ( REAL4 ) win );
      winFncRMS += win * win;
    }
    for ( k = kh; k <= N; k++ ) {
      /*dataSingle.data->data[k-1]=dataSingle.data->data[k-1]*0.5*( (REAL4)( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ) );*/
      win = 0.5 * ( 1.0 + cos( LAL_TWOPI / r - LAL_TWOPI / r * ( k - 1 ) / ( N - 1 ) - LAL_PI ) );
      dataSingle.data->data[k - 1] *= ( ( REAL4 ) win );
      winFncRMS += win * win;
    }

#if PRINTEXAMPLEDATA
    printf( "\nExample dataSingle values after windowing data in WindowData:\n" );
    printExampleDataSingle();
#endif
  } else {
    N = dataDouble.data->length;
    kl = r / 2 * ( N - 1 ) + 1;
    kh = N - r / 2 * ( N - 1 ) + 1;
    for ( k = 1; k < kl; k++ ) {
      /* dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ); */
      win = 0.5 * ( 1.0 + cos( LAL_TWOPI / r * ( k - 1 ) / ( N - 1 ) - LAL_PI ) );
      dataDouble.data->data[k - 1] *= win;
      winFncRMS += win * win;
    }
    for ( k = kh; k <= N; k++ ) {
      /* dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ); */
      win = 0.5 * ( 1.0 + cos( LAL_TWOPI / r - LAL_TWOPI / r * ( k - 1 ) / ( N - 1 ) - LAL_PI ) );
      dataDouble.data->data[k - 1] *= win;
      winFncRMS += win * win;
    }
#if PRINTEXAMPLEDATA
    printf( "\nExample dataDouble values after windowing data in WindowData:\n" );
    printExampleDataDouble();
#endif
  }

  /* Add to sum of squares of the window function the parts of window which are equal to 1, and then find RMS value*/
  winFncRMS += ( REAL8 )( kh - kl );
  winFncRMS = sqrt( ( winFncRMS / ( ( REAL8 ) N ) ) );

  return 0;
}

/* Same as window function given in lalapps/src/pulsar/make_sfts.c */
int WindowDataTukey2( struct CommandLineArgsTag CLA )
{
  /* Define the parameters to make the window */
  INT4 WINSTART = 4096;
  INT4 WINEND = 8192;
  INT4 WINLEN = ( WINEND - WINSTART );
  /* INT4 i; */
  INT4 i, N;
  REAL8 win;


  /* initialize the RMS of the window function */
  winFncRMS = 0.0;

  if ( CLA.useSingle ) {
    N = dataSingle.data->length;
    /* window data.  Off portion */
    for ( i = 0; i < WINSTART; i++ ) {
      dataSingle.data->data[i] = 0.0;
      dataSingle.data->data[dataSingle.data->length - 1 - i] = 0.0;
    }
    /* window data, smooth turn-on portion */
    for ( i = WINSTART; i < WINEND; i++ ) {
      win = ( ( sin( ( i - WINSTART ) * LAL_PI / ( WINLEN ) - LAL_PI_2 ) + 1.0 ) / 2.0 );
      dataSingle.data->data[i] *= win;
      dataSingle.data->data[dataSingle.data->length - 1 - i]  *= win;
      winFncRMS += 2.0 * win * win;
    }
#if PRINTEXAMPLEDATA
    printf( "\nExample dataSingle values after windowing data in WindowDataTukey2:\n" );
    printExampleDataSingle();
#endif
  } else {
    N = dataDouble.data->length;
    /* window data.  Off portion */
    for ( i = 0; i < WINSTART; i++ ) {
      dataDouble.data->data[i] = 0.0;
      dataDouble.data->data[dataDouble.data->length - 1 - i] = 0.0;
    }
    /* window data, smooth turn-on portion */
    for ( i = WINSTART; i < WINEND; i++ ) {
      win = ( ( sin( ( i - WINSTART ) * LAL_PI / ( WINLEN ) - LAL_PI_2 ) + 1.0 ) / 2.0 );
      dataDouble.data->data[i] *= win;
      dataDouble.data->data[dataDouble.data->length - 1 - i]  *= win;
      winFncRMS += 2.0 * win * win;
    }
#if PRINTEXAMPLEDATA
    printf( "\nExample dataDouble values after windowing data in WindowDataTukey2:\n" );
    printExampleDataDouble();
#endif
  }

  /* Add to sum of squares of the window function the parts of window which are equal to 1, and then find RMS value*/
  winFncRMS += ( REAL8 )( N - 2 * WINEND );
  winFncRMS = sqrt( ( winFncRMS / ( ( REAL8 ) N ) ) );

  return 0;
}

/* Hann window based on Matlab, but with C indexing: w[k] = 0.5*( 1 - cos(2*pi*k/(N-1)) ) k = 0, 1, 2,...N-1 */
int WindowDataHann( struct CommandLineArgsTag CLA )
{
  INT4 k;
  REAL8 win, N, Nm1;
  REAL8 real8TwoPi = 2.0 * ( ( REAL8 )( LAL_PI ) );

  /* initialize the RMS of the window function */
  winFncRMS = 0.0;

  if ( CLA.useSingle ) {
    N = ( ( REAL8 )dataSingle.data->length );
    Nm1 = N - 1;
    for ( k = 0; k < N; k++ ) {
      win = 0.5 * ( 1.0 - cos( real8TwoPi * ( ( REAL8 )( k ) ) / Nm1 ) );
      dataSingle.data->data[k] *= win;
      winFncRMS += win * win;
    }
#if PRINTEXAMPLEDATA
    printf( "\nExample dataSingle values after windowing data in WindowDataHann:\n" );
    printExampleDataSingle();
#endif
  } else {
    N = ( ( REAL8 )dataDouble.data->length );
    Nm1 = N - 1;
    for ( k = 0; k < N; k++ ) {
      win = 0.5 * ( 1.0 - cos( real8TwoPi * ( ( REAL8 )( k ) ) / Nm1 ) );
      dataDouble.data->data[k] *= win;
      winFncRMS += win * win;
    }
#if PRINTEXAMPLEDATA
    printf( "\nExample dataDouble values after windowing data in WindowDataHann:\n" );
    printExampleDataDouble();
#endif
  }

  /* Find RMS value; note that N is REAL8 in this function */
  winFncRMS = sqrt( ( winFncRMS / N ) );

  return 0;
}
