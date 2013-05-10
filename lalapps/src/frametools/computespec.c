/*
*  Copyright (C) 2007 Patrick Brady
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
 * \author Patrick R Brady
\file
*/

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>
#include <lal/LALMoment.h>
#include <lal/Units.h>
#include <lal/RealFFT.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/FrameCache.h>
#include <lal/LALFrameL.h>

#define TRUE       1
#define FALSE      0

#define FRAMEHEARSEEC_ENORM  0
#define FRAMEHEARSEEC_EARG   2

#define FRAMEHEARSEEC_MSGENORM  "Normal exit"
#define FRAMEHEARSEEC_MSGEARG   "Error parsing arguments"

/* Usage format string. */
#define USAGE "Usage: %s [options]\n"\
              "    Computes an averaged power spectrum according to convention\n"\
              "    in the Conventions Document\n"\
              "    --help                Print this help message\n" \
              "    --channel name        Name of frame channel\n" \
              "    --epoch sec nsec      The starting epoch\n"\
              "    --framedir dirname    Directory containing frame files\n"\
              "    --lowfreq lfreq       High-pass filter parameters \n"\
              "    --numpts npoints      Points per graph displayed\n"\
              "    --numavg navg         Number of segments to average\n"

INT4 lalDebugLevel = LALMSGLVL3;

#include <config.h>
#ifndef HAVE_LIBLALFRAME
int main( void )
{
  fputs( "Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr );
  return 77;
}
#else

int main( int argc, char *argv[] )
{
    static LALStatus  status;
    FrStream         *stream = NULL;
    FrChanIn          channelIn;
    REAL4             lowfreq, highfreq=0, norm;
    REAL4Vector      *spectrum = NULL;
    INT4              i, j, numPoints=4096, inarg = 1, points = 0, navg = 0, nmax;
    CHAR             *dirname;
    REAL4TimeSeries   series;
    LIGOTimeGPS       epoch = {0,0};
    BOOLEAN           epochSet = FALSE;
    BOOLEAN           highpass = FALSE;
    PassBandParamStruc highpassParam;
    FrCache           *frameCache = NULL;

    /* test files are version 4 frames */
    if ( FRAMELIB_VERSION < 4 )
        return 77;

    /* set default values for input */
    dirname = getenv( "LAL_FRAME_PATH" );

    /* Set up for a highpass filter */
    highpassParam.nMax = 4;
    highpassParam.f1 = -1.0;
    highpassParam.a1 = -1.0;

    /*******************************************************************
    * PARSE ARGUMENTS (arg stores the current position)               *
    *******************************************************************/

    if (argc <= 1){
        LALPrintError( USAGE, *argv );
        return 0;
    }

    while ( inarg < argc ) {
        /* Parse output file option. */
        if ( !strcmp( argv[inarg], "--help" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                fprintf(stderr,"Should be a usage message\n");
                exit(0);
            }else{
                LALPrintError( USAGE , *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--epoch" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                epoch.gpsSeconds = atoi( argv[inarg++] );
                epoch.gpsNanoSeconds = atoi( argv[inarg++] );
                epochSet = TRUE;
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--numpts" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                numPoints = atoi( argv[inarg++] );
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--navg" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                navg = atoi( argv[inarg++] );
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--lowfreq" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                highpass = TRUE;
                lowfreq = atof( argv[inarg++] );
                highpassParam.f2 = lowfreq/2.0;
                highpassParam.a2 = 0.1;
                if (argv[inarg] && (argv[inarg][0] != '-'))
                {
                    highfreq = atof( argv[inarg++] );
                }
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--channel" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                channelIn.name = argv[inarg++];
                channelIn.type = ADCDataChannel;
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--framedir" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                dirname = argv[inarg++];
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        /* Check for unrecognized options. */
        else if ( argv[inarg][0] == '-' ) {
            LALPrintError( USAGE, *argv );
            return FRAMEHEARSEEC_EARG;
        }
    } /* End of argument parsing loop. */


    /* Import the frame cache file */
    LALFrCacheImport( &status, &frameCache, dirname);

    /* Open frame stream */
    LALFrCacheOpen( &status, &stream, frameCache );

    /* Determine information about the channel */
    series.data = NULL;
    LALFrGetREAL4TimeSeries( &status, &series, &channelIn, stream);
    if (epochSet ){
        series.epoch.gpsSeconds = epoch.gpsSeconds;
        series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
    }
    LALFrSeek(&status, &(series.epoch), stream);

    /* allocate time series */
    points = (numPoints/2)*(navg+1);
    LALCreateVector( &status, &series.data, points);
    LALSCreateVector( &status, &spectrum, numPoints/2 + 1 );

    /* get the data */
    LALFrGetREAL4TimeSeries( &status, &series, &channelIn, stream);

    {
        RealFFTPlan *pfwd = NULL;   /* FFTW uses a plan to assess best FFT method */
        REAL4Vector hvec;
        COMPLEX8Vector *Hvec = NULL;

        /* Create an FFTW plan for forward REAL FFT */
        LALCreateForwardRealFFTPlan( &status, &pfwd, numPoints, 0);

        /* Create C (complex) vector of length n/2+1 to hold FFT */
        LALCCreateVector( &status, &Hvec, numPoints/2 + 1 );

        /* initialize everything */
        for (i=0 ; i<numPoints/2+1 ;i++){
            spectrum->data[i]=0.0;
        }
        hvec.length=numPoints;
        hvec.data = series.data->data;

        /* filter it if the highpass parameters were set */
        if (highpass){
            fprintf(stderr,"Doing higpass\n");
            LALButterworthREAL4TimeSeries(&status, &series, &highpassParam);
        }

        /* compute the average spectrum */
        for (i=0 ;  i<navg ; i++){
            REAL4 re, im;

            /* do a forward FFT and print results to file dft.dat */
            LALForwardRealFFT( &status, Hvec, &hvec, pfwd );

            for( j=0 ; j<numPoints/2+1; j++){
                re = crealf(Hvec->data[j]);
                im = cimagf(Hvec->data[j]);
                spectrum->data[j] += 2.0*(re*re+im*im)/((REAL4)navg);
            }
        }
        /* Get rid of the FFT plans */
        LALDestroyRealFFTPlan( &status, &pfwd );

        /* get rid of the vectors */
        LALCDestroyVector( &status, &Hvec );

        /* normalize the spectrum */
        norm = series.deltaT/((REAL4)numPoints);
        for( j=0 ; j<numPoints/2+1; j++){
            spectrum->data[j] *= norm;
        }
    }

    /* close the frame stream */
    LALFrClose( &status, &stream );

    /* print out the spectrum */
    if (highpass) {
        nmax=(INT4)(highfreq*((REAL4)numPoints*series.deltaT));
    }
    else
    {
        nmax=numPoints/2+1;
    }
        fprintf(stderr,"nmax = %i\n",nmax);
    for (i=0 ; i<nmax ; i++){
        fprintf(stdout,"%e %e\n",i/((REAL4)numPoints*series.deltaT),spectrum->data[i]);
    }

    /* clean up */
    LALDestroyFrCache( &status, &frameCache);
    LALSDestroyVector( &status, &spectrum );
    LALSDestroyVector( &status, &(series.data) );
    LALCheckMemoryLeaks();
    return 0;
}

#endif
