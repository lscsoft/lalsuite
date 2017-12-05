/*
*  Copyright (C) 2007 Jolien Creighton, Patrick Brady
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
 * \ingroup lalapps_frametools
 *
 * <dl>
 *
 * <dt>Name</dt><dd>
 * \c lal_animate --- produces an animated display showing the time
 * series output of a selected channel in a lower window, and a simultaneously
 * calculated FFT power spectrum in the upper window</dd>
 *
 * <dt>Synopsis</dt><dd>
 * \c lal_animate [<tt>--help</tt>] [<tt>--channel</tt> \e name]
 * [<tt>--duration</tt> \e secs] [<tt>--epoch</tt> \e sec
 * \e nsec] [<tt>--framedir</tt> \e dirname]
 * [<tt>--highpass</tt> \e freq \e attenuation]
 * [<tt>--numpts</tt> \e npoints] \f$|\f$ xmgr -pipe</dd>
 *
 * <dt>Description</dt><dd>
 * \c lal_animate produces an animated display showing the time
 * series output of a selected channel in a lower window, and a simultaneously
 * calculated FFT power spectrum in the upper window.  The output from
 * this program must be piped into the graphing program \c xmgr.</dd>
 *
 * <dt>Options</dt><dd>
 * <ul>
 * <li><tt>--help</tt>:
 * Print a help message.
 * <li><tt>--channel</tt> \e name:
 * Name of frame channel
 * <li><tt>--duration</tt> \e secs:
 * How many seconds to look at
 * <li><tt>--epoch</tt> \e sec \e nsec:
 * Starting epoch
 * <li><tt>--framedir</tt> \e dirname:
 * Directory containing frame files
 * <li><tt>--highpass</tt> \e freq \e attenuation:
 * High-pass filter parameters
 * <li><tt>--numpts</tt> \e npoints:
 * Points per graph to display
 * </ul></dd>
 *
 * <dt>Example</dt><dd>
 * To run the program,  type:
 * \code
 * lalapps_animate --channel H2:LSC-AS_Q --framedir ./h1 --numpts 16384 \
 * --epoch 693768272 0 --duration 1 --highpass 300 0.01 | xmgr -pipe
 * \endcode
 * This will search in directory <tt>./h1</tt> for frame files containing
 * the channel <tt>H2:LSC-AS_Q</tt> and pipe the data starting at 693768272
 * GPS seconds and 0 GPS nanonseconds to xmgr in segments containing 16384
 * points until 1 seconds of data has been reviewed.  The data is highpass
 * filtered to above 300 Hz with an attenuation of \f$0.1\f$;  the output is
 * shown in \ref f_animate "this figure".
 *
 * \anchor f_animate
 * \image html animate.png "Example of output from \c lalapps_animate program"
 *
 * <dt>Author</dt><dd>
 * Bruce Allen and Patrick Brady
 * </dd>
 * </dl>
 */

#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALFrStream.h>
#include <lal/LALMoment.h>
#include <lal/Units.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALFrameL.h>

#define TRUE       1
#define FALSE      0

#define FRAMEHEARSEEC_ENORM  0
#define FRAMEHEARSEEC_EARG   2

#define FRAMEHEARSEEC_MSGENORM  "Normal exit"
#define FRAMEHEARSEEC_MSGEARG   "Error parsing arguments"

/* Usage format string. */
#define USAGE "Usage: %s [options] | xmgr -pipe\n"\
              "    --help                Print this help message\n" \
              "    --channel name        Name of frame channel\n" \
              "   [--duration secs]      How many seconds to look at\n"\
              "   [--epoch sec nsec]     The starting epoch\n"\
              "   [--framedir dirname]   Directory containing frame files\n"\
              "   [--highpass freq attenuation]  High-pass filter parameters \n"\
              "   [--numpts npoints]     Points per graph displayed\n"


#include <config.h>
#ifndef HAVE_LIBLALFRAME
int main( void )
{
  fputs( "Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr );
  return 77;
}
#else

/* This routine is pipes output into the xmgr graphing program */
static void graphout(float x1,float x2,int thistime, int last) {
   static int count=0;
   printf("&\n");                            /* end of set marker             */
   /* first time we draw the plot */
   if (count==0) {
      printf("@doublebuffer true\n");       /* keeps display from flashing    */
      printf("@s0 color 3\n");              /* IFO graph is green             */
      printf("@view 0.1, 0.1, 0.9, 0.45\n"); /* set the viewport for IFO       */
      printf("@with g1\n");                 /* reset the current graph to FFT */
      printf("@view 0.1, 0.6, 0.9, 0.95\n");/* set the viewport FFT           */
      printf("@with g0\n");                 /* reset the current graph to IFO */
      printf("@world xmin %f\n",x1);        /* set min x                      */
      printf("@world xmax %f\n",x2);        /* set max x                      */
      printf("@autoscale\n");               /* autoscale first time through   */
      printf("@focus off\n");               /* turn off the focus markers     */
      printf("@xaxis label \"t (sec)\"\n"); /* IFO axis label                 */
      printf("@fft(s0, 1)\n");              /* compute the spectrum           */
      printf("@s1 color 2\n");              /* FFT is red                     */
      printf("@move g0.s1 to g1.s0\n");     /* move FFT to graph 1            */
      printf("@with g1\n");                 /* set the focus on FFT           */
      printf("@g1 type logy\n");            /* set FFT to log freq axis       */
      printf("@autoscale\n");               /* autoscale FFT                  */
      printf("@subtitle \"Spectrum\"\n");   /* set the subtitle               */
      printf("@xaxis label \"f (Hz)\"\n");  /* FFT axis label                 */
      printf("@with g0\n");                 /* reset the current graph IFO    */
      printf("@subtitle \"IFO output %d\"\n",thistime);/* set IFO subtitle       */
      count++;/* set IFO subtitle       */
      if (!last) printf("@kill s0\n");      /* kill IFO; ready to read again  */
   }
   else {
      /* other times we redraw the plot */
      printf("@s0 color 3\n");              /* set IFO green                   */
      printf("@fft(s0, 1)\n");              /* FFT it                          */
      printf("@s1 color 2\n");              /* set FFT red                     */
      printf("@move g0.s1 to g1.s0\n");     /* move FFT to graph 1             */
      printf("@subtitle \"IFO output %d\"\n",thistime);/* set IFO subtitle        */
      count++;
      printf("@world xmin %f\n",x1);        /* set min x                       */
      printf("@world xmax %f\n",x2);        /* set max x                       */
      printf("@autoscale yaxes\n");         /* autoscale IFO                   */
      printf("@clear stack\n");             /* clear the stack                 */
      if (!last) printf("@kill s0\n");      /* kill IFO data                   */
      printf("@with g1\n");                 /* switch to FFT                   */
      printf("@g1 type logy\n");            /* set FFT to log freq axis       */
      printf("@clear stack\n");             /* clear stack                     */
      if (!last) printf("@kill s0\n");      /* kill FFT                        */
      printf("@with g0\n");                 /* ready to read IFO again         */
   }
   return;
}



int main( int argc, char *argv[] )
{
    static LALStatus  status;
    LALFrStream         *stream = NULL;
    FrChanIn          channelIn;
    REAL4             itime, numSeconds=0;
    INT4              i, numPoints=4096, inarg = 1;
    CHAR             *dirname;
    REAL4TimeSeries   series;
    LIGOTimeGPS       epoch = {0,0};
    BOOLEAN           epochSet = FALSE;
    BOOLEAN           highpass = FALSE;
    PassBandParamStruc highpassParam;

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
        else if ( !strcmp( argv[inarg], "--highpass" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                highpass = TRUE;
                highpassParam.f2 = atof( argv[inarg++] );
                highpassParam.a2 = atof( argv[inarg++] );
            }else{
                LALPrintError( USAGE, *argv );
                return FRAMEHEARSEEC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--duration" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                numSeconds = atof( argv[inarg++] );
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

    /* Open frame stream */
    LALFrOpen( &status, &stream, dirname, "*.gwf" );

    /* Determine information about the channel */
    series.data = NULL;
    LALFrGetREAL4TimeSeries( &status, &series, &channelIn, stream);
    if (epochSet ){
        series.epoch.gpsSeconds = epoch.gpsSeconds;
        series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
    }
    LALFrSeek(&status, &(series.epoch), stream);

    /* allocate time series */
    LALCreateVector( &status, &series.data, numPoints);

    /* get the data */
    LALFrGetREAL4TimeSeries( &status, &series, &channelIn, stream);

    while ( !(status.statusCode) &&
            (series.epoch.gpsSeconds < epoch.gpsSeconds + (INT4)(numSeconds))){

        /* filter it if the highpass parameters were set */
        if (highpass){
            LALButterworthREAL4TimeSeries(&status, &series, &highpassParam);
        }

        /* print out the contents */
        for (i=0 ; i<32 ; i++ ){
            itime= i * series.deltaT;
            printf("%e\t%e\n", itime,0.0);
        }
        for (i=32 ; i<numPoints ; i++) {
            itime= i * series.deltaT;
            printf("%e\t%e\n", itime, series.data->data[i]);
        }

        /* put out information for the graphing program */
        graphout(0, 0+numPoints * series.deltaT, series.epoch.gpsSeconds,
                (status.statusCode == FRAMESTREAMH_EDONE ||
                 (series.epoch.gpsSeconds + (INT4)(numPoints*series.deltaT)
                  >= epoch.gpsSeconds + (INT4)(numSeconds)) ));

        /* get the data */
        LALFrGetREAL4TimeSeries( &status, &series, &channelIn, stream);
    }

    /* close the frame stream */
    LALFrClose( &status, &stream );

    /*
    * output the sound file
    * sprintf(fname, "%s", channelIn.name);
    * LALTimeSeriesToSound( &status, tSeries, fname, 1);
     */

    /* clean up */
    LALSDestroyVector( &status, &(series.data) );
    LALCheckMemoryLeaks();
    return 0;
}

#endif
