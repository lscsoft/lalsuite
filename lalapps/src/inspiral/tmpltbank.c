/*
*  Copyright (C) 2007 Alexander Dietz, Duncan Brown, Eirini Messaritaki, Gareth Jones, Benjamin Owen, Patrick Brady, Robert Adam Mercer, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer, Evan Ochsner
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

/*-----------------------------------------------------------------------
 *
 * File Name: tmpltbank.c
 *
 * Author: Brown, D. A.
 *
 *
 *-----------------------------------------------------------------------
 */

/**
 * \file
 * \ingroup lalapps_inspiral
 *
 * <dl>
 *
 * <dt>Name</dt><dd>
 * \c lalapps_tmpltbank --- program to generate inspiral template banks.</dd>
 *
 * <dt>Synopsis</dt><dd>
 * <tt>lalapps_tmpltbank</tt>
 * [<tt>--help</tt>]
 * [<tt>--verbose</tt>]
 * [<tt>--version</tt>]
 * [<tt>--user-tag</tt> <i>usertag</i>]
 * [<tt>--comment</tt> <i>comment</i>]
 * <tt>--gps-start-time</tt> <i>gps_start</i>
 * <tt>--gps-end-time</tt> <i>gps_end</i>
 * [<tt>--pad-data</tt> <i>time_pad</i>]
 * [<tt>--glob-frame-data</tt>]
 * [<tt>--frame-type</tt> <i>type</i>]
 * [<tt>--frame-cache</tt> <i>cache_file</i>]
 * <tt>--calibration-cache</tt> <i>cal_file</i>
 * <tt>--glob-calibration-data</tt>
 * <tt>--channel-name</tt> <i>channel</i>
 * [<tt>--calibrated-data</tt> <i>cal_type</i>]
 * [<tt>--geo-high-pass-freq</tt> <i>geo_freq</i>]
 * [<tt>--geo-high-pass-order</tt> <i>geo_order</i>]
 * [<tt>--geo-high-pass-atten</tt> <i>geo_atten</i>]
 * <tt>--sample-rate</tt> <i>sample_freq</i>
 * <tt>--resample-filter</tt> <i>filter_type</i>
 * [<tt>--disable-high-pass</tt>]
 * [<tt>--enable-high-pass</tt> <i>high_freq</i>]
 * [<tt>--high-pass-order</tt> <i>high_order</i>]
 * [<tt>--high-pass-attenuation</tt> <i>high_atten</i>]
 * <tt>--spectrum-type</tt> <i>spectype</i>
 * [<tt>--dynamic-range-exponent</tt> <i>exp</i>]
 * <tt>--segment-length</tt> <i>seglen</i>
 * [<tt>--number-of-segments</tt> <i>segnum</i>]
 * [<tt>--standard-candle</tt>]
 * [<tt>--candle-snr</tt> <i>candle_snr</i>]
 * [<tt>--candle-mass1</tt> <i>candle_mass1</i>]
 * [<tt>--candle-mass2</tt> <i>candle_mass2</i>]
 * <tt>--low-frequency-cutoff</tt> <i>cutlow</i>
 * <tt>--high-frequency-cutoff</tt> <i>cuthigh</i>
 * [<tt>--minimum-mass</tt> <i>minmass</i>]
 * [<tt>--maximum-mass</tt> <i>maxmass</i>]
 * [<tt>--minimum-psi0</tt> <i>psi0min</i>]
 * [<tt>--maximum-psi0</tt> <i>psi0max</i>]
 * [<tt>--minimum-psi3</tt> <i>psi3min</i>]
 * [<tt>--maximum-psi3</tt> <i>psi3max</i>]
 * [<tt>--maximum-fcut-tmplts</tt> <i>maxTemp</i>]
 * [<tt>--alpha</tt> <i>alpha</i>]
 * <tt>--minimal-match</tt> <i>match</i>
 * <tt>--order</tt> <i>order</i>
 * <tt>--approximant</tt> <i>approx</i>
 * <tt>--space</tt> <i>space</i>
 * [<tt>--write-raw-data</tt>]
 * [<tt>--write-response</tt>]
 * [<tt>--write-spectrum</tt>]
 * [<tt>--write-strain-spectrum</tt>]</dd>
 *
 * <dt>Options</dt><dd>
 * The following command line arguments are available when running tmpltbank.c
 * \\</dd>
 *
 * <dt><tt>--alpha</tt> <i> alpha</i></dt><dd>
 * Set BCV amplitude correction to <i>alpha</i>.</dd>
 *
 * <dt><tt>--approximant</tt> <i> approx</i></dt><dd>
 * Sets the approximant of the waveform to <i>approx</i>. <tt>TaylorT2</tt> is the standard stationary phase frequency domain chirp used in the BNS search. Available parameters: <tt>TaylorT1</tt>, <tt>TaylorT2</tt>, <tt>TaylorT3</tt>, <tt>TaylorF1</tt>, <tt>TaylorF2</tt>, <tt>PadeT1</tt>, <tt>PadeF1</tt>, <tt>EOB</tt>, <tt>BCV</tt>, <tt>SpinTaylorT3</tt>, <tt>BCVSpin</tt>.</dd>
 *
 * <dt><tt>--calibrated-data</tt> <i>type</i></dt><dd>
 * Calibrated data of <i>type</i> <tt>real_4</tt> or <tt>real_8</tt>.</dd>
 *
 * <dt><tt>--calibration-cache</tt> <i>cal_file</i></dt><dd>
 * Obtain calibration from LAL frame cache <i>cal_file</i>.</dd>
 *
 * <dt><tt>--candle-mass1</tt> <i>candle_mass1</i></dt><dd>
 * Mass  <i>candle_mass1</i> of first component in candle binary.  Must be specified is the option <tt>--standard-candle</tt> is set.</dd>
 *
 * <dt><tt>--candle-mass2</tt> <i>candle_mass2</i></dt><dd>
 * Mass <i>candle_mass2</i> of second component in candle binary. Must be specified is the option <tt>--standard-candle</tt> is set.</dd>
 *
 * <dt><tt>--candle-snr</tt> <i>candle_snr</i></dt><dd>
 * Set the signal-to-noise ratio of standard candle to <i>candle_snr</i>. Must be specified is the option <tt>--standard-candle</tt> is set.</dd>
 *
 * <dt><tt>--channel-name</tt> <i>channel</i></dt><dd>
 * Read data from interferometer channel <i>channel</i>.</dd>
 *
 * <dt><tt>--comment</tt> <i>comment</i></dt><dd>
 * Set the process table comment to <i>comment</i>.</dd>
 *
 * <dt><tt>--disable-high-pass</tt></dt><dd>
 * Turn off the IIR highpass filter.  This is an optimistc option. Someday the data will be so good we won't need high pass filtering!  </dd>
 *
 * <dt><tt>--dynamic-range-exponent</tt> <i>exp</i></dt><dd>
 * Set dynamic range scaling to \f${2}^{exp}\f$.</dd>
 *
 * <dt><tt>--enable-high-pass</tt> <i>high_freq</i></dt><dd>
 * High pass data above <i>high_freq</i> Hz using an IIR filter.</dd>
 *
 * <dt><tt>--frame-cache</tt> <i>cache_file</i></dt><dd>
 * This option is used instead of <tt>--glob-frame-data</tt> to read frame data from a frame cache file <i>cache_file</i>. </dd>
 *
 * <dt><tt>--frame-type</tt> <i>type</i></dt><dd>
 * This option specified the type of frames containing the input data. This option must be specified with the <tt>--glob-frame-data</tt> option.???????????</dd>
 *
 * <dt><tt>--geo-high-pass-atten</tt> <i>geo_atten</i></dt><dd>
 * Set the attenuation of the high pass filter to <i>geo_atten</i>. Only if <tt>--calibrated-data</tt> is set to <tt>real_8</tt>.</dd>
 *
 * <dt><tt>--geo-high-pass-freq</tt> <i>geo_freq</i></dt><dd>
 * This sets the high pass filter frequency for GEO data above <i>geo_freq</i> Hz using an IIR filter. Only if <tt>--calibrated-data</tt> is set to <tt>real_8</tt>.</dd>
 *
 * <dt><tt>--geo-high-pass-order</tt> <i>geo_order</i></dt><dd>
 * Set the order of the GEO high pass filter to <i>geo_order</i>. Only if <tt>--calibrated-data</tt> is set to <tt>real_8</tt>.</dd>
 *
 * <dt><tt>--glob-calibration-data</tt></dt><dd>
 * Is this option is specified, the calibration is obtained by globbing in the working directory.?????????</dd>
 *
 * <dt><tt>--glob-frame-data</tt></dt><dd>
 * This option along with <tt>--frame-type</tt>
 * can be used instead of <tt>--frame-cache</tt> to read data stored locally in
 * the working directory.  It finds files of the specified frame type with a *.gwf
 * extension. </dd>
 *
 * <dt><tt>--gps-end-time</tt> <i>gps_end</i></dt><dd>
 * Set the integer part of the GPS time <i>gps_end</i> you want
 * to stop reading data. </dd>
 *
 * <dt><tt>--gps-start-time</tt> <i>gps_start</i></dt><dd>
 * Set the integer part of the GPS time <i>gps_start</i> from which you wish to begin reading data.</dd>
 *
 * <dt><tt>--help</tt></dt><dd> display the help message which gives brief explanations
 * of the command arguments.  </dd>
 *
 * <dt><tt>--high-frequency-cutoff</tt> <i>cuthigh</i></dt><dd>
 * Do not filter above <i>cuthigh</i> Hz.</dd>
 *
 * <dt><tt>--high-pass-attenuation</tt> <i>high_atten</i></dt><dd>
 * Set the attenuation of the high pass filter to <i>high_atten</i>.</dd>
 *
 * <dt><tt>--high-pass-order</tt> <i>high_order</i></dt><dd>
 * Set the order of the high pass filter to <i>high_order</i>.</dd>
 *
 * <dt><tt>--low-frequency-cutoff</tt> <i>cutlow</i></dt><dd>
 * Do not filter below <i>cutlow</i> Hz.</dd>
 *
 * <dt><tt>-maximum-fcut-tmplts</tt> <i> maxTemp</i></dt><dd>
 * Set the maximum number of templates in fcut direction to <i>maxTemp</i>.</dd>
 *
 * <dt><tt>--maximum-mass</tt> <i>maxmass</i></dt><dd>
 * Set maximum component mass of bank to <i>maxmass</i>.</dd>
 *
 * <dt><tt>--maximum-psi0</tt> <i>psi0max</i></dt><dd>
 * Set maximum range of BCV parameter psi0 to <i> psi0max</i>.</dd>
 *
 * <dt><tt>--maximum-psi3</tt> <i>psi3max</i></dt><dd>
 * Set maximum range of BCV parameter psi3 to  <i> psi3max</i>.</dd>
 *
 * <dt><tt>--minimal-match</tt> <i>match</i></dt><dd>
 * Specifies the minimal match <i>match</i> between templates in the
 * bank and all possible signals in the parameter space.</dd>
 *
 * <dt><tt>--minimum-mass</tt> <i>minmass</i></dt><dd>
 * Set minimum component mass of bank to <i>minmass</i>.</dd>
 *
 * <dt><tt>--minimum-psi0</tt> <i>psi0min</i></dt><dd>
 * Set minimum range of BCV parameter psi0 to <i> psi0min</i>.</dd>
 *
 * <dt><tt>--minimum-psi3</tt> <i>psi3min</i></dt><dd>
 * Set minimum range of BCV parameter psi3 to <i> psi3min</i>.</dd>
 *
 * <dt><tt>--number-of-segments</tt> <i>segnum</i></dt><dd>
 * Set number of data segments to <i>segnum</i>.</dd>
 *
 * <dt><tt>--order</tt> <i>order</i></dt><dd>
 * This sets the order of the waveform to <i>order</i>. Usually it is set to <tt>twoPN</tt> (second order post newtonian). Available parameters: <tt>newtonian</tt>, <tt>oneHalfPN</tt>, <tt>onePN</tt>, <tt>onePointFivePN</tt>, <tt>twoPN</tt>, <tt>twoPointFivePN</tt>, <tt>threePN</tt>, <tt>threePointFivePN</tt>.</dd>
 *
 * <dt><tt>--pad-data</tt> <i>time_pad</i></dt><dd>
 * This flag specifies an amount of time <i>time_pad</i> to add to
 * the beginning and end of the input time series data.  Padding the data is
 * necessary because resampling and filtering corrupts these portions.
 * 8 seconds is the accepted choice for this paramenter.  See LAL documentation
 * for a description of resampling and high pass filtering.  </dd>
 *
 * <dt><tt>--resample-filter</tt> <i>filter_type</i></dt><dd>
 * Set resample filter <i>filter_type</i> to <tt>ldas</tt> or <tt>butterworth</tt>. In the normal case the <i>ldas</i> filter is used.</dd>
 *
 * <dt><tt>--sample-rate</tt> <i>sample_freq</i></dt><dd>
 * Specifies the sampling frequency <i>sample_freq</i> at which you
 * want to filter the data downsampling if necessary.</dd>
 *
 * <dt><tt>--segment-length</tt> <i>seglen</i></dt><dd>
 * Set data segment length to <i>seglen</i> points.</dd>
 *
 * <dt><tt>--space</tt> <i> space</i></dt><dd>
 * In order to make the template bank coordinates nice and friendly these parameters are used instead of masses.
 * Usually \c Tau0Tau3 is used. Available parameters:
 * <tt>Tau0Tau2</tt>, <tt>Tau0Tau3</tt>, <tt>Psi0Psi3</tt>.</dd>
 *
 * <dt><tt>--spectrum-type</tt> <i>spec_type</i></dt><dd>
 * Use PSD estimator <i>spec_type</i> <tt>mean</tt> or <tt>median</tt> to choose how the average is calculated. Since the median average is less affected by a loud glitch <tt>median</tt> is used generally.</dd>
 *
 * <dt><tt>--standard-candle</tt></dt><dd>
 * Compute a standard candle from the PSD. In that case the arguments <tt>candle-mass1</tt>, <tt>candle-mass2</tt> and  <tt>candle-snr</tt> must also be specified.</dd>
 *
 * <dt><tt>--verbose</tt></dt><dd> print progress information as the code executes.</dd>
 *
 * <dt><tt>--version</tt></dt><dd> print version information and exit without running
 * the tmpltbank code. </dd>
 *
 * <dt><tt>--user-tag</tt> <i>usertag</i></dt><dd>
 * Set the user tag to the string <i>usertag</i>.
 * This string must not contain spaces or dashes ("-").  This string will appear
 * in the name of the file to which output information is written, and is recorded
 * in the various XML tables within the file.</dd>
 *
 * <dt><tt>--write-raw-data</tt></dt><dd>
 * Write raw data to a frame file.</dd>
 *
 * <dt><tt>--write-response</tt></dt><dd>
 * Write the computed response function to a frame.</dd>
 *
 * <dt><tt>--write-spectrum</tt></dt><dd>
 * Write the uncalibrated psd to a frame.</dd>
 *
 * <dt><tt>--write-strain-spectrum</tt></dt><dd>
 * Write the calibrated strain psd to a text file.
 *
 *  </dd>
 *
 * <dt>Description</dt><dd>
 * \c lalapps_tmpltbank is a stand alone code for generating inspiral
 * template banks for LIGO or GEO data with the LAL bank package.  The code
 * generates a calibrated power spectrum at the specified time for the
 * requested channel and uses this to compute the template bank.
 * The number of templates and the
 * values of the bank parameters in the bank also depend on the minimal
 * match, the
 * minimum and maximum values of mass1 and mass2 (for the BNS search) or the
 * minimum and maximum values of psi0, psi3, the bank-alpha and the number of
 * fcut values (for the BCV search), which are all command-line arguments.
 * Other necessary pieces of information are the approximant and its order and
 * the space that the template bank will be laid on. The output of the code is
 * an xml file and the bank is contained in a \c sngl_inspiral table. The code has
 * also the capability of outputing the raw data, the response function and the
 * calibrated and unclibrated power spectra to frame files.
 * See the LAL bank package
 * documentation for detailed information on the algorithms used to generate the
 * template banks.</dd>
 *
 * <dt>Example</dt><dd>
 *
 * \code
 * lalapps_tmpltbank \
 * --gps-start-time 734357353 --gps-end-time 734358377 \
 * --frame-cache cache/L-734357345-734361107.cache \
 * --segment-length 1048576 --number-of-segments 7 \
 * --pad-data 7 --sample-rate 4096 --resample-filter ldas \
 * --enable-high-pass 5.000000e+01 --spectrum-type median
 * --low-frequency-cutoff 7.000000e+01 --high-frequency-cutoff 2.048000e+03 \
 * --minimum-mass 1.000000e+00  --maximum-mass 3.000000e+00 \
 * --minimal-match 9.700000e-01 --calibration-cache  \
 * /ldas_outgoing/calibration/cache_files/L1-CAL-V03-729273600-734367600.cache \
 * --space Tau0Tau3 --approximant TaylorT1 --order twoPN \
 * --channel-name L1:LSC-AS_Q
 *
 * \endcode</dd>
 *
 * <dt>Author</dt><dd>
 * Duncan Brown and Alexander Dietz</dd>
 * </dl>
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <series.h>
#include <lalappsfrutils.h>

#ifdef LALAPPS_CUDA_ENABLED
#include <cuda_runtime_api.h>
#endif

#include <lal/LALConfig.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALFrStream.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOLwXMLlegacy.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALSimNoise.h>

#include <LALAppsVCSInfo.h>

#include "inspiral.h"

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "tmpltbank"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

/* type of data to analyze */
enum
{
  undefined,
  real_4,
  real_8
} calData = undefined;

/*
 * which type of PSD to use - 'simulated' refers to LALSimNoise PSD functions
 *
 */
enum
{
  specType_mean,
  specType_median,
  specType_simulated,
  specType_undefined
} specType = specType_undefined;

/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbosity of lal function    */

/* parameters used to generate calibrated power spectrum */
LIGOTimeGPS gpsStartTime = { 0, 0 };    /* input data GPS start time    */
LIGOTimeGPS gpsEndTime = { 0, 0 };      /* input data GPS end time      */
INT4  padData = 0;                      /* saftety margin on input data */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
INT4  globFrameData     = 0;            /* glob *.gwf to get frame data */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
CHAR  *frInType         = NULL;         /* type of data frames          */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
CHAR  ifo[3];                           /* two character ifo code       */
CHAR *channelName = NULL;               /* channel string               */
INT4  inputDataLength = 0;              /* number of points in input    */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq     = 0;            /* high pass frequency          */
INT4   highPassOrder    = -1;           /* order of the td iir filter   */
REAL4  highPassAtten    = -1;           /* attenuation of the td filter */
REAL4  fLow             = -1;           /* low frequency cutoff         */
CHAR  *calCacheName     = NULL;         /* location of calibration data */
INT4   globCalData      = 0;            /* glob for calibration frames  */
INT4   pointCal         = 0;            /* don't average cal over chunk */
REAL4  dynRangeExponent = 0;            /* exponent of dynamic range    */
REAL4  strainHighPassFreq = -1;         /* h(t) high pass frequency     */
INT4   strainHighPassOrder = -1;        /* h(t) high pass filter order  */
REAL4  strainHighPassAtten = -1;        /* h(t) high pass attenuation   */
REAL8 (*specFunc)(REAL8) = NULL;        /* pointer to simPSD functions  */

/* template bank generation parameters */
REAL4   minMass         = -1;           /* minimum component mass       */
REAL4   maxMass         = -1;           /* maximum component mass       */
REAL4   minTotalMass    = -1;           /* minimum total mass           */
REAL4   maxTotalMass    = -1;           /* maximum total mass           */
REAL4   chirpMassCutoff = -1;           /* maximum chirp mass to keep   */
REAL4   etaMinCutoff    = -1;           /* minimum eta to keep          */
REAL4   etaMaxCutoff    = -1;           /* maximum eta to keep          */
REAL4   psi0Min         = 0;            /* minimum value of psi0        */
REAL4   psi0Max         = 0;            /* maximum value of psi0        */
REAL4   psi3Min         = 0;            /* minimum value of psi3        */
REAL4   psi3Max         = 0;            /* maximum value of psi3        */
REAL4   alpha           = 0;            /* BCV amplitude correction     */
REAL4   betaMin         = 0;            /* minimum BCV spin parameter   */
REAL4   betaMax         = 0;            /* maximum BCV spin parameter   */
INT4    maxFcutTmplts   = -1;           /* num tmplts in fcut direction */
REAL4   minMatch        = -1;           /* minimum requested match      */
REAL4   fUpper          = -1;           /* upper frequency cutoff       */
REAL4   chiMin          = 0.0;          /* minimum value of chi for PTF */
REAL4   chiMax          = 1.0;          /* maximum value of chi for PTF */
REAL4   kappaMin        = -1.0;         /* minimum value of kappa for PTF */
REAL4   kappaMax        = 1.0;          /* maximum value of kappa for PTF */
INT4    nPointsChi      = 3;            /* PTF template bank density    */
INT4    nPointsKappa    = 5;            /* PTF templated bank density   */
LALPNOrder order;                       /* post-Newtonian order         */
Approximant approximant;                /* approximation method         */
CoordinateSpace space;                  /* coordinate space used        */
INT4    haveGridSpacing = 0;            /* flag to indicate gridspacing */
INT4    computeMoments  = 1;
FreqCut maxFreqCut;                     /* Max. upper frequency cutoff  */
FreqCut minFreqCut;                     /* Min. upper frequency cutoff  */
INT4    numFreqCut      = 0;            /* # of upper freq. cuts to use */

GridSpacing gridSpacing = SquareNotOriented; /* grid spacing (square or hexa)*/
int     polygonFit      = 1;            /* fit a polygon around BCV bank */

/* standard candle parameters */
INT4    computeCandle = 0;              /* should we compute a candle?  */
REAL4   candleSnr     = -1;             /* candle signal to noise ratio */
REAL4   candleMinMass = -1;             /* standard candle mass (solar) */
REAL4   candleMaxMass = 50;             /* standard candle mass (solar) */

/* TD follow up filenames */
CHAR **tdFileNames = NULL;
INT4  numTDFiles = 0;

/* output parameters */
CHAR  *userTag          = NULL;
CHAR  *ifoTag           = NULL;
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeResponse    = 0;            /* write response function used */
int    writeSpectrum    = 0;            /* write computed psd to file   */
int    writeStrainSpec  = 0;            /* write computed strain spec   */
INT4   outCompress      = 0;

/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main ( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             XLAL_INIT_DECL(status);

  /* frame input data */
  LALCache     *frInCache = NULL;
  LALCache     *frGlobCache = NULL;
  LALFrStream     *frStream = NULL;
  FrChanIn      frChan;

  /* frame output data */
  struct FrFile *frOutFile = NULL;
  struct FrameH *outFrame  = NULL;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL8TimeSeries               strainChan;
  REAL4FrequencySeries          spec;
  COMPLEX8FrequencySeries       resp;

  /* structures for preconditioning */
  ResampleTSParams              resampleParams;

  /* templates */
  InspiralCoarseBankIn          bankIn;
  SnglInspiralTable            *tmplt  = NULL;
  INT4                          numCoarse = 0;

  /* output data */
  MetadataTable         templateBank;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;


  /* counters and other variables */
  UINT4 cut, i, j, k;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  CHAR  fname[256];
  REAL8 respRe, respIm;
  REAL8 shf;
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;
  REAL8 tsLength;
  REAL8 dynRange = 0;

  /* TD follow-up variables */
  int numTdFollow;           /* Number of events to follow up in this time */
  SnglInspiralTable *tdFollowUp   = NULL;
  SnglInspiralTable *thisTdFollow = NULL;

  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *)
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, lalAppsVCSIdentInfo.vcsId,
      lalAppsVCSIdentInfo.vcsStatus, lalAppsVCSIdentInfo.vcsDate, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* can use LALMalloc() / LALCalloc() from here */

  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        " " );
  }
  else
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* make sure the pointer to the first template is null */
  templateBank.snglInspiralTable = NULL;

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;


  /* If we're doing a td follow-up we dont want to generate a bank
   * if there was no BCV trigger in this time
   */
  if ( tdFileNames )
  {
    if ( vrbflg )
    {
      fprintf( stdout, "We are doing a TD follow-up\n" );
      fprintf( stdout, "Following up %d files...\n", numTDFiles );
    }

    numTdFollow = 0;

    for (i = 0; i < (UINT4)numTDFiles; i++ )
    {
      INT4 thisTDNum = 0;
      if ( !tdFollowUp )
      {
        thisTDNum = LALSnglInspiralTableFromLIGOLw(&tdFollowUp,
          tdFileNames[i], 0, -1);
        thisTdFollow = tdFollowUp;
      }
      else
      {
        thisTDNum = LALSnglInspiralTableFromLIGOLw(&(thisTdFollow->next),
          tdFileNames[i], 0, -1);
      }
      if ( thisTDNum < 0 )
      {
        fprintf( stderr, "Error reading file %s\n", tdFileNames[i] );
        exit( 1 );
      }
      numTdFollow += thisTDNum;

      while ( thisTdFollow->next )
      {
        thisTdFollow = thisTdFollow->next;
      }
    }

    if (numTdFollow <= 0) goto cleanExit;

    tdFollowUp  = XLALIfoCutSingleInspiral( &tdFollowUp, ifo );
    if ( tdFollowUp )
       tdFollowUp = XLALTimeCutSingleInspiral( tdFollowUp, &gpsStartTime,
           &gpsEndTime );

    /* If there are no events to follow up, we just exit */
    if ( !tdFollowUp ) goto cleanExit;

    /* Free the follow-up events */
    while (tdFollowUp)
    {
      thisTdFollow = tdFollowUp;
      tdFollowUp = tdFollowUp->next;
      XLALFreeSnglInspiral(&thisTdFollow);
    }
  }


  if ( dynRangeExponent )
  {
    /* compute the dynamic range scaling for the psd computation */
    dynRange = (REAL8) pow( 2.0, dynRangeExponent );
  }
  else
  {
    dynRange = 1.0;
  }
  if ( vrbflg )
    fprintf( stdout, "using dynamic range scaling %e\n", dynRange );



  if ( (specType==specType_mean) || (specType==specType_median) )
  {

    /*
     *
     * read in the input data channel
     *
     */

    /* set the time series parameters of the input data and resample params */
    memset( &resampleParams, 0, sizeof(ResampleTSParams) );
    resampleParams.deltaT = 1.0 / (REAL8) sampleRate;

    /* set the params of the input data time series */
    memset( &chan, 0, sizeof(REAL4TimeSeries) );
    memset( &strainChan, 0, sizeof(REAL8TimeSeries) );
    chan.epoch = gpsStartTime;
    chan.epoch.gpsSeconds -= padData; /* subtract pad seconds from start */

    /* copy the start time into the REAL8 h(t) time series */
    strainChan.epoch = chan.epoch;

    if ( globFrameData )
    {
      CHAR ifoRegExPattern[6];

      if ( vrbflg ) fprintf( stdout, "globbing for *.gwf frame files from %c "
          "of type %s in current directory\n", fqChanName[0], frInType );

      frGlobCache = NULL;

      /* create a frame cache by globbing all *.gwf files in the pwd */
      frGlobCache = XLALCacheGlob(NULL, NULL);

      /* check we globbed at least one frame file */
      if ( ! frGlobCache->length )
      {
        fprintf( stderr, "error: no frame files of type %s found\n",
            frInType );
        exit( 1 );
      }

      /* sieve out the requested data type */
      snprintf( ifoRegExPattern,
          XLAL_NUM_ELEM(ifoRegExPattern), ".*%c.*",
          fqChanName[0] );
      frInCache = XLALCacheDuplicate(frGlobCache);
      XLALCacheSieve(frInCache, 0, 0, ifoRegExPattern, frInType, NULL);

      /* check we got at least one frame file back after the sieve */
      if ( ! frInCache->length )
      {
        fprintf( stderr, "error: no frame files of type %s globbed as input\n",
            frInType );
        exit( 1 );
      }

      XLALDestroyCache( frGlobCache );
    }
    else
    {
      if ( vrbflg ) fprintf( stdout,
          "reading frame file locations from cache file: %s\n", frInCacheName );

      /* read a frame cache from the specified file */
      frInCache = XLALCacheImport(frInCacheName);
    }

    /* open the input data frame stream from the frame cache */
    LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );

    /* set the mode of the frame stream to fail on gaps or time errors */
    frStream->mode = LAL_FR_STREAM_VERBOSE_MODE;

    /* enable frame-file checksum checking */
    XLALFrStreamSetMode( frStream, frStream->mode | LAL_FR_STREAM_CHECKSUM_MODE );

    /* seek to required epoch and set chan name */
    LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
    frChan.name = fqChanName;

    if ( calData == real_8 )
    {
      /* determine the sample rate of the raw data */
      LAL_CALL( LALFrGetREAL8TimeSeries( &status, &strainChan, &frChan,
          frStream ), &status );

      /* copy the data parameters from the h(t) channel to input data channel */
      snprintf( chan.name, LALNameLength, "%s", strainChan.name );
      chan.epoch          = strainChan.epoch;
      chan.deltaT         = strainChan.deltaT;
      chan.f0             = strainChan.f0;
      chan.sampleUnits    = strainChan.sampleUnits;
    }
    else
    {
      /* determine the sample rate of the raw data and allocate enough memory */
      LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
          &status );
    }

    /* determine if we need to resample the channel */
    if ( vrbflg )
    {
      fprintf( stdout, "resampleParams.deltaT = %e\n", resampleParams.deltaT );
      fprintf( stdout, "chan.deltaT = %e\n", chan.deltaT );
    }
    if ( ! ( fabs( resampleParams.deltaT - chan.deltaT ) < epsilon ) )
    {
      resampleChan = 1;
      if ( vrbflg )
        fprintf( stdout, "input channel will be resampled\n" );

      if ( resampFiltType == 0 )
      {
        resampleParams.filterType = LDASfirLP;
      }
      else if ( resampFiltType == 1 )
      {
        resampleParams.filterType = defaultButterworth;
      }
    }

    /* determine the number of points to get and create storage for the data */
    inputLengthNS = (REAL8) ( LAL_INT8_C(1000000000) *
        ( gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds + 2 * padData ) );
    chan.deltaT *= 1.0e9;
    numInputPoints = (UINT4) floor( inputLengthNS / chan.deltaT + 0.5 );
    if ( calData == real_8 )
    {
      /* create storage for the REAL8 h(t) input data */
      LAL_CALL( LALDCreateVector( &status, &(strainChan.data), numInputPoints ),
          &status );
    }
    LAL_CALL( LALSCreateVector( &status, &(chan.data), numInputPoints ),
        &status );

    if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
        "(deltaT) = %e\nreading %d points from frame stream\n", fqChanName,
        chan.deltaT / 1.0e9, numInputPoints );

    if ( calData == real_8 )
    {
      /* read in the REAL8 h(t) data here */
      PassBandParamStruc strainHighpassParam;

      /* read the REAL8 h(t) data from the time series into strainChan      */
      /* which already has the correct amount of memory allocated */
      if ( vrbflg ) fprintf( stdout, "reading REAL8 h(t) data from frames... " );

      LAL_CALL( LALFrGetREAL8TimeSeries( &status, &strainChan, &frChan,
          frStream ), &status);

      if ( vrbflg ) fprintf( stdout, "done\n" );

      /* high pass the h(t) data using the parameters specified on the cmd line*/
      strainHighpassParam.nMax = strainHighPassOrder;
      strainHighpassParam.f1 = -1.0;
      strainHighpassParam.f2 = (REAL8) strainHighPassFreq;
      strainHighpassParam.a1 = -1.0;
      strainHighpassParam.a2 = (REAL8)(1.0 - strainHighPassAtten);
      if ( vrbflg ) fprintf( stdout,
          "applying %d order high pass to REAL8 h(t) data: "
          "%3.2f of signal passes at %4.2f Hz\n",
          strainHighpassParam.nMax, strainHighpassParam.a2,
          strainHighpassParam.f2 );

      LAL_CALL( LALButterworthREAL8TimeSeries( &status, &strainChan,
            &strainHighpassParam ), &status );

      /* cast the REAL8 h(t) data to REAL4 in the chan time series       */
      /* which already has the correct amount of memory allocated */
      for ( j = 0 ; j < numInputPoints ; ++j )
      {
        chan.data->data[j] = (REAL4) ( strainChan.data->data[j] * dynRange );
      }

      /* re-copy the data parameters from h(t) channel to input data channel */
      snprintf( chan.name, LALNameLength, "%s", strainChan.name );
      chan.epoch          = strainChan.epoch;
      chan.deltaT         = strainChan.deltaT;
      chan.f0             = strainChan.f0;
      chan.sampleUnits    = strainChan.sampleUnits;

      /* free the REAL8 h(t) input data */
      LAL_CALL( LALDDestroyVector( &status, &(strainChan.data) ), &status );
      strainChan.data = NULL;
    }
    else
    {
      /* read the data channel time series from frames */
      LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
          &status );

      if ( calData == real_4 )
      {
        /* multiply the input data by dynRange */
        for ( j = 0 ; j < numInputPoints ; ++j )
        {
          chan.data->data[j] *= dynRange;
        }
      }
    }
    memcpy( &(chan.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );

    /* store the start and end time of the raw channel in the search summary */
    /* FIXME:  loss of precision;  consider
    searchsumm.searchSummaryTable->in_start_time = searchsumm.searchSummaryTable->in_end_time = chan.epoch;
    XLALGPSAdd(&searchsumm.searchSummaryTable->in_end_time, chan.deltaT * (REAL8) chan.data->length);
    */
    searchsumm.searchSummaryTable->in_start_time = chan.epoch;
    tsLength = XLALGPSGetREAL8(&(chan.epoch) );
    tsLength += chan.deltaT * (REAL8) chan.data->length;
    XLALGPSSetREAL8( &(searchsumm.searchSummaryTable->in_end_time), tsLength );

    /* close the frame file stream and destroy the cache */
    LAL_CALL( LALFrClose( &status, &frStream ), &status );
    XLALDestroyCache( frInCache );

    /* write the raw channel data as read in from the frame files */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
        &chan, "ct", "RAW" );

    if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
        "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n",
        chan.name, chan.data->length, chan.deltaT,
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

    /* resample the input data */
    if ( resampleChan )
    {
      if (vrbflg) fprintf( stdout, "resampling input data from %e to %e\n",
          chan.deltaT, resampleParams.deltaT );

      LAL_CALL( LALResampleREAL4TimeSeries( &status, &chan, &resampleParams ),
          &status );

      if ( vrbflg ) fprintf( stdout, "channel %s resampled:\n"
          "%d points with deltaT %e\nstarting at GPS time %d sec %d ns\n",
          chan.name, chan.data->length, chan.deltaT,
          chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

      /* write the resampled channel data as read in from the frame files */
      if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
          &chan, "ct", "RAW_RESAMP" );
    }
  }

  /*
   *
   * compute a calibrated strain spectrum
   *
   */

  /* create storage for the response and spectrum */
  memset( &spec, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( &status, &(spec.data), numPoints / 2 + 1 ),
      &status );
  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ),
      &status );
  resp.epoch = spec.epoch = gpsStartTime;

  if ( (specType==specType_mean) || (specType==specType_median) )
  {
    /* iir filter to remove low frequencies from data channel */
    if ( highPass )
    {
      PassBandParamStruc highpassParam;
      highpassParam.nMax = highPassOrder;
      highpassParam.f1 = -1.0;
      highpassParam.f2 = (REAL8) highPassFreq;
      highpassParam.a1 = -1.0;
      highpassParam.a2 = (REAL8)(1.0 - highPassAtten); /* a2 is not attenuation */

      if ( vrbflg ) fprintf( stdout, "applying %d order high pass: "
          "%3.2f of signal passes at %4.2f Hz\n",
          highpassParam.nMax, highpassParam.a2, highpassParam.f2 );

      LAL_CALL( LALDButterworthREAL4TimeSeries( &status, &chan, &highpassParam ),
          &status );
    }

    /* remove pad from requested data from start and end of time series */
    memmove( chan.data->data, chan.data->data + padData * sampleRate,
        (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
    XLALRealloc( chan.data->data,
        (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
    chan.data->length -= 2 * padData * sampleRate;
    chan.epoch.gpsSeconds += padData;

    if ( vrbflg ) fprintf( stdout, "after removal of %d second padding at "
        "start and end:\ndata channel sample interval (deltaT) = %e\n"
        "data channel length = %d\nstarting at %d sec %d ns\n",
        padData , chan.deltaT , chan.data->length,
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

    /* store the start and end time of the filter channel in the search summ */
    /* FIXME:  loss of precision;  consider
    searchsumm.searchSummaryTable->out_start_time = chan.epoch;
    XLALGPSAdd(&searchsumm.searchSummaryTable->out_start_time, chan.deltaT * (REAL8) chan.data->length);
    */
    searchsumm.searchSummaryTable->out_start_time = chan.epoch;
    tsLength = XLALGPSGetREAL8( &(chan.epoch) );
    tsLength += chan.deltaT * (REAL8) chan.data->length;
    XLALGPSSetREAL8( &(searchsumm.searchSummaryTable->out_end_time), tsLength );

    /* compute the windowed power spectrum for the data channel */
    {	/* Kipp:  ported this block from LALStatus 20190806 */
    REAL4Window *window = XLALCreateHannREAL4Window( numPoints );
    REAL4FFTPlan *plan = XLALCreateForwardREAL4FFTPlan( numPoints, 0 );
    switch ( specType )
    {
      case specType_mean:
        if ( vrbflg ) fprintf( stdout, "computing mean psd with overlap %d\n", numPoints / 2 );
        XLALREAL4AverageSpectrumWelch( &spec, &chan, numPoints, numPoints / 2, window, plan );
        break;
      case specType_median:
        if ( vrbflg ) fprintf( stdout, "computing median psd with overlap %d\n", numPoints / 2 );
        XLALREAL4AverageSpectrumMedian( &spec, &chan, numPoints, numPoints / 2, window, plan );
        break;
      case specType_simulated:
        fprintf( stderr, "This should never happen! Exiting..." );
        exit( 1 );
      default:
        fprintf( stderr, "unknown spectrum type %d\n", specType );
        exit( 1 );
    }

    XLALDestroyREAL4Window( window );
    XLALDestroyREAL4FFTPlan( plan );
    }
    LAL_CALL( LALSDestroyVector( &status, &(chan.data) ), &status );
  }

  if ( specType == specType_simulated )
  /* initialize spectrum frequency bins */
  {
    spec.f0 = 0.0;
    /* frequency bin (Hz) is 1/duration of segment */
    /* = sample rate/number of points */
    spec.deltaF = ( (REAL8) sampleRate ) / ( (REAL8) numPoints );

    /* evaluate the simulated PSD */
    for ( k = 0; k < spec.data->length; ++k )
    {
      REAL8 sim_psd_freq = (REAL8) k * spec.deltaF;

      /* PSD can be a very small number, rescale it before casting to REAL4 */
      spec.data->data[k] = (REAL4) (dynRange * dynRange * specFunc( sim_psd_freq ));
    }
  }

  /* write the spectrum data to a file */
  if ( writeSpectrum )
  {
    strcpy( spec.name, chan.name );
    outFrame = fr_add_proc_REAL4FrequencySeries( outFrame,
        &spec, "ct/sqrtHz", "PSD" );
  }

  /* set the parameters of the response to match the data and spectrum */
  resp.deltaF = spec.deltaF;
  resp.f0 = spec.f0;
  resp.sampleUnits = strainPerCount;

  if ( calData || (specType==specType_simulated) )
  {
    /* if we are using calibrated data or a design PSD set the response to unity */
    /* i.e. 1 over the dynamic range scaling factor */
    if ( vrbflg ) fprintf( stdout, "generating unity response function\n" );
    for( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k] = (REAL4) (1.0 / dynRange);
    }
  }
  else
  {
    fprintf( stderr, "uncalibrated data no longer supported" );
    exit( 1 );
  }

  /* write the calibration data to a file */
  if ( writeResponse )
  {
    strcpy( resp.name, chan.name );
    outFrame = fr_add_proc_COMPLEX8FrequencySeries( outFrame,
        &resp, "strain/ct", "RESPONSE" );
  }

  /* set low frequency cutoff of power spectrum */
  cut = fLow / spec.deltaF > 1 ?  fLow / spec.deltaF : 1;

  /* compute a calibrated strain power spectrum */
  bankIn.shf.epoch = spec.epoch;
  memcpy( bankIn.shf.name, spec.name, LALNameLength * sizeof(CHAR) );
  bankIn.shf.deltaF = spec.deltaF;
  bankIn.shf.f0 = spec.f0;
  bankIn.shf.data = NULL;
  if (XLALUnitMultiply( &(bankIn.shf.sampleUnits), &(spec.sampleUnits), &(resp.sampleUnits) ) == NULL) {
    return LAL_EXLAL;
  }
  LAL_CALL( LALDCreateVector( &status, &(bankIn.shf.data), spec.data->length ),
      &status );
  memset( bankIn.shf.data->data, 0,
      bankIn.shf.data->length * sizeof(REAL8) );

  shf = spec.data->data[cut] *
    ( crealf(resp.data->data[cut]) * crealf(resp.data->data[cut]) +
      cimagf(resp.data->data[cut]) * cimagf(resp.data->data[cut]) );
  for ( k = 1; k < cut ; ++k )
  {
    bankIn.shf.data->data[k] = shf;
  }
  for ( k = cut; k < bankIn.shf.data->length; ++k )
  {
    respRe = (REAL8) crealf(resp.data->data[k]);
    respIm = (REAL8) cimagf(resp.data->data[k]);
    bankIn.shf.data->data[k] = (REAL8) spec.data->data[k] *
      ( respRe * respRe + respIm * respIm );
  }

  /* write the scaled strain spectrum data to a file */
  if ( writeStrainSpec )
  {
#if 0
    strcpy( spec.name, chan.name );
    outFrame = fr_add_proc_REAL8FrequencySeries( outFrame,
        &(bankIn.shf), "strain/sqrtHz", "STRAIN_PSD" );
#endif
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.strainspec.txt",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
    LALDPrintFrequencySeries( &(bankIn.shf), fname );
  }


  /*
   *
   * compute the standard candle distance
   *
   */

  if ( computeCandle )
  {
    while ( candleMinMass <= candleMaxMass )
    {
      if ( vrbflg ) fprintf( stdout, "maximum distance for (%3.2f,%3.2f) "
          "at signal-to-noise %3.2f = ", candleMinMass, candleMinMass, candleSnr );

      candleMinMass = candleMinMass + 1.0;
    }
  }


  /*
   *
   * compute the template bank
   *
   */


  /* bank generation parameters */
  bankIn.mMin          = (REAL8) minMass;
  bankIn.mMax          = (REAL8) maxMass;
  if (maxTotalMass > 0)
  {
    bankIn.MMax        = (REAL8) maxTotalMass;
    bankIn.mMax        = bankIn.MMax - bankIn.mMin;
    if (minTotalMass > 0)
    {
      bankIn.massRange = MinMaxComponentTotalMass;
      bankIn.MMin      = (REAL8) minTotalMass;
    }
    else
    {
      bankIn.massRange   = MinComponentMassMaxTotalMass;
    }
  }
  else
  {
    bankIn.massRange     = MinMaxComponentMass;
    bankIn.MMax          = bankIn.mMax * 2.0;
  }
  bankIn.psi0Min          = (REAL8) psi0Min;
  bankIn.psi0Max          = (REAL8) psi0Max;
  bankIn.psi3Min          = (REAL8) psi3Min;
  bankIn.psi3Max          = (REAL8) psi3Max;
  bankIn.numFcutTemplates = (UINT4) maxFcutTmplts;
  bankIn.alpha            = (REAL8) alpha;
  bankIn.betaMin          = (REAL8) betaMin;
  bankIn.betaMax          = (REAL8) betaMax;
  bankIn.chiMin           = (REAL8) chiMin;
  bankIn.chiMax           = (REAL8) chiMax;
  bankIn.kappaMin         = (REAL8) kappaMin;
  bankIn.kappaMax         = (REAL8) kappaMax;
  bankIn.nPointsChi       = nPointsChi;
  bankIn.nPointsKappa     = nPointsKappa;
  bankIn.mmCoarse         = (REAL8) minMatch;
  bankIn.mmFine           = 0.99; /* doesn't matter since no fine bank yet */
  bankIn.fLower           = (REAL8) fLow;
  bankIn.fUpper           = (REAL8) fUpper;
  bankIn.iflso            = 0; /* currently not implemented */
  bankIn.tSampling        = (REAL8) sampleRate;
  bankIn.order            = order;
  bankIn.approximant      = approximant;
  bankIn.gridSpacing      = gridSpacing;
  bankIn.space            = space;
  bankIn.insidePolygon    = polygonFit; /*by default it uses a polygon fitting. */
  bankIn.etamin           = bankIn.mMin * ( bankIn.MMax - bankIn.mMin) /
    ( bankIn.MMax * bankIn.MMax );
  bankIn.LowGM            = -4.;
  bankIn.HighGM           = 6.;
  bankIn.computeMoments   = computeMoments; /* by default, gammas/moments are recomputed */
  bankIn.maxFreqCut       = maxFreqCut;
  bankIn.minFreqCut       = minFreqCut;
  bankIn.numFreqCut       = numFreqCut;

  /* generate the template bank */
  if ( vrbflg )
  {
    fprintf( stdout, "generating template bank parameters... " );
    fflush( stdout );
  }
  LAL_CALL( LALInspiralBankGeneration( &status, &bankIn, &tmplt, &numCoarse),
      &status );

  /* Do chirp mass cut */
  if ( chirpMassCutoff > 0 )
  {
    tmplt = XLALMassCut( tmplt, "mchirp", 0, chirpMassCutoff, -1, -1 );
    /* count the remaining tmplts */
    numCoarse = XLALCountSnglInspiral( tmplt );
  }

  /* Do eta cut */
  if ( etaMinCutoff >= 0 || etaMaxCutoff > 0 )
  {
    tmplt = XLALMassCut( tmplt, "eta", etaMinCutoff, etaMaxCutoff, -1, -1 );
    /* count the remaining tmplts */
    numCoarse = XLALCountSnglInspiral( tmplt );
  }

  if ( vrbflg )
  {
    fprintf( stdout, "done. Got %d templates\n", numCoarse );
    fflush( stdout );
  }

  if ( numCoarse )
  {
    templateBank.snglInspiralTable = tmplt;
    snprintf( tmplt->ifo, LIGOMETA_IFO_MAX, "%s", ifo );
    snprintf( tmplt->search, LIGOMETA_SEARCH_MAX, "tmpltbank" );
    snprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX,
        "%s", channelName );
    while( (tmplt = tmplt->next) )
    {
      snprintf( tmplt->ifo, LIGOMETA_IFO_MAX, "%s", ifo );
      snprintf( tmplt->search, LIGOMETA_SEARCH_MAX, "tmpltbank" );
      snprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX,
          "%s", channelName );
    }
  }

  /* save the number of templates in the search summary table */
  searchsumm.searchSummaryTable->nevents = numCoarse;


  /*
   *
   * free the data storage
   *
   */


  LAL_CALL( LALDDestroyVector( &status, &(bankIn.shf.data) ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(spec.data) ), &status );
  LAL_CALL( LALCDestroyVector( &status, &(resp.data) ), &status );


  /*
   *
   * write the results to disk
   *
   */


  /* write the output frame */
  if ( writeRawData || writeResponse || writeSpectrum )
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.gwf",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
    frOutFile = FrFileONew( fname, 0 );
    FrameWrite( outFrame, frOutFile );
    FrFileOEnd( frOutFile );
  }

cleanExit:
  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if ( userTag && ifoTag && !outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s_%s-%d-%d.xml",
        ifo, ifoTag, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if (userTag && !ifoTag && !outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s-%d-%d.xml",
        ifo, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if (!userTag && ifoTag && !outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s-%d-%d.xml",
        ifo, ifoTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( userTag && ifoTag && outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s_%s-%d-%d.xml.gz",
        ifo, ifoTag, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if (userTag && !ifoTag && outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s-%d-%d.xml.gz",
        ifo, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if (!userTag && ifoTag && outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s-%d-%d.xml.gz",
        ifo, ifoTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if (!userTag && !ifoTag && outCompress)
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.xml.gz",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    snprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.xml",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname ), &status );

  /* write the process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFO_MAX, "%s", ifo );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable,
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
  }

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams,
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* write the search summary table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results,
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm,
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  /* write the template bank to the file */
  if ( templateBank.snglInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, templateBank,
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while ( templateBank.snglInspiralTable )
  {
    tmplt = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( tmplt );
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( tdFileNames ) free( tdFileNames );
  if ( calCacheName ) free( calCacheName );
  if ( frInCacheName ) free( frInCacheName );
  if ( frInType ) free( frInType );
  if ( channelName ) free( channelName );
  if ( fqChanName ) free( fqChanName );
  LALCheckMemoryLeaks();

#ifdef LALAPPS_CUDA_ENABLED
  cudaDeviceReset();
#endif

  exit( 0 );
}

/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) );\
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
   PROGRAM_NAME );\
   snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
     long_options[option_index].name );\
     snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype );\
     snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define USAGE( a ) \
fprintf(a, "  --help                       display this message\n");\
fprintf(a, "  --verbose                    print progress information\n");\
fprintf(a, "  --version                    print version information and exit\n");\
fprintf(a, "  --user-tag STRING            set the process_params usertag to STRING\n");\
fprintf(a, "  --ifo-tag STRING             set the ifotag to STRING - for file naming\n");\
fprintf(a, "  --comment STRING             set the process table comment to STRING\n");\
fprintf(a, "  --write-compress             write a compressed xml file\n");\
fprintf(a, "\n");\
fprintf(a, "  NOTE: Data-related options not required when using a simulated PSD are labelled @\n");\
fprintf(a, "\n");\
fprintf(a, "  --gps-start-time SEC         GPS second of data start time\n");\
fprintf(a, "  --gps-end-time SEC           GPS second of data end time\n");\
fprintf(a, "@ --pad-data T                 pad the data start and end time by T seconds\n");\
fprintf(a, "\n");\
fprintf(a, "@ --glob-frame-data            glob *.gwf files in the pwd to obtain frame data\n");\
fprintf(a, "@ --frame-type TAG             input data is contained in frames of type TAG\n");\
fprintf(a, "@ --frame-cache                obtain frame data from LAL frame cache FILE\n");\
fprintf(a, "@ --calibration-cache FILE     obtain calibration from LAL frame cache FILE\n");\
fprintf(a, "@ --glob-calibration-data      obtain calibration by globbing in working dir\n");\
fprintf(a, "\n");\
fprintf(a, "@ --channel-name CHAN          read data from interferometer channel CHAN\n");\
fprintf(a, "@ --calibrated-data TYPE       calibrated data of TYPE real_4 or real_8\n");\
fprintf(a, "@ --strain-high-pass-freq F    high pass REAL8 h(t) data above F Hz\n");\
fprintf(a, "@ --strain-high-pass-order O   set the order of the h(t) high pass filter to O\n");\
fprintf(a, "@ --strain-high-pass-atten A   set the attenuation of the high pass filter to A\n");\
fprintf(a, "@ --point-calibration          use the first point in the chunk to calibrate\n");\
fprintf(a, "\n");\
fprintf(a, "  --sample-rate F              filter data at F Hz, downsampling if necessary\n");\
fprintf(a, "@ --resample-filter TYPE       set resample filter to TYPE [ldas|butterworth]\n");\
fprintf(a, "\n");\
fprintf(a, "  --disable-high-pass          turn off the IIR highpass filter\n");\
fprintf(a, "@ --enable-high-pass F         high pass data above F Hz using an IIR filter\n");\
fprintf(a, "@ --high-pass-order O          set the order of the high pass filter to O\n");\
fprintf(a, "@ --high-pass-attenuation A    set the attenuation of the high pass filter to A\n");\
fprintf(a, "  --spectrum-type TYPE         use PSD estimator TYPE \n");\
fprintf(a, "                               (mean|median|iLIGOSRD|eLIGOModel|GEOModel|\n");\
fprintf(a, "                               |aLIGONoSRMLoP|aLIGONoSRMHiP|aLIGOZDLoP|aLIGOZDHiP|\n");\
fprintf(a, "                               |iVirgoModel|aVirgoModel|KAGRAModel)\n");\
fprintf(a, "  --dynamic-range-exponent X   set dynamic range scaling to 2^X (eg X=69.0)\n");\
fprintf(a, "\n");\
fprintf(a, "  --segment-length N           set data segment length to N points\n");\
fprintf(a, "@ --number-of-segments N       set number of data segments to N\n");\
fprintf(a, "\n");\
fprintf(a, "  --td-follow-up FILE          follow up BCV events contained in FILE\n");\
fprintf(a, "\n");\
fprintf(a, "  --standard-candle            compute a standard candle from the PSD\n");\
fprintf(a, "  --candle-snr SNR             signal-to-noise ratio of standard candle\n");\
fprintf(a, "  --candle-minmass M           minimum component mass for (equal-mass) candle binary\n");\
fprintf(a, "  --candle-maxmass M           maximum component mass for candle binary, default=50\n");\
fprintf(a, "\n");\
fprintf(a, "  --low-frequency-cutoff F     do not filter below F Hz\n");\
fprintf(a, "  --high-frequency-cutoff F    upper frequency cutoff in Hz\n");\
fprintf(a, "  --disable-compute-moments    do not recompute the moments stored in the template bank. \n");\
fprintf(a, "\n");\
fprintf(a, "  --minimum-mass MASS          set minimum component mass of bank to MASS: required\n");\
fprintf(a, "  --maximum-mass MASS          set maximum component mass of bank to MASS\n");\
fprintf(a, "  --max-total-mass MASS        set maximum total mass of the bank to MASS. Will override --maximum-mass option\n");\
fprintf(a, "  --min-total-mass MASS        set minimum total mass of the bank to MASS: --max-total-mass must also be given\n");\
fprintf(a, "  --chirp-mass-cutoff MASS     set chirp mass cutoff to MASS\n");\
fprintf(a, "  --max-eta ETA                set maximum symmetric mass ratio of the bank to ETA\n");\
fprintf(a, "  --min-eta ETA                set minimum symmetric mass ratio of the bank to ETA\n");\
fprintf(a, "\n");\
fprintf(a, "  --minimum-psi0 PSI0          set minimum range of BCV parameter psi0 to PSI0\n");\
fprintf(a, "  --maximum-psi0 PSI0          set maximum range of BCV parameter psi0 to PSI0\n");\
fprintf(a, "  --minimum-psi3 PSI3          set minimum range of BCV parameter psi3 to PSI3\n");\
fprintf(a, "  --maximum-psi3 PSI3          set maximum range of BCV parameter psi3 to PSI3\n");\
fprintf(a, "  --maximum-fcut-tmplts N      maximum number of tmplts in fcut direction is N\n");\
fprintf(a, "  --disable-polygon-fit        disable the polygon fitting for BCV bank\n");\
fprintf(a, "  --alpha ALPHA                set alpha for the BCV bank generation\n");\
fprintf(a, "  --minimum-beta BETA          set minimum BCV spin parameter beta to BETA\n");\
fprintf(a, "  --maximum-beta BETA          set maximum BCV spin parameter beta to BETA\n");\
fprintf(a, "\n");\
fprintf(a, "  --minimum-spin1 SPIN1_MIN    set minimum value of chi for PTF to SPIN1_MIN (0.0)\n");\
fprintf(a, "  --maximum-spin1 SPIN1_MAX    set maximum value of chi for PTF to SPIN1_MAX (1.0)\n");\
fprintf(a, "  --minimum-kappa1 KAPPA1_MIN  set minimum value of kappa for PTF to KAPPA1_MIN (-1.0)\n");\
fprintf(a, "  --maximum-kappa1 KAPPA1_MAX  set maximum value of kappa for PTF to KAPPA1_MAX (1.0)\n");\
fprintf(a, "  --npoints-chi N-CHI          set number of points in the Chi direction for PTF template bank to N-CHI (3)\n");\
fprintf(a, "  --npoints-kappa N-KAPPA      set number of points in the Kappa direction for PTF template bank to N-KAPPA (5)\n");\
fprintf(a, "\n");\
fprintf(a, "  --minimal-match M            generate bank with minimal match M\n");\
fprintf(a, "\n");\
fprintf(a, "  --order ORDER                set post-Newtonian order of the waveform to ORDER\n");\
fprintf(a, "                                 (newtonian|oneHalfPN|onePN|onePointFivePN|\n");\
fprintf(a, "                                 twoPN|twoPointFive|threePN|threePointFivePN)\n");\
fprintf(a, "  --approximant APPROX         set approximant of the waveform to APPROX\n");\
fprintf(a, "                                 (TaylorT1|TaylorT2|TaylorT3|TaylorF1|TaylorF2|\n");\
fprintf(a, "                                 PadeT1|PadeT2|EOB|EOBNR|BCV|SpinTaylorT3|BCVSpin)\n");\
fprintf(a, "  --num-freq-cutoffs Ncut       create a template bank with Ncut different upper \n");\
fprintf(a, "                                 frequency cutoffs (must be a positive integer) \n");\
fprintf(a, "  --max-high-freq-cutoff MAX    formula to compute the largest high freq. cutoff\n");\
fprintf(a, "                                 possible choices in ascending order: (SchwarzISCO|BKLISCO|LightRing|FRD|ERD|LRD)\n");\
fprintf(a, "  --min-high-freq-cutoff MIN    formula to compute the smallest high freq. cutoff\n");\
fprintf(a, "                                 possible choices in ascending order: (SchwarzISCO|BKLISCO|LightRing|FRD|ERD|LRD)\n");\
fprintf(a, "  --space SPACE                grid up template bank with mass parameters SPACE\n");\
fprintf(a, "                                 (Tau0Tau2|Tau0Tau3|Psi0Psi3)\n");\
fprintf(a, "  --grid-spacing GRIDSPACING   grid up template bank with GRIDSPACING\n");\
fprintf(a, "                                 (Hexagonal|SquareNotOriented)\n");\
fprintf(a, "\n");\
fprintf(a, "  --write-response             write the computed response function to a frame\n");\
fprintf(a, "  --write-spectrum             write the uncalibrated psd to a frame\n");\
fprintf(a, "  --write-strain-spectrum      write the calibrated strain psd to a text file\n");


int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* LALgetopt arguments */
  struct LALoption long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"disable-high-pass",       no_argument,       &highPass,         0 },
    {"standard-candle",         no_argument,       &computeCandle,    1 },
    {"glob-frame-data",         no_argument,       &globFrameData,    1 },
    {"glob-calibration-data",   no_argument,       &globCalData,      1 },
    {"point-calibration",       no_argument,       &pointCal,         1 },
    /* parameters used to generate calibrated power spectrum */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"channel-name",            required_argument, 0,                'c'},
    {"segment-length",          required_argument, 0,                'd'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"sample-rate",             required_argument, 0,                'g'},
#ifdef LALAPPS_CUDA_ENABLED
    {"gpu-device-id",           required_argument, 0,                '+'},
#endif
    {"calibrated-data",         required_argument, 0,                'M'},
    {"strain-high-pass-freq",   required_argument, 0,                'J'},
    {"strain-high-pass-order",  required_argument, 0,                'K'},
    {"strain-high-pass-atten",  required_argument, 0,                'L'},
    {"help",                    no_argument,       0,                'h'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"dynamic-range-exponent",  required_argument, 0,                'f'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass",        required_argument, 0,                't'},
    {"high-pass-order",         required_argument, 0,                'H'},
    {"high-pass-attenuation",   required_argument, 0,                'I'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"frame-type",              required_argument, 0,                'n'},
    {"pad-data",                required_argument, 0,                'x'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"ifo-tag",                 required_argument, 0,                'Y'},
    {"version",                 no_argument,       0,                'V'},
    {"resample-filter",         required_argument, 0,                'r'},
    /* template bank generation parameters */
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"minimum-psi0",            required_argument, 0,                'P'},
    {"maximum-psi0",            required_argument, 0,                'Q'},
    {"minimum-psi3",            required_argument, 0,                'R'},
    {"maximum-psi3",            required_argument, 0,                'S'},
    {"maximum-fcut-tmplts",     required_argument, 0,                'U'},
    {"minimum-beta",            required_argument, 0,                'o'},
    {"maximum-beta",            required_argument, 0,                'O'},
    {"alpha",                   required_argument, 0,                'T'},
    {"minimum-spin1",           required_argument, 0,                '4'},
    {"maximum-spin1",           required_argument, 0,                '5'},
    {"minimum-kappa1",          required_argument, 0,                '6'},
    {"maximum-kappa1",          required_argument, 0,                '7'},
    {"npoints-chi",             required_argument, 0,                '8'},
    {"npoints-kappa",           required_argument, 0,                '9'},
    {"minimal-match",           required_argument, 0,                'C'},
    {"high-frequency-cutoff",   required_argument, 0,                'D'},
    {"order",                   required_argument, 0,                'E'},
    {"approximant",             required_argument, 0,                'F'},
    {"num-freq-cutoffs",        required_argument, 0,                '1'},
    {"max-high-freq-cutoff",    required_argument, 0,                '2'},
    {"min-high-freq-cutoff",    required_argument, 0,                '3'},
    {"space",                   required_argument, 0,                'G'},
    {"grid-spacing",            required_argument, 0,                'v'},
    {"max-total-mass",          required_argument, 0,                'y'},
    {"min-total-mass",          required_argument, 0,                'W'},
    {"chirp-mass-cutoff",       required_argument, 0,                'q'},
    {"max-eta",                 required_argument, 0,                '0'},
    {"min-eta",                 required_argument, 0,                'X'},
    {"disable-polygon-fit",     no_argument,            &polygonFit,       0 },
    {"disable-compute-moments", no_argument,            &computeMoments,   0 },
    /* standard candle parameters */
    {"candle-snr",              required_argument, 0,                'k'},
    {"candle-minmass",          required_argument, 0,                'l'},
    {"candle-maxmass",          required_argument, 0,                'm'},
    /* frame writing options */
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-response",          no_argument,       &writeResponse,    1 },
    {"write-spectrum",          no_argument,       &writeSpectrum,    1 },
    {"write-strain-spectrum",   no_argument,       &writeStrainSpec,  1 },
    /* td follow up file */
    {"td-follow-up",            required_argument, 0,                'w'},
    {0, 0, 0, 0}
  };

  int c;
#ifdef LALAPPS_CUDA_ENABLED
  INT4 gpuDeviceID = 0;
  cudaError_t cudaError = cudaSuccess;
#endif
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;
  UINT4   haveOrder       = 0;
  UINT4   haveApprox      = 0;
  UINT4   haveSpace       = 0;
  UINT4   havePsi0Min     = 0;
  UINT4   havePsi0Max     = 0;
  UINT4   havePsi3Min     = 0;
  UINT4   havePsi3Max     = 0;
  UINT4   haveAlpha       = 0;
  UINT4   haveNumFcut     = 0;
  UINT4   haveMaxFcut     = 0;
  UINT4   haveMinFcut     = 0;

  /*
   *
   * parse command line arguments
   *
   */

  while ( 1 )
  {
    /* LALgetopt_long stores long option here */
    int option_index = 0;
    size_t LALoptarg_len;

    c = LALgetopt_long_only( argc, argv,
#ifdef LALAPPS_CUDA_ENABLED
        "a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:r:s:t:u:v:x:yX:0:"
        "A:B:C:D:E:F:G:H:I:J:K:L:M:O:P:Q:R:S:T:U:VZ:1:2:3:4:5:6:7:8:9:+:",
#else
        "a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:r:s:t:u:v:x:yX:0:"
        "A:B:C:D:E:F:G:H:I:J:K:L:M:O:P:Q:R:S:T:U:VZ:1:2:3:4:5:6:7:8:9:",
#endif
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'a':
        {
          long int gstartt = atol( LALoptarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is prior to "
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          gpsStartTime.gpsSeconds = (INT4) gstartt;
          gpsStartTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'b':
        {
          long int gendt = atol( LALoptarg );
          if ( gendt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is prior to "
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          gpsEndTime.gpsSeconds = (INT4) gendt;
          gpsEndTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'c':
        {
          /* create storage for the channel name and copy it */
          char *channamptr = NULL;
          LALoptarg_len = strlen( LALoptarg ) + 1;
          fqChanName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
          memcpy( fqChanName, LALoptarg, LALoptarg_len );
          ADD_PROCESS_PARAM( "string", "%s", LALoptarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName, ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, LALoptarg );
            exit( 1 );
          }
          LALoptarg_len = strlen( ++channamptr ) + 1;
          channelName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
          memcpy( channelName, channamptr, LALoptarg_len );

          /* copy the first two characters to ifo */
          memset( ifo, 0, sizeof(ifo) );
          memcpy( ifo, LALoptarg, sizeof(ifo) - 1 );
        }
        break;

      case 'd':
        numPoints = (INT4) atoi( LALoptarg );
        if ( numPoints < 2 || numPoints % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of points must be a non-zero power of 2: "
              "(%d specified) \n",
              long_options[option_index].name, numPoints );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numPoints );
        break;

      case 'e':
        numSegments = (INT4) atoi( LALoptarg );
        if ( numSegments < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of data segment must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, numSegments );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSegments );
        break;

      case 'g':
        sampleRate = (INT4) atoi( LALoptarg );
        if ( sampleRate < 2 || sampleRate > 16384 || sampleRate % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "rate must be power of 2 between 2 and 16384 inclusive: "
              "(%d specified)\n",
              long_options[option_index].name, sampleRate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", sampleRate );
        break;

#ifdef LALAPPS_CUDA_ENABLED
      case '+':
        gpuDeviceID = (INT4) atoi( LALoptarg );
        cudaError = cudaSetDevice( gpuDeviceID );
        if ( cudaError != cudaSuccess )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
                   "could not associate thread to GPU %d\n"
                   "CudaError: %s\n",
                   long_options[option_index].name, gpuDeviceID,
                   cudaGetErrorString(cudaError));
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", gpuDeviceID );
        break;
#endif

      case 'M':
        /* specify which type of calibrated data */
        {
          if ( ! strcmp( "real_4", LALoptarg ) )
          {
            calData = real_4;
          }
          else if ( ! strcmp( "real_8", LALoptarg ) )
          {
            calData = real_8;
          }
          else
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "unknown data type specified;\n"
                "%s (must be one of: real_4, real_8)\n",
                long_options[option_index].name, LALoptarg);
          }
          ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        }
        break;

      case 'J':
        strainHighPassFreq = (REAL4) atof( LALoptarg );
        if ( strainHighPassFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "REAL8 h(t) high pass filter frequency must be greater than 0 Hz:"
              "(%f Hz specified)\n",
              long_options[option_index].name, strainHighPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", strainHighPassFreq );
        break;

      case 'K':
        strainHighPassOrder = (INT4) atoi( LALoptarg );
        if ( strainHighPassOrder <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "REAL8 h(t) high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, strainHighPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", strainHighPassOrder );
        break;

      case 'L':
        strainHighPassAtten = (REAL4) atof( LALoptarg );
        if ( strainHighPassAtten < 0.0 || strainHighPassAtten > 1.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "REAL8 h(t) high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, strainHighPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", strainHighPassAtten );
        break;

      case 'h':
        USAGE( stdout );
        exit( 0 );
        break;

      case 'i':
        fLow = (REAL4) atof( LALoptarg );
        if ( fLow < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow );
        break;

      case 'j':
        if ( ! strcmp( "mean", LALoptarg ) )
        {
          specType = specType_mean;
        }
        else if ( ! strcmp( "median", LALoptarg ) )
        {
          specType = specType_median;
        }
        else if ( ! strcmp( "iLIGOSRD", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDiLIGOSRD; }
        else if ( ! strcmp( "eLIGOModel", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDeLIGOModel; }
        else if ( ! strcmp( "GEOModel", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDGEO; }
        else if ( ! strcmp( "aLIGONoSRMLoP", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDaLIGONoSRMLowPower; }
        else if ( ! strcmp( "aLIGONoSRMHiP", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDaLIGONoSRMHighPower; }
        else if ( ! strcmp( "aLIGOZDLoP", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDaLIGOZeroDetLowPower; }
        else if ( ! strcmp( "aLIGOZDHiP", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDaLIGOZeroDetHighPower; }
        else if ( ! strcmp( "iVirgoModel", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDVirgo; }
        else if ( ! strcmp( "aVirgoModel", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDAdvVirgo; }
        else if ( ! strcmp( "KAGRAModel", LALoptarg ) )
        { specType = specType_simulated;
          specFunc = XLALSimNoisePSDKAGRA; }

        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown power spectrum type: %s\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'f':
        dynRangeExponent = (REAL4) atof( LALoptarg );
        ADD_PROCESS_PARAM( "float", "%e", dynRangeExponent );
        break;

      case 'p':
        /* create storage for the calibration frame cache name */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        calCacheName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR));
        memcpy( calCacheName, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'r':
        if ( ! strcmp( "ldas", LALoptarg ) )
        {
          resampFiltType = 0;
        }
        else if ( ! strcmp( "butterworth", LALoptarg ) )
        {
          resampFiltType = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown resampling filter type: "
              "%s (must be ldas or butterworth)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 's':
        if ( strlen( LALoptarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", LALoptarg);
        }
        break;

      case 't':
        highPass = 1;
        highPassFreq = (REAL4) atof( LALoptarg );
        if ( highPassFreq < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, highPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassFreq );
        break;

      case 'H':
        highPassOrder = (INT4) atoi( LALoptarg );
        if ( highPassOrder <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, highPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", highPassOrder );
        break;

      case 'I':
        highPassAtten = (REAL4) atof( LALoptarg );
        if ( highPassAtten < 0.0 || highPassAtten > 1.0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, highPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassAtten );
        break;

      case 'u':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        frInCacheName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( frInCacheName, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'n':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        frInType = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( frInType, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'x':
        padData = (UINT4) atoi( LALoptarg );
        if ( padData < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of seconds to pad from input data"
              "must be greater than 0: (%d specified)\n",
              long_options[option_index].name, padData );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", padData );
        break;

      case 'Z':
        /* create storage for the usertag */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        userTag = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( userTag, LALoptarg, LALoptarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--user-tag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );
        break;

      case 'A':
        minMass = (REAL4) atof( LALoptarg );
        if ( minMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMass );
        break;

      case 'B':
        maxMass = (REAL4) atof( LALoptarg );
        if ( maxMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maximum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxMass );
        break;

      case 'P':
        psi0Min = (REAL4) atof( LALoptarg );
        if ( psi0Min <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimum value of psi0 must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, psi0Min );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", psi0Min );
        havePsi0Min = 1;
        break;

      case 'Q':
        psi0Max = (REAL4) atof( LALoptarg );
        if ( psi0Max <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maximum value of psi0 must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, psi0Max );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", psi0Max );
        havePsi0Max = 1;
        break;

      case 'R':
        psi3Min = (REAL4) atof( LALoptarg );
        if ( psi3Min >= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum value of psi3 must be < 0: "
              "(%f specified)\n",
              long_options[option_index].name, psi3Min );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", psi3Min );
        havePsi3Min = 1;
        break;

      case 'S':
        psi3Max = (REAL4) atof( LALoptarg );
        ADD_PROCESS_PARAM( "float", "%e", psi3Max );
        havePsi3Max = 1;
        break;

      case 'o':
        betaMin = (REAL4) atof( LALoptarg );
        ADD_PROCESS_PARAM( "float", "%e", betaMin );
        break;

      case 'O':
        betaMax = (REAL4) atof( LALoptarg );
        ADD_PROCESS_PARAM( "float", "%e", betaMax );
        break;

      case 'U':
        maxFcutTmplts = (INT4) atof( LALoptarg );
        if ( maxFcutTmplts < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "number of templates in f_final direction must be >= 0"
              "(%d specified)\n",
              long_options[option_index].name, maxFcutTmplts );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", maxFcutTmplts );
        break;

      case 'T':
        alpha = (REAL4) atof( LALoptarg );
        if ( alpha < -1 || alpha > 1 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "value of alpha must be the range [0:1]"
              "(%f specified)\n",
              long_options[option_index].name, alpha );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", alpha );
        haveAlpha = 1;
        break;

      case 'C':
        minMatch = (REAL4) atof( LALoptarg );
        if ( minMatch <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimum match of bank must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minMatch );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMatch );
        break;

      case 'D':
        fUpper = (REAL4) atof( LALoptarg );
        if ( fUpper <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumu component mass must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, fUpper );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fUpper );
        break;

      case 'E':
        if ( ! strcmp( "newtonian", LALoptarg ) )
        {
          order = LAL_PNORDER_NEWTONIAN;
        }
        else if ( ! strcmp( "oneHalfPN", LALoptarg ) )
        {
          order = LAL_PNORDER_HALF;
        }
        else if ( ! strcmp( "onePN", LALoptarg ) )
        {
          order = LAL_PNORDER_ONE;
        }
        else if ( ! strcmp( "onePointFivePN", LALoptarg ) )
        {
          order = LAL_PNORDER_ONE_POINT_FIVE;
        }
        else if ( ! strcmp( "twoPN", LALoptarg ) )
        {
          order = LAL_PNORDER_TWO;
        }
        else if ( ! strcmp( "twoPointFive", LALoptarg ) )
        {
          order = LAL_PNORDER_TWO_POINT_FIVE;
        }
        else if ( ! strcmp( "threePN", LALoptarg ) )
        {
          order = LAL_PNORDER_THREE;
        }
        else if ( ! strcmp( "threePointFivePN", LALoptarg ) )
        {
          order = LAL_PNORDER_THREE_POINT_FIVE;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of: newtonian, oneHalfPN, onePN,\n"
              "onePointFivePN, twoPN, twoPointFivePN, threePN or\n"
              "threePointFivePN)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        haveOrder = 1;
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'F':
        if ( ! strcmp( "TaylorT1", LALoptarg ) )
        {
          approximant = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", LALoptarg ) )
        {
          approximant = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", LALoptarg ) )
        {
          approximant = TaylorT3;
        }
        else if ( ! strcmp( "TaylorF1", LALoptarg ) )
        {
          approximant = TaylorF1;
        }
        else if ( ! strcmp( "TaylorF2", LALoptarg ) )
        {
          approximant = TaylorF2;
        }
        else if ( ! strcmp( "PadeT1", LALoptarg ) )
        {
          approximant = PadeT1;
        }
        else if ( ! strcmp( "PadeF1", LALoptarg ) )
        {
          approximant = PadeF1;
        }
        else if ( ! strcmp( "EOB", LALoptarg ) )
        {
          approximant = EOB;
        }
        else if ( ! strcmp( "EOBNR", LALoptarg ) )
        {
          approximant = EOBNR;
        }
        else if ( ! strcmp( "EOBNRv2", LALoptarg ) )
        {
          approximant = EOBNRv2;
        }
        else if ( ! strcmp( "IMRPhenomA", LALoptarg ) )
        {
          approximant = IMRPhenomA;
        }
        else if ( ! strcmp( "IMRPhenomB", LALoptarg ) )
        {
          approximant = IMRPhenomB;
        }
        else if ( ! strcmp( "BCV", LALoptarg ) )
        {
          approximant = BCV;
        }
        else if ( ! strcmp( "SpinTaylorT3", LALoptarg ) )
        {
          approximant = SpinTaylorT3;
        }
        else if ( ! strcmp( "BCVSpin", LALoptarg ) )
        {
          approximant = BCVSpin;
        }
        else if ( ! strcmp( "FindChirpPTF", LALoptarg ) )
        {
          approximant = FindChirpPTF;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of: TaylorT1, TaylorT2, TaylorT3, TaylorF1,\n"
              "TaylorF2, PadeT1, PadeF1, EOB, EOBNR, EOBNRv2, IMRPhenomA,\n"
              "IMRPhenomB, BCV, SpinTaylorT3, BCVSpin\n"
              "or FindChirpPTF)\n", long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        haveApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'G':
        if ( ! strcmp( "Tau0Tau2", LALoptarg ) )
        {
          space = Tau0Tau2;
        }
        else if ( ! strcmp( "Tau0Tau3", LALoptarg ) )
        {
          space = Tau0Tau3;
        }
        else if ( ! strcmp( "Psi0Psi3", LALoptarg ) )
        {
          space = Psi0Psi3;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown space specified: "
              "%s (must be one of: Tau0Tau2, Tau0Tau3 or Psi0Psi3)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        haveSpace = 1;
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'v':
        if ( ! strcmp( "Hexagonal", LALoptarg) )
        {
          haveGridSpacing = 1;
          gridSpacing = Hexagonal;
        }
        else if ( ! strcmp( "SquareNotOriented", LALoptarg) )
        {
          haveGridSpacing = 1;
          gridSpacing = SquareNotOriented;
        }
        else
        {
          fprintf(stderr, "invalid argument to --%s:\n"
              "unknown grid spacing specified: "
              "%s (must be one of  Hexagonal, SquareNotOriented )\n",
              long_options[option_index].name, LALoptarg );
          exit(1);
        }
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'y':
        maxTotalMass = (REAL4) atof( LALoptarg );
        if ( maxTotalMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maximum total mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxTotalMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxTotalMass );
        break;

      case 'W':
        minTotalMass = (REAL4) atof( LALoptarg );
        if ( minTotalMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimum total mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minTotalMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minTotalMass );
        break;

      case 'q':
        chirpMassCutoff = (REAL4) atof( LALoptarg );
        if ( chirpMassCutoff <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "chirp mass cutoff must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, chirpMassCutoff );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", chirpMassCutoff );
        break;

      case 'X':
        etaMinCutoff = (REAL4) atof( LALoptarg );
        if ( etaMinCutoff < 0 || etaMinCutoff >= 0.25 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimum eta must be >= 0 and < 0.25: "
              "(%f specified)\n",
              long_options[option_index].name, etaMinCutoff);
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", etaMinCutoff);
        break;

      case '0':
        etaMaxCutoff = (REAL4) atof( LALoptarg );
        if ( etaMaxCutoff <= 0 || etaMaxCutoff > 0.25 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maximum eta must be > 0 and <= 0.25: "
              "(%f specified)\n",
              long_options[option_index].name, etaMaxCutoff );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", etaMaxCutoff );
        break;

      case 'k':
        candleSnr = (REAL4) atof( LALoptarg );
        if ( candleSnr <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "standard candle signal-to-noise ratio must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, candleSnr );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", candleSnr );
        break;

      case 'l':
        candleMinMass = (REAL4) atof( LALoptarg );
        if ( candleMinMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "standard candle minimum component mass must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, candleMinMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", candleMinMass );
        break;

      case 'm':
        candleMaxMass = (REAL4) atof( LALoptarg );
        if ( ( candleMaxMass <= 0 ) || ( candleMaxMass < candleMinMass ) )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "standard candle maximum component mass must be > 0 and greater than min mass: "
              "(%f specified)\n",
              long_options[option_index].name, candleMaxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", candleMaxMass );
        break;

      case 'w':
        numTDFiles = 1;
        /* Count the number of thinca files to follow up */
        while ( !strstr( argv[LALoptind], "--" ) )
        {
          numTDFiles++;
          LALoptind++;
        }
        LALoptind = LALoptind - numTDFiles;

        /* Set pointers to the relevant filenames */
        tdFileNames = (CHAR **) calloc( numTDFiles, sizeof(CHAR *));
        numTDFiles = 0;

        while ( !strstr( argv[LALoptind], "--" ) )
        {
          tdFileNames[numTDFiles++] = argv[LALoptind];
          ADD_PROCESS_PARAM( "string", "%s", argv[LALoptind] );
          LALoptind++;

        }
        break;

      case 'Y':
        /* create storage for the ifo-tag */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        ifoTag = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( ifoTag, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "LIGO/LSC Standalone Inspiral Template Bank Code\n"
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n");
        XLALOutputVCSInfo(stderr, lalAppsVCSInfoList, 0, "%% ");
        exit( 0 );
        break;

      case '?':
        USAGE( stderr );
        exit( 1 );
        break;

      case '1':
        numFreqCut = (INT4) atof( LALoptarg );
        if( numFreqCut < 1 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Value must be a positive integer "
              "(%d specified)\n",
              long_options[option_index].name, numFreqCut );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numFreqCut );
        haveNumFcut = 1;
        break;

      case '2':
        if ( ! strcmp( "SchwarzISCO", LALoptarg ) )
        {
          maxFreqCut = FreqCut_SchwarzISCO;
        }
        else if( ! strcmp( "BKLISCO", LALoptarg ) )
        {
          maxFreqCut = FreqCut_BKLISCO;
        }
        else if ( ! strcmp( "LightRing", LALoptarg ) )
        {
          maxFreqCut = FreqCut_LightRing;
        }
        else if ( ! strcmp( "FRD", LALoptarg ) )
        {
          maxFreqCut = FreqCut_FRD;
        }
        else if ( ! strcmp( "ERD", LALoptarg ) )
        {
          maxFreqCut = FreqCut_ERD;
        }
        else if ( ! strcmp( "LRD", LALoptarg ) )
        {
          maxFreqCut = FreqCut_LRD;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown cutoff frequency specified: "
              "%s (must be one of: SchwarzISCO, BKLISCO, LightRing, FRD, ERD or LRD)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        haveMaxFcut = 1;
        break;

      case '3':
        if ( ! strcmp( "SchwarzISCO", LALoptarg ) )
        {
          minFreqCut = FreqCut_SchwarzISCO;
        }
        else if ( ! strcmp( "BKLISCO", LALoptarg ) )
        {
          minFreqCut = FreqCut_BKLISCO;
        }
        else if ( ! strcmp( "LightRing", LALoptarg ) )
        {
          minFreqCut = FreqCut_LightRing;
        }
        else if ( ! strcmp( "FRD", LALoptarg ) )
        {
          minFreqCut = FreqCut_FRD;
        }
        else if ( ! strcmp( "ERD", LALoptarg ) )
        {
          minFreqCut = FreqCut_ERD;
        }
        else if ( ! strcmp( "LRD", LALoptarg ) )
        {
          minFreqCut = FreqCut_LRD;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown cutoff frequency specified: "
              "%s (must be one of: SchwarzISCO, BKLISCO, LightRing, FRD, ERD, or LRD)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        haveMinFcut = 1;
        break;

      case '4':
        chiMin = atof( LALoptarg );
        if ( chiMin < 0. )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Spin magnitude can only take values between 0 and 1. : "
              "(%f specified)\n",
              long_options[option_index].name, chiMin );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", chiMin );
        break;

      case '5':
        chiMax = atof( LALoptarg );
        if ( chiMax > 1. )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Spin magnitude can only take values between 0 and 1. : "
              "(%f specified)\n",
              long_options[option_index].name, chiMax );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", chiMax );
        break;

      case '6':
        kappaMin = atof( LALoptarg );
        if ( kappaMin < -1. )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Kappa can only take values between -1. and 1. : "
              "(%f specified)\n",
              long_options[option_index].name, kappaMin );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", kappaMin );
        break;

      case '7':
        kappaMax = atof( LALoptarg );
        if ( kappaMax > 1. )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Kappa can only take values between -1. and 1. : "
              "(%f specified)\n",
              long_options[option_index].name, kappaMax );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", kappaMax );
        break;

      case '8':
        nPointsChi = atof( LALoptarg );
        if ( nPointsChi < 1 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Number of points in the Chi direction must be greater than 0 : "
              "(%d specified)\n",
              long_options[option_index].name, nPointsChi );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", nPointsChi );
        break;

      case '9':
        nPointsKappa = atof( LALoptarg );
        if ( nPointsKappa < 1 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "Number of points in the Kappa direction must be greater than 0 : "
              "(%d specified)\n",
              long_options[option_index].name, nPointsKappa );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", nPointsKappa );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        USAGE( stderr );
        exit( 1 );
    }
  }




  if ( LALoptind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( LALoptind < argc )
    {
      fprintf ( stderr, "%s\n", argv[LALoptind++] );
    }
    exit( 1 );
  }

  /* add option without arguments into the process param table */
  if (vrbflg==1)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--verbose" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if (outCompress==1)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--write-compress" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if (computeMoments==0)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--disable-compute-moments" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if (polygonFit==0)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--disable-polygon-fit" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if (globFrameData==1)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--glob-frame-data" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if (globCalData==1)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--glob-calibration-data" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if (pointCal==1)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--point-calibration" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /*
   *
   * check validity of arguments
   *
   */

  /* if using a simulated PSD, use the first two characters of */
  /* ifoTag as an ifo prefix */
  if ( specType == specType_simulated )
  {
    if ( ! ifoTag )
    {
      fprintf( stderr, "--ifo-tag must be specified if using a simulated PSD\n" );
      exit( 1 );
    }
    else
    {
      memset( ifo, 0, sizeof(ifo) );
      memcpy( ifo, ifoTag, sizeof(ifo) - 1 );
    }
  }

  /* check validity of input data time */
  if ( ! gpsStartTime.gpsSeconds )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }
  if ( ! gpsEndTime.gpsSeconds )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  if ( gpsEndTime.gpsSeconds <= gpsStartTime.gpsSeconds )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }

  /* check validity of data length parameters */
  if ( numPoints < 0 )
  {
    fprintf( stderr, "--segment-length must be specified\n" );
    exit( 1 );
  }
  if ( ( specType==specType_mean || specType==specType_median ) && numSegments < 0 )
  {
    fprintf( stderr, "--number-of-segments must be specified if using frama data\n" );
    exit( 1 );
  }

  /* check sample rate has been given */
  if ( sampleRate < 0 )
  {
    fprintf( stderr, "--sample-rate must be specified\n" );
    exit( 1 );
  }

  /* check high pass option has been given */
  if ( highPass < 0 )
  {
    fprintf( stderr, "--disable-high-pass or --enable-high-pass (freq)"
        " must be specified\n" );
    exit( 1 );
  }
  else if ( ! highPass )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--disable-high-pass" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  else
  {
    /* check that all the high pass parameters have been specified */
    if ( highPassOrder < 0 )
    {
      fprintf( stderr, "--high-pass-order must be specified\n" );
      exit( 1 );
    }
    if ( highPassAtten < 0 )
    {
      fprintf( stderr, "--high-pass-attenuation must be specified\n" );
      exit( 1 );
    }
  }

  if ( ( specType==specType_mean || specType==specType_median ) && calData == real_8 )
  {
    /* check that strain high pass parameters have been specified */
    if ( strainHighPassFreq < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-freq must be specified for REAL8 h(t) data\n" );
      exit( 1 );
    }
    if ( strainHighPassOrder < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-order must be specified for REAL8 h(t) data\n");
      exit( 1 );
    }
    if ( strainHighPassAtten < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-atten must be specified for REAL8 h(t) data\n");
      exit( 1 );
    }
  }

  /* check validity of input data length */
  if ( specType==specType_mean || specType==specType_median )
  {
    inputDataLength = numPoints * numSegments - ( numSegments - 1 ) *
      (numPoints / 2);
    {
      UINT8 gpsChanIntervalNS = gpsEndTime.gpsSeconds * LAL_INT8_C(1000000000) -
        gpsStartTime.gpsSeconds * LAL_INT8_C(1000000000);
      UINT8 inputDataLengthNS = (UINT8) inputDataLength * LAL_INT8_C(1000000000)
        / (UINT8) sampleRate;

      if ( inputDataLengthNS != gpsChanIntervalNS )
      {
        fprintf( stderr, "length of input data and data chunk do not match\n" );
        fprintf( stderr, "start time: %d, end time %d\n",
            gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
        fprintf( stderr, "gps channel time interval: %" LAL_UINT8_FORMAT " ns\n"
            "computed input data length: %" LAL_UINT8_FORMAT " ns\n",
            gpsChanIntervalNS, inputDataLengthNS );
        exit( 1 );
      }
    }
  }

  /* check standard candle arguments */
  if ( computeCandle )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--standard-candle" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );

    if ( candleSnr < 0 )
    {
      fprintf( stderr,
          "--candle-snr must be specified if --standard-candle is given\n" );
      exit( 1 );
    }
    if ( candleMinMass < 0 )
    {
      fprintf( stderr,
          "--candle-minmass must be specified if --standard-candle is given\n" );
      exit( 1 );
    }
    if ( candleMaxMass < 0 )
    {
      fprintf( stderr,
          "--candle-maxmass must be specified if --standard-candle is given\n" );
      exit( 1 );
    }
    if ( ( specType == specType_simulated ) && numPoints < 0 )
    {
      fprintf( stderr,
          "--num-points must be specified if --standard-candle is given\n"
          "when using a simulated PSD\n" );
      exit( 1 );
    }
  }

  /* check that the spectrum generation parameters have been given */
  if ( fLow < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff must be specified\n" );
    exit( 1 );
  }
  if ( ( specType==specType_mean || specType==specType_median ) && resampFiltType < 0 )
  {
    fprintf( stderr, "--resample-filter must be specified for frame data\n" );
    exit( 1 );
  }
  if ( specType == specType_undefined )
  {
    fprintf( stderr, "--spectrum-type must be specified\n" );
    exit( 1 );
  }

  /* check for potential underflows in the simulated spectrum */
  if ( ( specType == specType_simulated ) && dynRangeExponent < 10 )
  {
    fprintf( stderr, "If using a simulated PSD, a suitable dynamic \n"
        "range exponent must be given, eg 69.0. Exiting...\n" );
    exit( 1 );
  }

  /* check that a channel has been requested and fill the ifo */
  if ( ( specType==specType_mean || specType==specType_median ) && ! fqChanName )
  {
    fprintf( stderr, "--channel-name must be specified for frame data\n" );
    exit( 1 );
  }

  /* check that we can correctly obtain the input frame data */
  if ( specType==specType_mean || specType==specType_median )
  {
    if ( globFrameData )
    {
      if ( frInCacheName )
      {
        fprintf( stderr,
            "--frame-cache must not be specified when globbing frame data\n" );
        exit( 1 );
      }

      if ( ! frInType )
      {
        fprintf( stderr,
            "--frame-type must be specified when globbing frame data\n" );
        exit( 1 );
      }
    }
    else
    {
      if ( ! frInCacheName )
      {
        fprintf( stderr,
            "--frame-cache must be specified when not globbing frame data\n" );
        exit( 1 );
      }

      if ( frInType )
      {
        fprintf( stderr, "--frame-type must not be specified when obtaining "
            "frame data from a cache file\n" );
        exit( 1 );
      }
    }
  }

  /* record the glob frame data option in the process params */
  if ( globFrameData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--glob-frame-data" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store point calibration option */
  if ( pointCal )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--point-calibration" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* check we can calibrate the data if it's not h(t) */
  if ( ( specType==specType_mean || specType==specType_median ) && ! calData )
  {
    if ( ! ( calCacheName || globCalData ) )
    {
      fprintf( stderr, "either --calibration-cache or "
          "--glob-calibration-data must be specified\n" );
      exit( 1 );
    }
    else if ( calCacheName && globCalData )
    {
      fprintf( stderr, "only one of --calibration-cache or "
          "--glob-calibration-data can be specified\n" );
      exit( 1 );
    }
  }
  else
  {
    if ( calCacheName || globCalData )
    {
      fprintf( stderr, "neither --calibration-cache nor --glob-calibration-data\n"
          "should be given when using calibrated data or a simulated PSD\n" );
      exit( 1 );
    }
  }

  /* record the glob calibration data option in the process params */
  if ( globCalData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--glob-calibration-data" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* check that the bank type has been specified */
  if ( ! haveOrder )
  {
    fprintf( stderr, "--order must be specified\n" );
    exit( 1 );
  }
  if ( ! haveApprox )
  {
    fprintf( stderr, "--approximant must be specified\n" );
    exit( 1 );
  }
  if ( ! haveSpace )
  {
    fprintf( stderr, "--space must be specified\n" );
    exit( 1 );
  }

  /* check validity of grid spacing with respect to approximant */
  if ( ! haveGridSpacing )
  {
    fprintf( stderr, "--grid-spacing must be specified\n" );
    exit( 1 );
  }
  if (gridSpacing != SquareNotOriented && gridSpacing != Hexagonal)
  {
    fprintf( stderr, "--grid-spacing must be either SquareNotOriented or Hexagonal\n" );
    exit( 1 );
  }

  /* check that the correct range parameters have been given for the bank */
  if ( approximant == BCV )
  {
    if ( ! havePsi0Min )
    {
      fprintf( stderr, "--minimum-psi0 must be specified\n" );
      exit( 1 );
    }
    if ( ! havePsi0Max )
    {
      fprintf( stderr, "--maximum-psi0 must be specified\n" );
      exit( 1 );
    }
    if ( ! havePsi3Min )
    {
      fprintf( stderr, "--minimum-psi3 must be specified\n" );
      exit( 1 );
    }
    if ( ! havePsi3Max )
    {
      fprintf( stderr, "--maximum-psi3 must be specified\n" );
      exit( 1 );
    }
    if ( ! haveAlpha )
    {
      fprintf( stderr, "--alpha must be specified\n" );
      exit( 1 );
    }
    if ( maxFcutTmplts < 0 )
    {
      fprintf( stderr, "--maximum-fcut-tmplts must be specified\n" );
      exit( 1 );
    }

    if ( psi3Max <= psi3Min )
    {
      fprintf( stdout, "invalid argument to --maximum-psi3:\n"
          "maximum value of psi3 must be greater than minimum value of psi3: "
          "(%f specified)\n",
          psi3Max );
      exit( 1 );
    }

    minMass = maxMass = 0;
  }
  else if ( approximant == BCVSpin )
  {
    if ( minMass < 0 && ( ! havePsi0Min || ! havePsi3Min ) )
    {
      fprintf( stderr, "--minimum-mass or --minimum-psi0 and --minimum-psi3 "
                       "must be specified\n" );
      exit( 1 );
    }
    if ( maxMass < 0 && ( ! havePsi0Max || ! havePsi3Max ) )
    {
      fprintf( stderr, "--maximum-mass or --maximum-psi0 and --maximum-psi3 "
                       "must be specified\n" );
      exit( 1 );
    }
  }
  else
  {
    if ( minMass < 0 )
    {
      fprintf( stderr, "--minimum-mass must be specified\n" );
      exit( 1 );
    }
    if ( maxMass < 0 && maxTotalMass < 0 )
    {
      fprintf( stderr, "Either --maximum-mass or --max-total-mass must be specified\n" );
      exit( 1 );
    }
    if ( minTotalMass > 0 && maxTotalMass < 0 )
    {
      fprintf( stderr, "--max-total-mass must be specified with --min-total-mass\n" );
      exit( 1 );
    }
  }

  /* check that the bank parameters have been specified */
  if ( minMatch < 0 )
  {
    fprintf( stderr, "--minimal-match must be specified\n" );
    exit( 1 );
  }
  if ( fUpper < 0 )
  {
    fprintf( stderr, "--high-frequency-cutoff must be specified\n" );
    exit( 1 );
  }

  /* Check that multiple cutoff freq. options specified */
  if ( ! haveNumFcut )
    {
      fprintf( stderr, "must specify --num-freq-cutoffs\n" );
      exit( 1 );
    }
  if ( ! haveMaxFcut )
    {
      fprintf( stderr, "must specify --max-high-freq-cutoff\n" );
      exit( 1 );
    }
  if ( ! haveMinFcut )
    {
      fprintf( stderr, "must specify --min-high-freq-cutoff\n" );
      exit( 1 );
    }
  /* Check Min and Max upper freq. cuts are the same if NumFreqCut = 1 */
  if ( numFreqCut == 1 )
  {
    if( maxFreqCut < minFreqCut || maxFreqCut > minFreqCut )
    {
      fprintf(stderr, "--max-high-freq-cutoff must equal --min-high-freq-cutoff when --num-freq-cutoffs = 1\n" );
      exit( 1 );
    }
  }

  if ( approximant==FindChirpPTF )
  {
    /* check max and mins are the correct way around */
    if (chiMin > chiMax )
    {
      fprintf( stderr,
          "Error: argument to --minimum-spin1 must be less than --maximum-spin1 .\n" );
      exit( 1 );
    }

    /* check that kappa min-max are set correctly */
    if (kappaMin > kappaMax)
    {
      fprintf( stderr,
          "Error: argument to --minimum-kappa1 must be less than --maximum-kappa1 .\n" );
      exit( 1 );
    }
  }

  if( etaMaxCutoff > 0 || etaMinCutoff >= 0 )
  {
    if( etaMaxCutoff <= 0 )
    {
      fprintf( stderr,
          "Error: argument --max-eta must be given if --min-eta is given\n");
      exit(1);
    }

    if( etaMinCutoff < 0 )
    {
      fprintf( stderr,
          "Error: argument --min-eta must be given if --max-eta is given\n");
      exit(1);
    }

    if( etaMaxCutoff < etaMinCutoff )
    {
      fprintf( stderr,
            "Error: value for --max-eta must be greater than or equal to value for --min-eta\n");
      exit(1);
    }
  }

  return 0;
}

#undef ADD_PROCESS_PARAM

