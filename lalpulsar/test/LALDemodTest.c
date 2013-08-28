/*
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Teviet Creighton
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
 * \author Berukoff, S.J., Papa, M.A.,
 * \file
 * \ingroup pulsarTODO
 *
 * \brief Performs required tests of LALDemod().
 *
 * \heading{Usage}
 * \code
 * LALDemodTest -i <input data file> [-d <gap>] [-n] [-o]
 * \endcode
 *
 * \heading{Description}
 *
 * This routine performs tests on the routine <tt>LALDemod()</tt>.
 * Options:
 * <ul>
 * <li> <tt>-i</tt> -- the input data file (default is 'in.data'; an example is included, format below)</li>
 * <li> <tt>-n</tt> -- add zero-mean Gaussian noise to the signal</li>
 * <li> <tt>-d \<gap\></tt> -- simulate gaps in the data.  The number <tt>\<gaps\></tt> refers to the integral number of SFT timescales between adjacent timestamps.</li>
 * <li> <tt>-o</tt> -- print out result data files</li>
 * </ul>
 *
 * Structure:
 * In more detail, let us begin with a discussion of the structure of the test
 * code, which is composed of several modules.
 * <ul>
 * <li> The first module reads in data from an input parameter data file.
 * The parameters must be listed in the input file in the following order,
 * with the corresponding format:
 * \code
 * total observation time -- float
 * coherent search time -- float
 * factor by which to modify SFT timescale -- float
 * Cross amplitude -- float
 * Plus amplitude -- float
 * DeFT frequency band (centered by default around f0) -- float
 * f0, intrinsic frequency of the signal at the beginning of the
 * observation -- float
 * maximum order of signal spindown parameters -- int
 * signal spindown parameter 1 --  scientific notation
 * signal spindown parameter 2 --  scientific notation
 * signal spindown parameter 3 --  scientific notation
 * signal spindown parameter 4 --  scientific notation
 * signal spindown parameter 5 --  scientific notation
 * signal source right ascension (alpha) -- float (value in DEGREES)
 * signal source declination (delta) -- float (value in DEGREES)
 * maximum order of template spindown parameters -- int
 * template spindown parameter 1 --  scientific notation (NOT
 * template spindown parameter 2 --  scientific notation (NOT scaled by f0)
 * template spindown parameter 3 --  scientific notation (NOT scaled by f0)
 * template spindown parameter 4 --  scientific notation (NOT scaled by f0)
 * template spindown parameter 5 --  scientific notation (NOT scaled by f0)
 * template source right ascension (alpha) -- float (value in DEGREES)
 * template source declination (delta) -- float (value in DEGREES)
 * \endcode
 * Note: Above, the *signal* spindown parameters are scaled by the intrinsic frequency, while the
 * template* spindown parameters are not.  This is due to the difference in definitions between the
 * PulsarSimulateCoherentGW() package, which generates the signal, and this package.
 *
 * </li><li> The next module in the test code, which is optionally executed with
 * the '<tt>-n</tt>' switch, creates noise using LALs
 * <tt>LALNormalDeviates()</tt> routine.  By design, the noise is created in
 * single
 * precision, and is zero-mean and Gaussian.  This noise is added,
 * datum-by-datum, to the time series created in
 * the next module, after the amplitude of the time series has been
 * changed by a
 * factor of \c SNR, which is specified in the input data file.
 *
 * </li><li> The next module to be invoked creates a time series, according to
 * the standard model for pulsars with spindown.  This is done by using the
 * <tt>LALGenerateTaylorCW()</tt> and <tt>LALPulsarSimulateCoherentGW()</tt> functions.  This time series
 * undergoes an FFT, and this transformed data then constitutes the SFT data to be
 * input to the demodulation code.  The fake signal data is characterized
 * by an intrinsic frequency at the beginning of the observation plus some
 * other source parameters.  The DeFT is produced in a band \c f0Band
 * (as specified as an input parameter) and centered at this frequency.
 * The width of the band (plus some extra width of \f$2\cdot 10^{-4}f0 \f$ Hz) determines the
 * sampling frequency of the time series (Nyquist theorem).  In practice this
 * would be the inverse FFT of a data set that has been band-passed around
 * \c f0 and then appropriately down-sampled (e.g. with a lock-in).  The
 * normalization rule for FFT data is the following: if sinusoidal data over a
 * time \f$T\f$ and with amplitude \f$A\f$ is FFT-ed, the sum of the square amplitude of
 * the output of the FFT (power) is equal to \f${A^2 T}\f$.  Thus, the power peak at
 * the sinusoids frequency should be expected to be \f$\sim\f$ \f$\frac{A^2}{2} T\f$,
 * within a factor of 2.  The same normalization rule applies to the DeFT data.
 * Thus by piecing together \f$N\f$ SFTs we expect a DeFT power peak \f$\sim N\f$ higher
 * than that of the SFTs - at least in the case of perfect signal-template match.
 *
 * Let us now spend a few words on the choice of the SFT time baseline.  Given an
 * intrinsic search frequency one can compute the longest time baseline which is
 * still compatible with the requirement that the instantaneous signal frequency
 * during such time baseline does not shift by more than a frequency bin.  This
 * is the default choice for the SFT length, having assumed that the modulation
 * is due to the spin of the Earth and having taken a simple epicyclic model to
 * evaluate the magnitude of this effect.  It is possible to choose a different
 * time baseline by specifying a value for
 * the variable \c gap other than 1.  Note that the SFT time baseline is
 * approximated to the nearest value such that the number of SFT samples is a
 * power of two.  This is also well documented in the code.
 *
 * The set of SFTs does not necessarily come from contiguous data sets: a set of
 * time stamps is created that defines the time of the first sample of each SFT
 * data chunk.  The timestamps which are required in many parts of the code are
 * generated in a small subroutine <tt>times2()</tt>.  This routine takes as input
 * the SFT timescale \c tSFT, the number of SFTs which will be created,
 * \c mObsSFT, and a switch which lets the code know whether to make even
 * timestamps, or timestamps with gaps (see below for more on this).  The
 * subroutine then writes the times to the \c LIGOTimeGPS vector containing
 * the timestamps for the entire test code, and returns this vector.  Note that
 * each datum of the  \c LIGOTimeGPS vector is comprised of two fields; if
 * accessing the \f$i^{th}\f$ datum, the seconds part of the timestamp vector
 * \c ts is <tt>ts[i].gpsSeconds</tt> and the nanoseconds part is
 * <tt>ts[i].gpsNanoSeconds</tt>.  These are the fields which are written in this
 * <tt>times()</tt>.
 *
 * As an important side note, let us discuss the effect that a vector of
 * timestamps with gaps has on the resulting transformed data.  Since each of the
 * timestamps refers to the first datum of each SFT, the existence of the gaps
 * means that instead of transforming a continuous set of data, we are reduced to
 * transforming a piecewise continuous set.  Since we can envision these gaps as
 * simply replacing real data with zeros, we correspondingly should see a power
 * loss in the resulting FFTs signal bins and a broadening of the power spectrum.
 * Since real detectors will clearly have gaps in the data, this effect is
 * obviously something we seek to minimize or eliminate if possible.  This work
 * continues to be under development.
 *
 * The total observation time determines how many SFTs and how many DeFTs are
 * created.  The actual time baseline for both the DeFTs and the total
 * observation time might differ from the ones defined in the input file, the
 * reason being that they are rounded to the nearest multiple of the SFT time
 * baseline.
 *
 * Note that use is made of the <tt>LALBarycenter()</tt> routine (see
 * \ref LALBarycenter.h), which (among other things) provides,  at any given
 * time, the actual instantaneous position and velocity of a  detector at any
 * specified location of the Earth with respect to the SSB.
 *
 * </li><li> Following the creation of a short chunk of time series data, an FFT is
 * performed with the internal FFTW routines.  This outputs a frequency domain
 * chunk which is placed into the \c SFTData array of structures.  This will
 * contain all of the SFT data we need to demodulate, and in the future, will be
 * the storage area for the real data.
 *
 * </li><li> The next module begins the demodulation process.  First, the parameters
 * for the demodulation routine are assigned from values previously calculated in
 * the test code.  Similarly, parameters for the <tt>LALComputeSky()</tt> routine are
 * assigned.  This routine computes the coefficients \f$A_{s\alpha}\f$ and
 * \f$B_{s\alpha}\f$ (see \ref ComputeSky.h) of the spindown parameters
 * for the phase model weve assumed.  These coefficients are used within the
 * <tt>LALDemod()</tt> routine itself. Since they only depend on the template sky
 * position, in a search over many different spin-down parameters they are
 * reused, thus one needs compute them only once.  Then, the <tt>LALComputeAM()</tt>
 * routine is called, to calculate the amplitude modulation filter information.  Finally, at last, the
 * demodulation routine itself is called, and, if the command line option
 * '<tt>-o</tt>' is used,  output are several data files containing demodulated
 * data (these are by default named '<tt>xhat_#</tt>').  These output files have two columns, one
 * for the value of the periodogram and one for the frequency.
 *
 * </li></ul>
 *
 * \heading{Exit codes}
 *
 * \heading{Uses}
 * \code
 * lalDebugLevel
 * LALMalloc()
 * LALFopen()
 * LALFclose()
 * LALSCreateVector()
 * LALCreateRandomParams()
 * LALNormalDeviates()
 * LALDestroyRandomParams()
 * LALSDestroyVector()
 * LALCCreateVector()
 * LALCreateForwardRealFFTPlan()
 * LALREAL4VectorFFT()
 * LALCDestroyVector()
 * LALDestroyRealFFTPlan()
 * LALGenerateTaylorCW()
 * LALPulsarSimulateCoherentGW()
 * LALComputeSky()
 * LALFree()
 * LALDemod()
 * LALBarycenter()
 * LALComputeAM()
 * \endcode
 *
 * \heading{Notes}
 * The implementation of the code here is intended to give a general outline of
 * what the demodulation code needs to work.  Most of this test function performs
 * steps (e.g., noise, time- and frequency-series generation) that will be already
 * present in the data.
 *
 */

#ifndef LALDEMODTEST_C
#define LALDEMODTEST_C
#endif

/* Usage */
#define USAGE "Usage: %s [-i basicInputsFile] [-n] [-d] [-o]\n"


/* Error macro, taken from ResampleTest.c in LAL's pulsar/test package */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,		\
		 __LINE__, "$Id$", statement ? statement : "",\
		 (msg) );                                            \
}                                                                    \
else (void)(0)

#include <lal/LALDemod.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SeqFactories.h>
#include <lal/FileIO.h>
#include <unistd.h>
#include <sys/types.h>


static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);

static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f);

static void times(REAL8 , INT4, LIGOTimeGPS *, INT4 );

static void times2(REAL8 tSFT, INT4 howMany, LIGOTimeGPS **ts, INT4 **sftPerCoh, INT4 sw, INT4 mCohSFT);


static PulsarCoherentGW emptySignal;

int main(int argc, char **argv)
{
  static LALStatus status;

  /***** VARIABLE DECLARATION *****/

  ParameterSet *signalParams;
  ParameterSet *templateParams;
  char *basicInputsFile;
  FILE *bif;
  REAL8 tObs, tCoh, tSFT;
  REAL8 oneOverSqrtTSFT;
  REAL8 aCross, aPlus, SNR;
  REAL8 f0;

  INT4 mCohSFT, mObsCoh, mObsSFT;
  REAL8 dfSFT, dt;
  INT4 if0Min, if0Max, ifMin, ifMax;
  REAL8 f0Min, f0Max, fMin, f0Band, fWing;
  INT4 nDeltaF,ntermsdivbytwo;
  LALFstat Fstat;

  LIGOTimeGPS *timeStamps;

  REAL4Vector *temp = NULL;
  REAL4Vector *noise = NULL;
  static RandomParams *params;
  INT4 seed=0;

  INT4 a,i,k;

  FFT **SFTData;
  RealFFTPlan *pfwd = NULL;
  COMPLEX8Vector *fvec = NULL;

  DemodPar *demParams;
  CSParams *csParams;
  INT4 iSkyCoh;

  DeFTPeriodogram **xHat;

  const CHAR *noi=NULL;
  INT4 deletions=0;
  const CHAR *output=NULL;

  REAL8 factor;
  INT2 arg;

  /* file name if needed for SFT output files */
  CHAR *sftoutname=NULL;

  /* Quantities for use in LALBarycenter package */
  BarycenterInput baryinput;
  EmissionTime emit;
  EarthState earth;
  LALDetector cachedDetector;

  EphemerisData *edat=NULL;
  char earthEphemeris[]="earth98.dat";
  char sunEphemeris[]="sun98.dat";

  REAL4TimeSeries *timeSeries = NULL;

  /*
   *  Variables for AM correction
   */
  PulsarCoherentGW cgwOutput = emptySignal;
  TaylorCWParamStruc genTayParams;
  PulsarDetectorResponse cwDetector;
  AMCoeffsParams *amParams;
  AMCoeffs amc;
  INT4 *sftPerCoh;
  INT4 totalSFT=0;

#define DEBUG 1
#if (0)

  INT2 tmp=0;
  REAL8 dtSFT, dfCoh;
  REAL8 fMax;
  INT4 nSFT;

#endif

  /***** END VARIABLE DECLARATION *****/


  /***** PARSE COMMAND LINE OPTIONS *****/
  basicInputsFile=(char *)LALMalloc(50*sizeof(char));
  snprintf(basicInputsFile, 9, "in.data\0");
  arg=1;
  while(arg<argc) {

    /* the input file */
    if(!strcmp(argv[arg],"-i"))
      {
	strcpy(basicInputsFile, argv[++arg]);
	arg++;
	if(LALOpenDataFile(basicInputsFile)==NULL)
	  {
	    ERROR(LALDEMODH_ENOFILE, LALDEMODH_MSGENOFILE, 0);
	    XLALPrintError(USAGE, *argv);
	    return LALDEMODH_ENOFILE;
	  }
      }

    /* turn noise on? if string is not NULL, it will be */
    else if(!strcmp(argv[arg],"-n"))
      {
	noi="n";
	arg++;
      }

    /* make SFT file output (such that the driver code can read in as input data) */
    else if (!strcmp(argv[arg],"-m"))
      {
      sftoutname=argv[++arg];
      arg++;
      }
    /* turn timestamp deletions on? if string is not NULL, it will be */
    else if(!strcmp(argv[arg],"-d"))
      {
	deletions=1;
	arg++;
      }

    /* turn on output? if string is not NULL, it will be */
    else if(!strcmp(argv[arg],"-o"))
      {
	output="o";
	arg++;
      }

    /* default: no input file specified */
    else if(basicInputsFile==NULL)
      {
	bif=LALOpenDataFile("in.data");
      }

    /* erroneous command line argument */
    else if(basicInputsFile!=NULL && arg<argc)
      {
	ERROR(LALDEMODH_EBADARG, LALDEMODH_MSGEBADARG, 0);
	XLALPrintError(USAGE, *argv);
	arg=argc;
	return LALDEMODH_EBADARG;
      }
  }

  /***** END COMMAND LINE PARSING *****/


  /***** INITIALIZATION OF SIGNAL AND TEMPLATE *****/

  /* Allocate space for signal parameters */
  signalParams=LALMalloc(sizeof(ParameterSet));
  signalParams->spind=LALMalloc(sizeof(Spindown));
  signalParams->skyP=LALMalloc(2*sizeof(SkyPos));
  signalParams->spind->spParams=LALMalloc(5*sizeof(REAL8));

  /* Allocate space for template parameters */
  templateParams=LALMalloc(sizeof(ParameterSet));
  templateParams->spind=LALMalloc(sizeof(Spindown));
  templateParams->skyP=LALMalloc(2*sizeof(SkyPos));
  templateParams->spind->spParams=LALMalloc(5*sizeof(REAL8));

  /***** END INITIALIZATION *****/


  /***** GET INPUTS FROM FILES *****/

  bif=LALOpenDataFile(basicInputsFile);

  fscanf(bif, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n%le\n%le\n%le\n%le\n%le\n%lf\n"
	 " %lf\n%d\n%le\n%le\n%le\n%le\n%le\n%lf\n%lf\n",
	 &tObs, &tCoh, &factor, &aCross, &aPlus, &SNR, &f0Band, &f0,
	 &signalParams->spind->m,
	 &signalParams->spind->spParams[0], &signalParams->spind->spParams[1],
	 &signalParams->spind->spParams[2], &signalParams->spind->spParams[3],
	 &signalParams->spind->spParams[4],
	 &signalParams->skyP->alpha, &signalParams->skyP->delta,
	 &templateParams->spind->m,
	 &templateParams->spind->spParams[0], &templateParams->spind->spParams[1],
	 &templateParams->spind->spParams[2], &templateParams->spind->spParams[3],
	 &templateParams->spind->spParams[4],
	 &templateParams->skyP->alpha, &templateParams->skyP->delta);

  LALFclose(bif);

  /***** END FILE INPUT *****/


  /***** CALCULATE USEFUL QUANTITIES *****/

  /*	2^n gives the number of samples of the SFT. This will be a signal in a 	*/
  /*	band fSample/2 with fSample=2*(2.e-4 f0+f0Band). This makes sure that	*/
  /*	that the sampling frequency enables us to produce the DeFT frequency 	*/
  /*	band that we want (f0Band) and that there's enought wings of SFT data	*/
  /*	(4.e-4*f0) to produce such DeFT band. 					*/
  /* 	The main formula used here is that the gives the maximum allowed	*/
  /*	time-baseline (Tmax) such that a signal of frequency f0 Doppler 	*/
  /*	modulated due to Earth spin is seen as monochromatic during such	*/
  /*	observation time: 							*/
  /* 	  	Tmax < 9.6e4/sqrt(f0). 						*/
  /*	We have thus taken Tmax=9.4e4/sqrt(f0). A different criterion can be	*/
  /*	chosen by setting the variable "factor" to a suitable value. We have	*/
  /*	rounded the resulting number of samples to the nearest smallest power 	*/
  /*	of two. 							       	*/

  {
    REAL8 tempf0Band;

    ntermsdivbytwo=64;
    /* compute size of SFT wings */
    fWing = 2.0e-4 * f0;
    /* Adjust search band to include wings */
    f0Band += 2.0*fWing;
    tempf0Band = f0Band;
    tSFT=1.0e3;/*9.6e4/sqrt(f0);*/
    nDeltaF = 2*ntermsdivbytwo+(INT4)(ceil(f0Band*tSFT));
      /* tSFT = (REAL8)nDeltaF/f0Band;*/
    oneOverSqrtTSFT = 1.0/sqrt(tSFT);
    /* Number of data in time series for 1 SFT */
    dfSFT=1.0/tSFT;
    /* Get rid of the wings */
    f0Band -= 2.0*fWing;
    /* the index of the right side of the band, NO WINGS */
    if0Max=ceil((f0+f0Band/2.0)/dfSFT);
    /* the index of the left side of the band, NO WINGS */
    if0Min=floor((f0-f0Band/2.0)/dfSFT);
    /* frequency of the right side of the band, NO WINGS */
    f0Max=dfSFT*if0Max;
    /* frequency of the left side of the band, NO WINGS */
    f0Min=dfSFT*if0Min;

    f0Band = tempf0Band;
  }
  /* the hard-wired 64 is for the necessary NTERMS_COH_DIV_TWO*/
  /* index of the left side of the band, WITH WINGS */

  ifMin=floor(f0/dfSFT)-nDeltaF/2;
  /* indexof the right side of the band, WITH WINGS */
  ifMax=ifMin+nDeltaF;
  /* frequency of the left side of the band, WITH WINGS */
  fMin=dfSFT*ifMin;
  nDeltaF=2*ntermsdivbytwo+(INT4)(ceil(f0Band/dfSFT));
  /* Number of SFTs which make one DeFT */
  mCohSFT=ceil(tCoh/tSFT);
  if (mCohSFT==0) mCohSFT=1;
  /* Coherent search time baseline */
  tCoh=tSFT*mCohSFT;
  /* Number of coherent timescales in total obs. time */
  mObsCoh=floor(tObs/tCoh);
  if (mObsCoh ==0) mObsCoh=1;
  /* Total observation time */
  tObs=tCoh*mObsCoh;
  /* Number of SFTs needed during observation time */
  mObsSFT=mCohSFT*mObsCoh;

  printf("%d %d %d\n",mObsSFT,mObsCoh,mCohSFT);

  /* convert input angles from degrees to radians */
  signalParams->skyP->alpha=signalParams->skyP->alpha*LAL_TWOPI/360.0;
  signalParams->skyP->delta=signalParams->skyP->delta*LAL_TWOPI/360.0;
  templateParams->skyP->alpha=templateParams->skyP->alpha*LAL_TWOPI/360.0;
  templateParams->skyP->delta=templateParams->skyP->delta*LAL_TWOPI/360.0;

  /* Convert signal spindown to conform with GenerateTaylorCW() definition */
  for(i=0;i<signalParams->spind->m;i++){signalParams->spind->spParams[i] /= f0;}

  /***** END USEFUL QUANTITIES *****/


  /***** CALL ROUTINE TO GENERATE TIMESTAMPS *****/

  /*times(tSFT, mObsSFT, timeStamps, deletions);*/
  i=0;
  times2(tSFT, mObsCoh, &timeStamps, &sftPerCoh, deletions, mCohSFT);
  /*
    In order to tell many of the loops in this code how many times to iterate,
     we set a variable which is equal to the last entry in the sftPerCoh array.
     Note, that for no gaps, this number is simply mObsSFT; with gaps, it's
     something else.  Regardless, this is the value that the subroutines know
     as mObsSFT.
  */
  totalSFT = sftPerCoh[mObsCoh];
    /***** END TIMESTAMPS *****/


  /***** CREATE NOISE *****/

  if(noi!=NULL)
    {
      LALCreateVector(&status, &noise, (UINT4)nDeltaF*mObsSFT);
      LALCreateVector(&status, &temp, (UINT4)nDeltaF*mObsSFT);
      LALCreateRandomParams(&status, &params, seed);

      /* fill temp vector with normal deviates*/
      LALNormalDeviates(&status, temp, params);

      /* rewrite into COMPLEX8 VECTOR *noise */
      i=0;
      while(i<nDeltaF*mObsSFT)
	{
	  noise->data[i]=temp->data[i];
	  i++;
	}

      /* Destroy structures that were used here, if we can */
      LALDestroyRandomParams(&status, &params);
      LALSDestroyVector(&status, &temp);
    }

  /***** END CREATE NOISE *****/

  /* Quantities computed for barycentering */
  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = earthEphemeris;
  (*edat).ephiles.sunEphemeris = sunEphemeris;

  /* Read in ephemerides */
  LALInitBarycenter(&status, edat);
  /*Getting detector coords from DetectorSite module of tools package */

  /* Cached options are:
     LALDetectorIndexLHODIFF, LALDetectorIndexLLODIFF,
     LALDetectorIndexVIRGODIFF, LALDetectorIndexGEO600DIFF,
     LALDetectorIndexTAMA300DIFF,LALDetectorIndexCIT40DIFF
  */
  cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  baryinput.site.location[0]=cachedDetector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=cachedDetector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=cachedDetector.location[2]/LAL_C_SI;
  baryinput.alpha=signalParams->skyP->alpha;
  baryinput.delta=signalParams->skyP->delta;
  baryinput.dInv=0.e0;

  /***** CREATE SIGNAL *****/

  /* Set up parameter structure */
  /* Source position */
  genTayParams.position.latitude  = signalParams->skyP->alpha;
  genTayParams.position.longitude = signalParams->skyP->delta;
  genTayParams.position.system = COORDINATESYSTEM_EQUATORIAL;
  /* Source polarization angle */
  /* Note, that the way we compute the statistic, we don't care what this value is. */
  genTayParams.psi = 0.0;
  /* Source initial phase */
  genTayParams.phi0 = 0.0;
  /* Source polarization amplitudes */
  genTayParams.aPlus = aPlus;
  genTayParams.aCross = aCross;
  /* Source intrinsic frequency */
  genTayParams.f0 = f0;
  /* Source spindown parameters: create vector and assign values */
  genTayParams.f = NULL;
  LALDCreateVector(&status, &(genTayParams.f), signalParams->spind->m);
  for(i=0;i<signalParams->spind->m;i++){genTayParams.f->data[i] = signalParams->spind->spParams[i];}
  /* Time resolution for source */
  /* Note, that this needs to be sampled only very coarsely! */
  genTayParams.deltaT = 100.0;
  /* Length should include fudge factor to allow for barycentring */
  genTayParams.length = (tObs+1200.0)/genTayParams.deltaT;
  memset(&cwDetector, 0, sizeof(PulsarDetectorResponse));
  /* The ephemerides */
  cwDetector.ephemerides = edat;
  /* Specifying the detector site (set above) */
  cwDetector.site = &cachedDetector;
  /* The transfer function.
   * Note, this xfer function has only two points */
  cwDetector.transfer = (COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
  memset(cwDetector.transfer, 0, sizeof(COMPLEX8FrequencySeries));
  cwDetector.transfer->epoch = timeStamps[0];
  cwDetector.transfer->f0 = 0.0;
  cwDetector.transfer->deltaF = 16384.0;
  cwDetector.transfer->data = NULL;
  LALCCreateVector(&status, &(cwDetector.transfer->data), 2);
  /* Allocating space for the eventual time series */
  timeSeries = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
  timeSeries->data = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));

  /* The length of the TS will be plenty long; below, we have it set so that
   * it is sampled at just better than twice the Nyquist frequency, then
   * increased that to make it a power of two (to facilitate a quick FFT)
   */
  {
    INT4 lenPwr;
    lenPwr = ceil(log(tSFT*2.01*f0)/log(2.0));
    timeSeries->data->length = pow(2.0,lenPwr);
  }

  timeSeries->data->data = (REAL4 *)LALMalloc(timeSeries->data->length*sizeof(REAL4));
  dt = timeSeries->deltaT = tSFT/timeSeries->data->length;

  /* unit response function */
  cwDetector.transfer->data->data[0].re = 1.0;
  cwDetector.transfer->data->data[1].re = 1.0;
  cwDetector.transfer->data->data[0].im = 0.0;
  cwDetector.transfer->data->data[1].im = 0.0;
  /* again, the note about the barycentring */
  genTayParams.epoch.gpsSeconds = timeStamps[0].gpsSeconds - 600;
  genTayParams.epoch.gpsNanoSeconds = timeStamps[0].gpsNanoSeconds;
  memset(&cgwOutput, 0, sizeof(PulsarCoherentGW));

  /*
   *
   * OKAY, GENERATE THE SIGNAL @ THE SOURCE
   *
   */
  LALGenerateTaylorCW(&status, &cgwOutput, &genTayParams);

  {
    INT4 len,len2;

    len=timeSeries->data->length;
    len2=len/2+1;

    /* Create vector to hold frequency series */
    LALCCreateVector(&status, &fvec, (UINT4)len2);

    /* Compute measured plan for FFTW */
    LALCreateForwardRealFFTPlan(&status, &pfwd, (UINT4)len, 0);

    /* Allocate memory for the SFTData structure */
    /*
     * Note that the length allocated for the data
     * array is the width of the band we're interested in,
     * plus f0*4E-4 for the 'wings', plus a bit of wiggle room.
     */
    SFTData=(FFT **)LALMalloc(totalSFT*sizeof(FFT *));
    for(i=0;i<totalSFT;i++){
      SFTData[i]=(FFT *)LALMalloc(sizeof(FFT));
      SFTData[i]->fft=(COMPLEX8FrequencySeries *)
	LALMalloc(sizeof(COMPLEX8FrequencySeries));
      SFTData[i]->fft->data=(COMPLEX8Vector *)
	LALMalloc(sizeof(COMPLEX8Vector));
      SFTData[i]->fft->data->data=(COMPLEX8 *)
	LALMalloc((nDeltaF+1)*sizeof(COMPLEX8));
    }

    /* Lots of debugging stuff.  If you want to use these, go ahead, but
     * be warned that these files may be HUGE (10s-100s of GB).
     */

    {
      REAL4Vector *tempTS = NULL;

      LALSCreateVector(&status, &tempTS, (UINT4)len);

      for(a=0;a<totalSFT;a++)
	{
	  REAL4Vector *ts = timeSeries->data;
	  /*
	   *  Note that this epoch is different from the epoch of GenTay.
	   *  The difference is seconds, a little bit more than the
	   *  maximum propagation delay of a signal between the Earth
	   *  and the solar system barycentre.
	   */
	  timeSeries->epoch = timeStamps[a];
	  /*
	   *
	   *  OKAY, CONVERT SOURCE SIGNAL TO WHAT SOMEBODY SEES AT THE DETECTOR OUTPUT
	   *
	   */
	  LALPulsarSimulateCoherentGW(&status, timeSeries, &cgwOutput, &cwDetector);

	  /* Write time series of correct size to temp Array, for FFT */
	  for(i=0;i<len;i++)
	    {
	      tempTS->data[i] = ts->data[i];
	    }

	  /* Perform FFTW-LAL Fast Fourier Transform */
	  LALForwardRealFFT(&status, fvec, tempTS, pfwd);

	  {
	    INT4 fL = ifMin;
	    INT4 cnt = 0;
	    while(fL < ifMin+nDeltaF+1)
	      {
		COMPLEX8 *tempSFT=SFTData[a]->fft->data->data;
		COMPLEX8 *fTemp=fvec->data;

		/* Also normalise */
		tempSFT[cnt].re = fTemp[fL].re * oneOverSqrtTSFT;
		tempSFT[cnt].im = fTemp[fL].im * oneOverSqrtTSFT;

		cnt++; fL++;

	      }

	  }
	  /* assign particulars to each SFT */
	  SFTData[a]->fft->data->length = nDeltaF+1;
	  SFTData[a]->fft->epoch = timeStamps[a];
	  SFTData[a]->fft->f0 = f0Min; /* this is the frequency of the first freq in the band */
	  SFTData[a]->fft->deltaF = dfSFT;
	  printf("Created SFT number %d of %d\n.",a, mObsSFT);

	}
      LALSDestroyVector(&status, &tempTS);
    }

    /*
     * Note, we have to destroy the memory that GenTay allocates.
     * This is currently (04.02) not documented!
     */

    if(&(cgwOutput.a) !=NULL) {
      LALSDestroyVectorSequence(&status, &(cgwOutput.a->data));
      LALFree(cgwOutput.a);
    }
    if(&(cgwOutput.f) !=NULL) {
      LALSDestroyVector(&status, &(cgwOutput.f->data));
      LALFree(cgwOutput.f);
    }
    if(&(cgwOutput.phi) !=NULL) {
      LALDDestroyVector(&status, &(cgwOutput.phi->data));
      LALFree(cgwOutput.phi);
    }

    if(noi!=NULL) {LALDestroyVector(&status, &noise);}
    LALCDestroyVector(&status, &fvec);
    LALFree(timeSeries->data->data);
    LALFree(timeSeries->data);
    LALFree(timeSeries);
    LALDestroyRealFFTPlan(&status, &pfwd);
  }
  LALDDestroyVector(&status,&(genTayParams.f));
  LALCDestroyVector(&status, &(cwDetector.transfer->data));
  LALFree(cwDetector.transfer);

  /***** END CREATE SIGNAL *****/


 /* BEGIN AMPLITUDE MODULATION */

  /* Allocate space for amParams stucture */
  /* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  /* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth;
  amParams->edat = edat;
  amParams->das->pDetector = &cachedDetector;
  amParams->das->pSource->equatorialCoords.latitude = signalParams->skyP->alpha;
  amParams->das->pSource->equatorialCoords.longitude = signalParams->skyP->delta;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = genTayParams.psi;
  amParams->mObsSFT = totalSFT;
 /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  LALSCreateVector(&status, &(amc.a), (UINT4)totalSFT);
  LALSCreateVector(&status, &(amc.b), (UINT4)totalSFT);

 /*
  * Compute timestamps for middle of each ts chunk
  * Note, we have decided to use the values of 'a' and 'b'
  * that correspond to the midpoint of the timeSeries for
  * which the SFT denotes.  In practice, this can be a flaw
  * since 'a' is a sinusoid with T=45ksec and T_b=90ksec; both
  * of these have amplitude about 0.4.  Thus, for 'a', on a timescale
  * of 10ksec, the value changes significantly.
  */
  {
    LIGOTimeGPS *midTS;
	REAL8 t;

	midTS = (LIGOTimeGPS *)LALCalloc(totalSFT,sizeof(LIGOTimeGPS));
    /* compute timestamps of end and beg */
    amParams->ts0 = timeStamps[0];
    TimeToFloat(&t, &(timeStamps[0]));
    t+=tObs;
    FloatToTime(&(amParams->tsEnd), &t);

    for(k=0; k<totalSFT; k++)
      {
	REAL8 teemp=0.0;

	TimeToFloat(&teemp, &(timeStamps[k]));
	teemp += 0.5*tSFT;
	FloatToTime(&(midTS[k]), &teemp);
      }
    /* Compute the AM coefficients */
    LALComputeAM(&status, &amc, midTS, amParams);
    LALFree(midTS);
  }

  /***** DEMODULATE SIGNAL *****/

 /* Allocate space and set quantity values for demodulation parameters */
  demParams=(DemodPar *)LALMalloc(sizeof(DemodPar));
  demParams->skyConst=(REAL8 *)LALMalloc((2*templateParams->spind->m *
					  (totalSFT+1)+2*totalSFT+3)*sizeof(REAL8));
  demParams->spinDwn=(REAL8 *)LALMalloc(templateParams->spind->m*sizeof(REAL8));
  demParams->if0Max=if0Max;
  demParams->if0Min=if0Min;
  demParams->mCohSFT=mCohSFT;
  demParams->mObsCoh=mObsCoh;
  demParams->ifMin=ifMin;
  demParams->spinDwnOrder=templateParams->spind->m;
  demParams->amcoe = &amc;
  demParams->sftPerCoh = sftPerCoh;
  for(i=0;i<signalParams->spind->m;i++)
    {
      demParams->spinDwn[i]=templateParams->spind->spParams[i];
    }

  /* Allocate space and set quantities for call to LALComputeSky() */
  csParams=(CSParams *)LALMalloc(sizeof(CSParams));
  csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csParams->skyPos[0]=templateParams->skyP->alpha;
 csParams->skyPos[1]=templateParams->skyP->delta;
 csParams->tGPS=timeStamps;
 csParams->spinDwnOrder=templateParams->spind->m;
 csParams->mObsSFT=totalSFT;
 csParams->tSFT=tSFT;
 csParams->edat=edat;
 csParams->emit=&emit;
 csParams->earth=&earth;
 csParams->baryinput=&baryinput;

 iSkyCoh=0;

  /* Call COMPUTESKY() */
  LALComputeSky(&status, demParams->skyConst, iSkyCoh, csParams);

  /* Deallocate space for ComputeSky parameters */
  LALFree(csParams->skyPos);
  LALFree(csParams);

  /* Allocate memory for demodulated data */
  xHat=(DeFTPeriodogram **)LALMalloc(mObsCoh*sizeof(DeFTPeriodogram *));
  for(i=0; i<mObsCoh; i++)
    {
      xHat[i]=(DeFTPeriodogram *)LALMalloc(sizeof(DeFTPeriodogram));
      xHat[i]->fft=(REAL8FrequencySeries *)
	LALMalloc(sizeof(REAL8FrequencySeries));
      xHat[i]->fA=(COMPLEX16FrequencySeries *)LALMalloc(sizeof(COMPLEX16FrequencySeries));
      xHat[i]->fB=(COMPLEX16FrequencySeries *)LALMalloc(sizeof(COMPLEX16FrequencySeries));
      xHat[i]->fft->data=(REAL8Vector *)LALMalloc(sizeof(REAL8Vector));
      xHat[i]->fft->data->data=(REAL8 *)LALMalloc((UINT4)((if0Max-if0Min+1)*mCohSFT)*sizeof(REAL8));
      xHat[i]->fft->data->length=(UINT4)((if0Max-if0Min+1)*mCohSFT);
      xHat[i]->fA->data=(COMPLEX16Vector *)LALMalloc(sizeof(COMPLEX16Vector));
      xHat[i]->fA->data->data=(COMPLEX16 *)LALMalloc((UINT4)((if0Max-if0Min+1)*mCohSFT)*sizeof(COMPLEX16));
      xHat[i]->fA->data->length=(UINT4)((if0Max-if0Min+1)*mCohSFT);
      xHat[i]->fB->data=(COMPLEX16Vector *)LALMalloc(sizeof(COMPLEX16Vector));
      xHat[i]->fB->data->data=(COMPLEX16 *)LALMalloc((UINT4)((if0Max-if0Min+1)*mCohSFT)*sizeof(COMPLEX16));
      xHat[i]->fB->data->length=(UINT4)((if0Max-if0Min+1)*mCohSFT);
    }

  for(k=0; k<mObsCoh; k++)
    {
      demParams->iCoh=k;

      /**************************/
      /*       DEMODULATE       */
      /**************************/
      Fstat.F = xHat+k;
      /*      LALDemod(&status, *(xHat+k), SFTData, demParams); */
      LALDemod(&status, &Fstat, SFTData, demParams);

    }

/***** END DEMODULATION *****/


  /***** DEALLOCATION *****/

  /* Deallocate AM  */
  LALSDestroyVector(&status, &(amc.a));
  LALSDestroyVector(&status, &(amc.b));
  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);

  /* Deallocate SFTData structure, since we don't need it anymore */
  for(i=0;i<totalSFT;i++)
    {
      LALFree(SFTData[i]->fft->data->data);
      LALFree(SFTData[i]->fft->data);
      LALFree(SFTData[i]->fft);
      LALFree(SFTData[i]);
    }
  LALFree(SFTData);

  /* Deallocate memory used by demodulated data */
  for(i=0; i<mObsCoh; i++)
    {
      LALFree(xHat[i]->fft->data->data);
      LALFree(xHat[i]->fft->data);
      LALFree(xHat[i]->fft);
      LALFree(xHat[i]->fA->data->data);
      LALFree(xHat[i]->fA->data);
      LALFree(xHat[i]->fA);
      LALFree(xHat[i]->fB->data->data);
      LALFree(xHat[i]->fB->data);
      LALFree(xHat[i]->fB);
      LALFree(xHat[i]);
    }
  LALFree(xHat);

  /* Deallocate template and signal params */
  LALFree(templateParams->spind->spParams);
  LALFree(templateParams->spind);
  LALFree(templateParams->skyP);
  LALFree(templateParams);

  LALFree(signalParams->spind->spParams);
  LALFree(signalParams->spind);
  LALFree(signalParams->skyP);
  LALFree(signalParams);

  /* Deallocate demodulation params */
  LALFree(demParams->skyConst);
  LALFree(demParams->spinDwn);
  LALFree(demParams);

  LALFree(basicInputsFile);
  /* Anything else */

  /* Deallocate timestamps */
  LALFree(timeStamps);
  LALFree(sftPerCoh);

  LALFree(edat->ephemS);
  LALFree(edat->ephemE);
  LALFree(edat);
  LALCheckMemoryLeaks();
  return 0;
}


/***** This is the routine which computes the timestamps *****/
static void times(REAL8 tSFT, INT4 howMany, LIGOTimeGPS *ts, INT4 sw)
{
  int i=0, j=0;
  int temp1=0, temp2=0;

  while(i<sw*howMany)
    {
      temp1=floor(tSFT*(double)i);
      temp2=(int)((tSFT*(double)i-temp1)*1E9);
      /* This is Jan 1 1998 + 30 days, roughly */
      ts[j].gpsSeconds=temp1+567648000+86400*30;
      ts[j].gpsNanoSeconds=temp2;
      i=i+sw;
      j++;
    }
}

/* This routine can generate timestamps with random dropouts.
   Memory is allocated inside this routine, and it must be
   freed by the user at the end of its usefulness.  The allocated
   memory is for the timeStamps and sftPerCoh vectors.
*/
/* THIS IS HOW THIS ROUTINE WORKS */
/* This routine can generate timestamps with random dropouts.
   This is done by applying the '-d' switch on the command line.  The t.s. are
   made inside two nested while loops.  The outer loop loops over the order
   of DeFT to create; the inner loop, the SFTs which create a DeFT.  Thus,
   in order to create gaps, the random number generator decides inside the
   inner loop whether to create a gap or not.  If it decides yes (create
   a gap) the index over the SFTs is incremented, but no timestamp is
   computed.  If no, the timestamp is computed, and two counters are
   incremented to signify this fact.  One counter counts the number of
   timestamps that have been made for this tCoh, and inserts this number
   into the output array sftPerCoh.  This array contains the value of alpha
   that each demodulation needs to begin with.  The other counter counts the
   number of SFTs generated.
*/
static void times2(REAL8 tSFT, INT4 mObsCoh, LIGOTimeGPS **ts, INT4 **sftPerCoh, INT4 sw, INT4 mCohSFT)
{
  int i=0, j=0, k=0, m=0;
  int temp2=0;
  double temp1=0;

  LIGOTimeGPS *tempTS;
  INT4 *tempPC;

  tempTS=(LIGOTimeGPS *)LALCalloc((mObsCoh+1)*mCohSFT, sizeof(LIGOTimeGPS));
  tempPC=(INT4 *)LALCalloc(mObsCoh+1, sizeof(INT4));
  /* this is because first alpha of demod needs 0 index */
  tempPC[0]=0;

  if(sw!=0){
    int seed;

    seed=getpid();
    srand(seed);

    while(i<mObsCoh)
      {
	k=0;
	while(k<mCohSFT)
	  {
	    if( (rand()/(RAND_MAX+1.0)) <= 0.70 )
	      {
		int temp=i*mCohSFT+k;
		temp1=floor(tSFT*(double)temp);
		temp2=(int)((tSFT*(double)temp-temp1)*1.0E9);
		/* This is Jan 1 1998 + 30 days, roughly */
		tempTS[j].gpsSeconds=temp1+567648000+86400*30;
		tempTS[j].gpsNanoSeconds=temp2;
		j++;m++;
	      }
	    k++;
	  }
	if( m == 0 ) k--;
	else
	  {
	    tempPC[i+1] = m;
	    i++;
	  }
      }

    /* Now write to output arrays */
    *ts=(LIGOTimeGPS *)LALCalloc(j, sizeof(LIGOTimeGPS));
    *sftPerCoh=(INT4 *)LALCalloc(mObsCoh+1, sizeof(INT4));
    i=0;
    while(i<mObsCoh+1)
      {
	(*sftPerCoh)[i] = tempPC[i];
	i++;
      }
    i=0;
    while(i<j)
      {
	(*ts)[i].gpsSeconds = tempTS[i].gpsSeconds;
	(*ts)[i].gpsNanoSeconds = tempTS[i].gpsNanoSeconds;
	i++;
      }
  }

  else
    {
      *ts=(LIGOTimeGPS *)LALCalloc(mObsCoh*mCohSFT, sizeof(LIGOTimeGPS));
      *sftPerCoh=(INT4 *)LALCalloc(mObsCoh+1, sizeof(INT4));
      (*sftPerCoh)[0]=0;
      while(i<mObsCoh)
	{
	  j=0;
	  (*sftPerCoh)[i+1] += (i+1)*mCohSFT;
	  while(j<mCohSFT)
	    {
	      int x = i*mCohSFT+j;
	      temp1=(int)(floor(tSFT*(double)x));
	      temp2=(int)((tSFT*(double)x-temp1)*1.E9);
	      /* This is Jan 1 1998 + 30 days, roughly */
	      (*ts)[x].gpsSeconds = temp1+567648000+86400*30;
	      (*ts)[x].gpsNanoSeconds = temp2;
	      j++;
	    }
	  i++;
	}
    }

    /* Free up local memory */
    LALFree(tempTS);
    LALFree(tempPC);
}


static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps)
{
  INT4 x, y;

  x=tgps->gpsSeconds;
  y=tgps->gpsNanoSeconds;
  *f=(REAL8)x+(REAL8)y*1.e-9;
}


static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;

  temp0 = floor(*f);     /* this is tgps.S */
  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2);
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
}





