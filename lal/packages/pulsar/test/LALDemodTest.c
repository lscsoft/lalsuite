/* <lalVerbatim file="LALDemodTestCV">
Author: Berukoff, S.J., Papa, M.A., $Id$
 </lalVerbatim> */

#if 0
<lalLaTeX>

\subsection{Program \texttt{LALDemodTest.c}}
\label{s:LALDemodTest.c}

Performs required tests of \verb@LALDemod()@.

\subsubsection*{Usage}
\begin{verbatim}
LALDemodTest -i <input data file> [-g <gap>] [-m] [-n] [-o] 
\end{verbatim}

\subsubsection*{Description}

\noindent This routine performs tests on the routine \verb@LALDemod()@ as per
LAL standards.  It is a lengthy code, in comparison with the demodulation
routine itself.  In brief, this is an overview of the options: an input data
file may be specified on the command line with '\verb@-i@', and contains parameters relevant to the search (search band, timescales, spindown parameters, sky positions, etc.).  A fake signal is then produced, which can be superimposed on zero-mean, Gaussian noise using the '\verb@-n@' switch.  This signal can be modulated (by adding the '-m' switch on the command line), and one can simulate 'gaps' in the data by adding the '\verb@-g <gap>@' option, where \verb@<gap>@ refers to the integral number SFT timescales between adjacent timestamps.  Finally, to print output files containing the DeFT data, one uses the '\verb@-o@' option.  Please note that an example input file, \verb@in.data@, is included with this distribution.  This file contains the default search values, if none are specified in a different data file.

In more detail, let us begin with a discussion of the structure of the test
code, which is composed of several modules.  
\begin{itemize}
\item The first module reads in data from an input parameter data file.  The parameters must be listed in the input file in the following order, with the corresponding format:
\begin{verbatim}
total observation time -- float
coherent search time -- float
factor by which to modify SFT timescale -- float
amplitude SNR -- float
DeFT frequency band (centered by default around f0) -- float
f0, intrinsic frequency of the signal at the beginning of the observation -- float
maximum order of signal spindown parameters -- int
signal spindown parameter 1 --  scientific notation (value in Hz^2)
signal spindown parameter 2 --  scientific notation (value in Hz^3)
signal spindown parameter 3 --  scientific notation (value in Hz^4)
signal spindown parameter 4 --  scientific notation (value in Hz^5)
signal spindown parameter 5 --  scientific notation (value in Hz^6)
signal source right ascension (alpha) -- float (value in DEGREES)
signal source declination (delta) -- float (value in DEGREES)
maximum order of template spindown parameters -- int
template spindown parameter 1 --  scientific notation (value in Hz^2)
template spindown parameter 2 --  scientific notation (value in Hz^3)
template spindown parameter 3 --  scientific notation (value in Hz^4)
template spindown parameter 4 --  scientific notation (value in Hz^5)
template spindown parameter 5 --  scientific notation (value in Hz^6)
template source right ascension (alpha) -- float (value in DEGREES)
template source declination (delta) -- float (value in DEGREES)
\end{verbatim}


\item The next module in this test code creates noise using LALs
\verb@LALNormalDeviates()@ routine.  By design, the noise is created in single
precision, and is zero-mean and Gaussian.  This noise is added, datum-by-datum, to the time series created in
the next module, after the amplitude of the time series has been changed by a
factor of \verb@SNR@, which is specified in the input data file.

\item The next module to be invoked creates a time series, which undergoes a FFT; this transformed data constitutes the SFT data to be input to the demodulation code. This is the fake signal data set that we need to generate in order to test the demodulation code.  The fake signal data is characterized by an intrinsic frequency at the beginning of the observation plus some other source parameters.  The DeFT is produced in a band \verb@f0Band@ centered at this frequency.  The band, as explained above, is an input parameter. The width of this band (plus some extra width of $2\cdot 10^{-4}f0 $ Hz) determines the
sampling frequency of the time series (Nyquist theorem).  In practice this
would be the inverse FFT of a data set that has been band-passed around
\verb@f0@ and then appropriately down-sampled (e.g. with a lock-in).  The
normalization rule for FFT data is the following: if sinusoidal data over a
time $T$ and with amplitude $A$ is FFT-ed, the sum of the square amplitude of
the output of the FFT (power) is equal to ${A^2 T}$.  Thus, the power peak at
the sinusoids frequency should be expected to be $\sim$ $\frac{A^2}{2} T$,
within a factor of 2.  The same normalization rule applies to the DeFT data.
Thus by piecing together $N$ SFTs we expect a DeFT power peak $\sim N$ higher
than that of the SFTs - in the case of perfect signal-template match.


Let us now spend a few words on the choice of the SFT time baseline.  Given an
intrinsic search frequency one can compute the longest time baseline which is
still compatible with the requirement that the instantaneous signal frequency
during such time baseline does not shift by more than a frequency bin.  This
is the default choice for the SFT length, having assumed that the modulation
is due to the spin of the Earth and having taken a simple epicyclic model to
evaluate the magnitude of this effect. 


It is possible to choose a different time baseline by specifying a value for
the variable \verb@gap@ other than 1.  Note that the SFT time baseline is
approximated to the nearest value such that the number of SFT samples is a
power of two.  This is also well documented in the code.  


The set of SFTs does not necessarily come from contiguous data sets: a set of
time stamps is created that defines the time of the first sample of each SFT
data chunk.  The timestamps which are required in many parts of the code are
generated in a small subroutine \verb@times()@.  This routine take as input
the SFT timescale \verb@tSFT@, the number of SFTs which will be created,
\verb@mObsSFT@, and a switch which lets the code know whether to make even
timestamps, or timestamps with gaps (see below for more on this).  The
subroutine then writes the times to the \verb@LIGOTimeGPS@ vector containing
the timestamps for the entire test code, and returns this vector.  Note that
each datum of the  \verb@LIGOTimeGPS@ vector is comprised of two fields; if
accessing the $i^{th}$ datum, the seconds part of the timestamp vector
\verb@ts@ is \verb@ts[i].gpsSeconds@ and the nanoseconds part is
\verb@ts[i].gpsNanoSeconds@.  These are the fields which are written in this
\verb@times()@.


As an important side note, let us discuss the effect that a vector of
timestamps with gaps has on the resulting transformed data.  Since each of the
timestamps refers to the first datum of each SFT, the existence of the gaps
means that instead of transforming a continuous set of data, we are reduced to
transforming a piecewise continuous set.  Since we can envision these gaps as
simply replacing real data with zeros, we correspondingly should see a power
loss in the resulting FFTs signal bins and a broadening of the power spectrum.
Since real detectors will clearly have gaps in the data, this effect is
obviously something we seek to minimize or eliminate if possible.  This work
continues to be under development.


The total observation time determines how many SFTs and how many DeFTs are
created.  The actual time baseline for both the DeFTs and the total
observation time might differ from the ones defined in the input file, the
reason being that they are rounded to the nearest multiple of the SFT time
baseline.

Note that use is made of an internal function named \verb@tdb()@
written by C. Cutler.  This function is a basic routine which simulates the
modulation effects due to the Earths motion.  Currently it uses an epicyclic
motion model for the Earth spin and the Earths orbital motion.  In time it
will be replaced by a routine submitted to LAL by C. Cutler. This routine will
compute at any given time the actual instantaneous position and velocity of a
detector at any specified location of the Earth with respect to the SSB.  For now, the user of the test routine may supply his or her own routine to accomplish the same task, provided that the function prototype matches that of \verb@tdb()@.

\item Following the creation of a short chunk of time series data, an FFT is
performed with the internal FFTW routines.  This outputs a frequency domain
chunk which is placed into the \verb@SFTData@ array of structures.  This will
contain all of the SFT data we need to demodulate, and in the future, will be
the storage area for the real data.

\item The next module begins the demodulation process.  First, the parameters
for the demodulation routine are assigned from values previously calculated in
the test code.  Similarly, parameters for the \verb@ComputeSky()@ routine are
assigned.  This routine computes the coefficients $A_{s\alpha}$ and
$B_{s\alpha}$ (see section \ref{s:ComputeSky.h}) of the spindown parameters
for the phase model weve assumed.  These coefficients are used within the
\verb@LALDemod()@ routine itself. Since they only depend on the template sky
position, in a search over many different spin-down parameters they are
reused, thus one needs compute them only once.  Finally, at last, the
demodulation routine itself is called, and, if the command line option
'\verb@-o@' is used,  output are several data files containing demodulated
data (these are by default named '\verb@xhat_#@').  The files have three
columns: one for the frequency, and two for the real and imaginary parts of
the data.

\item There is an option within the code which allows one to print out SFT
data, information about the peaks of the SFTs, and information about the
DeFTs' peak position (whether they various DeFTs' maximal power peaks are in
the same bins as all the others or not).  This option is available to
those interested by uncommenting some code in this test suite.  The comments
look like \begin{verbatim}
/*****
*****
	
	Code and/or commentary

*****
*****/
\end{verbatim}
Simply remove the necessary characters.  Note that in some places, the
commentary should be commented out (so the code can compile!)  These snippets
of code provide someone with the ability to see what the intermediate results
of the test code look like.
\end{itemize}


\subsubsection*{Exit codes}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()
LALOpenDataFile()
LALFclose()
LALSCreateVector()
LALCreateRandomParams()
LALNormalDeviates()
LALDestroyRandomParams()
LALSDestroyVector()
LALCCreateVector()
LALCreateForwardComplexFFTPlan()
LALCOMPLEX8VectorFFT()
LALCDestroyVector()
LALDestroyComplexFFTPlan()
ComputeSky()
LALFree()
LALDemod()
tdb() 
\end{verbatim}

\subsubsection*{Notes}
The implementation of the code here is intended to give a general outline of what the demodulation code needs to work.  Most of this test function performs steps (e.g., noise, time- and frequency-series generation) that will be already present in the data.  

\vfill{\footnotesize\input{LALDemodTestCV}}

</lalLaTeX>
#endif /* autodoc block */
 
#ifndef LALDEMODTEST_C
#define LALDEMODTEST_C
#endif

/* Usage */
#define USAGE "Usage: %s [-i basicInputsFile] [-n] [-d] [-o]\n"


/* Error macro, taken from ResampleTest.c in LAL's pulsar/test package */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,		\
		 __LINE__, LALDEMODTESTC, statement ? statement : "",\
		 (msg) );                                            \
}                                                                    \
else (void)(0)

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALConstants.h>
#include <lal/LALDemod.h>

NRCSID(LALDEMODTESTC, "$Id$");

void tdb(REAL8 alpha, REAL8 delta, REAL8 t_GPS, REAL8 *T, REAL8 *Tdot);

int lalDebugLevel =3;

int main(int argc, char **argv)
{
	static LALStatus status;
	
/***** VARIABLE DECLARATION *****/

	ParameterSet *signalParams;
	ParameterSet *templateParams;
	CHAR *basicInputsFile = NULL;
	FILE *bif;
	REAL8 tObs, tCoh, tSFT, fSample;
	REAL8 SNR;
	REAL8 f0;
	
	INT4 mCohSFT, mObsCoh, mObsSFT;
	REAL8 dfSFT, dtEFF;
	INT4 if0Min, if0Max, ifMin, ifMax;
	REAL8 f0Min, f0Max, fMin, f0Band;
	INT4 nDeltaF, n;

	LIGOTimeGPS *timeStamps;
	
	REAL4Vector *temp = NULL;
	REAL4Vector *noise = NULL;
	static RandomParams *params;
	INT4 seed=0;
	
	REAL8 baseTbary0, tDot, t, tBary;
	REAL8 dTbary, dTbary2, dTbary3, dTbary4, dTbary5, dTbary6;
	INT4 i, k;
	COMPLEX8Vector *tvec = NULL;
	
	FFT **SFTData;
	ComplexFFTPlan *pfwd = NULL;
	COMPLEX8Vector *fvec = NULL;
	
	DemodPar *demParams;
	CSParams *csParams;
	INT4 iSkyCoh;
	
	FFT **xHat;
	
	FILE *fp;
	REAL8 ts=0.0, ts0=0.0;
	REAL8 tn=0.0;
	
	const CHAR *modulation=NULL;
	const CHAR *noi=NULL;
	INT4 deletions=1;
	const CHAR *output=NULL;
	
	CHAR filename[13];
	REAL8 factor;
	INT2 arg;
	void (*funcName)(REAL8, REAL8, REAL8, REAL8 *, REAL8 *);
	
/* Comment out the following when printing the 'extra' output */
/*****
*****/
	
	/* INT4 *maxArray; 	*/
	/* REAL8 pw, pwMax; 	*/
	/* INT4 ipwMax; 		*/
	/* INT2 tmp=0; 		*/
	/* FILE *fp1; 		*/
	/* INT4 j; 			*/
	/* REAL8 dtSFT, dfCoh 	*/
	/* REAL8 fMax; 		*/
	/* INT4 nSFT; 		*/
/******
*******/
	
/***** END VARIABLE DECLARATION *****/

	
/***** PARSE COMMAND LINE OPTIONS *****/	
	basicInputsFile=(CHAR *)LALMalloc(50*sizeof(CHAR));
	sprintf(basicInputsFile, "in.data");
	arg=1;
	while(arg<argc) {
		
		/* the input file */
		if(!strcmp(argv[arg],"-i"))
		{
				basicInputsFile=argv[++arg];
				arg++;
			if(LALOpenDataFile(basicInputsFile)==NULL)
			{
				ERROR(LALDEMODH_ENOFILE, LALDEMODH_MSGENOFILE, 0);
				LALPrintError(USAGE, *argv);
				return LALDEMODH_ENOFILE;
			}
		}
		
		/* turn noise on? if string is not NULL, it will be */
		else if(!strcmp(argv[arg],"-n"))
		{
			noi="n";
			arg++;
		}
		
		/* turn timestamp deletions on? if string is not NULL, it will be */
		else if(!strcmp(argv[arg],"-g"))
		{
			deletions=atoi(argv[++arg]);
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
			LALPrintError(USAGE, *argv);
			arg=argc;
			return LALDEMODH_EBADARG;
		}
		
	}
	
/***** END COMMAND LINE PARSING *****/

	
/***** SET NAME OF TIMING ROUTINE *****/

/* Here you should set the name of your timing routine to 		*/
/* whatever it is.  Of course you'll also need to compile it in! 	*/
funcName=tdb;

/*** END TIMING ROUTINE NAME ***/
	
	
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
	
	fscanf(bif, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n%le\n%le\n%le\n%le\n%le\n%lf\n%lf\n%d\n%le\n%le\n%le\n%le\n%le\n%lf\n%lf\n",
	&tObs, &tCoh, &factor, &SNR, &f0Band, &f0,
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

/*	2^n gives the number of samples of the SFT. This will be a signal in a 	*/	/*	band fSample/2 with fSample=2*(2.e-4 f0+f0Band). This makes sure that	*/ /*	that the sampling frequency enables us to produce the DeFT frequency 	*/ /*	band that we want (f0Band) and that there's enought wings of SFT data	*/ /*	(2.e-4*f0) to produce such DeFT band. 						*/

	fSample=2.0*(2.e-4*f0+f0Band);

/* 	The main formula used here is that the gives the maximum allowed		*/ /*	time-baseline (Tmax) such that a signal of frequency f0 Doppler 		*/ /*	modulated due to Earth spin is seen as monochromatic during such		*/ /*	observation time: 									*/
/* 	  	Tmax < 9.6e4/sqrt(f0). 								*/
/*	We have thus taken Tmax=9.4e4/sqrt(f0). A different criterion can be	*/ /*	chosen by setting the variable "factor" to a suitable value. We have	*/ /*	rounded the resulting number of samples to the nearest smallest power 	*/ /*	of two. 											*/

 	n=floor((log(factor*94000.0*fSample/sqrt(f0))/log(2.0))); 
	/* Number of SFT time series data points */
	nDeltaF=ldexp(1., n);
	tSFT=(REAL8)nDeltaF/fSample;
	dfSFT=1/tSFT;
	/* Size of time bin */
	dtEFF=1.0/fSample;
	ifMin=floor(f0/dfSFT)-nDeltaF/4.0;
	ifMax=ifMin+nDeltaF/2.0;
	fMin=dfSFT*ifMin;
	/* Maximum search band frequency, in Hz */
	if0Max=ifMax-ceil((f0Band/2.0)/dfSFT);
	f0Max=dfSFT*if0Max;
	/* Minimum search band frequency, in Hz*/
	if0Min=ifMin+floor((f0Band/2.0)/dfSFT);
	f0Min=dfSFT*if0Min;
	/* DeFT band, in Hz */
	f0Band=f0Max-f0Min;

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

	/* convert input angles from degrees to radians */
	signalParams->skyP->alpha=signalParams->skyP->alpha*LAL_TWOPI/360.0;
	signalParams->skyP->delta=signalParams->skyP->delta*LAL_TWOPI/360.0;
	templateParams->skyP->alpha=templateParams->skyP->alpha*LAL_TWOPI/360.0;
	templateParams->skyP->delta=templateParams->skyP->delta*LAL_TWOPI/360.0;

/***** END USEFUL QUANTITIES *****/	


/***** CALL ROUTINE TO GENERATE TIMESTAMPS *****/
	
	timeStamps=(LIGOTimeGPS *)LALMalloc(mObsSFT*sizeof(LIGOTimeGPS));
	times(tSFT, mObsSFT, timeStamps, deletions);
	
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


/***** CREATE SIGNAL *****/

	/* Create vector to hold frequency series */
	LALCCreateVector(&status, &fvec, (UINT4)nDeltaF);
	
	/* Compute plan for FFTW */
	LALCreateForwardComplexFFTPlan(&status, &pfwd, (UINT4)nDeltaF, 1);
	
	/* Create vector to hold time series */
	LALCCreateVector(&status, &tvec, (UINT4)nDeltaF);
	
	/* Allocate memory for the SFTData structure */
	SFTData=(FFT **)LALMalloc(mObsSFT*sizeof(FFT *));
	for(i=0;i<mObsSFT;i++)
	{
		SFTData[i]=(FFT*)LALMalloc(sizeof(FFT));
		SFTData[i]->fft=(COMPLEX16FrequencySeries 
			*)LALMalloc(sizeof(COMPLEX16FrequencySeries));
		SFTData[i]->fft->data=(COMPLEX16Vector 
			*)LALMalloc(sizeof(COMPLEX16Vector));
		SFTData[i]->fft->data->data=(COMPLEX16 
			*)LALMalloc(nDeltaF*sizeof(COMPLEX16));
	}
	
	/* conversion: detector->barycenter time,  t=GPS_0 */ 
	ts0=(REAL8)(timeStamps[0].gpsSeconds)+	 
		(REAL8)(timeStamps[0].gpsNanoSeconds)*1.0E-9;
	(*funcName)((REAL8)signalParams->skyP->alpha,(REAL8)signalParams->skyP->delta,	 
		ts0, &baseTbary0, &tDot);

	/*****
	*****
		
		Below I've commented out file pointers.  These open the file descriptors for files which can contain output data; one is 'sft.data', which 	contains all of the data from all of the SFTs, 'peaks.data', which contains the peaks of the SFTs.  These are used in conjunction with areas commented out below.
	
	fp=LALFopen("sft.data","w");
	fp1=LALFopen("peaks.data","w");	
	maxArray=LALMalloc(mObsSFT*sizeof(INT4));	
	
	*****
	*****/
	
	for(k=0; k<mObsSFT; k++) 
	{
		ts=(REAL8)(timeStamps[k].gpsSeconds)*1.00;
		tn=(REAL8)(timeStamps[k].gpsNanoSeconds)*1.00E-9;
		
		for(i=0;i<nDeltaF;i++)
		{
			t=ts+tn+(REAL8)i*dtEFF;
			/* conversion again */
			(*funcName)((REAL8)signalParams->skyP->alpha,
				(REAL8)signalParams->skyP->delta, (REAL8)t, &tBary,  							&tDot);
			
			/* pipeline for better performance */
			dTbary=tBary-baseTbary0;
			dTbary2=dTbary*dTbary;
			dTbary3=dTbary2*dTbary;
			dTbary4=dTbary3*dTbary;
			dTbary5=dTbary4*dTbary;
			dTbary6=dTbary5*dTbary;
			
			/* create time series according to phase model */
			
			if(noi!=NULL) /* first, if noise is selected */
			{
			tvec->data[i].re=SNR*cos(LAL_TWOPI*((REAL8)f0*dTbary-
				(REAL8)fMin*(t-ts0)+                            
				0.5*signalParams->spind->spParams[0]*dTbary2+   
				1/3.0*signalParams->spind->spParams[1]*dTbary3+ 
				0.25*signalParams->spind->spParams[2]*dTbary4+  
				0.20*signalParams->spind->spParams[3]*dTbary5+
				1/6.0*signalParams->spind->spParams[4]*dTbary6))
				+noise->data[k*nDeltaF+i]; 
			tvec->data[i].im=0.0;
			}
			else /* if no noise is selected */
			{
			tvec->data[i].re=SNR*cos(LAL_TWOPI*((REAL8)f0*dTbary-
				(REAL8)fMin*(t-ts0)+                            
				0.5*signalParams->spind->spParams[0]*dTbary2+   
				1/3.0*signalParams->spind->spParams[1]*dTbary3+ 
				0.25*signalParams->spind->spParams[2]*dTbary4+  
				0.20*signalParams->spind->spParams[3]*dTbary5+
				1/6.0*signalParams->spind->spParams[4]*dTbary6)); 
			tvec->data[i].im=0.0;
	

			}
		}
		
		/* Perform FFTW-LAL Fast Fourier Transform */
		LALCOMPLEX8VectorFFT(&status, fvec, tvec, pfwd); 
		
		/* write the SIGNAL+NOISE to the SFTData structure and normalize */
		for(i=0;i<(int)fvec->length;i++)
		{
			SFTData[k]->fft->data->data[i].re=fvec->data[i].re*dtEFF/sqrt(tSFT);
			SFTData[k]->fft->data->data[i].im=fvec->data[i].im*dtEFF/sqrt(tSFT);
		}

       	/***** Uncomment this to enable file output of SFTs and peaks. *****/
		/*****
		
			pwMax=0;
	
		for(i=0;i<fvec->length;i++)
		{	pw=SFTData[k]->fft->data->data[i].re*SFTData[k]->fft->data->data[i].re+SFTData[k]->fft->data->data[i].im*SFTData[k]->fft->data->data[i].im;

			if(pwMax<pw)
			{
				pwMax=pw;
				ipwMax=i;
			}
		}
 
		fprintf(fp1,"%d\t%d\n",ipwMax,k);
		
		for(i=0;i<fvec->length/2;i++)
		{
			fprintf(fp,"%lf\t%lf\n",SFTData[k]->fft->data->data[i].re*SFTData[k]->fft->data->data[i].re+SFTData[k]->fft->data->data[i].im*SFTData[k]->fft->data->data[i].im,(ifMin+i)*dfSFT); 
		}	
		
		*****
		*****/
		
		/* assign particulars to each SFT */
		SFTData[k]->fft->epoch=timeStamps[k];
		SFTData[k]->fft->f0=fMin;
		SFTData[k]->fft->deltaF=dfSFT;
		SFTData[k]->fft->data->length=nDeltaF;
		
	} 
	
	/***** Uncomment this too.
	*****
	
	LALFree(maxArray);
	fclose(fp);
	fclose(fp1);
	
	*****
	*****/
	
	/* Free entities used here */
	if(noi!=NULL) {LALDestroyVector(&status, &noise);} 
	LALCDestroyVector(&status, &fvec);
	LALCDestroyVector(&status, &tvec);
	LALDestroyComplexFFTPlan(&status, &pfwd);
	
/***** END CREATE SIGNAL *****/
	
	
/***** DEMODULATE SIGNAL *****/

	/* Allocate space and set quantity values for demodulation parameters */
	demParams=(DemodPar *)LALMalloc(sizeof(DemodPar));
	demParams->skyConst=(REAL8 *)LALMalloc((2*templateParams->spind->m *				(mObsSFT+1)+2*mObsSFT+3)*sizeof(REAL8));
	demParams->spinDwn=(REAL8 *)LALMalloc(templateParams->spind->m*sizeof(REAL8));
	demParams->if0Max=if0Max;
	demParams->if0Min=if0Min;
	demParams->mCohSFT=mCohSFT;
	demParams->mObsCoh=mObsCoh;
	demParams->ifMin=ifMin;
	demParams->spinDwnOrder=templateParams->spind->m;
		
	for(i=0;i<signalParams->spind->m;i++)
	{
		demParams->spinDwn[i]=templateParams->spind->spParams[i];
	}

	/* Allocate space and set quantities for call to ComputeSky() */
	csParams=(CSParams *)LALMalloc(sizeof(CSParams));
	csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
	csParams->skyPos[0]=templateParams->skyP->alpha;
	csParams->skyPos[1]=templateParams->skyP->delta;
	csParams->tGPS=timeStamps;
	csParams->spinDwnOrder=templateParams->spind->m;
	csParams->mObsSFT=mObsSFT;
	csParams->tSFT=tSFT;
	csParams->funcName=funcName;
	iSkyCoh=0;

	/* Call COMPUTESKY() */
	ComputeSky(&status, demParams->skyConst, iSkyCoh, csParams);

	/* Deallocate space for ComputeSky parameters */
	LALFree(csParams->skyPos);
	LALFree(csParams);

	/* Allocate memory for demodulated data */
	xHat=(FFT **)LALMalloc(mObsCoh*sizeof(FFT *));
	for(i=0; i<mObsCoh; i++)
	{
		xHat[i]=(FFT *)LALMalloc(sizeof(FFT));
		xHat[i]->fft=(COMPLEX16FrequencySeries *) 	
		LALMalloc(sizeof(COMPLEX16FrequencySeries));
		xHat[i]->fft->data=(COMPLEX16Vector *)LALMalloc(sizeof(COMPLEX16Vector));
	xHat[i]->fft->data->data=(COMPLEX16*) LALMalloc(((if0Max-if0Min)*mCohSFT)*				sizeof(COMPLEX16));
	}
	
	/***** Uncomment this
	*****
	
		maxArray=LALMalloc(mObsCoh*sizeof(INT4));
	
	*****
	*****/
	
	for(k=0; k<mObsCoh; k++)
	{
		demParams->iCoh=k;
		
		/**************************/
		/*       DEMODULATE       */
		/**************************/

		LALDemod(&status, xHat, SFTData, demParams);
		if(output!=NULL)
		{
			sprintf(filename,"xhat_%d.data",k);
			fp=LALFopen(filename,"w");
			printf("Dumping demodulated data to disk: xhat_%d.data  \n",k);
			for(i=0;i<(if0Max-if0Min)*mCohSFT;i++) 
			{
			fprintf(fp,"%24.16f\t%24.16f\t%24.16f\n",f0Min+(REAL8)i/tCoh, 	
				xHat[k]->fft->data->data[i].re, xHat[k]->fft->data->data[i].im); 			}
			LALFclose(fp);
		}
	
		/***** Uncomment this to enable printout of the position of maximal 				*****		power DeFT frequency bins.
		
		pwMax=0;
		for(i=0;i<(if0Max-if0Min)*mCohSFT;i++)
		{
			pw=xHat[k]->fft->data->data[i].re*xHat[k]->fft->data->data[i].re+
				xHat[k]->fft->data->data[i].im*xHat[k]->fft->data->data[i].im;
			if(pwMax<pw)
			{
				pwMax=pw;
				ipwMax=i;
			}
		}

		maxArray[k]=ipwMax; */
	}
	
	/*

	for(i=0; i<mObsCoh; i++)
	{
		for(k=mObsCoh-1; k>=0; k--)
		{
			if(maxArray[i]!=maxArray[k]) 
				{
				printf("Max bin of DeFT %d (%d) is different from DeFT %d (%d). \n",i ,maxArray[i], k, maxArray[k]);
				tmp=1;
				}	
		}
	}
	if(tmp==0) printf("All max power DeFT bins are identical.\n");

		*****
		*****/
		
/***** END DEMODULATION *****/
			
		
/***** DEALLOCATION *****/
	
	/* Deallocate SFTData structure, since we don't need it anymore */
	for(i=0;i<mObsSFT;i++)
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
	
	LALFree(timeStamps);
	LALFree(basicInputsFile);
	/* Anything else */
	
	/***** 
	*****
	
	LALFree(maxArray);
	
	*****
	*****/
	
	LALCheckMemoryLeaks(); 
        return 0;
}


/* internal routine to convert detector time to barycentric time */
/* Note: this routine will be replaced by Barycenter(), written  */
/* by C. Cutler, Albert-Einstein-Institut.  			 */

void tdb(REAL8 alpha, REAL8 delta, REAL8 t_GPS, REAL8 *T, REAL8 *Tdot) 
{

	/*RA alpha and DEC delta have their usual meanings, but units 
  are radians, not hours or degrees */

	REAL8 tdiff;
	REAL8 tdiff_day;
	REAL8 tdiff_yr; 
	REAL8 beta;
	REAL8 x;
	REAL8 y; 
	REAL8 lam;
	REAL8 phi0;
	REAL8 psi0; 
	REAL8 phi;
	REAL8 psi;
	REAL8 cose;
	REAL8 sine;
	REAL8 sind; 
	REAL8 cosd;
	REAL8 sina;
	REAL8 cosa;
	REAL8 cosb;

	/*phi0 and psi0 given location and orientation, resp., of Earth at t_GPS =0 */
	REAL8 T_bary, Tdot_bary;
	/* T_bary (what we want to know) is pulse emission time in TDB coords.*/
	/* The lab measures GPS time, but we set t_GPS =0 at first instant of 1998*/
	REAL8 AU = 1.49597870e13/2.99792458e10;
	REAL8 RE = 6.37814e8/2.99792458e10;
	REAL8 eps = (23.0 + 26.0/60 + 21.488/3600)*LAL_PI/180;
	/*    printf("eps is %lf \n",eps); */
	phi0=2.86e0;
	psi0=4.12e0;
	/* Those values of phi0 and psi0 are completely made up */

	/* Now converting from equatorial to ecliptic coords. */
	sine=sin(eps);
	cose=cos(eps);
	sind=sin(delta);
	cosd=cos(delta);
	sina=sin(alpha);
	cosa=cos(alpha);
	beta=cose*sind-sine*cosd*sina;
	cosb=sqrt(1.0-beta*beta);
	x=cosd*cosa/cosb;
	y=(sine*sind+cose*cosd*sina)/cosb;
	lam=atan2(y,x);
	/*  printf("beta, lam are (rad) %lf  %lf \n", beta, lam);  */
	phi=phi0+2*LAL_PI*t_GPS/3.15576e7;
	psi=psi0+2*LAL_PI*t_GPS/8.64e4;
	tdiff_day=RE*cosd*cos(psi-alpha);
	tdiff_yr=AU*cosb*cos(phi-lam);
	tdiff=tdiff_day+tdiff_yr;
	T_bary=t_GPS+tdiff;
	Tdot_bary=1.e0-RE*cosd*sin(psi-alpha)*2*LAL_PI/8.64e4 		
		-AU*cosb*sin(phi-lam)*2*LAL_PI/3.15576e7;
	
	/* switch to turn mod off/on */
	/* note that 'off' negates the function of the code! */

	*T=T_bary; 
 	*Tdot=Tdot_bary;
		
	}

void times(REAL8 tSFT, int howMany, LIGOTimeGPS *ts, INT4 sw)
{
	int i=0, j=0;
	int temp1, temp2;
	FILE *fp;
	
	fp=fopen("times.blah","w");
	while(i<sw*howMany)
	{
		temp1=floor(tSFT*i);
		temp2=(int)((tSFT*(double)i-temp1)*1E9);
	
		if(i==0) {
			ts[j].gpsSeconds=0;
			ts[j].gpsNanoSeconds=0;
			}
		else {
			ts[j].gpsSeconds=temp1;
			ts[j].gpsNanoSeconds=temp2;
			}
		
		fprintf(fp,"%d %d\n",ts[j].gpsSeconds, ts[j].gpsNanoSeconds);
		i=i+sw;
		j++;
	}fclose(fp);
}

