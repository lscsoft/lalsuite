/*
 * Copyright (C) 2009 Chris Messenger, Pinkesh Patel
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

/** \author R. Prix, J. T. Whelan
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space,
 * following the equations in \ref JKS98.
 *
 * This code is partly a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/FindRoot.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include "ComputeFstat_RS.h"
#include <lal/ComplexAM.h>

NRCSID( COMPUTEFSTATRSC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


#define LD_SMALL4       (2.0e-4)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

/*----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define SQ(x) ( (x) * (x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* empty initializers  */
static const LALStatus empty_status;
static const AMCoeffs empty_AMCoeffs;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFParams empty_ComputeFParams;
const ComputeFBuffer empty_ComputeFBuffer;
const EarthState empty_EarthState;

static REAL8 p,q,r;          /* binary time delay coefficients (need to be global so that the LAL root finding procedure can see them) */

/*---------- internal prototypes ----------*/
int finite(double x);

/*==================== FUNCTION DEFINITIONS ====================*/


/** Function to compute a vector of Fstatistic values for a number of frequency bins.
    The output, i.e. fstatVector must be properly allocated
    before this function is called.  The values of the start frequency, the step size
    in the frequency and the number of frequency values for which the Fstatistic is
    to be calculated are read from fstatVector.  The other parameters are not checked and
    they must be correctly set outside this function.
*/ 
void ComputeFStatFreqBand_RS ( LALStatus *status,
			       REAL4FrequencySeries *fstatVector, 		/**< [out] Vector of Fstat values */
			       const PulsarDopplerParams *doppler,		/**< parameter-space point to compute F for */
			       const MultiSFTVector *multiSFTs, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
			       const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
			       const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
			       const ComputeFParams *params		/**< addition computational params */
			       )
{

  UINT4 numDetectors; 
  ComputeFBuffer *cfBuffer = NULL;
  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiCOMPLEX8TimeSeries *multiTimeseries = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv;
  SkyPosition skypos;
  UINT4 i,j,k;
  UINT4 nspins = PULSAR_MAX_SPINS - 1;
  REAL8 deltaref;
  COMPLEX8Vector *Faf = NULL;
  COMPLEX8Vector *Fbf = NULL;
  REAL8Vector *Fstat = NULL;
  UINT4 N;

  INITSTATUS( status, "ComputeFStatFreqBand_RS", COMPUTEFSTATRSC );
  ATTATCHSTATUSPTR (status);

  /* check that the input data and parameters structures don't point to NULL */
  ASSERT ( multiSFTs, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );

  numDetectors = multiSFTs->length;               /* set the number of detectors to the number of sets of SFTs */
  printf("numDetectors = %d\n",numDetectors);
  printf("alpha = %f delta = %f\n",doppler->Alpha,doppler->Delta);
  printf("doppler ref time = %f\n",XLALGPSGetREAL8(&(doppler->refTime)));

  /* check that the pre-allocated output vector doesn't point to NULL */ 
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  ASSERT ( fstatVector, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data->data, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data->length > 0, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );

  /* check that the multidetector states have the same length as the multiSFTs */
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  }

  /* check whether the multi-detector SSB times and detector states are already buffered */
  cfBuffer = params->buffer;
  if ( cfBuffer  
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta )
       && cfBuffer->multiSSB )
    { /* if already buffered then we reuse the values */
      multiSSB = cfBuffer->multiSSB;

      /* re-use (LWL) AM coefficients whenever available */
      if ( cfBuffer->multiAMcoef )
	multiAMcoef = cfBuffer->multiAMcoef;

    } 
  /* otherwise we need to compute the SSB times */
  else {
    skypos.system =   COORDINATESYSTEM_EQUATORIAL;
    skypos.longitude = doppler->Alpha;
    skypos.latitude  = doppler->Delta;

    TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, skypos, doppler->refTime, params->SSBprec ), status );
    
    if ( cfBuffer ) {
      XLALDestroyMultiSSBtimes ( cfBuffer->multiSSB );
      cfBuffer->multiSSB = multiSSB;
      cfBuffer->Alpha = doppler->Alpha;
      cfBuffer->Delta = doppler->Delta;
      cfBuffer->multiDetStates = multiDetStates ;
    } /* buffer new SSB times */
    
  } /* could not reuse previously buffered quantites */
  
  /* compute AM parameters if not buffered */
  if ( !multiAMcoef )
    {
      /* compute new AM-coefficients */
      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      for (i=0;i<10;i++) printf("a = %6.12f b = %6.12f\n",multiAMcoef->data[0]->a->data[i],multiAMcoef->data[0]->b->data[i]);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
      }

      for (i=0;i<10;i++) printf("noise weights = %6.12f\n",multiWeights->data[0]->data[i]);
      for (i=0;i<10;i++) printf("a = %6.12f b = %6.12f\n",multiAMcoef->data[0]->a->data[i],multiAMcoef->data[0]->b->data[i]);

      /* store these in buffer if available */
      if ( cfBuffer ) {
	XLALDestroyMultiAMCoeffs ( cfBuffer->multiAMcoef );
	cfBuffer->multiAMcoef = multiAMcoef;
      }

    } /* if LWL AM coefficient needed to be computed */

  /* store AM coefficient integrals in local variables */
  if ( multiAMcoef ) {
    Ad = multiAMcoef->Mmunu.Ad;
    Bd = multiAMcoef->Mmunu.Bd;
    Cd = multiAMcoef->Mmunu.Cd;
    Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
  }
  else {
    LALPrintError ( "Programming error: 'multiAMcoef' not available!\n");
    ABORT ( status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  }
  printf("hello\n");
  /* Generate weighted resampled timeseries */
  /* check whether the multi-detector resampled timeseries are already buffered */
  if ( cfBuffer && ( cfBuffer->multiTimeseries ) ) multiTimeseries = cfBuffer->multiTimeseries;
  else {
    TRY ( ResampleMultiSFTs ( status->statusPtr, &multiTimeseries, fstatVector->deltaF, multiAMcoef, multiSSB, multiSFTs ), status );
  }
  N = multiTimeseries->Fat[0]->length;
  printf("N = %d\n",N);
  
  /* determine number of spin down's */
  while (doppler->fkdot[nspins]==0.0) nspins--;
  printf("number of spin downs = %d\n",nspins); 
  
  /* compute the time differnce between timeseries epoch and reference time */
  deltaref = XLALGPSGetREAL8(&(multiTimeseries->epoch)) - XLALGPSGetREAL8(&(doppler->refTime));
  
  /* apply spin derivitive correction to resampled timeseries */
  /* loop over spin derivitives (nspins = 1 means first derivitive, = 2 means second derivitive etc.. )*/
  for (j=1;j<=nspins;j++) {
    
    /* loop over time samples */
    for (k=0;k<multiTimeseries->Fat[0]->length;k++) {
      
      REAL8 phase = (1.0/LAL_TWOPI)*doppler->fkdot[j]*pow(deltaref + k*multiTimeseries->deltaT,(REAL8)j);
      REAL8 cosphase = cos(phase);
      REAL8 sinphase = sin(phase);
      
      /* loop over detectors */
      for (i=0;i<numDetectors;i++) {
	
	REAL8 Fare = multiTimeseries->Fat[0]->data[i].re*cosphase - multiTimeseries->Fat[0]->data[i].im*sinphase;
	REAL8 Faim = multiTimeseries->Fat[0]->data[i].re*sinphase + multiTimeseries->Fat[0]->data[i].im*cosphase;
	REAL8 Fbre = multiTimeseries->Fbt[0]->data[i].re*cosphase - multiTimeseries->Fbt[0]->data[i].im*sinphase;
	REAL8 Fbim = multiTimeseries->Fbt[0]->data[i].re*sinphase + multiTimeseries->Fbt[0]->data[i].im*cosphase;
	
	multiTimeseries->Fat[0]->data[i].re = Fare;
	multiTimeseries->Fat[0]->data[i].im = Faim;
	multiTimeseries->Fbt[0]->data[i].re = Fbre;
	multiTimeseries->Fbt[0]->data[i].im = Fbim;
	
      }
    
    }
  
  }

  {
    /* TESTING */
    FILE *fp = NULL;
    fp = fopen("/home/chmess/test/out/timeseries.txt","w");
    REAL8 SSBstart = XLALGPSGetREAL8(&(multiTimeseries->epoch));
    for (i=0;i<N;i++) fprintf(fp,"%6.12f %6.12f %6.12f %6.12f %6.12f\n",
			      SSBstart + i*multiTimeseries->deltaT,
			      multiTimeseries->Fat[0]->data[i].re,multiTimeseries->Fat[0]->data[i].im,
			      multiTimeseries->Fbt[0]->data[i].re,multiTimeseries->Fbt[0]->data[i].im);
    fclose(fp);
  }

  /* allocate memory for Fa(f) and Fb(f) and initialise */
  Faf = XLALCreateCOMPLEX8Vector(N);
  Fbf = XLALCreateCOMPLEX8Vector(N);
  Fstat = XLALCreateREAL8Vector(N);
  memset(Faf->data,0,N*sizeof(COMPLEX8));
  memset(Fbf->data,0,N*sizeof(COMPLEX8));

  /* compute Fourier transforms to get Fa(f) and Fb(f) */
  /* loop over detectors */
  for (i=0;i<numDetectors;i++) {

    ComplexFFTPlan *pfwd = NULL;
    COMPLEX8Vector *ina = NULL;
    COMPLEX8Vector *inb = NULL;
    COMPLEX8Vector *outa = NULL;
    COMPLEX8Vector *outb = NULL;
    REAL8 startminusreftime = XLALGPSGetREAL8(&(multiTimeseries->epoch)) - XLALGPSGetREAL8(&(multiTimeseries->refTime));
    printf("startminusreftime = %6.12f\n",startminusreftime);
    printf("SSBstart = %d %d\n",multiTimeseries->epoch.gpsSeconds,multiTimeseries->epoch.gpsNanoSeconds);
    printf("SSBrefTime = %d %d\n",multiTimeseries->refTime.gpsSeconds,multiTimeseries->refTime.gpsNanoSeconds);

    /* point the input to the timeseries */
    ina = multiTimeseries->Fat[i];
    inb = multiTimeseries->Fbt[i];
    
    /* allocate memory for the outputs */
    outa = XLALCreateCOMPLEX8Vector(N);
    outb = XLALCreateCOMPLEX8Vector(N);
    {
      INT4 err = xlalErrno;
      if ( err != XLAL_SUCCESS ) {
	ABORT ( status, err, "XLALCreateCOMPLEX8Vector() failed!\n");
      }
    }

    /* make forwards FFT plan */
    pfwd = XLALCreateCOMPLEX8FFTPlan(N,1,0);

    /* Fourier transform Fa(t)  */
    if (XLALCOMPLEX8VectorFFT(outa,ina,pfwd)!= XLAL_SUCCESS) 
      ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");

    /* Fourier transform Fb(t)  */
    if (XLALCOMPLEX8VectorFFT(outb,inb,pfwd)!= XLAL_SUCCESS) 
      ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
    
    /* correct for reference time effects */
    /* loop over freq samples */
    for (k=0;k<outa->length;k++) {
      
      REAL8 phase = -LAL_TWOPI*(k*fstatVector->deltaF)*startminusreftime;
      REAL8 cosphase = cos(phase);
      REAL8 sinphase = sin(phase);
     /*  printf("f = %6.12f phase = %6.12f\n",k*fstatVector->deltaF,phase); */
      
      REAL8 Fare = outa->data[k].re*cosphase - outa->data[k].im*sinphase;
      REAL8 Faim = outa->data[k].re*sinphase + outa->data[k].im*cosphase;
      REAL8 Fbre = outb->data[k].re*cosphase - outb->data[k].im*sinphase;
      REAL8 Fbim = outb->data[k].re*sinphase + outb->data[k].im*cosphase;
      
      outa->data[k].re = Fare;
      outa->data[k].im = Faim;
      outb->data[k].re = Fbre;
      outb->data[k].im = Fbim;
      
    }
    
  
    {
      /* TESTING */
      FILE *fp = NULL;
      fp = fopen("/home/chmess/test/out/freqseries.txt","w");
      REAL8 f0 = multiTimeseries->f0;
      for (i=0;i<N;i++) fprintf(fp,"%6.12f %6.12f %6.12f %6.12f %6.12f\n",
				f0 + i*fstatVector->deltaF,
				outa->data[i].re,outa->data[i].im,
				outb->data[i].re,outb->data[i].im);
      fclose(fp);
    }
    
    /* add to summed Faf and Fbf */
    for (j=0;j<N;j++) {
      Faf->data[j].re += outa->data[j].re;
      Faf->data[j].im += outa->data[j].im;
      Fbf->data[j].re += outb->data[j].re;
      Fbf->data[j].im += outb->data[j].im;
    }
    
  }
  
  printf("length of fstatvector = %d\n",fstatVector->data->length);
  
  /* loop over frequency and construct 2F */
  for (k=0;k<N;k++) {
    
    /* ----- compute final Fstatistic-value ----- */
    
    /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
     * therefore there is a factor of 2 difference with respect to the equations in JKS, which
     * where based on the single-sided PSD.
     */
    Fstat->data[k] = Dd_inv * (  Bd * (SQ(Faf->data[k].re) + SQ(Faf->data[k].im) )
				 + Ad * ( SQ(Fbf->data[k].re) + SQ(Fbf->data[k].im) )
				 - 2.0 * Cd *( Faf->data[k].re * Fbf->data[k].re + Faf->data[k].im * Fbf->data[k].im )
				 );
    
  }
  
  {
    /* TESTING */
    FILE *fp = NULL;
    fp = fopen("/home/chmess/test/out/fstat.txt","w");
    REAL8 f0 = multiTimeseries->f0;
    for (i=0;i<N;i++) fprintf(fp,"%6.12f %6.12f\n",
			      f0 + i*fstatVector->deltaF,
			      Fstat->data[i]);
    fclose(fp);
  }

  /* free memory if no buffer was available */
  if ( !cfBuffer ) {
    XLALDestroyMultiSSBtimes ( multiSSB );
    XLALDestroyMultiAMCoeffs ( multiAMcoef );
  } /* if !cfBuffer */
  
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStatFreqBand_RS() */


/** Function to compute a resampled timeseries from a multidetector set of SFTs */
void ResampleMultiSFTs ( LALStatus *status, 
			 MultiCOMPLEX8TimeSeries **multitimeseries,      	/**< [out] Array of resampled timeseries */
			 REAL8 deltaF,                                          /**< [in] user requested frequency resolution (determines zero padding) */
			 const MultiAMCoeffs *multiAMcoef,                      /**< [in] the multidetector AM coefficients */
			 const MultiSSBtimes *multiSSB,                         /**< [in] the SSB times ascociated with the MID-POINT of each SFT */
			 const MultiSFTVector *multiSFTs                        /**< [in] the multi-detector SFT data */
			 ) 
{

  UINT4 i,j,k;
  REAL8 SSBstart,SSBend,DETstart,DETend;
  LIGOTimeGPS SSBstart_GPS,SSBend_GPS;
  REAL8 Tprime,T;
  UINT4 numDetectors = multiSFTs->length;
  UINT4 N;
  REAL8 deltaT;
  ComplexFFTPlan *prev = NULL;
  COMPLEX8Vector *out = NULL;
  REAL8 SSBrefTime = XLALGPSGetREAL8(&(multiSSB->data[0]->refTime));
  LIGOTimeGPS SSBrefTime_GPS = multiSSB->data[0]->refTime;

  /* define the SFT parameters */
  REAL8 fhet = multiSFTs->data[0]->data[0].f0;                      /* the start frequency of each SFT (we assume that this is the same for all SFTS) */ 
  REAL8 SFTdeltaF = multiSFTs->data[0]->data[0].deltaF;             /* the frequency spacing in each SFT (again we assume this is equal for all SFTs) */
  REAL8 tSFT = 1.0/SFTdeltaF;                                       /* the length of the SFT (try not to use this in case of rounding) */ 
  UINT4 SFTbins = multiSFTs->data[0]->data[0].data[0].length;       /* the number of complex frequency bins in each SFT (we assume that this is the same for all SFTS) */ 
  REAL8 deltaTSFT = 1.0/(SFTdeltaF*SFTbins);                        /* the time spacing of the inverse-FFTed SFT data */
  printf("fhet = %6.12f SFTdeltaF = %6.12f tSFT = %6.12f SFTbins = %d\n",fhet,SFTdeltaF,tSFT,SFTbins);

  INITSTATUS( status, "ResampleMultiSFTs", COMPUTEFSTATRSC );
  ATTATCHSTATUSPTR (status);

  /* initialise the start and end SSB times for the SFTs */
  printf("multiSFTs->data[0]->data[0].epoch = %d %d\n",multiSFTs->data[0]->data[0].epoch.gpsSeconds,multiSFTs->data[0]->data[0].epoch.gpsNanoSeconds);
  SSBstart = XLALGPSGetREAL8(&(multiSSB->data[0]->refTime)) + multiSSB->data[0]->DeltaT->data[0] - 0.5*tSFT*multiSSB->data[0]->Tdot->data[0];
  SSBend = SSBstart;
  XLALGPSSetREAL8(&SSBstart_GPS,SSBstart);
  XLALGPSSetREAL8(&SSBend_GPS,SSBend);

  /* loop over all SFTs to determine the earliest SFT midpoint of the input data in the SSB frame */
  for (i=0;i<multiSSB->length;i++) {

    /* loop over the SFTs for the current detector */
    for (j=0;j<multiSSB->data[i]->DeltaT->length;j++) {
      
      /* compute the MID-POINT of this SFT in the SSB frame */
      double tempt = SSBrefTime + multiSSB->data[i]->DeltaT->data[j];
      
      /*  record it if earlier or later than the current earliest and latest */
      if (tempt<SSBstart) {
        SSBstart = tempt - 0.5*tSFT*multiSSB->data[i]->Tdot->data[j];     /* we approximate the SFT start time in the SSB by T_mid - 0.5*tSFT*(dt_SSB/dt_DET) */
	DETstart = tempt - 0.5*tSFT;
	XLALGPSSetREAL8(&SSBstart_GPS,SSBstart);
      }
      if (tempt>SSBend) {
	SSBend = tempt + 0.5*tSFT*multiSSB->data[i]->Tdot->data[j];       /* we approximate the SFT end time in the SSB by T_mid + 0.5*tSFT*(dt_SSB/dt_DET) */
	DETend = tempt + 0.5*tSFT;
        XLALGPSSetREAL8(&SSBend_GPS,SSBend);
      }
	
    }
    
  }
  printf("final SSBstart = %6.12f SSBend = %6.12f\n",SSBstart,SSBend);

  /* compute resampling parameters */        
  Tprime = SSBend - SSBstart;              /* the total time span of data in the SSB */ 
  T = 1.0/deltaF;                          /* the effective observation time given the user requested frequency resolution */
  printf("Tprime = %6.12f\n",Tprime);
  printf("T = %6.12f\n",T);

  /* if the requested deltaF is too large for the data span then we exit */
  if (Tprime>T) {
    LALPrintError("\nResampleMultiSFTs() failed because the users requested frequency resolution was too coarse for the observation span.\n\n");
    ABORT ( status->statusPtr, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
  }
    
  N = ceil(T*SFTbins*SFTdeltaF);                              /* the number of time bins to use in the resampled timeseries */
  deltaT = T/(double)N;                                       /* the sampling time in the resampled timeseries (these is the same for all detectors) */
  printf("N = %d\n",N);
  printf("deltaT = %6.12f deltaTSFT = %6.12f\n",deltaT,deltaTSFT);


  /* alocate memory for the resampled timeseries */
  *multitimeseries = (MultiCOMPLEX8TimeSeries *)calloc(1,sizeof(MultiCOMPLEX8TimeSeries));
  (*multitimeseries)->length = numDetectors;
  (*multitimeseries)->f0 = fhet;
  (*multitimeseries)->deltaT = deltaT;
  (*multitimeseries)->epoch = SSBstart_GPS;
  (*multitimeseries)->refTime = SSBrefTime_GPS;
  (*multitimeseries)->Fat = (COMPLEX8Vector **)LALCalloc(numDetectors,sizeof(COMPLEX8Vector *));
  (*multitimeseries)->Fbt = (COMPLEX8Vector **)LALCalloc(numDetectors,sizeof(COMPLEX8Vector *));
  if (((*multitimeseries)->Fat == NULL)||((*multitimeseries)->Fbt == NULL)) {
    LALPrintError("\nComputeFStatFreqBand_RS() was unable to allocate memory for the resampled timeseries.\n\n");
    ABORT ( status->statusPtr, COMPUTEFSTATRSC_EMEM, COMPUTEFSTATRSC_MSGEMEM );
  }
  for (i=0;i<numDetectors;i++) {
    (*multitimeseries)->Fat[i] = XLALCreateCOMPLEX8Vector(N);
    (*multitimeseries)->Fbt[i] = XLALCreateCOMPLEX8Vector(N);
    {
      INT4 err = xlalErrno;
      if ( err != XLAL_SUCCESS ) {
	ABORT ( status, err, "XLALGPSLeapSeconds() failed!\n");
      }
    }
  }
     
  /* allocate memory for the temporary complex time domain vector and create FFT plan */
  out = XLALCreateCOMPLEX8Vector(SFTbins);
  prev = XLALCreateCOMPLEX8FFTPlan(SFTbins,0,0); 
  {
    INT4 err = xlalErrno;
    if ( err != XLAL_SUCCESS ) {
      ABORT ( status, err, "XLALGPSLeapSeconds() failed!\n");
    }
  }
  
  /* loop over each detector */
  for (i=0;i<numDetectors;i++) {
    
    REAL8 tdiffstart = SSBstart - DETstart;

    /* initialise output timeseries vector to zeros */
    memset((*multitimeseries)->Fat[i]->data,0,N*sizeof(COMPLEX8));
    memset((*multitimeseries)->Fbt[i]->data,0,N*sizeof(COMPLEX8));

    /* loop over the SFTs for this detector */
    for (j=0;j<multiSFTs->data[i]->length;j++) {
      
      COMPLEX8Vector *in = NULL;                                                      /* pointer for the input to the FFT function (points to SFT data) */
      REAL8 a = (REAL8)multiAMcoef->data[i]->a->data[j];                              /* value of the antenna pattern a(t) at the MID-POINT of the SFT */
      REAL8 b = (REAL8)multiAMcoef->data[i]->b->data[j];                              /* value of the antenna pattern b(t) at the MID-POINT of the SFT */

      REAL8 SFTstartDET = XLALGPSGetREAL8(&(multiSFTs->data[i]->data[j].epoch));              /* START time of the SFT at the detector */
      REAL8 SFTmidDET = SFTstartDET + 0.5*tSFT;                                       /* MID-POINT time of the SFT at the detector */      
      REAL8 SFTmidSSB = SSBrefTime + multiSSB->data[i]->DeltaT->data[j];              /* MID-POINT time of the SFT at the SSB */  
      REAL8 Tdot = multiSSB->data[i]->Tdot->data[j];                                  /* the instantaneous time derivitive dt_SSB/dt_DET at the MID-POINT of the SFT */
      
      REAL8 tSFT_SSB = tSFT*Tdot;                                                     /* the approximated span of the SFT in the SSB -> tSFT_SSB = tSFT_DET*(dt_SSB/dt_DET) */
      
      REAL8 SFTstartSSB = SFTmidSSB - 0.5*tSFT*Tdot;                                  /* the approximated start time of the SFT at the SSB */
      REAL8 SFTendSSB = SFTmidSSB + (0.5*tSFT - deltaTSFT)*Tdot;                      /* the approximated time of the last SFT datum in the SSB */
      UINT4 start = floor(0.5 + (SFTstartSSB - SSBstart)/deltaT);                     /* the index of the resampled timeseries corresponding to the first datum of the SFT */ 
      UINT4 end = floor(0.5 + (SFTendSSB - SSBstart)/deltaT);                         /* the index of the resampled timeseries corresponding to the last datum of the SFT */ 
      REAL8 avg = 0.0;
      REAL8 tdiffstartSFT = SFTstartSSB - SFTstartDET;

      printf("SFT data = %6.12f %6.12f\n",multiSFTs->data[i]->data[j].data->data[0].re,multiSFTs->data[i]->data[j].data->data[0].im);
      printf("SFTmidDET = %6.12f SFTmidSSB = %6.12f\n",SFTmidDET,SFTmidSSB);
      printf("SFTstartDET = %6.12f SFTstartSSB = %6.12f\n",SFTstartDET,SFTstartSSB);
      printf("a = %6.12f b = %6.12f\n",a,b);
      printf("start = %d end = %d (SFTbins = %d)\n",start,end,SFTbins);
      printf("tSFT_SSB = %6.12f\n",tSFT_SSB);
      printf("Tdot = %6.12f\n",Tdot);

      /* inverse Fourier transform this SFT into the complex time domain */
      in = multiSFTs->data[i]->data[j].data;
      if (XLALCOMPLEX8VectorFFT(out,in,prev)!= XLAL_SUCCESS) 
	ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");

      /* window timeseries */
      /* for (k=0;k<10;k++) { */
/* 	out->data[k].re *= 0.5*(1.0 - cos(k*LAL_PI/10)); */
/* 	out->data[k].im *= 0.5*(1.0 - cos(k*LAL_PI/10)); */
/*       } */
/*       for (k=0;k<10;k++) { */
/* 	out->data[SFTbins-1-k].re *= 0.5*(1.0 - cos(k*LAL_PI/10)); */
/* 	out->data[SFTbins-1-k].im *= 0.5*(1.0 - cos(k*LAL_PI/10)); */
/*       } */

      /* Loop Over SSB timesamples for this SFT and for each uniformly spaced SSB sample we find the corresponding closest DET sample */
      for (k=start;k<=end;k++) {

	REAL8 tSSB = SSBstart + k*deltaT;                                  /* the SSB time of the current resampled time sample */
	REAL8 tDET = SFTmidDET + (tSSB - SFTmidSSB)/Tdot;                  /* the approximated DET time of the current resampled time sample */
 	REAL8 deltatdiff = SSBrefTime - SFTstartDET - (tSSB - tDET);                       /* the difference t_SSB - t_DET */
 	REAL8 hetphase = LAL_TWOPI*deltatdiff*(*multitimeseries)->f0;     /* the phase for the heterodyne correction -> 2*pi*f_het*(t_SSB-t_DET) */
 	REAL8 cosphase = cos(hetphase);                                    /* the cosine of the heterodyne phase (real part) */
 	REAL8 sinphase = sin(hetphase);                                    /* the sine of the heterodyne phase (imaginary part) */
	UINT4 SFTidx = floor(0.5 + (tDET - SFTstartDET)/deltaTSFT);        /* the SFT time index corresponding to the current SSB time sample */
	/* if (SFTidx<SFTbins) break; */

	(*multitimeseries)->Fat[i]->data[k].re = a*(out->data[SFTidx].re*cosphase - out->data[SFTidx].im*sinphase);
	(*multitimeseries)->Fat[i]->data[k].im = a*(out->data[SFTidx].im*cosphase + out->data[SFTidx].re*sinphase);

	(*multitimeseries)->Fbt[i]->data[k].re = b*(out->data[SFTidx].re*cosphase - out->data[SFTidx].im*sinphase);
	(*multitimeseries)->Fbt[i]->data[k].im = b*(out->data[SFTidx].im*cosphase + out->data[SFTidx].re*sinphase);
	/* printf("k = %d tSSB = %6.6f tDET = %6.6f tdiff = %6.6f hetphase = %6.6f SFTidx = %d data = %6.6f %6.6f\n",k,tSSB,tDET,tdiff,hetphase,SFTidx,out->data[SFTidx].re,out->data[SFTidx].im);  */
	if ((k==end)||(k==start)) printf("k = %d SFTidx = %d SFTbins = %d Fa(t) = %6.12f\n",k,SFTidx,SFTbins,(*multitimeseries)->Fat[i]->data[k].re);
	avg += sqrt(out->data[SFTidx].re*out->data[SFTidx].re + out->data[SFTidx].im*out->data[SFTidx].im);
      }
      printf("\navg = %6.12f\n\n",avg/(end-start));
      /* reinitialise the FFT output */
      memset(out->data,0,SFTbins*sizeof(COMPLEX8));

    }
    
  }
  
  /* free memory */
  XLALDestroyCOMPLEX8FFTPlan(prev);
  XLALDestroyCOMPLEX8Vector(out);
  {
    INT4 err = xlalErrno;
    if ( err != XLAL_SUCCESS ) {
      ABORT ( status, err, "XLALGPSLeapSeconds() failed!\n");
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ResampleMultiSFTs() */
