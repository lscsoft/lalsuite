/*  Copyright (C) 2010 Chris Messenger
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

/** \author C.Messenger
 * \ingroup lalapps_pulsar_Xray
 * \file
 * \brief
 * This code is designed to compute the Bayes factor for a semi-coherent analysis
 * of input SFT data specific to searching for continuous signals in a binary system.
 *
 * It generates likelihood samples from a coarse grid of templates placed on each SFT and
 * combines them using a fine binary template band.  The parameter space is integrated over
 * and a Bayes factor is produced.
 *
 */

/***********************************************************************************************/
/* includes */
#include "config.h"
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>




#include <gsl/gsl_interp.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_spline.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_rng.h>           /* for random number generation */
#include <gsl/gsl_randist.h>       /* for random number generation */
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_log.h>        /* for log computation */
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/ComplexFFT.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>

#include "SemiCoherent.h"


/***********************************************************************************************/
/* some useful macros */

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/** this function initialises the gsl random number generation
 *
 * If the input seed is zero then a random seed is drawn from
 * /dev/urandom.
 *
 */
int XLALInitgslrand(gsl_rng **gslrnd,     /**< [out] the gsl random number generator */
		    INT8 seed             /**< [in] the random number generator seed */
		    )
{
  FILE *devrandom = NULL;      /* pointer to the /dev/urandom file */

  /* if the seed is 0 then we draw a random seed from /dev/urandom */
  if (seed == 0) {

    /* open /dev/urandom */
    if ((devrandom=fopen("/dev/urandom","r")) == NULL)  {
      LogPrintf(LOG_CRITICAL,"%s: Error, unable to open device /dev/random\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* read a random seed */
    if (fread((void*)&seed,sizeof(INT8),1,devrandom) != 1) {
      LogPrintf(LOG_CRITICAL,"%s: Error, unable to read /dev/random\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }
    fclose(devrandom);

  }

  /* setup gsl random number generation */
  *gslrnd = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(*gslrnd,seed);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Free the memory allocated within a REAL4DemodulatedPowerVector structure
 *
 */
int XLALFreeREAL4DemodulatedPowerVector(REAL4DemodulatedPowerVector *power            /**< [in] the data to be freed */
					)
{
  UINT4 i;                     /* counter */

  /* free each segment */
  for (i=0;i<power->length;i++) {
    XLALDestroyREAL4Vector(power->segment[i]->data);
    XLALFree(power->segment[i]);
  }
  XLALFree(power->segment);
  XLALFree(power);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Free the memory allocated within a ParameterSpace structure
 *
 */
int XLALFreeParameterSpace(ParameterSpace *pspace            /**< [in] the parameter space to be freed */
			   )
{

  /* free parameter space */
  XLALFree(pspace->space->data);
  XLALFree(pspace->space);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}




/** Compute the instantaneous frequency derivitives for a given binary template and segment
 *
 */
int XLALComputeBinaryFreqDerivitives(Template *fdots,                        /**< [out] the frequency derivitives */
				     Template *bintemp,                      /**< [in] the binary template */
				     REAL8 tmid                              /**< [in] the midpoint time of the segment */
				     )
{
  UINT4 n;                               /* counters */
  REAL8 nu = bintemp->x[0];              /* define nu */
  REAL8 asini = bintemp->x[1];           /* define asini */
  REAL8 tasc = bintemp->x[2];            /* define tasc */
  REAL8 omega = bintemp->x[3];           /* define omega */

  /* precompute repeated quantities */
  REAL8 nuasiniomega = nu*asini*omega;
  REAL8 orbphase = omega*(tmid-tasc);
  REAL8 omegan = 1;

  /* the instantanous frequency is therefore f0 = nu - a*nu*W*cos(W*(t-tasc) ) */
  fdots->x[0] = nu - nuasiniomega*cos(orbphase);

  /* the instantanous nth frequency derivitive is therefore fn = - a * nu * W^(n+1) * cos ( W*(t-tasc) + n*pi/2 ) */
  for (n=1;n<fdots->ndim;n++) {
    omegan *= omega;
    fdots->x[n] = (-1.0)*nuasiniomega*omegan*cos(orbphase + 0.5*n*LAL_PI);
  }

  return XLAL_SUCCESS;

}

/** Compute the next template in a grid
 *
 * the templates are generated sequentially from the grid parameters file.  The n-dimensional
 * virtual indices of each dimension are also generated
 *
 */
int XLALGetNextRandomBinaryTemplate(Template **temp,                        /**< [out] the signal model template parameters */
				    GridParameters *gridparams,             /**< [in] the parameter space grid params */
				    ParameterSpace *space,		/**< UNDOCUMENTED */
				    void *r		/**< UNDOCUMENTED */
				    )
{

  /* if the input template is null then we allocate memory and assume we are on the first template */
  if ((*temp) == NULL) {
    if ( ((*temp) = XLALCalloc(1,sizeof(Template))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( ((*temp)->x = XLALCalloc(gridparams->ndim,sizeof(REAL8))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( ((*temp)->idx = XLALCalloc(gridparams->ndim,sizeof(UINT4))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    (*temp)->currentidx = 0;
    (*temp)->ndim = gridparams->ndim;

  }
  else if ((*temp)->currentidx == (UINT4)gridparams->Nr - 1) {

    /* free binary template memory */
    XLALFree((*temp)->x);
    XLALFree((*temp)->idx);
    XLALFree(*temp);

    LogPrintf(LOG_DEBUG,"%s: at last template.\n",__func__);
    return 0;
  }
  else (*temp)->currentidx++;

  REAL8 n1 = space->space->data[0].min;
  REAL8 n2 = space->space->data[0].max;
  REAL8 a1 = space->space->data[1].min;
  REAL8 a2 = space->space->data[1].max;
  REAL8 t1 = space->space->data[2].min;
  REAL8 t2 = space->space->data[2].max;
  REAL8 O1 = space->space->data[3].min;
  REAL8 O2 = space->space->data[3].max;
  REAL8 nu;
  REAL8 Om;
  REAL8 a;
  REAL8 tasc;

  /* while we haven't selected a template in the space */
  INT4 flag = 0;
  while (!flag) {

    REAL8 temp1 = gsl_ran_flat((gsl_rng*)r,0,1);
    nu = n2*pow((pow(n1/n2,4.0) + (1.0 - pow(n1/n2,4.0))*temp1),1.0/4.0);
    REAL8 temp2 = gsl_ran_flat((gsl_rng*)r,0,1);
    a = a2*pow((pow(a1/a2,3.0) + (1.0 - pow(a1/a2,3.0))*temp2),1.0/3.0);
    REAL8 temp3 = gsl_ran_flat((gsl_rng*)r,0,1);
    tasc = t1 + (t2-t1)*temp3;
    REAL8 temp4 = gsl_ran_flat((gsl_rng*)r,0,1);
    Om = O2*pow((pow(O1/O2,5.0) + (1.0 - pow(O1/O2,5.0))*temp4),1.0/5.0);

    /* fprintf(stdout,"%f (%f %f) %f (%f %f) %f (%f %f) %f (%f %f)\n",nu,n1,n2,a,a1,a2,tasc,t1,t2,Om,O1,O2); */

    if ((nu>n1)&&(nu<n2)&&(a>a1)&&(a<a2)&&(tasc>t1)&&(tasc<t2)&&(Om>O1)&&(Om<O2)) flag = 1;

  }

  (*temp)->x[0] = nu;
  (*temp)->x[1] = a;
  (*temp)->x[2] = tasc;
  (*temp)->x[3] = Om;

  return 1;

}

/** Compute the next template in a grid
 *
 * the templates are generated sequentially from the grid parameters file.  The n-dimensional
 * virtual indices of each dimension are also generated
 *
 */
int XLALGetNextTemplate(Template **temp,                        /**< [out] the signal model template parameters */
			GridParameters *gridparams,              /**< [in] the parameter space grid params */
			ParameterSpace *space,		/**< UNDOCUMENTED */
			UNUSED void *r		/**< UNDOCUMENTED */
			)
{
  UINT4 idx;                             /* the index variable */
  INT4 j;                                /* counters */

  /* check input */
  if (space != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input ParameterSpace structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* if the input template is null then we allocate memory and assume we are on the first template */
  if ((*temp) == NULL) {
    if ( ((*temp) = XLALCalloc(1,sizeof(Template))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( ((*temp)->x = XLALCalloc(gridparams->ndim,sizeof(REAL8))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( ((*temp)->idx = XLALCalloc(gridparams->ndim,sizeof(UINT4))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    (*temp)->currentidx = 0;
    (*temp)->ndim = gridparams->ndim;

  }
  else if ((*temp)->currentidx == gridparams->max - 1) {

    /* free binary template memory */
    XLALFree((*temp)->x);
    XLALFree((*temp)->idx);
    XLALFree(*temp);

    LogPrintf(LOG_DEBUG,"%s: at last template.\n",__func__);
    return 0;
  }
  else (*temp)->currentidx++;

  /* initialise index */
  idx = (*temp)->currentidx;

  /* loop over each dimension and obtain the index for that dimension (store both) */
  for (j=gridparams->ndim-1;j>=0;j--) {

    /* compute the index for the j'th dimension and compute the actual value */
    UINT4 q = idx/gridparams->prod[j];
    (*temp)->x[j] = gridparams->grid[j].min + q*gridparams->grid[j].delta;
    (*temp)->idx[j] = q;

    /* update the index variable for the next dimension */
    idx = idx - q*gridparams->prod[j];

  }

  /* update index */
 /*  (*temp)->currentidx++; */


  return 1;

}

/** Compute the demodulated power for a single SFT
 *
 * This involves taking the downsampled SFT timeseries and multiplying by the
 * timeseries spin-derivitive templates in turn.  Then we inverse FFT the result
 * and square to obtain the power at a given set of freq derivitive parameters.
 *
 */
int XLALApplyPhaseCorrection(COMPLEX8TimeSeries **outts,            /**< [out] the output complex time series */
			     COMPLEX8TimeSeries *ints,          /**< [in] the input complex time series */
			     Template *fn                       /**< [in] the spin derivitives */
			     )
{

  UINT4 j,k;

  /* validate input arguments */
  if ((*outts) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8Vector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (ints == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8 structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* apply time domain phase correction - first loop over time and then spin derivitive */
  for (j=0;j<ints->data->length;j++) {

    /* compute phase correction including heterodyne to shift frequencies to match up with grid */
    const REAL8 t = j*ints->deltaT - 0.5*ints->deltaT*ints->data->length;
    REAL8 tn = t;
    REAL8 arg = 0.0;
    UINT4 fac = 1;

    /* loop over each spin derivitive and add to phase contribution for current time sample */
    for (k=0;k<fn->ndim;k++) {
      tn *= t;
      fac *= k+2;
      arg += (-1.0)*LAL_TWOPI*fn->x[k]*tn/fac;
    }

    /* compute real and imaginary parts of phase correction timeseries */
    REAL8 xr = cos(arg);
    REAL8 xi = sin(arg);
    COMPLEX8 x = xr + I*xi;

    /* multiply data by phase correction - leave the zero-padding */
    (*outts)->data->data[j] = ints->data->data[j]*x;

  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Compute the demodulated power for a single SFT
 *
 * This involves taking the downsampled SFT timeseries and multiplying by the
 * timeseries spin-derivitive templates in turn.  Then we inverse FFT the result
 * and square to obtain the power at a given set of freq derivitive parameters.
 *
 */
int XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(COMPLEX8FrequencySeries **fs,     /**< [out] the over-resolved frequency series */
						    const COMPLEX8TimeSeries *ts,     /**< [in] the downsampled SFT data */
						    GridParameters **gridparams        /**< [in/out] the spin derivitive gridding parameters */
						    )
{

  COMPLEX8FFTPlan *plan = NULL;           /* plan for the inverse FFT */
  COMPLEX8Vector *temp_input = NULL;      /* the temporary input of the inverse FFT = data*template */
  COMPLEX8Vector *temp_output = NULL;     /* the temporary output of the inverse FFT = data*template */
  UINT4 j;

  /* validate input arguments */
  if ((*fs) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8FrequencySeries structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8TimeSeries structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (gridparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParameters structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  UINT4 N = ts->data->length;                         /* the original length of the time series */
  REAL8 deltaT = ts->deltaT;                          /* the fixed time sampling of the input data */
  REAL8 T = N*deltaT;                                 /* the intrinsic duration */
  REAL8 deltaF = 1.0/(ts->deltaT*ts->data->length);   /* the intrinsic deltaF */
  /* define new deltaF accounting for the fixed deltaT */
  REAL8 newdeltaF = (*gridparams)->grid[0].delta;
  UINT4 newN = ceil((T*deltaF/newdeltaF)/deltaT);
  REAL8 newT = deltaT*newN;
  newdeltaF = 1/newT;

  /* allocate memory for the temporary zero-padded input data */
  if ( (temp_input = XLALCreateCOMPLEX8Vector(newN)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d.\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  /* allocate memory for the temporary output data */
  if ( (temp_output = XLALCreateCOMPLEX8Vector(newN)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d.\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }

  /* create a forward complex fft plan */
  if ((plan = XLALCreateForwardCOMPLEX8FFTPlan(newN,0)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateForwardCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }

  /* initialise the input data with zeros */
  memset(temp_input->data,0.0,temp_input->length*sizeof(COMPLEX8));

  /* put the input data into the temporary input structure and normalise it with deltaT */
  for (j=0;j<N;j++) {
    temp_input->data[j] = ts->data->data[j]*deltaT;
  }

  /* FFT the data */
  if (XLALCOMPLEX8VectorFFT(temp_output,temp_input,plan)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX8VectorFFT() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the FFT\n",__func__);

  /* and modify grid based on new resolution  - just shift the lowest frequency to match up with existing bins */
  /* also update the resolution BUT keep the same number of bins */
  INT4 binoffset = (INT4)floor(((*gridparams)->grid[0].min - ts->f0)/newdeltaF);
  if (binoffset<0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, we have a negative binoffset.\n",__func__);
    return XLAL_EINVAL;
  }
  (*gridparams)->grid[0].min = ts->f0 + newdeltaF*binoffset;
  (*gridparams)->grid[0].delta = newdeltaF;

  /* allocate memory for output */
  if ( ((*fs) = XLALCreateCOMPLEX8FrequencySeries("FS",&(ts->epoch),(*gridparams)->grid[0].min,(*gridparams)->grid[0].delta,&lalDimensionlessUnit,(*gridparams)->grid[0].length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL4Vector() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }

  /*  extract desired frequencies */
  for (j=0;j<(*gridparams)->grid[0].length;j++) {
    INT4 k = j + binoffset;
    (*fs)->data->data[j] = temp_output->data[k];
  }

  /* free memory */
  XLALDestroyCOMPLEX8Vector(temp_input);
  XLALDestroyCOMPLEX8Vector(temp_output);
  XLALDestroyCOMPLEX8FFTPlan(plan);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Compute the demodulated power for all downsampled timeseries
 *
 * This function is simply a wrapper for XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries()
 *
 */
int XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(REAL4DemodulatedPowerVector **power,     /**< [out] the spin derivitive demodulated power */
						  COMPLEX8TimeSeriesArray *dsdata,         /**< [in] the downsampled SFT data */
						  GridParametersVector *gridparams,        /**< [in/out] the spin derivitive gridding parameters */
	                                          FILE *fp		/**< UNDOCUMENTED */
						)
{
  UINT4 i;
  REAL8 maxspintemp[NFREQMAX];

  /* validate input */
  if ((*power) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output REAL4DemodulatedPowerVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (dsdata == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8TimeSeriesArray structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (gridparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParametersVector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (dsdata->length != gridparams->length) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, length of downsampled data vector and grid parameters vector not equal.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* allocate memory */
  if ( ((*power) = (REAL4DemodulatedPowerVector*)XLALCalloc(1,sizeof(REAL4DemodulatedPowerVector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for REAL4DemodulatedPowerVector structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*power)->segment = (REAL4DemodulatedPower**)XLALCalloc(dsdata->length,sizeof(REAL4DemodulatedPower*))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for array of REAL4Vector structures.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*power)->length = dsdata->length;

  /* loop over each segment */
  for (i=0;i<dsdata->length;i++) {

    COMPLEX8TimeSeries *ts = dsdata->data[i];
    COMPLEX8TimeSeries *temp_ts = NULL;
    COMPLEX8FrequencySeries *fs = NULL;
    UINT4 j;
    UINT4 idx = 0;
    REAL8 maxpower = 0.0;
    UINT4 maxdim = 0;

    /* allocate memory for all templates in the frequency derivitive space for this segment */
    if ( ((*power)->segment[i] = XLALCalloc(1,sizeof(REAL4Vector*))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    if ( ((*power)->segment[i]->data = XLALCreateREAL4Vector(gridparams->segment[i]->max)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL4Vector() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }

    /* allocate memory for temporary complex timeseries */
    if ( (temp_ts = XLALCreateCOMPLEX8TimeSeries("DS",&(ts->epoch),ts->f0,ts->deltaT,&lalDimensionlessUnit,ts->data->length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8TimeSeries() failed to allocate memory for inverse FFT output.\n",__func__);
      return XLAL_ENOMEM;
    }

    /* redefine grid on spin derivitives only */
    GridParameters tempgrid;
    Template *spintemp = NULL;
    tempgrid.ndim = gridparams->segment[i]->ndim - 1;
    tempgrid.grid = XLALCalloc(tempgrid.ndim,sizeof(Grid));
    tempgrid.prod = XLALCalloc(tempgrid.ndim,sizeof(UINT4));
    tempgrid.max = 1;
    for (j=0;j<tempgrid.ndim;j++) {
      tempgrid.grid[j].min = gridparams->segment[i]->grid[j+1].min;
      tempgrid.grid[j].delta = gridparams->segment[i]->grid[j+1].delta;
      tempgrid.grid[j].length = gridparams->segment[i]->grid[j+1].length;
      tempgrid.max *= tempgrid.grid[j].length;
    }

    tempgrid.prod[0] = 1;
    for (j=1;j<tempgrid.ndim;j++) {
      tempgrid.prod[j] = tempgrid.prod[j-1]*tempgrid.grid[j-1].length;
    }

    /* loop over spin derivitive values - not over frequency */
    ParameterSpace *temp = NULL;
    while (XLALGetNextTemplate(&spintemp,&tempgrid,temp,NULL)) {

      /* reinitilaise the reused complex timeseries */
      fs = NULL;

      /* apply phase correction to complex timeseries */
      if (XLALApplyPhaseCorrection(&temp_ts,ts,spintemp)) {
	LogPrintf(LOG_CRITICAL,"%s : XLALApplyPhseCorrection() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s : applied phase correction for template index %d/%d on segment %d/%d\n",__func__,spintemp->currentidx,tempgrid.max,i,dsdata->length);

      /* convert to the complex frequency domain - on the frequency grid specified */
      if (XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(&fs,temp_ts,&(gridparams->segment[i]))) {
	LogPrintf(LOG_CRITICAL,"%s : XLALSFTCOMPLEX8TimeseriesToCOMPLEX8FrequencySeries() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s : computed demodulated frequency series for SFT %d/%d\n",__func__,i+1,dsdata->length);

      /* compute power and store it */
      for (j=0;j<fs->data->length;j++) {
	(*power)->segment[i]->data->data[idx] = crealf(fs->data->data[j])*crealf(fs->data->data[j]) + cimagf(fs->data->data[j])*cimagf(fs->data->data[j]);

        /* find maximum for this segment */
        if (fp!=NULL) {
          if ((*power)->segment[i]->data->data[idx]>maxpower) {
            maxpower = (*power)->segment[i]->data->data[idx];
            maxdim = spintemp->ndim + 1;
            UINT4 k;
            maxspintemp[0] = fs->f0 + j*fs->deltaF;
            for (k=0;k<spintemp->ndim;k++) maxspintemp[k+1] = spintemp->x[k];
          }
        }
	idx++;
      }
      (*power)->segment[i]->epoch.gpsSeconds = ts->epoch.gpsSeconds;
      (*power)->segment[i]->epoch.gpsNanoSeconds = ts->epoch.gpsNanoSeconds;

      /* free memory */
      XLALDestroyCOMPLEX8FrequencySeries(fs);

    } /* end loop over spin derivitive templates */

    /* output loudest segment candidate - normalise to be a chi-squared variable */
    if (fp!=NULL) {
      fprintf(fp,"%d\t%d\t%d\t",ts->epoch.gpsSeconds,ts->epoch.gpsNanoSeconds,(INT4)floor(ts->deltaT*ts->data->length+0.5));
      for (j=0;j<maxdim;j++) fprintf(fp,"%.12f\t",maxspintemp[j]);
      for (j=maxdim;j<NFREQMAX;j++) fprintf(fp,"0.0\t");
      fprintf(fp,"%.12f\n",2.0*maxpower); /*/background->data[i]); */
    }

    /* free memory */
    XLALDestroyCOMPLEX8TimeSeries(temp_ts);
    XLALFree(tempgrid.grid);
    XLALFree(tempgrid.prod);

  } /* end the loop over segments */

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Compute the gridding parameters on spin derivitives
 *
 * The circular orbit binary phase model is phi = 2*pi*nu*( (t-tref) - a*sin( W*(t-tasc) )
 * from which we compute the min and maximum instantaneous spin derivitives.
 *
 */
int XLALComputeFreqGridParams(GridParameters **gridparams,              /**< [out] the gridding parameters */
			      REAL8Space *space,                        /**< [in] the orbital parameter space */
			      REAL8 tmid,                               /**< [in] the segment mid point */
			      REAL8 Tseg,                               /**< [in] the segment length */
			      REAL8 mu                                  /**< [in] the required mismatch */
			      )
{
  UINT4 i,j,k,l;                         /* counters */
  INT4 n;                                /* counter */
  REAL8 fnmin[NFREQMAX],fnmax[NFREQMAX]; /* min and max values of spin derivitives */
  INT4 dim[NFREQMAX];                    /* flag indicating whether a dimension has width */
  INT4 ndim = -1;                        /* the number of spin derivitive dimensions required */
  Template fdots;                        /* template for an instance of spin parameters */
  Template bintemp;                      /* template for instance of binary parameters */
  UINT4 ngrid = 100;                     /* the number of grid points per omega and tasc to search for finding true fdot spans per sft */

  /* validate input arguments */
  if ((*gridparams) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GridParams structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (space == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input ParameterSpace structure == NULL.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }
   if (tmid < 0) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GPS time < 0.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }
   if (Tseg < 0) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input Tseg parameter < 0.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }
   if ((mu < 0)||(mu>1)) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input mismatch parameter, not in range 0 -> 1.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }

   /* allocte memory */
   if (((*gridparams) = (GridParameters*)XLALCalloc(1,sizeof(GridParameters))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

   /* allocate memory for the fdots */
   if ((fdots.x = XLALCalloc(NFREQMAX,sizeof(REAL8))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
     XLAL_ERROR(XLAL_ENOMEM);
   }
   if ((bintemp.x = XLALCalloc(NBINMAX,sizeof(REAL8))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
     XLAL_ERROR(XLAL_ENOMEM);
   }
   fdots.ndim = NFREQMAX;
   bintemp.ndim = NBINMAX;

   /* initialise the min and max spin derivitives */
   for (n=0;n<NFREQMAX;n++) {
     fnmin[n] = 1e38;
     fnmax[n] = -1e38;
   }

   /* loop over each parameter in turn and compute the spin derivitives at the corners of the parameter space */
   for (i=0;i<2;i++) {     /* nu */
     if (i==0) bintemp.x[0] = space->data[0].min;
     else bintemp.x[0] = space->data[0].max;

     for (j=0;j<2;j++) {    /* a */
       if (j==0) bintemp.x[1] = space->data[1].min;
       else bintemp.x[1] = space->data[1].max;

       /* tasc and omega are the problematic ones so we'll perform a fine grid search over them */
       for (k=0;k<ngrid;k++) {   /* tasc */
	 bintemp.x[2] = space->data[2].min + k*(space->data[2].max-space->data[2].min)/(ngrid-1);

	 for (l=0;l<ngrid;l++) {  /* omega */
	   bintemp.x[3] = space->data[3].min + l*(space->data[3].max-space->data[3].min)/(ngrid-1);

	   if (XLALComputeBinaryFreqDerivitives(&fdots,&bintemp,tmid)) {
	     LogPrintf(LOG_CRITICAL,"%s : XLALComputeBinaryFreqDerivitives() failed with error = %d\n",__func__,xlalErrno);
	     XLAL_ERROR(XLAL_EFAULT);
	   }

	   /* find min and max values */
	   for (n=0;n<NFREQMAX;n++) {
	     if (fdots.x[n] < fnmin[n]) fnmin[n] = fdots.x[n];
	     if (fdots.x[n] > fnmax[n]) fnmax[n] = fdots.x[n];
	   }

	 }

       }

     }

   }
   for (n=0;n<NFREQMAX;n++) {
     LogPrintf(LOG_DEBUG,"%s : determined f%d range as [%6.12e -> %6.12e].\n",__func__,n,fnmin[n],fnmax[n]);
   }
   LogPrintf(LOG_DEBUG,"%s : midpoint epoch for this SFT is %6.12f\n",__func__,tmid);

   /* free templates */
   XLALFree(fdots.x);
   XLALFree(bintemp.x);

   /* compute the required dimensionality of the frequency derivitive grid */
   /* we check the width of a 1-D template across each dimension span */
   for (n=0;n<NFREQMAX;n++) {
     REAL8 gnn = pow(LAL_PI,2.0)*pow(Tseg,2*n+2)/(pow(2.0,2*n)*(2*n+3.0));
     REAL8 deltafn = 2.0*sqrt(mu/gnn);
     REAL8 span = fnmax[n] - fnmin[n];
     dim[n] = 0;
     if (span > deltafn) dim[n] = 1;
     LogPrintf(LOG_DEBUG,"%s : single template span for %d'th derivitive = %e.\n",__func__,n,deltafn);
   }
   n = NFREQMAX-1;
   while ( (n>=0) && (ndim == -1) ) {
     if (dim[n] > 0) ndim = n+1;
     n--;
   }
   if (ndim < 0) {
      LogPrintf(LOG_CRITICAL,"%s: dimensionality of frequency space < 0.  No templates required.\n",__func__);
      return XLAL_EINVAL;
   }
   LogPrintf(LOG_DEBUG,"%s : determined dimensionality of frequency space = %d.\n",__func__,ndim);

   /* allocate memory to the output */
   if ( ((*gridparams)->grid = XLALCalloc(ndim,sizeof(Grid))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
     return XLAL_ENOMEM;
   }
   (*gridparams)->ndim = ndim;
   LogPrintf(LOG_DEBUG,"%s : allocated memory for the output grid parameters.\n",__func__);

   /* Compute the grid spacing, grid start and span for each spin derivitive dimension */
   for (n=0;n<ndim;n++) {

     /* compute diagonal metric element and corresponding spacing */
     REAL8 gnn = pow(LAL_PI,2.0)*pow(Tseg,2*n+2)/(pow(2.0,2*n)*(2*n+3.0));
     REAL8 deltafn = 2.0*sqrt(mu/(ndim*gnn));

     /* compute number of grid points in this dimension and enforce a grid centered on the middle of the parameter space */
     INT4 length = (INT4)ceil((fnmax[n]-fnmin[n])/deltafn) + 2*NBINS;                      /* add bins at each end for safety */
     REAL8 minfn = 0.5*(fnmin[n]+fnmax[n]) - 0.5*(length-1)*deltafn;

     (*gridparams)->grid[n].delta = deltafn;
     (*gridparams)->grid[n].oneoverdelta = 1.0/deltafn;
     (*gridparams)->grid[n].length = length;
     (*gridparams)->grid[n].min = minfn;
     snprintf((*gridparams)->grid[n].name,LALNameLength,"f%d",n);

     LogPrintf(LOG_DEBUG,"%s : %s -> [%e - %e] (%e) %d grid points.\n",__func__,(*gridparams)->grid[n].name,(*gridparams)->grid[n].min,
	       (*gridparams)->grid[n].min+((*gridparams)->grid[n].length-1)*(*gridparams)->grid[n].delta,
	       (*gridparams)->grid[n].delta,(*gridparams)->grid[n].length);
   }
   LogPrintf(LOG_DEBUG,"%s : computed output grid parameters.\n",__func__);

   /* compute some internally required parameters for the grid */
   if ( ((*gridparams)->prod = XLALCalloc(ndim,sizeof(UINT4))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
     return XLAL_ENOMEM;
   }
   (*gridparams)->ndim = ndim;
   (*gridparams)->mismatch = mu;
   (*gridparams)->max = 1;
   for (k=0;k<(*gridparams)->ndim;k++) (*gridparams)->max *= (*gridparams)->grid[k].length;

   (*gridparams)->prod[0] = 1;
   for (k=1;k<(*gridparams)->ndim;k++) (*gridparams)->prod[k] = (*gridparams)->prod[k-1]*(*gridparams)->grid[k-1].length;

   LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
   return XLAL_SUCCESS;

 }

/*******************************************************************************/
/** Compute the gridding parameters on spin derivitives for all segments
 *
 * This is simply a wrapper for the single segment function
 *
 */
int XLALComputeFreqGridParamsVector(GridParametersVector **freqgridparams,    /**< [out] the gridding parameters */
				    REAL8Space *space,                        /**< [in] the orbital parameter space */
				    SFTVector *sftvec,                        /**< [in] the input SFTs */
				    REAL8 mu                                  /**< [in] the required mismatch */
				    )
{
  UINT4 i;                              /* counter */

  /* validate input arguments */
  if ((*freqgridparams) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GridParamsVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (space == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input REAL8Space structure == NULL.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((mu < 0)||(mu>1)) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input mismatch parameter, not in range 0 -> 1.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* allocate memory for each set of grid parameters */
  if (((*freqgridparams) = (GridParametersVector*)XLALCalloc(1,sizeof(GridParametersVector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() falied with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*freqgridparams)->segment = (GridParameters**)XLALCalloc(sftvec->length,sizeof(GridParameters*))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a COMPLEX8TimeSeriesArray structure\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*freqgridparams)->length = sftvec->length;

  /* loop over each SFT */
  for (i=0;i<sftvec->length;i++) {

    REAL8 t0 = XLALGPSGetREAL8(&(sftvec->data[i].epoch));
    REAL8 tsft = 1.0/sftvec->data[i].deltaF;
    REAL8 tmid = t0 + 0.5*tsft;

    if (XLALComputeFreqGridParams(&((*freqgridparams)->segment[i]),space,tmid,tsft,mu)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALComputeFreqGridParams() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s : computed frequency grid for SFT %d/%d\n",__func__,i+1,sftvec->length);

  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Inverse FFT all narrowband SFTs
 *
 * In order to apply the frequency derivitive corrections we must work in the time domain
 * so here we convert all SFTs to the complex time domain.
 *
 */
int XLALSFTToCOMPLEX8TimeSeries(COMPLEX8TimeSeries **ts,           /**< [out] the downsampled timeseries */
				COMPLEX8FrequencySeries *sft,      /**< [in] the input SFT vector */
				COMPLEX8FFTPlan **plan             /**< [in/out] the FFT plan */
				)
{

  UINT4 N;                               /* the length of the SFTs */

  /* validate input arguments */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8TimeSeries structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (sft == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFT structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  /* we check that all input SFTs are of identical length so we make a single plan */
  N = sft->data->length;

  /* make the reverse plan if not cached */
  if ((*plan)==NULL) {
    if (((*plan) = XLALCreateReverseCOMPLEX8FFTPlan(N,1)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateReverseCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_EINVAL;
    }
    LogPrintf(LOG_DEBUG,"%s : created the inverse FFT plan\n",__func__);
  }

  COMPLEX8Vector temp_output;
  REAL8 deltaF = sft->deltaF;
  REAL8 Tsft = 1.0/sft->deltaF;
  REAL8 deltaT = Tsft/N;
  UINT4 j;

  /* allocate output memory - create a COMPLEX8TimeSeries */
  if (((*ts) = XLALCreateCOMPLEX8TimeSeries("DS",&(sft->epoch),sft->f0,deltaT,&lalDimensionlessUnit,N)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8TimeSeries() failed to allocate memory for inverse FFT output.\n",__func__);
    return XLAL_ENOMEM;
  }

  /* point to input */
  COMPLEX8Vector temp_input;
  temp_input.length = sft->data->length;
  temp_input.data = (COMPLEX8*)sft->data->data;

  /* point temp output to timeseries */
  temp_output.length = N;
  temp_output.data = (COMPLEX8*)(*ts)->data->data;

  /* perform inverse FFT */
  if (XLALCOMPLEX8VectorFFT(&temp_output, &temp_input,(*plan))) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX8VectorFFT() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_EINVAL;
  }

  /* normalise outputs by multiplying by df */
  for (j=0;j<N;j++) {
    temp_output.data[j] *= deltaF;
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Inverse FFT all narrowband SFTs
 *
 * In order to apply the frequency derivitive corrections we must work in the time domain
 * so here we convert all SFTs to the complex time domain.
 *
 */
int XLALSFTVectorToCOMPLEX8TimeSeriesArray(COMPLEX8TimeSeriesArray **dstimevec,      /**< [out] the downsampled timeseries */
					   SFTVector *sftvec                         /**< [in] the input SFT vector */
					   )
{
  INT4 i;                                /* counter */
  COMPLEX8FFTPlan *plan = NULL;          /* inverse FFT plan */

  /* validate input arguments */
  if ((*dstimevec) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8TimeSeriesArray structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  /* allocate memory for output */
  if (((*dstimevec) = (COMPLEX8TimeSeriesArray*)XLALCalloc(1,sizeof(COMPLEX8TimeSeriesArray))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a COMPLEX8TimeSeriesArray structure\n",__func__);
    return XLAL_ENOMEM;
  }
  if (((*dstimevec)->data = (COMPLEX8TimeSeries**)XLALCalloc(sftvec->length,sizeof(COMPLEX8TimeSeries *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a vector of COMPLEX8TimeSeries pointers\n",__func__);
    return XLAL_ENOMEM;
  }
  (*dstimevec)->length = sftvec->length;                       /* the number of timeseries */
  LogPrintf(LOG_DEBUG,"%s : allocated memory for the output data structure\n",__func__);

  /* loop over each SFT */
  for (i=0;i<(INT4)sftvec->length;i++) {

    /* convert each SFT to complex timeseries */
    if (XLALSFTToCOMPLEX8TimeSeries(&((*dstimevec)->data[i]),&(sftvec->data[i]),&plan)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALSFTToCOMPLEX8TimeSeries() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_EINVAL;
    }

  }
  LogPrintf(LOG_DEBUG,"%s : performed inverse FFT on all %d SFTs\n",__func__,sftvec->length);

  /* free memory */
  XLALDestroyCOMPLEX8FFTPlan(plan);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Read in SFTs to an SFTVector
 *
 */
int XLALReadSFTs(SFTVector **sftvec,        /**< [out] the input SFT data */
		 CHAR *sftbasename,         /**< [in] the SFT file basename to read in */
		 REAL8 freq,                /**< [in] the starting frequency to read in */
		 REAL8 freqband,            /**< [in] the bandwidth to read */
		 INT4 start,                /**< [in] the min GPS time of the input data */
		 INT4 end,                  /**< [in] the max GPS time of the input data*/
		 INT4 tsft		/**< UNDOCUMENTED */
  		 )
{
  static SFTConstraints constraints;
  SFTCatalog *catalog = NULL;
  /* INT4 sft_check_result = 0; */
  REAL8 freqmin,freqmax;
  LIGOTimeGPS *dummy_gpsstart = NULL;
  LIGOTimeGPS *dummy_gpsend = NULL;
  LIGOTimeGPS gpsstart, gpsend;
  /* LALStatus status = blank_status; */             /* for use wih non-XLAL functions */

  /* validate input variables */
  if (*sftvec != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input SFTVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (sftbasename == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input SFT basename string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (freqband < 0 ) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, frequency band must be > 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((start > 0) && (end > 0)) {
    if (start - end >= 0) {
      LogPrintf(LOG_CRITICAL,"%s: Invalid input, the start time %d >= %d end time.\n",__func__,start,end);
      XLAL_ERROR(XLAL_EINVAL);
    }
  }

  /* get sft catalog */
  /* if the input gps times are negative i.e. not set, then we pass null pointers to LALLALSFTDataFind */
  if (start > 0) {
    XLALGPSSetREAL8(&gpsstart,(REAL8)start);
    dummy_gpsstart = &gpsstart;
  }
  if (end > 0) {
    XLALGPSSetREAL8(&gpsend,(REAL8)end);
    dummy_gpsend = &gpsend;
  }
  constraints.minStartTime = dummy_gpsstart;
  constraints.maxStartTime = dummy_gpsend;
  if ((catalog = XLALSFTdataFind(sftbasename, &constraints))==NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Null pointer returned from XLALSFTdataFind.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  LogPrintf(LOG_DEBUG,"%s : found %d SFTs\n",__func__,catalog->length);

  /* define actual frequency range to read in */
  freqmin = freq;
  freqmax = freqmin + freqband;

  /* check CRC sums of SFTs */
  /* LAL_CALL ( LALCheckSFTCatalog ( &status, &sft_check_result, catalog ), &status );
  if (sft_check_result) {
    LogPrintf(LOG_CRITICAL,"%s : LALCheckSFTCatalogSFT() validity check failed with error = %d\n", sft_check_result);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : checked the SFTs\n",__func__); */

  /* load the SFT-vectors */
  if ((*sftvec = XLALLoadSFTs (catalog, freqmin, freqmax))==NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Null pointer returned from XLALLoadSFTs.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  LogPrintf(LOG_DEBUG,"%s : loaded the sfts\n",__func__);

  /* we don't need the original catalog anymore */
  XLALDestroySFTCatalog(catalog);
  LogPrintf(LOG_DEBUG,"%s : destroyed the catalogue(s)\n",__func__);

  /* check if we found any SFTs */
  if ((*sftvec)->length == 0) {
    LogPrintf(LOG_CRITICAL,"%s : No SFTs found in specified frequency range.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* check if the SFTs had the expected length */
  if (fabs((*sftvec)->data[0].deltaF - 1.0/(REAL8)tsft)>1e-12) {
    LogPrintf(LOG_CRITICAL,"%s : SFT length not as specified.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}



/** Compute the grid on binary parameters based on the semi-coherent metric
 *
 * We use this grid to perform the integration of the posterior and to ultimately
 * compute the Bayes factor.
 *
 */
int XLALComputeBinaryGridParams(GridParameters **binarygridparams,  /**< [out] the binary parameter grid */
				REAL8Space *space,                  /**< [in] the signal parameter space */
				REAL8 T,                            /**< [in] the duration of the observation */
				REAL8 DT,                           /**< [in] the length of the coherent segments */
				REAL8 mu,                           /**< [in] the mismatch */
				REAL8 coverage		/**< UNDOCUMENTED */
				)
{
  REAL8 gnn[NBINMAX];                    /* stores the diagonal metric elements */
  INT4 ndim = 0;                         /* the number of actual search dimensions */
  INT4 n,k;                              /* counters */

  /* validate input arguments */
  if ((*binarygridparams) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GridParameters structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (space == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input ParameterSpace structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (T < 0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input T parameter < 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((DT < 0) || (DT > T)) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input DT parameter < 0 or < T.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ( (mu < 0) || (mu>1) ) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input mismatch parameter, not in range 0 -> 1.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* compute the semi-coherent binary metric diagonal elements */
  {
    REAL8 numax = space->data[0].max;
    REAL8 amax = space->data[1].max;
    REAL8 omegamax = space->data[3].max;
    gnn[0] = (pow(LAL_PI,2.0)/6.0)*pow(DT,2.0);                                  /* nu */
    gnn[1] = (pow(LAL_PI,2.0)/6.0)*pow(numax*omegamax*DT,2.0);                   /* a */
    gnn[2] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*omegamax*DT,2.0);     /* tasc */
    gnn[3] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*DT*T,2.0)/12.0;       /* W */

    /* add eccentricity parameters at some point */
    /*  gnn->data[4] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*DT,2.0);       /\* kappa *\/ */
    /*  gnn->data[5] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*DT,2.0);       /\* eta *\/ */
  }

  /* allocate memory to the output */
  if ( ((*binarygridparams) = (GridParameters*)XLALCalloc(1,sizeof(GridParameters))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*binarygridparams)->grid = (Grid*)XLALCalloc(NBINMAX,sizeof(Grid))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*binarygridparams)->prod = (UINT4*)XLALCalloc(NBINMAX,sizeof(UINT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*binarygridparams)->ndim = NBINMAX;
  LogPrintf(LOG_DEBUG,"%s : allocated memory for the output grid parameters.\n",__func__);

  /* we need to determine the true number of searchable dimensions */
  /* we check the width of a 1-D template across each dimension span */
  for (n=0;n<NBINMAX;n++) {
    REAL8 deltax = 2.0*sqrt(mu/gnn[n]);
    if (space->data[n].span > deltax) ndim++;
  }
  LogPrintf(LOG_DEBUG,"%s : determined true dimensionality of binary space = %d.\n",__func__,ndim);

  /* Compute the grid spacing, grid start and span for each spin derivitive dimension */
  for (n=0;n<NBINMAX;n++) {

    REAL8 deltax;
    UINT4 length;
    REAL8 xmin;

    /* only if we have a non-zero span */
    if (space->data[n].span > 0) {

      /* compute spacing for this parameter given the total number of true search dimensions and the mismatch */
      deltax = 2.0*sqrt(mu/(ndim*gnn[n]));

      /* compute number of grid points in this dimension and enforce a grid centered on the middle of the parameter space */
      length = MYMAX((UINT4)ceil((space->data[n].span)/deltax),1);
      xmin = 0.5*(space->data[n].min + space->data[n].max) - 0.5*(length-1)*deltax;

    }
    else {

      /* otherwise set the space boundaries accordingly */
      deltax = 0.0;
      length = 1;
      xmin = space->data[n].min;

    }

    (*binarygridparams)->grid[n].delta = deltax;
    (*binarygridparams)->grid[n].oneoverdelta = 1.0/deltax;
    (*binarygridparams)->grid[n].length = length;
    (*binarygridparams)->grid[n].min = xmin;
    strncpy((*binarygridparams)->grid[n].name,space->data[n].name,LALNameLength*sizeof(CHAR));

    LogPrintf(LOG_DEBUG,"%s : %s -> [%e - %e] (%e) %d grid points.\n",__func__,(*binarygridparams)->grid[n].name,(*binarygridparams)->grid[n].min,
	      (*binarygridparams)->grid[n].min+(*binarygridparams)->grid[n].length*(*binarygridparams)->grid[n].delta,
	      (*binarygridparams)->grid[n].delta,(*binarygridparams)->grid[n].length);
  }
  LogPrintf(LOG_DEBUG,"%s : computed output grid parameters.\n",__func__);

  /* compute some internally required parameters for the grid */
  (*binarygridparams)->mismatch = mu;
  (*binarygridparams)->max = 1;
  for (k=0;k<(INT4)(*binarygridparams)->ndim;k++) (*binarygridparams)->max *= (*binarygridparams)->grid[k].length;

  (*binarygridparams)->prod[0] = 1;
  for (k=1;k<(INT4)(*binarygridparams)->ndim;k++) (*binarygridparams)->prod[k] = (*binarygridparams)->prod[k-1]*(*binarygridparams)->grid[k-1].length;

  /* if we've specified a random template bank */
  if (coverage>0) {

    REAL8 Vn = pow(LAL_PI,ndim/2.0)/gsl_sf_gamma(1.0+ndim/2.0);

    REAL8 G11 = pow(LAL_PI,2.0)*DT*DT/3;
    REAL8 G22 = pow(LAL_PI,2.0)*DT*DT/6;
    REAL8 G33 = pow(LAL_PI,2.0)*DT*DT/6;
    REAL8 G44 = pow(LAL_PI,2.0)*DT*DT*T*T/72;

    REAL8 dVr = sqrt(G11*G22*G33*G44);
    REAL8 Vsr = dVr*space->data[2].span*(1.0/60.0)*(pow(space->data[3].max,5.0)-pow(space->data[3].min,5.0))*(pow(space->data[0].max,4.0)-pow(space->data[0].min,4.0))*(pow(space->data[1].max,3.0)-pow(space->data[1].min,3.0));

    (*binarygridparams)->Nr = (INT4)ceil((1.0/Vn)*log(1.0/(1.0-coverage))*(pow(mu,-ndim/2.0))*Vsr);
    LogPrintf(LOG_DEBUG,"%s : computed the number of random binary templates to be %d.\n",__func__,(*binarygridparams)->Nr);
    LogPrintf(LOG_DEBUG,"%s : to br compared to the total number of cubic templates %d (%.6f).\n",__func__,(*binarygridparams)->max,(REAL8)(*binarygridparams)->max/(REAL8)(*binarygridparams)->Nr);

  }
  else (*binarygridparams)->Nr = -1;

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Append the given SFTtype to the SFT-vector (no SFT-specific checks are done!) */
int XLALAppendSFT2Vector (SFTVector *vect,		/**< destinatino SFTVector to append to */
			  const SFTtype *sft            /**< the SFT to append */
			  )
{
  UINT4 oldlen = vect->length;

  if ( (vect->data = LALRealloc ( vect->data, (oldlen + 1)*sizeof( *vect->data ) )) == NULL ) {
     LogPrintf(LOG_CRITICAL,"%s: Error, unable to allocate memory\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
  }
  memset ( &(vect->data[oldlen]), 0, sizeof( vect->data[0] ) );
  vect->length ++;

  XLALCopySFT(&vect->data[oldlen], sft );

  return XLAL_SUCCESS;

} /* XLALAppendSFT2Vector() */

/** Convert an input binary file into an SFTVector
 *
 *
 */
int XLALBinaryToSFTVector(SFTVector **SFTvect, 	   /**< [out] copied SFT (needs to be allocated already) */
			  CHAR *filename,          /**< [in] the input filename */
			  LIGOTimeGPS *fileStart,  /**< [in] the file start time */
 			  BinaryToSFTparams *par,  /**< [in] the required parameters  */
			  INT8Vector **np,         /**< [out] the number of photons in each SFT */
			  REAL8Vector **R          /**< [out] the estimated photon rate per bin */
			 )
{

  INT4 i,k;
  static const LALUnit empty_LALUnit;
  fprintf(stdout,"---> working on file %s\n",filename);

  /* initialise results vector */
  if ((*SFTvect)==NULL) {
    if (((*SFTvect) = XLALCalloc(1,sizeof(SFTVector)))==NULL) {
      LogPrintf(LOG_CRITICAL,"%s: Unable to allocate memory for SFTVector\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
  }

  /**********************************************************************************/
  /* open the binary file and count the number of elements */
  FILE *binfp = NULL;                              /* file pointer for input file */
  REAL4 dummy;                                  /* dummy variable */
  if ((binfp = fopen(filename,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,filename);
    return 1;
  }
  i = 0;
  while (fread(&dummy,sizeof(REAL4),1,binfp)) i++;
  INT8 N = i;
  fclose(binfp);
  fprintf(stdout,"---> %s : counted %" LAL_UINT8_FORMAT " samples in the file.\n",__func__,N);

  /* return if file has no length */
  if (N==0) {
    LogPrintf(LOG_CRITICAL,"%s : found no time samples in file %s.\n",__func__,filename);
    return XLAL_SUCCESS;
  }

  /* open the file for reading */
  if ((binfp = fopen(filename,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,filename);
    return 1;
  }

  /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
  REAL8 dt = par->tsamp;
  LIGOTimeGPS startTimeGPS, endTimeGPS;                     /* the start and end time of the observation */
  fprintf(stdout,"---> the filestart time is %d %d\n",fileStart->gpsSeconds,fileStart->gpsNanoSeconds);
  memcpy(&startTimeGPS,fileStart,sizeof(LIGOTimeGPS));
  memcpy(&endTimeGPS,&startTimeGPS,sizeof(LIGOTimeGPS));
  XLALGPSAdd(&endTimeGPS,N*dt);
  fprintf(stdout,"---> input binary file has start and end [%d %d - %d %d]\n",startTimeGPS.gpsSeconds,startTimeGPS.gpsNanoSeconds,endTimeGPS.gpsSeconds,endTimeGPS.gpsNanoSeconds);
  fprintf(stdout,"---> input binary file has length %.12f sec\n",N*dt);
  REAL4TimeSeries *Tseries = XLALCreateREAL4TimeSeries ( "X1", &(startTimeGPS), 0, dt, &empty_LALUnit, N);
  fprintf(stdout,"---> made timeseries with length %f sec\n",Tseries->deltaT*Tseries->data->length);
  INT8 sum = 0;

  /* read in the data to the timeseries - max of MAXNTSERIES values */
  for (i=0;i<N;i++) {
    fread(&dummy,sizeof(REAL4),1,binfp);
    Tseries->data->data[i] = (REAL4)dummy;
    sum += (INT8)Tseries->data->data[i];
  }

  /* check for zero data */
  if (sum == 0) {
    LogPrintf(LOG_CRITICAL,"%s : found no photons in file %s.\n",__func__,filename);
    return XLAL_SUCCESS;
  }
  fprintf(stdout,"---> the total number of photons in the binary file = %" LAL_INT8_FORMAT "\n",sum);

  /* inject signal if requested */
  if (par->amp_inj>0) {
    fprintf(stdout,"%s : injecting signal into the data with fractional amplitude %f.\n",__func__,par->amp_inj);
    REAL8 Om = LAL_TWOPI/par->P_inj;
    REAL8 bg = (REAL8)sum/(REAL8)N;
    for (i=0;i<N;i++) {
      LIGOTimeGPS t;
      memcpy(&t,&startTimeGPS,sizeof(LIGOTimeGPS));
      XLALGPSAdd(&t,dt*i);                              /* the current sample GPS time */
      REAL8 tmtref = XLALGPSDiff(&t,&(par->tref));         /* t - t_ref */
      REAL8 tmtasc = XLALGPSDiff(&t,&(par->tasc_inj));     /* t - t_asc */
      REAL8 phase = par->phi_inj + LAL_TWOPI*par->f_inj*(tmtref - par->asini_inj*sin(Om*tmtasc));
      REAL8 x = bg*(1.0 + par->amp_inj*sin(phase));
      Tseries->data->data[i] = (REAL4)gsl_ran_poisson(par->r,x);
      /* fprintf(stdout,"%d %d %.12f %.12f %.12f %.12f %d %.2f\n",t.gpsSeconds,t.gpsNanoSeconds,tmtref,tmtasc,phase,x,y,Tseries->data->data[i]); */
    }
  }

  /**********************************************************************************/
  /* make timestamps */
  LIGOTimeGPSVector timestamps;
  fprintf(stdout,"---> requested starttime = %d %d endtime = %d %d\n",par->tstart.gpsSeconds,par->tstart.gpsNanoSeconds,endTimeGPS.gpsSeconds,endTimeGPS.gpsNanoSeconds);
  INT4 tspan = (INT4)XLALGPSDiff(&endTimeGPS,&(par->tstart));
  fprintf(stdout,"---> requested time span = %d sec\n",tspan);
  if (tspan<0) {
    LogPrintf(LOG_CRITICAL,"%s : requested sft start time after end of data!\n",__func__);
    return XLAL_SUCCESS;
  }
  INT4 maxsft = (INT4)floor((REAL8)tspan/par->tsft);        /* compute the max number of SFTs */
  timestamps.length = maxsft;
  fprintf(stdout,"---> initial estimate of %d SFTs will be made\n",maxsft);
  *np = XLALCreateINT8Vector(maxsft);

  /* if we have one or more SFTs to make */
  if (maxsft>0) {

    INT4 cnt = 0;
    INT4 Nsft = 0;
    INT4 Ndt = (INT4)(par->tsft/par->tsamp);
    fprintf(stdout,"---> number of samples per SFT = %d\n",Ndt);
    timestamps.data = XLALCalloc(maxsft,sizeof(LIGOTimeGPS));
    for (i=0;i<maxsft;i++) {

      /* test if the SFT has photons */
      sum = 0;
      INT4 j;
      for (j=0;j<Ndt;j++) sum += (INT8)Tseries->data->data[cnt++];
      fprintf(stdout,"---> the photon sum for the current SFT is %" LAL_INT8_FORMAT "\n",sum);
      if (sum > 0) {
        memcpy(&(timestamps.data[Nsft]),&(par->tstart),sizeof(LIGOTimeGPS));
        fprintf(stdout,"---> about to add %f to %d %d\n",(REAL8)(i*par->tsft),timestamps.data[Nsft].gpsSeconds,timestamps.data[Nsft].gpsNanoSeconds);
        XLALGPSAdd(&(timestamps.data[Nsft]),(REAL8)(i*par->tsft));
	(*np)->data[Nsft] = sum;
        fprintf(stdout,"np->data[%d] = %" LAL_INT8_FORMAT "\n",Nsft,(*np)->data[Nsft]);
        Nsft++;
      }

    }

    /* if after testing the individual SFT photon counts we still have valid SFTs */
    if (Nsft>0) {

      timestamps.data = LALRealloc(timestamps.data,Nsft*sizeof(LIGOTimeGPS));
      timestamps.length = Nsft;
      (*np)->data = XLALRealloc((*np)->data,Nsft*sizeof(INT8));
      (*np)->length = Nsft;
      *R = XLALCreateREAL8Vector(Nsft);
      fprintf(stdout,"---> start SFT [%d %d - %d %d]\n",timestamps.data[0].gpsSeconds,timestamps.data[0].gpsNanoSeconds,timestamps.data[0].gpsSeconds+par->tsft,timestamps.data[0].gpsNanoSeconds);
      fprintf(stdout,"---> end SFT [%d %d - %d %d]\n",timestamps.data[Nsft-1].gpsSeconds,timestamps.data[Nsft-1].gpsNanoSeconds,timestamps.data[Nsft-1].gpsSeconds+par->tsft,timestamps.data[Nsft-1].gpsNanoSeconds);

      /**********************************************************************************/
      /* define the high pass filter params and high pass filter the data */
      PassBandParamStruc filterpar;
      char tmpname[] = "Butterworth High Pass";
      filterpar.name  = tmpname;
      filterpar.nMax  = 10;
      filterpar.f2    = par->highpassf;
      filterpar.a2    = 0.5;
      filterpar.f1    = -1.0;
      filterpar.a1    = -1.0;
      XLALButterworthREAL4TimeSeries(Tseries, &filterpar);
      fprintf(stdout,"---> Filtered time-series now has length %f sec\n",Tseries->deltaT*Tseries->data->length);

      /**********************************************************************************/
      /* compute SFTs from timeseries */
      SFTParams XLAL_INIT_DECL(sftParams);
      sftParams.Tsft = par->tsft;
      sftParams.timestamps = &timestamps;
      sftParams.noiseSFTs = NULL;       // not used here any more!
      sftParams.window = NULL;
      SFTVector *sftvect = NULL;
      fprintf(stdout,"---> timeseries prior to making SFTs has length %f sec\n",Tseries->deltaT*Tseries->data->length);
      XLAL_CHECK ( (sftvect = XLALSignalToSFTs (Tseries, &sftParams)) != NULL, XLAL_EFUNC );

      /**********************************************************************************/
      /* extract effective band from this, if neccessary (ie if faster-sampled output SFTs) */
      SFTVector *tempSFTvect = XLALExtractBandFromSFTVector ( sftvect, par->freq, par->freqband-1.0/par->tsft );    /* the last bin has to be avoided */
      XLALDestroySFTVector(sftvect);

      /**********************************************************************************/
      /* append these SFTs to the full list of SFTs */
      for (k=0;k<Nsft;k++) {
        XLALAppendSFT2Vector((*SFTvect),&(tempSFTvect->data[k]));
        XLALNormalizeSFTMedian(&(tempSFTvect->data[k]),&((*R)->data[k]),0);
      }
      fprintf(stdout,"---> appended SFTs\n");

      /**********************************************************************************/

      XLALDestroySFTVector(tempSFTvect);

    }

    XLALFree(timestamps.data);
    XLALDestroyREAL4TimeSeries(Tseries);

  }

  fclose(binfp);
  fprintf(stdout,"---> completed analysis of file %s\n",filename);

  return XLAL_SUCCESS;

}

/**
  * Normalize an sft based on median estimate PSD
  */
  int
  XLALNormalizeSFTMedian (SFTtype *sft,           /**< SFT to be normalized */
                          REAL8 *mean,            /**< [out] the mean computed from the median */
                          INT4 flag               /**< [in] flag for normalising, 0 = no normalising */
			)
  {
    /* check input argments */
    XLAL_CHECK (sft && sft->data && sft->data->data && sft->data->length > 0, XLAL_EINVAL, "Invalid NULL or zero-length input in 'sft'" );
    UINT4 length = sft->data->length;

    /* allocate memory for power */
    gsl_vector *P = gsl_vector_alloc(length);
    REAL8 *s = XLALCalloc(length,sizeof(REAL8));

    /* compute power */
    for (UINT4 j = 0; j < length; j++) {
      gsl_vector_set(P,j,crealf(sft->data->data[j])*crealf(sft->data->data[j]) + cimagf(sft->data->data[j])*cimagf(sft->data->data[j]));
    }

    /* sort the power vector */
    gsl_sort_vector(P);
    for (UINT4 j = 0; j < length; j++) s[j] = gsl_vector_get(P,j);

    /* compute the median and the median bias */
    REAL8 med = gsl_stats_median_from_sorted_data (s,1,length);
    REAL8 medianBias = XLALRngMedBias (length);
    *mean = med/medianBias;

    /* free memory */
    gsl_vector_free(P);
    XLALFree(s);

    /* loop over sft and normalize by the median */
    if (flag) {
      for (UINT4 j = 0; j < length; j++) {
        sft->data->data[j] *= (REAL4)(sqrt(medianBias/med));
      } // for j < length
    }
    return XLAL_SUCCESS;


  } /* XLALNormalizeSFTMedian() */


 /**
  * Function for normalizing a vector of SFTs.
  */
  int
  XLALNormalizeSFTVectMedian ( SFTVector  *sftVect,            /**< [in/out] pointer to a vector of SFTs which will be normalized */
                               REAL8Vector *mean,              /**< [out] a vector of mean estimates based on the median */
			       INT4 flag                       /**< [in] flag for normalising, 0 = no normalising */
			     )
  {
    /* check input argments */
    XLAL_CHECK ( sftVect && sftVect->data && sftVect->length > 0, XLAL_EINVAL, "Invalid NULL or zero-length input in 'sftVect'");
    XLAL_CHECK ( mean && mean->length > 0, XLAL_EINVAL, "Invalid NULL or zero-length input in 'mean'");

    /* loop over sfts and normalize them */
    for (UINT4 j = 0; j < sftVect->length; j++)
      {
        SFTtype *sft = &sftVect->data[j];

        /* call sft normalization function */
        XLAL_CHECK ( XLALNormalizeSFTMedian ( sft, &(mean->data[j]), flag ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALNormalizeSFTMedian() failed." );

      } /* for j < sftVect->length */

    return XLAL_SUCCESS;

  } /* XLALNormalizeSFTVectMedian() */


/**
  * Normalize an sft based on mean estimate PSD
  */
  int
  XLALNormalizeSFTMean (SFTtype *sft,          /**< SFT to be normalized */
                        REAL8 *mean,           /**< [out] the mean power */
			INT4 flag              /**< [in] flag for normalising, 0 = no normalising */
		       )
  {
    /* check input argments */
    XLAL_CHECK (sft && sft->data && sft->data->data && sft->data->length > 0, XLAL_EINVAL, "Invalid NULL or zero-length input in 'sft'" );
    UINT4 length = sft->data->length;

    /* find mean psd */
    REAL4 sum = 0.0;
    for (UINT4 j = 0; j < length; j++) {
      sum += crealf(sft->data->data[j])*crealf(sft->data->data[j]) + cimagf(sft->data->data[j])*cimagf(sft->data->data[j]);
    }
    sum /= (REAL4)length;
    *mean = sum;

    /* loop over sft and normalize */
    if (flag) {
      for (UINT4 j = 0; j < length; j++) {
        sft->data->data[j] /= ((REAL4)sqrt(sum));
      } // for j < length
    }
    return XLAL_SUCCESS;

  } /* XLALNormalizeSFTMean() */


 /**
  * Function for normalizing a vector of SFTs.
  */
  int
  XLALNormalizeSFTVectMean ( SFTVector  *sftVect,            /**< [in/out] pointer to a vector of SFTs which will be normalized */
                             REAL8Vector *mean,              /**< [out] a vector of means estimating the noise floor */
 		             INT4 flag                       /**< [in] flag for normalising, 0 = no normalising */
                           )
  {
    /* check input argments */
    XLAL_CHECK ( sftVect && sftVect->data && sftVect->length > 0, XLAL_EINVAL, "Invalid NULL or zero-length input in 'sftVect'");
    XLAL_CHECK ( mean && mean->length > 0, XLAL_EINVAL, "Invalid NULL or zero-length input in 'mean'");

    /* loop over sfts and normalize them */
    for (UINT4 j = 0; j < sftVect->length; j++)
      {
        SFTtype *sft = &sftVect->data[j];

        /* call sft normalization function */
        XLAL_CHECK ( XLALNormalizeSFTMean ( sft, &(mean->data[j]), flag ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALNormalizeSFT() failed." );

      } /* for j < sftVect->length */

    return XLAL_SUCCESS;

  } /* XLALNormalizeSFTVectMean() */
