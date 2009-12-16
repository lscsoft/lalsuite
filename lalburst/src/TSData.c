 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
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

/*-----------------------------------------------------------------------
 *
 * File Name: TSData.c
 *
 * Author: Torres Cristina V. (LLO)
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <lal/TSData.h>
#include <lal/TSSearch.h>
#include <lal/LALStdlib.h>

/* macro to "use" unused function parameters */
#define UNUSED(expr) do { (void)(expr); } while(0)

NRCSID (TSDATAC,"$Id$");

/*
 * Local static functions
 */
static REAL4 localGauss2(REAL4 lineWidth)
{
  REAL8 limit1,limit2; /* the limits over which to average the function */
  INT4 numberInt=20;   /* number of bins;*/
  REAL8 interval;      /* the bin width */
  REAL8 loop;          /* the loop variable */
  INT4 i;              /* loop variable */
  REAL8 sum;           /* the sum */
  REAL4 result;        /* the value to return */
  REAL8 norm;          /* normalizing constant */
  REAL8 temp;          /* temporary variable */
  REAL8 sigma=0;

  sigma = lineWidth/(2*sqrt(3));
  limit1 = -(lineWidth/2);
  limit2 =  (lineWidth/2);
  interval= (limit2-limit1)/numberInt;
  sum=0.0;
  norm = 1.0/(sigma*sigma*sigma*sqrt(LAL_PI*2.0));
  for(i=0;i<=numberInt;i++){
    loop = limit1 + i*interval;
    temp = ((loop*loop)/(sigma*sigma));
    sum += norm * exp(-0.5*temp) * (temp - 1.0);
  }
  sum /= (numberInt+1);
  result = sum;
  return result;
}
/*End local gauss func*/

void
LALTracksearchFindLambdaMean(
			     LALStatus                *status,
			     TimeFreqRep               map,
			     TSSearchParams           *searchParams
			     )
{
  INT4            i=0;/* Counter in F axis */
  INT4            j=0;/* Counter in T axis */
  REAL8        sumX=0;/* Sum of h differences */
  REAL8     sumXSqr=0;/* Sum of h^2 differences */
  INT8      counter=0;/* Number of h diff measurements taken */
  REAL8     current=0;/* Current h value */
  REAL8       meanH=0;/* Mean feature step height */
  REAL8        stdH=0;/* Variance on the mean step height */
  REAL8 upperThresh=0;/* Auto discovered upper curvatuve threshold */
  REAL8 lowerThresh=0;/* Value set from upperThresh */
  REAL8 myGaussian=0;
#if 0
  REAL8 myFloor=10e20;/*Lowest value in TFR */
#endif
  INT4        lowerFRow=-1;
  INT4        upperFRow=-1;
  REAL8       binPerHz=0;
  INITSTATUS(status,"LALTracksearchFindLambdaMean", TSDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT(searchParams,status, TSDATA_ENNUL, TSDATA_MSGENNUL);
  counter=0;
  sumX=0;
  sumXSqr=0;
  current=0;
  /*
   * If there was band passing determine regions of F to
   * calculate the Lambda values inside of band passes region
   */
  binPerHz=((map.fRow/2+1)/(searchParams->SamplingRate/2.0));
  if (searchParams->highPass > 0)
    {
      lowerFRow=floor(binPerHz*searchParams->highPass);
    }
  else
    {
      lowerFRow=0;
    }
  if (searchParams->lowPass > 0)
    {
      upperFRow=ceil(binPerHz*searchParams->lowPass);
    }
  else
    {
      upperFRow=map.fRow/2+1;
    }

  for (i = 0;i < map.tCol;i++)
    {
      for (j = 0;j <(map.fRow/2+1);j++)
	{
	  current=0;
	  current=map.map[i][j];
	  sumX=sumX+current;
	  sumXSqr=sumXSqr+(current*current);
          counter=counter+1;
	};
    };
  /* Determine mean H value */
  meanH=(sumX/counter);
  /* Determine STD of H value */
  stdH=sqrt((sumXSqr-(sumX*meanH))/(counter-1));
  myGaussian=localGauss2(searchParams->LineWidth);
  /* Given desired contrast Z score */
  upperThresh=(meanH+(searchParams->autoThresh*stdH))*myGaussian;
  lowerThresh=upperThresh*searchParams->relaThresh;
  searchParams->StartThresh=fabs(upperThresh);
  searchParams->LinePThresh=fabs(lowerThresh);
  if (searchParams->verbosity >= verbose)
    {
      fprintf(stdout,"Auto lambda invoked\n");
      fprintf(stdout,"Lh %e \t Ll %e \n 2nd D Gauss %f \n Mean h:  %10.20e \n Std  h: %10.20e \n",searchParams->StartThresh,searchParams->LinePThresh,myGaussian,meanH,stdH);
    }
  /* Need to throw an error if we get a value Lh <= 0 */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}/*
  * End LALTracksearchFindLambda
  */

static int lambdaSortOpt(const void *a,const void *b)
{
  const REAL8 *A = ((const REAL8 *) a);
  const REAL8 *B = ((const REAL8 *) b);

  if (*A<*B)
    return -1;
  else if (*A==*B)
    return 0;
  else
    return 1;
}

void
LALTracksearchFindLambdaMedian(
			       LALStatus                *status,
			       TimeFreqRep               map,
			       TSSearchParams           *searchParams
			       )
{
  REAL8       medianH=0;/*Median value calculated*/
  REAL8      *vector=NULL;/*Pointer to array of REAL8s for sorting*/
  REAL8       upperThresh=0;
  REAL8       lowerThresh=0;
  REAL8       myGaussian=0;
  REAL8       binPerHz=0;
  INT4        count=0;/*Count the number of elements inserted into
			vector*/
  INT4        i=0;/*Index*/
  INT4        j=0;/*Index*/
  INT4        k=0;/*Index*/
  INT4        lowerFRow=-1;
  INT4        upperFRow=-1;

  INITSTATUS(status,"LALTracksearchFindLambdaMedian", TSDATAC);
  ATTATCHSTATUSPTR (status);
  /*
   * If there was band passing determine regions of F to
   * calculate the Lambda values inside of band passes region
   */
  binPerHz=((map.fRow/2+1)/(searchParams->SamplingRate/2.0));
  if (searchParams->highPass > 0)
    {
      lowerFRow=floor(binPerHz*searchParams->highPass);
    }
  else
    {
      lowerFRow=0;
    }
  if (searchParams->lowPass > 0)
    {
      upperFRow=ceil(binPerHz*searchParams->lowPass);
    }
  else
    {
      upperFRow=map.fRow/2+1;
    }

  /*Read out pixel h from map into vector for sorting*/
  count=((map.tCol)*(map.fRow/2+1));
  count=(map.tCol*(upperFRow-lowerFRow));
  if (count < map.tCol)
    ABORTXLAL(status);

  vector=(REAL8*)LALMalloc(count*sizeof(REAL8));
  for (i = 0;i < map.tCol;i++)
    {
      for (j = lowerFRow;j <(upperFRow);j++)
	{
	  /*fprintf(stdout,"%i %i %i %i %f\n",i,j,k,count,map.map[i][j]);*/
	  vector[k]=map.map[i][j];
	  k=k+1;
	}
    }
  qsort(vector,count,sizeof(REAL8),lambdaSortOpt);
  if (k%2 == 0)
    medianH=(vector[k/2]+vector[(k/2)+1])/2;
  else
    medianH=(vector[(k+1)/2]);
  LALFree(vector);
  /* Given desired contrast Z score */
  /*Valid interval of autoThresh(x) are {x|-1<x<=\inf} */
  if (searchParams->autoThresh <= -1)
    searchParams->autoThresh = -0.99999999;
  myGaussian=localGauss2(searchParams->LineWidth);
  upperThresh=myGaussian*medianH*(1+searchParams->autoThresh);
  lowerThresh=upperThresh*searchParams->relaThresh;
  searchParams->StartThresh=fabs(upperThresh);
  searchParams->LinePThresh=fabs(lowerThresh);
  if (searchParams->verbosity >= verbose)
    {
      fprintf(stdout,"Auto lambda invoked\n");
      fprintf(stdout,"Lh %e \t Ll %e \n 2nd D Gauss %f \n Median h: %10.20e \n Pixel Count: %i \n",searchParams->StartThresh,
	      searchParams->LinePThresh,
	      myGaussian,
	      medianH,
	      k);
    }
  /* Need to throw an error if we get a value Lh <= 0 */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}/*
  * End LALTracksearchFindLambda
  */


void
LALCreateTSDataSegmentVector (
			      LALStatus                    *status,
			      TSSegmentVector             **vector,
			      TSCreateParams               *params
			      )
{
  INT4                           i;
  UINT4                  segmentLength=0;
  const LIGOTimeGPS      gps_zero = LIGOTIMEGPSZERO;

  INITSTATUS (status, "LALCreateTSSegmentVector", TSDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT (!*vector, status, TSDATA_ENNUL, TSDATA_MSGENNUL);
  ASSERT (params, status, TSDATA_ENULL, TSDATA_MSGENULL);
  ASSERT (params->numberDataSegments > 0,
	  status, TSDATA_ESEGZ, TSDATA_MSGESEGZ);
  ASSERT (params->dataSegmentPoints > 0,
	  status, TSDATA_ENUMZ, TSDATA_MSGENUMZ);

  *vector = (TSSegmentVector *) LALMalloc(sizeof(TSSegmentVector));


  (*vector)->length = params->numberDataSegments;
  (*vector)->SegBufferPoints=params->SegBufferPoints;
  (*vector)->dataSeg  = (REAL4TimeSeries **)
    LALMalloc(params->numberDataSegments*sizeof(REAL4TimeSeries*));
  if ( !((*vector)->dataSeg) )
    {
      LALFree(*vector);
      ABORT( status, TSDATA_EALOC, TSDATA_MSGEALOC );
    }
  /*
   * Intialize structure to empty state include space for
   * SegBufferPoints also.
   */
  segmentLength=params->dataSegmentPoints+(2*params->SegBufferPoints);
  for (i = 0; i < (INT4)((*vector)->length); i++)
    {
      (*vector)->dataSeg[i] = XLALCreateREAL4TimeSeries("Uninitialized",
			       &gps_zero,
			       0,
			       1,
			       &lalDimensionlessUnit,
			       segmentLength);
      if(!(*vector)->dataSeg[i])
        ABORTXLAL (status);
    }
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
/*END LALCreateTSDataSegmentVector*/


void
LALDestroyTSDataSegmentVector (
			       LALStatus                  *status,
			       TSSegmentVector            *vector
			       )
{
  INT4                    i;

  INITSTATUS (status, "LALDestroyTSDataSegmentVector", TSDATAC);
  ATTATCHSTATUSPTR (status);
  ASSERT (vector,
	  status,
	  TSDATA_ENULL,
	  TSDATA_MSGENULL);
  for (i = 0; i < (INT4)(vector->length) ; i++)
      XLALDestroyREAL4TimeSeries(vector->dataSeg[i]);
  if (vector->dataSeg)
    LALFree(vector->dataSeg);
  if (vector)
    LALFree(vector);
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
/*END LALDestroyTSDataSegmentVector */


void
LALTrackSearchConnectSigma(
			   LALStatus                     *status,
			   TrackSearchOut                *curveinfo,
			   TimeFreqRep                    map,
			   TrackSearchParams              params
			   )
{
  INT4 i,j;
  Curve *curveA;
  Curve *curveB;

  INITSTATUS(status,"LALTrackSearchConnectSigma",TSDATAC);
  ATTATCHSTATUSPTR (status);


  /* Outer Loop trys to establish joints */
  for (i=0;i<curveinfo->numberOfCurves;i++)
    {
      /*Ignore curves labeled with marker 'D' */
      if (curveinfo->curves[i].trash != 'D' )
	{
	  for (j=0;j<curveinfo->numberOfCurves;j++)
	    {
	      /*Ignore 'D' label curves*/
	      if (curveinfo->curves[j].trash != 'D' )
		{
		  /*Inside Loop test this track to all other*/
		  if (i!=j)
		    {
		      curveA=&(curveinfo->curves[i]);
		      curveB=&(curveinfo->curves[j]);
		      if ((abs(curveA->row[curveA->n]-curveB->row[0]) <
			   params.sigma+1)
			  &&
			  (abs(curveA->col[curveA->n]-curveB->row[0]) < params.sigma+1))
			{
			  /* Following function is iffy */
			  /* Ok they touch join them */
			  connect2Segments(map,curveA,curveB);
			  /* Reset interior loop to try and check new longer track*/
			  j=0;
			}
		    }
		}
	    }
	}
    }
  /* Done joing useful curve candidates */
  /* Cleanup linked list */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALTrackSearchApplyThreshold(
			     LALStatus         *status,
			     TrackSearchOut    *curveinfo,
			     TrackSearchOut    *dataProduct,
			     TSSearchParams     params
			     )
{
  INT4 UsefulCurves;
  INT4 i;
  INT4 j;
  INT4 cicn;

  INITSTATUS(status,"LALTrackSearchApplyThreshold",TSDATAC);
  ATTATCHSTATUSPTR (status);
  ASSERT(dataProduct->curves==NULL,status,TS_NON_NULL_POINTER,TS_MSGNON_NULL_POINTER);
  /*Trackstore struct field not copied! */
  UsefulCurves=0; /*Array numbering starts with zero */
  dataProduct->numberOfCurves=0;/*Default value*/
  if (curveinfo->numberOfCurves > 0)
    {
      /* Count up number of useful curves to keep */
      for(i=0;i<curveinfo->numberOfCurves;i++)
	{
	  if ((curveinfo->curves[i].n >= params.MinLength) &&
	      (curveinfo->curves[i].totalPower >= params.MinPower) &&
	      (curveinfo->curves[i].snrEstimate >= params.MinSNR))
	    {
	      /*UsefulCurve Var is not count but curve index so add 1*/
	      /* Expand the structure to take in another curve */
	      UsefulCurves=dataProduct->numberOfCurves;
	      dataProduct->curves = (Curve*)LALRealloc(dataProduct->curves,(sizeof(Curve) * (UsefulCurves+1)));

	      /*Length of curve to copy*/
	      cicn=curveinfo->curves[i].n;
	      dataProduct->curves[UsefulCurves].row=(INT4*)LALMalloc(sizeof(INT4) * cicn);
	      dataProduct->curves[UsefulCurves].col=(INT4*)LALMalloc(sizeof(INT4) * cicn);
	      dataProduct->curves[UsefulCurves].depth=(REAL4*)LALMalloc(sizeof(REAL4) * cicn);
	      dataProduct->curves[UsefulCurves].fBinHz=(REAL4*)LALMalloc(sizeof(REAL4) * cicn);
	      dataProduct->curves[UsefulCurves].gpsStamp=(LIGOTimeGPS*)LALMalloc(sizeof(LIGOTimeGPS) * cicn);
	      dataProduct->curves[UsefulCurves].n = cicn;
	      dataProduct->curves[UsefulCurves].totalPower = curveinfo->curves[i].totalPower;
	      dataProduct->curves[UsefulCurves].snrEstimate = curveinfo->curves[i].snrEstimate;
	      /* Copy the data over */
	      for(j=0;j<dataProduct->curves[UsefulCurves].n;j++)
		{
		  dataProduct->curves[UsefulCurves].row[j]=curveinfo->curves[i].row[j];
		  dataProduct->curves[UsefulCurves].col[j]=curveinfo->curves[i].col[j];
		  dataProduct->curves[UsefulCurves].depth[j]=
		    curveinfo->curves[i].depth[j];
		  dataProduct->curves[UsefulCurves].fBinHz[j]=
		    curveinfo->curves[i].fBinHz[j];
		  dataProduct->curves[UsefulCurves].gpsStamp[j].gpsSeconds=
		    curveinfo->curves[i].gpsStamp[j].gpsSeconds;
		  dataProduct->curves[UsefulCurves].gpsStamp[j].gpsNanoSeconds=
		    curveinfo->curves[i].gpsStamp[j].gpsNanoSeconds;
		}
	      /*Increase index recall count is +1 of index*/
	      dataProduct->numberOfCurves++;
	      /* Move on to check next curve */
	    }
	}
      /*If nothing passes the threshold deallocate the pointer*/
      if ((dataProduct->numberOfCurves == 0) &&
	  (dataProduct->curves != NULL))
	{
	  LALFree(dataProduct->curves);
	  dataProduct->curves=NULL;
	}
    }
  /* Done moving useful curve candidates */
  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* End LALTrackSearchApplyThreshold */


/* Routine to whiten the data stream */
/* We will overwhiten by calling 2 times in order */
/* Currently unused and not thoroughly tested */
void
LALTrackSearchWhitenREAL4TimeSeries(
				    LALStatus              *status,
				    REAL4TimeSeries        *signalvec,
				    REAL4FrequencySeries   *signalPSD,
				    TSWhitenParams          params
				    )
{
  UINT4                      i=0;
  RealFFTPlan               *forwardPlan=NULL;
  RealFFTPlan               *reversePlan=NULL;
  COMPLEX8FrequencySeries   *signalFFT=NULL;
  const LIGOTimeGPS              gps_zero = LIGOTIMEGPSZERO;
  UINT4                      planLength=0;
  REAL8                     factor=0;
  LALUnit                   tmpUnit1=lalDimensionlessUnit;
  LALUnitPair               tmpUnitPair;
  RAT4                      exponent;
  INITSTATUS(status,"LALTrackSearchWhitenREAL4TimeSeries",TSDATAC);
  ATTATCHSTATUSPTR (status);

  /* params is unused in this function */
  UNUSED(params);

  /*Setup FFT Plans*/
  planLength=signalvec->data->length;
  LALCreateForwardREAL4FFTPlan(status->statusPtr,
			       &forwardPlan,
			       planLength,
			       0);
  CHECKSTATUSPTR (status);
  LALCreateReverseREAL4FFTPlan(status->statusPtr,
			       &reversePlan,
			       planLength,
			       0);
  CHECKSTATUSPTR (status);

  /* Allocate space for FFT */
  signalFFT = XLALCreateCOMPLEX8FrequencySeries("tmpSegPSD",
				   &gps_zero,
				   0,
				   1/(signalvec->deltaT*signalvec->data->length),
				   &lalDimensionlessUnit,
				   planLength/2+1);
  /* FFT the time series */
  LALForwardREAL4FFT(status->statusPtr,
		     signalFFT->data,
		     signalvec->data,
		     forwardPlan);
  CHECKSTATUSPTR (status);

  /*
   * Perform whitening
   * Look at Tech Doc T010095-00  Sec3
   */
  for (i=0;i<signalFFT->data->length;i++)
    {
      if ( signalPSD->data->data[i] == 0.0 )
	factor=0;
      else
	/*
	 * This whitening filter pulled from EPsearch code
	 */
	factor=2*sqrt(signalFFT->deltaF/signalPSD->data->data[i]);

      signalFFT->data->data[i].re = signalFFT->data->data[i].re * factor;
      signalFFT->data->data[i].im = signalFFT->data->data[i].im * factor;
    }
  /*
   * Manipulate the LALUnits structure to reflect above operation
   */
  exponent.numerator=-1;
  exponent.denominatorMinusOne=1;/*2*/
  LALUnitRaise(status->statusPtr,
	       &tmpUnit1,
	       &(signalPSD->sampleUnits),
	       &exponent);
  CHECKSTATUSPTR (status);
  tmpUnitPair.unitOne=&tmpUnit1;
  tmpUnitPair.unitTwo=&(signalFFT->sampleUnits);
  LALUnitMultiply(status->statusPtr,
		  &(signalFFT->sampleUnits),
		  &tmpUnitPair);
  CHECKSTATUSPTR (status);
#if 0
  /*
   * Diagnostic code
   */
  print_complex8fseries(signalFFT,"dataFFTComplexPOST.txt");
#endif

  /*
   * Transform back to time domain
   */
  LALReverseREAL4FFT(status->statusPtr,
		     signalvec->data,
		     signalFFT->data,
		     reversePlan);
  CHECKSTATUSPTR (status);
  /*
   * The 1/n factor need to be applied
   * See lsd-5 p 259 10.1
   */
  for (i=0;i<signalvec->data->length;i++)
    signalvec->data->data[i]= signalvec->data->data[i]/signalvec->data->length;

#if 0
  /*
   * Diagnostic code
   */
  print_real4tseries(signalvec,"dataSegi.txt");
#endif

  /*
   *Release the temporary memory
   */
  XLALDestroyCOMPLEX8FrequencySeries(signalFFT);
  if (forwardPlan)
    {
      LALDestroyREAL4FFTPlan(status->statusPtr,&forwardPlan);
      CHECKSTATUSPTR (status);
    }
  if (reversePlan)
    {
      LALDestroyREAL4FFTPlan(status->statusPtr,&reversePlan);
      CHECKSTATUSPTR (status);
    }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* End whiten routine */


/* Begin Fourier Domain whitening routine */
void
LALTrackSearchWhitenCOMPLEX8FrequencySeries(
					    LALStatus                *status,
					    COMPLEX8FrequencySeries  *fSeries,
					    REAL4FrequencySeries     *PSD,
					    UINT4                     level
					    )
{
  UINT4         i=0;
  REAL8         factor=0;
  LALUnit       tmpUnit1=lalDimensionlessUnit;
  LALUnitPair   tmpUnitPair;
  RAT4          exponent;

  INITSTATUS(status,"LALTrackSearchWhitenCOMPLEX8FrequencySeries",TSDATAC);
  ATTATCHSTATUSPTR (status);
  /*
   * Error checking
   */

  /*
   * Need to add error check to see that deltaF for PSD matches
   * fSeries deltaF with TOL if not issue warning and match the PSD
   * to the fSeries otherwise throw error.
   * I.E.  Run if fSeries.deltaF>=PSD.deltaF matching as needed!
   */
  ASSERT(level > 0, status,TSDATA_EINVA,TSDATA_MSGEINVA);

  for (i=0;i<fSeries->data->length;i++)
    {
      if ( PSD->data->data[i] == 0.0 )
	factor=0;
      else
	/*
	 * This whitening filter pulled from EPsearch code
	 */
	factor=2*sqrt(fSeries->deltaF/PSD->data->data[i]);

      fSeries->data->data[i].re = fSeries->data->data[i].re * factor;
      fSeries->data->data[i].im = fSeries->data->data[i].im * factor;
    }
 /*
  * LALUnits manipulation
  */
   exponent.numerator=-1;
  exponent.denominatorMinusOne=1;/*2*/
  LALUnitRaise(status->statusPtr,
	       &tmpUnit1,
	       &(PSD->sampleUnits),
	       &exponent);
  CHECKSTATUSPTR (status);
  tmpUnitPair.unitOne=&tmpUnit1;
  tmpUnitPair.unitTwo=&(fSeries->sampleUnits);
  LALUnitMultiply(status->statusPtr,
		  &(fSeries->sampleUnits),
		  &tmpUnitPair);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/* End Fourier domain whitening */

/* Begin calibration routine */
void
LALTrackSearchCalibrateREAL4TimeSeries(LALStatus               *status,
				       REAL4TimeSeries         *signalvec,
				       COMPLEX8FrequencySeries *response)
{
  UINT4                      i=0;
  RealFFTPlan               *forwardPlan=NULL;
  RealFFTPlan               *reversePlan=NULL;
  COMPLEX8FrequencySeries   *signalFFT=NULL;
  UINT4                      planLength=0;

  INITSTATUS(status,"LALTrackSearchCalibrateREAL4TimeSeries",TSDATAC);
  ATTATCHSTATUSPTR (status);
  /* Need consistency checks for inputs so that the df of each match */
  /*
   * Setup FFT plans for FFTing the data segment
   */
  planLength=signalvec->data->length;
  LALCreateForwardREAL4FFTPlan(status->statusPtr,
			       &forwardPlan,
			       planLength,
			       0);
  CHECKSTATUSPTR (status);
  LALCreateReverseREAL4FFTPlan(status->statusPtr,
			       &reversePlan,
			       planLength,
			       0);
  CHECKSTATUSPTR (status);
  /*
   * Allocate RAM for temp Freq series
   */
  signalFFT = XLALCreateCOMPLEX8FrequencySeries("tmpSignalFFT",
				   &signalvec->epoch,
				   0,
				   1/signalvec->deltaT,
				   &lalDimensionlessUnit,
				   planLength/2+1);
  /*
   * FFT the data segment for calibration
   */
  LALForwardREAL4FFT(status->statusPtr,
		     signalFFT->data,
		     signalvec->data,
		     forwardPlan);
  CHECKSTATUSPTR (status);
  /*
   * Perform the frequency basis calibration as defined in
   * LSD Conventions Eq 23.1 p 601
   */
  for (i=0;i<signalvec->data->length;i++)
    {
      signalFFT->data->data[i].re=
	response->data->data[i].re*signalFFT->data->data[i].re;
      signalFFT->data->data[i].im=
	response->data->data[i].im*signalFFT->data->data[i].im;
    }
  /*
   * Bring this back to the time domain
   * this is the calibrated data set
   */
  LALReverseREAL4FFT(status->statusPtr,
		     signalvec->data,
		     signalFFT->data,
		     reversePlan);
  CHECKSTATUSPTR (status);
  /*
   * The 1/n factor need to be applied
   * See lsd-5 p 259 10.1
   */
  for (i=0;i<signalvec->data->length;i++)
    signalvec->data->data[i]= signalvec->data->data[i]/signalvec->data->length;

  /*
   * Destroy signalFFT Temp variable
   */
  XLALDestroyCOMPLEX8FrequencySeries(signalFFT);
  /*
   * Destroy the FFT plans
   */
  if (forwardPlan)
    {
      LALDestroyREAL4FFTPlan(status->statusPtr,&forwardPlan);
      CHECKSTATUSPTR (status);
    }
  if (reversePlan)
    {
      LALDestroyREAL4FFTPlan(status->statusPtr,&reversePlan);
      CHECKSTATUSPTR (status);
    }
  DETATCHSTATUSPTR(status);
  RETURN(status);

}
/* End calibration routine */


/* Begin Fourier Domain calibration routine */
void
LALTrackSearchCalibrateCOMPLEX8FrequencySeries(
					       LALStatus                 *status,
					       COMPLEX8FrequencySeries   *fSeries,
					       COMPLEX8FrequencySeries   *response
					       )
{
  UINT4          i=0;
  LALUnitPair    tmpUnitPair;
  REAL4          a=0;
  REAL4          b=0;
  REAL4          c=0;
  REAL4          d=0;

  INITSTATUS(status,"LALTrackSearchCalibrateCOMPLEX8FrequencySeries",TSDATAC);
  ATTATCHSTATUSPTR (status);
  /*
   * Error checking
   */
  ASSERT(fSeries != NULL,status,TSDATA_ENULL, TSDATA_MSGENULL);
  ASSERT(response != NULL,status,TSDATA_ENULL, TSDATA_MSGENULL);
  /*
   * Calibration is done via applying expression
   * 23.1 Conventions
   * s(f) = R(f;t)v(f)
   * Unit field is adjust appropriately and we return a
   * calibrated data solution
   */
  for(i=0;i<fSeries->data->length;i++)
    {
      a=fSeries->data->data[i].re;
      b=fSeries->data->data[i].im;
      c=response->data->data[i].re;
      d=response->data->data[i].im;
      /*(a+bi)*(c+di)*/
      fSeries->data->data[i].re=(a*c - b*d);
      fSeries->data->data[i].im=(a*d + b*c);
    }
  /*
   * Unit manipulation
   */
  tmpUnitPair.unitOne=&(fSeries->sampleUnits);
  tmpUnitPair.unitTwo=&(response->sampleUnits);
  LALUnitMultiply(status->statusPtr,
		  &(fSeries->sampleUnits),
		  &tmpUnitPair);
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/* End Fourier Domain calibration routine */

/*
 * This is the function to break up long stretch of input data
 * into the requested chunks accounting for overlap it should also
 * account for implicit overlaps via the segmenting buffer to remove
 * TFR edge effects.
 */
void
LALTrackSearchDataSegmenter(
			    LALStatus           *status,
			    REAL4TimeSeries     *TSSearchData,
			    TSSegmentVector     *PreparedData,
			    TSSearchParams       params)
{
  UINT4         k=0;
  UINT4         l=0;
  UINT4         j=0;
  REAL8         kTime;
  LIGOTimeGPS   timeInterval;

  INITSTATUS (status, "LALTrackSearchDataSegmenter", TSDATAC);

  /*
   * Error checking
   */
  ASSERT(PreparedData != NULL,status,TSDATA_ENULL,TSDATA_MSGENULL);
  ASSERT(TSSearchData != NULL,status,TSDATA_ENULL,TSDATA_MSGENULL);
  ASSERT((params.SegLengthPoints+(2*params.SegBufferPoints) == PreparedData->dataSeg[0]->data->length),
	 status,
	 TSDATA_ENUMZ,
	 TSDATA_MSGENUMZ);
  ASSERT(PreparedData->length == params.NumSeg,
	 status,
	 TSDATA_ESEGZ,
	 TSDATA_MSGESEGZ);
  ASSERT((params.SegBufferPoints == PreparedData->SegBufferPoints),
	 status,
	 TSDATA_ENUMZ,
	 TSDATA_MSGENUMZ);
  /*
   * We want to fill up our TSDataVector structure accounding for
   * desired number of segments and overlaps, we also need to account
   * implicitly specified buffer around each segment.
   */
  ATTATCHSTATUSPTR (status);
  for (l=0;l<PreparedData->length;l++)
    {
      /*Determine Segment Epoch*/
      kTime=TSSearchData->deltaT*k;
      XLALGPSSetREAL8(&(timeInterval), kTime);
      for (j=0;j<PreparedData->dataSeg[l]->data->length;j++)
	{
	  PreparedData->dataSeg[l]->data->data[j]=TSSearchData->data->data[k];
	  k++;
	};
      /*Ajust for segment overlap and buffer point implicit overlap*/
      k = k - (params.overlapFlag + 2*PreparedData->SegBufferPoints);
      /*Changing element count may cause memory leak*/
      /*PreparedData->dataSeg[l]->data->length = params.SegLengthPoints;*/
      PreparedData->dataSeg[l]->deltaT=TSSearchData->deltaT;
      PreparedData->dataSeg[l]->sampleUnits=TSSearchData->sampleUnits;
      PreparedData->dataSeg[l]->epoch.gpsSeconds=
	TSSearchData->epoch.gpsSeconds+timeInterval.gpsSeconds;
      PreparedData->dataSeg[l]->epoch.gpsNanoSeconds=
	TSSearchData->epoch.gpsNanoSeconds+timeInterval.gpsNanoSeconds;
      sprintf(PreparedData->dataSeg[l]->name,"%s","Initialized");
    };

  /*
   * End segment setup
   */
 DETATCHSTATUSPTR (status);
 RETURN (status);
}
/* End the data segmenter */

void
LALSVectorPolynomialInterpolation(
				  LALStatus         *status,
				  REAL4Sequence     *newDomain,
				  REAL4Sequence     *newRange,
				  REAL4Sequence     *domain,
				  REAL4Sequence     *range
				  )
{
  SInterpolatePar               interpolateParams;
  SInterpolateOut               interpolateResult;
  REAL4Vector                   *xfrag=NULL;
  REAL4Vector                   *yfrag=NULL;
  UINT4                         i=0;
  INT4                          j=0;
  INT4                          k=0;
  INT4                         bottomElement,topElement,currentElement;
  BOOLEAN                       cut;

  INITSTATUS(status,"LALSVectorPolynomialInterpolation",TSDATAC);
  ATTATCHSTATUSPTR (status);
  ASSERT(range->length > 5,status,TSDATA_EINTP,TSDATA_MSGEINTP);
  LALCreateVector(status->statusPtr,&xfrag,5);
  CHECKSTATUSPTR(status);
  LALCreateVector(status->statusPtr,&yfrag,5);
  CHECKSTATUSPTR(status);
  for (i=0;i<newRange->length;i++)
    {
      /*
       *Setup params
       */

      /*
       * Domain is ordered so use bisection technique to extract domain
       * point of interest.
       * Since we are interpolating using 5 points from the orignal range
       */
      currentElement=((UINT4) (domain->length/2));
      topElement=domain->length;
      bottomElement=0;
      cut=1;
      while (cut)
	{
	  if (newDomain->data[i] >= domain->data[currentElement])
	    {
	      bottomElement=currentElement;
	      currentElement=((UINT4) (topElement-bottomElement)/2)+bottomElement;
	    }
	  else
	    {
	      topElement=currentElement;
	      currentElement=topElement-((UINT4) (topElement-bottomElement)/2);
	    }
	  if ((topElement-bottomElement) == 1) cut=0;
	}
      if (currentElement < 2) currentElement=2;
      if ((domain->length - currentElement) < 2)
	currentElement=domain->length - 2;
      for (j=0,k=currentElement-2;k<currentElement+3;k++,j++)
	{
	  xfrag->data[j]=domain->data[k];
	  yfrag->data[j]=range->data[k];
	}
      interpolateParams.n=5;
      interpolateParams.x=xfrag->data;
      interpolateParams.y=yfrag->data;
      LALSPolynomialInterpolation(status->statusPtr,
				  &interpolateResult,
				  newDomain->data[i],
				  &interpolateParams);
      CHECKSTATUSPTR(status);
      newRange->data[i]=interpolateResult.y;
    }
  if (xfrag)
    {
      LALDestroyVector(status->statusPtr,&xfrag);
      CHECKSTATUSPTR(status);
    }
  if (yfrag)
    {
      LALDestroyVector(status->statusPtr,&yfrag);
      CHECKSTATUSPTR(status);
    }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}/* End interpolate vector like Matlab interp1 */



/*
 * As written this routine most likely wont work
 * I think there is variable scope problems
 */
void
connect2Segments(
		 TimeFreqRep    map,
		 Curve          *curveA,
		 Curve          *curveB
		 )
{
  INT4         i,j,k,m;
  INT4         deltaY;
  INT4         deltaX;
  INT4         newLength;
  INT4         intermediateSegment;
  Curve        tempCurve;
  /* Determine DeltaY */
  /* Determine DeltaX */
  /* Connect via greater of two */

  /* Determine points between head-tail to include */
  /* Calculate straight line segement between two points */

  deltaY=curveA->row[curveA->n]-curveB->row[0];
  deltaX=curveA->row[curveA->n]-curveB->row[0];
  newLength=0;
  intermediateSegment=(int)floor(sqrt(deltaY*deltaY+deltaX*deltaX));

  newLength=curveA->n+curveB->n+intermediateSegment;
  /* Allocate concatenated curve */
  tempCurve.row=LALMalloc(sizeof(INT4) * newLength);
  tempCurve.col=LALMalloc(sizeof(INT4) * newLength);
  tempCurve.fBinHz=LALMalloc(sizeof(REAL4) * newLength);
  tempCurve.gpsStamp=LALMalloc(sizeof(LIGOTimeGPS) * newLength);
  tempCurve.depth=LALMalloc(sizeof(REAL4) * newLength);

  /* Copy in segments information use map also */
  i=0;
  /* Copy in curveA */
  for (j=0;j<curveA->n;j++)
    {
      tempCurve.row[j]=curveA->row[j];
      tempCurve.col[j]=curveA->col[j];
      tempCurve.fBinHz[j]=curveA->fBinHz[j];
      tempCurve.gpsStamp[j].gpsSeconds=curveA->gpsStamp[j].gpsSeconds;
      tempCurve.gpsStamp[j].gpsNanoSeconds=curveA->gpsStamp[j].gpsNanoSeconds;
      tempCurve.depth[j]=curveA->depth[j];
      i++;
    };
  /* Copy intermediate curve information */
  /* X direction loop */
  m=0;
  for (j=curveA->n;j<(curveA->n+deltaX);j++)
    {
      for (k=0;k<deltaY;k++)
	{
	  tempCurve.row[curveA->n+m]=curveA->row[curveA->n]+j*k;
	  tempCurve.col[curveA->n+m]=curveA->col[curveA->n]+k;
	  tempCurve.fBinHz[curveA->n+m]=0;
	  tempCurve.gpsStamp[curveA->n+m].gpsSeconds=0;
	  tempCurve.gpsStamp[curveA->n+m].gpsNanoSeconds=0;
	  tempCurve.depth[curveA->n+m]=map.map[tempCurve.col[curveA->n+m]][tempCurve.row[curveA->n+m]];
	  i++;
	  if ( m >= intermediateSegment)
	    {
	      k=deltaY+1;
	      j=deltaX+1;
	    }
	}
    }
  /* Copy End segment information */
  for (k=curveA->n+deltaX;k<curveA->n+curveB->n;k++)
    {
      tempCurve.row[k]=curveB->row[k];
      tempCurve.col[k]=curveB->col[k];
      tempCurve.fBinHz[k]=curveB->fBinHz[k];
      tempCurve.gpsStamp[k].gpsSeconds=curveB->gpsStamp[k].gpsSeconds;
      tempCurve.gpsStamp[k].gpsNanoSeconds=curveB->gpsStamp[k].gpsNanoSeconds;
      tempCurve.depth[k]=curveB->depth[k];
      i++;
      k++;
    };
  /*Index i tracking elements of tempCurve*/
  tempCurve.n=i;
  /*Set total sum of power field */
  for (k=0;i<tempCurve.n;k++)
    {
      tempCurve.totalPower=tempCurve.totalPower+tempCurve.depth[k];
    };
  /* Set Used Marker on second curve when returned to calling function
   */
  curveB->trash='D';
  /* Set 'K'-keep Marker on tempCurve so we know is a joined track */
  tempCurve.trash='K';
  /* Set curveA to tempCurve free original curveA struct first!*/
  /* Freeing curveA */
  LALFree(curveA->row);
  LALFree(curveA->col);
  LALFree(curveA->gpsStamp);
  LALFree(curveA->fBinHz);
  LALFree(curveA->depth);
  LALFree(curveA);
  /*Set pointer to curveA to show joined product */
  curveA=&tempCurve;
  /*Deallocate curveB memory also*/
  LALFree(curveB->row);
  LALFree(curveB->col);
  LALFree(curveB->gpsStamp);
  LALFree(curveB->fBinHz);
  LALFree(curveB->depth);
  curveB->n=0;
  curveB->totalPower=0;
  return;
}

/*
 * Meant as companion function to connect track function
 */
void cleanLinkedList(
		     TrackSearchOut      *inList,
		     TrackSearchOut      *outList
		     )
{/*Will be used to clean up linked list structure and organize it */
  INT4 UsefulCurves;
  INT4 i;
  INT4 j;
  INT4 noc2;

  UsefulCurves=0; /*Array numbering starts with zero */
  outList->curves = LALMalloc(sizeof(outList->curves));
  /* Count up number of useful curves to keep */
  for(i=0;i<inList->numberOfCurves;i++)
    {
      if ((inList->curves[i].trash = 'K'))
	{
	  UsefulCurves++;
	  /* Expand the structure to take in another curve */
	  outList->curves = LALRealloc(outList->curves,
				       (sizeof(Curve) * (UsefulCurves)));
	  outList->numberOfCurves = UsefulCurves;
	  /* Numbering starts at zero */
	  noc2 = outList->numberOfCurves-1;
	  outList->curves[noc2].row=LALMalloc(sizeof(INT4) *
					      inList->curves[i].n);
	  outList->curves[noc2].col=LALMalloc(sizeof(INT4) *
					      inList->curves[i].n);
	  outList->curves[noc2].depth=LALMalloc(sizeof(REAL4) *
						inList->curves[i].n);
	  outList->curves[noc2].n = inList->curves[i].n;
	  outList->curves[noc2].totalPower = inList->curves[i].totalPower;
	  /* Copy the data over */
	  for(j=0;j<outList->curves[noc2].n;j++)
	    {
	      outList->curves[noc2].row[j]=inList->curves[i].row[j];
	      outList->curves[noc2].col[j]=inList->curves[i].col[j];
	      outList->curves[noc2].depth[j]=inList->curves[i].depth[j];
	    }
	  outList->curves[noc2].trash='K';
	  /* Move on to check next curve */
	}
    }
  return;
}
