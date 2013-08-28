/*
*  Copyright (C) 2007 Jolien Creighton
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

/* ****** INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

/* ****** INCLUDE ANY LDAS LIBRARY HEADERS ************/

/* ****** INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALCorrelation.h>
#include <math.h>

/* ****** DEFINE LOCAL CONSTANTS AND MACROS ************/

/* ****** DECLARE LOCAL (static) FUNCTIONS ************/
/* (definitions can go here or at the end of the file) */

static REAL4 findMean( REAL4 data[], UINT4 length );

static REAL4 findVariance( REAL4 data[], REAL4 mean, UINT4 length );

static REAL4 findCrossProduct( REAL4 data1[], REAL4 mean1, REAL4 data2[], REAL4 mean2, UINT4 length );


/* ****** DEFINE GLOBAL FUNCTIONS ************/

/**
 * \brief LALCorrelation() is designed to compute a time shifted correlation between two time series given in <tt>input-\>one</tt> and <tt>input-\>two</tt>.
 * \author Yakushin, Igor
 *
 * The maximum time shift in nanoseconds is given in <tt>params->maxTimeShiftNan</tt>. The output consists of a correlation for each
 * time shift in the range <tt>out->timeShiftedCorrelation</tt>, maximum and minimum values of correlations and corresponding time shifts.
 * The original intention is to use this function to test coincendence bursts found in two detectors for correlation.
 * For this to work one must apply a response function to the raw time series in order to get rid of hardware specific contributions
 * to each time series. The signature of the coincendence event is a clear maximum above some threshold in the graph of correlation vs time shift
 * (no more than 10 ms).
 *
 * One might, of course, try to use the code to search for any correlations in the data caused by any kind of gravitational waves but
 * that seems to be too computationally expensive.
 *
 * ### Algorithm ###
 *
 * Just a straightforward computation of correlation for different time shifts. This computation is applied to time series of the length
 * <tt>originalLength - maxShift</tt>.
 *
 * ### Notes ###
 *
 * One must figure out how to prefilter the raw data, what length of time series s appropriate to use, what
 * threshold on the maximum correlation value should be applied to declare a good correlation.
 */
void
LALCorrelation( LALStatus                      *status,
		OutputCorrelation              **out,
		const InputCorrelation         *input,
		const CorrelationParams        *params)


{
  /* ****** DECLARE VARIABLES; for example: ************/

  REAL4 mean1, mean2, var1, var2, cov12, cor;
  INT4 shift;
  INT4 i;
  UINT4 length;
  REAL4 *data1, *data2;
  OutputCorrelation *output;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* ****** CHECK VALIDITY OF ARGUMENTS; for example ************/


  if ( input == NULL || params == NULL ) {
    ABORT( status, LALCORRELATIONH_ENULLP, LALCORRELATIONH_MSGENULLP );
  }

  if( input->one == NULL || input->two == NULL ){
    ABORT( status, LALCORRELATIONH_ENULLP, LALCORRELATIONH_MSGENULLP );
  }

  if (input->one->epoch.gpsSeconds != input->two->epoch.gpsSeconds ||
      input->one->epoch.gpsNanoSeconds != input->two->epoch.gpsNanoSeconds) {
    ABORT( status, LALCORRELATIONH_ESTART, LALCORRELATIONH_MSGESTART);
  }

  if(input->one->deltaT != input->two->deltaT ||
     input->one->data->length!=input->two->data->length){
    ABORT( status, LALCORRELATIONH_ESAMPLING, LALCORRELATIONH_MSGESAMPLING);
  }


  output=(OutputCorrelation*)LALMalloc(sizeof(OutputCorrelation));

  /* ****** EXTRACT INPUTS AND PARAMETERS ************/

/*    printf("LALCorrelation: maxTimeShiftNan=%g, deltaT=%g\n",params->maxTimeShiftNan, input->one->deltaT); */

  shift=(INT4)(params->maxTimeShiftNan/pow(10.0,9.0)/input->one->deltaT);

  output->timeShiftedCorrelation=LALMalloc(sizeof(REAL4)*(2*shift+1));

  length = input->one->data->length - shift;

  output->maxCorrelationValue=-2.0;
  output->maxCorrelationTimeShift=length+2;
  output->minCorrelationValue=2.0;
  output->minCorrelationTimeShift=length+2;

  output->shift=shift;

  /* ****** DO ANALYSIS ************/

  for(i=-shift;i<=shift;i++){
    data1=input->one->data->data;
    data2=input->two->data->data;

    if(i<0){
      data2+=-i;
    }
    else if(i>0){
      data1+=i;
    }


    mean1=findMean(data1,length);
    mean2=findMean(data2,length);

    var1=findVariance(data1,mean1,length);
    var2=findVariance(data2,mean2,length);

    cov12=findCrossProduct(data1,mean1,data2,mean2,length);
    cor=cov12/var1/var2;

    output->timeShiftedCorrelation[shift+i]=cor;

    if(cor>output->maxCorrelationValue){
      output->maxCorrelationValue=cor;
      output->maxCorrelationTimeShift=i;
    }
    if(cor<output->minCorrelationValue){
      output->minCorrelationValue=cor;
      output->minCorrelationTimeShift=i;
    }
  }





  /* ****** CONSTRUCT OUTPUT ************/

  output->start = input->one->epoch;
  output->length = input->one->data->length;
  output->deltaT = input->one->deltaT;

  *out=output;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

static REAL4 findMean( REAL4 data[], UINT4 length ) {
  UINT4     i;
  REAL4     mean;

  mean = 0.0;

  for ( i = 0 ; i < length ; ++i )
  {
    mean += data[i];
  }

  return mean / length;
}

static REAL4 findVariance( REAL4 data[], REAL4 mean, UINT4 length ) {
  UINT4     i;
  REAL4     var;

  var = 0.0;

  for ( i = 0 ; i < length ; ++i )
  {
    var += (data[i]-mean)*(data[i]-mean) ;
  }

  return sqrt( (double)(var / (length - 1)) );
}

static REAL4 findCrossProduct( REAL4 data1[], REAL4 mean1, REAL4 data2[], REAL4 mean2, UINT4 length ){
  UINT4 i;
  REAL4 result=0.0;

  for(i=0;i<length;i++){
    result+=(data1[i]-mean1)*(data2[i]-mean2);
  }

  result/=length-1;

  return result;

}

