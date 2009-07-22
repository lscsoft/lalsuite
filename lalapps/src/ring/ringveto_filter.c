/*
*  Copyright (C) 2007 Chad Hanna
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

#include <math.h>
#include <string.h>
#include <limits.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Ring.h>
#include <lal/Date.h>

#include "lalapps.h"
#include "errutil.h"
#include "gpstime.h"
#include "ringveto.h"

RCSID( "$Id$" );

static REAL8 compute_template_variance(
    COMPLEX8FrequencySeries  *stilde,
    REAL4FrequencySeries     *invspec,
    REAL8                     dynRangeFac
    );

static int filter_segment_template(
    REAL4TimeSeries          *result,
    COMPLEX8FrequencySeries  *rtilde,
    COMPLEX8FrequencySeries  *stilde,
    COMPLEX8FrequencySeries  *segment,
    REAL4FFTPlan             *plan
    );

static SnglRingdownTable * find_events(
    SnglRingdownTable  *events,
    UINT4              *numEvents,
    REAL4TimeSeries    *result,
    REAL4               tmpltSigma,
    REAL4               tmpltFrequency,
    REAL4               tmpltQuality,
    RingVetoResults    *Result,
    RingVetoResults    *firstResult,
    RingVetoCC         *thisCC,
    UINT4              sC,
    struct ring_params *params
    );


static REAL4Vector *snipInvSpectrum(UINT4 snipLength,
                                    REAL4FrequencySeries *invspec);

static RingVetoCC computeCC( RingVetoResults *Result, 
                             REAL4FrequencySeries *snipSpec,
                             COMPLEX8FFTPlan *fwdplan,
                             COMPLEX8FFTPlan *revplan );


static REAL4 computeChisq(RingVetoResults *thisResult,
                                     RingVetoResults *firstResult,
                                     RingVetoCC *RVCC,
                                     struct ring_params *params,
                                     UINT4 i,
                                     UINT4 j );



struct tagRV_BT{
  REAL4 beta;
  INT4 tau;
};
typedef struct tagRV_BT RV_BT;

static RV_BT getBT(COMPLEX8FrequencySeries *A,
                   COMPLEX8FrequencySeries *B,
                   COMPLEX8FFTPlan *revplan);

static void normSnip( COMPLEX8Vector *vec );

SnglRingdownTable * ringveto_filter(
    RingDataSegments         *segments,
    RingVetoTemplateBank     *vetobank,
    REAL4FrequencySeries     *invSpectrum,
    REAL4FFTPlan             *fwdPlan,
    REAL4FFTPlan             *revPlan,
    struct ring_params       *params
    )
{
  SnglRingdownTable       *events = NULL; /* head of linked list of events */
  RingTemplateBank        *bank = NULL; /* the smaller banks */
  RingVetoTemplateBank    *firstvetobank = NULL;
  REAL4TimeSeries          signal;
  /*REAL4TimeSeries          result;*/
  RingVetoResults          *Result = NULL;
  RingVetoResults          *firstResult = NULL;
  RingVetoResults          *thisResult = NULL;
  RingVetoCC		   thisCC;
  COMPLEX8FFTPlan  	   *fwdSnipPlan = NULL;
  COMPLEX8FFTPlan             *revSnipPlan = NULL;
  COMPLEX8FrequencySeries  stilde;
  COMPLEX8FrequencySeries  rtilde;
  REAL4FrequencySeries     snipSpec; /*this will actually store the amplitude */
                                     /*spectrum NOT PSD! */
  REAL4TimeSeries          chisqSeries;
  UINT4 segmentLength;
  UINT4 sgmnt;
  UINT4 tmplt;
  UINT4 subCounter = 0;
  UINT4 sC = 0;
  const REAL4 snipFac = 2.0; /*template snip length in seconds should be a 
                               power of two */
  SnglRingdownTable *thisTmplt = NULL;
  REAL8 sigma = 0;
  UINT4 snipLength;
  if ( ! params->doFilter )
    return NULL;
  
  segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
  snipLength = floor( snipFac * params->sampleRate + 0.5 );
  snipSpec.deltaF = 1.0/snipLength;
  snipSpec.data = snipInvSpectrum( snipLength, invSpectrum );
  snipSpec.sampleUnits = invSpectrum->sampleUnits;
  fwdSnipPlan = XLALCreateCOMPLEX8FFTPlan( snipLength, 1, 0 );
  revSnipPlan = XLALCreateCOMPLEX8FFTPlan( snipLength, 0, 0 );
  memset( &signal, 0, sizeof( signal ) );
  /* changed to the structure Result not just the timeseries result */
  /* memset( Result, 0, sizeof( Result ) );*/ /* not necessary */
  memset( &stilde, 0, sizeof( stilde ) );
  memset( &rtilde, 0, sizeof( rtilde ) );

  signal.deltaT = 1.0/params->sampleRate;
  signal.sampleUnits = lalStrainUnit;
  chisqSeries.deltaT = 1.0/params->sampleRate;
  rtilde.deltaF = 1.0 / params->segmentDuration;
  signal.data = XLALCreateREAL4Vector( segmentLength );
  stilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  rtilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  chisqSeries.data = XLALCreateREAL4Vector( segmentLength );

  firstvetobank = vetobank; /* just remember the first one */
  /* loop over all elements in the template bank */
  UINT4 numEvents = 0;
  while (vetobank->next && vetobank->bank){
    subCounter++;
    bank = vetobank->bank;
    for ( sgmnt = 0; sgmnt < segments->numSgmnt; ++sgmnt )
    {
      /* set up the Results structure to store all the ~10 templates */
      if (1){/*(sgmnt == 0){allocate memory for the first segment (TRUE) */
        Result = (RingVetoResults *) calloc(1, sizeof(RingVetoResults));
        firstResult = Result; /* remember the first one */
        }
      else Result = firstResult;
      
      for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt ){
        /* changed to Result from result - extra layer of pointers... */
        if (1){
          Result->result.data = XLALCreateREAL4Vector( segmentLength );
          Result->tmpSnip.data = XLALCreateCOMPLEX8Vector( snipLength );
          Result->tmpSnipTilde.data = XLALCreateCOMPLEX8Vector( snipLength );
          Result->tmpSnip.deltaT = 1.0/params->sampleRate;
          Result->tmpSnipTilde.deltaF = snipSpec.deltaF;
          Result->tmpSnip.sampleUnits = lalStrainUnit;
          Result->numResults = bank->numTmplt;
          thisTmplt = bank->tmplt + tmplt;
          verbose( "creating template %d for segment %d in subbank %d\n",
                  tmplt, sgmnt, subCounter );
          }
        /* make template and fft it */
        /* a two second snip is sufficient for the cross correlations*/
        XLALComputeRingTemplate( &signal, thisTmplt );
        if (params->vetoThresh){
          /* snip the template for later WILL THIS WORK?*/
          for(sC=0;sC<snipLength;sC++){ 
             Result->tmpSnip.data->data[sC].re = signal.data->data[sC];
             Result->tmpSnip.data->data[sC].im = 0;
             }
          normSnip( Result->tmpSnip.data ); /* normalize the template snip */
          }          
        snprintf( signal.name, sizeof(signal.name), "BANK_%u_TMPLT_%u", 
                  subCounter, tmplt );
        /*write_REAL4TimeSeries( &signal );*/
        XLALREAL4TimeFreqFFT( &stilde, &signal, fwdPlan );
        snprintf( stilde.name, sizeof(stilde.name), "BANK_%u_TMPLT_%u_FFT", 
                  subCounter, tmplt );
        /*    write_COMPLEX8FrequencySeries( &stilde );    */

        /* compute sigma for this template MAKE SURE AND SAVE THIS*/
        /* to memory so that all the subbanks variances are available */
        sigma = sqrt( compute_template_variance( &stilde, invSpectrum,
            params->dynRangeFac ) );
        Result->sigma = sigma; /* store this value for all these tmplts */
        verbose( "  filtering segment %d against template %d in subbank %d\n", 
                sgmnt, tmplt, subCounter );

        /* filter the segment with the template */
        filter_segment_template( &(Result->result), &rtilde, &stilde,
            &segments->sgmnt[sgmnt], revPlan );

        /* write filter output if requested */
        if ( params->writeFilterOutput )
        { /* guess we better normalize it so it is SNR-like... */
          REAL4 snrFactor = 2 * params->dynRangeFac / sigma;
          UINT4 k;
          snprintf( Result->result.name, sizeof(Result->result.name), 
                      "SNR_BANK_%u_TMPLT_%u_SEGMENT_%u",
                  subCounter, tmplt, sgmnt );
          for ( k = 0; k < Result->result.data->length; ++k )
            Result->result.data->data[k] *= snrFactor;
          write_REAL4TimeSeries( &(Result->result) );
        }
        /* save some additional template parameters and add a link */
        if (1){ /*only needed for first segment SET TO TRUE IF IT BREAKS*/
          Result->quality = thisTmplt->quality;
          Result->frequency = thisTmplt->frequency;
          Result = Result->next = 
                 (RingVetoResults *) calloc(1, sizeof(RingVetoResults));
          }
        else Result = Result->next; /*just write over results*/
                
      } /* end loop over templates */
      Result->next = NULL; /* terminate the linked list */
      Result = firstResult;
      /* Add a for loop here that computes the cross correlation */
      /* remember to make a time series which is shorter than the */
      /* 256 second signal - that is way too much time to pad and */
      /* fft */
      if (sgmnt==0 && params->vetoThresh) thisCC = computeCC( Result, &snipSpec, fwdSnipPlan, revSnipPlan );
      /* search through results for threshold crossings and record events */
      /* REWRITE this to loop over ALL RESULTS!!! */
      numEvents = 0;
      sC = 0;
      while(Result->next){
        /*if (params->vetoThresh) computeChisqVec(chisqSeries.data,Result,firstResult,&thisCC,params,sC);*/
        events = find_events( events, &numEvents, &Result->result, 
            Result->sigma, Result->frequency, Result->quality,
            Result, firstResult, &thisCC, sC, params );
        if ( params->writeFilterOutput ){
          snprintf( chisqSeries.name, sizeof(chisqSeries.name),
                      "CHISQ_BANK_%u_TMPLT_%u_SEGMENT_%u",
                       subCounter, sC, sgmnt );
          write_REAL4TimeSeries( &(chisqSeries) );
          }
        params->numEvents += numEvents;
        verbose( "found %u event%s in subbank %u and segment%d\n", numEvents,
          numEvents == 1 ? "" : "s", subCounter, sgmnt );
        Result = Result->next;
        numEvents = 0;
        sC++;
        }
      Result=firstResult; /* go to head */
      /* clean up all the results if done with sub bank*/
      if (sgmnt == (segments->numSgmnt-1) && params->vetoThresh){
        XLALDestroyREAL4Vector(thisCC.Beta);
        XLALDestroyINT4Vector(thisCC.Tau);
        }
 
      if (1){
        Result = firstResult;
        while (Result->next){
          thisResult = Result;
          Result = Result->next;
          XLALDestroyREAL4Vector(thisResult->result.data);
          XLALDestroyCOMPLEX8Vector(thisResult->tmpSnip.data);
          XLALDestroyCOMPLEX8Vector(thisResult->tmpSnipTilde.data);
          free(thisResult);
          }
        }
    } /* end loop over segments */
  vetobank = vetobank->next;
  } /* end loop over sub banks */

  /* tell the world what we saw */
  verbose( "found %u event%s\n", params->numEvents,
      params->numEvents == 1 ? "" : "s" );

  /* cleanup */
  /*XLALDestroyREAL4Vector( result.data );*/
  XLALDestroyCOMPLEX8Vector( rtilde.data );
  XLALDestroyCOMPLEX8Vector( stilde.data );
  XLALDestroyREAL4Vector( signal.data );
  vetobank = firstvetobank; /* make sure you send back the head */
  return events;
}

static REAL8 compute_template_variance(
    COMPLEX8FrequencySeries  *stilde,
    REAL4FrequencySeries     *invspec,
    REAL8                     dynRangeFac
    )
{
  UINT4 k;
  REAL8 var;

  var = 0;
  for ( k = 0; k < stilde->data->length; ++k )
  {
    REAL8 re = stilde->data->data[k].re;
    REAL8 im = stilde->data->data[k].im;
    var += invspec->data->data[k] * (re*re + im*im);
  }

  var *= 4.0 * dynRangeFac * dynRangeFac * stilde->deltaF;
  return var;
}

static int filter_segment_template(
    REAL4TimeSeries          *result,
    COMPLEX8FrequencySeries  *rtilde,
    COMPLEX8FrequencySeries  *stilde,
    COMPLEX8FrequencySeries  *segment,
    REAL4FFTPlan             *plan
    )
{
  char *s;

  /* name of rtilde */
  snprintf( rtilde->name, sizeof( rtilde->name ), "%s_%s",
      segment->name, stilde->name );
  /* name of result is the same but without the _FFT */
  strncpy( result->name, rtilde->name, sizeof( result->name ) - 1 );
  /* make sure that the string ends with _FFT */
  s = result->name + strlen( result->name ) - strlen( "_FFT" );
  if ( strcmp( s, "_FFT" ) == 0 )
    *s = 0; /* it does: terminate here */

  /* multiply segment by filter and store in fft of result */
  XLALCCVectorMultiplyConjugate( rtilde->data, segment->data, stilde->data );
  XLALUnitMultiply( &rtilde->sampleUnits, &segment->sampleUnits, &stilde->sampleUnits );
  rtilde->epoch = segment->epoch;
  /*write_COMPLEX8FrequencySeries( rtilde );*/

  /* inverse fft to obtain result */
  XLALREAL4FreqTimeFFT( result, rtilde, plan );
  /*this works without the FFT so the FFT must be raising the units error */
  /*this works without the FFT so the FFT must be raising the units error */
 /*write_REAL4TimeSeries( result );*/
  result->epoch = rtilde->epoch;
  return 0;
}



/* NOTE: numEvents must be number of events _FROM_CURRENT_TEMPLATE_ so far. */
/* It must be set to zero when filtering against a new template is started. */
static SnglRingdownTable * find_events(
    SnglRingdownTable  *events,
    UINT4              *numEvents,
    REAL4TimeSeries    *result,
    REAL4               tmpltSigma,
    REAL4               tmpltFrequency,
    REAL4               tmpltQuality,
    RingVetoResults    *Result,
    RingVetoResults    *firstResult,
    RingVetoCC         *thisCC,
    UINT4              sC,
    struct ring_params *params
    )
{
  const REAL4 efolds = 10.0; /* number of efolds of ringdown in template */
  SnglRingdownTable *thisEvent = events; /* the current event */
  REAL4 snrFactor; /* factor to convert from filter result to snr */
  REAL4 threshold; /* modified threshold on filter result (rather than snr) */
  REAL4 vetoThresh;
  REAL4 filterDuration;
  INT8  lastTimeNS;
  INT8  gapTimeNS;
  UINT4 segmentStride;
  UINT4 eventCount = 0;
  UINT4 jmin;
  UINT4 jmax;
  UINT4 j;
/*  LALStatus             status = blank_status;*/
 /* LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };*/

  /* compute filter duration: sum of rindown duration and spec trunc duration */
  filterDuration  = efolds * tmpltQuality / (LAL_PI * tmpltFrequency);
  filterDuration += params->truncateDuration;
  
  /* gap time: filter duration unless explicitly specified */
  if ( params->maximizeEventDuration < 0 ) /* maximize based on duration*/
    gapTimeNS = sec_to_ns( filterDuration );
  else /* maximize with specified duration */
    gapTimeNS = sec_to_ns( params->maximizeEventDuration );

  /* determine if we are maximizing over current event or not */
  /* set lastTimeNS to be time of last event or before any possible */
  /* event in this segment depending on whether there was a previous */
  /* event from this template */
  if ( *numEvents ) /* continue maximizing over the most recent event */
    lastTimeNS = epoch_to_ns( &thisEvent->start_time );
  else /* set last time to a time before any possible event in this segment */
    lastTimeNS = epoch_to_ns( &result->epoch ) - gapTimeNS - (INT8)(1000000000);

  /* compute modified threshold on filter output rather than snr */
  snrFactor = 2 * params->dynRangeFac / tmpltSigma;
  threshold = params->threshold / snrFactor;
  vetoThresh = params->vetoThresh;
  /* compute start and stop index for scanning time series */
  segmentStride = floor( params->strideDuration / result->deltaT + 0.5 );
  jmin = segmentStride/2;
  jmax = jmin + segmentStride;

  for ( j = jmin; j < jmax; ++j )
    if ( (fabs( result->data->data[j] ) > threshold) && 
               (computeChisq(Result,firstResult,thisCC,params,sC,j) 
                < vetoThresh)) /* threshold crossing */
    {
      REAL4 snr;
      INT8  timeNS;
      REAL4 sigma;
      REAL4 amp;
      REAL4 chisq;

      snr     = fabs( result->data->data[j] ) * snrFactor;
      chisq   = computeChisq(Result,firstResult,thisCC,params,sC,j);
      timeNS  = epoch_to_ns( &result->epoch );
      timeNS += sec_to_ns( j * result->deltaT );

      if ( timeNS > lastTimeNS + gapTimeNS ) /* new distinct event */
      {
        /* create a new node on the linked list */
        thisEvent = LALCalloc( 1, sizeof( *thisEvent ) );
        thisEvent->next = events;
        events = thisEvent;

        /* copy general information about the filter */
        strncpy( thisEvent->ifo, params->ifoName, sizeof( thisEvent->ifo ) );
        strncpy( thisEvent->channel, strchr( params->channel, ':') + 1, 
            sizeof( thisEvent->channel ) );
/*       LAL_CALL( LALGPStoGMST1( &status, &(thisEvent->start_time_gmst), 
                        &(thisEvent->start_time), &gmstUnits ), &status); */ 
        /*this isnt working properly, will just leave awhile, as it is not
         * needed */

        
        thisEvent->frequency = tmpltFrequency;
        thisEvent->quality = tmpltQuality;
        thisEvent->mass = ( 1.0 - 0.63 * pow( (2.0 / thisEvent->quality) , (2.0 / 3.0) ) ) *
            pow( LAL_C_SI, 3.0) / LAL_G_SI / LAL_TWOPI / LAL_MSUN_SI / thisEvent->frequency; 
        thisEvent->spin = 1.0 - pow( (2.0 / thisEvent->quality) , (20.0 / 9.0) );
        
        /* specific information about this threshold crossing */
        ns_to_epoch( &thisEvent->start_time, timeNS );
        thisEvent->snr = snr;
        thisEvent->epsilon = chisq;

        
        amp = sqrt( 5.0 / 2.0 * 0.01 )  * ( LAL_G_SI * thisEvent->mass * LAL_MSUN_SI / 
            pow( LAL_C_SI, 2) * 2.0 / LAL_PC_SI /1000000.0 ) * pow( thisEvent->quality, -0.5 ) * 
            pow( 1.0 + 7.0 / 24.0 / pow( thisEvent->quality, 2.0), -0.5 ) *
            pow(  1.0 - 0.63 * pow( 1.0 - thisEvent->spin,0.3 ), -0.5);
        sigma=tmpltSigma * amp;
        thisEvent->sigma_sq = pow(sigma, 2.0);
        thisEvent->eff_dist = sigma / thisEvent->snr;
        thisEvent->amplitude=amp;       
        ++eventCount;
      }
      else if ( snr > thisEvent->snr ) /* maximize within a set of crossings */
      {
        /* update to specific information about this threshold crossing */
        ns_to_epoch( &thisEvent->start_time, timeNS );
        thisEvent->snr        = snr;
        thisEvent->epsilon    = chisq;

        amp = sqrt( 5.0 / 2.0 * 0.01 )  * ( LAL_G_SI * thisEvent->mass * LAL_MSUN_SI/
            pow( LAL_C_SI, 2) * 2.0 / LAL_PC_SI /1000000.0 ) * pow( thisEvent->quality, -0.5 ) *
          pow( 1.0 + 7.0 / 24.0 / pow( thisEvent->quality, 2.0), -0.5 ) *
          pow(  1.0 - 0.63 * pow( 1.0 - thisEvent->spin,0.3 ), -0.5);
        sigma=tmpltSigma * amp;
        thisEvent->eff_dist = sigma / thisEvent->snr;        
        thisEvent->amplitude = amp;
      }
      
      /* update last threshold crossing time */
      lastTimeNS = timeNS;
    }

  *numEvents += eventCount;
  return events;
}

static REAL4Vector *snipInvSpectrum(UINT4 snipLength, 
                                      REAL4FrequencySeries *invspec) {
  FILE *FP;
  UINT4 i=0;
  UINT4 cnt=0;
  REAL4FrequencySeries smallspec;
  REAL4FrequencySeries doublesmallspec;
  UINT4 invSpecLength = snipLength/2 + 1;
  smallspec.data = XLALCreateREAL4Vector( invSpecLength );
  doublesmallspec.data = XLALCreateREAL4Vector(snipLength);
  int factor = floor( invspec->data->length/invSpecLength+0.5 );
  /*this is really bad if both aren't a power of two */
  smallspec.data->data[0] = invspec->data->data[0];
  for(i=1;i<invspec->data->length;i++){
    if( !(i%factor) && cnt <invSpecLength){
      smallspec.data->data[cnt] = 
         (invspec->data->data[i-1] +  
         invspec->data->data[i] +  
         invspec->data->data[i+1])/3.0;
      cnt++;
      }
    }
  cnt = invSpecLength;
  for(i=0;i<smallspec.data->length;i++) 
     smallspec.data->data[i] = sqrt(smallspec.data->data[i]);
  for(i=0;i<doublesmallspec.data->length/2;i++){
     doublesmallspec.data->data[i] = smallspec.data->data[cnt];
     cnt--;
     }
  cnt = 0;
  for(i=doublesmallspec.data->length/2;i<doublesmallspec.data->length;i++){
     doublesmallspec.data->data[i] = smallspec.data->data[cnt];
     cnt++;
     }

  /*FP = fopen("miniasd.dat","w");
  for(i=0;i<doublesmallspec.data->length;i++) fprintf(FP,"%f\n",doublesmallspec.data->data[i]);*/

  return doublesmallspec.data;
  }

static RingVetoCC computeCC( RingVetoResults *Result,
                             REAL4FrequencySeries *snipSpec,
			     COMPLEX8FFTPlan *fwdplan,
			     COMPLEX8FFTPlan *revplan){

  RingVetoResults *first = NULL;
  RingVetoResults *thisResult = NULL;
  REAL4TimeSeries thisCCresult; 
  REAL4Vector *Beta = NULL;
  INT4Vector *Tau = NULL; 
  FILE *FP = NULL;
  RingVetoCC finalCC;
  UINT4 cnt = 0;
  UINT4 i=0;
  UINT4 j=0;
  RV_BT BT;
  first = Result;
  thisCCresult.data = XLALCreateREAL4Vector( Result->tmpSnip.data->length );
 
  while(Result->next){
    cnt++;
    Result=Result->next;
    }

  Beta  = XLALCreateREAL4Vector(cnt*cnt);
  Tau  = XLALCreateINT4Vector(cnt*cnt);
  Result = first;
  /*I will go ahead and FFT all the templates now */
   
  while(Result->next){
    XLALCOMPLEX8TimeFreqFFT( &Result->tmpSnipTilde, 
                          &Result->tmpSnip, fwdplan );
    /* whiten the templates with the amplitude spectral density */
    XLALSCVectorMultiply( Result->tmpSnipTilde.data, 
                          snipSpec->data, 
                          Result->tmpSnipTilde.data);
    Result->numResults = cnt; /*store the number of results */
    normSnip( Result->tmpSnipTilde.data ); /* normalize whitened temps */
    Result=Result->next;
    }
  Result = first;
 /*FP = fopen("template.dat","w");
 for(i=0;i<Result->tmpSnip.data->length;i++)
   fprintf(FP,"%f\n",Result->tmpSnip.data->data[i].re);  
  */
  i=0;
  /* Now I will evaluate the cc matrix */
  Result = first;
  while(Result->next){
    thisResult = Result;
    BT = getBT(&thisResult->tmpSnipTilde,&thisResult->tmpSnipTilde,revplan);
    Beta->data[cnt*i+i] = BT.beta;
    Tau->data[cnt*i+i] = BT.tau;
    verbose("BETA[%d,%d] = %e \n",i,i,BT.beta);
    i++;
    Result=Result->next;
    }
  Result = first;
  i=0;
  while(Result->next){
    thisResult = Result;
    j=i+1;
    while(thisResult->next && j<thisResult->numResults){ /*also check j */
      /* one pair */
      BT = getBT(&Result->tmpSnipTilde,&thisResult->next->tmpSnipTilde,
                 revplan);
      verbose("norm BETA[%d,%d] = %e \n",j,i,BT.beta/Beta->data[cnt*i+i]);
      Beta->data[cnt*i+j] = BT.beta; /*/sqrt(Beta->data[cnt*i+i]*
                                         Beta->data[cnt*j+j]);*/
      /*if(BT.tau > Result->tmpSnipTilde.data->length/2) 
         BT.tau-=Result->tmpSnipTilde.data->length;*/
      Tau->data[cnt*i+j] = BT.tau;
      /* conjugate pair */
      BT = getBT(&thisResult->next->tmpSnipTilde,&Result->tmpSnipTilde,
                 revplan);
      Beta->data[cnt*j+i] = BT.beta; /*/sqrt(Beta->data[cnt*i+i]*
                                         Beta->data[cnt*j+j]);*/
      /*if(BT.tau > Result->tmpSnipTilde.data->length/2)          
         BT.tau-=Result->tmpSnipTilde.data->length;*/
      Tau->data[cnt*j+i] = BT.tau;
      thisResult=thisResult->next;
      j++;
      }
    i++;
    Result=Result->next;
    }

  /*for(i=0;i<cnt;i++) Beta->data[cnt*i+i] = 1;  normalize */

  Result = first;
  XLALDestroyREAL4Vector( thisCCresult.data );
  finalCC.Beta = Beta;
  finalCC.Tau = Tau;
  return finalCC;
  }


/* this has to be switched to a COMPLEX FFT */
static RV_BT getBT(COMPLEX8FrequencySeries *A, 
                   COMPLEX8FrequencySeries *B, 
                   COMPLEX8FFTPlan *revplan){
  
  COMPLEX8FrequencySeries *C=NULL;
  COMPLEX8TimeSeries    *OUT=NULL;
  RV_BT BT;
  UINT4 i = 0;
  BT.beta = 0;
  BT.tau = 0;
  FILE *FP = NULL;
  OUT = calloc(1,sizeof(COMPLEX8FrequencySeries));
  /*C = calloc(1,sizeof(COMPLEX8FrequencySeries));*/
  OUT->data = XLALCreateCOMPLEX8Vector(A->data->length);
  /*C->data = XLALCreateCOMPLEX8Vector(A->data->length);*/
  OUT->sampleUnits = lalDimensionlessUnit; /* fixed error with FFT? */
  /*C->sampleUnits = lalDimensionlessUnit;
  C->deltaF = A->deltaF;*/
  XLALCCVectorMultiplyConjugate(OUT->data,B->data,A->data);
  /*XLALCOMPLEX8FreqTimeFFT( OUT, C, revplan );*/ 
  
  for(i=0; i < OUT->data->length; i++){
    BT.beta += A->data->data[i].re*B->data->data[i].re+
                 B->data->data[i].im*A->data->data[i].im;
    BT.tau = 0;
    /*if (sqrt(OUT->data->data[i].re*OUT->data->data[i].re 
             + OUT->data->data[i].im*OUT->data->data[i].im) > BT.beta){
      BT.beta = sqrt(OUT->data->data[i].re*OUT->data->data[i].re+ 
                 OUT->data->data[i].im*OUT->data->data[i].im);
      BT.tau = i;*/
    };
  /*FP = fopen("OUT.dat","w");
  for(i=0;i<OUT->data->length;i++) fprintf(FP,"%e\n",sqrt(OUT->data->data[i].re*OUT->data->data[i].re + OUT->data->data[i].im*OUT->data->data[i].im));*/
  XLALDestroyCOMPLEX8Vector(OUT->data);
  /*XLALDestroyCOMPLEX8Vector(C->data);*/
  free(OUT);
  /*free(C);*/
  return BT;
  }

static void normSnip( COMPLEX8Vector *vec){

  UINT4 i = 0;
  REAL4 vecSum = 0;
  for(i=0;i<vec->length;i++) vecSum+= 
     vec->data[i].re*vec->data[i].re + vec->data[i].im*vec->data[i].im;
  for(i=0;i<vec->length;i++) {
    vec->data[i].re/=sqrt(vecSum);
    vec->data[i].im/=sqrt(vecSum);
    }
  }


static REAL4 computeChisq (RingVetoResults *thisResult,
                                     RingVetoResults *firstResult,
                                     RingVetoCC *RVCC,
                                     struct ring_params *params,
                                     UINT4 i,
                                     UINT4 dataIX ) {

  FILE *FP = NULL;
  UINT4 j = 0;
  UINT4 k = 0;
  REAL4 deltasq = 0.1; /*THIS SHOULD BE A PARAMETER?! */
  REAL4 betaFac = 0;
  REAL4 normFac = 0;
  REAL4 betaii = 0;
  REAL4 snrFactor = 0;
  REAL4 chisqVal = 0;
  RingVetoResults *Result = NULL;
  REAL4Vector *betaji = NULL;
  INT4Vector *tau = NULL;
  REAL4Vector *betajj = NULL;
  betaji = XLALCreateREAL4Vector(thisResult->numResults);
  betajj = XLALCreateREAL4Vector(thisResult->numResults);
  tau = XLALCreateINT4Vector(thisResult->numResults);

  Result = firstResult;

  for(j=0; j < thisResult->numResults; j++){
    betaji->data[j] = RVCC->Beta->data[(thisResult->numResults)*j+i];
    betajj->data[j] = RVCC->Beta->data[(thisResult->numResults)*j+j];
    tau->data[j] = RVCC->Tau->data[(thisResult->numResults)*j+i];
    /*verbose("beta[%d] %f, tau[%d] %d\n",j,beta->data[j],j,tau->data[j]);*/
    }
  betaii = betaji->data[i];

  Result = firstResult;

  chisqVal = 0.0;
  Result = firstResult;
  j = 0;
  snrFactor = 1.0 * params->dynRangeFac / thisResult->sigma;
  for (j=0;j < (thisResult->numResults);j++){
    REAL4 otherSNRfactor  = 1.0 * params->dynRangeFac / Result->sigma;
    betaFac = sqrt(betajj->data[j]*betaii)/betaji->data[j];
    if (i != j){
      normFac = betaji->data[j]/betaii - betajj->data[j]/betaji->data[j];
      chisqVal += 1.0*
        (fabs(thisResult->result.data->data[dataIX])*snrFactor-
        1.0*betaFac*fabs(Result->result.data->data[dataIX+tau->data[j]])
        *otherSNRfactor)/normFac*
        (fabs(thisResult->result.data->data[dataIX])*snrFactor-
        1.0*betaFac*fabs(Result->result.data->data[dataIX+tau->data[j]])
        *otherSNRfactor)/normFac/
        (1.0 + 2.0*deltasq*fabs(Result->result.data->data[dataIX+tau->data[j]])
        *otherSNRfactor + 
        deltasq*deltasq*fabs(Result->result.data->data[dataIX+tau->data[j]])
        *otherSNRfactor*fabs(Result->result.data->data[dataIX+tau->data[j]])
        *otherSNRfactor) ;
      }
    Result=Result->next;
    }
  chisqVal /= (thisResult->numResults-1)/(1.0+2.0*deltasq+deltasq*deltasq); 
  /* Divide by DOF */  
  XLALDestroyREAL4Vector(betajj);
  XLALDestroyREAL4Vector(betaji);
  XLALDestroyINT4Vector(tau);
  return chisqVal;
  }

 
