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
#include "ring.h"

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
    SnglRingdownTable        *events,
    UINT4                    *numEvents,
    REAL4TimeSeries          *result,
    REAL4                     tmpltSigma,
    REAL4                     tmpltFrequency,
    REAL4                     tmpltQuality,
    struct ring_params       *params
    );


SnglRingdownTable * ring_filter(
    RingDataSegments         *segments,
    RingTemplateBank         *bank,
    REAL4FrequencySeries     *invSpectrum,
    REAL4FFTPlan             *fwdPlan,
    REAL4FFTPlan             *revPlan,
    struct ring_params       *params
    )
{
  SnglRingdownTable       *events = NULL; /* head of linked list of events */
  REAL4TimeSeries          signal;
  REAL4TimeSeries          result;
  COMPLEX8FrequencySeries  stilde;
  COMPLEX8FrequencySeries  rtilde;

  UINT4 segmentLength;
  UINT4 sgmnt;
  UINT4 tmplt;

  if ( ! params->doFilter )
    return NULL;

  segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );

  memset( &signal, 0, sizeof( signal ) );
  memset( &result, 0, sizeof( result ) );
  memset( &stilde, 0, sizeof( stilde ) );
  memset( &rtilde, 0, sizeof( rtilde ) );

  signal.deltaT = 1.0/params->sampleRate;
  signal.sampleUnits = lalStrainUnit;
  rtilde.deltaF = 1.0 / params->segmentDuration;
  signal.data = XLALCreateREAL4Vector( segmentLength );
  stilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  rtilde.data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );
  result.data = XLALCreateREAL4Vector( segmentLength );

  /* loop over all elements in the template bank */
  for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt )
  {
    RingTemplateInput *thisTmplt = bank->tmplt + tmplt;
    UINT4 numEvents = 0;
    REAL8 sigma;
    
    verbose( "creating template %d\n", tmplt );

    /* make template and fft it */
    XLALComputeRingTemplate( &signal, thisTmplt );
    LALSnprintf( signal.name, sizeof(signal.name), "TMPLT_%u", tmplt );
    XLALREAL4TimeFreqFFT( &stilde, &signal, fwdPlan );
    LALSnprintf( stilde.name, sizeof(stilde.name), "TMPLT_%u_FFT", tmplt );

    /* compute sigma for this template */
    sigma = sqrt( compute_template_variance( &stilde, invSpectrum,
          params->dynRangeFac ) );

    /* loop over segments */
    for ( sgmnt = 0; sgmnt < segments->numSgmnt; ++sgmnt )
    {
      verbose( "  filtering segment %d against template %d\n", sgmnt, tmplt );

      /* filter the segment with the template */
      filter_segment_template( &result, &rtilde, &stilde,
          &segments->sgmnt[sgmnt], revPlan );

      /* search through results for threshold crossings and record events */
      events = find_events( events, &numEvents, &result, sigma,
          thisTmplt->frequency, thisTmplt->quality, params );

      /* write filter output if requested */
      if ( params->writeFilterOutput )
      { /* guess we better normalize it so it is SNR-like... */
        REAL4 snrFactor = 2 * params->dynRangeFac / sigma;
        UINT4 k;
        for ( k = 0; k < result.data->length; ++k )
          result.data->data[k] *= snrFactor;
        write_REAL4TimeSeries( &result );
      }
    }
    params->numEvents += numEvents;
    verbose( "found %u event%s in template %u\n", numEvents,
        numEvents == 1 ? "" : "s", tmplt );
  }

  verbose( "found %u event%s\n", params->numEvents,
      params->numEvents == 1 ? "" : "s" );

  /* cleanup */
  XLALDestroyREAL4Vector( result.data );
  XLALDestroyCOMPLEX8Vector( rtilde.data );
  XLALDestroyCOMPLEX8Vector( stilde.data );
  XLALDestroyREAL4Vector( signal.data );

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
  LALSnprintf( rtilde->name, sizeof( rtilde->name ), "%s_%s",
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

  /* inverse fft to obtain result */
  XLALREAL4FreqTimeFFT( result, rtilde, plan );
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
    struct ring_params *params
    )
{
  const REAL4 efolds = 10.0; /* number of efolds of ringdown in template */
  SnglRingdownTable *thisEvent = events; /* the current event */
  REAL4 snrFactor; /* factor to convert from filter result to snr */
  REAL4 threshold; /* modified threshold on filter result (rather than snr) */
  REAL4 filterDuration;
  INT8  lastTimeNS;
  INT8  gapTimeNS;
  UINT4 segmentStride;
  UINT4 eventCount = 0;
  UINT4 jmin;
  UINT4 jmax;
  UINT4 j;
  LALStatus             status = blank_status;
  LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };

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

  /* compute start and stop index for scanning time series */
  segmentStride = floor( params->strideDuration / result->deltaT + 0.5 );
  jmin = segmentStride/2;
  jmax = jmin + segmentStride;

  for ( j = jmin; j < jmax; ++j )
    if ( fabs( result->data->data[j] ) > threshold ) /* threshold crossing */
    {
      REAL4 snr;
      INT8  timeNS;
      REAL4 sigma;
      REAL4 amp;

      snr     = fabs( result->data->data[j] ) * snrFactor;
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
        strncpy( thisEvent->search, "ring", sizeof( thisEvent->search ) );
        strncpy( thisEvent->channel, params->channel, sizeof( thisEvent->channel ) );
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
#if 0        
        amp = 2.415e-21; /* factor given in PRD 022001 (1999) */
        amp *= sqrt( 2.0 * LAL_PI ); /* convert NEW conventions to OLD conventions */
        amp *= pow( LAL_C_SI, 3.0) / LAL_G_SI / LAL_MSUN_SI / 2.0 / LAL_PI; 
        amp *= 1.0 / sqrt( tmpltQuality ) * sqrt( 1 - 0.63 * pow( (tmpltQuality / 2.0), (-2.0/3.0) ) ) / tmpltFrequency ;
#endif
                               /* aplitude for ringdown at a distance of 1Mpc with epsilon=0.01 */
        
        amp = sqrt(5.0*0.01 / 2.0) / LAL_TWOPI * LAL_C_SI * 
          sqrt( 1.0-0.63*pow( (tmpltQuality / 2.0), (-2.0/3.0) )) / tmpltFrequency /
          1e6 / LAL_PC_SI / sqrt( tmpltQuality ) * pow (1+7/24/pow(tmpltQuality,2),-0.5);
          
        sigma=tmpltSigma * amp;
        thisEvent->sigma_sq = pow(sigma, 2.0);
        thisEvent->eff_dist = sigma / thisEvent->snr;
        
        ++eventCount;
      }
      else if ( snr > thisEvent->snr ) /* maximize within a set of crossings */
      {
        /* update to specific information about this threshold crossing */
        ns_to_epoch( &thisEvent->start_time, timeNS );
        thisEvent->snr        = snr;

#if 0        
        amp = 2.415e-21; /* factor given in PRD 022001 (1999) */
        amp *= sqrt( 2.0 * LAL_PI ); /* convert NEW conventions to OLD conventions */
        amp *= pow( LAL_C_SI, 3.0) / LAL_G_SI / LAL_MSUN_SI / 2.0 / LAL_PI;
        amp *= 1.0 / sqrt( tmpltQuality ) * sqrt( 1.0 - 0.63 * pow( (tmpltQuality / 2.0), (-2.0/3.0) ) ) / tmpltFrequency ;
#endif
        amp = sqrt(5.0*0.01 / 2.0) / LAL_TWOPI * LAL_C_SI *
          sqrt( 1.0-0.63*pow( (tmpltQuality / 2.0), (-2.0/3.0) )) / tmpltFrequency /
          1e6 / LAL_PC_SI / sqrt( tmpltQuality ) * pow (1+7/24/pow(tmpltQuality,2),-0.5);
        
        sigma=tmpltSigma * amp;
        thisEvent->eff_dist = sigma / thisEvent->snr;        
      }

      /* update last threshold crossing time */
      lastTimeNS = timeNS;
    }

  *numEvents += eventCount;
  return events;
}
