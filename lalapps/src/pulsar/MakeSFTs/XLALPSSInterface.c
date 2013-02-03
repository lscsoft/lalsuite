#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <lal/XLALError.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include "XLALPSSInterface.h"
#include "pss_serv.h"

/* variables in PSS that aren't exported in the headers */
extern FILE*LOG_INFO;

static const char* captionString = "count";

FILE* XLALPSSOpenLog(const char*name) {
  if (strcmp(name, "-") == 0)
    LOG_INFO = stderr;
  else
    LOG_INFO = logfile_open(name);
  return LOG_INFO;
}

void XLALPSSCloseLog(void) {
  if(LOG_INFO != stderr) 
    logfile_close(LOG_INFO);
}

PSSEventParams *XLALCreatePSSEventParams(UINT4 length, XLALPSSParamSet setpar) {
  PSSEventParams*ep = NULL;
  if ( length == 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EINVAL ); 
  ep = crea_evenparam(length,0);
  if ( ep == NULL )
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EFAULT );

  /* values that we might set differently from the defaults in PSS' crea_evenpar() */
  // PSS default // Paola used /* PSS paramter comment */

  // 1.0   // 0.0f           /* 1.0f = use absolute values. Else use signed values */
  if(setpar.set & XLALPSS_SET_ABS) {
    fprintf(stderr, "overriding PSS default %f for abs: %.13f\n", ep->absvalue, setpar.abs);
    ep->absvalue = setpar.abs;
  }

  // 600   // 20.0f          /* memory time of the autoregressive average */
  if(setpar.set & XLALPSS_SET_TAU) {
    fprintf(stderr, "overriding PSS default %f for tau: %.13f\n", ep->tau, setpar.tau);
    ep->tau = setpar.tau;
  }

  // 20    // ep->tau/(6.103515625e-05f); /* e.g. 10 or 20 : when re-evaluate mean and std */
  if(setpar.set & XLALPSS_SET_FACT) {
    fprintf(stderr, "overriding PSS default %f for tau: %.13f\n", ep->factor, setpar.fact);
    ep->factor = setpar.fact;
  }

  // 6.0   // 5.0f;          /* critical ratio of the threshold */
  if(setpar.set & XLALPSS_SET_CR) {
    fprintf(stderr, "overriding PSS default %f for cr: %.13f\n", ep->cr, setpar.cr);
    ep->cr = setpar.cr;
  }

  // 0.15  // 0.00061035f;   /* how many seconds around (before and after) the event have to be "purged" */
  if(setpar.set & XLALPSS_SET_EDGE) {
    fprintf(stderr, "overriding PSS default %f for edge: %.13f\n", ep->edge, setpar.edge);
    ep->edge = setpar.edge;
  }

  // 0.00001 // 1e-25f       /* value added to denominator to avoid division by zero */
  ep->notzero = 1e-25f;
  ep->w_norm=1.0f;           /* to be sure, might already be set in crea_evenparam() */

  return ep;
}

void XLALDestroyPSSEventParams(PSSEventParams *ts) {
  /* not using (X)LALFree because PSS isn't using (X)LALMalloc */
  if ( !ts )
    XLAL_ERROR_VOID( "XLALDestroyPSSEventParams", XLAL_EFAULT );
  if(ts->xamed)
    free(ts->xamed);
  if(ts->xastd)
    free(ts->xastd);
  if(ts->begin)
    free(ts->begin);
  if(ts->duration)
    free(ts->duration);
  if(ts->imax)
    free(ts->imax);
  if(ts->crmax)
    free(ts->crmax);
  if(ts->ener)
    free(ts->ener);
  free(ts);
}

REAL4TimeSeries *XLALConvertREAL8TimeSeriesToREAL4TimeSeries(REAL4TimeSeries *r4ts, REAL8TimeSeries *r8ts) {
  UINT4 i;
  if (r4ts->data->length != r8ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertREAL8TimeSeriesToREAL4TimeSeries", XLAL_EINVAL ); 
  for(i = 0; i < r8ts->data->length; i++)
    r4ts->data->data[i] = r8ts->data->data[i];
  strncpy(r4ts->name, r8ts->name, LALNameLength);
  r4ts->epoch       = r8ts->epoch;
  r4ts->deltaT      = r8ts->deltaT;
  r4ts->f0          = r8ts->f0;
  r4ts->sampleUnits = r8ts->sampleUnits;
  return r4ts;
}

PSSTimeseries *XLALCreatePSSTimeseries(UINT4 length) {
  PSSTimeseries *pssGD;
  if ( length == 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EINVAL );
  pssGD = crea_gd(length,0,0,captionString); /* unit of a timeseries is 'count' */
  if( !pssGD )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
  return(pssGD);
}

void XLALDestroyPSSTimeseries(PSSTimeseries *ts) {
  /* not using (X)LALFree because PSS isn't using (X)LALMalloc */
  if ( !ts )
    XLAL_ERROR_VOID( "XLALDestroyPSSTimeseries", XLAL_EFAULT );
  if(ts->x)
    free(ts->x);
  if(ts->y)
    free(ts->y);
  free(ts);
}

PSSHeaderParams* XLALPSSInitializeHeaderParams(PSSHeaderParams* headerParams, REAL8 deltaT) {
  /*
    we don't actually use the header params here to write SFTs,
    we only need them as an argument to the PSS functions for
    time-domain cleaning. Note that this initialization 
    probably doesn't work for any other application.
    
    The PSS function to allocate and initialize the header parameters is crea_sfdbheader()
  */

  memset(headerParams,0,sizeof(PSSHeaderParams));
  /*
    this implies
    headerParams.typ = 0;
    headerParams.nfft = 0;
  */
  headerParams->tsamplu = deltaT;

  return(headerParams);
}

PSSTimeseries *XLALPrintPSSTimeseriesToFile(PSSTimeseries *tsPSS, const char*name, UINT4 numToPrint) {
  UINT4 i,n,len;
  FILE*fp;
  LIGOTimeGPS gpsepoch;

  if ( !tsPSS || !name )
    XLAL_ERROR_NULL( "XLALPrintPSSTimeseriesToFile", XLAL_EFAULT );

  fp = fopen(name,"w");
  if ( !fp ) {
    fprintf(stderr,"XLALPrintPSSTimeseriesToFile: can't open file '%s'(%d)\n",name,errno);
    XLAL_ERROR_NULL( "XLALPrintPSSTimeseriesToFile", XLAL_EFAULT );
  }

  XLALGPSSetREAL8(&gpsepoch, tsPSS->ini);
  if (xlalErrno) {
    fprintf(stderr,"XLALPrintPSSTimeseriesToFile: unhandled XLAL Error in XLALGPSSetREAL8 %s,%d\n",__FILE__,__LINE__);
    XLAL_ERROR_NULL( "XLALPrintPSSTimeseriesToFile", XLAL_EFAULT );
  }

  len = tsPSS->n;
  fprintf(fp,"%% Length: %d\n",len);
  fprintf(fp,"%% deltaT: %f\n",tsPSS->dx);
  fprintf(fp,"%% gpsSeconds: %i\n",gpsepoch.gpsSeconds);
  fprintf(fp,"%% gpsNanoSeconds: %i\n",gpsepoch.gpsNanoSeconds);
  if(tsPSS->name)
    fprintf(fp,"%% Name: '%s'\n",tsPSS->name);
  if(tsPSS->capt)
    fprintf(fp,"%% Capt: '%s'\n",tsPSS->capt);
  n = numToPrint;
  if ((n == 0) || (2*n > len))
    n = len;
  fprintf(fp,"%% First %d values:\n", n);
  for(i = 0; i < n; i++)
    fprintf(fp,"%23.16e\n",tsPSS->y[i]);
  if(i + n < len) {
    fprintf(fp,"%% Last %d values:\n", n);
    for(i = len - n; i < len; i++)
      fprintf(fp,"%23.16e\n",tsPSS->y[i]);
  }
  fclose(fp);

  xlalErrno = 0;
  return tsPSS;
}

REAL8TimeSeries *XLALPrintREAL8TimeSeriesToFile(REAL8TimeSeries *ts, const char*name, UINT4 numToPrint, BOOLEAN scaleToREAL4) {
  UINT4 i,n,len;
  FILE*fp;
  char unit[LALNameLength];

  if ( !ts || !name )
    XLAL_ERROR_NULL( "XLALPrintREAL8TimeSeriesToFile", XLAL_EFAULT );

  fp = fopen(name,"w");
  if ( !fp ) {
    fprintf(stderr,"XLALPrintREAL8TimeSeriesToFile: can't open file '%s'(%d)\n",name,errno);
    XLAL_ERROR_NULL( "XLALPrintREAL8TimeSeriesToFile", XLAL_EFAULT );
  }

  len = ts->data->length;
  fprintf(fp,"%% Length: %d\n",len);
  fprintf(fp,"%% deltaT: %f\n",ts->deltaT);
  fprintf(fp,"%% gpsSeconds: %i\n",ts->epoch.gpsSeconds);
  fprintf(fp,"%% gpsNanoSeconds: %i\n",ts->epoch.gpsNanoSeconds);
  if(ts->name)
    fprintf(fp,"%% Name: '%s'\n",ts->name);
  if( XLALUnitAsString( unit, LALNameLength, &(ts->sampleUnits) ) )
    fprintf(fp,"%% Unit: '%s'\n",unit);
  n = numToPrint;
  if ((n == 0) || (2*n > len))
    n = len;
  fprintf(fp,"%% First %d values:\n", n);
  for(i = 0; i < n; i++)
    if(scaleToREAL4) {
      REAL4 val = ts->data->data[i];
      fprintf(fp,"%23.16e\n",val);
    } else
      fprintf(fp,"%23.16e\n",ts->data->data[i]);

  if(i + n < len) {
    fprintf(fp,"%% Last %d values:\n", n);
    for(i = len - n; i < len; i++)
      if(scaleToREAL4) {
	REAL4 val = ts->data->data[i];
	fprintf(fp,"%23.16e\n",val);
      } else
	fprintf(fp,"%23.16e\n",ts->data->data[i]);
  }
  fclose(fp);

  xlalErrno = 0;
  return ts;
}

REAL4TimeSeries *XLALPrintREAL4TimeSeriesToFile(REAL4TimeSeries *ts, const char*name, UINT4 numToPrint) {
  UINT4 i,n,len;
  FILE*fp;
  char unit[LALNameLength];

  if ( !ts || !name )
    XLAL_ERROR_NULL( "XLALPrintREAL4TimeSeriesToFile", XLAL_EFAULT );

  fp = fopen(name,"w");
  if ( !fp ) {
    fprintf(stderr,"XLALPrintREAL4TimeSeriesToFile: can't open file '%s'(%d)\n",name,errno);
    XLAL_ERROR_NULL( "XLALPrintREAL4TimeSeriesToFile", XLAL_EFAULT );
  }

  len = ts->data->length;
  fprintf(fp,"%% Length: %d\n",len);
  fprintf(fp,"%% deltaT: %f\n",ts->deltaT);
  fprintf(fp,"%% gpsSeconds: %i\n",ts->epoch.gpsSeconds);
  fprintf(fp,"%% gpsNanoSeconds: %i\n",ts->epoch.gpsNanoSeconds);
  if(ts->name)
    fprintf(fp,"%% Name: '%s'\n",ts->name);
  if( XLALUnitAsString( unit, LALNameLength, &(ts->sampleUnits) ) )
    fprintf(fp,"%% Unit: '%s'\n",unit);
  n = numToPrint;
  if ((n == 0) || (2*n > len))
    n = len;
  fprintf(fp,"%% First %d values:\n", n);
  for(i = 0; i < n; i++)
    fprintf(fp,"%23.16e\n",ts->data->data[i]);

  if(i + n < len) {
    fprintf(fp,"%% Last %d values:\n", n);
    for(i = len - n; i < len; i++)
      fprintf(fp,"%23.16e\n",ts->data->data[i]);
  }
  fclose(fp);

  xlalErrno = 0;
  return ts;
}

REAL8TimeSeries *XLALConvertPSSTimeseriesToREAL8Timeseries(REAL8TimeSeries *ts, PSSTimeseries *tsPSS) {
  UINT4 i;

  if ( !tsPSS || !ts )
    XLAL_ERROR_NULL( "XLALConvertPSSTimeseriesToREAL8Timeseries", XLAL_EFAULT );
  if (tsPSS->n != (int)ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertPSSTimeseriesToREAL8Timeseries", XLAL_EINVAL ); 

  for(i = 0; i < ts->data->length; i++)
    ts->data->data[i] = tsPSS->y[i];

  strncpy(ts->name, tsPSS->name, LALNameLength);
  XLALGPSSetREAL8(&(ts->epoch), tsPSS->ini);
  if (xlalErrno)
    fprintf(stderr,"XLALConvertPSSTimeseriesToREAL8Timeseries: unhandled XLAL Error in XLALGPSSetREAL8 %s,%d\n",__FILE__,__LINE__);
  ts->deltaT = tsPSS->dx;
  ts->f0     = 0.0; /* no heterodyning */
  XLALParseUnitString( &(ts->sampleUnits), tsPSS->capt );
  if (xlalErrno) {
    fprintf(stderr,"XLALConvertPSSTimeseriesToREAL8Timeseries: unhandled XLAL Error in XLALParseUnitString:%d %s,%d\n",xlalErrno,__FILE__,__LINE__);
    if (tsPSS->capt)
      fprintf(stderr,"[DEBUG] PSS caption: '%s'\n", tsPSS->capt); 
    xlalErrno = 0;
  }

  return ts;
}

PSSTimeseries *XLALConvertREAL8TimeseriesToPSSTimeseries(PSSTimeseries *tsPSS, REAL8TimeSeries *ts) {
  UINT4 i;

  /* input sanity checking */
  if ( !tsPSS || !ts )
    XLAL_ERROR_NULL( "XLALConvertREAL8TimeseriesToPSSTimeseries", XLAL_EFAULT );
  if (tsPSS->nall != (int)ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertREAL8TimeseriesToPSSTimeseries", XLAL_EINVAL ); 

  /* convert the actual data */
  for(i = 0; i < ts->data->length; i++)
    tsPSS->y[i] = ts->data->data[i];
  tsPSS->n = ts->data->length;

  /* other parameters */
  strncpy(tsPSS->name, ts->name, sizeof(tsPSS->name)); /* sizeof(tsPSS->name) = 20 */
  tsPSS->ini = XLALGPSGetREAL8(&(ts->epoch));
  tsPSS->dx = ts->deltaT;
  tsPSS->capt = captionString;

  return tsPSS;
}

REAL4TimeSeries *XLALConvertPSSTimeseriesToREAL4Timeseries(REAL4TimeSeries *ts, PSSTimeseries *tsPSS) {
  UINT4 i;

  if ( !tsPSS || !ts )
    XLAL_ERROR_NULL( "XLALConvertPSSTimeseriesToREAL4Timeseries", XLAL_EFAULT );
  if (tsPSS->n != (int)ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertPSSTimeseriesToREAL4Timeseries", XLAL_EINVAL ); 

  for(i = 0; i < ts->data->length; i++)
    ts->data->data[i] = tsPSS->y[i];

  strncpy(ts->name, tsPSS->name, LALNameLength);
  XLALGPSSetREAL8(&(ts->epoch), tsPSS->ini);
  ts->deltaT = tsPSS->dx;
  ts->f0     = 0.0; /* no heterodyning */
  XLALParseUnitString( &(ts->sampleUnits), tsPSS->capt );

  return ts;
}

PSSTimeseries *XLALConvertREAL4TimeseriesToPSSTimeseries(PSSTimeseries *tsPSS, REAL4TimeSeries *ts) {
  UINT4 i;

  /* input sanity checking */
  if ( !tsPSS || !ts )
    XLAL_ERROR_NULL( "XLALConvertREAL4TimeseriesToPSSTimeseries", XLAL_EFAULT );
  if (tsPSS->nall != (int)ts->data->length)
    XLAL_ERROR_NULL( "XLALConvertREAL4TimeseriesToPSSTimeseries", XLAL_EINVAL ); 

  /* convert the actual data */
  for(i = 0; i < ts->data->length; i++)
    tsPSS->y[i] = ts->data->data[i];
  tsPSS->n = ts->data->length;

  /* other parameters */
  strncpy(tsPSS->name, ts->name, sizeof(tsPSS->name)); /* sizeof(tsPSS->name) = 20 */
  tsPSS->ini = XLALGPSGetREAL8(&(ts->epoch));
  tsPSS->dx = ts->deltaT;
  tsPSS->capt = captionString;

  return tsPSS;
}

PSSTimeseries *XLALPSSHighpassData(PSSTimeseries *tsout, PSSTimeseries *tsin, PSSHeaderParams* hp, REAL4 f) {
  if ( !tsout || !tsin || !hp )
    XLAL_ERROR_NULL( "XLALPSSHighpassData", XLAL_EFAULT );
  fprintf(stderr, "[DEBUG] highpass_data_bil: %d\n", highpass_data_bil(tsout,tsin,hp,f) );
  return tsout;
}

PSSEventParams *XLALPSSComputeARMeanAndStdev(PSSEventParams *events, PSSTimeseries *ts, PSSHeaderParams* hp) {
  long ret;
  if ( !events || !ts || !hp)
    XLAL_ERROR_NULL( "XLALPSSComputeARMeanAndStdev", XLAL_EFAULT );
  if ((ret = sn_medsig(ts,events,hp)) != ts->n ) {
    fprintf(stderr, "[DEBUG] sn_medsig: %ld\n", ret);
    XLAL_ERROR_NULL( "XLALPSSComputeARMeanAndStdev", XLAL_EFUNC );
  }
  return events;
}

PSSEventParams *XLALPSSComputeExtARMeanAndStdev(PSSEventParams *events, PSSTimeseries *ts, PSSHeaderParams* hp) {
  long ret;
  unsigned long len, newlen, pre, i;
  float *xamed, *xastd, *tsdata;

  if ( !events || !ts || !hp)
    XLAL_ERROR_NULL( "XLALPSSComputeExtARMeanAndStdev", XLAL_EFAULT );

  pre = ceil(events->tau / ts->dx);
  len = ts->n;
  newlen = len+pre;

  /* allocate extended data, save original data in local pointers */

  xamed = events->xamed;
  events->xamed = (float*)XLALMalloc(newlen*sizeof(float));
  if( !events->xastd ) {
    events->xamed = xamed;
    XLAL_ERROR_NULL( "XLALPSSComputeExtARMeanAndStdev", XLAL_ENOMEM );
  }

  xastd = events->xastd;
  events->xastd = (float*)XLALMalloc(newlen*sizeof(float));
  if( !events->xastd ) {
    events->xastd = xastd;
    XLALFree(events->xamed);
    events->xamed = xamed;
    XLAL_ERROR_NULL( "XLALPSSComputeExtARMeanAndStdev", XLAL_ENOMEM );
  }

  tsdata = ts->y;
  ts->y = (float*)XLALMalloc(newlen*sizeof(float));
  if( !ts->y ) {
    ts->y = tsdata;
    XLALFree(events->xamed);
    events->xamed = xamed;
    XLALFree(events->xastd);
    events->xastd = xastd;
    XLAL_ERROR_NULL( "XLALPSSComputeExtARMeanAndStdev", XLAL_ENOMEM );
  }
  ts->n = newlen;

  /* fill extended timeseries */
  for(i=0;i<pre;i++)
    ts->y[i] = tsdata[pre-i-1];
  for(i=0;i<len;i++)
    ts->y[i+pre] = tsdata[i];

  /* call original PSS routine */
  if ((ret = sn_medsig(ts,events,hp)) != ts->n ) {
    fprintf(stderr, "[DEBUG] sn_medsig: %ld\n", ret);
    /* cleanup extended timeseries */
    XLALFree(events->xamed);
    events->xamed = xamed;
    XLALFree(events->xastd);
    events->xastd = xastd;
    ts->n = len;
    XLALFree(ts->y);
    ts->y = tsdata;
    XLAL_ERROR_NULL( "XLALPSSComputeExtARMeanAndStdev", XLAL_EFUNC );
  }

  /* get results back into original structures */
  for(i=0;i<len;i++) {
    xamed[i] = events->xamed[i+pre];
    xastd[i] = events->xastd[i+pre];
  }

  /* free extended arrays and restore original pointers */
  XLALFree(events->xamed);
  events->xamed = xamed;
  XLALFree(events->xastd);
  events->xastd = xastd;
  ts->n = len;
  XLALFree(ts->y);
  ts->y = tsdata;

  return events;
}

PSSEventParams *XLALIdentifyPSSCleaningEvents(PSSEventParams *events, PSSTimeseries *ts) {
  long ret;
  if ( !events || !ts )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFAULT );
  if ((ret = even_anst(ts,events)) != ts->n ) {
    fprintf(stderr, "[DEBUG] even_anst: %ld\n", ret );
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFUNC );
  }
  return events;
}

PSSTimeseries *XLALSubstractPSSCleaningEvents(PSSTimeseries *tsout, PSSTimeseries *tsin,
					      PSSTimeseries *tshp, PSSEventParams *events,
					      PSSHeaderParams* hp) {
  if ( !tsout || !tsin || !tshp || !events )
    XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFAULT );
  fprintf(stderr, "[DEBUG] data_subtract: %d\n", data_subtract(tsout,tsin,tshp,events,hp) );
  /* XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFUNC ); */
  return tsout;
}

PSSTimeseries *XLALPSSPurgeEvents(PSSTimeseries *tsout, PSSTimeseries *tsin,
				  PSSTimeseries *tshp, PSSEventParams *events,
				  PSSHeaderParams* hp) {
  if ( !tsout || !tsin || !tshp || !events )
    XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFAULT );
  fprintf(stderr, "[DEBUG] purge_data_subtract: %d\n", purge_data_subtract(tsout,tsin,tshp,events,hp) );
  /* XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFUNC ); */
  return tsout;
}
