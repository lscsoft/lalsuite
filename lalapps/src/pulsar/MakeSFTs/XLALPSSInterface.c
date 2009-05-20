#include <stdio.h>
#include "XLALPSSInterface.h"

#define UINT4 unsigned long
#define UINT8 unsigned long long
#define REAL4 float
#define XLAL_EFAULT -1
#define XLAL_ERROR_NULL(a,b) fprintf(stderr,"XLAL Error: %s:%d\n",a,b);

/* as this code isn't used to actually write SFTs, the SFT header
   information that some PSS functions work on is kept private in this module */
PSSHeaderParams headerParams;

PSSEventParams *XLALCreatePSSEventParams(UINT4 length) { 
  PSSEventParams*ep = NULL;
  if ( length == 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EFAULT ); 
  ep = crea_evenparam(length);
  if ( ep == NULL )
    XLAL_ERROR_NULL( "XLALCreatePSSEventParams", XLAL_EFAULT ); 
  return ep;
} 


PSSTimeseries *XLALCreatePSSTimeseries(UINT4 length) {
  if ( length == 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
}

void XLALDestroyPSSTimeseries(PSSTimeseries *ts) {
  if ( ts == NULL )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
}

PSSTimeseries *XLALConvertPSSTimeseriesToREAL4Timeseries(void) {
}
PSSTimeseries *XLALConvertREAL4TimeseriesToPSSTimeseries(void) {
}

PSSTimeseries *XLALPSSHighpassData(PSSTimeseries *tsout, PSSTimeseries *tsin, REAL4 f) {
  int ret;
  ret = highpass_data_bil(tsout,tsin,&headerParams,f);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALCreatePSSTimeseries", XLAL_EFAULT );
  return tsout;
}

PSSEventParams *XLALIdentifyPSSCleaningEvents(PSSEventParams *events, PSSTimeseries *ts) {
  int ret;
  if ( !events || !ts )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFAULT );
  ret = sn_medsig(ts,events,&headerParams);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFAULT );
  ret = even_anst(ts,events);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALIdentifyPSSCleaningEvents", XLAL_EFAULT );
  return events;
}

PSSTimeseries *XLALSubstractPSSCleaningEvents(PSSTimeseries *tsout, PSSTimeseries *tsin,
					      PSSTimeseries *tshp, PSSEventParams *events) {
  int ret;
  if ( !tsout || !tsin || !tshp || !events )
    XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFAULT );
  ret = purge_data_subtract(tsout,tsin,tshp,events,&headerParams);
  if ( ret != 0 )
    XLAL_ERROR_NULL( "XLALSubstractPSSCleaningEvents", XLAL_EFAULT );
  return tsout;
}

main() {
  GD *gd=NULL, *gdhp=NULL, *gdc=NULL;
  EVEN_PARAM *ep=NULL;
  int res;

  ep = crea_evenparam(100);/*(long len)*/
  res = highpass_data_bil(gdhp,gd,&headerParams,1000.0);/*float freqc*/
  res = sn_medsig(gdhp,ep,&headerParams);
  res = even_anst(gdhp,ep);
  res = purge_data_subtract(gdc,gd,gdhp,ep,&headerParams);
}
