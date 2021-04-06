/* double inclusion protection */
#ifndef XLALPSSINTERFACE_H
#define XLALPSSINTERFACE_H

/* PSS has some strange internal dependencies between its headers,
   better don't change the following order */
#include <stdio.h>
#include "pss_math.h"
#include "pss_sfc.h"
#include "pss_snag.h"
#include "pss_sfdb.h"
#include <lal/LALDatatypes.h>

#define PSS_PI PIG

#define XLALPSS_SET_ABS  (1<<0)
#define XLALPSS_SET_TAU  (1<<1)
#define XLALPSS_SET_FACT (1<<2)
#define XLALPSS_SET_CR   (1<<3)
#define XLALPSS_SET_EDGE (1<<4)

typedef struct {
  UINT4 set; /* specify which PSS parameters to modify in XLALCreatePSSEventParams(),
		ORed list of XLALPSS_SET flags */
  REAL4 abs, tau, fact, cr, edge; /* new values for the PSS parameters */
} XLALPSSParamSet;

/* PSS interface datatypes based on PSS datatypes */
typedef GD           PSSTimeseries;
typedef EVEN_PARAM   PSSEventParams;
typedef HEADER_PARAM PSSHeaderParams;

/* creator & destructor functions */
extern PSSEventParams *XLALCreatePSSEventParams(UINT4 length, XLALPSSParamSet setpar);
extern PSSTimeseries  *XLALCreatePSSTimeseries(UINT4 length);
extern void XLALDestroyPSSTimeseries(PSSTimeseries *ts);
extern void XLALDestroyPSSEventParams(PSSEventParams *ts);
extern PSSHeaderParams* XLALPSSInitializeHeaderParams(PSSHeaderParams* headerParams, REAL8 deltaT);

/* open and close a log file */
/* the PSS functions definitely need a working file pointer they can log to */
/* XLALPSSOpenLog("-") logs to stderr */
extern FILE* XLALPSSOpenLog(const const char*name);
extern void XLALPSSCloseLog(void);

/* LAL <-> PSS timeseries conversion functions */
extern REAL8TimeSeries
*XLALConvertPSSTimeseriesToREAL8Timeseries
(REAL8TimeSeries *ts,
 PSSTimeseries *tsPSS);

extern PSSTimeseries
*XLALConvertREAL8TimeseriesToPSSTimeseries
(PSSTimeseries *tsPSS,
 REAL8TimeSeries *ts);

extern REAL4TimeSeries
*XLALConvertPSSTimeseriesToREAL4Timeseries
(REAL4TimeSeries *ts,
 PSSTimeseries *tsPSS);

extern PSSTimeseries
*XLALConvertREAL4TimeseriesToPSSTimeseries
(PSSTimeseries *tsPSS,
 REAL4TimeSeries *ts);

/* debug: convert a LAL REAL8 Timeseries to a REAL4 one */
REAL4TimeSeries
*XLALConvertREAL8TimeSeriesToREAL4TimeSeries
(REAL4TimeSeries *r4ts,
 REAL8TimeSeries *r8ts);

/* debug: write timeseries data to a file */
extern PSSTimeseries
*XLALPrintPSSTimeseriesToFile
(PSSTimeseries *tsPSS,
 const char*name,
 UINT4 numToPrint);

extern REAL8TimeSeries
*XLALPrintREAL8TimeSeriesToFile
(REAL8TimeSeries *ts,
 const char*name,
 UINT4 numToPrint,
 BOOLEAN scaleToREAL4);

extern REAL4TimeSeries
*XLALPrintREAL4TimeSeriesToFile
(REAL4TimeSeries *ts,
 const char*name,
 UINT4 numToPrint);

/* functions for time domain cleaning */
extern PSSTimeseries
*XLALPSSHighpassData
(PSSTimeseries *tsout,
 PSSTimeseries *tsin,
 PSSHeaderParams* hp,
 REAL4 f);

extern PSSEventParams
*XLALPSSComputeARMeanAndStdev
(PSSEventParams *events,
 PSSTimeseries *ts,
 PSSHeaderParams* hp);

extern PSSEventParams
*XLALPSSComputeExtARMeanAndStdev
(PSSEventParams *events,
 PSSTimeseries *ts,
 PSSHeaderParams* hp);

extern PSSEventParams
*XLALIdentifyPSSCleaningEvents
(PSSEventParams *events,
 PSSTimeseries *ts);

extern PSSTimeseries
*XLALSubstractPSSCleaningEvents
(PSSTimeseries *tsout,
 PSSTimeseries *tsin,
 PSSTimeseries *tshp,
 PSSEventParams *events,
 PSSHeaderParams* hp);

/* debug: just for testing backwards compatibility.
   This function should not be used anymore, instead use the sequence
   XLALPSSComputeARMeanAndStdev();
   XLALIdentifyPSSCleaningEvents();
   XLALSubstractPSSCleaningEvents();
*/
extern PSSTimeseries
*XLALPSSPurgeEvents
(PSSTimeseries *tsout,
 PSSTimeseries *tsin,
 PSSTimeseries *tshp,
 PSSEventParams *events,
 PSSHeaderParams* hp);

#endif /*< XLALPSSINTERFACE_H double inclusion protection */
