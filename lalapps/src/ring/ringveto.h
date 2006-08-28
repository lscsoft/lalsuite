#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>

enum { write_frame, write_ascii };

struct ring_params {
  char        *programName;
  char        *cvsRevision;
  char        *cvsSource;
  char        *cvsDate;
  char         ifoName[3];
  INT4         randomSeed;
  LIGOTimeGPS  startTime;
  LIGOTimeGPS  endTime;
  INT8         trigStartTimeNS;
  INT8         trigEndTimeNS;
  LIGOTimeGPS  frameDataStartTime;
  REAL8        frameDataDuration;
  REAL8        duration;
  const char  *channel;
  const char  *calibCache;
  const char  *dataCache;
  const char  *injectFile;
  REAL8        sampleRate;
  REAL8        padData;
  REAL8        segmentDuration;
  REAL8        strideDuration;
  REAL8        truncateDuration;
  UINT4        numOverlapSegments;
  REAL4        dynRangeFac;
  REAL4        lowCutoffFrequency;
  REAL4        highpassFrequency;
  REAL4        invSpecLen;
  char         bankFile[256];
  RingTemplateBankInput bankParams;
  REAL4        bankMinFrequency;
  REAL4        bankMaxFrequency;
  REAL4        bankMinQuality;
  REAL4        bankMaxQuality;
  REAL4        bankMaxMismatch;
  REAL4        bankTemplatePhase;
  REAL4        threshold;
  REAL4        maximizeEventDuration;
  const char  *segmentsToDoList;
  const char  *templatesToDoList;
  UINT4        numEvents;
  char         outputFile[256];
  int          geoData;
  REAL8        geoHighpassFrequency;
  REAL8        geoScale;
  /* flags */
  int          strainData;
  int          simData;
  int          zeroData;
  int          whiteSpectrum;
  int          bankOnly;
  int          getData;
  int          getResponse;
  int          getSpectrum;
  int          getBank;
  int          doFilter;
  /* write intermediate result flags */
  int          writeRawData;
  int          writeProcessedData;
  int          writeResponse;
  int          writeSpectrum;
  int          writeInvSpectrum;
  int          writeSegment;
  int          writeFilterOutput;
};

typedef struct tagRingDataSegments
{
  UINT4                    numSgmnt;
  COMPLEX8FrequencySeries *sgmnt;
}
RingDataSegments;

/* Ring Veto Structure */
struct tagRingVetoTemplateBank {
  RingTemplateBank *bank;
  struct tagRingVetoTemplateBank *next;
}; 
typedef struct tagRingVetoTemplateBank RingVetoTemplateBank; 

/* Ring veto results structure - I need to store ~10 result vectors */
/* and the template variances in one structure in order to compute */
/* the cross correlations and construct the chisq bank time series */

struct tagRingVetoResults {
  REAL4TimeSeries            result;
  COMPLEX8TimeSeries            tmpSnip;
  COMPLEX8FrequencySeries    tmpSnipTilde;
  REAL8			     sigma;
  REAL8			     quality;
  REAL8			     frequency;
  UINT4		             numResults;			
  struct tagRingVetoResults  *next;
};
typedef struct tagRingVetoResults RingVetoResults;

struct tagRingVetoCC {
  REAL4Vector                 Beta;
  INT4Vector		      Tau;
};
typedef struct tagRingVetoCC RingVetoCC;

/* routines in ring_option */
int ring_parse_options( struct ring_params *params, int argc, char **argv );
int ring_params_sanity_check( struct ring_params *params );

/* routines in ring_output */
ProcessParamsTable * create_process_params( int argc, char **argv,
    const char *program );
int ring_output_events_xml( 
    char               *outputFile,
    SnglRingdownTable  *events,
    ProcessParamsTable *processParamsTable,
    struct ring_params *params
    );

/* routines to write intermediate results in ring_output */
int write_REAL4TimeSeries( REAL4TimeSeries *series );
int write_REAL4FrequencySeries( REAL4FrequencySeries *series );
int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series );
int write_bank( RingTemplateBank *bank );

/* routines in ring_filter */
SnglRingdownTable * ringveto_filter(
    RingDataSegments         *segments,
    RingVetoTemplateBank         *bank,
    REAL4FrequencySeries     *invSpectrum,
    REAL4FFTPlan             *fwdPlan,
    REAL4FFTPlan             *revPlan,
    struct ring_params       *params
    );


