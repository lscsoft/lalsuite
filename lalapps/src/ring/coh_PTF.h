/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>

#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/Units.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/RealFFT.h>
#include <lal/FrameStream.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirpDatatypes.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpPTF.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LIGOLwXMLInspiralHeaders.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/TimeSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/FindChirpPTF.h>
#include <lal/RingUtils.h>
#include <LALAppsVCSInfo.h>
#include <lal/SkyCoordinates.h>
#include <lal/XLALError.h>

#include "lalapps.h"
#include "getdata.h"
#include "injsgnl.h"
#include "getresp.h"
#include "spectrm.h"
#include "segment.h"
#include "errutil.h"
#include "processtable.h"
#include "gpstime.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <lal/GSLSupport.h>

#define BUFFER_SIZE 256
#define FILENAME_SIZE 256
#define MAXIFO 4

/* macro for testing validity of a condition that prints an error if invalid */
#define sanity_check( condition ) \
  ( condition ? 0 : ( fputs( #condition " not satisfied\n", stderr ), error( "sanity check failed\n" ) ) )

/* macro is_option() (and friends): determine if string is an option */
/* i.e., does it start with "-[a-zA-Z]" or "--[a-zA-Z]" */
#define is_long_option(s) \
  ( strlen(s) > 2 && (s)[0] == '-' && (s)[1] == '-' && isalpha( s[2] ) )
#define is_short_option(s) \
  ( strlen(s) > 1 && (s)[0] == '-' && isalpha( s[1] ) )
#define is_option(s) ( is_long_option(s) || is_short_option(s) )

/* Supress unused variable warnings */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

extern int vrbflg;

enum { write_frame, write_ascii };

/* The coh PTF params structure */

struct coh_PTF_params {
  char        *programName;
  char        *cvsRevision;
  char        *cvsSource;
  char        *cvsDate;
  char         ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  UINT4        numIFO;
  INT4         randomSeed;
  INT4         haveTrig[LAL_NUM_IFO];
  LIGOTimeGPS  startTime;
  LIGOTimeGPS  endTime;
  LIGOTimeGPS  trigTime;
  INT8         trigStartTimeNS;
  INT8         trigEndTimeNS;
  LIGOTimeGPS  frameDataStartTime;
  REAL8        frameDataDuration;
  REAL8        duration;
  const char  *channel[LAL_NUM_IFO];
  const char  *dataCache[LAL_NUM_IFO];
  const char  *injectFile;
  REAL8        sampleRate;
  REAL8        padData;
  REAL8        segmentDuration;
  REAL8        psdSegmentDuration;
  REAL8        strideDuration;
  REAL8        psdStrideDuration;
  REAL8        truncateDuration;
  UINT4        numTimePoints;
  UINT4        numFreqPoints;
  UINT4        analStartPoint;
  REAL8        analStartTime;
  UINT4        analEndPoint;
  REAL8        analEndTime;
  UINT4        analStartPointBuf;
  UINT4        analEndPointBuf;
  UINT4        numAnalPoints;
  UINT4        numAnalPointsBuf;
  UINT4        numBufferPoints;
  REAL4        maxTempLength;
  UINT4        numOverlapSegments;
  REAL4        dynRangeFac;
  REAL4        lowTemplateFrequency;
  REAL4        lowFilterFrequency;
  REAL4        highFilterFrequency;
  REAL4        highpassFrequency;
  Approximant  approximant;
  LALPNOrder   order;
  REAL4        invSpecLen;
  REAL4        threshold;
  REAL4        spinThreshold;
  REAL4        snglSNRThreshold;
  REAL4        timeWindow;
  REAL4        spinSNR2threshold;
  REAL4        nonspinSNR2threshold;
  REAL4        rightAscension;
  REAL4        declination;
  REAL4        skyError;
  REAL4        timingAccuracy;
  const char  *skyPositionsFile;
  UINT4        skyLooping;
  const char  *bankFile;
  const char  *segmentsToDoList;
  const char  *templatesToDoList;
  UINT4        numEvents;
  UINT4        BVsubBankSize;
  const char  *bankVetoBankName;
  UINT4        numAutoPoints;
  REAL4        autoVetoTimeStep;
  REAL4        bankVeton;
  REAL4        bankVetoq;
  REAL4        autoVeton;
  REAL4        autoVetoq;
  REAL4        nullStatThreshold;
  REAL4        nullStatGradOn;
  REAL4        nullStatGradient;
  UINT4        numChiSquareBins;
  REAL4        chiSquareCalcThreshold;
  char         outputFile[256];
  UINT4        spinBank;
  char         spinBankName[256];
  UINT4        noSpinBank;
  char         noSpinBankName[256];
  char         userTag[256];
  char         ifoTag[256];
  UINT4        slideSegments[LAL_NUM_IFO+1];
  REAL4        shortSlideOffset;
  UINT4        numShortSlides;
  UINT4        fftLevel;
  UINT4        simDataType;
  REAL4        clusterWindow;
  SimInspiralTable *injectList;
  REAL4        injSearchWindow;
  REAL4        injMchirpWindow;
  /* flags */
  int          strainData;
  int          ligoDoubleData;
  int          virgoDoubleData;
  int          simData;
  int          zeroData;
  int          whiteSpectrum;
  int          getData;
  int          getSpectrum;
  int          getBank;
  int          doFilter;
  int          analyzeInjSegsOnly;
  int          doNullStream;
  int          doTraceSNR;
  int          doBankVeto;
  int          doAutoVeto;
  int          doChiSquare;
  int          doSnglChiSquared;
  int          singlePolFlag;
  int          clusterFlag;
  int          faceOnAnalysis;
  int          faceAwayAnalysis;
  int          dynTempLength;
  int          storeAmpParams;
  int          analSegmentEnd;
  int          doShortSlides;
  int          writeSnglInspiralTable;
  /* write intermediate result flags */
  int          writeRawData;
  int          writeProcessedData;
  int          writeInvSpectrum;
  int          writeSegment;
  int          writeFilterOutput;
  LIGOTimeGPS  jobStartTime;
  int          faceOnStatistic;
};

/* Other structures */

struct bankTemplateOverlaps {
  REAL8Array  *PTFM[LAL_NUM_IFO+1];
};

struct bankComplexTemplateOverlaps {
  COMPLEX8Array  *PTFM[LAL_NUM_IFO];
  INT4           timeStepDelay;
};

struct bankDataOverlaps {
  COMPLEX8VectorSequence *PTFqVec[LAL_NUM_IFO+1];
};

struct bankCohTemplateOverlaps {
  gsl_matrix *rotReOverlaps;
  gsl_matrix *rotImOverlaps;
};

typedef struct tagRingDataSegments
{
  UINT4                    numSgmnt;
  COMPLEX8FrequencySeries *sgmnt;
}
RingDataSegments;

typedef struct tagCohPTFSkyPositions
{
  UINT4       numPoints;
  SkyPosition *data;
}
CohPTFSkyPositions; 

typedef struct tagTimeSlideVectorList
{
  REAL4       timeSlideVectors[LAL_NUM_IFO];
  INT8        timeSlideID;
  UINT4        analStartPoint;
  UINT4        analEndPoint;
}
TimeSlideVectorList;

/* ENUM for sky location looping */

typedef enum
{
  SINGLE_SKY_POINT,
  TWO_DET_SKY_PATCH,
  SKY_PATCH,
  ALL_SKY,
  TWO_DET_ALL_SKY
}
Skyloopingtype;

/* Function declarations for coh_PTF_inspiral */

void coh_PTF_statistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    REAL4                   *timeOffsets,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *autoVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *chiSquare[LAL_NUM_IFO+1],
    UINT4                   subBankSize,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankTemplateOverlaps *bankNormOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    RingDataSegments        *segment[LAL_NUM_IFO+1],
    COMPLEX8FFTPlan         *invPlan,
    struct bankDataOverlaps **chisqOverlapsP,
    struct bankDataOverlaps **chisqSnglOverlapsP,
    REAL4 *frequencyRangesPlus[LAL_NUM_IFO+1],
    REAL4 *frequencyRangesCross[LAL_NUM_IFO+1],
    struct timeval          startTime,
    UINT4                   segStartPoint,
    UINT4                   segEndPoint
);

UINT8 coh_PTF_add_triggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    UINT4                   spinTrigger,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *autoVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *chiSquare[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4                   rightAscension,
    REAL4                   declination,
    INT8                    slideId,
    REAL4                   *timeOffsets,
    UINT4                   startPoint,
    UINT4                   endPoint
);
void coh_PTF_cluster_triggers(
  struct coh_PTF_params   *params,
  MultiInspiralTable      **eventList,
  MultiInspiralTable      **thisEvent
);

UINT4 coh_PTF_accept_trig_check(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      thisEvent
);

SnglInspiralTable* coh_PTF_create_sngl_event(
    struct coh_PTF_params   *params,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   *eventId,
    REAL4TimeSeries         **pValues,
    REAL4TimeSeries         **bankVeto,
    REAL4TimeSeries         **autoVeto,
    REAL4TimeSeries         **chiSquare,
    REAL8Array              **PTFM,
    UINT4                   currPos
);

UINT8 coh_PTF_add_sngl_triggers(
    struct coh_PTF_params   *params,
    SnglInspiralTable       **eventList,
    SnglInspiralTable       **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    REAL4TimeSeries         **pValues,
    REAL4TimeSeries         **bankVeto,
    REAL4TimeSeries         **autoVeto,
    REAL4TimeSeries         **chiSquare,
    REAL8Array              **PTFM,
    UINT4                   startPoint,
    UINT4                   endPoint
);

UINT4 coh_PTF_accept_sngl_trig_check(
    struct coh_PTF_params   *params,
    SnglInspiralTable      **eventList,
    SnglInspiralTable      thisEvent
);

void coh_PTF_cluster_sngl_triggers(
    struct coh_PTF_params   *params,
    SnglInspiralTable      **eventList,
    SnglInspiralTable      **thisEvent
);


/* Function declarations for coh_PTF_spin_checker */

int coh_PTF_spin_checker(
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL8Array              *PTFN[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   singleDetector,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    INT4                    segmentNumber
);

/* Function declarations for coh_PTF_utils */

INT4 coh_PTF_data_condition(
              struct coh_PTF_params *params,
              REAL4TimeSeries          **channel,
              REAL4FrequencySeries     **invspec,
              RingDataSegments         **segments,
              REAL4FFTPlan             *fwdplan,
              REAL4FFTPlan             *psdplan,
              REAL4FFTPlan             *revplan,
              REAL4                    **timeSlideVectors,
              struct timeval           startTime
);

REAL4TimeSeries *coh_PTF_get_data( 
    struct coh_PTF_params *params,
    const char *ifoChannel,
    const char *dataCache,
    UINT4 ifoNumber
);

void coh_PTF_setup_null_stream(
    struct coh_PTF_params   *params,
    REAL4TimeSeries         **channel,
    REAL4FrequencySeries    **invspec,
    RingDataSegments        **segments,
    REAL4                   *Fplustrig,
    REAL4                   *Fcrosstrig,
    REAL4                   *timeOffsets,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    REAL4FFTPlan            *psdplan,
    REAL4                   *timeSlideVectors,
    struct timeval           startTime
);

int coh_PTF_get_null_stream(
    struct coh_PTF_params *params,
    REAL4TimeSeries *channel[LAL_NUM_IFO + 1],
    REAL4 *Fplus,
    REAL4 *Fcross,
    REAL4 *timeOffsets
);

REAL4FrequencySeries *coh_PTF_get_invspec(
    REAL4TimeSeries         *channel,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    REAL4FFTPlan            *psdplan,
    struct coh_PTF_params   *params
);

void coh_PTF_rescale_data (
    REAL4TimeSeries *channel,
    REAL8 rescaleFactor
);

RingDataSegments *coh_PTF_get_segments(
    REAL4TimeSeries         *channel,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    InterferometerNumber     NumberIFO,
    REAL4                   *timeSlideVectors,
    struct coh_PTF_params   *params
);

void coh_PTF_create_time_slide_table(
  struct coh_PTF_params   *params,
  INT8                    *slideIDList,
  RingDataSegments        **segments,
  TimeSlide               **time_slide_headP,
  TimeSlideSegmentMapTable **time_slide_map_headP,
  SegmentTable            **segment_table_headP,
  TimeSlideVectorList     **longTimeSlideListP,
  TimeSlideVectorList     **shortTimeSlideListP,
  REAL4                   *timeSlideVectors,
  INT4                    numSegments
);

void coh_PTF_initialize_structures(
  struct coh_PTF_params    *params,
  FindChirpInitParams      **fcInitParamsP,
  FindChirpTemplate        **fcTmpltP,
  FindChirpTmpltParams     **fcTmpltParamsP,
  REAL8Array               **PTFM,
  REAL8Array               **PTFN,
  COMPLEX8VectorSequence   **PTFqVec,
  REAL4FFTPlan             *fwdplan
);

void coh_PTF_initialize_time_series(
  struct coh_PTF_params    *params,
  LIGOTimeGPS              segStartTime,
  REAL8                    fLower,
  REAL4TimeSeries          **cohSNRP,
  REAL4TimeSeries          **nullSNRP,
  REAL4TimeSeries          **traceSNRP,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  REAL4TimeSeries          **chiSquare,
  REAL4TimeSeries          **snrComps,
  REAL4TimeSeries          **pValues,
  REAL4TimeSeries          **gammaBeta,
  UINT4                    spinTemplates
);

void coh_PTF_reset_time_series(
  struct coh_PTF_params    *params,
  LIGOTimeGPS              segStartTime,
  REAL4TimeSeries          *cohSNR,
  REAL4TimeSeries          *nullSNR,
  REAL4TimeSeries          *traceSNR,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  REAL4TimeSeries          **chiSquare,
  REAL4TimeSeries          **snrComps,
  REAL4TimeSeries          **pValues,
  REAL4TimeSeries          **gammaBeta,
  UINT4                    spinTemplates
);

void coh_PTF_destroy_time_series(
  REAL4TimeSeries          *cohSNR,
  REAL4TimeSeries          *nullSNR,
  REAL4TimeSeries          *traceSNR,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  REAL4TimeSeries          **chiSquare,
  REAL4TimeSeries          **pValues,
  REAL4TimeSeries          **gammaBeta,
  REAL4TimeSeries          **snrComps
);

void coh_PTF_calculate_single_detector_filters(
  struct coh_PTF_params      *params,
  FindChirpTemplate          *fcTmplt,
  REAL4FrequencySeries       **invspec,
  REAL8Array                 **PTFM,
  COMPLEX8VectorSequence     **PTFqVec,
  REAL4TimeSeries            **snrComps,
  RingDataSegments           **segments,
  COMPLEX8FFTPlan            *invPlan,
  UINT4                      spinTemplate,
  UINT4                      segNum
);

void coh_PTF_calculate_single_det_spin_snr(
  struct coh_PTF_params      *params,
  REAL8Array                 **PTFM,
  COMPLEX8VectorSequence     **PTFqVec,
  REAL4TimeSeries            **snrComps,
  UINT4                      ifoNumber
);

REAL4 coh_PTF_get_spin_SNR(
  REAL4 *v1p,
  REAL4 *v2p,
  UINT4 vecLengthTwo
);

void coh_PTF_calculate_coherent_SNR(
  struct coh_PTF_params      *params,
  REAL4                      *snrData,
  REAL4TimeSeries            **pValues,
  REAL4TimeSeries            **snrComps,
  INT4                       *timeOffsetPoints,
  COMPLEX8VectorSequence     **PTFqVec,
  REAL4                      *Fplus,
  REAL4                      *Fcross,
  gsl_matrix                 *eigenvecs,
  gsl_vector                 *eigenvals,
  UINT4                      segStartPoint,
  UINT4                      segEndPoint,
  UINT4                      vecLength,
  UINT4                      vecLengthTwo,
  UINT4                      spinTemplate
);

UINT4 coh_PTF_template_time_series_cluster(
  REAL4TimeSeries *cohSNR,
  UINT4 *acceptPoints,
  INT4 numPointCheck,
  UINT4 startPoint,
  UINT4 endPoint
);

UINT4 coh_PTF_test_veto_vals(
  struct coh_PTF_params    *params,
  REAL4TimeSeries          *cohSNR,
  REAL4TimeSeries          *nullSNR,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  UINT4                    currPointLoc
);

void coh_PTF_calculate_null_stream_filters(
  struct coh_PTF_params      *params,
  FindChirpTemplate          *fcTmplt,
  REAL4FrequencySeries       **invspec,
  REAL8Array                 **PTFM,
  COMPLEX8VectorSequence     **PTFqVec,
  RingDataSegments           **segments,
  COMPLEX8FFTPlan            *invPlan,
  UINT4                      spinTemplate,
  UINT4                      segNum
);

void coh_PTF_calculate_null_stream_norms(
  UINT4 vecLength,
  gsl_matrix *eigenvecsNull,
  gsl_vector *eigenvalsNull,
  REAL8Array *PTFM[LAL_NUM_IFO+1]
);

void coh_PTF_calculate_null_stream_snr(
  struct coh_PTF_params   *params,
  REAL4TimeSeries         *nullSNR,
  COMPLEX8VectorSequence  **PTFqVec,
  gsl_matrix              *eigenvecsNull,
  gsl_vector              *eigenvalsNull,
  UINT4                   spinTemplate,
  UINT4                   vecLength,
  UINT4                   vecLoc,
  UINT4                   snrLoc
);

void coh_PTF_calculate_trace_snr(
  struct coh_PTF_params   *params,
  REAL4TimeSeries         *traceSNR,
  COMPLEX8VectorSequence  **PTFqVec,
  gsl_matrix              *eigenvecs,
  gsl_vector              *eigenvals,
  REAL4                   *Fplus,
  REAL4                   *Fcross,
  INT4                    *timeOffsetPoints,
  UINT4                   spinTemplate,
  UINT4                   vecLength,
  UINT4                   vecLengthTwo,
  UINT4                   vecLoc,
  UINT4                   snrLoc
);

UINT4 coh_PTF_initialize_bank_veto(
  struct coh_PTF_params   *params,
  struct bankTemplateOverlaps **bankNormOverlapsP,
  struct bankComplexTemplateOverlaps **bankOverlapsP,
  struct bankDataOverlaps **dataOverlapsP,
  FindChirpTemplate        **bankFcTmpltsP,
  FindChirpTemplate        *fcTmplt,
  FindChirpTmpltParams     *fcTmpltParams,
  REAL4FrequencySeries     **invspec,
  struct timeval           startTime
);

void coh_PTF_bank_veto_segment_setup(
  struct coh_PTF_params   *params,
  struct bankDataOverlaps *dataOverlaps,
  FindChirpTemplate       *bankFcTmplts,
  RingDataSegments        **segments,
  COMPLEX8VectorSequence  **PTFqVec,
  COMPLEX8FFTPlan         *invplan,
  INT4                    segmentNum,
  struct timeval           startTime
);

void coh_PTF_bank_veto_coh_setup(
  struct coh_PTF_params   *params,
  gsl_matrix              **Bankeigenvecs,
  gsl_vector              **Bankeigenvals,
  struct bankCohTemplateOverlaps **bankCohOverlapsP,
  struct bankComplexTemplateOverlaps *bankOverlaps,
  REAL4                   *Fplus,
  REAL4                   *Fcross,
  REAL8Array              **PTFM,
  struct bankTemplateOverlaps *bankNormOverlaps,
  UINT4                   csVecLength,
  UINT4                   csVecLengthTwo,
  UINT4                   vecLength
);

UINT4 coh_PTF_initialize_auto_veto(
  struct coh_PTF_params   *params,
  struct bankComplexTemplateOverlaps **autoTempOverlapsP,
  struct timeval           startTime
);

void coh_PTF_calculate_bank_veto_template_filters(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *bankFcTmplts,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    **invspec,
    struct bankComplexTemplateOverlaps *bankOverlaps
);

void coh_PTF_calculate_auto_veto_template_filters(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    REAL4FrequencySeries    **invspec,
    COMPLEX8FFTPlan         *invplan,
    UINT4                   timeStepPoints
);

void coh_PTF_auto_veto_coh_setup(
  struct coh_PTF_params   *params,
  gsl_matrix              **AutoeigenvecsP,
  gsl_vector              **AutoeigenvalsP,
  struct bankCohTemplateOverlaps **autoCohOverlapsP,
  struct bankComplexTemplateOverlaps *autoTempOverlaps,
  REAL4                   *Fplus,
  REAL4                   *Fcross,
  REAL8Array              **PTFM,
  UINT4                   csVecLength,
  UINT4                   csVecLengthTwo,
  UINT4                   vecLength
);

void coh_PTF_chi_square_coh_setup(
  struct coh_PTF_params   *params,
  gsl_matrix              **AutoeigenvecsP,
  gsl_vector              **AutoeigenvalsP,
  REAL4                   **frequencyRangesPlus,
  REAL4                   **frequencyRangesCross,
  REAL4                   **powerBinsPlus,
  REAL4                   **powerBinsCross,
  struct bankDataOverlaps **chisqOverlapsP,
  FindChirpTemplate       *fcTmplt,
  REAL4FrequencySeries    **invspec,
  RingDataSegments        **segments,
  REAL4                   *Fplus,
  REAL4                   *Fcross,
  REAL8Array              **PTFM,
  COMPLEX8FFTPlan         *invPlan,
  INT4                    segmentNumber,
  UINT4                   csVecLength,
  UINT4                   csVecLengthTwo,
  UINT4                   vecLength
);

void coh_PTF_chi_square_sngl_setup(
  struct coh_PTF_params   *params,
  REAL4                   **frequencyRangesPlus,
  REAL4                   **frequencyRangesCross,
  REAL4                   **powerBinsPlus,
  REAL4                   **powerBinsCross,
  struct bankDataOverlaps **chisqSnglOverlapsP,
  FindChirpTemplate       *fcTmplt,
  REAL4FrequencySeries    **invspec,
  RingDataSegments        **segments,
  REAL8Array              **PTFM,
  COMPLEX8FFTPlan         *invPlan,
  INT4                    segmentNumber
);

void coh_PTF_calculate_det_stuff(
  struct coh_PTF_params   *params,
  LALDetector             **detectors,
  REAL4                   *timeOffsets,
  REAL4                   *Fplustrig,
  REAL4                   *Fcrosstrig,
  CohPTFSkyPositions      *skyPoints,
  UINT4                   skyPointNum
);

void coh_PTF_convert_time_offsets_to_points(
  struct coh_PTF_params   *params,
  REAL4                   *timeOffsets,
  INT4                    *timeOffsetPoints
);

void coh_PTF_calculate_bmatrix(
  struct coh_PTF_params   *params,
  gsl_matrix *eigenvecs,
  gsl_vector *eigenvals,
  REAL4 Fplus[LAL_NUM_IFO],
  REAL4 Fpcross[LAL_NUM_IFO],
  REAL8Array              *PTFM[LAL_NUM_IFO+1],
  UINT4 vecLength,
  UINT4 vecLengthTwo,
  UINT4 PTFMlen
);

void coh_PTF_calculate_rotated_vectors(
    struct coh_PTF_params   *params,
    COMPLEX8VectorSequence  **PTFqVec,
    REAL4 *u1,
    REAL4 *u2,
    REAL4 *Fplus,
    REAL4 *Fcross,
    INT4  *timeOffsetPoints,
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    UINT4 numPoints,
    UINT4 position,
    UINT4 vecLength,
    UINT4 vecLengthTwo,
    UINT4 detectorNum
);

MultiInspiralTable* coh_PTF_create_multi_event(
    struct coh_PTF_params   *params,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   *eventId,
    UINT4                   spinTrigger,
    REAL4TimeSeries         **pValues,
    REAL4TimeSeries         **gammaBeta,
    REAL4TimeSeries         **snrComps,
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         **bankVeto,
    REAL4TimeSeries         **autoVeto,
    REAL4TimeSeries         **chiSquare,
    REAL8Array              **PTFM,
    REAL4                   rightAscension,
    REAL4                   declination,
    INT8                    slideId,
    INT4                   *timeOffsetPoints,
    UINT4                   currPos
);

void coh_PTF_cleanup(
    struct coh_PTF_params   *params,
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *psdplan,
    REAL4FFTPlan            *revplan,
    COMPLEX8FFTPlan         *invPlan,
    REAL4TimeSeries         **channel,
    REAL4FrequencySeries    **invspec,
    RingDataSegments        **segments,
    MultiInspiralTable      *events,
    SnglInspiralTable       *snglEvents,
    InspiralTemplate        *PTFbankhead,
    FindChirpTemplate       *fcTmplt,
    FindChirpTmpltParams    *fcTmpltParams,
    FindChirpInitParams     *fcInitParams,
    REAL8Array              **PTFM,
    REAL8Array              **PTFN,
    COMPLEX8VectorSequence  **PTFqVec,
    REAL4                   *timeOffsets,
    REAL4                   *slidTimeOffsets,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    REAL4                   *Fplustrig,
    REAL4                   *Fcrosstrig,
    CohPTFSkyPositions      *skyPoints,
    TimeSlide               *time_slide_head,
    TimeSlideVectorList     *longTimeSlideList,
    TimeSlideVectorList     *shortTimeSlideList,
    REAL4                   *timeSlideVectors,
    LALDetector             **detectors,
    INT8                    *slideIDList,
    TimeSlideSegmentMapTable *time_slide_map_head,
    SegmentTable            *segment_table_head
);

REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params );

REAL4FFTPlan *coh_PTF_get_fft_psdplan( struct coh_PTF_params *params );

REAL4FFTPlan *coh_PTF_get_fft_revplan( struct coh_PTF_params *params );

COMPLEX8FFTPlan *coh_PTF_get_fft_invplan( struct coh_PTF_params *params );

int is_in_list( int i, const char *list );

SnglInspiralTable *conv_insp_tmpl_to_sngl_table(
    InspiralTemplate        *template,
    UINT4                   eventNumber
);

long int timeval_subtract(struct timeval *t1);

void timeval_print(struct timeval *tv);

/* Function declarations for coh_PTF_template.c */

void coh_PTF_template (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
);

void coh_PTF_template_PTF (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
);

void coh_PTF_normalize(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt,
    REAL4FrequencySeries       *invspec,
    REAL8Array                 *PTFM,
    REAL8Array                 *PTFN,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan,
    UINT4                      spinTemplate
);

void coh_PTF_template_overlaps(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt1,
    FindChirpTemplate          *fcTmplt2,
    REAL4FrequencySeries       *invspec,
    UINT4                      spinBank,
    REAL8Array                 *PTFM
);

void coh_PTF_complex_template_overlaps(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt1,
    FindChirpTemplate          *fcTmplt2,
    REAL4FrequencySeries       *invspec,
    UINT4                      spinBank,
    COMPLEX8Array                 *PTFM
);

void coh_PTF_bank_filters(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt,
    UINT4                      spinBank,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invBankPlan,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8VectorSequence     *PTFBankqVec,
    REAL8                      f_min,
    REAL8                      fFinal
);

void coh_PTF_auto_veto_overlaps(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    REAL4FrequencySeries       *invspec,
    COMPLEX8FFTPlan            *invBankPlan,
    UINT4                      spinBank,
    UINT4                      numAutoPoints,
    UINT4                      timeStepPoints,
    UINT4                      ifoNumber
);

/* Function declarations in coh_PTF_chi_squared.c*/

UINT4 coh_PTF_read_sub_bank(
    struct coh_PTF_params   *params,
    InspiralTemplate        **PTFBankTemplates
);

void coh_PTF_initialise_sub_bank(
    struct coh_PTF_params   *params,
    InspiralTemplate        *PTFBankTemplates,
    FindChirpTemplate       *bankFcTmplts,
    UINT4                    subBankSize,
    UINT4                    numPoints
);

REAL4 coh_PTF_calculate_bank_veto(
    UINT4           numPoints,
    UINT4           position,
    UINT4           subBankSize,
    REAL4           Fplus[LAL_NUM_IFO],
    REAL4           Fcross[LAL_NUM_IFO],
    struct coh_PTF_params      *params,
    struct bankCohTemplateOverlaps *cohBankOverlaps,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    struct bankTemplateOverlaps *bankNormOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL8Array      *PTFM[LAL_NUM_IFO+1],
    INT4            timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix **Bankeigenvecs,
    gsl_vector **Bankeigenvals,
    UINT4       detectorNum,
    UINT4       vecLength,
    UINT4       vecLengthTwo
);

REAL4 coh_PTF_calculate_auto_veto(
    UINT4           numPoints,
    UINT4           position,
    REAL4           Fplus[LAL_NUM_IFO],
    REAL4           Fcross[LAL_NUM_IFO],
    struct coh_PTF_params      *params,
    struct bankCohTemplateOverlaps *cohAutoOverlaps,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL8Array      *PTFM[LAL_NUM_IFO+1],
    INT4            timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *Autoeigenvecs,
    gsl_vector *Autoeigenvals,
    UINT4       detectorNum,
    UINT4       vecLength,
    UINT4       vecLengthTwo
);

void coh_PTF_free_veto_memory(
  struct coh_PTF_params   *params,
  struct bankTemplateOverlaps *bankNormOverlaps,
  FindChirpTemplate       *bankFcTmplts,
  struct bankComplexTemplateOverlaps *bankOverlaps,
  struct bankDataOverlaps *dataOverlaps,
  struct bankComplexTemplateOverlaps *autoTempOverlaps
);

void coh_PTF_calculate_coherent_bank_overlaps(
    struct coh_PTF_params   *params,
    struct bankComplexTemplateOverlaps bankOverlaps,
    struct bankCohTemplateOverlaps cohBankOverlaps,
    REAL4           Fplus[LAL_NUM_IFO],
    REAL4           Fcross[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    gsl_matrix *Bankeigenvecs,
    gsl_vector *Bankeigenvals,
    UINT4 vecLength,
    UINT4 vecLengthTwo
);

void coh_PTF_calculate_standard_chisq_freq_ranges(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4 Fplus[LAL_NUM_IFO],
    REAL4 Fcross[LAL_NUM_IFO],
    REAL4 *frequencyRangesPlus,
    REAL4 *frequencyRangesCross,
    gsl_matrix *eigenvecs,
    UINT4 detectorNum,
    UINT4 singlePolFlag
);

void coh_PTF_calculate_standard_chisq_power_bins(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4 Fplus[LAL_NUM_IFO],
    REAL4 Fcross[LAL_NUM_IFO],
    REAL4 *frequencyRangesPlus,
    REAL4 *frequencyRangesCross,
    REAL4 *powerBinsPlus,
    REAL4 *powerBinsCross,
    gsl_matrix *eigenvecs,
    UINT4 detectorNum,
    UINT4 singlePolFlag
);

REAL4 coh_PTF_calculate_chi_square(
    struct coh_PTF_params   *params,
    UINT4           position,
    struct bankDataOverlaps *chisqOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL8Array      *PTFM[LAL_NUM_IFO+1],
    REAL4           Fplus[LAL_NUM_IFO],
    REAL4           Fcross[LAL_NUM_IFO],
    INT4            timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    REAL4 *powerBinsPlus,
    REAL4 *powerBinsCross,
    UINT4 detectorNum,
    UINT4 vecLength,
    UINT4 vecLengthTwo
);

/* routines in coh_PTF_option */

int coh_PTF_parse_options(
    struct coh_PTF_params *params,
    int argc,
    char **argv
);

int coh_PTF_default_params( struct coh_PTF_params *params );

int coh_PTF_params_sanity_check( struct coh_PTF_params *params );

int coh_PTF_params_inspiral_sanity_check( struct coh_PTF_params *params );

int coh_PTF_params_spin_checker_sanity_check( struct coh_PTF_params *params );

int coh_PTF_usage( const char *program );

/* routines in coh_PTF_output */
/* There is some duplication with the ring_output */

ProcessParamsTable * create_process_params(
    int argc,
    char **argv,
    const char *program
);

int coh_PTF_output_events_xml(
    char               *outputFile,
    MultiInspiralTable  *events,
    SnglInspiralTable *snglEvents,
    SimInspiralTable *injections,
    ProcessParamsTable *processParamsTable,
    TimeSlide          *time_slide_head,
    TimeSlideSegmentMapTable *time_slide_map_head,
    SegmentTable       *segment_table_head,
    struct coh_PTF_params *params
    );

int coh_PTF_output_tmpltbank(
    char               *outputFile,
    SnglInspiralTable   *tmplts,
    ProcessParamsTable *processParamsTable,
    struct coh_PTF_params *params
);

ProcessTable * coh_PTF_create_process_table(
    struct coh_PTF_params *params
);   

SearchSummaryTable *coh_PTF_create_search_summary(
    struct coh_PTF_params *params
);

int write_REAL4TimeSeries( REAL4TimeSeries *series );

int write_REAL4FrequencySeries( REAL4FrequencySeries *series );

int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series );

int generate_file_name(
    char *fname,
    size_t size,
    const char *sname,
    int t,
    int dt
);

/* generate array of sky points based on RA, Dec, error radius, and timing
 * accuracy */

CohPTFSkyPositions *coh_PTF_generate_sky_points(
    struct coh_PTF_params *params
);

CohPTFSkyPositions *coh_PTF_generate_sky_grid(
    struct coh_PTF_params *params
);

CohPTFSkyPositions *coh_PTF_circular_grid(
    REAL4 angularResolution,
    REAL4 skyError
);

CohPTFSkyPositions *coh_PTF_parse_time_delays(
    CohPTFSkyPositions *skyPoints,
    struct coh_PTF_params *params
);

CohPTFSkyPositions *coh_PTF_read_grid_from_file(
    const char *fname,
    UINT4      raColumn,
    UINT4      decColumn
);

void coh_PTF_rotate_skyPoints(
    CohPTFSkyPositions *skyPoints,
    gsl_vector *axis,
    REAL8 angle
);

void coh_PTF_rotate_SkyPosition(
    SkyPosition *skyPoint,
    gsl_matrix  *matrix
);

CohPTFSkyPositions *coh_PTF_two_det_sky_grid(
    struct coh_PTF_params *params
);

CohPTFSkyPositions *coh_PTF_three_det_sky_grid(
    struct coh_PTF_params *params
);

void normalise(
    gsl_vector *vec
);

void cross_product(
    gsl_vector *product,
    const gsl_vector *u,
    const gsl_vector *v
);

void rotation_matrix(
    gsl_matrix *matrix,
    gsl_vector *axis,
    REAL8 angle
);

void REALToGSLVector(
    const REAL8 *input,
    gsl_vector  *output,
    size_t      size
);

void findInjectionSegment(
    UINT4 *start,
    UINT4 *end,
    LIGOTimeGPS *epoch,
    struct coh_PTF_params *params
);

UINT4 coh_PTF_trig_time_check(
    struct coh_PTF_params *params,
    LIGOTimeGPS segStartTime,
    LIGOTimeGPS segEndTime);

UINT4 checkInjectionMchirp(
    struct coh_PTF_params *params,
    InspiralTemplate *tmplt,
    LIGOTimeGPS *epoch
);

void coh_PTF_set_null_input_REAL4TimeSeries(
  REAL4TimeSeries** timeSeries,
  UINT4 length
);

void coh_PTF_set_null_input_REAL4FrequencySeries(
  REAL4FrequencySeries** freqSeries,
  UINT4 length
);

void coh_PTF_set_null_input_RingDataSegments(
  RingDataSegments** segment,
  UINT4 length
);

void coh_PTF_set_null_input_REAL8Array(
  REAL8Array** array,
  UINT4 length
);

void coh_PTF_set_null_input_COMPLEX8VectorSequence(
  COMPLEX8VectorSequence** vecSeq,
  UINT4 length
);

void coh_PTF_set_null_input_REAL4(
  REAL4** array,
  UINT4 length
);

void coh_PTF_set_null_input_LALDetector(
  LALDetector** detector,
  UINT4 length
);
