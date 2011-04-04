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

#define BUFFER_SIZE 256
#define FILENAME_SIZE 256

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
  char         ifoName[3];
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
  REAL8        strideDuration;
  REAL8        truncateDuration;
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
  REAL4        timeWindow;
  REAL4        spinSNR2threshold;
  REAL4        nonspinSNR2threshold;
  REAL4        rightAscension;
  REAL4        declination;
  REAL4        skyError;
  REAL4        timingAccuracy;
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
  const char  *spinBank;
  const char  *noSpinBank;
  char         userTag[256];
  char         ifoTag[256];
  /* flags */
  int          strainData;
  int          doubleData;
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
  /* write intermediate result flags */
  int          writeRawData;
  int          writeProcessedData;
  int          writeInvSpectrum;
  int          writeSegment;
  int          writeFilterOutput;
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

struct coh_PTF_skyPoints {
  UINT4 numPoints;
  REAL8 *rightAscension;
  REAL8 *declination;
};

/* ENUM for sky location looping */

typedef enum
{
  SINGLE_SKY_POINT,
  TWO_DET_SKY_POINT_ERROR,
  SKY_POINT_ERROR,
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
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    REAL8                   *Fplustrig,
    REAL8                   *Fcrosstrig,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto,
    REAL4TimeSeries         *autoVeto,
    REAL4TimeSeries         *chiSquare,
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
    REAL4 **frequencyRangesPlusP,
    REAL4 **frequencyRangesCrossP,
    struct timeval          startTime
);

UINT8 coh_PTF_add_triggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    UINT4                   spinTrigger,
    UINT4                   singleDetector,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto,
    REAL4TimeSeries         *autoVeto,
    REAL4TimeSeries         *chiSquare,
    REAL8Array              *PTFM[LAL_NUM_IFO+1]
);
void coh_PTF_cluster_triggers(
  MultiInspiralTable      **eventList,
  MultiInspiralTable      **thisEvent
);

/* Function declarations for coh_PTF_spin_checker */

int coh_PTF_spin_checker(
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL8Array              *PTFN[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   singleDetector,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber
);

/* Function declarations for coh_PTF_utils */

REAL4TimeSeries *coh_PTF_get_data( 
    struct coh_PTF_params *params,
    const char *ifoChannel,
    const char *dataCache,
    UINT4 ifoNumber
);

int coh_PTF_get_null_stream(
    struct coh_PTF_params *params,
    REAL4TimeSeries *channel[LAL_NUM_IFO + 1],
    REAL8 *Fplus,
    REAL8 *Fcross,
    REAL8 *timeOffsets
);

REAL4FrequencySeries *coh_PTF_get_invspec(
    REAL4TimeSeries         *channel,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
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
    struct coh_PTF_params      *params
);

void coh_PTF_calculate_bmatrix(
  struct coh_PTF_params   *params,
  gsl_matrix *eigenvecs,
  gsl_vector *eigenvals,
  REAL4 a[LAL_NUM_IFO],
  REAL4 b[LAL_NUM_IFO],
  REAL8Array              *PTFM[LAL_NUM_IFO+1],
  UINT4 vecLength,
  UINT4 vecLengthTwo,
  UINT4 PTFMlen
);

void coh_PTF_calculate_rotated_vectors(
    struct coh_PTF_params   *params,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL4 *u1,
    REAL4 *u2,
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    INT4  timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    UINT4 numPoints,
    UINT4 position,
    UINT4 vecLength,
    UINT4 vecLengthTwo
);

void coh_PTF_cleanup(
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    COMPLEX8FFTPlan         *invPlan,
    REAL4TimeSeries         *channel[LAL_NUM_IFO+1],
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    RingDataSegments        *segments[LAL_NUM_IFO+1],
    MultiInspiralTable      *events,
    InspiralTemplate        *PTFbankhead,
    FindChirpTemplate       *fcTmplt,
    FindChirpTmpltParams    *fcTmpltParams,
    FindChirpInitParams     *fcInitParams,
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL8Array              *PTFN[LAL_NUM_IFO+1],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    REAL8                   *Fplustrig,
    REAL8                   *Fcrosstrig
);

REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params );

REAL4FFTPlan *coh_PTF_get_fft_revplan( struct coh_PTF_params *params );

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

void coh_PTF_template_TaylorF2 (
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

REAL4 coh_PTF_calculate_bank_veto_max_phase(
    UINT4           numPoints,
    UINT4           position,
    UINT4           subBankSize,
    REAL8Array      *PTFM[LAL_NUM_IFO+1],
    struct coh_PTF_params      *params,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankTemplateOverlaps *bankNormOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    INT4            timeOffsetPoints[LAL_NUM_IFO]
);
    
REAL4 coh_PTF_calculate_bank_veto(
    UINT4           numPoints,
    UINT4           position,
    UINT4           subBankSize,
    REAL4           a[LAL_NUM_IFO],
    REAL4           b[LAL_NUM_IFO],
    struct coh_PTF_params      *params,
    struct bankCohTemplateOverlaps *cohBankOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    INT4            timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *Bankeigenvecs[50],
    gsl_vector *Bankeigenvals[50]
);

REAL4 coh_PTF_calculate_auto_veto(
    UINT4           numPoints,
    UINT4           position,
    REAL4           a[LAL_NUM_IFO],
    REAL4           b[LAL_NUM_IFO],
    struct coh_PTF_params      *params,
    struct bankCohTemplateOverlaps *cohAutoOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    INT4            timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *Autoeigenvecs,
    gsl_vector *Autoeigenvals
);

void coh_PTF_free_bank_veto_memory(
    struct bankTemplateOverlaps *bankNormOverlaps,
    InspiralTemplate        *PTFBankTemplates,
    FindChirpTemplate       *bankFcTmplts,
    UINT4 subBankSize,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankDataOverlaps *dataOverlaps
);

void coh_PTF_calculate_coherent_bank_overlaps(
    struct coh_PTF_params   *params,
    struct bankComplexTemplateOverlaps bankOverlaps,
    struct bankCohTemplateOverlaps cohBankOverlaps,
    REAL4           a[LAL_NUM_IFO],
    REAL4           b[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    gsl_matrix *Bankeigenvecs,
    gsl_vector *Bankeigenvals
);

void coh_PTF_calculate_standard_chisq_freq_ranges(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    REAL4 *frequencyRangesPlus,
    REAL4 *frequencyRangesCross,
    gsl_matrix *eigenvecs
);

void coh_PTF_calculate_standard_chisq_power_bins(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    REAL4 *frequencyRanges,
    REAL4 *powerBinsPlus,
    REAL4 *powerBinsCross,
    gsl_matrix *eigenvecs
);

REAL4 coh_PTF_calculate_chi_square(
    struct coh_PTF_params   *params,
    UINT4           numPoints,
    UINT4           position,
    struct bankDataOverlaps *chisqOverlaps,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL4           a[LAL_NUM_IFO],
    REAL4           b[LAL_NUM_IFO],
    INT4            timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    REAL4 *powerBinsPlus,
    REAL4 *powerBinsCross
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
    MultiInspiralTable *events,
    ProcessParamsTable *processParamsTable,
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

/* generate array of sky points based on RA, DEC and DECERROR */
/* should return linked list */
void coh_PTF_generate_sky_points(
    struct coh_PTF_skyPoints *skyPoints,
    struct coh_PTF_params *params
);

void coh_PTF_sky_grid(
    struct coh_PTF_skyPoints *skyPoints,
    struct coh_PTF_params *params
);

struct coh_PTF_skyPoints coh_PTF_circular_grid(
    REAL4 angularResolution,
    REAL4 skyError
);
