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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>


enum { write_frame, write_ascii };

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
  REAL4        lowCutoffFrequency;
  REAL4        highpassFrequency;
  REAL4        invSpecLen;
  REAL4        threshold;
  REAL4        timeWindow;
  REAL4        spinSNR2threshold;
  REAL4        nonspinSNR2threshold;
  REAL4        rightAscension;
  REAL4        declination;
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

/* routines in ring_option */
int coh_PTF_parse_options(struct coh_PTF_params *params,int argc,char **argv );
int coh_PTF_params_sanity_check( struct coh_PTF_params *params );
int coh_PTF_params_inspiral_sanity_check( struct coh_PTF_params *params );
int coh_PTF_params_spin_checker_sanity_check( struct coh_PTF_params *params );

/* routines in ring_output */
ProcessParamsTable * create_process_params( int argc, char **argv,
    const char *program );

/* routines to write intermediate results in ring_output */
int write_REAL4TimeSeries( REAL4TimeSeries *series );
int write_REAL4FrequencySeries( REAL4FrequencySeries *series );
int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series );

int cohPTF_output_events_xml(
    char               *outputFile,
    MultiInspiralTable *events,
    ProcessParamsTable *processParamsTable,
    struct coh_PTF_params *params
    );

int cohPTF_output_tmpltbank(
    char               *outputFile,
    SnglInspiralTable   *tmplts,
    ProcessParamsTable *processParamsTable,
    struct coh_PTF_params *params
    );

UINT4 read_sub_bank(
struct coh_PTF_params   *params,
InspiralTemplate        **PTFBankTemplates);

void initialise_sub_bank(
struct coh_PTF_params   *params,
InspiralTemplate        *PTFBankTemplates,
FindChirpTemplate       *bankFcTmplts,
UINT4                    subBankSize,
UINT4                    numPoints,
UINT4                    spinBank);

void
cohPTFTemplate (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    );

void
cohPTFNormalize(
    FindChirpTemplate          *fcTmplt,
    REAL4FrequencySeries       *invspec,
    REAL8Array                 *PTFM,
    REAL8Array                 *PTFN,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan,
    UINT4                      spinTemplate
    );

void cohPTFTemplateOverlaps(
    FindChirpTemplate          *fcTmplt1,
    FindChirpTemplate          *fcTmplt2,
    REAL4FrequencySeries       *invspec,
    UINT4                      spinBank,
    REAL8Array                 *PTFM);

void
cohPTFComplexTemplateOverlaps(
    FindChirpTemplate          *fcTmplt1,
    FindChirpTemplate          *fcTmplt2,
    REAL4FrequencySeries       *invspec,
    UINT4                      spinBank,
    COMPLEX8Array                 *PTFM
    );

void cohPTFBankFilters(
    FindChirpTemplate          *fcTmplt,
    UINT4                      spinBank,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invBankPlan,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8VectorSequence     *PTFBankqVec,
    REAL8                      f_min,
    REAL8                      fFinal);

REAL4 cohPTFDataNormalize(
    COMPLEX8FrequencySeries    *sgmnt,
    REAL4FrequencySeries       *invspec);

void autoVetoOverlaps(
    FindChirpTemplate          *fcTmplt,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    REAL4FrequencySeries       *invspec,
    COMPLEX8FFTPlan            *invBankPlan,
    UINT4                      spinBank,
    UINT4                      numAutoPoints,
    UINT4                      timeStepPoints,
    UINT4                      ifoNumber );

REAL4 calculate_bank_veto(
UINT4           numPoints,
UINT4           position,
UINT4           subBankSize,
UINT4           vecLength,
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
REAL4           SNR,
REAL8Array      *PTFM[LAL_NUM_IFO+1],
struct coh_PTF_params      *params,
struct bankComplexTemplateOverlaps *bankOverlaps,
struct bankTemplateOverlaps *bankNormOverlaps,
struct bankDataOverlaps *dataOverlaps,
REAL4TimeSeries         *pValues[10],
REAL4TimeSeries         *gammaBeta[2],
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO],
UINT4           singleDetector );

REAL4 calculate_bank_veto_max_phase(
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

REAL4 calculate_bank_veto_max_phase_coherent(
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
gsl_vector *Bankeigenvals[50]);

REAL4 calculate_auto_veto_max_phase_coherent(
UINT4           numPoints,
UINT4           position,
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
struct coh_PTF_params      *params,
struct bankCohTemplateOverlaps *cohAutoOverlaps,
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO],
gsl_matrix *Autoeigenvecs,
gsl_vector *Autoeigenvals);

void calculate_coherent_bank_overlaps(
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


void free_bank_veto_memory(
  struct bankTemplateOverlaps *bankNormOverlaps,
  InspiralTemplate        *PTFBankTemplates,
  FindChirpTemplate       *bankFcTmplts,
  UINT4 subBankSize,
  struct bankComplexTemplateOverlaps *bankOverlaps,
  struct bankDataOverlaps *dataOverlaps);

void calculate_bmatrix(
  struct coh_PTF_params   *params,
  gsl_matrix *eigenvecs,
  gsl_vector *eigenvals,
  REAL4 a[LAL_NUM_IFO],
  REAL4 b[LAL_NUM_IFO],
  REAL8Array              *PTFM[LAL_NUM_IFO+1],
  UINT4 vecLength,
  UINT4 vecLengthTwo,
  UINT4 PTFMlen);

void calculate_rotated_vectors(
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
    UINT4 vecLengthTwo);

void calculate_standard_chisq_freq_ranges(
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

void calculate_standard_chisq_power_bins(
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
REAL4 calculate_chi_square(
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

