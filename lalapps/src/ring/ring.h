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
  INT4         injectType;
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
  char         userTag[256];
  char         ifoTag[256];
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
SnglRingdownTable * ring_filter(
    RingDataSegments         *segments,
    RingTemplateBank         *bank,
    REAL4FrequencySeries     *invSpectrum,
    REAL4FFTPlan             *fwdPlan,
    REAL4FFTPlan             *revPlan,
    struct ring_params       *params
    );
