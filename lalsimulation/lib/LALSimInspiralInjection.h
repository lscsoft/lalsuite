#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LALDict.h>

struct tagLALSimInspiralInjectionSequence;
typedef struct tagLALSimInspiralInjectionSequence LALSimInspiralInjectionSequence;

void XLALSimInspiralDestroyInjectionSequence(LALSimInspiralInjectionSequence *sequence);
LALSimInspiralInjectionSequence * XLALSimInspiralCreateInjectionSequence(size_t length);
LALSimInspiralInjectionSequence * XLALSimInspiralCutInjectionSequence(const LALSimInspiralInjectionSequence *sequence, size_t first, size_t length);
LALSimInspiralInjectionSequence * XLALSimInspiralCopyInjectionSequence(const LALSimInspiralInjectionSequence *sequence);
void XLALSimInspiralShiftInjectionSequence(LALSimInspiralInjectionSequence *sequence, int count);
LALSimInspiralInjectionSequence * XLALSimInspiralResizeInjectionSequence(LALSimInspiralInjectionSequence *sequence, int first, size_t length);

size_t LALSimInspiralInjectionSequenceLength(LALSimInspiralInjectionSequence *sequence);
LALDict * LALSimInspiralInjectionSequenceGet(LALSimInspiralInjectionSequence *sequence, int pos);
int LALSimInspiralInjectionSequenceSet(LALSimInspiralInjectionSequence *sequence, LALDict *inparmas, int pos);

LALSimInspiralInjectionSequence * XLALSimInspiralInjectionSequenceFromH5File(const char *fname);
int XLALSimInspiralInjectionSequenceToH5File(const LALSimInspiralInjectionSequence *sequence, const char *fname);

LIGOTimeGPS * XLALSimInspiralInjectionEndTime(LIGOTimeGPS *epoch, LALDict *injparams);
LIGOTimeGPS * XLALSimInspiralInjectionStartTime(LIGOTimeGPS *epoch, LALDict *injparams);

int XLALSimInspiralInjectionSequenceIsEndTimeOrdered(LALSimInspiralInjectionSequence *sequence);
int XLALSimInspiralInjectionSequenceIsStartTimeOrdered(LALSimInspiralInjectionSequence *sequence);
int XLALSimInspiralInjectionSequenceOrderByEndTime(LALSimInspiralInjectionSequence *sequence);
int XLALSimInspiralInjectionSequenceOrderByStartTime(LALSimInspiralInjectionSequence *sequence);
LALSimInspiralInjectionSequence * XLALSimInspiralInjectionSequenceInInterval(const LALSimInspiralInjectionSequence *sequence, const LIGOTimeGPS *start, const LIGOTimeGPS *end);

int XLALSimInspiralInjectionTDWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *injparams, REAL8 deltaT);

REAL8TimeSeries * XLALSimInspiralInjectionStrain(LALDict *injparams, REAL8 deltaT, const LALDetector *detector);
