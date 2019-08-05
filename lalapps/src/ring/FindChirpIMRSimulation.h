#include <lal/FindChirp.h>

void
LALFindChirpInjectIMR (
    LALStatus                     *status,
    REAL4TimeSeries               *chan,
    SimInspiralTable              *events,
    SimRingdownTable              *ringdownevents,
    COMPLEX8FrequencySeries       *resp,
    INT4                           injectSignalType
    );
