#ifndef _LAL_SIM_INSPIRAL_GENERATOR_PRIVATE_H
#define _LAL_SIM_INSPIRAL_GENERATOR_PRIVATE_H

#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LALSimSphHarmSeries.h>
#include "LALSimInspiral.h"

struct tagLALSimInspiralGenerator {

    const char *name;

    int (*initialize) (LALSimInspiralGenerator * myself, LALDict *params);

    int (*finalize) (LALSimInspiralGenerator * myself);

    int (*generate_td_modes) (
        SphHarmTimeSeries **hlm,
        LALDict *params,
        LALSimInspiralGenerator *myself
    );

    int (*generate_td_waveform) (
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        LALDict *params,
        LALSimInspiralGenerator *myself
    );

    int (*generate_fd_modes) (
        SphHarmFrequencySeries **hlm,
        LALDict *params,
        LALSimInspiralGenerator *myself
    );

    int (*generate_fd_waveform) (
        COMPLEX16FrequencySeries **hplus,
        COMPLEX16FrequencySeries **hcross,
        LALDict *params,
        LALSimInspiralGenerator *myself
    );

    /* ... */
    void *internal_data;
};

#endif
