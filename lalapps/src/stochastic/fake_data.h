/*
 * fake_data.h - SGWB Standalone Analysis Pipeline
 *             - Fake Data Function Prototypes
 * 
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/SimulateSB.h>

#include <lalapps.h>

/* generate random noise */
REAL4TimeSeries *generate_random_noise(LALStatus *status,
    INT4 duration,
    INT4 sample_rate);

/* generate fake detector output */
SSSimStochBGOutput *generate_fake_detector_output(LALStatus *status,
    REAL4TimeSeries *noise_one,
    REAL4TimeSeries *noise_two,
    REAL8 deltaF,
    REAL8 f_min,
    REAL8 f_max);

/*
 * vim: et
 */
