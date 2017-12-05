/*
 * fake_data.h - SGWB Standalone Analysis Pipeline
 *             - Fake Data Function Prototypes
 *
 * Copyright (C) 2002-2006 Adam Mercer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 */

#ifndef FAKE_DATA_H
#define FAKE_DATA_H

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

#endif /* FAKE_DATA_H */

/*
 * vim: et
 */
