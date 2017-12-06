/*
*  Copyright (C) 2007 Sukanta Bose, Duncan Brown, Stephen Fairhurst
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

/*-----------------------------------------------------------------------
 *
 * File Name: lalappsfrutils.h
 *
 * Author: Brown, D. A.
 *
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALFrameL.h>

#include "series.h"

FrameH *fr_add_proc_REAL4TimeSeries (
    FrameH          *frame,
    REAL4TimeSeries *chan,
    const char      *unit,
    const char      *suffix
    );
FrameH *fr_add_proc_REAL8TimeSeries (
    FrameH          *frame,
    REAL8TimeSeries *chan,
    const char      *unit,
    const char      *suffix
    );
FrameH *fr_add_proc_REAL4FrequencySeries (
    FrameH               *frame,
    REAL4FrequencySeries *chan,
    const char           *unit,
    const char           *suffix
    );
FrameH *fr_add_proc_REAL8FrequencySeries (
    FrameH               *frame,
    REAL8FrequencySeries *chan,
    const char           *unit,
    const char           *suffix
    );
FrameH *fr_add_proc_COMPLEX8FrequencySeries (
    FrameH                  *frame,
    COMPLEX8FrequencySeries *chan,
    const char              *unit,
    const char              *suffix
    );
FrameH *fr_add_proc_COMPLEX8TimeSeries (
    FrameH                        *frame,
    COMPLEX8TimeSeries            *chan,
    const char                    *unit,
    const char                    *suffix
    );

