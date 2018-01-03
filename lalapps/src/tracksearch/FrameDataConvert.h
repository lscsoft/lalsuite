 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
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
  *
  */

 /*
  * Author: Cristina Valeria Torres (LIGO Livingston)
  */

#ifndef _FRAMEDATACONVERTER_H
#define _FRAMEDATACONVERTER_H

#include <stdio.h>
#include <lal/LALDatatypes.h>
#include <lal/LALFrStream.h>

int XLALFrGetREAL8FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, LALFrStream *stream);
int XLALFrGetREAL4FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, LALFrStream *stream);
int XLALFrGetINT2FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, LALFrStream *stream);
int XLALFrGetINT4FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, LALFrStream *stream);
int XLALFrGetINT8FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, LALFrStream *stream);

#endif
