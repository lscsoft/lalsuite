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
 * File Name: lalappsfrutils.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "lalappsfrutils.h"

/*
 *
 * functions to write data to frame files, since lal doesn't have this yet
 *
 */


FrameH *fr_add_proc_REAL4TimeSeries ( 
    FrameH                     *frame, 
    REAL4TimeSeries            *chan,
    const char                 *unit,
    const char                 *suffix
    )
{
  char          chname[256];
  struct        series fdata;
  
  if ( suffix )
  {
    snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
  }
  else
  {
    snprintf( chname, sizeof(chname), "%s", chan->name );
  } 
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &(chan->epoch), 
        chan->deltaT * (REAL8) chan->data->length );
    fdata.dom = Time;
    fdata.type = FR_VECT_4R;
    fdata.step = (double) chan->deltaT;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.data = (float *) chan->data->data;
    return fr_add_proc_data( frame, &fdata );
}

FrameH *fr_add_proc_REAL8TimeSeries ( 
    FrameH                     *frame, 
    REAL8TimeSeries            *chan,
    const char                 *unit,
    const char                 *suffix
    )
{
  char          chname[256];
  struct        series fdata;
  
  if ( suffix )
  {
    snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
  }
  else
  {
    snprintf( chname, sizeof(chname), "%s", chan->name );
  } 
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &(chan->epoch), 
        chan->deltaT * (REAL8) chan->data->length );
    fdata.dom = Time;
    fdata.type = FR_VECT_8R;
    fdata.step = (double) chan->deltaT;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.ddata = (double *) chan->data->data;
    return fr_add_proc_data( frame, &fdata );
}


FrameH *fr_add_proc_REAL4FrequencySeries ( 
    FrameH                     *frame, 
    REAL4FrequencySeries       *chan,
    const char                 *unit,
    const char                 *suffix
    )
{
  char          chname[256];
  struct        series fdata;

  snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &chan->epoch, 
        (chan->data->length - 1) / (chan->deltaF * chan->data->length) );
    fdata.dom = Freq;
    fdata.type = FR_VECT_4R;
    fdata.step = (double) chan->deltaF;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.data = (float *) chan->data->data;
    return fr_add_proc_data( frame, &fdata );
}

FrameH *fr_add_proc_COMPLEX8FrequencySeries ( 
    FrameH                        *frame, 
    COMPLEX8FrequencySeries       *chan,
    const char                    *unit,
    const char                    *suffix
    )
{
  char          chname[256];
  struct        series fdata;

  snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &chan->epoch, 
        (chan->data->length - 1) / (chan->deltaF * chan->data->length) );
    fdata.dom = Freq;
    fdata.type = FR_VECT_8C;
    fdata.step = (double) chan->deltaF;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.data = (float *) chan->data->data;
    return fr_add_proc_data( frame, &fdata );
}

FrameH *fr_add_proc_COMPLEX8TimeSeries (
    FrameH                        *frame,
    COMPLEX8TimeSeries            *chan,
    const char                    *unit,
    const char                    *suffix
    ) 
{
  char          chname[256];
  struct        series fdata;

  snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &chan->epoch, 
	chan->deltaT * (REAL8) chan->data->length );
    fdata.dom = Time;
    fdata.type = FR_VECT_8C;
    fdata.step = (float) chan->deltaT;
    fdata.f0 = (float) chan->f0;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.data = (float *) chan->data->data;
    return fr_add_proc_data( frame, &fdata );
}



FrameH *fr_add_proc_REAL8FrequencySeries ( 
    FrameH                     *frame, 
    REAL8FrequencySeries       *chan,
    const char                 *unit,
    const char                 *suffix
    )
{
  char          chname[256];
  struct        series fdata;
  char comment[] = "Generated by $Id$";
  char seconds[] = "s";
  char hertz[]   = "Hz";
  struct FrVect     *vect;
  struct FrProcData *proc;
  size_t i;
  const char *channel = fdata.name;

  snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &chan->epoch, 
        (chan->data->length - 1) / (chan->deltaF * chan->data->length) );
    fdata.dom = Freq;
    fdata.type = FR_VECT_8R;
    fdata.step = (double) chan->deltaF;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;

  if ( ! frame )
  {
    char src[2];
    src[0] = fdata.name[0];
    src[1] = 0;
    frame = FrameHNew( src );
    frame->run    = 1;
    frame->frame  = 1;
    frame->GTimeS = (int) fdata.tbeg.gpsSeconds;
    frame->GTimeN = (int) fdata.tbeg.gpsNanoSeconds;
    frame->dt     = epoch_diff( &fdata.tend, &fdata.tbeg );
  }

  vect = FrVectNew1D( channel, fdata.type, fdata.size, fdata.step,
      IS_TIME( fdata.dom) ? seconds : hertz, fdata.unit );
  proc = calloc( 1, sizeof( *proc ) );
  proc->classe     = FrProcDataDef();
#if defined FR_VERS && FR_VERS < 5000
  proc->sampleRate = IS_TIME( fdata.dom) ? 1.0 / fdata.step : -1;
#endif
  proc->fShift     = 0;
  proc->data       = vect;
  proc->next       = frame->procData;
  frame->procData  = proc;
  FrStrCpy( &proc->name, channel );
  FrStrCpy( &proc->comment, comment );
  for ( i = 0; i < fdata.size; ++i )
  {
    vect->dataF[i] = fdata.data[i];
  }

  return frame;
}
