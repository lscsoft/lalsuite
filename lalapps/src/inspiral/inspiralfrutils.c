/*----------------------------------------------------------------------- 
 * 
 * File Name: inspiralfrutils.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


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

  snprintf( chname, sizeof(chname), "%s_%s", chan->name, suffix );
    fdata.name = chname;
    fdata.tbeg = chan->epoch;
    memset( &fdata.tend, 0, sizeof(LIGOTimeGPS) );
    epoch_add( &fdata.tend, &(chan->epoch), 
        chan->deltaT * (REAL8) chan->data->length );
    fdata.dom = Time;
    fdata.type = FR_VECT_4R;
    fdata.step = (float) chan->deltaT;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.data = (float *) chan->data->data;
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
    fdata.step = (float) chan->deltaF;
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
    fdata.step = (float) chan->deltaF;
    fdata.unit = unit;
    fdata.size = (size_t) chan->data->length;
    fdata.data = (float *) chan->data->data;
    return fr_add_proc_data( frame, &fdata );
}
