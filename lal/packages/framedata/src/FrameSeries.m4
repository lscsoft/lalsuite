/**** <lalVerbatim file="FrameSeriesCV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{FrameSeries.c}}
 *
 * These routines are used to read/write time series data from/to frame files.
 *
 * \subsubsection*{Prototypes}
 * \input{FrameSeriesCP}
 * \idx{LALFrGetINT2TimeSeries}
 * \idx{LALFrGetINT4TimeSeries}
 * \idx{LALFrGetINT8TimeSeries}
 * \idx{LALFrGetREAL4TimeSeries}
 * \idx{LALFrGetREAL8TimeSeries}
 * \idx{LALFrWriteREAL4TimeSeries}
 * \idx{LALFrGetTimeSeriesType}
 *
 * \subsubsection*{Description}
 *
 * The routines
 * \texttt{LALFrGet}$\langle\mbox{datatype}\rangle$\texttt{TimeSeries()}
 * search the frame for a specified channel.  If the time series supplied has
 * data storage allocated, then the specified amount of data is filled from
 * the frame stream.  If no space has been allocated (so that the data field
 * is \texttt{NULL}), then only the channel information is returned in the
 * time series (e.g., the start time of the next data and the time step size).
 *
 * The routine \texttt{LALFrGetTimeSeriesType} returns the type of the time
 * series corresponding to the specified channel.  This is needed if it is not
 * known in advance what type of data is stored in the channel.
 *
 * The routines
 * \texttt{LALFrWrite}$\langle\mbox{datatype}\rangle$\texttt{TimeSeries()}
 * outputs a given time series as a new frame file with the filename
 * $\langle\mbox{source}\rangle$\verb+-+$\langle\mbox{description}\rangle$%
 * \verb+-+$\langle\mbox{GPS-seconds}\rangle$\verb+-+%
 * $\langle\mbox{duration}\rangle$\verb+.gwf+
 * where source and description are the specified frame source and description
 * identifiers, GPS-seconds is the start time of the data in
 * seconds since 0h UTC 6 Jan 1980 (GPS time origin), or the nearest second
 * before the start of the data, and duration is the number of seconds between
 * value of GPS-seconds and the GPS end time of the data, rounded up.
 *
 * \vfill{\footnotesize\input{FrameSeriesCV}}
 *
 **** </lalLaTeX> */ 



#include <FrameL.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <FrameStreamDef.h>

NRCSID( FRAMESERIESC, "$Id$" );


/*
 *
 * Routine to load a FrVect associated with a given channel from a given frame.
 *
 */
static struct FrVect *loadFrVect( struct FrameH *frame, char *name, int chtype )
{
  struct FrVect *vect = NULL;
  if ( chtype == ProcDataChannel )
  {
    struct FrProcData *proc = FrProcDataFind( frame, name );
    vect = proc ? proc->data : NULL;
  }
  else if ( chtype == ADCDataChannel )
  {
    struct FrAdcData *adc = FrAdcDataFind( frame, name );
    vect = adc ? adc->data : NULL;
  }
  else if ( chtype == SimDataChannel )
  {
    struct FrSimData *sim = FrSimDataFind( frame, name );
    vect = sim ? sim->data : NULL;
  }
  return vect;
}


/*
 *
 * Routine to make a 1D FrVect structure and attach it to a frame as the
 * appropriate channel type.
 *
 */
static struct FrVect *makeFrVect1D( struct FrameH *frame, int chtype,
    char *name, char *comment, char *unitx, char *unity, int datatype,
    double rate, double fshift, double dx, unsigned int npts )
{
  struct FrVect *vect = FrVectNew1D( name, datatype, npts, dx, unitx, unity );
  if ( ! vect ) return NULL;
  if ( chtype == ProcDataChannel )
  {
    struct FrProcData *proc = calloc( 1, sizeof( *proc ) );
    if ( ! proc )
    {
      FrVectFree( vect );
      return NULL;
    }
    proc->classe = FrProcDataDef();
    proc->fShift = fshift;
    proc->data = vect;
    proc->next = frame->procData;
    frame->procData = proc;
    if ( ! FrStrCpy( &proc->name, name ) ) return NULL;
    if ( ! FrStrCpy( &proc->comment, comment ) ) return NULL;
  }
  else if ( chtype == ADCDataChannel )
  {
    struct FrAdcData *adc = calloc( 1, sizeof( *adc ) );
    struct FrRawData *raw = frame->rawData ? frame->rawData :
      FrRawDataNew( frame );
    if ( ! adc || ! raw )
    {
      FrVectFree( vect );
      if ( adc ) free( adc );
      return NULL;
    }
    adc->classe = FrAdcDataDef();
    adc->sampleRate = rate;
    adc->fShift = fshift;
    adc->data = vect;
    adc->next = raw->firstAdc;
    raw->firstAdc = adc;
    if ( ! FrStrCpy( &adc->name, name ) ) return NULL;
    if ( ! FrStrCpy( &adc->comment, comment ) ) return NULL;
    if ( ! FrStrCpy( &adc->units, unity ) ) return NULL;
  }
  else if ( chtype == SimDataChannel )
  {
    struct FrSimData *sim = calloc( 1, sizeof( *sim ) );
    if ( ! sim )
    {
      FrVectFree( vect );
      return NULL;
    }
    sim->classe = FrSimDataDef();
    sim->sampleRate = rate;
    sim->data = vect;
    sim->next = frame->simData;
    frame->simData = sim;
    if ( ! FrStrCpy( &sim->name, name ) ) return NULL;
    if ( ! FrStrCpy( &sim->comment, comment ) ) return NULL;
  }
  else
  {
    FrVectFree( vect );
    return NULL;
  }
  return vect;
}


/* <lalVerbatim file="FrameSeriesCP"> */
void
LALFrGetTimeSeriesType(
    LALStatus   *status,
    LALTYPECODE *output,
    FrChanIn    *chanin,
    FrStream    *stream
    )
{ /* </lalVerbatim> */
  INT4 type = -1;
  INITSTATUS( status, "FUNC", FRAMESERIESC );  

  ASSERT( output, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( chanin, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( chanin->name, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->err )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->end )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  if ( chanin->type == ProcDataChannel )
  {
    FrProcData *proc = stream->frame ? stream->frame->procData : NULL;
    while ( proc && proc->name && strcmp( proc->name, chanin->name ) )
      proc = proc->next;
    type = ( proc && proc->data ) ? proc->data->type : -1;
  }
  else if ( chanin->type == ADCDataChannel )
  {
    FrAdcData *adc = ( stream->frame && stream->frame->rawData )
      ? stream->frame->rawData->firstAdc : NULL;
    while ( adc && adc->name && strcmp( adc->name, chanin->name ) )
      adc = adc->next;
    type = ( adc && adc->data ) ? adc->data->type : -1;
  }
  else if ( chanin->type == SimDataChannel )
  {
    FrSimData *sim = stream->frame ? stream->frame->simData : NULL;
    while ( sim && sim->name && strcmp( sim->name, chanin->name ) )
      sim = sim->next;
    type = ( sim && sim->data ) ? sim->data->type : -1;
  }

  switch ( type )
  {
    case FR_VECT_C:
      *output = LAL_CHAR_TYPE_CODE;
      break;
    case FR_VECT_2S:
      *output = LAL_I2_TYPE_CODE;
      break;
    case FR_VECT_4S:
      *output = LAL_I4_TYPE_CODE;
      break;
    case FR_VECT_8S:
      *output = LAL_I8_TYPE_CODE;
      break;
    case FR_VECT_1U:
      *output = LAL_UCHAR_TYPE_CODE;
      break;
    case FR_VECT_2U:
      *output = LAL_U2_TYPE_CODE;
      break;
    case FR_VECT_4U:
      *output = LAL_U4_TYPE_CODE;
      break;
    case FR_VECT_8U:
      *output = LAL_U8_TYPE_CODE;
      break;
    case FR_VECT_4R:
      *output = LAL_S_TYPE_CODE;
      break;
    case FR_VECT_8R:
      *output = LAL_D_TYPE_CODE;
      break;
    case FR_VECT_C8:
      *output = LAL_C_TYPE_CODE;
      break;
    case FR_VECT_C16:
      *output = LAL_Z_TYPE_CODE;
      break;
    default:
      ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );      
  }

  RETURN( status );
}


define(`TYPE',`REAL8')
include(`FrameSeriesRead.m4')
include(`FrameSeriesWrite.m4')

define(`TYPE',`REAL4')
include(`FrameSeriesRead.m4')
include(`FrameSeriesWrite.m4')

define(`TYPE',`INT8')
include(`FrameSeriesRead.m4')
include(`FrameSeriesWrite.m4')

define(`TYPE',`INT4')
include(`FrameSeriesRead.m4')
include(`FrameSeriesWrite.m4')

define(`TYPE',`INT2')
include(`FrameSeriesRead.m4')
include(`FrameSeriesWrite.m4')
