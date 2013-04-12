/*
*  Copyright (C) 2007 Jolien Creighton
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

int isinf( double );
int isnan( double );
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdio.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include <lal/LALString.h>
#include "LALASCIIFileRead.h"


#ifndef LINE_MAX
#define LINE_MAX 2048
#endif


int XLALDataFileNameParse( LALDataFileNameFields *fields, const char *fname )
{
  const char *fmt = "%[^-]-%[^-]-%u-%u.%s";
  const char *basename;
  UINT4 t0;
  UINT4 dt;
  int c;
  if ( ! fields || ! fname )
    XLAL_ERROR( XLAL_EFAULT );
  basename = strrchr( fname, '/' ); /* find the last slash */
  basename = basename ? basename + 1 : fname;
  if ( strlen( basename ) > FILENAME_MAX )
  {
    XLALPrintError( "XLAL Error - %s: Filename %s too long\n", __func__, basename );
    XLAL_ERROR( XLAL_EBADLEN );
  }
  c = sscanf( basename, fmt, fields->site, fields->description, &t0, &dt, fields->extension );
  if ( c != 5 )
  {
    XLALPrintError( "XLAL Error - %s: Could not parse basename %s\n\tinto <site>-<description>-<tstart>-<duration>.<extension> format\n", __func__, basename );
    XLAL_ERROR( XLAL_EINVAL );
  }
  fields->tstart   = t0;
  fields->duration = dt;
  return 0;
}

int XLALCalRefFileNameDescriptionParse( LALCalRefFileNameDescriptionFields *fields, const char *description )
{
  char dsc[FILENAME_MAX];
  char *dsc_channel_start;
  char *dsc_version_start;
  char *dsc_run_start;
  char *dsc_ifo_start;
  int c;

  if ( ! fields || ! description )
    XLAL_ERROR( XLAL_EFAULT );
  if ( strlen( description ) > FILENAME_MAX )
  {
    XLALPrintError( "XLAL Error - %s: Description string %s too long\n", __func__, description );
    XLAL_ERROR( XLAL_EBADLEN );
  }

  XLALStringCopy( dsc, description, sizeof( dsc ) );
  dsc_ifo_start = dsc;

  /* find start of channel string */
  dsc_channel_start = strchr( dsc, '_' );
  if ( ! dsc_channel_start )
  {
    XLALPrintError( "XLAL Error - %s: Could not parse description field %s\n\tinto <ifo>_CAL_REF_<channelpostfix>_<run>_<version> format\n", __func__, description );
    XLAL_ERROR( XLAL_EINVAL );
  }
  *dsc_channel_start++ = 0;
  /* now skip the CAL_REF_ */
  if ( dsc_channel_start != strstr( dsc_channel_start, "CAL_REF_" ) )
  {
    XLALPrintError( "XLAL Error - %s: Could not parse description field %s\n\tinto <ifo>_CAL_REF_<channelpostfix>_<run>_<version> format\n", __func__, description );
    XLAL_ERROR( XLAL_EINVAL );
  }
  dsc_channel_start += strlen( "CAL_REF" );
  *dsc_channel_start++ = 0;

  /* find start of version string */
  dsc_version_start = strrchr( dsc_channel_start, '_' );
  if ( ! dsc_version_start )
  {
    XLALPrintError( "XLAL Error - %s: Could not parse description field %s\n\tinto <ifo>_CAL_REF_<channelpostfix>_<run>_<version> format\n", __func__, description );
    XLAL_ERROR( XLAL_EINVAL );
  }
  *dsc_version_start++ = 0;

  /* find start of run string */
  dsc_run_start = strrchr( dsc_channel_start, '_' );
  if ( ! dsc_run_start )
  {
    XLALPrintError( "XLAL Error - %s: Could not parse description field %s\n\tinto <ifo>_CAL_REF_<channelpostfix>_<run>_<version> format\n", __func__, description );
    XLAL_ERROR( XLAL_EINVAL );
  }
  *dsc_run_start++ = 0;

  XLALStringCopy( fields->ifo, dsc_ifo_start, sizeof( fields->ifo ) );
  XLALStringCopy( fields->channelPostfix, dsc_channel_start, sizeof( fields->channelPostfix ) );
  XLALStringCopy( fields->run, dsc_run_start, sizeof( fields->run ) );
  c = sscanf( dsc_version_start, "V%d", &fields->version );
  if ( c != 1 )
  {
    if ( ! strcmp( dsc_version_start, "VU" ) )
      fields->version = -1; /* means unity */
    else
    {
      XLALPrintError( "XLAL Error - %s: Could not parse description field %s\n\tinto <ifo>_CAL_REF_<channelpostfix>_<run>_<version> format\n", __func__, description );
      XLAL_ERROR( XLAL_EINVAL );
    }
  }
  return 0;
}

int XLALCalFacFileNameDescriptionParse( LALCalFacFileNameDescriptionFields *fields, const char *description )
{
  const char *fmt1 = "%[^_]_CAL_FAC_%[^_]_V%d_%d";
  const char *fmt2 = "%[^_]_CAL_FAC_%[^_]_VU_%d";
  int c;
  if ( ! fields || ! description )
    XLAL_ERROR( XLAL_EFAULT );
  if ( strlen( description ) > FILENAME_MAX )
  {
    XLALPrintError( "XLAL Error - %s: Description string %s too long\n", __func__, description );
    XLAL_ERROR( XLAL_EBADLEN );
  }
  c = sscanf( description, fmt1, fields->ifo, fields->run, &fields->version, &fields->deltaT );
  if ( c != 4 )
  {
    c = sscanf( description, fmt2, fields->ifo, fields->run, fields->deltaT );
    if ( c != 3 )
    {
      XLALPrintError( "XLAL Error - %s: Could not parse description field %s\n\tinto <ifo>_CAL_FAC_<run>_<version>_<deltat> format\n", __func__, description );
      XLAL_ERROR( XLAL_EINVAL );
    }
    fields->version = -1; /* means unity */
  }
  return 0;
}


int XLALASCIIFileCountRows( const char *fname )
{
  char line[LINE_MAX];
  int nline;
  int nrow;
  FILE *fp;

  fp = fopen( fname, "r" );
  if ( ! fp )
    XLAL_ERROR( XLAL_EIO );

  nline = 0;
  nrow = 0;
  while ( 1 )
    if ( fgets( line, sizeof( line ), fp ) )
    {
      ++nline;
      /* check that the line is not too long */
      if ( strlen( line ) >= sizeof( line ) - 1 )
      {
        fclose( fp );
        XLALPrintError( "XLAL Error - %s: line %d too long\n\tfile: %s\n", __func__, nline, fname );
        XLAL_ERROR( XLAL_EBADLEN );
      }
      /* check to see if this line is a comment line */
      if ( line[0] == '#' || line[0] == '%' )
        continue; /* ignore these */
      ++nrow;
    }
    else if ( feof( fp ) ) /* we're done */
      break;
    else /* must have been a file reading error */
    {
      fclose( fp );
      XLAL_ERROR( XLAL_EIO );
    }

  fclose( fp );
  return nrow;
}


REAL8VectorSequence * XLALASCIIFileReadColumns( INT4 ncol, const char *fname )
{
  char line[LINE_MAX];
  REAL8VectorSequence *data;
  int nline;
  int nrow;
  int row;
  int col;
  FILE *fp;

  if ( ncol < 0 )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  /* count rows */
  nrow = XLALASCIIFileCountRows( fname );
  if ( nrow < 0 )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  /* allocate memory for data */
  /* column 0 will contain line number of input file */
  data = XLALCreateREAL8VectorSequence( nrow, ncol + 1 );
  if ( ! data )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  /* open file */
  fp = fopen( fname, "r" );
  if ( ! fp )
  {
    XLALDestroyREAL8VectorSequence( data );
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  nline = 0;
  row = 0;
  while ( row < nrow )
    if ( fgets( line, sizeof( line ), fp ) )
    {
      char *p = line;
      ++nline;
      if ( strlen( line ) >= sizeof( line ) - 1 )
      {
        XLALPrintError( "XLAL Error - %s: line %d too long\n\tfile: %s\n", __func__, nline, fname );
        XLALDestroyREAL8VectorSequence( data );
        fclose( fp );
        XLAL_ERROR_NULL( XLAL_EBADLEN );
      }
      if ( line[0] == '#' || line[0] == '%' )
        continue;
      data->data[row*data->vectorLength] = nline;
      for ( col = 1; col <= ncol; ++col )
      {
        char *endp;
        REAL8 val;
        val = strtod( p, &endp );
        if ( p == endp ) /* no conversion */
        {
          XLALPrintError( "XLAL Error - %s: unable to parse line %d\n\tfile: %s\n", __func__, nline, fname );
          XLALDestroyREAL8VectorSequence( data );
          fclose( fp );
          XLAL_ERROR_NULL( XLAL_EFAILED );
        }
        if ( isnan( val ) )
        {
          XLALPrintWarning( "XLAL Warning - %s: invalid data (nan) in column %d of line %d\n\tfile: %s\n", __func__, col, nline, fname );
          val = 0;
        }
        if ( isinf( val ) )
        {
          XLALPrintWarning( "XLAL Warning - %s: invalid data (inf) in column %d of line %d\n\tfile: %s\n", __func__, col, nline, fname );
          val = 0;
        }
        data->data[row*data->vectorLength + col] = val;
        p = endp;
      }
      ++row;
    }
    else
    {
      XLALDestroyREAL8VectorSequence( data );
      XLAL_ERROR_NULL( XLAL_EIO );
    }

  fclose( fp );
  return data;
}


REAL4 XLALASCIIFileReadCalFacHeader( const char *fname )
{
  static const char *fmt1 = "%% deltaT = %f\n"; /* format of line 1 */
  static const char *fmt2 = "%% ""$""Name: %[^$]""$""\n"; /* format of line 2 */
  char line[LINE_MAX];
  char cvsname[LINE_MAX];
  REAL4 deltaT;
  FILE *fp;
  char *rc;
  int c;

  fp = fopen( fname, "r" );
  if ( ! fp )
    XLAL_ERROR_REAL4( XLAL_EIO );

  rc = fgets( line, sizeof( line ), fp );
  if ( rc == NULL )
  {
    XLALPrintError( "XLAL Error - %s: unable to read file: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }

  c = sscanf( line, fmt1, &deltaT );
  if ( c != 1 ) /* wrong number of conversions */
  {
    XLALPrintError( "XLAL Error - %s: incorrect first header line\n\tfile: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }

  rc = fgets( line, sizeof( line ), fp );
  if ( rc == NULL )
  {
    XLALPrintError( "XLAL Error - %s: unable to read file: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }

  c = sscanf( line, fmt2, cvsname );
  if ( c != 1 ) /* wrong number of conversions */
  {
    XLALPrintError( "XLAL Error - %s: incorrect second header line\n\tfile: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }

  XLALPrintInfo( "XLAL Info - %s: CVS tag used to produce factors data file: \"%s\"\n", __func__, cvsname );

  fclose( fp );
  return deltaT;
}


REAL4 XLALASCIIFileReadCalRefHeader( const char *fname )
{
  static const char *fmt1 = "%% deltaF = %f\n"; /* format of line 1 */
  static const char *fmt2 = "%% ""$""Name: %[^$]""$""\n"; /* format of line 2 */
  char line[LINE_MAX];
  REAL4 deltaF;
  char cvsname[LINE_MAX];
  FILE *fp;
  char *rc;
  int c;

  fp = fopen( fname, "r" );
  if ( ! fp )
    XLAL_ERROR_REAL4( XLAL_EIO );

  rc = fgets( line, sizeof( line ), fp );
  if ( rc == NULL )
  {
    XLALPrintError( "XLAL Error - %s: unable to read file: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }

  c = sscanf( line, fmt1, &deltaF );
  if ( c != 1 ) /* wrong number of conversions */
  {
    XLALPrintError( "XLAL Error - %s: incorrect first header line\n\tfile: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }

  rc = fgets( line, sizeof( line ), fp );
  if ( rc == NULL )
  {
    XLALPrintError( "XLAL Error - %s: unable to read file: %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EFAILED );
  }
  c = sscanf( line, fmt2, cvsname );
  if ( c != 1 ) /* wrong number of conversions */
  {
    XLALPrintError( "XLAL Error - %s: incorrect second header line\n\tfile: %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAILED );
  }

  XLALPrintInfo( "XLAL Info - %s: CVS tag used to produce factors data file: \"%s\"\n", __func__, cvsname, fname );

  fclose( fp );
  return deltaF;
}

#define ELEM(data,row,col) ((data)->data[(row)*(data)->vectorLength+(col)])

int XLALASCIIFileReadCalFac( REAL4TimeSeries **alpha, REAL4TimeSeries **lal_gamma, const char *fname )
{
  const REAL8 fuzzfactor = 1e-3; /* fraction of a sample of fuzziness */
  REAL8 fuzz; /* possible discrepancies in times */
  char alphaName[] = "Xn:CAL-CAV_FAC";
  char gammaName[] = "Xn:CAL-OLOOP_FAC";
  LALDataFileNameFields              fileFields;
  LALCalFacFileNameDescriptionFields descFields;
  REAL8VectorSequence *data;
  LIGOTimeGPS epoch;
  REAL8 tstart;
  REAL8 tend;
  REAL4 deltaT;
  int ndat;
  int nrow;
  int row;
  int dat;

  if ( ! alpha || ! lal_gamma )
    XLAL_ERROR( XLAL_EFAULT );
  if ( *alpha || *lal_gamma )
    XLAL_ERROR( XLAL_EINVAL );

  if ( XLALDataFileNameParse( &fileFields, fname ) < 0 )
  {
    XLALPrintError( "XLAL Error - %s: invalid file name %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EINVAL );
  }

  if ( XLALCalFacFileNameDescriptionParse( &descFields, fileFields.description ) < 0 )
  {
    XLALPrintError( "XLAL Error - %s: invalid description part of file name %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EINVAL );
  }

  XLALPrintInfo( "XLAL Info - %s: Reading calibration factors from file %s\n", __func__, fname );

  /* setup channel names */
  memcpy( alphaName, descFields.ifo, 2 );
  memcpy( gammaName, descFields.ifo, 2 );

  deltaT = XLALASCIIFileReadCalFacHeader( fname );
  if ( XLAL_IS_REAL4_FAIL_NAN( deltaT ) )
    XLAL_ERROR( XLAL_EFUNC );

  if ( deltaT <= 0 ) /* bad time step */
    XLAL_ERROR( XLAL_EDATA );

  /* check consistency of time steps */
  if ( fabs( descFields.deltaT - deltaT ) > 1 ) /* step in file name might only be accurate to one second */
    XLAL_ERROR( XLAL_EDATA );

  /* read three columns */
  data = XLALASCIIFileReadColumns( 3, fname );
  if ( ! data )
    XLAL_ERROR( XLAL_EFUNC );
  if ( ! data->length )
  {
    XLALPrintError( "XLAL Error - %s: no rows of data\n", __func__ );
    XLAL_ERROR( XLAL_EFAILED );
  }
  nrow = data->length;

  /* sanity checks on time stamps */
  /* check for mismatch in start or end times compared to file name */
  tstart = ELEM(data,0,1);
  tend   = ELEM(data,nrow-1,1);
  ndat   = 1 + (int)floor( (tend - tstart)/deltaT + 0.5 );
  if ( ( fileFields.tstart != (int)floor( tstart ) )
    || ( fileFields.tstart + fileFields.duration != (int)ceil( tend + deltaT ) ) )
  {
    XLALPrintError( "XLAL Error - %s: filename start time and duration not consistent with contents\n", __func__ );
    XLALDestroyREAL8VectorSequence( data );
    XLAL_ERROR( XLAL_EDATA );
  }
  if ( nrow > ndat ) /* this should be impossible */
  {
    XLALDestroyREAL8VectorSequence( data );
    XLAL_ERROR( XLAL_EDATA );
  }

  XLALGPSSetREAL8( &epoch, tstart );

  *alpha = XLALCreateREAL4TimeSeries( alphaName, &epoch, 0.0, deltaT, &lalDimensionlessUnit, ndat );
  *lal_gamma = XLALCreateREAL4TimeSeries( gammaName, &epoch, 0.0, deltaT, &lalDimensionlessUnit, ndat );
  if ( ! *alpha || ! *lal_gamma )
  {
    XLALDestroyREAL4TimeSeries( *lal_gamma );
    XLALDestroyREAL4TimeSeries( *alpha );
    XLALDestroyREAL8VectorSequence( data );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* clear the data memory */
  memset( (*alpha)->data->data, 0, (*alpha)->data->length * sizeof( *(*alpha)->data->data ) );
  memset( (*lal_gamma)->data->data, 0, (*lal_gamma)->data->length * sizeof( *(*lal_gamma)->data->data ) );

  /* IMPORTANT: SPECIFICATION SAYS THAT ALPHA IS COL 2 AND GAMMA IS COL 3 */
  fuzz = fuzzfactor * deltaT;
  dat = -1;
  for ( row = 0; row < nrow; ++row )
  {
    int   line     = ELEM(data,row,0);
    REAL8 trow     = ELEM(data,row,1);
    REAL4 alphaval = ELEM(data,row,2);
    REAL4 gammaval = ELEM(data,row,3);
    int thisdat = (int)floor( (trow - tstart) / deltaT + 0.5 );
    if ( thisdat <= dat ) /* rows must be monotonically increasing in time */
    {
      XLALPrintError( "XLAL Error - %s: error on line %d of file %s\n\trows must be monotonically increasing in time\n", __func__, line, fname );
      XLALDestroyREAL4TimeSeries( *lal_gamma );
      XLALDestroyREAL4TimeSeries( *alpha );
      XLALDestroyREAL8VectorSequence( data );
      XLAL_ERROR( XLAL_EDATA );
    }
    dat = thisdat;
    if ( fabs( (tstart + dat * deltaT) - trow ) > fuzz ) /* time between rows must be an integral multiple of deltaT */
    {
      XLALPrintError( "XLAL Error - %s: error on line %d of file %s\n\ttimes must be integral multiples of deltaT\n", __func__, line, fname );
      XLALDestroyREAL4TimeSeries( *lal_gamma );
      XLALDestroyREAL4TimeSeries( *alpha );
      XLALDestroyREAL8VectorSequence( data );
      XLAL_ERROR( XLAL_EDATA );
    }
    if ( dat >= ndat ) /* beyond length of array */
    {
      printf( "%d\t%d\n", dat, ndat );
      XLALPrintError( "XLAL Error - %s: error on line %d of file %s\n\ttime beyond end time\n", __func__, line, fname );
      XLALDestroyREAL4TimeSeries( *lal_gamma );
      XLALDestroyREAL4TimeSeries( *alpha );
      XLALDestroyREAL8VectorSequence( data );
      XLAL_ERROR( XLAL_EDATA );
    }
    (*alpha)->data->data[dat] = alphaval;
    (*lal_gamma)->data->data[dat] = gammaval;
  }

  XLALDestroyREAL8VectorSequence( data );
  return 0;
}


int XLALASCIIFileReadCalRef( COMPLEX8FrequencySeries **series, REAL8 *duration, const char *fname )
{
  const REAL8 fuzzfactor = 1e-3; /* fraction of a sample of fuzziness */
  REAL8 fuzz; /* possible discrepancies in times */
  LALDataFileNameFields              fileFields;
  LALCalRefFileNameDescriptionFields descFields;
  REAL8VectorSequence *data;
  LIGOTimeGPS epoch;
  LALUnit unit;
  char channel[FILENAME_MAX];
  REAL4 deltaF;
  REAL8 fstart;
  REAL8 fend;
  REAL8 df;
  int ndat;
  int nrow;
  int row;

  if ( ! series )
    XLAL_ERROR( XLAL_EFAULT );
  if ( *series )
    XLAL_ERROR( XLAL_EINVAL );

  if ( XLALDataFileNameParse( &fileFields, fname ) < 0 )
  {
    XLALPrintError( "XLAL Error - %s: invalid file name %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EINVAL );
  }
  *duration = fileFields.duration;

  if ( XLALCalRefFileNameDescriptionParse( &descFields, fileFields.description ) < 0 )
  {
    XLALPrintError( "XLAL Error - %s: invalid description part of file name %s\n", __func__, fname );
    XLAL_ERROR( XLAL_EINVAL );
  }

  XLALPrintInfo( "XLAL Info - %s: Reading calibration factors from file %s\n", __func__, fname );

  /* setup units */
  if ( strstr( descFields.channelPostfix, "RESPONSE" ) )
    XLALUnitDivide( &unit, &lalStrainUnit, &lalADCCountUnit );
  else if ( strstr( descFields.channelPostfix, "CAV_GAIN" ) )
    XLALUnitDivide( &unit, &lalADCCountUnit, &lalStrainUnit );
  else if ( strstr( descFields.channelPostfix, "OLOOP_GAIN" ) )
    unit = lalDimensionlessUnit;
  else if ( strstr( descFields.channelPostfix, "ACTUATION" ) )
    XLALUnitDivide( &unit, &lalADCCountUnit, &lalStrainUnit );
  else if ( strstr( descFields.channelPostfix, "DIGFLT" ) )
    unit = lalDimensionlessUnit;
  else
  {
    XLALPrintError( "XLAL Error - %s: invalid channel %s\n", __func__, descFields.channelPostfix );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* setup channel names */
  XLALStringCopy( channel, descFields.ifo, sizeof( channel ) );
  XLALStringConcatenate( channel, ":CAL-", sizeof( channel ) );
  XLALStringConcatenate( channel, descFields.channelPostfix, sizeof( channel ) );

  /* read header */
  deltaF = XLALASCIIFileReadCalRefHeader( fname );
  if ( XLAL_IS_REAL4_FAIL_NAN( deltaF ) )
    XLAL_ERROR( XLAL_EFUNC );

  if ( deltaF <= 0 ) /* bad time step */
    XLAL_ERROR( XLAL_EDATA );


  /* read three columns */
  data = XLALASCIIFileReadColumns( 3, fname );
  if ( ! data )
    XLAL_ERROR( XLAL_EFUNC );
  if ( ! data->length )
  {
    XLALPrintError( "XLAL Error - %s: no rows of data\n", __func__ );
    XLAL_ERROR( XLAL_EFAILED );
  }
  nrow = data->length;

  /* sanity checks on time stamps */
  fstart = ELEM(data,0,1);
  fend   = ELEM(data,nrow-1,1);
  ndat   = 1 + (int)floor( (fend - fstart)/deltaF + 0.5 );
  df     = (fend - fstart) / (ndat - 1);

  if ( nrow != ndat ) /* this should be impossible */
  {
    XLALDestroyREAL8VectorSequence( data );
    XLAL_ERROR( XLAL_EDATA );
  }

  XLALGPSSetREAL8( &epoch, fileFields.tstart );

  *series = XLALCreateCOMPLEX8FrequencySeries( channel, &epoch, fstart, df, &unit, ndat );
  if ( ! *series )
  {
    XLALDestroyREAL8VectorSequence( data );
    XLAL_ERROR( XLAL_EFUNC );
  }

  fuzz = fuzzfactor * deltaF;
  for ( row = 0; row < nrow; ++row )
  {
    int   line = ELEM(data,row,0);
    REAL8 freq = ELEM(data,row,1);
    REAL8 mod  = ELEM(data,row,2);
    REAL8 arg  = ELEM(data,row,3);
    REAL8 fexp = fstart + row * df;
    if ( fabs( freq - fexp ) / fexp > fuzz )
    {
      XLALPrintError( "XLAL Error - %s: error on line %d of file %s\n\tunexpected frequency\n", __func__, line, fname );
      XLALDestroyCOMPLEX8FrequencySeries( *series );
      XLALDestroyREAL8VectorSequence( data );
      XLAL_ERROR( XLAL_EDATA );
    }
    (*series)->data->data[row].realf_FIXME = mod * cos( arg );
    (*series)->data->data[row].im = mod * sin( arg );
  }

  XLALDestroyREAL8VectorSequence( data );
  return 0;
}
