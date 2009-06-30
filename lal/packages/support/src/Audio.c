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

/* vim: set noet ts=4 sw=4: */
/** \file
 * \ingroup support
 * \author Creighton, J. D. E.
 * \date $Date$
 * \brief Routines for exporting time series as sound files.
 *
 * $Id$
 */

#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Audio.h>

#define LAL_SOUND_MAX 32760

NRCSID( AUDIOC, "$Id$" );

/* if file pointer is NULL then create a new file based on series name */
static FILE * fopen_if_null( FILE *fp, const char *name, const char *extn )
{
	if ( ! fp ) {
		char fname[FILENAME_MAX];
		LALSnprintf( fname, sizeof( fname ), "%s.%s", name, extn );
		fp = fopen( fname, "w" );
	}
	return fp;
}

static int output_wav_hdr( FILE *fp, INT4 samplerate, UINT4 datasize )
{
	UINT4 totalsize = datasize + 36; /* for header */
	INT4 byterate = samplerate * sizeof( UINT2 );

	/*
	 *
	 * Write RIFF Header
	 *
	 */

	/* magic */
	fprintf( fp, "RIFF" );

	/* size of rest of file (bytes) */
	fputc( (totalsize & 0x000000ff)       , fp );
	fputc( (totalsize & 0x0000ff00) >> 8  , fp );
	fputc( (totalsize & 0x00ff0000) >> 16 , fp );
	fputc( (totalsize & 0xff000000) >> 24 , fp );

	/* file type "WAVE" */
	fprintf( fp, "WAVE" );


	/*
	 *
	 * Write Format Chunk
	 *
	 */


	/* magic */
	fprintf( fp, "fmt " );

	/* chunk size */
	fputc( 16, fp );
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 0, fp );

	/* format tag (1=uncompressed) */
	fputc( 1, fp );
	fputc( 0, fp );

	/* number of channels */
	fputc( 1, fp );
	fputc( 0, fp );

	/* samples per second */
	fputc( (samplerate & 0x000000ff)      , fp );
	fputc( (samplerate & 0x0000ff00) >> 8 , fp );
	fputc( (samplerate & 0x00ff0000) >> 16, fp );
	fputc( (samplerate & 0xff000000) >> 24, fp );

	/* average bytes per second */
	fputc( (byterate & 0x000000ff)      , fp );
	fputc( (byterate & 0x0000ff00) >> 8 , fp );
	fputc( (byterate & 0x00ff0000) >> 16, fp );
	fputc( (byterate & 0xff000000) >> 24, fp );

	/* block alignment (bytes of data per time step) */
	/* i.e., number of channels times sample width in bytes */
	fputc( 2, fp );
	fputc( 0, fp );

	/* bits per sample */
	fputc( 16, fp );
	fputc( 0, fp );


	/*
	 *
	 * Write Data Chunk
	 *
	 */


	/* magic */
	fprintf( fp, "data" );

	/* size of this chunk in bytes */
	fputc( (datasize & 0x000000ff),       fp );
	fputc( (datasize & 0x0000ff00) >> 8,  fp );
	fputc( (datasize & 0x00ff0000) >> 16, fp );
	fputc( (datasize & 0xff000000) >> 24, fp );

	return 0;
}

static int output_au_hdr( FILE *fp, INT4 samplerate, UINT4 datasize )
{
	/*
	 *
	 * Write header
	 *
	 */

	/* magic number */
	fprintf( fp, ".snd" );

	/* header length */
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 24, fp );

	/* data size (bytes) */
	fputc( (datasize & 0xff000000) >> 24, fp );
	fputc( (datasize & 0x00ff0000) >> 16, fp );
	fputc( (datasize & 0x0000ff00) >> 8,  fp );
	fputc( (datasize & 0x000000ff),       fp );

	/* encoding (16 bit pcm) */
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 3, fp );

	/* sample rate */
	fputc( (samplerate & 0xff000000) >> 24, fp );
	fputc( (samplerate & 0x00ff0000) >> 16, fp );
	fputc( (samplerate & 0x0000ff00) >> 8,  fp );
	fputc( (samplerate & 0x000000ff),       fp );

	/* number of channels */
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 0, fp );
	fputc( 1, fp );

	return 0;
}

static int output_REAL4Vector( FILE *fp, REAL4Vector *data, int wavfmt )
{
	REAL4 maxval, minval, midval, scale;
	UINT4 i;

	/* compute range */
	maxval = minval = data->data[0];
	for ( i = 1; i < data->length; ++i ) {
		if ( data->data[i] > maxval )
			maxval = data->data[i];
		if ( data->data[i] < minval )
			minval = data->data[i];
	}
	midval  = 0.5*(maxval + minval);
	maxval -= midval;
	minval -= midval;
	if ( fabs( minval ) > fabs( maxval ) )
		maxval = fabs( minval );
	scale = LAL_SOUND_MAX / maxval;

	/* output data */
	if ( wavfmt == 1 ) { /* wav format */
		for ( i = 0; i < data->length; ++i ) {
			UINT2 val;
			val = (UINT2)( scale * (data->data[i] - midval) );
			fputc( (val & 0x00ff),      fp );
			fputc( (val & 0xff00) >> 8, fp );
		}
	} else {  /* au format */
		for ( i = 0; i < data->length; ++i ) {
			UINT2 val;
			val = (UINT2)(scale*data->data[i]);
			fputc( (val & 0xff00) >> 8, fp );
			fputc( (val & 0x00ff),      fp );
		}
	}
	return 0;
}

static int output_REAL8Vector( FILE *fp, REAL8Vector *data, int wavfmt )
{
	REAL4 maxval, minval, midval, scale;
	UINT4 i;

	/* compute range */
	maxval = minval = data->data[0];
	for ( i = 1; i < data->length; ++i ) {
		if ( data->data[i] > maxval )
			maxval = data->data[i];
		if ( data->data[i] < minval )
			minval = data->data[i];
	}
	midval  = 0.5*(maxval + minval);
	maxval -= midval;
	minval -= midval;
	if ( fabs( minval ) > fabs( maxval ) )
		maxval = fabs( minval );
	scale = LAL_SOUND_MAX / maxval;

	/* output data */
	if ( wavfmt == 1 ) { /* wav format */
		for ( i = 0; i < data->length; ++i ) {
			UINT2 val;
			val = (UINT2)( scale * (data->data[i] - midval) );
			fputc( (val & 0x00ff),      fp );
			fputc( (val & 0xff00) >> 8, fp );
		}
	} else {  /* au format */
		for ( i = 0; i < data->length; ++i ) {
			UINT2 val;
			val = (UINT2)(scale*data->data[i]);
			fputc( (val & 0xff00) >> 8, fp );
			fputc( (val & 0x00ff),      fp );
		}
	}

	return 0;
}

/** Records a time series as a .wav audio file */
int XLALAudioWAVRecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series )
{
	static const char * func = "XLALAudioWAVRecordREAL4TimeSeries";
	INT4  samplerate;
	UINT4 datasize;
	FILE *fpout;

	if ( ! series )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! series->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( ! series->data->length )
		XLAL_ERROR( func, XLAL_EBADLEN );

	datasize  = series->data->length * sizeof( UINT2 );
	samplerate = (INT4)(1.0/series->deltaT);
	if ( samplerate < 1 )
		XLAL_ERROR( func, XLAL_EINVAL );

	fpout = fopen_if_null( fp, series->name, "wav" );
	if ( ! fpout )
		XLAL_ERROR( func, XLAL_EIO );

	/* write header */
	output_wav_hdr( fpout, samplerate, datasize );

	/* write data */
	output_REAL4Vector( fpout, series->data, 1 );

	if ( !fp )
		fclose( fpout );

	return 0;
}

/** Records a time series as a .wav audio file */
int XLALAudioWAVRecordREAL8TimeSeries( FILE *fp, REAL8TimeSeries *series )
{
	static const char * func = "XLALAudioWAVRecordREAL8TimeSeries";
	INT4  samplerate;
	UINT4 datasize;
	FILE *fpout;

	if ( ! series )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! series->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( ! series->data->length )
		XLAL_ERROR( func, XLAL_EBADLEN );

	datasize  = series->data->length * sizeof( UINT2 );
	samplerate = (INT4)(1.0/series->deltaT);
	if ( samplerate < 1 )
		XLAL_ERROR( func, XLAL_EINVAL );

	fpout = fopen_if_null( fp, series->name, "wav" );
	if ( ! fpout )
		XLAL_ERROR( func, XLAL_EIO );

	/* write header */
	output_wav_hdr( fpout, samplerate, datasize );

	/* write data */
	output_REAL8Vector( fpout, series->data, 1 );

	if ( !fp )
		fclose( fpout );

	return 0;
}


/** Records a time series as a .au audio file */
int XLALAudioAURecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series )
{
	static const char * func = "XLALAudioAURecordREAL4TimeSeries";
	INT4  samplerate;
	UINT4 datasize;
	FILE *fpout;

	if ( ! series )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! series->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( ! series->data->length )
		XLAL_ERROR( func, XLAL_EBADLEN );

	datasize  = series->data->length * sizeof( UINT2 );
	samplerate = (INT4)(1.0/series->deltaT);
	if ( samplerate < 1 )
		XLAL_ERROR( func, XLAL_EINVAL );

	fpout = fopen_if_null( fp, series->name, "au" );
	if ( ! fpout )
		XLAL_ERROR( func, XLAL_EIO );

	/* write header */
	output_au_hdr( fpout, samplerate, datasize );

	/* write data */
	output_REAL4Vector( fpout, series->data, 0 );

	if ( !fp )
		fclose( fpout );

	return 0;
}

/** Records a time series as a .au audio file */
int XLALAudioAURecordREAL8TimeSeries( FILE *fp, REAL8TimeSeries *series )
{
	static const char * func = "XLALAudioAURecordREAL8TimeSeries";
	INT4  samplerate;
	UINT4 datasize;
	FILE *fpout;

	if ( ! series )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! series->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( ! series->data->length )
		XLAL_ERROR( func, XLAL_EBADLEN );

	datasize  = series->data->length * sizeof( UINT2 );
	samplerate = (INT4)(1.0/series->deltaT);
	if ( samplerate < 1 )
		XLAL_ERROR( func, XLAL_EINVAL );

	fpout = fopen_if_null( fp, series->name, "au" );
	if ( ! fpout )
		XLAL_ERROR( func, XLAL_EIO );

	/* write header */
	output_au_hdr( fpout, samplerate, datasize );

	/* write data */
	output_REAL8Vector( fpout, series->data, 0 );

	if ( !fp )
		fclose( fpout );

	return 0;
}
