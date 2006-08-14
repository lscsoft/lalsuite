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

/** Records a time series as a .wav audio file */
int XLALAudioWAVRecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series )
{
	static const char * func = "XLALAudioWAVRecordREAL4TimeSeries";
	REAL4 maxval, minval, midval, scale;
	INT4  samplerate, byterate;
	UINT4 totalsize, datasize;
	UINT4 i;
	int fclose_required = 0;

	if ( ! series )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! series->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( ! series->data->length )
		XLAL_ERROR( func, XLAL_EBADLEN );

	datasize  = series->data->length * sizeof( UINT2 ); 
	totalsize = datasize + 36; /* for header */
	samplerate = (INT4)(1.0/series->deltaT);
	byterate = samplerate * sizeof( UINT2 );
	if ( samplerate < 1 )
		XLAL_ERROR( func, XLAL_EINVAL );

	/* compute range */
	maxval = minval = series->data->data[0];
	for ( i = 1; i < series->data->length; ++i ) {
		if ( series->data->data[i] > maxval )
			maxval = series->data->data[i];
		if ( series->data->data[i] < minval )
			minval = series->data->data[i];
	}
	midval  = 0.5*(maxval + minval);
	maxval -= midval;
	minval -= midval;
	if ( fabs( minval ) > fabs( maxval ) )
		maxval = fabs( minval );
	scale = LAL_SOUND_MAX / maxval;

	/* if file pointer is NULL then create a new file based on series name */
	if ( ! fp ) {
		char fname[FILENAME_MAX];
		LALSnprintf( fname, sizeof( fname ), "%s.wav", series->name );
		fp = fopen( fname, "w" );
		fclose_required = 1;
	}
	if ( ! fp )
		XLAL_ERROR( func, XLAL_EIO );


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

	/* write data */
	for ( i = 0; i < series->data->length; ++i ) {
		UINT2 val;
		val = (UINT2)( scale * (series->data->data[i] - midval) );
		fputc( (val & 0x00ff),      fp );
		fputc( (val & 0xff00) >> 8, fp );
	}

	if ( fclose_required )
		fclose( fp );

	return 0;
}


/** Records a time series as a .au audio file */
int XLALAudioAURecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series )
{
	static const char * func = "XLALAudioAURecordREAL4TimeSeries";
	REAL4 maxval, minval, midval, scale;
	INT4  samplerate;
	UINT4 datasize;
	UINT4 i;
	int fclose_required = 0;

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

	/* compute range */
	maxval = minval = series->data->data[0];
	for ( i = 1; i < series->data->length; ++i ) {
		if ( series->data->data[i] > maxval )
			maxval = series->data->data[i];
		if ( series->data->data[i] < minval )
			minval = series->data->data[i];
	}
	midval  = 0.5*(maxval + minval);
	maxval -= midval;
	minval -= midval;
	if ( fabs( minval ) > fabs( maxval ) )
		maxval = fabs( minval );
	scale = LAL_SOUND_MAX / maxval;

	/* if file pointer is NULL then create a new file based on series name */
	if ( ! fp ) {
		char fname[FILENAME_MAX];
		LALSnprintf( fname, sizeof( fname ), "%s.wav", series->name );
		fp = fopen( fname, "w" );
		fclose_required = 1;
	}
	if ( ! fp )
		XLAL_ERROR( func, XLAL_EIO );

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

	/*
	 *
	 * Write Data
	 *
	 */

	for ( i = 0; i < series->data->length; ++i ) {
		UINT2 val;
		val = (UINT2)(scale*series->data->data[i]);
		fputc( (val & 0xff00) >> 8, fp );
		fputc( (val & 0x00ff),      fp );
	}

	if ( fclose_required )
		fclose( fp );

	return 0;
}
