/*
 * Copyright (C) 2007 Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


/** \file
 *
 * \ingroup LIGOLwXMLArray
 * \author Kipp Cannon
 *
 * \brief I/O library for LIGO Light Weight XML Array document fragments.
 *
 * The following example code constructs a nonsense frequency series and
 * writes it to a LIGO Light Weight XML file.
 *
 * #include <math.h>
 * #include <stdlib.h>
 * #include <lal/FrequencySeries.h>
 * #include <lal/LIGOLwXML.h>
 * #include <lal/LIGOLwXMLArray.h>
 * #include <lal/Units.h>
 *
 * int main()
 * {
 * 	LALStatus stat;
 * 	LIGOLwXMLStream xml;
 * 	REAL8FrequencySeries *series;
 * 	LIGOTimeGPS epoch = {772847291, 500000000};
 * 	int i;
 *
 * 	memset(&stat, 0, sizeof(stat));
 *
 * 	series = XLALCreateREAL8FrequencySeries("PSD", &epoch, 0.0, 1.0 / 16384, &lalVoltUnit, 16384);
 * 	series->data->data[0] = 1.0;
 * 	for(i = 1; i < series->data->length; i++)
 * 		series->data->data[i] = sin(3.14 * i / 1000) / (3.14 * i / 1000);
 *
 * 	LALOpenLIGOLwXMLFile(&stat, &xml, "fseries.xml");
 * 	XLALWriteLIGOLwXMLArrayREAL8FrequencySeries(&xml, "a sin(x)/x function", series);
 * 	LALCloseLIGOLwXMLFile(&stat, &xml);
 *
 * 	return 0;
 * }
 *
 *
 * $Id$
 */


#include <stdarg.h>
#include <stdio.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LIGOLwXML.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <lal/LIGOLwXMLArray.h>

#define FILE LALFILE
#define fprintf XLALFilePrintf


NRCSID(LIGOLWXMLARRAYC, "$Id$");


/*
 * ============================================================================
 *
 *                            Common Writing Code
 *
 * ============================================================================
 */


/*
 * Common code for writing series meta data.
 */


static int WriteLIGOLwXMLArrayMeta(
	FILE *stream,
	const char *indent,
	const LIGOTimeGPS *epoch,
	REAL8 f0
)
{
	int error = 0;
	char *s = XLALGPSToStr(NULL, epoch);

	if(!s)
		return -1;

	error |= fprintf(stream, "%s<Time Name=\"epoch\" Type=\"GPS\">%s</Time>\n", indent, s) < 0;
	error |= fprintf(stream, "%s<Param Name=\"f0\" Unit=\"s^-1\" Type=\"real_8\">%.16g</Param>\n", indent, f0) < 0;

	XLALFree(s);

	return error ? -1 : 0;
}


/*
 * Common code for writing Array header and footer.
 */


static int WriteLIGOLwXMLArrayHeader(
	FILE *stream,
	const char *indent,
	const char *name,
	const char *type,
	const LALUnit *units,
	int length,
	double delta,
	int is_tseries,
	int is_complex
)
{
	const char *dim1 = is_tseries ? "Time" : "Frequency";
	const char *units1 = is_tseries ? "s" : "s^-1";
	const char *dim2 = is_complex ? "Real,Imaginary" : "Real";
	char units2[100];
	int error = 0;

	if(!XLALUnitAsString(units2, 100, units))
		return -1;

	error |= fprintf(stream, "%s<Array Name=\"%s\" Type=\"%s\" Unit=\"%s\">\n", indent, name, type, units2) < 0;
	error |= fprintf(stream, "%s\t<Dim Name=\"%s\" Unit=\"%s\" Start=\"0\" Scale=\"%.16g\">%d</Dim>\n", indent, dim1, units1, delta, length) < 0;
	error |= fprintf(stream, "%s\t<Dim Name=\"%s,%s\">%d</Dim>\n", indent, dim1, dim2, is_complex ? 3 : 2) < 0;

	return error ? -1 : 0;
}


static int WriteLIGOLwXMLArrayFooter(
	FILE *stream,
	const char *indent
)
{
	return fprintf(stream, "%s</Array>\n", indent) < 0 ? -1 : 0;
}


/*
 * Common code for writing stream contents.  The loops are complicated by
 * the need to omit the delimiter at the end of the last line.
 */


static int WriteLIGOLwXMLArrayREAL4Stream(
	FILE *stream,
	const char *indent,
	int length,
	double delta,
	const REAL4 *data
)
{
	int error = 0;
	int i;

	error |= fprintf(stream, "%s<Stream Type=\"Local\" Delimiter=\",\">\n", indent) < 0;
	if(length) {
		for(i = 0; i < length - 1; i++)
			error |= fprintf(stream, "%s\t%.16g,%.8g,\n", indent, i * delta, (double) data[i]) < 0;
		error |= fprintf(stream, "%s\t%.16g,%.8g\n", indent, i * delta, (double) data[i]) < 0;
	}
	error |= fprintf(stream, "%s</Stream>\n", indent) < 0;

	return error ? -1 : 0;
}


static int WriteLIGOLwXMLArrayREAL8Stream(
	FILE *stream,
	const char *indent,
	int length,
	double delta,
	const REAL8 *data
)
{
	int error = 0;
	int i;

	error |= fprintf(stream, "%s<Stream Type=\"Local\" Delimiter=\",\">\n", indent) < 0;
	if(length) {
		for(i = 0; i < length - 1; i++)
			error |= fprintf(stream, "%s\t%.16g,%.16g,\n", indent, i * delta, data[i]) < 0;
		error |= fprintf(stream, "%s\t%.16g,%.16g\n", indent, i * delta, data[i]) < 0;
	}
	error |= fprintf(stream, "%s</Stream>\n", indent) < 0;

	return error ? -1 : 0;
}


static int WriteLIGOLwXMLArrayCOMPLEX8Stream(
	FILE *stream,
	const char *indent,
	int length,
	double delta,
	const COMPLEX8 *data
)
{
	int error = 0;
	int i;

	error |= fprintf(stream, "%s<Stream Type=\"Local\" Delimiter=\",\">\n", indent) < 0;
	if(length) {
		for(i = 0; i < length - 1; i++)
			error |= fprintf(stream, "%s\t%.16g,%.8g,%.8g,\n", indent, i * delta, (double) data[i].re, (double) data[i].im) < 0;
		error |= fprintf(stream, "%s\t%.16g,%.8g,%.8g\n", indent, i * delta, (double) data[i].re, (double) data[i].im) < 0;
	}
	error |= fprintf(stream, "%s</Stream>\n", indent) < 0;

	return error ? -1 : 0;
}


static int WriteLIGOLwXMLArrayCOMPLEX16Stream(
	FILE *stream,
	const char *indent,
	int length,
	double delta,
	const COMPLEX16 *data
)
{
	int error = 0;
	int i;

	error |= fprintf(stream, "%s<Stream Type=\"Local\" Delimiter=\",\">\n", indent) < 0;
	if(length) {
		for(i = 0; i < length - 1; i++)
			error |= fprintf(stream, "%s\t%.16g,%.16g,%.16g,\n", indent, i * delta, data[i].re, data[i].im) < 0;
		error |= fprintf(stream, "%s\t%.16g,%.16g,%.16g\n", indent, i * delta, data[i].re, data[i].im) < 0;
	}
	error |= fprintf(stream, "%s</Stream>\n", indent) < 0;

	return error ? -1 : 0;
}


/*
 * Common code for writing series structure, called by type-specific
 * wrappers.
 */


static int WriteLIGOLwXMLArray(
	FILE *stream,
	const char *comment,
	const char *name,
	const LIGOTimeGPS *epoch,
	const LALUnit *units,
	double f0,
	double delta,
	int is_tseries,
	int is_complex,
	int is_real4,
	int length,
	const void *data
)
{
	int error = 0;

	error |= fprintf(stream, "\t<LIGO_LW Name=\"%s%s\">\n", is_complex ? (is_real4 ? "COMPLEX8" : "COMPLEX16") : (is_real4 ? "REAL4" : "REAL8"), is_tseries ? "TimeSeries" : "FrequencySeries") < 0;
	if(comment)
		error |= fprintf(stream, "\t\t<Comment>%s</Comment>\n", comment) < 0;
	error |= WriteLIGOLwXMLArrayMeta(stream, "\t\t", epoch, f0) < 0;
	error |= WriteLIGOLwXMLArrayHeader(stream, "\t\t", name, is_real4 ? "real_4" : "real_8", units, length, delta, is_tseries, is_complex) < 0;

	if(is_real4) {
		if(is_complex)
			error |= WriteLIGOLwXMLArrayCOMPLEX8Stream(stream, "\t\t\t", length, delta, data) < 0;
		else
			error |= WriteLIGOLwXMLArrayREAL4Stream(stream, "\t\t\t", length, delta, data) < 0;
	} else {
		if(is_complex)
			error |= WriteLIGOLwXMLArrayCOMPLEX16Stream(stream, "\t\t\t", length, delta, data) < 0;
		else
			error |= WriteLIGOLwXMLArrayREAL8Stream(stream, "\t\t\t", length, delta, data) < 0;
	}

	error |= WriteLIGOLwXMLArrayFooter(stream, "\t\t") < 0;
	error |= fprintf(stream, "\t</LIGO_LW>\n") < 0;

	return error ? -1 : 0;
}


/*
 * ============================================================================
 *
 *                     Exported Series Writing Functions
 *
 * ============================================================================
 */


/** Write a \c REAL4TimeSeries to a \c LIGOLwXMLStream.  Returns 0 on
 * success, less than 0 on failure.  If \a comment is not NULL, it will be
 * added to the output as the string in a Comment element.
 */


int XLALWriteLIGOLwXMLArrayREAL4TimeSeries(
	LIGOLwXMLStream *xml,	/*< LIGOLwXMLStream target */
	const char *comment,	/*< Optional comment string */
	const REAL4TimeSeries *series	/*< REAL4TimeSeries to write */
)
{
	const char func[] = "XLALWriteLIGOLwXMLArrayREAL4TimeSeries";

	if(WriteLIGOLwXMLArray(xml->fp, comment, series->name, &series->epoch, &series->sampleUnits, series->f0, series->deltaT, 1, 0, 1, series->data->length, series->data->data) < 0)
		XLAL_ERROR(func, XLAL_EIO);
	return 0;
}


/** Write a \c REAL8TimeSeries to a \c LIGOLwXMLStream.  Returns 0 on
 * success, less than 0 on failure.  If \a comment is not NULL, it will be
 * added to the output as the string in a Comment element.
 */


int XLALWriteLIGOLwXMLArrayREAL8TimeSeries(
	LIGOLwXMLStream *xml,	/*< LIGOLwXMLStream target */
	const char *comment,	/*< Optional comment string */
	const REAL8TimeSeries *series	/*< REAL8TimeSeries to write */
)
{
	const char func[] = "XLALWriteLIGOLwXMLArrayREAL8TimeSeries";

	if(WriteLIGOLwXMLArray(xml->fp, comment, series->name, &series->epoch, &series->sampleUnits, series->f0, series->deltaT, 1, 0, 0, series->data->length, series->data->data) < 0)
		XLAL_ERROR(func, XLAL_EIO);
	return 0;
}


/** Write a \c REAL4FrequencySeries to a \c LIGOLwXMLStream.  Returns 0 on
 * success, less than 0 on failure.  If \a comment is not NULL, it will be
 * added to the output as the string in a Comment element.
 */


int XLALWriteLIGOLwXMLArrayREAL4FrequencySeries(
	LIGOLwXMLStream *xml,	/*< LIGOLwXMLStream target */
	const char *comment,	/*< Optional comment string */
	const REAL4FrequencySeries *series	/*< REAL4TimeSeries to write */
)
{
	const char func[] = "XLALWriteLIGOLwXMLArrayREAL4FrequencySeries";

	if(WriteLIGOLwXMLArray(xml->fp, comment, series->name, &series->epoch, &series->sampleUnits, series->f0, series->deltaF, 0, 0, 1, series->data->length, series->data->data) < 0)
		XLAL_ERROR(func, XLAL_EIO);
	return 0;
}


/** Write a \c REAL8FrequencySeries to a \c LIGOLwXMLStream.  Returns 0 on
 * success, less than 0 on failure.  If \a comment is not NULL, it will be
 * added to the output as the string in a Comment element.
 */


int XLALWriteLIGOLwXMLArrayREAL8FrequencySeries(
	LIGOLwXMLStream *xml,	/*< LIGOLwXMLStream target */
	const char *comment,	/*< Optional comment string */
	const REAL8FrequencySeries *series	/*< REAL8FrequencySeries to write */
)
{
	const char func[] = "XLALWriteLIGOLwXMLArrayREAL8FrequencySeries";

	if(WriteLIGOLwXMLArray(xml->fp, comment, series->name, &series->epoch, &series->sampleUnits, series->f0, series->deltaF, 0, 0, 0, series->data->length, series->data->data) < 0)
		XLAL_ERROR(func, XLAL_EIO);
	return 0;
}


/** Write a \c COMPLEX8FrequencySeries to a \c LIGOLwXMLStream.  Returns 0
 * on success, less than 0 on failure.  If \a comment is not NULL, it will
 * be added to the output as the string in a Comment element.
 */


int XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries(
	LIGOLwXMLStream *xml,	/*< LIGOLwXMLStream target */
	const char *comment,	/*< Optional comment string */
	const COMPLEX8FrequencySeries *series	/*< COMPLEX8FrequencySeries to write */
)
{
	const char func[] = "XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries";

	if(WriteLIGOLwXMLArray(xml->fp, comment, series->name, &series->epoch, &series->sampleUnits, series->f0, series->deltaF, 0, 1, 1, series->data->length, series->data->data) < 0)
		XLAL_ERROR(func, XLAL_EIO);
	return 0;
}


/** Write a \c COMPLEX16FrequencySeries to a \c LIGOLwXMLStream.  Returns 0
 * on success, less than 0 on failure.  If \a comment is not NULL, it will
 * be added to the output as the string in a Comment element.
 */


int XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries(
	LIGOLwXMLStream *xml,	/*< LIGOLwXMLStream target */
	const char *comment,	/*< Optional comment string */
	const COMPLEX16FrequencySeries *series	/*< COMPLEX16FrequencySeries to write */
)
{
	const char func[] = "XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries";

	if(WriteLIGOLwXMLArray(xml->fp, comment, series->name, &series->epoch, &series->sampleUnits, series->f0, series->deltaF, 0, 1, 0, series->data->length, series->data->data) < 0)
		XLAL_ERROR(func, XLAL_EIO);
	return 0;
}

