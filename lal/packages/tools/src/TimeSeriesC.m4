dnl $Id$
dnl
dnl Copyright (C) 2007  Kipp Cannon
dnl
dnl This program is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by the
dnl Free Software Foundation; either version 2 of the License, or (at your
dnl option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
dnl Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License along
dnl with this program; if not, write to the Free Software Foundation, Inc.,
dnl 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

define(`SERIESTYPE',DATATYPE`TimeSeries')
define(`SEQUENCETYPE',DATATYPE`Sequence')
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
)
{
	if(series)
		`XLALDestroy'SEQUENCETYPE (series->data);
	XLALFree(series);
}


SERIESTYPE *`XLALCreate'SERIESTYPE (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
)
{
	static const char func[] = "`XLALCreate'SERIESTYPE";
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	new = XLALMalloc(sizeof(*new));
	sequence = `XLALCreate'SEQUENCETYPE (length);
	if(!new || !sequence) {
		XLALFree(new);
		`XLALDestroy'SEQUENCETYPE (sequence);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	if(name)
		strncpy(new->name, name, LALNameLength);
	else
		new->name[0] = '\0';
	new->epoch = *epoch;
	new->f0 = f0;
	new->deltaT = deltaT;
	new->sampleUnits = *sampleUnits;
	new->data = sequence;

	return new;
}


SERIESTYPE *`XLALCut'SERIESTYPE (
	const SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	static const char func[] = "`XLALCut'SERIESTYPE";
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	new = XLALMalloc(sizeof(*new));
	sequence = `XLALCut'SEQUENCETYPE (series->data, first, length);
	if(!new || !sequence) {
		XLALFree(new);
		`XLALDestroy'SEQUENCETYPE (sequence);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	*new = *series;
	new->data = sequence;
	XLALGPSAdd(&new->epoch, first * new->deltaT);

	return new;
}


SERIESTYPE *`XLALResize'SERIESTYPE (
	SERIESTYPE *series,
	int first,
	size_t length
)
{
	static const char func[] = "`XLALResize'SERIESTYPE";

	if(!`XLALResize'SEQUENCETYPE (series->data, first, length))
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	XLALGPSAdd(&series->epoch, first * series->deltaT);

	return series;
}


SERIESTYPE *`XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	static const char func[] = "`XLALShrink'SERIESTYPE";

	if(!`XLALResize'SERIESTYPE (series, first, length))
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	return series;
}


SERIESTYPE *`XLALAdd'SERIESTYPE (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
)
{
	static const char func[] = "`XLALAdd'SERIESTYPE";
	REAL8 Delta_epoch = XLALGPSDiff(&arg2->epoch, &arg1->epoch);
	/* number of arg1 units per arg2 unit.  XLALUnitRatio() returns the
	 * number one obtains when one divides 1 of the first argument by 1
	 * of the second argument, for example if arg2 is in m and arg1 is
	 * in cm then unit_ratio = 100.0 */
	REAL8 unit_ratio = XLALUnitRatio(&arg2->sampleUnits, &arg1->sampleUnits);
	unsigned i, j;

	/* make sure arguments are compatible */
	if(XLALIsREAL8FailNaN(unit_ratio)) {
		XLALPrintError("%s(): incompatible units\n", func);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	if(arg1->f0 != arg2->f0) {
		XLALPrintError("%s(): incompatible heterodyne frequencies\n", func);
		XLAL_ERROR_NULL(func, XLAL_EDATA);
	}
	if(arg1->deltaT != arg2->deltaT) {
		XLALPrintError("%s(): incompatible sample periods\n", func);
		XLAL_ERROR_NULL(func, XLAL_EDATA);
	}

	/* set start indexes */
	if(Delta_epoch >= 0) {
		i = floor(Delta_epoch / arg1->deltaT + 0.5);
		j = 0;
	} else {
		i = 0;
		j = floor(-Delta_epoch / arg2->deltaT + 0.5);
	}

	/* add arg2 to arg1, adjusting the units */
	for(; i < arg1->data->length && j < arg2->data->length; i++, j++) {
		ifelse(DATATYPE, COMPLEX8,
		arg1->data->data[i].re += arg2->data->data[j].re * unit_ratio;
		arg1->data->data[i].im += arg2->data->data[j].im * unit_ratio;
		, DATATYPE, COMPLEX16,
		arg1->data->data[i].re += arg2->data->data[j].re * unit_ratio;
		arg1->data->data[i].im += arg2->data->data[j].im * unit_ratio;
		, 
		arg1->data->data[i] += arg2->data->data[j] * unit_ratio;)
	}

	return arg1;
}

