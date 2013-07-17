/*
*  Copyright (C) 2013 Jolien Creighton
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

#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,FrequencySeries)
#define READSERIES CONCAT2(XLALFrFileRead,STYPE)
#define STREAMGETSERIES CONCAT2(XLALFrStreamGet,STYPE)
#define STREAMREADSERIES CONCAT2(XLALFrStreamRead,STYPE)

int STREAMGETSERIES(STYPE * series, LALFrStream * stream)
{
    STYPE *tmpser;

    /* seek to the relevant point in the stream */
    if (XLALFrStreamSeek(stream, &series->epoch))
        XLAL_ERROR(XLAL_EFUNC);

    tmpser = READSERIES(stream->file, series->name, stream->pos);
    if (!tmpser)
        XLAL_ERROR(XLAL_EFUNC);

    *series = *tmpser;
    LALFree(tmpser);
    return 0;
}

STYPE *STREAMREADSERIES(LALFrStream * stream, const char *chname,
    const LIGOTimeGPS * epoch)
{
    STYPE *series;

    /* seek to the relevant point in the stream */
    if (XLALFrStreamSeek(stream, epoch))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    series = READSERIES(stream->file, chname, stream->pos);
    if (!series)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return series;
}

#undef STYPE
#undef READSERIES
#undef STREAMREADSERIES

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3
#undef STRING
