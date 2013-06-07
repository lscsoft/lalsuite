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
