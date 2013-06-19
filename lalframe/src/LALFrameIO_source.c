#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#if DOM != TDOM && DOM != FDOM
#error time/frequency domain not set
#endif

#define VTYPE CONCAT2(LAL_FRAMEU_FR_VECT_, VEXT)
#if DOM == TDOM
#define STYPE CONCAT2(TYPE,TimeSeries)
#elif DOM == FDOM
#define STYPE CONCAT2(TYPE,FrequencySeries)
#endif

#define CFUNC CONCAT2(XLALCreate,STYPE)
#define RFUNC CONCAT2(XLALFrFileRead,STYPE)
#define MFUNC CONCAT3(XLALFrFileRead,STYPE,Metadata)
#define FUNC_ CONCAT3(XLALFrFileRead,STYPE,_)

static STYPE *FUNC_(LALFrFile * stream, const char *name, size_t pos,
    int load)
{
    STYPE *series;
    LALFrameUFrChan *channel;
#   if 0
    const char *unitX;
#   endif
    const char *unitY;
    LALUnit sampleUnits;
    LIGOTimeGPS epoch;
    double deltaX;
    size_t length;
    size_t bytes;
    void *data;
    int errnum;

    channel = XLALFrameUFrChanRead(stream->file, name, pos);
    if (!channel)
        XLAL_ERROR_NULL(XLAL_ENAME);

    /* make sure it is 1d */
    if (XLALFrameUFrChanVectorQueryNDim(channel) != 1) {
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EDIMS);
    }

    /* check type */
    if (XLALFrameUFrChanVectorQueryType(channel) != VTYPE) {
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }

    /* make sure units are correct and figure out unitY */
#   if 0
    unitX = XLALFrameUFrChanVectorQueryUnitX(channel, 0);
#   endif
    unitY = XLALFrameUFrChanVectorQueryUnitY(channel);
#   if 0        /* don't do unit checking - frames likely have unparsable units */
#   if DOM == TDOM
    if (strcmp(unitX, "s") && strcmp(unitX, "time")) {
        /* doesn't seem to be a tseries */
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EUNIT);
    }
#   elif DOM == FDOM
    if (strcmp(unitX, "s^-1") && strcmp(unitX, "Hz")) {
        /* doesn't seem to be a fseries */
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EUNIT);
    }
#   endif
#   endif
    XLAL_TRY(XLALParseUnitString(&sampleUnits, unitY), errnum);
    if (errnum) {
        XLAL_PRINT_WARNING("Could not parse unit string %s\n", unitY);
        sampleUnits = lalDimensionlessUnit;
    }

    XLALFrFileQueryGTime(&epoch, stream, pos);
    XLALGPSAdd(&epoch, XLALFrameUFrChanQueryTimeOffset(channel));
#   if DOM == TDOM
    XLALGPSAdd(&epoch, XLALFrameUFrChanVectorQueryStartX(channel, 0));
#   endif
    deltaX = XLALFrameUFrChanVectorQueryDx(channel, 0);
    length = XLALFrameUFrChanVectorQueryNData(channel);

    if (!load) {
        /* not expected to load the data vector
         * so exit now with a zero-length vector */
        XLALFrameUFrChanFree(channel);
        series = CFUNC(name, &epoch, 0.0, deltaX, &sampleUnits, 0);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
        return series;
    }

    XLALFrameUFrChanVectorExpand(channel);
    data = XLALFrameUFrChanVectorQueryData(channel);
    if (!data) {
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EDATA);
    }
    bytes = XLALFrameUFrChanVectorQueryNBytes(channel);
    /* make sure bytes, type, and length are sane */
    if (bytes != length * sizeof(TYPE)) {
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EBADLEN);
    }

    series = CFUNC(name, &epoch, 0.0, deltaX, &sampleUnits, length);
    if (!series) {
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }
    memcpy(series->data->data, data, bytes);

    XLALFrameUFrChanFree(channel);
    return series;
}

STYPE *MFUNC(LALFrFile * stream, const char *chname, size_t pos)
{
    return FUNC_(stream, chname, pos, 0);
}

STYPE *RFUNC(LALFrFile * stream, const char *chname, size_t pos)
{
    return FUNC_(stream, chname, pos, 1);
}

#undef FUNC_
#undef MFUNC
#undef RFUNC
#undef CFUNC
#undef STYPE
#undef VTYPE
#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3
#undef STRING
