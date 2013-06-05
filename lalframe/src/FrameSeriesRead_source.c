#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,TimeSeries)
#define FSTYPE CONCAT2(TYPE,FrequencySeries)

#define XFUNC CONCAT2(XLALFrStreamGet,STYPE)
#define XFUNCM CONCAT3(XLALFrStreamGet,STYPE,Metadata)
#define FUNC CONCAT2(LALFrGet,STYPE)
#define FUNCM CONCAT3(LALFrGet,STYPE,Metadata)
#define XFSFUNC CONCAT2(XLALFrStreamGet,FSTYPE)
#define FSFUNC CONCAT2(LALFrGet,FSTYPE)
#define XRFUNC CONCAT2(XLALFrStreamRead,STYPE)
#define XCFUNC CONCAT2(XLALCreate,STYPE)
#define XDFUNC CONCAT2(XLALDestroy,STYPE)
#define XREFUNC CONCAT2(XLALResize,STYPE)



int XFSFUNC(FSTYPE * series, LALFrStream * stream)
{
    struct FrVect *vect;

    if (stream->state & LAL_FR_STREAM_ERR)
        XLAL_ERROR(XLAL_EIO);
    if (stream->state & LAL_FR_STREAM_END)
        XLAL_ERROR(XLAL_EIO);

    vect = loadFrVect(stream, series->name);
    if (!vect || !vect->data)
        XLAL_ERROR(XLAL_ENAME); /* couldn't find channel */
    if (vect->type != FRTYPE)
        XLAL_ERROR(XLAL_ETYPE); /* data has wrong type */

#if defined FR_VERS && FR_VERS >= 5000
    series->epoch.gpsSeconds = floor(vect->GTime);
    series->epoch.gpsNanoSeconds =
        floor(1e9 * (vect->GTime - floor(vect->GTime)));
#else
    series->epoch.gpsSeconds = vect->GTimeS;
    series->epoch.gpsNanoSeconds = vect->GTimeN;
#endif
    series->deltaF = vect->dx[0];
    series->f0 = 0;     /* FIXME: should get correct value... */
    series->sampleUnits = lalADCCountUnit;
    series->data = LALCalloc(1, sizeof(*series->data));
    if (!series->data) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_ENOMEM);
    }
    series->data->length = vect->nData;
    series->data->data =
        LALMalloc(series->data->length * sizeof(*series->data->data));
    if (!series->data->data) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_ENOMEM);
    }
    memcpy(series->data->data, vect->FRDATA,
           series->data->length * sizeof(*series->data->data));

    FrVectFree(vect);

    return 0;
}



void
FSFUNC(LALStatus * status,
       FSTYPE * series, FrChanIn * chanin, LALFrStream * stream)
{
    struct FrVect *vect;
    INITSTATUS(status);

    ASSERT(series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(!series->data, status, FRAMESTREAMH_ENNUL,
           FRAMESTREAMH_MSGENNUL);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);

    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }
    if (stream->state & LAL_FR_STREAM_END) {
        ABORT(status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE);
    }

    strncpy(series->name, chanin->name, sizeof(series->name));
    vect = loadFrVect(stream, series->name);
    if (!vect || !vect->data) {
        ABORT(status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN);
    }
    if (vect->type != FRTYPE) {
        ABORT(status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE);
    }
#if defined FR_VERS && FR_VERS >= 5000
    series->epoch.gpsSeconds = floor(vect->GTime);
    series->epoch.gpsNanoSeconds =
        floor(1e9 * (vect->GTime - floor(vect->GTime)));
#else
    series->epoch.gpsSeconds = vect->GTimeS;
    series->epoch.gpsNanoSeconds = vect->GTimeN;
#endif
    series->deltaF = vect->dx[0];
    series->f0 = 0;     /* FIXME: should get correct value... */
    series->sampleUnits = lalADCCountUnit;
    series->data = LALCalloc(1, sizeof(*series->data));
    if (!series->data) {
        ABORT(status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC);
    }
    series->data->length = vect->nData;
    series->data->data =
        LALMalloc(series->data->length * sizeof(*series->data->data));
    if (!series->data->data) {
        ABORT(status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC);
    }
    memcpy(series->data->data, vect->FRDATA,
           series->data->length * sizeof(*series->data->data));

    FrVectFree(vect);

    vect = NULL;

    RETURN(status);
}


int XFUNCM(STYPE * series, LALFrStream * stream)
{
    const REAL8 fuzz = 0.1 / 16384.0;   /* smallest discernable unit of time */
    struct FrVect *vect;
    UINT4 noff;
    INT8 tnow;
    INT8 tbeg;
    REAL8 rate;

    if (stream->state & LAL_FR_STREAM_ERR)
        XLAL_ERROR(XLAL_EIO);
    if (stream->state & LAL_FR_STREAM_END)
        XLAL_ERROR(XLAL_EIO);

    vect = loadFrVect(stream, series->name);
    if (!vect || !vect->data)
        XLAL_ERROR(XLAL_ENAME); /* couldn't find channel */
    if (vect->type != FRTYPE)
        XLAL_ERROR(XLAL_ETYPE); /* data has wrong type */

    tnow = EPOCH_TO_I8TIME(stream->epoch);
#if defined FR_VERS && FR_VERS >= 5000
    tbeg = 1e9 * vect->GTime;
#else
    tbeg = SECNAN_TO_I8TIME(vect->GTimeS, vect->GTimeN);
#endif
    if (tnow + 1000 < tbeg) {   /* added 1000 ns to account for double precision */
        FrVectFree(vect);
        XLAL_ERROR(XLAL_ETIME); /* invalid time offset */
    }

    /* compute number of points offset very carefully:
     * if current time is within fuzz of a sample, get that sample;
     * otherwise get the sample just after the requested time */
    rate = vect->dx[0] ? 1.0 / vect->dx[0] : 0.0;
    noff = ceil((1e-9 * (tnow - tbeg) - fuzz) * rate);

    /* adjust current time to be exactly the first sample
     * (rounded to nearest nanosecond) */
    tnow = tbeg + floor(1e9 * noff * vect->dx[0] + 0.5);

    SET_EPOCH(&series->epoch, tnow);
    series->deltaT = vect->dx[0];
    series->sampleUnits = lalADCCountUnit;

    FrVectFree(vect);
    return 0;
}


int XFUNC(STYPE * series, LALFrStream * stream)
{
    const REAL8 fuzz = 0.1 / 16384.0;   /* smallest discernable unit of time */
    struct FrVect *vect;
    UINT4 need;
    UINT4 noff;
    UINT4 mult;
    UINT4 ncpy;
    TYPE *dest;
    INT8 tnow;
    INT8 tbeg;
    INT8 tend;
    REAL8 rate;
    INT4 gap = 0;

    if (stream->state & LAL_FR_STREAM_ERR)
        XLAL_ERROR(XLAL_EIO);
    if (stream->state & LAL_FR_STREAM_END)
        XLAL_ERROR(XLAL_EIO);

    vect = loadFrVect(stream, series->name);
    if (!vect || !vect->data)
        XLAL_ERROR(XLAL_ENAME); /* couldn't find channel */
    if (vect->type != FRTYPE)
        XLAL_ERROR(XLAL_ETYPE); /* data has wrong type */

    tnow = EPOCH_TO_I8TIME(stream->epoch);
#if defined FR_VERS && FR_VERS >= 5000
    tbeg = 1e9 * vect->GTime;
#else
    tbeg = SECNAN_TO_I8TIME(vect->GTimeS, vect->GTimeN);
#endif
    if (tnow + 1000 < tbeg) {   /* added 1000 ns to account for double precision */
        FrVectFree(vect);
        XLAL_ERROR(XLAL_ETIME); /* invalid time offset */
    }

    /* compute number of points offset very carefully:
     * if current time is within fuzz of a sample, get that sample;
     * otherwise get the sample just after the requested time */
    rate = vect->dx[0] ? 1.0 / vect->dx[0] : 0.0;
    noff = ceil((1e-9 * (tnow - tbeg) - fuzz) * rate);

    /* adjust current time to be exactly the first sample
     * (rounded to nearest nanosecond) */
    tnow = tbeg + floor(1e9 * noff * vect->dx[0] + 0.5);

    SET_EPOCH(&series->epoch, tnow);
    series->deltaT = vect->dx[0];
    series->sampleUnits = lalADCCountUnit;

    if (!series->data) {        /* end here... just metadata was requested */
        FrVectFree(vect);
        return 0;
    }

    /* check to see if data vector is ok */
    if (!series->data->data)
        XLAL_ERROR(XLAL_EFAULT);
    if (!series->data->length)
        XLAL_ERROR(XLAL_EBADLEN);

    /* mult is two if output series is complex */
    mult = sizeof(*series->data->data) / sizeof(*vect->FRDATA);
    dest = series->data->data;
    need = series->data->length;
    if (noff > vect->nData) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_ETIME); /* invalid time offset */
    }

    /* number of points to copy */
    ncpy = (vect->nData - noff < need) ? (vect->nData - noff) : need;
    memcpy(dest, vect->FRDATA + noff * mult,
           ncpy * sizeof(*series->data->data));

    FrVectFree(vect);
    vect = NULL;

    dest += ncpy;
    need -= ncpy;

    /* if still data remaining */
    while (need) {
        if (XLALFrStreamNext(stream) < 0) {
            if (vect)
                FrVectFree(vect);
            memset(dest, 0, need * sizeof(*series->data->data));
            XLAL_ERROR(XLAL_EFUNC);
        }
        if (stream->state & LAL_FR_STREAM_END) {
            if (vect)
                FrVectFree(vect);
            memset(dest, 0, need * sizeof(*series->data->data));
            XLAL_ERROR(XLAL_EIO);
        }

        /* load more data */
        vect = loadFrVect(stream, series->name);
        if (!vect || !vect->data) {
            memset(dest, 0, need * sizeof(*series->data->data));
            XLAL_ERROR(XLAL_ENAME);     /* now channel is missing ... */
        }
        if (vect->type != FRTYPE) {
            FrVectFree(vect);
            memset(dest, 0, need * sizeof(*series->data->data));
            XLAL_ERROR(XLAL_ETYPE);     /* now type is wrong ... */
        }

        if (stream->state & LAL_FR_STREAM_GAP) {        /* gap in data */
            dest = series->data->data;
            need = series->data->length;
#if defined FR_VERS && FR_VERS >= 5000
            series->epoch.gpsSeconds = floor(vect->GTime);
            series->epoch.gpsNanoSeconds =
                (INT8) (1e9 * vect->GTime) % (INT8) 1000000000;
#else
            series->epoch.gpsSeconds = vect->GTimeS;
            series->epoch.gpsNanoSeconds = vect->GTimeN;
#endif
            gap = 1;    /* need to record this because next FrNext will erase! */
        }

        /* copy data */
        ncpy = vect->nData < need ? vect->nData : need;
        memcpy(dest, vect->FRDATA, ncpy * sizeof(*series->data->data));


        FrVectFree(vect);
        vect = NULL;


        dest += ncpy;
        need -= ncpy;
    }

    /* update stream start time very carefully:
     * start time must be the exact time of the next sample, rounded to the
     * nearest nanosecond */
    SET_EPOCH(&stream->epoch, EPOCH_TO_I8TIME(series->epoch)
              + (INT8) floor(1e9 * series->data->length * series->deltaT +
                             0.5));

    /* check to see that we are still within current frame */
    tend = SECNAN_TO_I8TIME(stream->file->toc->GTimeS[stream->pos],
                            stream->file->toc->GTimeN[stream->pos]);
    tend += (INT8) floor(1e9 * stream->file->toc->dt[stream->pos]);
    if (tend <= EPOCH_TO_I8TIME(stream->epoch)) {
        int keepmode = stream->mode;
        LIGOTimeGPS keep;
        keep = stream->epoch;
        /* advance a frame */
        /* failure is benign, so we return results */
        stream->mode |= LAL_FR_STREAM_IGNOREGAP_MODE;
        if (XLALFrStreamNext(stream) < 0) {
            stream->mode = keepmode;
            XLAL_ERROR(XLAL_EFUNC);
        }
        if (!(stream->state & LAL_FR_STREAM_GAP))
            stream->epoch = keep;
        stream->mode = keepmode;
    }

    if (gap)    /* there was a gap in the data */
        stream->state = (FrState) (stream->state | LAL_FR_STREAM_GAP);
    /* FIXME: does this need to cause a failure if mode is set to fail on gaps? */

    if (stream->state & LAL_FR_STREAM_ERR)
        XLAL_ERROR(XLAL_EIO);

    return 0;
}



void
FUNC(LALStatus * status,
     STYPE * series, FrChanIn * chanin, LALFrStream * stream)
{
    const REAL8 fuzz = 0.1 / 16384.0;   /* smallest discernable unit of time */
    struct FrVect *vect;
    UINT4 need;
    UINT4 noff;
    UINT4 mult;
    UINT4 ncpy;
    TYPE *dest;
    INT8 tnow;
    INT8 tbeg;
    INT8 tend;
    REAL8 rate;
    INT4 gap = 0;

    INITSTATUS(status);

    ASSERT(series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);

    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }
    if (stream->state & LAL_FR_STREAM_END) {
        ABORT(status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE);
    }

    strncpy(series->name, chanin->name, sizeof(series->name));
    vect = loadFrVect(stream, series->name);
    if (!vect || !vect->data) {
        ABORT(status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN);
    }
    if (vect->type != FRTYPE) {
        FrVectFree(vect);
        ABORT(status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE);
    }

    tnow = EPOCH_TO_I8TIME(stream->epoch);
#if defined FR_VERS && FR_VERS >= 5000
    tbeg = 1e9 * vect->GTime;
#else
    tbeg = SECNAN_TO_I8TIME(vect->GTimeS, vect->GTimeN);
#endif
    if (tnow + 1000 < tbeg) {   /* added 1000 ns to account for double precision */
        FrVectFree(vect);
        ABORT(status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME);
    }

    /* compute number of points offset very carefully:
     * if current time is within fuzz of a sample, get that sample;
     * otherwise get the sample just after the requested time */
    rate = vect->dx[0] ? 1.0 / vect->dx[0] : 0.0;
    noff = ceil((1e-9 * (tnow - tbeg) - fuzz) * rate);

    /* adjust current time to be exactly the first sample
     * (rounded to nearest nanosecond) */
    tnow = tbeg + floor(1e9 * noff * vect->dx[0] + 0.5);


    SET_EPOCH(&series->epoch, tnow);
    series->deltaT = vect->dx[0];
    series->sampleUnits = lalADCCountUnit;

    if (!series->data) {        /* no data requested: return now */
        FrVectFree(vect);
        RETURN(status);
    }
    ASSERT(series->data->data, status, FRAMESTREAMH_ENULL,
           FRAMESTREAMH_MSGENULL);
    ASSERT(series->data->length > 0, status, FRAMESTREAMH_ESIZE,
           FRAMESTREAMH_MSGESIZE);

    ATTATCHSTATUSPTR(status);

    /* mult is two if output series is complex */
    mult = sizeof(*series->data->data) / sizeof(*vect->FRDATA);
    dest = series->data->data;
    need = series->data->length;
    if (noff > vect->nData) {
        FrVectFree(vect);
        ABORT(status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME);
    }

    /* number of points to copy */
    ncpy = (vect->nData - noff < need) ? (vect->nData - noff) : need;
    memcpy(dest, vect->FRDATA + noff * mult,
           ncpy * sizeof(*series->data->data));

    FrVectFree(vect);
    vect = NULL;

    dest += ncpy;
    need -= ncpy;


    /* if still data remaining */
    while (need) {
        LALFrNext(status->statusPtr, stream);
        BEGINFAIL(status) {
            if (vect)
                FrVectFree(vect);
            memset(dest, 0, need * sizeof(*series->data->data));
        }
        ENDFAIL(status);
        if (stream->state & LAL_FR_STREAM_END) {
            if (vect)
                FrVectFree(vect);
            memset(dest, 0, need * sizeof(*series->data->data));
            ABORT(status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE);
        }

        /* load more data */
        vect = loadFrVect(stream, series->name);
        if (!vect || !vect->data) {
            memset(dest, 0, need * sizeof(*series->data->data));
            ABORT(status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN);
        }
        if (vect->type != FRTYPE) {
            FrVectFree(vect);
            memset(dest, 0, need * sizeof(*series->data->data));
            ABORT(status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE);
        }

        if (stream->state & LAL_FR_STREAM_GAP) {        /* gap in data */
            dest = series->data->data;
            need = series->data->length;
#if defined FR_VERS && FR_VERS >= 5000
            series->epoch.gpsSeconds = floor(vect->GTime);
            series->epoch.gpsNanoSeconds =
                (INT8) (1e9 * vect->GTime) % (INT8) 1000000000;
#else
            series->epoch.gpsSeconds = vect->GTimeS;
            series->epoch.gpsNanoSeconds = vect->GTimeN;
#endif
            gap = 1;    /* need to record this because next FrNext will erase! */
        }

        /* copy data */
        ncpy = vect->nData < need ? vect->nData : need;
        memcpy(dest, vect->FRDATA, ncpy * sizeof(*series->data->data));


        FrVectFree(vect);
        vect = NULL;


        dest += ncpy;
        need -= ncpy;
    }

    /* update stream start time very carefully:
     * start time must be the exact time of the next sample, rounded to the
     * nearest nanosecond */
    SET_EPOCH(&stream->epoch, EPOCH_TO_I8TIME(series->epoch)
              + (INT8) floor(1e9 * series->data->length * series->deltaT +
                             0.5));

    /* check to see that we are still within current frame */
    tend =
        SECNAN_TO_I8TIME(stream->file->toc->GTimeS[stream->pos],
                         stream->file->toc->GTimeN[stream->pos]);
    tend += (INT8) floor(1e9 * stream->file->toc->dt[stream->pos]);
    if (tend <= EPOCH_TO_I8TIME(stream->epoch)) {
        int keepmode = stream->mode;
        LIGOTimeGPS keep;
        keep = stream->epoch;
        /* advance a frame */
        /* failure is benign, so we return results */
        stream->mode |= LAL_FR_STREAM_IGNOREGAP_MODE;
        TRY(LALFrNext(status->statusPtr, stream), status);
        if (!(stream->state & LAL_FR_STREAM_GAP)) {
            stream->epoch = keep;
        }
        stream->mode = keepmode;
    }

    if (gap) {  /* there was a gap in the data */
        stream->state = (FrState) (stream->state | LAL_FR_STREAM_GAP);
    }

    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }

    /* remove this: the error will be reported on the *next* call! */
    /*
       if ( stream->state & LAL_FR_STREAM_END )
       {
       ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
       }
     */

    DETATCHSTATUSPTR(status);
    RETURN(status);
}



void
FUNCM(LALStatus * status,
      STYPE * series, FrChanIn * chanin, LALFrStream * stream)
{
    void *sequence;

    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    ASSERT(series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);

    /* save the sequence address, then wipe the series structure */
    sequence = series->data;
    memset(series, 0, sizeof(*series));

    /* call FUNC to populate the series' metadata */
    FUNC(status->statusPtr, series, chanin, stream);
    CHECKSTATUSPTR(status);

    /* restore the sequence address */
    series->data = sequence;

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


STYPE *XRFUNC(LALFrStream * stream,
              const char *chname,
              const LIGOTimeGPS * start,
              REAL8 duration, size_t lengthlimit)
{
    STYPE *series;
    size_t length;

    /* create and initialize a zero-length time series vector */
    series = XCFUNC(chname, start, 0.0, 0.0, &lalADCCountUnit, 0);
    if (!series)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    /* get the time series meta-data */
    if (XFUNCM(series, stream)) {
        XDFUNC(series);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* resize the time series to the correct number of samples */
    length = duration / series->deltaT;
    if (lengthlimit && (lengthlimit < length))
        length = lengthlimit;
    if (!XREFUNC(series, 0, length)) {
        XDFUNC(series);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* read the data */
    if (XLALFrStreamSeek(stream, start) || XFUNC(series, stream)) {
        XDFUNC(series);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    return (series);
}
