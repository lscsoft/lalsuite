REAL4 compute_candle_distance(
    REAL4 candleM1, 
    REAL4 candleM2,
    REAL4 snr, 
    REAL8 chanDeltaT, 
    INT4 nPoints, 
    REAL8FrequencySeries *spec, 
    UINT4 cut);

SummValueTable **add_summvalue_table(
    SummValueTable **newTable,
    LIGOTimeGPS gpsStartTime, 
    LIGOTimeGPS gpsEndTime, 
    const CHAR *programName, 
    const CHAR *ifoName, 
    const CHAR *summValueName, 
    const CHAR *comment, 
    REAL8 value
    );

