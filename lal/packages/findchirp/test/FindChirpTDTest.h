#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>

int Start(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    FindChirpInitParams     *initParams
    );

int Stop(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    UINT4 numChisqBins
    );

int SPInit(
    FindChirpTmpltParams **spTmpltParams,
    FindChirpDataParams  **spDataParams,
    FindChirpInitParams     *initParams,
    REAL4 srate,
    REAL4 fmin,
    REAL4 dynRange,
    UINT4 trunc
    );

int SPFini(
    FindChirpTmpltParams **spTmpltParams,
    FindChirpDataParams  **spDataParams
    );

int TDInit(
    FindChirpDataParams **tdDataParams,
    FindChirpInitParams    *initParams,
    REAL4 srate,
    REAL4 fmin,
    REAL4 dynRange,
    UINT4 trunc
    );

int TDFini(
    FindChirpDataParams  **tdDataParams
    );

int SPFilter(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    FindChirpFilterInput *filterInput,
    FindChirpFilterParams *filterParams,
    FindChirpSegmentVector *fcSegVec,
    FindChirpTmpltParams *tmpltParams,
    FindChirpDataParams *dataParams
    );

int TDFilter(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 fmax,
    FindChirpFilterInput *filterInput,
    FindChirpFilterParams *filterParams,
    FindChirpSegmentVector *fcSegVec,
    FindChirpDataParams *dataParams
    );

int MakeData(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 srate,
    REAL4 fmin,
    REAL4 fmax
    );
