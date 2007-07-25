#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
REAL8Sequence *XLALMoveREAL8Sequence( REAL8Sequence *destination, const REAL8Sequence *source, int first );
REAL8TimeSeries *XLALMoveREAL8TimeSeries( REAL8TimeSeries *destination, const REAL8TimeSeries *source, int first );
REAL8TimeSeries * XLALLeonorWhitenREAL8TimeSeries( REAL8TimeSeries *unwhitened, REAL8FFTPlan *fwdplan, REAL8FFTPlan *revplan, UINT4 seglen, COMPLEX8FrequencySeries *response );
