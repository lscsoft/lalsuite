#include <lal/LALDatatypes.h>
REAL8 XLALPearsonCrossCorrelation( REAL8 *x, REAL8 *y, UINT4 n );
REAL8 XLALPearsonCrossCorrelationSlow( REAL8 *x, REAL8 *y, UINT4 n );
REAL8Vector * XLALPearsonCrossCorrelationREAL8Vector( REAL8Vector *x, REAL8Vector *y, UINT4 seglen );
REAL8 XLALPearsonMaxCrossCorrelationREAL8Vector( REAL8Vector *x, REAL8Vector *y, UINT4 seglen );
