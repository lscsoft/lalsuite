#include <math.h>
#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include "crosscorr.h"

REAL8 XLALPearsonCrossCorrelation( REAL8 *x, REAL8 *y, UINT4 n )
{
	REAL8 xssq = 0.0;
	REAL8 yssq = 0.0;
	REAL8 xyssq = 0.0;
	REAL8 xave = x[0];
	REAL8 yave = x[0];
	REAL8 xstddev;
	REAL8 ystddev;
	REAL8 xycov;
	REAL8 r;
	UINT4 i;
	for ( i = 0; i < n; ++i ) {
		REAL8 fac = 1.0/(i + 1.0);
		REAL8 sweep = fac * i;
		REAL8 dx = x[i] - xave;
		REAL8 dy = y[i] - yave;
		xssq  += dx * dx * sweep;
		yssq  += dy * dy * sweep;
		xyssq += dx * dy * sweep;
		xave  += fac * dx;
		yave  += fac * dy;
	}
	xstddev = sqrt( xssq / n );
	ystddev = sqrt( yssq / n );
	xycov = xyssq / n;
	r = xycov / (xstddev * ystddev);
	return r;
}

REAL8 XLALPearsonCrossCorrelationSlow( REAL8 *x, REAL8 *y, UINT4 n )
{
	REAL8 xave = 0;
	REAL8 yave = 0;
	REAL8 xssq = 0;
	REAL8 yssq = 0;
	REAL8 xycor = 0;
	UINT4 i;

	for ( i = 0; i < n; ++i ) {
		xave += x[i];
		yave += y[i];
	}
	xave /= n;
	yave /= n;

	for ( i = 0; i < n; ++i ) {
		REAL8 dx = x[i] - xave;
		REAL8 dy = y[i] - yave;
		xssq  += dx * dx;
		yssq  += dy * dy;
		xycor += dx * dy;
	}
	xycor /= sqrt( xssq * yssq );
	return xycor;
}

REAL8Vector * XLALPearsonCrossCorrelationREAL8Vector( REAL8Vector *x, REAL8Vector *y, UINT4 seglen )
{
	REAL8Vector *r;
	UINT4 stride = seglen/2;
	UINT4 reclen = x->length;
	UINT4 numseg = 1 + (reclen - seglen)/stride;
	UINT4 seg;
	r = XLALCreateREAL8Vector( numseg );
	for ( seg = 0; seg < numseg; ++seg )
		r->data[seg] = XLALPearsonCrossCorrelation( x->data + seg*stride, y->data + seg*stride, seglen );
	return r;
}

REAL8 XLALPearsonMaxCrossCorrelationREAL8Vector( REAL8Vector *x, REAL8Vector *y, UINT4 seglen )
{
	UINT4 stride = seglen/2;
	UINT4 reclen = x->length;
	UINT4 numseg = 1 + (reclen - seglen)/stride;
	UINT4 seg;
	REAL8 rmax = -1.0;
	REAL8 r;
	for ( seg = 0; seg < numseg; ++seg )
		if ( rmax < ( r = XLALPearsonCrossCorrelation( x->data + seg*stride, y->data + seg*stride, seglen ) ) )
			rmax = r;
	return rmax;
}
