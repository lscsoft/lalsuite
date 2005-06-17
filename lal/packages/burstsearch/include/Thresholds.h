/********************************** <lalVerbatim file="ThresholdsHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _THRESHOLDS_H
#define _THRESHOLDS_H

#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (THRESHOLDSH, "$Id$");

REAL8
XLALChisqCdf(
	REAL8 chi2,
	REAL8 dof
);

REAL8
XLALOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
);

REAL8
XLALNoncChisqCdf (
	REAL8 chi2,
	REAL8 dof,
	REAL8 nonCentral
);

REAL8
XLALNoncChisqCdfNonSafe (
	REAL8 chi2,
	REAL8 dof,
	REAL8 nonCentral
);

REAL8 XLALChi2Threshold(
	REAL8 dof,
	REAL8 falseAlarm
);

REAL8 XLALRhoThreshold(
	REAL8 chi2,
	REAL8 dof,
	REAL8 falseDismissal
);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
