/**** <lalVerbatim file="FrameCalibrationHV">
 * Author: Brown, D. A.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{FrameCalibration.h}}
 *
 * High-level routines for exctracting calibration data from frames
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/FrameCalibration.h>
 * \end{verbatim}
 *
 * Provides a high level interface for building a transfer or response 
 * functions from raw calibration data provided by the calibration team.
 *
 **** </lalLaTeX> */

#include <lal/LALDatatypes.h>

#ifndef _FRAMECALIBRATION_H
#define _FRAMECALIBRATION_H

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( FRAMECALIBRATIONH, "$Id$" );

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define FRAMECALIBRATIONH_ENULL 00001
#define FRAMECALIBRATIONH_ENNUL 00002

#define FRAMECALIBRATIONH_MSGENULL "Null pointer"
#define FRAMECALIBRATIONH_MSGENNUL "Non-null pointer"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 *
 **** </lalLaTeX> */

/**** <lalLaTeX>
 * 
 * \vfill{\footnotesize\input{FrameCalibrationHV}}
 * \newpage\input{FrameCalibrationC}
 *
 **** </lalLaTeX> */

void
LALExtractFrameResponse(
    LALStatus               *status,
    COMPLEX8FrequencySeries *output,
    const CHAR              *catalog,
    const CHAR              *ifo
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _FRAMECALIBRATION_H */
