/**** <lalVerbatim file="CalibrationHV">
 * Author: P. R. Brady, J. D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{Calibration.h}}
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Calibration.h>
 * \end{verbatim}
 *
 **** </lalLaTeX> */

#ifndef _CALIBRATION_H
#define _CALIBRATION_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma } /** to match the previous brace **/
#endif

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define CALIBRATIONH_ENULL 001
#define CALIBRATIONH_ESIZE 002
#define CALIBRATIONH_ESZMM 004
#define CALIBRATIONH_EZERO 010
#define CALIBRATIONH_ETIME 020
#define CALIBRATIONH_EUNIT 040

#define CALIBRATIONH_MSGENULL "Null pointer"
#define CALIBRATIONH_MSGESIZE "Invalid size"
#define CALIBRATIONH_MSGESZMM "Size mismatch"
#define CALIBRATIONH_MSGEZERO "Zero factor"
#define CALIBRATIONH_MSGETIME "Time out of range"
#define CALIBRATIONH_MSGEUNIT "Incompatible units"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{CalibrationType}
 * \idx[Type]{CalibrationRecord}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef enum
{
  CalibrationAmplitude = 001,
  CalibrationOffset    = 002,
  CalibrationDelay     = 004,
  CalibrationTransfer  = 010,
  CalibrationZPG       = 020
}
CalibrationType;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * Document \verb+CalibrationType+
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalibrationRecord
{
  CHAR                     name[LALNameLength];
  LIGOTimeGPS              epoch;
  REAL8                    duration;
  CHAR                     reference[LALNameLength];
  LALUnit                  units;
  UINT4                    type;
  REAL8                    conversion;
  REAL8                    offset;
  REAL8                    delay;
  COMPLEX8FrequencySeries *transfer;
  REAL8Vector             *zeros;
  REAL8Vector             *poles;
  REAL8                    gain;
}
CalibrationRecord;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * Document \verb+CalibrationRecord+
 *
 * \vfill{\footnotesize\input{CalibrationHV}}
 * \newpage\input{ComputeTransferC}
 *
 **** </lalLaTeX> */


typedef struct
tagCalibrationFunctions
{
  COMPLEX8FrequencySeries *responseFunction;
  COMPLEX8FrequencySeries *sensingFunction;
}
CalibrationFunctions;

typedef struct
tagCalibrationUpdateParams
{
  LIGOTimeGPS epoch;
  COMPLEX8TimeSeries *openLoopFactor;
  COMPLEX8TimeSeries *sensingFactor;
}
CalibrationUpdateParams;




void LALComputeTransfer( LALStatus *status, CalibrationRecord *calrec );


void
LALUpdateCalibration(
    LALStatus               *status,
    CalibrationFunctions    *output,
    CalibrationFunctions    *input,
    CalibrationUpdateParams *params
    );

void
LALResponseConvert(
    LALStatus               *status,
    COMPLEX8FrequencySeries *output,
    COMPLEX8FrequencySeries *input
    );




#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _CALIBRATION_H */
