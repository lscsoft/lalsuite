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
#define CALIBRATIONH_MSGENULL "Null pointer"
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
 * \newpage\input{CalibrationC}
 * \newpage\input{CalibrationTestC}
 *
 **** </lalLaTeX> */

void LALComputeTransfer( LALStatus *status, CalibrationRecord *calrec );

#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _CALIBRATION_H */
