/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Xavier Siemens
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
#include <lal/BandPassTimeSeries.h>

NRCSID (CALIBRATIONH,"$Id$");

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
 * \idx[Type]{CalibrationFunctions}
 * \idx[Type]{CalibrationUpdateParams}
 *
 * \subsubsection*{Type \texttt{CalibrationType}}
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
 * \subsubsection*{Type \texttt{CalFactors}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalFactors
{
  COMPLEX16 alpha;
  COMPLEX16 alphabeta;
  COMPLEX16 beta;
  COMPLEX16 exc;
  COMPLEX16 asq;
  COMPLEX16 darm;
}
CalFactors;
/**** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * Document \verb+CalibrationType+
 *
 * \subsubsection*{Type \texttt{UpdateFactorsParams}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagUpdateFactorsParams
{
   REAL8 lineFrequency;
   COMPLEX16 openloop;
   COMPLEX16 digital;
   COMPLEX16 whitener;
   REAL4TimeSeries *darmCtrl;
   REAL4TimeSeries *asQ;
   REAL4TimeSeries *exc;
}
UpdateFactorsParams;
/**** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * Document \verb+CalibrationType+
 *
 * \subsubsection*{Type \texttt{CalibrationRecord}}
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
 * \subsubsection*{Type \texttt{CalibrationFunctions}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalibrationFunctions
{
  COMPLEX8FrequencySeries *responseFunction;
  COMPLEX8FrequencySeries *sensingFunction;
}
CalibrationFunctions;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * The type \texttt{CalibrationFunctions} contains two calibration functions,
 * the sensing function $C(f)$ and the response function $R(f)$.  While the
 * response function is the function that is most often wanted, the sensing
 * function is needed in updating calibration from one epoch to another.
 *
 * \subsubsection*{Type \texttt{CalibrationUpdateParams}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalibrationUpdateParams
{
  LIGOTimeGPS epoch;
  LIGOTimeGPS duration;
  CHAR *ifo;
  COMPLEX8 alpha;
  COMPLEX8 alphabeta;
  COMPLEX8TimeSeries *openLoopFactor;
  COMPLEX8TimeSeries *sensingFactor;
}
CalibrationUpdateParams;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * The type \texttt{CalibrationUpdateParams} contains two time series
 * representing an overall gain factor for the open-loop gain function $H(f)$
 * and the sensing function $C(f)$.  These transfer functions are known to
 * change (by an overall factor) with time, and these two factors can be
 * tracked using the injected calibration lines.  The factors are stored
 * in this structure as (very-slowly varying) time series, and are to be
 * used in updating the calibration functions described previously.
 * (The response function can be computed from the open-loop gain and the
 * sensing function.  It is simply $R(f)=[1+H(f)]/C(f)$.)  In addition, this
 * structure contains the present epoch and the duration of the data to be
 * calibrated to identify the particular set of
 * factors (from those recorded in the time series) to use.
 *
 * \vfill{\footnotesize\input{CalibrationHV}}
 * \newpage\input{ComputeTransferC}
 *
 **** </lalLaTeX> */

typedef
struct StrainOutTag {
  REAL8TimeSeries h;         /* timeseries containing h(t) */
  REAL8TimeSeries hC;         /* timeseries containing the control signal */
  REAL8TimeSeries hR;         /* timeseries containing the residual signal */
  COMPLEX16TimeSeries alpha; /* alpha time series */
  COMPLEX16TimeSeries beta;  /* beta time series */
  COMPLEX16TimeSeries alphabeta; /* alpha time series */
  INT2TimeSeries science_mode;   /* flag = 1 for science mode, 0 otherwise */
} StrainOut;

typedef
struct StrainInTag {
  REAL4TimeSeries AS_Q ;   /* timeseries containing ASQ */
  REAL4TimeSeries DARM_ERR;/* timeseries containing DARM_ERR */
  REAL4TimeSeries DARM ;   /* timeseries containing DARM_CTRL */
  REAL4TimeSeries EXC ;    /* timeseries containing the excitation */
  COMPLEX16 Do;            /* digital filter at cal line frequency */
  COMPLEX16 Go;            /* OLG at cal line frequency */
  COMPLEX16 Wo;            /* Whitening filter at cal line frequency */
  REAL8 f;                 /* calibration line frequency */
  REAL4 To;                /* factors integration time */
  REAL8IIRFilter *Cinv;    /* Filters for inverse of sensing function */
  INT4 CinvUSF;            /* Upsampling factor for sensing function */
  INT4 CinvDelay;          /* Overall inverse sensing function delay */
  REAL8IIRFilter *A;      /* Filters for analog actuation function */
  REAL8IIRFilter *D;      /* Filters for analog actuation function */
  REAL8IIRFilter *AW;      /* Filters for analog actuation function */
  REAL8 gamma_fudgefactor;
  INT4 delta;
  INT4 usefactors;
  INT4 wings;               /* size of wings in seconds */
  INT4 fftconv;
  INT4 outalphas;
  INT4 darmctrl;
  /* Stuff needed to run old IIR code */
  REAL8IIRFilter *AA;      /* Filters for analog actuation function */
  INT4 AADelay;            /* Overall analog actuation function delay */
  REAL8IIRFilter *AX;      /* Digital filters for x arm actuation function */
  REAL8IIRFilter *AY;      /* Digital filters for y arm actuation function */
  INT4 NCinv;              /* Numbers of filters of each type */
  INT4 ND;
  INT4 NAA;
  INT4 NAX;
  INT4 NAY;

} StrainIn;

typedef
struct MyIIRFilter {
  INT4 yOrder;
  INT4 xOrder;
  REAL8 a[20];
  REAL8 b[20];
  REAL8 yhist[20];
  REAL8 xhist[20];
} MyIIRFilter;



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

INT4
XLALResponseConvert(
    COMPLEX8FrequencySeries *output,
    COMPLEX8FrequencySeries *input
    );

void LALComputeCalibrationFactors(
    LALStatus              *status,
    CalFactors             *output,
    UpdateFactorsParams    *input
    );

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,
    StrainIn               *input
    );

void LALComputeStrainDMT(
    LALStatus              *status,
    StrainOut              *output,
    StrainIn               *input
    );

void LALGetFactors(
    LALStatus              *status,
    StrainOut              *output,
    StrainIn               *input
    );

int XLALFIRFilter(
    REAL8TimeSeries        *tseries,
    REAL8IIRFilter         *FIR
    );

void LALMakeFIRLP(LALStatus *status,
		  REAL8IIRFilter *G,
		  int USF);

void LALMakeFIRLPALPHAS(LALStatus *status,
		  REAL8IIRFilter *G);

void LALMakeFIRHP(LALStatus *status,
		  REAL8IIRFilter *G);

int XALFIRFilter(REAL8TimeSeries *tseries,
		  REAL8IIRFilter *FIR);

void LALFFTFIRFilter(LALStatus *status,
		     REAL8TimeSeries *tseries,
		     REAL8IIRFilter *FIR);

void LALFIRFilter(LALStatus *status,
		     REAL8TimeSeries *tseries,
		     REAL8IIRFilter *FIR);

void LALFreeFilter(LALStatus *status, REAL8IIRFilter *F2, int ORDER);
void LALCopyFilter(LALStatus *status, REAL8IIRFilter **F2, REAL8IIRFilter *F1, int ORDER);

int XLALDivideTimeSeries(REAL8TimeSeries *hR, REAL8TimeSeries *ALPHAS);
int XLALUpsample(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor);
int XLALUpsampleLinear(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor);

#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _CALIBRATION_H */
