/*-----------------------------------------------------------------------
 *
 * File Name: Thresholds.h
 *
 * Author: Eanna Flanagan
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Thresholds.h
 *
 * SYNOPSIS
 * #include "Thresholds.h"
 *
 * DESCRIPTION
 * Error codes, typedefs, and protypes for the functions related
 * to thresholds for chi-squared and non-central chi-squared
 * distributions
 * 
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */


#ifndef _THRESHOLDS_H
#define _THRESHOLDS_H


#include "LALStdlib.h"
#include "LALDatatypes.h"
#include "LALRCSID.h"

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (THRESHOLDSH, "$Id$");


#define THRESHOLDS_ENULLP    1
#define THRESHOLDS_EPOSARG   2
#define THRESHOLDS_EMXIT     4
#define THRESHOLDS_EBADPROB  8
#define THRESHOLDS_ERANGE    16


#define THRESHOLDS_MSGENULLP    "Null pointer"
#define THRESHOLDS_MSGEPOSARG   "Arguments must be non-negative"
#define THRESHOLDS_MSGEMXIT     "Maximum iterations exceeded"
#define THRESHOLDS_MSGEBADPROB  "Supplied probability must be between 0 and 1"
#define THRESHOLDS_MSGERANGE    "Arguments too large, cannot obtain finite probability"


typedef struct tagChisqCdfIn
{
  REAL8     chi2;               /* value of chi squared          */
  REAL8     dof;                /* number of degrees of freedom  */
  REAL8     nonCentral;         /* non-centrality parameter      */
}
ChisqCdfIn;

typedef struct tagChi2ThresholdIn
{
  REAL8     dof;                /* number of degrees of freedom  */
  REAL8     falseAlarm;         /* false alarm probability       */
}
Chi2ThresholdIn;

typedef struct tagRhoThresholdIn
{
  REAL8     chi2;               /* value of chi squared          */
  REAL8     dof;                /* number of degrees of freedom  */
  REAL8     falseDismissal;     /* false dismissal probability   */
}
RhoThresholdIn;



void
ChisqCdf (
          Status                        *status,
          REAL8                         *prob,
          ChisqCdfIn                    *input
          );

void
OneMinusChisqCdf (
                  Status                *status,
                  REAL8                 *prob,
                  ChisqCdfIn            *input
                  );

void
NoncChisqCdf (
              Status                    *status,
              REAL8                     *prob,
              ChisqCdfIn                *input
              );

void
Chi2Threshold (
              Status                    *status,
              REAL8                     *chi2,
              Chi2ThresholdIn           *input
               );

void
RhoThreshold (
              Status                    *status,
              REAL8                     *rho,
              RhoThresholdIn            *input
              );



#ifdef  __cplusplus
}
#endif  /* C++ protection. */


#endif






