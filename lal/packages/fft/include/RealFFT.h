/*-----------------------------------------------------------------------
 *
 * File Name: RealFFT.h
 *
 * Author: Creighton, J. D. E.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * RealFFT.h
 *
 * SYNOPSIS
 * #include "RealFFT.h"
 *
 * DESCRIPTION
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _REALFFT_H
#define _REALFFT_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (REALFFTH, "$Id$");

#define REALFFT_ENULL 1
#define REALFFT_ENNUL 2
#define REALFFT_ESIZE 4
#define REALFFT_ESZMM 8
#define REALFFT_ESLEN 16
#define REALFFT_ESAME 32
#define REALFFT_ESIGN 64
#define REALFFT_EDATA 128

#define REALFFT_MSGENULL "Null pointer"
#define REALFFT_MSGENNUL "Non-null pointer"
#define REALFFT_MSGESIZE "Invalid input size"
#define REALFFT_MSGESZMM "Size mismatch"
#define REALFFT_MSGESLEN "Invalid/mismatched sequence lengths"
#define REALFFT_MSGESAME "Input/Output data vectors are the same"
#define REALFFT_MSGESIGN "Incorrect plan sign"
#define REALFFT_MSGEDATA "Bad input data: DC/Nyquist should be real"

typedef struct
tagRealFFTPlan
{
  INT4   sign;
  UINT4  size;
  void  *plan;
}
RealFFTPlan;


void
EstimateFwdRealFFTPlan (
    Status       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
EstimateInvRealFFTPlan (
    Status       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
MeasureFwdRealFFTPlan (
    Status       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
MeasureInvRealFFTPlan (
    Status       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
DestroyRealFFTPlan (
    Status       *stat,
    RealFFTPlan **plan
    );


void
REAL4VectorFFT (
    Status      *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    );

void
REAL4VectorSequenceFFT (
    Status              *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    );


void
FwdRealFFT (
    Status         *stat,
    COMPLEX8Vector *vout,
    REAL4Vector    *vinp,
    RealFFTPlan    *plan
    );

void
InvRealFFT (
    Status         *stat,
    REAL4Vector    *vout,
    COMPLEX8Vector *vinp,
    RealFFTPlan    *plan
    );

void
RealPowerSpectrum (
    Status      *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    );


void
FwdRealSequenceFFT (
    Status                 *stat,
    COMPLEX8VectorSequence *vout,
    REAL4VectorSequence    *vinp,
    RealFFTPlan            *plan
    );

void
InvRealSequenceFFT (
    Status                 *stat,
    REAL4VectorSequence    *vout,
    COMPLEX8VectorSequence *vinp,
    RealFFTPlan            *plan
    );

void
RealSequencePowerSpectrum (
    Status              *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    );


#ifdef  __cplusplus
}
#endif

#endif /* _REALFFT_H */
