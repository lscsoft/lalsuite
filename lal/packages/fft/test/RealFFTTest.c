/*----------------------------------------------------------------------- 
 * 
 * File Name: RealFFTTest.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 *   main()
 *
 * SYNOPSIS 
 * 
 * DESCRIPTION 
 *   Tests real FFT functions.
 * 
 * DIAGNOSTICS
 * 
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include "LALStdlib.h"
#include "SeqFactories.h"
#include "RealFFT.h"
#include "VectorOps.h"


NRCSID (MAIN, "$Id$");

int lalDebugLevel = 1;

int main()
{
  const INT4 m = 4;   /* example length of sequence of vectors */
  const INT4 n = 32;  /* example vector length */

  static LALStatus status; 

  RealFFTPlan            *pfwd = NULL;
  RealFFTPlan            *pinv = NULL;
  REAL4Vector            *hvec = NULL;
  COMPLEX8Vector         *Hvec = NULL;
  REAL4Vector            *Pvec = NULL;
  REAL4VectorSequence    *hseq = NULL;
  COMPLEX8VectorSequence *Hseq = NULL;
  REAL4VectorSequence    *Pseq = NULL;
  CreateVectorSequenceIn  seqinp;
  CreateVectorSequenceIn  seqout;

  INT4 i;
  INT4 j;

  /* create vectors and sequences */

  seqinp.length       = m;
  seqinp.vectorLength = n;
  seqout.length       = m;
  seqout.vectorLength = n/2 + 1;

  LALSCreateVector         (&status, &hvec, n);
  REPORTSTATUS (&status);

  LALCCreateVector         (&status, &Hvec, n/2 + 1);
  REPORTSTATUS (&status);

  LALSCreateVector         (&status, &Pvec, n/2 + 1);
  REPORTSTATUS (&status);

  LALSCreateVectorSequence (&status, &hseq, &seqinp);
  REPORTSTATUS (&status);

  LALCCreateVectorSequence (&status, &Hseq, &seqout);
  REPORTSTATUS (&status);

  LALSCreateVectorSequence (&status, &Pseq, &seqout);
  REPORTSTATUS (&status);

  /* create plans */

  LALEstimateFwdRealFFTPlan (&status, &pfwd, n);
  REPORTSTATUS (&status);

  LALEstimateInvRealFFTPlan (&status, &pinv, n);
  REPORTSTATUS (&status);

  LALDestroyRealFFTPlan     (&status, &pfwd);
  REPORTSTATUS (&status);

  LALDestroyRealFFTPlan     (&status, &pinv);
  REPORTSTATUS (&status);

  LALMeasureFwdRealFFTPlan  (&status, &pfwd, n);
  REPORTSTATUS (&status);

  LALMeasureInvRealFFTPlan  (&status, &pinv, n);
  REPORTSTATUS (&status);

  /* assign data ... */
  printf ("\nInitial data:\n");
  for (i = 0; i < n; ++i)
  {
    REAL4 data = cos(7*i) - 1;
    hvec->data[i] = data;
    for (j = 0; j < m; ++j)
    {
      hseq->data[i + j*n] = data;
      printf ("% 9.3f\t", data);
      data += 1;
    }
    printf ("\n");
  }

  /* perform FFTs */

  printf ("\nSingle Forward FFT:\n");
  LALFwdRealFFT (&status, Hvec, hvec, pfwd);
  REPORTSTATUS (&status);
  for (i = 0; i < Hvec->length; ++i)
    printf ("(% 9.3f, % 9.3f)\n", Hvec->data[i].re, Hvec->data[i].im);

  printf ("\nSingle Forward FFT Power:\n");
  LALRealPowerSpectrum (&status, Pvec, hvec, pfwd);
  REPORTSTATUS (&status);
  for (i = 0; i < Pvec->length; ++i)
    printf ("%12.3f\n", Pvec->data[i]);

  printf ("\nSingle Inverse FFT:\n");
  LALInvRealFFT (&status, hvec, Hvec, pinv);
  REPORTSTATUS (&status);
  for (i = 0; i < hvec->length; ++i)
    printf ("% 9.3f\n", hvec->data[i]/n);

  printf ("\nMultiple Forward FFT:\n");
  LALFwdRealSequenceFFT (&status, Hseq, hseq, pfwd);
  REPORTSTATUS (&status);
  for (i = 0; i < Hseq->vectorLength; ++i)
  {
    for (j = 0; j < Hseq->length; ++j)
      printf ("(% 9.3f, % 9.3f)\t",
                    Hseq->data[i + j*Hseq->vectorLength].re,
                    Hseq->data[i + j*Hseq->vectorLength].im);
    printf ("\n");
  }

  printf ("\nMultiple Forward FFT Power:\n");
  LALRealSequencePowerSpectrum (&status, Pseq, hseq, pfwd);
  REPORTSTATUS (&status);
  for (i = 0; i < Pseq->vectorLength; ++i)
  {
    for (j = 0; j < Pseq->length; ++j)
      printf ("%12.3f\t", Pseq->data[i + j*Pseq->vectorLength]);
    printf ("\n");
  }

  printf ("\nMultiple Inverse FFT:\n");
  LALInvRealSequenceFFT (&status, hseq, Hseq, pinv);
  REPORTSTATUS (&status);
  for (i = 0; i < hseq->vectorLength; ++i)
  {
    for (j = 0; j < hseq->length; ++j)
      printf ("% 9.3f\t", hseq->data[i + j*hseq->vectorLength]/n);
    printf ("\n");
  }

  /* destroy plans, vectors, and sequences */
  LALDestroyRealFFTPlan     (&status, &pfwd);
  REPORTSTATUS (&status);

  LALDestroyRealFFTPlan     (&status, &pinv);
  REPORTSTATUS (&status);

  LALSDestroyVector         (&status, &hvec);
  REPORTSTATUS (&status);

  LALCDestroyVector         (&status, &Hvec);
  REPORTSTATUS (&status);

  LALSDestroyVector         (&status, &Pvec);
  REPORTSTATUS (&status);

  LALSDestroyVectorSequence (&status, &hseq);
  REPORTSTATUS (&status);

  LALCDestroyVectorSequence (&status, &Hseq);
  REPORTSTATUS (&status);

  LALSDestroyVectorSequence (&status, &Pseq);
  REPORTSTATUS (&status);

  LALCheckMemoryLeaks ();
  return 0;
}
