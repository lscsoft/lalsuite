/*----------------------------------------------------------------------- 
 * 
 * File Name: RealPowerSpectrumTest.c
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
 *   Constructs real power spectra of data with a variety of widows.
 * 
 * DIAGNOSTICS
 * 
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _SEQFACTORIES_H
#include "SeqFactories.h"
#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H
#endif
#endif

#ifndef _REALFFT_H
#include "RealFFT.h"
#ifndef _REALFFT_H
#define _REALFFT_H
#endif
#endif

#ifndef _WINDOW_H
#include "Window.h"
#ifndef _WINDOW_H
#define _WINDOW_H
#endif
#endif

#ifndef _VECTOROPS_H
#include "VectorOps.h"
#ifndef _VECTOROPS_H
#define _VECTOROPS_H
#endif
#endif

NRCSID (MAIN, "$Id$");

int debuglevel = 2;

int main()
{
  const INT4 m = NUMBERWINDOWTYPES;
  const INT4 n = 65536;

  static Status status;

  RealFFTPlan            *plan = NULL;
  REAL4Vector            *wss  = NULL;
  REAL4Vector            *hvec = NULL;
  REAL4VectorSequence    *hseq = NULL;
  REAL4VectorSequence    *Pseq = NULL;
  CreateVectorSequenceIn  hinp;
  CreateVectorSequenceIn  Pinp;

  FILE *fp = fopen ("RealPowerSpectrumTest.out", "w");
  
  INT4 j;
  INT4 k;

  /* create vectors and sequences */
  hinp.length       = m;
  hinp.vectorLength = n;
  Pinp.length       = m;
  Pinp.vectorLength = n/2 + 1;

  SCreateVector          (&status, &hvec, n);
  SCreateVectorSequence  (&status, &hseq, &hinp);
  SCreateVectorSequence  (&status, &Pseq, &Pinp);
  EstimateFwdRealFFTPlan (&status, &plan, n);

  /* initialize raw data vector */
  for (k = 0; k < hvec->length; ++k)
    hvec->data[k] = sin(0.001*k);

  /* create window sum-of-squares vector */
  SCreateVector (&status, &wss, m);

  for (j = 0; j < hseq->length; ++j)
  {
    REAL4Vector   dum;
    REAL4Vector  *win  = NULL;
    WindowParams  winp;

    dum.length  = hseq->vectorLength;
    dum.data    = hseq->data + j*hseq->vectorLength;

    winp.length = hseq->vectorLength;
    winp.type   = (WindowType) j;

    /* create window */
    SCreateVector (&status, &win, winp.length);
    Window (&status, win, &winp);
    wss->data[j] = winp.sumofsquares;

    /* apply window to data */
    SSVectorMultiply (&status, &dum, hvec, win);

    /* destroy window */
    SDestroyVector (&status, &win);
  }

  /* compute power spectra */
  RealSequencePowerSpectrum (&status, Pseq, hseq, plan);

  /* print power spectra omitting DC component */
  for (k = 1; k < Pseq->vectorLength; ++k)
  {
    fprintf ( fp, "%d", k );
    for (j = 0; j < Pseq->length; ++j)
    {
      REAL4 out = Pseq->data[k + j*Pseq->vectorLength]/wss->data[j];
      fprintf ( fp, "\t%e", ( out == 0 ? 1e-30 : out ) );
    }
    fprintf ( fp, "\n" );
  }
  fclose (fp);

  SDestroyVector         (&status, &wss );
  SDestroyVector         (&status, &hvec);
  SDestroyVectorSequence (&status, &hseq);
  SDestroyVectorSequence (&status, &Pseq);
  DestroyRealFFTPlan     (&status, &plan);
  return 0;
}
