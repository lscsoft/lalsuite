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

#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/VectorOps.h>

NRCSID (MAIN, "$Id$");

int lalDebugLevel = 2;

int main( void )
{
  const INT4 m = NumberWindowTypes;
  const INT4 n = 65536;

  static LALStatus status;

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

  LALSCreateVector          (&status, &hvec, n);
  LALSCreateVectorSequence  (&status, &hseq, &hinp);
  LALSCreateVectorSequence  (&status, &Pseq, &Pinp);
  LALEstimateFwdRealFFTPlan (&status, &plan, n);

  /* initialize raw data vector */
  for (k = 0; k < (INT4) hvec->length; ++k)
    hvec->data[k] = sin(0.001*k);

  /* create window sum-of-squares vector */
  LALSCreateVector (&status, &wss, m);

  for (j = 0; j < (INT4) hseq->length; ++j)
  {
    REAL4Vector   dum;
    REAL4Vector  *win  = NULL;
    LALWindowParams  winp;

    dum.length  = hseq->vectorLength;
    dum.data    = hseq->data + j*hseq->vectorLength;

    winp.length = hseq->vectorLength;
    winp.type   = (WindowType) j;

    /* create window */
    LALSCreateVector (&status, &win, winp.length);
    LALWindow (&status, win, &winp);
    wss->data[j] = winp.sumofsquares;

    /* apply window to data */
    LALSSVectorMultiply (&status, &dum, hvec, win);

    /* destroy window */
    LALSDestroyVector (&status, &win);
  }

  /* compute power spectra */
  LALRealSequencePowerSpectrum (&status, Pseq, hseq, plan);

  /* print power spectra omitting DC component */
  for (k = 1; k < (INT4) Pseq->vectorLength; ++k)
  {
    fprintf ( fp, "%d", k );
    for (j = 0; j < (INT4) Pseq->length; ++j)
    {
      REAL4 out = Pseq->data[k + j*Pseq->vectorLength]/wss->data[j];
      fprintf ( fp, "\t%e", ( out == 0 ? 1e-30 : out ) );
    }
    fprintf ( fp, "\n" );
  }
  fclose (fp);

  LALSDestroyVector         (&status, &wss );
  LALSDestroyVector         (&status, &hvec);
  LALSDestroyVectorSequence (&status, &hseq);
  LALSDestroyVectorSequence (&status, &Pseq);
  LALDestroyRealFFTPlan     (&status, &plan);
  return 0;
}
