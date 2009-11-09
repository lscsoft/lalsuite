/*
*  Copyright (C) 2007 Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: CLRTest.c
 * Author: Sintes, A. M.
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
 *   Test suite for CLR  operations
 *
 * DIAGNOSTICS
 *
 * CALLS
 *   LALI4CreateVector
 *   LALCCreateVector
 *   LALSCreateVector
 *   LALMalloc
 *   LALFree
 *   LALCreateForwardRealFFTPlan
 *   LALRealPowerSpectrum
 *   LALForwardRealFFT
 *   LALDestroyRealFFTPlan
 *   LALI4DestroyVector
 *   LALCDestroyVector
 *   LALSDestroyVector
 *   LALHarmonicFinder
 *   LALRefInterference
 *   LALCleanAll
 *   LALCreateForwardComplexFFTPlan
 *   LALCOMPLEX8VectorFFT
 *   LALDestroyComplexFFTPlan
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


/************************************ <lalVerbatim file="CLRTestCV">
Author: Sintes, A. M.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>

\subsection{Program \texttt{CLRTest.c}}
\label{s:CLRTest.c}

 Test for CLR  operations.

\subsubsection*{Usage}
\begin{verbatim}
CLRTest
\end{verbatim}

\subsubsection*{Description}
This program is just an example of the usage of the different prototypes.

 The program reads some data from the file
\texttt{CLRindata.asc}, finds
the position of several harmonics, builds a reference signal,
cleans the initial data of all interference harmonics and
writes the clean data into the file
\texttt{CLRoutdata.asc}.

\subsubsection*{Exit codes}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALI4CreateVector()
LALCCreateVector()
LALSCreateVector()
LALMalloc()
LALFree()
LALCreateForwardRealFFTPlan()
LALRealPowerSpectrum()
LALForwardRealFFT()
LALDestroyRealFFTPlan()
LALI4DestroyVector()
LALCDestroyVector()
LALSDestroyVector()
LALHarmonicFinder()
LALRefInterference()
LALCleanAll()
LALCreateReverseComplexFFTPlan()
LALCOMPLEX8VectorFFT()
LALDestroyComplexFFTPlan()
\end{verbatim}

\subsubsection*{Notes}
Take this program just as an example,  build
your own one  and feed it  with the data of your interest.
The CLR functions work on stretches of data from a few seconds up to a
couple of minutes.

\vfill{\footnotesize\input{CLRTestCV}}

</lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/CLR.h>

NRCSID (MAIN, "$Id$");

INT4  lalDebugLevel = 2;

int main ( void )
{
    /* example of operation instructions */

  static LALStatus  status;
  const  UINT4   n = 64*4096; /* example vector length */
  const  UINT4   l =  7; /* number of harmonics to build the ref. signa l*/

  FILE   *in_file, *out_file;
  char   filename[100];
  int    number;

  INT4            i;

  INT4Vector     *hk   = NULL;   /* the harmonic index (l) */
  INT4Vector     *hkff = NULL;  /* harmonic index and bin location (3*l) */

  RealFFTPlan    *pfwd = NULL;

  REAL4TVectorCLR    *xt = NULL;  /* x(t), data + information */
  REAL4FVectorCLR    *xp = NULL;  /* |x(f)|^2, data + information */

  COMPLEX8Vector     *xf = NULL;   /* x(nu), size n/2+1 */
  COMPLEX8Vector     *mt = NULL;   /* m(t), size n */

  REAL4Vector    *xtclean = NULL; /* clean data x(t), size n */
  REAL4Vector    *x       = NULL;       /* data x(t), size n */
  REAL4Vector    *Pvec    = NULL; /* Power spectrum, size n/2+1 */

  REAL4    dummy;


  /* -------------------------------------- */
  /* create data vectors, plans... */

  LALI4CreateVector(&status, &hk, l);
  LALI4CreateVector(&status, &hkff, 3*l);

  LALCreateForwardRealFFTPlan(&status, &pfwd, n, 0);

  xt =  (REAL4TVectorCLR *)LALMalloc(sizeof(REAL4TVectorCLR));
  xp =  (REAL4FVectorCLR *)LALMalloc(sizeof(REAL4FVectorCLR));

  LALCCreateVector(&status, &xf, n/2+1);
  LALCCreateVector(&status, &mt, n);

  LALSCreateVector(&status, &xtclean, n);
  LALSCreateVector(&status, &x,  n);
  LALSCreateVector(&status, &Pvec, n/2+1);

  /* ---------------------------------------- */
  /* assign data */

  /* the harmonics to be considered to build the reference signal */
  hk->data[0] = 3;
  hk->data[1] = 5;
  hk->data[2] = 9;
  hk->data[3] = 11;
  hk->data[4] = 13;
  hk->data[5] = 15;
  hk->data[6] = 19;

  /* The  CLR Time Vector  */
  xt->length = n;
  xt->data = x->data;
  xt->deltaT = 1.0/4000.0; /* inverse of the sampling frequency */
  xt->fLine = 50.0;        /* or 60.0 Hz */

  /* The  CLR Frequency Vector  */
  xp->length = n/2+1;
  xp->data = Pvec->data;
  xp->deltaF = 1.0/( xt->length *  xt->deltaT );
  xp->fLine =  xt->fLine;

  /* ----------------------- */
  /* read data  x(t) */
  /* insert here your own data from a given file/frame */

  strcpy(filename,"CLRindata.asc\0");
  in_file = LALOpenDataFile(filename);
  for (i = 0; i < (int)n; ++i) {
    number = fscanf(in_file, "%f\n", &dummy );
    x->data[i] = dummy;
  }
  LALFclose(in_file);

  /* --------------------------------------------------- */
  /*          what the program should do                 */
  /* --------------------------------------------------- */

  /* compute Spectrum */

  LALRealPowerSpectrum(&status,Pvec,x,pfwd);
  /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
  {
    REAL4Vector *myvector = Pvec;
    UINT4 mybin;
    for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
      myvector->data[mybin] *= 0.5;
  }



  /* find the position of the harmonics considered */
  LALHarmonicFinder(&status,hkff,xp,hk);


  /* for debugging only */
  for (i = 0; i< 3*(int)l; ++i)
    printf(" %d \n", hkff->data[i]);

  /* --------------------------------------------------- */
  /* this information could be provided  as an input, e.g.: */

  /* hkff->data[1] =  9868; */
  /* hkff->data[2] =  9894; */

  /* hkff->data[4] =  16449; */
  /* hkff->data[5] = 16487; */

  /* hkff->data[7] = 29607; */
  /* hkff->data[8] = 29675; */

  /* hkff->data[10] =  36189 ; */
  /* hkff->data[11] =  36267; */

  /* hkff->data[13] =  42761; */
  /* hkff->data[14] = 42871; */

  /* hkff->data[16] = 49335 ; */
  /* hkff->data[17] =  49465; */

  /* hkff->data[19] =  62498; */
  /* hkff->data[20] = 62654; */

  /* ------------------------------- */

  /* perform fft */
  LALForwardRealFFT(&status,xf,x,pfwd);

  /* generate the reference signal */
  LALRefInterference(&status,mt,xf,hkff);

  /* clean the data of all harmonics */
  LALCleanAll(&status,xtclean,mt,xt);

  /* ------------------------------------------------- */
  /* write clean  data  x(t) */

  strcpy(filename,"CLRoutdata.asc\0");
  out_file = fopen(filename,"w");
  for (i = 0; i < (int)n; ++i) {
    fprintf(out_file, "%f\n", xtclean->data[i] );
    fflush(out_file);
  }
  fclose(out_file);


  /* -------------------------------------- */
  /* destroy data vectors, plans... if not done before */

  LALI4DestroyVector(&status, &hk);
  LALI4DestroyVector(&status, &hkff);

  LALDestroyRealFFTPlan (&status,&pfwd);

  LALCDestroyVector(&status, &xf);
  LALCDestroyVector(&status, &mt);

  LALSDestroyVector(&status, &x);
  LALSDestroyVector(&status, &xtclean);
  LALSDestroyVector(&status, &Pvec);

  LALFree(xt);
  LALFree(xp);

  /*--------------------------------*/
  return 0;
}
