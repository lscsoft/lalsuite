/**** <lalVerbatim file="TimeFreqFFTCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Module \texttt{TimeFreqFFT.c}}
 * \label{ss:TimeFreqFFT.c}
 * 
 * Functions for time to frequency Fourier transforms.
 *
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{TimeFreqFFTCP}
 * \idx{LALTimeFreqRealFFT()}
 * \idx{LALFreqTimeRealFFT()}
 * \idx{LALTimeFreqComplexFFT()}
 * \idx{LALFreqTimeComplexFFT()}
 *
 * \subsubsection*{Description}
 * 
 * The routines \verb+LALTimeFreqRealFFT()+ and \verb+LALTimeFreqComplexFFT()+
 * transform time series $h_j$, $0\le j<n$, into a frequency series
 * $\tilde{h}_k$.  For \verb+LALTimeFreqRealFFT()+,
 * \[
 *    \tilde{h}_k = \Delta t \times H_k \;
 *    \mbox{for $0\le k\le\lfloor n/2\rfloor$.}
 * \]
 * The packing covers the range from dc (inclusive) to Nyquist (inclusive if
 * $n$ is even).
 * For \verb+LALTimeFreqComplexFFT()+,
 * \[
 *    \tilde{h}_k = \Delta t \left\{
 *    \begin{array}{ll}
 *      H_{k+\lfloor(n+1)/2\rfloor} &
 *        \mbox{for $0\le k<\lfloor n/2\rfloor$}, \\
 *      H_{k-\lfloor n/2\rfloor} &
 *        \mbox{for $\lfloor n/2\rfloor\le k<n$}. \\
 *    \end{array}
 *    \right.
 * \]
 * The packing covers the range from negative Nyquist (inclusive if $n$ is
 * even) up to (but not including) positive Nyquist.
 * Here $H_k$ is the DFT of $h_j$:
 * \[
 *   H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}.
 * \]
 * The units of $\tilde{h}_k$ are equal to the units of $h_j$ times seconds.
 *
 * The routines \verb+LALFreqTimeRealFFT()+ and \verb+LALFreqTimeComplexFFT()+
 * perform the inverse transforms from $\tilde{h}_k$ back to $h_j$.  This is
 * done by shuffling the data, performing the reverse DFT, and multiplying by
 * $\Delta f$.
 * 
 * \subsubsection*{Operating Instructions}
 *
 * \begin{verbatim}
 * const UINT4 n  = 65536;
 * const REAL4 dt = 1.0 / 16384.0;
 * static LALStatus status;
 * static REAL4TimeSeries         x;
 * static COMPLEX8FrequencySeries X;
 * static COMPLEX8TimeSeries      z;
 * static COMPLEX8FrequencySeries Z;
 * RealFFTPlan    *fwdRealPlan    = NULL;
 * RealFFTPlan    *revRealPlan    = NULL;
 * ComplexFFTPlan *fwdComplexPlan = NULL;
 * ComplexFFTPlan *revComplexPlan = NULL;
 *
 * LALSCreateVector( &status, &x.data, n );
 * LALCCreateVector( &status, &X.data, n / 2 + 1 );
 * LALCCreateVector( &status, &z.data, n );
 * LALCCreateVector( &status, &Z.data, n );
 * LALCreateForwardRealFFTPlan( &status, &fwdRealPlan, n, 0 );
 * LALCreateReverseRealFFTPlan( &status, &revRealPlan, n, 0 );
 * LALCreateForwardComplexFFTPlan( &status, &fwdComplexPlan, n, 0 );
 * LALCreateReverseComplexFFTPlan( &status, &revComplexPlan, n, 0 );
 *
 * x.f0 = 0;
 * x.deltaT = dt;
 * x.sampleUnits = lalMeterUnit;
 * strncpy( x.name, "x", sizeof( x.name ) );
 *
 * z.f0 = 0;
 * z.deltaT = dt;
 * z.sampleUnits = lalVoltUnit;
 * strncpy( z.name, "z", sizeof( z.name ) );
 *
 * <assign data>
 *
 * LALTimeFreqRealFFT( &status, &X, &x, fwdRealPlan );
 * LALFreqTimeRealFFT( &status, &x, &X, revRealPlan );
 * LALTimeFreqComplexFFT( &status, &Z, &z, fwdComplexPlan );
 * LALFreqTimeComplexFFT( &status, &z, &Z, revComplexPlan );
 *
 * LALDestroyRealFFTPlan( &status, &fwdRealPlan );
 * LALDestroyRealFFTPlan( &status, &revRealPlan );
 * LALDestroyComplexFFTPlan( &status, &fwdComplexPlan );
 * LALDestroyComplexFFTPlan( &status, &revComplexPlan );
 * LALCDestroyVector( &status, &Z.data );
 * LALCDestroyVector( &status, &z.data );
 * LALCDestroyVector( &status, &X.data );
 * LALSDestroyVector( &status, &x.data );
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * \begin{enumerate}
 * \item The routines do not presently work properly with heterodyned data,
 * i.e., the original time series data should have \verb+f0+ equal to zero.
 * \end{enumerate}
 *
 * \vfill{\footnotesize\input{TimeFreqFFTCV}}
 * 
 **** </lalLaTeX> */


#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>

NRCSID( TIMEFREQFFTC, "$Id$" );


/* <lalVerbatim file="TimeFreqFFTCP"> */
void
LALTimeFreqRealFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *freq,
    REAL4TimeSeries         *time,
    RealFFTPlan             *plan
    )
{ /* </lalVerbatim> */
  LALUnitPair unitPair;
  UINT4 k;

  INITSTATUS( status, "LALTimeFreqRealFFT", TIMEFREQFFTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data->length, status,
      TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( time->deltaT > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );

  unitPair.unitOne = time->sampleUnits;
  unitPair.unitTwo = lalSecondUnit;

  /*
   *
   * The field f0 is not consistently defined between time and freq series.
   * Just do something....  Assume data is not heterodyned.
   *
   */
  if ( time->f0 != 0 )
  {
    LALWarning( status, "Frequency series may have incorrect f0." );
  }
  freq->f0     = 0; /* correct value for unheterodyned data */
  freq->epoch  = time->epoch;
  freq->deltaF = 1.0 / ( time->deltaT * time->data->length );
  TRY( LALUnitMultiply( status->statusPtr, &freq->sampleUnits, &unitPair ),
      status );
  TRY( LALForwardRealFFT( status->statusPtr, freq->data, time->data, plan ),
      status );
  for ( k = 0; k < freq->data->length; ++k )
  {
    freq->data->data[k].re *= time->deltaT;
    freq->data->data[k].im *= time->deltaT;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="TimeFreqFFTCP"> */
void
LALFreqTimeRealFFT(
    LALStatus               *status,
    REAL4TimeSeries         *time,
    COMPLEX8FrequencySeries *freq,
    RealFFTPlan             *plan
    )
{ /* </lalVerbatim> */
  LALUnitPair unitPair;
  UINT4 j;

  INITSTATUS( status, "LALFreqTimeRealFFT", TIMEFREQFFTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data->length, status,
      TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( freq->deltaF > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );

  unitPair.unitOne = freq->sampleUnits;
  unitPair.unitTwo = lalHertzUnit;

  /*
   *
   * The field f0 is not consistently defined between time and freq series.
   * Just do something....  Assume data is not heterodyned.
   *
   */
  if ( freq->f0 != 0 )
  {
    LALWarning( status, "Time series may have incorrect f0." );
  }
  time->f0     = 0; /* correct value for unheterodyned data */
  time->epoch  = freq->epoch;
  time->deltaT = 1.0 / ( freq->deltaF * time->data->length );
  TRY( LALUnitMultiply( status->statusPtr, &time->sampleUnits, &unitPair ),
      status );
  TRY( LALReverseRealFFT( status->statusPtr, time->data, freq->data, plan ),
      status );
  for ( j = 0; j < time->data->length; ++j )
  {
    time->data->data[j] *= freq->deltaF;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/*
 *
 * We need to define the ComplexFFTPlan structure so that we can check the FFT
 * sign.  The plan is not really a void*, but it is a pointer so a void* is
 * good enough (we don't need it).
 *
 */
struct tagComplexFFTPlan
{
  INT4  sign;
  UINT4 size;
  void *junk;
};


/* <lalVerbatim file="TimeFreqFFTCP"> */
void
LALTimeFreqComplexFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *freq,
    COMPLEX8TimeSeries      *time,
    ComplexFFTPlan          *plan
    )
{ /* </lalVerbatim> */
  COMPLEX8Vector *tmp = NULL;
  LALUnitPair unitPair;
  UINT4 n;
  UINT4 k;

  INITSTATUS( status, "LALTimeFreqComplexFFT", TIMEFREQFFTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  n = time->data->length;
  ASSERT( n, status, TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( time->deltaT > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );
  ASSERT( plan->sign == -1, status, TIMEFREQFFTH_ESIGN, TIMEFREQFFTH_MSGESIGN );

  unitPair.unitOne = time->sampleUnits;
  unitPair.unitTwo = lalSecondUnit;
  TRY( LALUnitMultiply( status->statusPtr, &freq->sampleUnits, &unitPair ),
      status );

  freq->epoch  = time->epoch;
  freq->deltaF = 1.0 / ( time->deltaT * n );
  freq->f0     = time->f0 - freq->deltaF * floor( n / 2 );

  TRY( LALUnitMultiply( status->statusPtr, &freq->sampleUnits, &unitPair ),
      status );

  TRY( LALCCreateVector( status->statusPtr, &tmp, n ), status );

  LALCOMPLEX8VectorFFT( status->statusPtr, tmp, time->data, plan );
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, &tmp ), status );
  }
  ENDFAIL( status );

  /*
   *
   * Unpack the frequency series and multiply by deltaT.
   *
   */
  for ( k = 0; k < n / 2; ++k )
  {
    UINT4 kk = k + ( n + 1 ) / 2;
    freq->data->data[k].re = time->deltaT * tmp->data[kk].re;
    freq->data->data[k].im = time->deltaT * tmp->data[kk].im;
  }
  for ( k = n / 2; k < n; ++k )
  {
    UINT4 kk = k - n / 2;
    freq->data->data[k].re = time->deltaT * tmp->data[kk].re;
    freq->data->data[k].im = time->deltaT * tmp->data[kk].im;
  }

  TRY( LALCDestroyVector( status->statusPtr, &tmp ), status );
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="TimeFreqFFTCP"> */
void
LALFreqTimeComplexFFT(
    LALStatus               *status,
    COMPLEX8TimeSeries      *time,
    COMPLEX8FrequencySeries *freq,
    ComplexFFTPlan          *plan
    )
{ /* </lalVerbatim> */
  COMPLEX8Vector *tmp = NULL;
  LALUnitPair unitPair;
  UINT4 n;
  UINT4 k;

  INITSTATUS( status, "LALFreqTimeComplexFFT", TIMEFREQFFTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  n = freq->data->length;
  ASSERT( n, status, TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( freq->deltaF > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );
  ASSERT( plan->sign == 1, status, TIMEFREQFFTH_ESIGN, TIMEFREQFFTH_MSGESIGN );

  unitPair.unitOne = freq->sampleUnits;
  unitPair.unitTwo = lalHertzUnit;
  TRY( LALUnitMultiply( status->statusPtr, &time->sampleUnits, &unitPair ),
      status );

  time->f0     = freq->f0 + freq->deltaF * floor( n / 2 );
  time->epoch  = freq->epoch;
  time->deltaT = 1.0 / ( freq->deltaF * n );

  TRY( LALUnitMultiply( status->statusPtr, &freq->sampleUnits, &unitPair ),
      status );

  TRY( LALCCreateVector( status->statusPtr, &tmp, n ), status );

  /*
   *
   * Pack the frequency series and multiply by deltaF.
   *
   */
  for ( k = 0; k < n / 2; ++k )
  {
    UINT4 kk = k + ( n + 1 ) / 2;
    tmp->data[kk].re = freq->deltaF * freq->data->data[k].re;
    tmp->data[kk].im = freq->deltaF * freq->data->data[k].im;
  }
  for ( k = n / 2; k < n; ++k )
  {
    UINT4 kk = k - n / 2;
    tmp->data[kk].re = freq->deltaF * freq->data->data[k].re;
    tmp->data[kk].im = freq->deltaF * freq->data->data[k].im;
  }

  LALCOMPLEX8VectorFFT( status->statusPtr, time->data, tmp, plan );
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, &tmp ), status );
  }
  ENDFAIL( status );

  TRY( LALCDestroyVector( status->statusPtr, &tmp ), status );
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
