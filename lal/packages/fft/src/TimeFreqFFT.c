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
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALConstants.h>

NRCSID( TIMEFREQFFTC, "$Id$" );


/*************************************************************/
/* calculation of median for power spectrum data 7/17/02 MSW */
/* implements an insert sort using temporary array           */
/* assumes inputs are positive, sorts in descending order    */
/*************************************************************/

static REAL4 EPMedian(REAL4 *p, INT4 j, INT4 flength, INT4 numSegs, LALStatus *status)
{
    /* p points to array of power spectra data over time slices */
    /* j is desired frequency offset into power spectra array   */
    /* flength is size of frequency series obtained from DFT    */
    /* numSegs is the number of time slices to be evaluated     */
    /* status points to LAL status struct passed into main      */
    /* returns the median value, over time slice at given freq. */

    /* TODO should check for valid or reasonable inputs */
    INT4  outer  = 0;       /* local loop counter */
    INT4  middle = 0;       /* local loop counter */
    INT4  inner  = 0;       /* local loop counter */
    REAL4 returnVal = 0.0;  /* holder for return value */

    /* allocate memory array for insert sort, test for success */
    REAL4 *s = (REAL4 *)LALMalloc(numSegs * sizeof(REAL4));
    ASSERT(s, status, TIMEFREQFFTH_EMALLOC, TIMEFREQFFTH_MSGEMALLOC);

    /* zero out the sort array */
    for (outer = 0; outer < numSegs; outer++)
    {
      s[outer] = 0.0;
    }

    /* scan time slices for a given frequency */
    for (outer = 0; outer < numSegs; outer++)
    {
      /* insert power value into sort array */
      REAL4 tmp = p[outer * flength + j]; /* obtain value to insert */
      for (middle = 0; middle < numSegs; middle++)
      {
        if (tmp > s[middle])
        {
          /* insert taking place of s[middle] */
          for (inner = numSegs - 1; inner > middle; inner--)
          {
            s[inner] = s [inner - 1];  /* move old values */
          }
          s[middle] = tmp;   /* insert new value */
          break;  /* terminate loop */
        }
      }
    }  /* done inserting into sort array */

    /* check for odd or even number of segments */
    if (numSegs % 2)
    {
      /* if odd number of segments, return median */
      returnVal = s[numSegs / 2];
    }
    else
    {
      /* if even number of segments, return average of two medians */
      returnVal = 0.5 * (s[numSegs/2] + s[(numSegs/2) - 1]);
    }

    /* free memory used for sort array */
    LALFree(s);

    return returnVal;
}

/******** <lalVerbatim file="CreateRealDFTParamsCP"> ********/
void
LALCreateRealDFTParams ( 
                     LALStatus                         *status, 
                     RealDFTParams                  **dftParams, 
                     LALWindowParams                   *params,
                     INT2                           sign
		     )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALCreateRealDFTParams", TIMEFREQFFTC);
  ATTATCHSTATUSPTR (status);

  /* 
   * Check return structure: dftParams should point to a valid pointer
   * which should not yet point to anything.
   *
   */
  ASSERT (dftParams, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL); 
  ASSERT (*dftParams == NULL, status, TIMEFREQFFTH_EALLOC, 
          TIMEFREQFFTH_MSGEALLOC);


  ASSERT (params, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL);  

  ASSERT (params->length > 0, status, TIMEFREQFFTH_EPOSARG, 
          TIMEFREQFFTH_MSGEPOSARG);

  ASSERT( (sign==1) || (sign==-1), status, TIMEFREQFFTH_EINCOMP,
          TIMEFREQFFTH_MSGEINCOMP);

  /*  Assign memory for *dftParams and check allocation */
  if ( !( *dftParams = (RealDFTParams *) LALMalloc(sizeof(RealDFTParams)) ) ){
    ABORT (status, TIMEFREQFFTH_EMALLOC, TIMEFREQFFTH_MSGEMALLOC);
  }
  
  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  if(sign==1)
    {
      /* Estimate the FFT plan */
      LALCreateForwardRealFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
    }
  else
    {
      /* Estimate the FFT plan */
      LALCreateReverseRealFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
    }
  CHECKSTATUSPTR (status);

  LALSCreateVector (status->statusPtr, &((*dftParams)->window), params->length);
  CHECKSTATUSPTR (status);

  LALWindow (status->statusPtr, ((*dftParams)->window), params);
  CHECKSTATUSPTR (status);

  (*dftParams)->sumofsquares = params->sumofsquares;
  (*dftParams)->windowType = params->type;
  
  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/******** <lalVerbatim file="DestroyRealDFTParamsCP"> ********/
void
LALDestroyRealDFTParams (
		      LALStatus                 *status, 
		      RealDFTParams          **dftParams
		      )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALDestroyRealDFTParams", TIMEFREQFFTC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not null */
  ASSERT (dftParams, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL);
  ASSERT (*dftParams, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL);

  /* make sure that data pointed to is non-null */
  ASSERT ((*dftParams)->plan, status, TIMEFREQFFTH_ENULL, 
          TIMEFREQFFTH_MSGENULL); 
  ASSERT ((*dftParams)->window, status, TIMEFREQFFTH_ENULL, 
          TIMEFREQFFTH_MSGENULL); 

  /* Ok, now let's free allocated storage */
  LALSDestroyVector (status->statusPtr, &((*dftParams)->window));
  CHECKSTATUSPTR (status);
  LALDestroyRealFFTPlan (status->statusPtr, &((*dftParams)->plan));
  CHECKSTATUSPTR (status);
  LALFree ( *dftParams );      /* free DFT parameters structure itself */

  *dftParams = NULL;	       /* make sure we don't point to freed struct */

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



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

  unitPair.unitOne = &(time->sampleUnits);
  unitPair.unitTwo = &(lalSecondUnit);

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

  unitPair.unitOne = &(freq->sampleUnits);
  unitPair.unitTwo = &(lalHertzUnit);

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
 * THIS IS THE AVERAGE SPECTRUM ROUTINE
 *
 */

void
LALRealAverageSpectrum(
    LALStatus            *status,
    REAL4FrequencySeries *fSeries,
    REAL4TimeSeries      *tSeries,
    RealDFTParams        *params,
    AvgSpecMethod         method
    )
{
    INT4   numSegs=0;                        /* number of overlapping segments */
    INT4   flength = fSeries->data->length;  /* length of output freq. series  */
    REAL4Vector     *dummyVec=NULL;
    COMPLEX8Vector  *dummyCVec=NULL;
    REAL4TimeSeries *dummyTimeSeries=NULL;
    INT4 n,j,i;
    REAL4 fac;

    INITSTATUS (status, "LALRealAverageSpectrum", TIMEFREQFFTC);
    ATTATCHSTATUSPTR (status);

    /* new convenience variables */
    numSegs = floor( tSeries->data->length / ( flength-1 )  - 1);
    n = params->window->length;
    fac = 1 / sqrt( params->sumofsquares);

    /* copy over data into frequency series structure */
    fSeries->epoch = tSeries->epoch;
    fSeries->f0    = tSeries->f0;
    fSeries->deltaF = 1/((REAL8)(n)*tSeries->deltaT); 

    /* create temporary vector */
    LALCreateVector (status->statusPtr, &dummyVec, n);
    CHECKSTATUSPTR (status);
    LALCCreateVector (status->statusPtr, &dummyCVec, flength);
    CHECKSTATUSPTR (status);

    
    if (method == useMedian)
    { 
        /* allocate array to hold power at each freq in each time slice */
        REAL4 *ptr = (REAL4 *)LALMalloc(flength * numSegs * sizeof(REAL4));
        ASSERT(ptr, status, TIMEFREQFFTH_EMALLOC, TIMEFREQFFTH_MSGEMALLOC);

        for (j = 0; j < flength; j++){
            for (i = 0; i < numSegs; i++){
                ptr[i * flength + j] = 0.0;
            }
        }

        /* obtain power spectrum in each overlapping time slice */
        for (i = 0; i < numSegs; i++)
        {
            INT4 offset=i*(flength-1);  /* This overlap = 0.5*n */

            /* compute windowed version of time series data */
            for (j = 0; j < n; j++)
            {
                dummyVec->data[j] = fac * tSeries->data->data[j+offset]
                    * params->window->data[j];
            };

            /* compute the DFT */
            LALForwardRealFFT (status->statusPtr, dummyCVec, dummyVec, params->plan);
            CHECKSTATUSPTR (status);

            /* copy modulus of DFT output into two dimensional array */
            for (j=0 ; j < flength ; j++)
            {
                REAL4 redummy = dummyCVec->data[j].re ;
                REAL4 imdummy = dummyCVec->data[j].im ;
                ptr[i * flength + j] = redummy*redummy + imdummy*imdummy;
            }
        }  /* done computing spectrum */

        /* find median value over time slices for each frequency */
        for (j = 0; j < flength; j++)
        {
            fSeries->data->data[j] = EPMedian(ptr, j, flength, numSegs, status);
            /* factors needed to give PSD of conventions doc T010095
             * and to rescale Gaussian white noise to mean estimator */
            fSeries->data->data[j] *= (2.0 * tSeries->deltaT / LAL_LN2); 
        }

        LALFree(ptr);     
    }  /* end of new code using median */

    else if (method == useMean)
    {  
        /* allocate array to hold power at each freq in each time slice */
        REAL4 *ptr = (REAL4 *)LALMalloc(flength * sizeof(REAL4));
        ASSERT(ptr, status, TIMEFREQFFTH_EMALLOC, TIMEFREQFFTH_MSGEMALLOC);

        /* zero out memory array */
        for (j = 0; j < flength; j++){
            ptr[j] = 0.0;
        }

        /* obtain power spectrum in each time slice */
        for (i = 0; i < numSegs; i++)
        {
            INT4 offset=i*(flength-1);

            /* compute windowed version of time series data */
            for (j = 0; j < n; j++)
            {
                dummyVec->data[j] = fac * tSeries->data->data[j+offset]
                    * params->window->data[j];
            };

            /* compute the DFT */
            LALForwardRealFFT (status->statusPtr, dummyCVec, dummyVec, params->plan);
            CHECKSTATUSPTR (status);

            /* copy modulus of DFT output into two dimensional array */
            for (j=0 ; j < flength ; j++)
            {
                REAL4 redummy = dummyCVec->data[j].re ;
                REAL4 imdummy = dummyCVec->data[j].im ;
                ptr[j] += redummy*redummy + imdummy*imdummy;
            }
        }  /* done computing spectrum */

        /* find median value over time slices for each frequency */
        for (j = 0; j < flength; j++)
        {
            /* compute normalized power spectrum consistent with PSD 
             * of conventions doc T010095 */
            fSeries->data->data[j] = ptr[j] * (2.0 * tSeries->deltaT) /
                (REAL4)(numSegs);
        }

        LALFree(ptr);     

    }  /* end of original code using mean */

    /* default case for unknown method */
    else
    {
        ABORT (status, TIMEFREQFFTH_EINCOMP, TIMEFREQFFTH_MSGEINCOMP );
    }

    LALDestroyVector (status->statusPtr, &dummyVec);
    CHECKSTATUSPTR (status);
    LALCDestroyVector (status->statusPtr, &dummyCVec);
    CHECKSTATUSPTR (status);


    /* normal exit */
    DETATCHSTATUSPTR (status);
    RETURN (status);
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

  unitPair.unitOne = &(time->sampleUnits);
  unitPair.unitTwo = &(lalSecondUnit);
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

  unitPair.unitOne = &(freq->sampleUnits);
  unitPair.unitTwo = &(lalHertzUnit);
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





