#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>

/*
 *
 * Compute a "modified periodogram," i.e., the power spectrum of a windowed
 * time series.
 *
 */
int XLALREAL4ModifiedPeriodogram(
    REAL4FrequencySeries        *periodogram,
    REAL4TimeSeries             *tseries,
    REAL4Window                 *window,
    REAL4FFTPlan                *plan
    )
{
  static const char *func = "XLALREAL4ModifiedPeriodogram";
  REAL4Vector *work;
  REAL4 normfac;
  UINT4 k;

  if ( ! periodogram || ! tseries || ! plan )
      XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! periodogram->data || ! tseries->data )
      XLAL_ERROR( func, XLAL_EINVAL );
  if ( tseries->deltaT <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );
  if ( periodogram->data->length != tseries->data->length/2 + 1 )
      XLAL_ERROR( func, XLAL_EBADLEN );
  if ( window )
  {
    if ( ! window->data )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->sumofsquares <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->data->length != tseries->data->length )
      XLAL_ERROR( func, XLAL_EBADLEN );
  }

  /* if the window has been specified, apply it to data */
  if ( window )
  {
    UINT4 j;

    work = XLALCreateREAL4Vector( tseries->data->length );
    if ( ! work )
      XLAL_ERROR( func, XLAL_EFUNC );

    /* apply windowing to data */
    for ( j = 0; j < tseries->data->length; ++j )
      work->data[j] = tseries->data->data[j] * window->data->data[j];
  }
  else /* otherwise just set work to the timeseries data */
  {
    work = tseries->data;
  }

  /* compute the power spectrum of the windowed timeseries */
  /* CHECKME: are DC and Nyquist right? */
  if ( XLALREAL4PowerSpectrum( periodogram->data, work, plan ) == XLAL_FAILURE )
  {
    if ( window ) /* need to free the workspace */
    {
      int saveErrno = xlalErrno;
      xlalErrno = 0;
      XLALDestroyREAL4Vector( work );
      xlalErrno = saveErrno;
    }
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* destroy the workspace if it was created */
  if ( window )
  {
    XLALDestroyREAL4Vector( work );
    if ( xlalErrno )
      XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* normalize power spectrum to give correct units */
  /* CHECKME: is this the right factor? */
  normfac = tseries->deltaT;
  if ( window )
    normfac /= window->sumofsquares;
  else
    normfac /= tseries->data->length;

  for ( k = 0; k < periodogram->data->length; ++k )
    periodogram->data->data[k] *= normfac;

  /* now set rest of metadata */
  periodogram->epoch  = tseries->epoch;
  periodogram->f0     = tseries->f0; /* FIXME: is this right? */
  periodogram->deltaF = 1.0 / ( tseries->data->length * tseries->deltaT );

  /* compute units */
  if ( ! XLALUnitSquare( &periodogram->sampleUnits, &tseries->sampleUnits ) )
    XLAL_ERROR( func, XLAL_EFUNC );
  if ( !  XLALUnitMultiply( &periodogram->sampleUnits,
        &periodogram->sampleUnits, &lalSecondUnit ) )
    XLAL_ERROR( func, XLAL_EFUNC );

  return 0;
}


/* 
 *
 * Use Welch's method to compute the average power spectrum of a time series.
 *
 * See: Peter D. Welch "The Use of Fast Fourier Transform for the Estimation
 * of Power Spectra: A Method Based on Time Averaging Over Short, Modified
 * Periodograms"  IEEE Transactions on Audio and Electroacoustics,
 * Vol. AU-15, No. 2, June 1967.
 *
 */
int XLALREAL4AverageSpectrumWelch(
    REAL4FrequencySeries        *spectrum,
    REAL4TimeSeries             *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    REAL4Window                 *window,
    REAL4FFTPlan                *plan
    )
{
  static const char *func = "XLALREAL4AverageSpectrumWelch";
  REAL4FrequencySeries work; /* workspace */
  UINT4 reclen; /* length of entire data record */
  UINT4 numseg;
  UINT4 seg;
  UINT4 k;

  if ( ! spectrum || ! tseries || ! plan )
      XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! spectrum->data || ! tseries->data )
      XLAL_ERROR( func, XLAL_EINVAL );
  if ( tseries->deltaT <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );

  reclen = tseries->data->length;
  numseg = 1 + (reclen - seglen)/stride;

  /* consistency check for lengths: make sure that the segments cover the
   * data record completely */
  if ( (numseg - 1)*stride + seglen != reclen )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( spectrum->data->length != seglen/2 + 1 )
      XLAL_ERROR( func, XLAL_EBADLEN );

  /* make sure window, if present, is appropriate */
  if ( window )
  {
    if ( ! window->data )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->sumofsquares <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->data->length != seglen )
      XLAL_ERROR( func, XLAL_EBADLEN );
  }

  /* clear spectrum data */
  memset( spectrum->data->data, 0,
      spectrum->data->length * sizeof( *spectrum->data->data ) );

  /* create frequency series data workspace */
  work.data = XLALCreateREAL4Vector( spectrum->data->length );
  if ( ! work.data )
    XLAL_ERROR( func, XLAL_EFUNC );

  for ( seg = 0; seg < numseg; ++seg )
  {
    REAL4Vector savevec; /* save the time series data vector */
    int code;

    /* save the time series data vector */
    savevec = *tseries->data;

    /* set the data vector to be appropriate for this segment */
    tseries->data->length  = seglen;
    tseries->data->data   += seg * stride;

    /* compute the modified periodogram */
    code = XLALREAL4ModifiedPeriodogram( &work, tseries, window, plan );
    
    /* restore the time series data vector to its original state */
    *tseries->data = savevec;

    /* now check for failure of the XLAL routine */
    if ( code == XLAL_FAILURE )
    {
      /* cleanup workspace data first */
      int saveErrno = xlalErrno;
      xlalErrno = 0;
      XLALDestroyREAL4Vector( work.data );
      xlalErrno = saveErrno;
      XLAL_ERROR( func, XLAL_EFUNC );
    }

    /* add the periodogram to the running sum */
    for ( k = 0; k < spectrum->data->length; ++k )
      spectrum->data->data[k] += work.data->data[k];
  }

  /* free the workspace data */
  XLALDestroyREAL4Vector( work.data );
  if ( xlalErrno )
    XLAL_ERROR( func, XLAL_EFUNC );

  /* set metadata */
  spectrum->epoch       = work.epoch;
  spectrum->f0          = work.f0;
  spectrum->deltaF      = work.deltaF;
  spectrum->sampleUnits = work.sampleUnits;

  /* divide spectrum data by the number of segments in average */
  for ( k = 0; k < spectrum->data->length; ++k )
    spectrum->data->data[k] /= numseg;

  return 0;
}


/* 
 *
 * Median Method: use median average rather than mean.  Note: this will
 * cause a bias if the segments overlap, i.e., if the stride is less than
 * the segment length -- even though the median bias for Gaussian noise
 * is accounted for -- because the segments are not independent and their
 * correlation is non-zero.
 *
 */

/* compute the median bias */
REAL8 XLALMedianBias( UINT4 nn )
{
  const UINT4 nmax = 1000;
  REAL8 ans = 1;
  UINT4 n = (nn - 1)/2;
  UINT4 i;

  if ( nn >= nmax )
    return LAL_LN2;

  for ( i = 1; i <= n; ++i )
  {
    ans -= 1.0/(2*i);
    ans += 1.0/(2*i + 1);
  }

  return ans;
}

/* cleanup temporary workspace... ignore xlal errors */
static void median_cleanup( REAL4FrequencySeries *work, UINT4 n )
{
  int saveErrno = xlalErrno;
  UINT4 i;
  for ( i = 0; i < n; ++i )
    if ( work[i].data )
      XLALDestroyREAL4Vector( work[i].data );
  LALFree( work );
  xlalErrno = saveErrno;
  return;
}

/* comparison for floating point numbers */
static int compare_float( const void *p1, const void *p2 )
{
  REAL4 x1 = * ( ( const REAL4 * )p1 );
  REAL4 x2 = * ( ( const REAL4 * )p2 );
  if ( x1 < x2 )
    return -1;
  if ( x1 > x2 )
    return 1;
  return 0;
}


/* here is the median method */
int XLALREAL4AverageSpectrumMedian(
    REAL4FrequencySeries        *spectrum,
    REAL4TimeSeries             *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    REAL4Window                 *window,
    REAL4FFTPlan                *plan
    )
{
  static const char *func = "XLALREAL4AverageSpectrumMedian";
  REAL4FrequencySeries *work; /* array of frequency series */
  REAL4 *bin; /* array of bin values */
  REAL4 biasfac; /* median bias factor */
  REAL4 normfac; /* normalization factor */
  UINT4 reclen; /* length of entire data record */
  UINT4 numseg;
  UINT4 seg;
  UINT4 k;

  if ( ! spectrum || ! tseries || ! plan )
      XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! spectrum->data || ! tseries->data )
      XLAL_ERROR( func, XLAL_EINVAL );
  if ( tseries->deltaT <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );

  reclen = tseries->data->length;
  numseg = 1 + (reclen - seglen)/stride;

  /* consistency check for lengths: make sure that the segments cover the
   * data record completely */
  if ( (numseg - 1)*stride + seglen != reclen )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( spectrum->data->length != seglen/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* make sure window, if present, is appropriate */
  if ( window )
  {
    if ( ! window->data )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->sumofsquares <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->data->length != seglen )
      XLAL_ERROR( func, XLAL_EBADLEN );
  }

  /* create frequency series data workspaces */
  work = LALCalloc( numseg, sizeof( *work ) );
  if ( ! work )
    XLAL_ERROR( func, XLAL_ENOMEM );
  for ( seg = 0; seg < numseg; ++seg )
  {
    work[seg].data = XLALCreateREAL4Vector( spectrum->data->length );
    if ( ! work[seg].data )
    {
      median_cleanup( work, numseg ); /* cleanup */
      XLAL_ERROR( func, XLAL_EFUNC );
    }
  }

  for ( seg = 0; seg < numseg; ++seg )
  {
    REAL4Vector savevec; /* save the time series data vector */
    int code;

    /* save the time series data vector */
    savevec = *tseries->data;

    /* set the data vector to be appropriate for the even segment */
    tseries->data->length  = seglen;
    tseries->data->data   += seg * stride;

    /* compute the modified periodogram for the even segment */
    code = XLALREAL4ModifiedPeriodogram( work + seg, tseries, window, plan );
    
    /* restore the time series data vector to its original state */
    *tseries->data = savevec;

    /* now check for failure of the XLAL routine */
    if ( code == XLAL_FAILURE )
    {
      median_cleanup( work, numseg ); /* cleanup */
      XLAL_ERROR( func, XLAL_EFUNC );
    }
  }

  /* create array to hold a particular frequency bin data */
  bin = LALMalloc( numseg * sizeof( *bin ) );
  if ( ! bin )
  {
    median_cleanup( work, numseg ); /* cleanup */
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* compute median bias factor */
  biasfac = XLALMedianBias( numseg );

  /* normaliztion takes into account bias */
  normfac = 1.0 / biasfac;

  /* now loop over frequency bins and compute the median-mean */
  for ( k = 0; k < spectrum->data->length; ++k )
  {
    /* assign array of even segment values to bin array for this freq bin */
    for ( seg = 0; seg < numseg; ++seg )
      bin[seg] = work[seg].data->data[k];

    /* sort them and find median */
    qsort( bin, numseg, sizeof( *bin ), compare_float );
    if ( numseg % 2 ) /* odd number of evens */
      spectrum->data->data[k] = bin[(numseg+1)/2];
    else /* even number... take average */
      spectrum->data->data[k] = 0.5*(bin[numseg/2-1] + bin[numseg/2]);

    /* remove median bias */
    spectrum->data->data[k] *= normfac;
  }

  /* set metadata */
  spectrum->epoch       = work->epoch;
  spectrum->f0          = work->f0;
  spectrum->deltaF      = work->deltaF;
  spectrum->sampleUnits = work->sampleUnits;

  /* free the workspace data */
  LALFree( bin );
  median_cleanup( work, numseg );

  return 0;
}


/* 
 *
 * Median-Mean Method: divide overlapping segments into "even" and "odd"
 * segments; compute the bin-by-bin median of the "even" segments and the
 * "odd" segments, and then take the bin-by-bin average of these two median
 * averages.
 *
 */


/* cleanup temporary workspace... ignore xlal errors */
static void median_mean_cleanup( REAL4FrequencySeries *even,
    REAL4FrequencySeries *odd, UINT4 n )
{
  int saveErrno = xlalErrno;
  UINT4 i;
  for ( i = 0; i < n; ++i )
  {
    if ( even[i].data )
      XLALDestroyREAL4Vector( even[i].data );
    if ( odd[i].data )
      XLALDestroyREAL4Vector( odd[i].data );
  }
  LALFree( even );
  LALFree( odd );
  xlalErrno = saveErrno;
  return;
}

/* here is the median-mean-method */
int XLALREAL4AverageSpectrumMedianMean(
    REAL4FrequencySeries        *spectrum,
    REAL4TimeSeries             *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    REAL4Window                 *window,
    REAL4FFTPlan                *plan
    )
{
  static const char *func = "XLALREAL4AverageSpectrumMedianMean";
  REAL4FrequencySeries *even; /* array of even frequency series */
  REAL4FrequencySeries *odd;  /* array of odd frequency series */
  REAL4 *bin; /* array of bin values */
  REAL4 biasfac; /* median bias factor */
  REAL4 normfac; /* normalization factor */
  UINT4 reclen; /* length of entire data record */
  UINT4 numseg;
  UINT4 halfnumseg;
  UINT4 seg;
  UINT4 k;

  if ( ! spectrum || ! tseries || ! plan )
      XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! spectrum->data || ! tseries->data )
      XLAL_ERROR( func, XLAL_EINVAL );
  if ( tseries->deltaT <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );

  reclen = tseries->data->length;
  numseg = 1 + (reclen - seglen)/stride;

  /* consistency check for lengths: make sure that the segments cover the
   * data record completely */
  if ( (numseg - 1)*stride + seglen != reclen )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( spectrum->data->length != seglen/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* for median-mean to work, the number of segments must be even and
   * the stride must be greater-than or equal-to half of the seglen */
  halfnumseg = numseg/2;
  if ( numseg%2 || stride < seglen/2 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* make sure window, if present, is appropriate */
  if ( window )
  {
    if ( ! window->data )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->sumofsquares <= 0.0 )
      XLAL_ERROR( func, XLAL_EINVAL );
    if ( window->data->length != seglen )
      XLAL_ERROR( func, XLAL_EBADLEN );
  }

  /* create frequency series data workspaces */
  even = LALCalloc( halfnumseg, sizeof( *even ) );
  if ( ! even )
    XLAL_ERROR( func, XLAL_ENOMEM );
  odd = LALCalloc( halfnumseg, sizeof( *odd ) );
  if ( ! odd )
    XLAL_ERROR( func, XLAL_ENOMEM );
  for ( seg = 0; seg < halfnumseg; ++seg )
  {
    even[seg].data = XLALCreateREAL4Vector( spectrum->data->length );
    odd[seg].data  = XLALCreateREAL4Vector( spectrum->data->length );
    if ( ! even[seg].data || ! odd[seg].data )
    {
      median_mean_cleanup( even, odd, halfnumseg ); /* cleanup */
      XLAL_ERROR( func, XLAL_EFUNC );
    }
  }

  for ( seg = 0; seg < halfnumseg; ++seg )
  {
    REAL4Vector savevec; /* save the time series data vector */
    int code;

    /* save the time series data vector */
    savevec = *tseries->data;

    /* set the data vector to be appropriate for the even segment */
    tseries->data->length  = seglen;
    tseries->data->data   += 2 * seg * stride;

    /* compute the modified periodogram for the even segment */
    code = XLALREAL4ModifiedPeriodogram( even + seg, tseries, window, plan );
    
    /* now check for failure of the XLAL routine */
    if ( code == XLAL_FAILURE )
    {
      *tseries->data = savevec;
      median_mean_cleanup( even, odd, halfnumseg ); /* cleanup */
      XLAL_ERROR( func, XLAL_EFUNC );
    }

    /* set the data vector to be appropriate for the odd segment */
    tseries->data->data += stride;

    /* compute the modified periodogram for the odd segment */
    code = XLALREAL4ModifiedPeriodogram( odd + seg, tseries, window, plan );

    /* now check for failure of the XLAL routine */
    if ( code == XLAL_FAILURE )
    {
      *tseries->data = savevec;
      median_mean_cleanup( even, odd, halfnumseg ); /* cleanup */
      XLAL_ERROR( func, XLAL_EFUNC );
    }

    /* restore the time series data vector to its original state */
    *tseries->data = savevec;
  }

  /* create array to hold a particular frequency bin data */
  bin = LALMalloc( halfnumseg * sizeof( *bin ) );
  if ( ! bin )
  {
    median_mean_cleanup( even, odd, halfnumseg ); /* cleanup */
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* compute median bias factor */
  biasfac = XLALMedianBias( halfnumseg );

  /* normaliztion takes into account bias and a factor of two from averaging
   * the even and the odd */
  normfac = 1.0 / ( 2.0 * biasfac );

  /* now loop over frequency bins and compute the median-mean */
  for ( k = 0; k < spectrum->data->length; ++k )
  {
    REAL4 evenmedian;
    REAL4 oddmedian;

    /* assign array of even segment values to bin array for this freq bin */
    for ( seg = 0; seg < halfnumseg; ++seg )
      bin[seg] = even[seg].data->data[k];

    /* sort them and find median */
    qsort( bin, halfnumseg, sizeof( *bin ), compare_float );
    if ( halfnumseg % 2 ) /* odd number of evens */
      evenmedian = bin[(halfnumseg+1)/2];
    else /* even number... take average */
      evenmedian = 0.5*(bin[halfnumseg/2-1] + bin[halfnumseg/2]);

    /* assign array of odd segment values to bin array for this freq bin */
    for ( seg = 0; seg < halfnumseg; ++seg )
      bin[seg] = odd[seg].data->data[k];

    /* sort them and find median */
    qsort( bin, halfnumseg, sizeof( *bin ), compare_float );
    if ( halfnumseg % 2 ) /* odd number of odds */
      oddmedian = bin[(halfnumseg+1)/2];
    else /* even number... take average */
      oddmedian = 0.5*(bin[halfnumseg/2-1] + bin[halfnumseg/2]);

    /* spectrum for this bin is the mean of the medians */
    spectrum->data->data[k] = normfac * (evenmedian + oddmedian);
  }

  /* set metadata */
  spectrum->epoch       = even->epoch;
  spectrum->f0          = even->f0;
  spectrum->deltaF      = even->deltaF;
  spectrum->sampleUnits = even->sampleUnits;

  /* free the workspace data */
  LALFree( bin );
  median_mean_cleanup( even, odd, halfnumseg );

  return 0;
}

int XLALREAL4SpectrumInvertTruncate(
    REAL4FrequencySeries        *spectrum,
    REAL4                        lowfreq,
    UINT4                        seglen,
    UINT4                        trunclen,
    REAL4FFTPlan                *fwdplan,
    REAL4FFTPlan                *revplan
    )
{
  static const char *func = "XLALREAL4SpectrumInvertTruncate";
  UINT4 cut;
  UINT4 k;

  if ( ! spectrum )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! spectrum->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( spectrum->deltaF <= 0.0 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( lowfreq < 0 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! seglen || spectrum->data->length != seglen/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( trunclen && ! revplan )
    XLAL_ERROR( func, XLAL_EFAULT );

  cut = lowfreq / spectrum->deltaF;
  if ( cut < 1 ) /* need to get rid of DC at least */
    cut = 1;

  if ( trunclen ) /* truncate while inverting */
  {
    COMPLEX8Vector *vtilde; /* frequency-domain vector workspace */
    REAL4Vector    *vector; /* time-domain vector workspace */
    REAL4           normfac;

    vector = XLALCreateREAL4Vector( seglen );
    vtilde = XLALCreateCOMPLEX8Vector( seglen / 2 + 1 );

    /* clear the low frequency components */
    memset( vtilde->data, 0, cut * sizeof( *vtilde->data ) );

    /* stick the root-inv-spectrum into the complex frequency-domain vector */
    for ( k = cut; k < vtilde->length - 1; ++k )
    {
      vtilde->data[k].re = 1.0 / sqrt( spectrum->data->data[k] );
      vtilde->data[k].im = 0.0; /* purely real */
    }

    /* no Nyquist */
    vtilde->data[vtilde->length - 1].re = 0.0;
    vtilde->data[vtilde->length - 1].im = 0.0;

    /* construct time-domain version of root-inv-spectrum */
    XLALREAL4ReverseFFT( vector, vtilde, revplan );

    /* now truncate it: zero time from trunclen/2 to length - trunclen/2 */
    memset( vector->data + trunclen/2, 0,
        (vector->length - trunclen) * sizeof( *vector->data ) );

    /* reconstruct spectrum */
    XLALREAL4PowerSpectrum( spectrum->data, vector, fwdplan );

    /* clear the low frequency components... again, and Nyquist too */
    memset( spectrum->data->data, 0, cut * sizeof( *spectrum->data->data ) );
    spectrum->data->data[spectrum->data->length - 1] = 0.0;

    /* rescale spectrum: 0.5 to undo factor of two from REAL4PowerSpectrum */
    /* seglen^2 for the reverse fft (squared because power spectrum) */
    /* Note: cast seglen to real so that seglen*seglen is a real rather */
    /* than a four-byte integer which might overflow... */
    normfac = 0.5 / ( ((REAL4)seglen) * ((REAL4)seglen) );
    for ( k = cut; k < spectrum->data->length - 1; ++k )
      spectrum->data->data[k] *= normfac;

    /* cleanup workspace */
    XLALDestroyCOMPLEX8Vector( vtilde );
    XLALDestroyREAL4Vector( vector );
  }
  else /* otherwise just invert the spectrum */
  {
    /* clear the low frequency components */
    memset( spectrum->data->data, 0, cut * sizeof( *spectrum->data->data ) );

    /* invert the high frequency (non-Nyquist) components */
    for ( k = cut; k < spectrum->data->length - 1; ++k )
      spectrum->data->data[k] = 1.0 / spectrum->data->data[k];

    /* zero Nyquist */
    spectrum->data->data[spectrum->data->length - 1] = 0.0;
  }

  /* now correct the units */
  XLALUnitRaiseINT2( &spectrum->sampleUnits, &spectrum->sampleUnits, -1 );

  return 0;
}
