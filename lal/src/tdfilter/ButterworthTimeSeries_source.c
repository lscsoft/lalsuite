#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define STRING(a) #a

#ifdef COMPLEX_DATA
#   define DBLDATATYPE COMPLEX16
#   ifdef SINGLE_PRECISION
#       define DATATYPE COMPLEX8
#   else
#       define DATATYPE COMPLEX16
#   endif
#else
#   define DBLDATATYPE REAL8
#   ifdef SINGLE_PRECISION
#       define DATATYPE REAL4
#   else
#       define DATATYPE REAL8
#   endif
#endif

#define SERIESTYPE CONCAT2(DATATYPE,TimeSeries)
#define VECTORTYPE CONCAT2(DATATYPE,Vector)
#define FILTERTYPE CONCAT2(DBLDATATYPE,IIRFilter)

#define BFUNC CONCAT2(XLALButterworth,SERIESTYPE)
#define LFUNC CONCAT2(XLALLowPass,SERIESTYPE)
#define HFUNC CONCAT2(XLALHighPass,SERIESTYPE)

#define CFUNC CONCAT2(XLALCreate,FILTERTYPE)
#define DFUNC CONCAT2(XLALDestroy,FILTERTYPE)

#define FFUNC CONCAT2(XLALIIRFilter,VECTORTYPE)
#define RFUNC CONCAT2(XLALIIRFilterReverse,VECTORTYPE)

int BFUNC(SERIESTYPE *series, PassBandParamStruc *params)
{
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  /* Make sure the input pointers are non-null. */
  if ( ! params || ! series || ! series->data || ! series->data->data )
    XLAL_ERROR( XLAL_EFAULT );

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=XLALParsePassBandParamStruc(params,&n,&wc,series->deltaT);
  if(type<0)
    XLAL_ERROR( XLAL_EINVAL );

  /* An order n Butterworth filter has n poles spaced evenly along a
     semicircle in the upper complex w-plane.  By pairing up poles
     symmetric across the imaginary axis, the filter gan be decomposed
     into [n/2] filters of order 2, plus perhaps an additional order 1
     filter.  The following loop pairs up poles and applies the
     filters with order 2. */
  for(i=0,j=n-1;i<j;i++,j--){
    REAL8 theta=LAL_PI*(i+0.5)/n;
    REAL8 ar=wc*cos(theta);
    REAL8 ai=wc*sin(theta);
    FILTERTYPE *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      zpgFilter = XLALCreateCOMPLEX16ZPGFilter(2,2);
      if ( ! zpgFilter )
        XLAL_ERROR( XLAL_EFUNC );
      zpgFilter->zeros->data[0]=0.0;
      zpgFilter->zeros->data[1]=0.0;
      zpgFilter->gain=1.0;
    }else{
      zpgFilter = XLALCreateCOMPLEX16ZPGFilter(0,2);
      if ( ! zpgFilter )
        XLAL_ERROR( XLAL_EFUNC );
      zpgFilter->gain=-wc*wc;
    }
    zpgFilter->poles->data[0]=ar;
    zpgFilter->poles->data[0]+=ai*I;
    zpgFilter->poles->data[1]=-ar;
    zpgFilter->poles->data[1]+=ai*I;

    /* Transform to the z-plane and create the IIR filter. */
    if (XLALWToZCOMPLEX16ZPGFilter(zpgFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }
    iirFilter = CFUNC(zpgFilter);
    if (!iirFilter)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Filter the data, once each way. */
    if (FFUNC(series->data,iirFilter)<0
        || RFUNC(series->data,iirFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      DFUNC(iirFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Free the filters. */
    XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
    DFUNC(iirFilter);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    FILTERTYPE *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      zpgFilter=XLALCreateCOMPLEX16ZPGFilter(1,1);
      if(!zpgFilter)
        XLAL_ERROR(XLAL_EFUNC);
      *zpgFilter->zeros->data=0.0;
      zpgFilter->gain=1.0;
    }else{
      zpgFilter=XLALCreateCOMPLEX16ZPGFilter(0,1);
      if(!zpgFilter)
        XLAL_ERROR(XLAL_EFUNC);
      zpgFilter->gain=-wc*I;
    }
    *zpgFilter->poles->data=wc*I;

    /* Transform to the z-plane and create the IIR filter. */
    if (XLALWToZCOMPLEX16ZPGFilter(zpgFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR(XLAL_EFUNC);
    }
    iirFilter=CFUNC(zpgFilter);
    if (!iirFilter)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR(XLAL_EFUNC);
    }

    /* Filter the data, once each way. */
    if (FFUNC(series->data,iirFilter)<0
        || RFUNC(series->data,iirFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      DFUNC(iirFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Free the filters. */
    XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
    DFUNC(iirFilter);
  }

  return 0;
}

int LFUNC(SERIESTYPE *series, REAL8 frequency, REAL8 amplitude, INT4 filtorder)
{
  PassBandParamStruc params;
  params.nMax = filtorder;
  params.f1   = frequency;
  params.a1   = amplitude;
  params.f2   = -1;
  params.a2   = -1;
  if (BFUNC(series, &params) < 0)
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

int HFUNC(SERIESTYPE *series, REAL8 frequency, REAL8 amplitude, INT4 filtorder)
{
  PassBandParamStruc params;
  params.nMax = filtorder;
  params.f2   = frequency;
  params.a2   = amplitude;
  params.f1   = -1;
  params.a1   = -1;
  if (BFUNC(series, &params) < 0)
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

#undef BFUNC
#undef LFUNC
#undef HFUNC
#undef CFUNC
#undef DFUNC
#undef FFUNC
#undef RFUNC
#undef SERIESTYPE
#undef VECTORTYPE
#undef FILTERTYPE
#undef DBLDATATYPE
#undef DATATYPE
#undef CONCAT2x
#undef CONCAT2
#undef STRING
