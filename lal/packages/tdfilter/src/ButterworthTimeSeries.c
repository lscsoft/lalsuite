/************************ <lalVerbatim file="ButterworthTimeSeriesCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{ButterworthTimeSeries.c}}
\label{ss:ButterworthTimeSeries.c}

Applies a low- or high-pass Butterworth filter to a time series.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ButterworthTimeSeriesD}
\idx{LALButterworthREAL4TimeSeries()}
\idx{LALButterworthREAL8TimeSeries()}
\idx{LALDButterworthREAL4TimeSeries()}

\subsubsection*{Description}

These routines perform an in-place time-domain band-pass filtering of
a data sequence \verb@*series@, using a Butterworth filter generated
from parameters \verb@*params@.  The routines construct a filter with
the square root of the desired amplitude response, which it then
applied to the data once forward and once in reverse.  This gives the
full amplitude response with little or no frequency-dependent phase
shift.

The routine \verb@LALDButterworthREAL4TimeSeries()@ applies a
double-precision filter to single-precision data, using
\verb@LALDIIRFilterREAL4Vector()@ and
\verb@LALDIIRFilterREAL4VectorR()@.

\subsubsection*{Algorithm}

The frequency response of a Butterworth low-pass filter is easiest to
express in terms of the transformed frequency variable $w=\tan(\pi
f\Delta t)$, where $\Delta t$ is the sampling interval (i.e.\
\verb@series->deltaT@).  In this parameter, then, the \emph{power}
response (attenuation) of the filter is:
$$
|R|^2 = \sqrt{a} = \frac{1}{1+(w/w_c)^{2n}} \; ,
$$
where $n$ is the filter order and $w_c$ is the characteristic
frequency.  We have written the attenuation as $\sqrt{a}$ to emphasize
that the full attenuation $a$ is achieved only after filtering twice
(once forward, once in reverse).  Similarly, a Butterworth high-pass
filter is given by
$$
|R|^2 = \sqrt{a} = \frac{1}{1+(w_c/w)^{2n}} \; .
$$
If one is given a filter order $n$, then the characteristic frequency
can be determined from the attenuation at some any given frequency.
Alternatively, $n$ and $w_c$ can both be computed given attenuations
at two different frequencies.

Frequencies in \verb@*params@ are assumed to be real frequencies $f$
given in the inverse of the units used for the sampling interval
\verb@series->deltaT@.  In order to be used, the pass band parameters
must lie in the ranges given below; if a parameter lies outside of its
range, then it is ignored and the filter is calculated from the
remaining parameters.  If too many parameters are missing, the routine
will fail.  The acceptable parameter ranges are:

\begin{description}
\item[\texttt{params->nMax}]   = 1, 2, $\ldots$
\item[\texttt{params->f1}, \texttt{f2}] $\in
  (0,\{2\times\verb@series->deltaT@\}^{-1}) $
\item[\texttt{params->a1}, \texttt{a2}] $\in (0,1) $
\end{description}

If both pairs of frequencies and amplitudes are given, then \verb@a1@,
\verb@a2@ specify the minimal requirements on the attenuation of the
filter at frequencies \verb@f1@, \verb@f2@.  Whether the filter is a
low- or high-pass filter is determined from the relative sizes of
these parameters.  In this case the \verb@nMax@ parameter is optional;
if given, it specifies an upper limit on the filter order.  If the
desired attenuations would require a higher order, then the routine
will sacrifice performance in the stop band in order to remain within
the specified \verb@nMax@.

If one of the frequency/attenuation pairs is missing, then the filter
is computed using the remaining pair and \verb@nMax@ (which must be
given).  The filter is taken to be a low-pass filter if \verb@f1@,
\verb@a1@ are given, and high-pass if \verb@f2@, \verb@a2@ are given.
If only one frequency and no corresponding attenuation is specified,
then it is taken to be the characteristic frequency (i.e. the
corresponding attenuation is assumed to be $\sqrt{a}=1/2$).  If none
of these conditions are met, the routine will return an error.

The following table summarizes the decision algorithm.  A $\bullet$
symbol indicates that the parameter is specified in the range given
above.  A $\circ$ symbol indicates that the parameter is ``not
given'', i.e.\ not specified in the valid range.
\begin{center}
\begin{tabular}{|cccccp{10cm}|}
\hline
\tt nMax & \tt f1 & \tt a1 & \tt f2 & \tt a2 & Procedure \\
\hline
 $\circ$  & $\bullet$ & $\bullet$ & $\bullet$ & $\bullet$ &
	Type of filter (low- or high-pass), $w_c$, and $n$ are
	computed from all four transition-band parameters. \\
$\bullet$ & $\bullet$ & $\bullet$ & $\bullet$ & $\bullet$ &
	Ditto, but if the resulting $n>$ \texttt{nMax}, $w_c$ is
	computed from \texttt{nMax} and the (\texttt{f},\texttt{a})
	pair with the \emph{larger} \texttt{a}. \\
$\bullet$ & $\bullet$ & $\bullet$ &  $\circ$  &  $\circ$  &
	Low-pass filter; $w_c$ is computed from \texttt{nMax},
	\texttt{f1}, and \texttt{a1}. \\
$\bullet$ & $\bullet$ & $\bullet$ &  $\circ$  & $\bullet$ &
	Ditto; \texttt{a2} is ignored. \\
$\bullet$ & $\bullet$ & $\bullet$ & $\bullet$ &  $\circ$  &
	Ditto; \texttt{f2} is ignored. \\
$\bullet$ & $\bullet$ &  $\circ$  &  $\circ$  &  $\circ$  &
	Low-pass filter; $w_c$ is computed as above with \texttt{a1}
	treated as 1/4. \\
$\bullet$ & $\bullet$ &  $\circ$  &  $\circ$  & $\bullet$ &
	Ditto; \texttt{a2} is ignored. \\
$\bullet$ &  $\circ$  &  $\circ$  & $\bullet$ & $\bullet$ &
	High-pass filter; $w_c$ is computed from \texttt{nMax},
	\texttt{f2}, and \texttt{a2}. \\
$\bullet$ &  $\circ$  & $\bullet$ & $\bullet$ & $\bullet$ &
	Ditto; \texttt{a1} is ignored. \\
$\bullet$ & $\bullet$ &  $\circ$  & $\bullet$ & $\bullet$ &
	Ditto; \texttt{f1} is ignored. \\
$\bullet$ &  $\circ$  &  $\circ$  & $\bullet$ &  $\circ$  &
	High-pass filter; $w_c$ is computed as above with \texttt{a2}
	treated as 1/4. \\
$\bullet$ &  $\circ$  & $\bullet$ & $\bullet$ &  $\circ$  &
	Ditto; \texttt{a1} is ignored. \\
\multicolumn{5}{|c}{Other} & Subroutine returns an error. \\
\hline
\end{tabular}
\end{center}

Once an order $n$ and characteristic frequency $w_c$ are known, the
zeros and poles of a ZPG filter are readily determined.  A stable,
physically realizable Butterworth filter will have $n$ poles evenly
spaced on the upper half of a circle of radius $w_c$; that is,
$$
R = \frac{(-iw_c)^n}{\prod_{k=0}^{n-1}(w - w_c e^{2\pi i(k+1/2)/n})}
$$
for a low-pass filter, and
$$
R = \frac{w^n}{\prod_{k=0}^{n-1}(w - w_c e^{2\pi i(k+1/2)/n})}
$$
for a high-pass filter.  By choosing only poles on the upper-half
plane, one ensures that after transforming to $z$ the poles will have
$|z|<1$.  Furthermore, the phase factor $(-i)^n$ in the numerator of
the low-pass filter is chosen so that the DC response is purely real;
this ensures that the response function in the $z$-plane will have a
real gain factor, and the resulting IIR filter will be physically
realizable.  The high-pass filter has a purely real response at
Nyquist ($w\rightarrow\infty$), which similarly gives a physical IIR
filter.

Although higher orders $n$ would appear to produce better (i.e.\
sharper) filter responses, one rapidly runs into numerical errors, as
one ends up adding and subtracting $n$ large numbers to obtain small
filter responses.  One way around this is to break the filter up into
several lower-order filters.  The routines in this module do just
that.  Poles are paired up across the imaginary axis, (and combined
with pairs of zeros at $w=0$ for high-pass filters,) to form $[n/2]$
second-order filters.  If $n$ is odd, there will be an additional
first-order filter, with one pole at $w=iw_c$ (and one zero at $w=0$
for a high-pass filter).

Each ZPG filter in the $w$-plane is first transformed to the $z$-plane
by a bilinear transformation, and is then used to construct a
time-domain IIR filter.  Each filter is then applied to the time
series.  As mentioned in the description above, the filters are
designed to give an overall amplitude response that is the square root
of the desired attenuation; however, each time-domain filter is
applied to the data stream twice: once in the normal sense, and once
in the time-reversed sense.  This gives the full attenuation with very
little frequency-dependent phase shift.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALWarning()
LALCreateREAL4IIRFilter()       LALCreateREAL8IIRFilter()
LALCreateCOMPLEX8ZPGFilter()    LALCreateCOMPLEX16ZPGFilter()
LALDestroyREAL4IIRFilter()      LALDestroyREAL8IIRFilter()
LALDestroyCOMPLEX8ZPGFilter()   LALDestroyCOMPLEX16ZPGFilter()
LALWToZCOMPLEX8ZPGFilter()      LALWToZCOMPLEX16ZPGFilter()
LALIIRFilterREAL4Vector()       LALIIRFilterREAL8Vector()
LALIIRFilterREAL4VectorR()      LALIIRFilterREAL8VectorR()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ButterworthTimeSeriesCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>

NRCSID(BUTTERWORTHTIMESERIESC,"$Id$");

extern INT4 lalDebugLevel;

static INT4
ParsePassBandParamStruc( LALStatus          *stat,
			 PassBandParamStruc *params,
			 INT4               *n,
			 REAL8              *wc,
			 REAL8              deltaT );
/* Prototype for a local input parsing routine. */


/* <lalVerbatim file="ButterworthTimeSeriesD"> */
void
LALButterworthREAL4TimeSeries( LALStatus          *stat,
			       REAL4TimeSeries    *series,
			       PassBandParamStruc *params )
{ /* </lalVerbatim> */
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  INITSTATUS(stat,"LALButterworthREAL4TimeSeries",BUTTERWORTHTIMESERIESC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure the input pointers are non-null. */
  ASSERT(params,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=ParsePassBandParamStruc(stat,params,&n,&wc,series->deltaT);
  if(type==0) {
    ABORT(stat,BANDPASSTIMESERIESH_EBAD,BANDPASSTIMESERIESH_MSGEBAD);
  }

  /* An order n Butterworth filter has n poles spaced evenly along a
     semicircle in the upper complex w-plane.  By pairing up poles
     symmetric across the imaginary axis, the filter gan be decomposed
     into [n/2] filters of order 2, plus perhaps an additional order 1
     filter.  The following loop pairs up poles and applies the
     filters with order 2. */
  for(i=0,j=n-1;i<j;i++,j--){
    REAL4 theta=LAL_PI*(i+0.5)/n;
    REAL4 ar=wc*cos(theta);
    REAL4 ai=wc*sin(theta);
    REAL4IIRFilter *iirFilter=NULL;
    COMPLEX8ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,2,2),
	  stat);
      zpgFilter->zeros->data[0].re=0.0;
      zpgFilter->zeros->data[0].im=0.0;
      zpgFilter->zeros->data[1].re=0.0;
      zpgFilter->zeros->data[1].im=0.0;
      zpgFilter->gain.re=1.0;
      zpgFilter->gain.im=0.0;
    }else{
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,0,2),
	  stat);
      zpgFilter->gain.re=-wc*wc;
      zpgFilter->gain.im=0.0;
    }
    zpgFilter->poles->data[0].re=ar;
    zpgFilter->poles->data[0].im=ai;
    zpgFilter->poles->data[1].re=-ar;
    zpgFilter->poles->data[1].im=ai;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX8ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL4IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALIIRFilterREAL4Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALIIRFilterREAL4VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    REAL4IIRFilter *iirFilter=NULL;
    COMPLEX8ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,1,1),
	  stat);
      zpgFilter->zeros->data->re=0.0;
      zpgFilter->zeros->data->im=0.0;
      zpgFilter->gain.re=1.0;
      zpgFilter->gain.im=0.0;
    }else{
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,0,1),
	  stat);
      zpgFilter->gain.re=0.0;
      zpgFilter->gain.im=-wc;
    }
    zpgFilter->poles->data->re=0.0;
    zpgFilter->poles->data->im=wc;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX8ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL4IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALIIRFilterREAL4Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALIIRFilterREAL4VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Normal exit. */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="ButterworthTimeSeriesD"> */
void
LALButterworthREAL8TimeSeries( LALStatus          *stat,
			       REAL8TimeSeries    *series,
			       PassBandParamStruc *params )
{ /* </lalVerbatim> */
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  INITSTATUS(stat,"LALButterworthREAL8TimeSeries",BUTTERWORTHTIMESERIESC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure the input pointers are non-null. */
  ASSERT(params,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=ParsePassBandParamStruc(stat,params,&n,&wc,series->deltaT);
  if(type==0) {
    ABORT(stat,BANDPASSTIMESERIESH_EBAD,BANDPASSTIMESERIESH_MSGEBAD);
  }

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
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,2,2),
	  stat);
      zpgFilter->zeros->data[0].re=0.0;
      zpgFilter->zeros->data[0].im=0.0;
      zpgFilter->zeros->data[1].re=0.0;
      zpgFilter->zeros->data[1].im=0.0;
      zpgFilter->gain.re=1.0;
      zpgFilter->gain.im=0.0;
    }else{
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,0,2),
	  stat);
      zpgFilter->gain.re=-wc*wc;
      zpgFilter->gain.im=0.0;
    }
    zpgFilter->poles->data[0].re=ar;
    zpgFilter->poles->data[0].im=ai;
    zpgFilter->poles->data[1].re=-ar;
    zpgFilter->poles->data[1].im=ai;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX16ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL8IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALIIRFilterREAL8Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALIIRFilterREAL8VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,1,1),
	  stat);
      zpgFilter->zeros->data->re=0.0;
      zpgFilter->zeros->data->im=0.0;
      zpgFilter->gain.re=1.0;
      zpgFilter->gain.im=0.0;
    }else{
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,0,1),
	  stat);
      zpgFilter->gain.re=0.0;
      zpgFilter->gain.im=-wc;
    }
    zpgFilter->poles->data->re=0.0;
    zpgFilter->poles->data->im=wc;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX16ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL8IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALIIRFilterREAL8Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALIIRFilterREAL8VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Normal exit. */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="ButterworthTimeSeriesD"> */
void
LALDButterworthREAL4TimeSeries( LALStatus          *stat,
				REAL4TimeSeries    *series,
				PassBandParamStruc *params )
{ /* </lalVerbatim> */
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  INITSTATUS(stat,"LALButterworthREAL8TimeSeries",BUTTERWORTHTIMESERIESC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure the input pointers are non-null. */
  ASSERT(params,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=ParsePassBandParamStruc(stat,params,&n,&wc,series->deltaT);
  if(type==0) {
    ABORT(stat,BANDPASSTIMESERIESH_EBAD,BANDPASSTIMESERIESH_MSGEBAD);
  }

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
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,2,2),
	  stat);
      zpgFilter->zeros->data[0].re=0.0;
      zpgFilter->zeros->data[0].im=0.0;
      zpgFilter->zeros->data[1].re=0.0;
      zpgFilter->zeros->data[1].im=0.0;
      zpgFilter->gain.re=1.0;
      zpgFilter->gain.im=0.0;
    }else{
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,0,2),
	  stat);
      zpgFilter->gain.re=-wc*wc;
      zpgFilter->gain.im=0.0;
    }
    zpgFilter->poles->data[0].re=ar;
    zpgFilter->poles->data[0].im=ai;
    zpgFilter->poles->data[1].re=-ar;
    zpgFilter->poles->data[1].im=ai;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX16ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL8IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALDIIRFilterREAL4Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALDIIRFilterREAL4VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,1,1),
	  stat);
      zpgFilter->zeros->data->re=0.0;
      zpgFilter->zeros->data->im=0.0;
      zpgFilter->gain.re=1.0;
      zpgFilter->gain.im=0.0;
    }else{
      TRY(LALCreateCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter,0,1),
	  stat);
      zpgFilter->gain.re=0.0;
      zpgFilter->gain.im=-wc;
    }
    zpgFilter->poles->data->re=0.0;
    zpgFilter->poles->data->im=wc;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX16ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL8IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALDIIRFilterREAL4Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALDIIRFilterREAL4VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL8IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX16ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Normal exit. */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


static INT4
ParsePassBandParamStruc( LALStatus          *stat,
			 PassBandParamStruc *params,
			 INT4               *n,
			 REAL8              *wc,
			 REAL8              deltaT )
     /* This local function parses the pass band parameters according
        to the rules given in the documentation to this module,
        computing the order and characteristic frequency (in the
        w-variable) of the desired filter.  The function returns 2 for
        a high-pass filter, 1 for a low-pass filter, and 0 if the
        parameters were poorly specified. */
{
  /* Okay, first define all the temporary variables, and figure out
     just which parameter combinations have been given. */
  REAL8 w1=deltaT*params->f1;    /* Dimensionless frequencies; */
  REAL8 w2=deltaT*params->f2;    /* later transformed to w-plane. */
  REAL8 wLow = w1<w2 ? w1 : w2;  /* Sorted frequencies. */
  REAL8 wHigh = w1>w2 ? w1 : w2;
  REAL8 a1=params->a1;           /* Attenuation factors; later a */
  REAL8 a2=params->a2;           /* square root will be taken. */
  REAL8 aLow = w1<w2 ? a1 : a2;  /* Sorted according to w's. */
  REAL8 aHigh = w1>w2 ? a1 : a2;
  BOOLEAN wLowGiven = (wLow>0.0)&&(wLow<0.5);
  BOOLEAN bothLowGiven = wLowGiven&&(aLow>0.0)&&(aLow<1.0);
  BOOLEAN wHighGiven = (wHigh>0.0)&&(wHigh<0.5);
  BOOLEAN bothHighGiven = wHighGiven&&(aHigh>0.0)&&(aHigh<1.0);

  /* If both frequencies and both attenuations have been given,
     compute the required filter order and characteristic frequency
     from them. */
  if(bothLowGiven && bothHighGiven){
    aLow=sqrt(aLow);
    aHigh=sqrt(aHigh);
    wLow=tan(LAL_PI*wLow);
    wHigh=tan(LAL_PI*wHigh);

    /* First make sure that two different frequencies and attenuations
       have been specified. */
    if(w1==w2){
      LALInfo(stat,"The two frequencies defining the transition band"
	      " are the same.");
      return 0;
    }
    if(a1==a2){
      LALInfo(stat,"The two attenuations across the transition band"
	      " are the same.");
      return 0;
    }

    /* Now compute the filter order. */
    if(aHigh>aLow) /* High-pass. */
      *n=(INT4)(0.5*log((1.0/aHigh-1.0)/(1.0/aLow-1.0))
		/log(wLow/wHigh))+1;
    else /* Low-pass. */
      *n=(INT4)(0.5*log((1.0/aHigh-1.0)/(1.0/aLow-1.0))
		/log(wHigh/wLow))+1;

    /* If a positive params->nMax less than *n has been specified,
       reduce to that order, with appropriate warnings. */
    if((params->nMax>0)&&(params->nMax<*n)){
      LALWarning(stat,"Filter order required to achieve requested"
		 " performance exceeds\n"
		 "\tspecified limit.");
      if(lalDebugLevel&LALWARNING)
	LALPrintError("\tRequired: %i  Limit: %i\n",*n,params->nMax);
      *n=params->nMax;
    }

    /* Compute the characteristic frequency of the filter using the
       filter order and the minimal attenuation in the pass-band. */
    if(aHigh>aLow){ /* High-pass */
      *wc=wHigh*pow(1.0/aHigh - 1.0,-0.5/((REAL8)(*n)));
      return 2;
    }else{ /* Low-pass */
      *wc=wLow*pow(1.0/aLow - 1.0,0.5/((REAL8)(*n)));
      return 1;
    }
  }

  /* If one end of the transition band has both frequency and
     attenuation specified, and the other does not, use only the end
     that does.  The filter order must also be specified for this and
     all future cases. */
  else{
    if(params->nMax<=0){
      LALInfo(stat,"The filter order must be given if one or more"
	      " frequencies or\n"
	      "\tattenuations of the transition band have not been"
	      " (validly) specified.");
      return 0;
    }
    if(bothLowGiven){
      wLow=tan(LAL_PI*wLow);
      aLow=sqrt(aLow);
      *n=params->nMax;
      *wc=wLow*pow(1.0/aLow - 1.0, -0.5/((REAL8)(*n)));
      return 1;
    }else if(bothHighGiven){
      wHigh=tan(LAL_PI*wHigh);
      aHigh=sqrt(aHigh);
      *n=params->nMax;
      *wc=wHigh*pow(1.0/aHigh - 1.0, 0.5/((REAL8)(*n)));
      return 2;
    }

    /* If no attenuations are given, then there had better be only one
       frequency given, otherwise we don't know whether to make a low-
       or a high-pass filter. */
    else if(wHighGiven && wLowGiven){
      LALInfo(stat,"Neither attenuation has been specified, so only"
	      " one frequency should be.");
      return 0;
    }

    /* Treat the given frequency as the characteristic frequency and
       the specified nMax as the filter order. */
    else{
      *n=params->nMax;
      if(wHighGiven){
	*wc=tan(LAL_PI*wHigh);
	return 2;
      }else if(wLowGiven){
	*wc=tan(LAL_PI*wLow);
	return 1;
      }else{
	LALInfo(stat,"No frequencies within the Nyquist band have"
		" been specified!");
	return 0;
      }
    }
  }
}
