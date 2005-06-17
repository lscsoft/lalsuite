/********************************** <lalVerbatim file="SlopeDetectorFilterCV">
Author: Daw, E. J.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{SlopeDetectorFilter.c}}

Functions in this code segment implement various time domain 
search algorithms for transients and bursts in the LIGO data.

\subsubsection*{Prototypes}
\input{SlopeDetectorFilterCP}
\idx{LALSlopeDetectorFilter()}
\input{SlopeLineFitFilterCP}
\idx{LALSlopeLineFitFilter()}
\input{SlopeConvolutionFilterCP}
\idx{LALSlopeConvolutionFilter()}

\subsubsection*{Description}

The three fine fuctions prototyped above implement time domain linear and nonlinear
filters on input data. The resulting output data is intended to enhance
sensitivity to burst sources. 

The first function is \texttt{void LALSlopeDetectorFilter()}, which
takes input data through its third argument
\texttt{const REAL4Vector* input\_data}. For
each element $i$ of the input data ntuple $x_i$, the best fit slope
the data points $i$ to $i+N-1$ is found, where $N$ is the fourth
argument \texttt{ntaps} passed to the function. This function 
provides no provision for history buffers. The output data is written
to the address \texttt{REAL4Vector* output\_data} passed as the second
argument.

The second function is \texttt{LALSlopeLineFitFilter()}, which also
takes a \texttt{const REAL4Vector* input\_data} pointer to input data as
its third argument. For each element $i$ of the
input data ntuple $x_i$, the best fit slope, $a_i$ , or the best
fit intercept to the data (vertical) axis, $b_i$, or a particular
nonlinear combination of $a_i$ and $b_i$ is found. Any of these
three quantities can be returned to the output address 
\texttt{REAL4Vector* output\_data}, depending on the argument passed
to \texttt{fparams.function\_select}. Acceptable values for this parameter
are shown in table \ref{slopetable:fitsettings}.

\begin{table}
\begin{center}
\begin{tabular}{|l|c|r|} \hline
\texttt{fparams.function\_select} & 
value & filter output \\ \hline\hline
\texttt{FILTER\_OUTPUT\_SLOPE} & 1 & fit to slope $a_i$ over N bins \\ \hline
\texttt{FILTER\_OUTPUT\_OFFSET} & 2 & fit to offset $b_i$ over N bins \\ \hline
\texttt{FILTER\_OUTPUT\_ALF} & 3 & ALF filter (see algorithms) \\ \hline
\end{tabular}
\caption{settings for \texttt{fparams.function\_select} used in conjunction
with function \texttt{LALSlopeLineFitFilter()}.}
\label{slopetable:fitsettings}
\end{center}
\end{table}

The third function is \texttt{LALSlopeConvolutionFilter()}.This
function convolves input data with time domain
finite impulse response (FIR) filters of four types. Before
running the filter on data, the filter taps must be set. Allocate
enough memory to hold the number of filter taps you want. Set
\texttt{fparams.forder} to be equal to the number of taps.
Set the \texttt{fparams.function\_select} to one of the values given
in table \ref{slopetable:convsettings}.

\begin{table}
\begin{center}
\begin{tabular}{|l|c|r|} \hline
\texttt{fparams.function\_select} & 
value & filter output \\ \hline\hline
\texttt{FILTER\_OUTPUT\_BOXCAR} & 4 & 
set taps to uniform window over N bins \\ \hline
\texttt{FILTER\_OUTPUT\_GAUSSIAN} 
& 5 & set taps to Gaussian over N bins \\ \hline
\texttt{FILTER\_OUTPUT\_SINE} & 6 & 
set taps to sine wave period over N bins \\ \hline
\texttt{FILTER\_OUTPUT\_USER} & 7 & 
set taps to user defined \\ \hline
\texttt{FILTER\_OUTPUT\_CONVOLVE} & 8 & 
convolve data with a filter already set up \\ \hline
\end{tabular}
\caption{settings for \texttt{fparams.function\_select} used in conjunction
with function \texttt{LALSlopeConvolutionFilter()}.}
\label{slopetable:convsettings}
\end{center}
\end{table}

The \texttt{BOXCAR} option creates N taps that are all set to 1/N.
The \texttt{GAUSSIAN} option creates taps that represent part of
a gaussian curve between $\pm 3\sigma$. The \texttt{SINE} option
sets the taps to be one period of a sine wave. 
Both the sine and Gaussian tap waveforms are of unity height. If the 
\texttt{USER} option is set, then the \texttt{input\_data} field
is used to set the taps. In the latter case, the 
\texttt{input\_data->length} field must be the same as
\texttt{fparams.forder} field, and equal the number of elements
in the input data array. In all cases, memory must be allocated 
by the calling function to the \texttt{fparams.taps\_set} field,
which must be initialized to zero. Once the function is called
to set the taps, the \texttt{fparams.taps\_set} field is changed to
one. 

Once the filter has been initialized, it can be convolved with
data. History buffering is supported, and works in the same way
as for the \texttt{LALSlopeLineFitFilter} function described
above. Note that the \texttt{fparams.function\_select} field
must be set to \texttt{FILTER\_OUTPUT\_CONVOLVE}. 

\subsubsection*{Algorithm}

The algorithms used for applying filters with the three functions
described above are given below. For all formulae, $x_i$ represents
an array of input data. The filters embedded in the functions
\texttt{LALSlopeDetectionFilter} and \texttt{LALSlopeLineFitFilter}
work by fitting each overlapping set of N successive data points
to a straight line $a_i x_{i+j} + b_i$. The fit parameters $a_i$
and $b_i$ are used as the basis to form discriminants. The
\texttt{LALSlopeConvolutionFilter} function convolves input
data either with a standard set of waveforms, or with waveforms
set by the user as described above. In all the equations below,
N is the number of successive data samples used to compute the
quantities, the filter `order'. 

For \texttt{LALSlopeDetectionFilter} the output data array is $y_i$:

\begin{equation}
y_i = \frac{
\frac{1}{N} \sum_{j=0}^{N-1} x_{i+j} t_{i+j} - 
\left( \frac{1}{N} \sum_{k=0}^{N-1} x_{i+k} \right)
\left( \frac{1}{N} \sum_{m=0}^{N-1} t_{i+m} \right)
}{
\frac{1}{N} \sum_{p=0}^{N-1} t_{i+p}^2 -
\left( \frac{1}{N} \sum_{q=0}^{N-1} t_{i+q} \right)^2
}
\label{eqn:slopedetectone}
\end{equation}

For
\texttt{LALSlopeLineFitFilter} with 
\texttt{fparams.function\_select} set to \texttt{FILTER\_OUTPUT\_SLOPE},
the output data $a_i$ is

\begin{equation}
a_i = \frac{12}{\tau N(N^2 - 1)}
\left( \sum_{j=0}^{N-1} jx_j - \frac{N-1}{2} \sum_{k=0}^{N-1} x_k \right),
\label{eqn:slopedgrad}
\end{equation}

where $\tau$ is the sampling period. Apart from a factor of $\tau$, equations
\ref{eqn:slopedetectone} and \ref{eqn:slopedgrad} are equivalent,
since we have the following identities:

\begin{eqnarray}
\sum_{j=0}^{N-1} j & = & \frac{N(N-1)}{2} \\ \nonumber
\sum_{j=0}^{N-1} j^2 & = & \frac{N}{6} (2N-1)(N-1) \\ \nonumber
\label{eqn:slopepartialsums}
\end{eqnarray}

Note that slope filters have a mimimum order of 2, since it takes at least 2 points
to estimate the slope.

For 
\texttt{LALSlopeLineFitFilter} with
\texttt{fparams.function\_select} set to \texttt{FILTER\_OUTPUT\_OFFSET},
$a_i$ is calculated and used to estimate the output $b_i$ given by

\begin{equation}
b_i = \sum_{j=0}^{N-1} x_j - \frac{a_i \tau (N-1)}{2}.
\label{eqn:slopedetecttwogradient}
\end{equation}

For
\texttt{LALSlopeLineFitFilter} with
\texttt{fparams.function\_select} set to \texttt{FILTER\_OUTPUT\_ALF},
a nonlinear combination of $a_i$ and $b_i$ is constructed.
Define $\sigma_a$ and $\sigma_b$:

\begin{eqnarray}
\sigma_a & = & \frac{12}{\tau N(N^2 - 1)} \\ \nonumber
\sigma_b & = & \frac{4N+2}{N(N-1)}. \\ \nonumber
\label{eqn:slopedetecttwoalfone}
\end{eqnarray}

Let $X_a = a_i / \sigma_a$ and $X_b = b_i / \sigma_b$. Then the
output $c_i$ of the \texttt{ALF} filter is 

\begin{equation}
c_i = \frac{
X_a^2 + X_b^2 - 2 \alpha X_a X_b}{
1 - \alpha^2},
\label{eqn:slopedetecttwoalftwo}
\end{equation}

where

\begin{equation}
\alpha = -\sqrt{
\frac{3(N+1)}{2(2N+1)}.
\label{eqn:slopedetecttwoalfthree}
}
\end{equation}

This nonlinear combination was found by the Orsay group to 
be particularly sensitive to bursts in the Zwerger Mueller catalog.
See \texttt{gr/qc-0010037} for a full description.

The remaining filters use a standard convolution algorithm to
convolve input data with the waveforms described above (boxcar,
gaussian, sinewave, user). The convolution of input data
$x_i$ with filter taps $c_j$ is

\begin{equation}
y_i = \sum_{j=0}^{N-1} x_{i+j} c_j
\label{eqn:slopedetectconv}
\end{equation}

\subsubsection*{Uses}

% List any external functions called by this function.
\begin{verbatim}
These functions do not use any other LAL routines.
\end{verbatim}

\subsubsection*{Notes}

All pointers passed to these functions must be preallocated by
the calling routine to have enough memory to hold the objects
written to them. In particular, the pointer to the output data
must be of size \texttt{input\_data\_length - filter\_order + 1},
where \texttt{filter\_order} is \texttt{ntaps} in the function
\texttt{LALSlopeDetectorFilter} and \texttt{fparams.forder} in
\texttt{LALSlopeLineFitFilter} and \texttt{LALSlopeConvolutionFilter}.
Unused parameters that are pointers can be set to \texttt{NULL} (for
example, the \texttt{fparams.tap} field can be null if calling
the \texttt{LALSlopeLineFitFilter} function.

The functions
\texttt{LALSlopeLineFitFilter} and \texttt{LALSlopeConvolutionFilter}
support the passing of a history buffer between successive applications
of a filter. The history buffer pointer needs to be allocated by the
calling routine when using these functions. Enough memory must be
allocated to it to hold \texttt{fparams.forder - 1} \texttt{REAL4}
elements. After running the filter, the history buffer will contain
the last \texttt{fparams.forder - 1} elements of the input data sent
to the filter. By adding these elements to the beginning of the next
set of input data sent to the filter, a continuous filtered data stream
with a seamless boundary may be maintained between successive 
applications of the filter. Note that currently the passing of history
buffers between successive MPI jobs is not supported under LDAS. 

Both \texttt{LALSlopeLineFitFilter} and \texttt{LALSlopeConvolutionFilter}
require the field \texttt{fparams.sampling\_period\_s} to contain
the reciprocal of the sampling frequency for the channel. 

The \texttt{waveform\_offset} field in \texttt{SLOPEFilterParam} is
there in case you need to offset the filter waveform by a fraction
of a frequency bin (eg, if you want the maximum of a gaussian filter
to fall equally between two bins, you would set
this field to $0.5$).

\vfill{\footnotesize\input{SlopeDetectorFilterCV}}

*** </lalLaTeX> */ 

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

#include <math.h>

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/SlopeDetectorFilter.h>

/******* DEFINE RCS ID STRING ************/
NRCSID( SLOPEDETECTORFILTERC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/********************************************/
/**** LALSlopeDetectorFilter ****************/
/********************************************/

static void
CreateGaussianTaps( UINT4 ntaps,
                    REAL4 binoffset,
                    REAL4* tap );

static void 
CreatePeriodSineTaps( UINT4 ntaps,
                      REAL4 binoffset,
                      REAL4* tap );



/* <lalVerbatim file="SlopeDetectorFilterCP"> */
void
LALSlopeDetectorFilter( LALStatus          *status,
			REAL4Vector*       output_data,
			const REAL4Vector* input_data,
			const UINT4        ntaps )
/* </lalVerbatim> */
{
  /******* DECLARE VARIABLES; for example: ************/

  REAL4              meanxt,meant,meanx,meantsquared;
  UINT4              i,j,datalength;
  REAL4*             datain;
  REAL4*             dataout;
  REAL4              tindex, nreal;

  INITSTATUS( status, "LALSlopeDetectorFilter", SLOPEDETECTORFILTERC );
  /* only needed if you call other LAL functions within your code */
  /* ATTATCHSTATUSPTR (status); */

  if ( input_data == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
  }
  if ( output_data == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
  }

  datalength = input_data->length;
  datain = input_data->data;
  dataout = output_data->data;

  /******* CHECK VALIDITY OF ARGUMENTS; for example ************/

  if ( dataout == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
  }
  if ( datain == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
  }
  if ( datalength < ntaps ) {
    ABORT( status, SLOPEDETECTORFILTERH_EDATATOOSHORT, 
	   SLOPEDETECTORFILTERH_MSGEDATATOOSHORT );
  }

  /* printf("CODE: At div by zero point, ntaps is %u\n",ntaps); */
  /* check that the number of taps is not zero */
  if ( ntaps == 0 ) {
    ABORT( status, SLOPEDETECTORFILTERH_EDIVBYZERO, 
	   SLOPEDETECTORFILTERH_MSGEDIVBYZERO );
  }

  /******* DO ANALYSIS ************/

  nreal=(REAL4)ntaps;
  for(i=0;i<(datalength-ntaps+1);++i) {
    meanx=0;
    meant=0; 
    meanxt=0; 
    meantsquared=0;
    for(j=0;j<ntaps;++j) {
      tindex = (REAL4)j;
      meanx += datain[i+j]/nreal;
      meant += tindex/nreal;
      meanxt += tindex*datain[i+j]/nreal;
      meantsquared += tindex*tindex/nreal;
    }

    /* check that the denominator is not zero */
    if ( (meantsquared - meant*meant) == 0 ) {
      ABORT( status, SLOPEDETECTORFILTERH_EDIVBYZERO, 
	     SLOPEDETECTORFILTERH_MSGEDIVBYZERO );
    }
    
    dataout[i] = (meanxt - meanx * meant) / ( meantsquared - meant * meant );

    /*these lines present for testing purposes only*/
    /*printf("%u\t%f\t%f\t%f\t%f\t%f\n",
      i,meanx,meant,meanxt,meantsquared,dataout[i]); */
    /*dataout[i] = datain[i];*/
  }
  
  /* DETATCHSTATUSPTR(status); */
  RETURN(status);
}

/**********************************************/
/*** LALSlopeLineFitFilter ********************/
/**********************************************/

/* <lalVerbatim file="SlopeLineFitFilterCP"> */
void
LALSlopeLineFitFilter( LALStatus                *status,
		       REAL4Vector*             output_data,
		       const REAL4Vector*       input_data,
		       const SLOPEFilterParams  fparams )
/* </lalVerbatim> */
{
  /******* DECLARE VARIABLES ***************************/

  REAL4              sumxt,sumx,meant,normalization;
  UINT4              i,j,datalength;
  REAL4*             datain;
  REAL4*             dataout;
  REAL4              calculated_slope,calculated_offset,sigmaa,sigmab,alpha; 
  REAL4              xnorma, xnormb;
  REAL4              nreal;

  /******* FUNCTION BODY ********************************/

  INITSTATUS( status, "LALSlopeLineFitFilter", SLOPEDETECTORFILTERC );
  /* only needed if you call other LAL functions within your code */
  /* ATTATCHSTATUSPTR (status); */

  if ( input_data == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
  }
  if ( output_data == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
  }
  if ( fparams.history == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EHISTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEHISTNULLP );
  }
  /* filter order needs to be at least two for slope detector filters */
  /* since you can't fit the slope off one data point */
  if ( (fparams.forder < 2) || 
       (fparams.forder > MAX_SLOPE_DETECTOR_FILTER_ORDER ) ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINVFILTLEN, 
	   SLOPEDETECTORFILTERH_MSGEINVFILTLEN );
  }
  if ( ( *(fparams.history_allocated) != 0) &&
       ( *(fparams.history_allocated) != 1) ) {
    ABORT( status, SLOPEDETECTORFILTERH_EHISTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEHISTNULLP );
  }

  datalength = input_data->length;
  datain = input_data->data;
  dataout = output_data->data;
  nreal = (REAL4)fparams.forder;

  if ( dataout == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
  }
  if ( datain == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
  }
  if ( datalength <= fparams.forder ) {
    ABORT( status, SLOPEDETECTORFILTERH_EDATATOOSHORT, 
	   SLOPEDETECTORFILTERH_MSGEDATATOOSHORT );
  }

  /* add checks that the input series is valid data, etc */

  /******* DO ANALYSIS ************/

  /* run filter */  

  /* calculate overall normalization to get output in counts per second */
  normalization = 12 / (nreal * 
			fparams.sampling_period_s *
			(nreal * nreal - 1));

  meant = fparams.sampling_period_s * ( nreal - 1 ) / 2;
  
  /* loop over input data */
  for(i=0;i<(datalength-fparams.forder+1);++i) {

    /* compute sums needed for calculating best fit slope */ 
    sumx=0; sumxt=0;
    for(j=0;j<fparams.forder;++j) {
      sumx += datain[i+j];
      sumxt += j*datain[i+j];
    }

    /* calculate slope, needed for all filter outputs */
    calculated_slope
      = normalization * (sumxt - meant * sumx / fparams.sampling_period_s );

    switch( fparams.function_select ) {
    case FILTER_OUTPUT_SLOPE:
      /* the output is the slope */
      dataout[i] = calculated_slope;
      break;
    case FILTER_OUTPUT_OFFSET:
      /* the output is the offset */
      dataout[i] = sumx / nreal - calculated_slope * meant;
      break;
    case FILTER_OUTPUT_ALF:
      /* output data is a nonlinear function of slope and offset */
      sigmaa = (1 / fparams.sampling_period_s) * (REAL4)sqrt 
	( 12 / (    (REAL8)nreal * ( (REAL8)nreal*(REAL8)nreal - 1) )    ) ;
      sigmab = (REAL4)sqrt ( (4 * (REAL8)nreal + 2) / 
		      ( (REAL8)nreal * ( (REAL8)nreal - 1 ) ) );
      calculated_offset = sumx / nreal - calculated_slope * meant;
      alpha = ( 0 - sqrt( 1.5 * ( (REAL8)nreal + 1) / ( 2*(REAL8)nreal + 1 ) ) );
      xnorma = calculated_slope / sigmaa;
      xnormb = calculated_offset / sigmab;
      dataout[i] = ( ( xnorma*xnorma + xnormb*xnormb - 2*alpha*xnorma*xnormb ) /
		     ( 1 - alpha * alpha ) );
      break;
    default:
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDACTION, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDACTION );

    }
  }
  
  /* set history buffer */
  for(i=0;i<(fparams.forder-1);++i) {
    fparams.history[i] = datain[datalength-(fparams.forder-1)+i];
  }
  *(fparams.history_allocated) = 1;

  /* DETATCHSTATUSPTR(status); */
  RETURN(status);

}

/*************************************************/
/*** LALSlopeConvolutionFilter *******************/
/*************************************************/

/* <lalVerbatim file="SlopeConvolutionFilterCP"> */
void
LALSlopeConvolutionFilter( LALStatus                *status,
			   REAL4Vector*             output_data,
			   const REAL4Vector*       input_data,
			   const SLOPEFilterParams  fparams )
/* </lalVerbatim> */
{
  /******* DECLARE VARIABLES; for example: ************/

  REAL4* datain;
  REAL4* dataout;
  REAL4 nreal;
  UINT4 datalength,i,j;

  INITSTATUS( status, "LALSlopeConvolutionFilter", SLOPEDETECTORFILTERC );
  /* only needed if you call other LAL functions within your code */
  /* ATTATCHSTATUSPTR (status);                                   */

  switch( fparams.function_select ) {
  case FILTER_OUTPUT_TAPS_SET_GAUSSIAN:
    if ( ( fparams.forder < 1) || 
	 ( fparams.forder > MAX_SLOPE_DETECTOR_FILTER_ORDER ) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVFILTLEN, 
	     SLOPEDETECTORFILTERH_MSGEINVFILTLEN );
    }
    if ( *(fparams.taps_set) != 0 ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDTAPSBIT, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDTAPSBIT );
    }
    if ( fparams.tap == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_ETAPSNULLP, 
	     SLOPEDETECTORFILTERH_MSGETAPSNULLP );
    }
    if ( (fparams.waveform_offset < 0) || 
	 (fparams.waveform_offset > 1) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EBINOFFINVALID, 
	     SLOPEDETECTORFILTERH_MSGEBINOFFINVALID );
    }
    /* create the filter taps */
    CreateGaussianTaps(fparams.forder,fparams.waveform_offset,fparams.tap);
    /* set the taps bit */
    *(fparams.taps_set) = 1;
    break;
  case FILTER_OUTPUT_TAPS_SET_SINE:
    if ( ( fparams.forder < 1 ) || 
	 ( fparams.forder > MAX_SLOPE_DETECTOR_FILTER_ORDER ) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVFILTLEN, 
	     SLOPEDETECTORFILTERH_MSGEINVFILTLEN );
    }
    if ( *(fparams.taps_set) != 0 ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDTAPSBIT, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDTAPSBIT );
    }
    if ( fparams.tap == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_ETAPSNULLP, 
	     SLOPEDETECTORFILTERH_MSGETAPSNULLP );
    }
    if ( (fparams.waveform_offset < 0) || 
	 (fparams.waveform_offset > 1) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EBINOFFINVALID, 
	     SLOPEDETECTORFILTERH_MSGEBINOFFINVALID );
    }
    /* create the filter taps */
    CreatePeriodSineTaps(fparams.forder,fparams.waveform_offset,fparams.tap);
    /* set the taps bit */
    *(fparams.taps_set) = 1;
    break;
  case FILTER_OUTPUT_TAPS_SET_BOXCAR:
    if ( ( fparams.forder < 1 ) || 
	 ( fparams.forder > MAX_SLOPE_DETECTOR_FILTER_ORDER ) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVFILTLEN, 
	     SLOPEDETECTORFILTERH_MSGEINVFILTLEN );
    }
    if ( *(fparams.taps_set) != 0 ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDTAPSBIT, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDTAPSBIT );
    }
    if ( fparams.tap == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_ETAPSNULLP, 
	     SLOPEDETECTORFILTERH_MSGETAPSNULLP );
    }
    if ( (fparams.waveform_offset < 0) || 
	 (fparams.waveform_offset > 1) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EBINOFFINVALID, 
	     SLOPEDETECTORFILTERH_MSGEBINOFFINVALID );
    }
    nreal = (REAL4)fparams.forder;
    /*printf("CODE: nreal = %f\n",nreal);*/
    /* create the filter taps */
    for(i=0;i<fparams.forder;++i) {
      *(fparams.tap + i) = 1/nreal;
      /*printf("CODE: tap %d is %f.\n",i,*(fparams.tap + i));*/
    }
    /* set the taps bit */
    *(fparams.taps_set) = 1;
    break;
  case FILTER_OUTPUT_TAPS_SET_USER:
    if ( input_data == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	     SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
    }
    datalength = input_data->length;
    datain = input_data->data;
    /* check that the data array is of a usable length */
    if ( ((INT4)datalength < 0) ||
	 (datalength > MAX_SLOPE_DETECTOR_FILTER_ORDER) ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVFILTLEN, 
	     SLOPEDETECTORFILTERH_MSGEINVFILTLEN );
    }
    /* check that the taps bit is set to zero */
    if ( *(fparams.taps_set) != 0 ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDTAPSBIT, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDTAPSBIT );
    }
    /* check that the taps pointer is not null */
    if ( fparams.tap == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_ETAPSNULLP, 
	     SLOPEDETECTORFILTERH_MSGETAPSNULLP );
    }
    /* copy the input data array to the filter taps */
    for(i=0;i<datalength;++i) {
      *(fparams.tap + i)=datain[i];
    }
    /* set taps bit */
    *(fparams.taps_set) = 1;
    break;
  case FILTER_OUTPUT_CONVOLVE:
    if ( input_data == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	     SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
    }
    if ( output_data == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	     SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
    }
    datalength = input_data->length;
    datain = input_data->data;
    dataout = output_data->data;
    /* check that the taps bit is set to one */
    if ( *(fparams.taps_set) != 1 ) {
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDTAPSBIT, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDTAPSBIT );
    }
    /* check that the taps pointer is not null */
    if ( fparams.tap == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_ETAPSNULLP, 
	     SLOPEDETECTORFILTERH_MSGETAPSNULLP );
    }
    /* check that the history buffer pointer is not null */
    if ( (fparams.history) == NULL ) {
      ABORT( status, SLOPEDETECTORFILTERH_EHISTNULLP, 
	     SLOPEDETECTORFILTERH_MSGEHISTNULLP );
    }
    /* convolve the data with the filter */
    for(i=0;i<(datalength-fparams.forder+1);++i) {

      /* calculate convolution for each input data point */
      dataout[i]=0;
      for(j=0;j<fparams.forder;++j) {
	dataout[i]+=fparams.tap[j]*datain[i+j];
      }

    }
    /* set the history buffer */
    for(i=0;i<(fparams.forder-1);++i) {
      fparams.history[i] = datain[datalength-(fparams.forder-1)+i];
    }
    *(fparams.history_allocated) = 1;
    break;
  default:
      ABORT( status, SLOPEDETECTORFILTERH_EINVALIDACTION, 
	     SLOPEDETECTORFILTERH_MSGEINVALIDACTION );
  }
  /* DETATCHSTATUSPTR(status); */
  RETURN(status);

}

/******** functions used interally by the slope detector library *******/

static void CreateGaussianTaps(UINT4 ntaps, REAL4 binoffset, REAL4* tap) {
  UINT4 i;
  REAL8 midbin,sigbins,offs;
  midbin = ((REAL8)ntaps - 1)/2;
  sigbins = midbin / (REAL8)GAUSSIAN_TAPS_NSIGMA_AT_EDGE;
  for(i=0;i<ntaps;++i) {
    offs = (REAL8)i;
    tap[i] = (REAL4)exp( - (offs - (midbin+(REAL8)binoffset)) * 
			 (offs - (midbin+(REAL8)binoffset)) /
			 (2 * sigbins * sigbins) );
  }
  return;
}

static void CreatePeriodSineTaps(UINT4 ntaps, REAL4 binoffset, REAL4* tap) {
  UINT4 i;
  REAL8 midbin/*,sigbins*/;
  REAL8 offs;
  midbin = ((REAL8)ntaps - 1)/2;
  for(i=0;i<ntaps;++i) {
    offs = (REAL8)i;
    tap[i] = (REAL4)sin( (offs - (midbin+(REAL8)binoffset))*SLOPE_PI/midbin );
  }
  return;
}
