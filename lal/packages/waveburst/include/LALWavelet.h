/********************************* <lalVerbatim file="LALWaveletHV">
Author: Klimenko, Sergey and Yakushin, Igor
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{LALWavelet.h}}

Provides structures for wavelets, clusters, pixels.
Declares prototypes of all the functions used by waveburst DSO.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALWavelet.h>
\end{verbatim}

This package provides all the structures and procedures needed to
work with wavelets, do percentile transform of wavelets, perform coincidence
between two channels, cluster the remaining pixels and select only those
events that satisfy certain criteria.

\subsection*{Error conditions}
\input{LALWaveletHE}

\subsection*{Structures}
\input{WaveburstStructs}





\vfill{\footnotesize\input{LALWaveletHV}}
\newpage\input{LALWaveletC}
%\newpage\input{clusterTestC.tex}
%\newpage\input{LALGetLayerWaveletTestC.tex}
%\newpage\input{setMaskTestC.tex}
********************************** </lalLaTeX> */

#ifndef _LALWAVELET_H
#define _LALWAVELET_H

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/Random.h>

/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALWAVELETH, "$Id$");


/******************************** <lalErrTable file="LALWaveletHE"> */

#define LALWAVELETH_ENULLP         1
#define LALWAVELETH_ETYPEMISMATCH  2
#define LALWAVELETH_EDIFF          3

#define LALWAVELETH_MSGENULLP        "Null pointer"
#define LALWAVELETH_MSGETYPEMISMATCH "Wavelet type mismatch"
#define LALWAVELETH_MSGEDIFF         "Difference between computed and expected values exceeds the threshold"

/************************************ </lalErrTable> */

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

#define TRUE 1
#define FALSE 0

/****** DEFINE NEW STRUCTURES AND TYPES ************/

enum BORDER {B_PAD_ZERO, 
	     B_CYCLE, 
	     B_MIRROR, 
	     B_PAD_EDGE, 
	     B_POLYNOM};


enum WAVETYPE {HAAR, 
	       BIORTHOGONAL,
	       DAUBECHIES, 
	       SYMLET,
	       DMEYER};

enum TREETYPE {DIADIC, BINARY};

  /*enum SLICETYPE {F,L};  */

  /*enum COINCIDENCETYPE {GG,GV};*/

typedef enum { ORIGINAL_CL, SWAPPED_CL, MIXED_CL } CLUSTER_TYPE;
typedef enum { NONE_CO=-1, CROSS_CO=0, RECTANGLE_CO=1 } COINCIDENCE_LEVEL;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{Slice}}
\index{\texttt{Slice}}
\noindent Slice structure is used to extract each \textbf{step}th element from an array starting from \textbf{start}.
\begin{description}
\item[\texttt{UINT4 start}] Index to start with.
\item[\texttt{UINT4 size}] Number of elements in the slice.
\item[\texttt{UINT4 step}] Skip $step-1$ elements when choosing the next one.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagSlice
{
  /*  enum SLICETYPE type;*/
  UINT4 start;
  UINT4 size;
  UINT4 step;
}
Slice;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{Wavelet}}
\index{\texttt{Wavelet}}
\noindent Wavelet structure contains all the wavelet data and metadata.
\begin{description}
\item[\texttt{enum WAVETYPE type}] Possible values: HAAR, BIORTHOGONAL, DAUBECHIES, SYMLET, DMEYER.
\item[\texttt{enum BORDER border}] Possible values: B\_PAD\_ZERO, B\_CYCLE,  B\_MIRROR,  B\_PAD\_EDGE, B\_POLYNOM.
\item[\texttt{enum TREETYPE treeType}] Possible values: DIADIC, BINARY.
\item[\texttt{UINT4 level}] Wavelet decomposition level.
\item[\texttt{HPFilterLength}] Highpass filter length.
\item[\texttt{LPFilterLength}] Lowpass filter length.
\item[\texttt{REAL4TimeSeries *data}] Wavelet coeficients.
%\item[\texttt{}] .

\end{description}
******************************************************* </lalLaTeX> */

typedef struct 
tagWavelet
{
  enum WAVETYPE type;
  enum BORDER border;
  enum TREETYPE treeType;
  UINT4 level;              
  UINT4 HPFilterLength;
  UINT4 LPFilterLength;
  REAL4TimeSeries *data;

  REAL8 *PForward;
  REAL8 *PInverse;
  REAL8 *UForward;
  REAL8 *UInverse;

  REAL8 *pLForward;
  REAL8 *pLInverse;
  REAL8 *pHForward;
  REAL8 *pHInverse;
}
Wavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{PixelWavelet}}
\index{\texttt{PixelWavelet}}
\noindent PixelWavelet describes pixel in two dimensional time frequency spectrogram.
\begin{description}
\item[\texttt{UINT4 time}] Time index.
\item[\texttt{UINT4 frequency}] Frequency index.
\item[\texttt{UINT4 clusterID}] Id of the cluster the pixel belongs to.
\item[\texttt{BOOLEAN core}] Defines if the pixel is in the core or in the halo of the cluster.
\item[\texttt{UINT4 neighbors[8]}] List of neighbors indexed the same way as in \texttt{pMask}.
\item[\texttt{UINT4 neighborsCount}] The number of neighbors.
\item[\texttt{REAL4 amplitude}] Percentile amplitude.
\item[\texttt{REAL4 amplitudeOriginal}] Original amplitude.
%\item[\texttt{}] .

\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagPixelWavelet
{
  UINT4 time;
  UINT4 frequency;
  UINT4 clusterID;
  BOOLEAN core;
  UINT4 neighbors[8];
  UINT4 neighborsCount;
  REAL4 amplitude;/*percentile*/
  REAL4 amplitudeOriginal;
}
PixelWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{ClusterBlobWavelet}}
\index{\texttt{ClusterBlobWavelet}}
\noindent ClusterBlobWavelet describes blobs of percentile and original amplitudes written to the database.
\begin{description}
\item[\texttt{UINT4 start\_time\_indx}]  Time index where the cluster starts.
\item[\texttt{UINT4 stop\_time\_indx}] Time index where the cluster ends.
\item[\texttt{UINT4 start\_freq\_indx}] Minimum frequency index of the cluster.
\item[\texttt{UINT4 stop\_freq\_indx}] Maximum frequency index of the cluster.
\item[\texttt{UINT4 time\_width}] $stop\_time\_indx-start\_time\_indx+1$.
\item[\texttt{UINT4 freq\_width}] $stop\_freq\_indx-start\_freq\_indx+1$.
\item[\texttt{REAL4 *pBlob}] Array of percentile amplitudes for the smallest rectangle containing the cluster.
\item[\texttt{REAL4 *oBlob}] Array of original amplitudes for the smallest rectangle containing the cluster.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagClusterBlobWavelet
{
  UINT4 start_time_indx;
  UINT4 stop_time_indx;
  UINT4 start_freq_indx;
  UINT4 stop_freq_indx;
  UINT4 time_width;
  UINT4 freq_width;
  REAL4 *pBlob;
  REAL4 *oBlob;
}
ClusterBlobWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{ClusterWavelet}}
\index{\texttt{ClusterBlobWavelet}}
\noindent ClusterWavelet organizes non-zero pixels into clusters and contains various characteristics of those clusters.
\begin{description}
\item[\texttt{Wavelet *wavelet}] The wavelet for which clustering is done.
\item[\texttt{REAL4 *medians}] Medians of the amplitudes on each layer.
\item[\texttt{UINT4 pMaskCount}] Number of non-zero pixels that survived after percentile transform and coincidence.
\item[\texttt{UINT4 clusterCount}] Number of clusters after percentile transform, coincidence and initial selection on cluster size.
\item[\texttt{UINT4 clusterCountFinal}] Final number of clusters that survived all the selection criteria.
\item[\texttt{INT4 clusterType}] Currently possible values are: ORIGINAL\_CL, PIXELMIXING\_CL, SWAP\_CL, SIMULATION\_CL. 
 At the moment this field is not yet actively used
\item[\texttt{REAL4 delta\_t}] Wavelet time step.
\item[\texttt{REAL4 delta\_f}] Wavelet frequency step.
\item[\texttt{PixelWavelet **pMask}] An array of all non-zero pixels (after percentilie transform and coincidence).
\item[\texttt{UINT4 *sCuts}] An array elements of which can take two values: 0 (the cluster has not passed some selection criteria) and 1.
\item[\texttt{UINT4 **cList}] An array of arrays that list pixels (in terms of pMask index) belonging to the corresponding cluster.
\item[\texttt{UINT4 *volumes}] An array of volumes ($coreSize + halo$) of clusters.
\item[\texttt{UINT4 *coreSize}] An array of core sizes of clusters.
\item[\texttt{REAL4 *correlation}] An array of correlations (asymmetries) for each cluster computed as follows: 
$\frac{NumberOfPositivePercentileAmplitudes-NumberOfNegativePercentileAmplitudes}{coreSize}$.
\item[\texttt{REAL4 *likelihood}] Likelihood of each cluster.
\item[\texttt{REAL4 *power}] Power of each cluster.
\item[\texttt{REAL4 *maxAmplitude}] Absolute maximum value of a pixel amplitude for each cluster.
\item[\texttt{REAL8 *relativeStartTime}] Start time of the cluster relative to the beginning of time series given as an input.
\item[\texttt{REAL8 *relativeStopTime}] Stop time of the cluster relative to the beginning of time series given as an input.
\item[\texttt{REAL8 *duration}] $relativeStopTime-relativeStartTime$.
\item[\texttt{LIGOTimeGPS *absoluteStartTime}] Absolute start time of the cluster expressed as GPS time.
\item[\texttt{LIGOTimeGPS *absoluteStopTime}] Absolute stop time of the cluster expressed as GPS time.
\item[\texttt{REAL4 *startFrequency}] Minimum frequency of the cluster.
\item[\texttt{REAL4 *stopFrequency}] Maximum frequency of the cluster.
\item[\texttt{REAL4 *bandwidth}] $stopFrequency-startFrequency$.
\item[\texttt{ClusterBlobWavelet *blobs}] The actual percentile or original amplitudes describing the cluster.
\item[\texttt{REAL4 nonZeroFractionAfterPercentile}] Fraction of non-zero pixels after percentile transform.
\item[\texttt{REAL4 nonZeroFractionAfterCoincidence}] Fraction of non-zero pixels after coincidence.
\item[\texttt{REAL4 nonZeroFractionAfterSetMask}] Fraction of non-zero pixels after setMask was applied and small cluster were zeroed out.
%\item[\texttt{}] .
%\item[\texttt{}] .
%\item[\texttt{}] .
\item[\texttt{BOOLEAN pixelSwapApplied}] The flag that indicates whether pixel swap was applied.
\item[\texttt{BOOLEAN pixelMixerApplied}] The flag that indicates whether pixel mixing was applied.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct 
tagClusterWavelet
{
  Wavelet *wavelet;
  Wavelet *original;

  REAL4FrequencySeries *psd;

  REAL4 *medians;
  REAL4 *norm50;
  REAL4 *avgPSD;

  REAL4 calibration_max_freq;

  UINT4 pMaskCount;
  UINT4 clusterCount;
  UINT4 clusterCountFinal;
  CLUSTER_TYPE clusterType;
  INT4 simulationType;
  REAL4 delta_t;
  REAL4 delta_f;

  PixelWavelet **pMask;

  UINT4 *sCuts;  
  UINT4 **cList;
  UINT4 *volumes;

  UINT4 *coreSize;
  REAL4 *correlation;
  REAL4 *likelihood;
  REAL4 *power;
  REAL4 *maxAmplitude;
  REAL8 *relativeStartTime;
  REAL8 *relativeStopTime;
  REAL8 *duration;
  LIGOTimeGPS *absoluteStartTime;
  LIGOTimeGPS *absoluteStopTime;
  REAL4 *startFrequency;
  REAL4 *stopFrequency;
  REAL4 *bandwidth;
  REAL4 *noise_rms;
  BOOLEAN noise_rms_flag;

  ClusterBlobWavelet *blobs;
  
  REAL4 nonZeroFractionAfterPercentile;
  REAL4 nonZeroFractionAfterCoincidence;
  REAL4 nonZeroFractionAfterSetMask;
  REAL4 nonZeroFractionAfterClustering;
  REAL4 nonZeroFractionAfterCuts;
  REAL4 nonZeroFractionAfterVetoes;

  BOOLEAN pixelSwapApplied;
  BOOLEAN pixelMixerApplied;
}
ClusterWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputLayerWavelet}}
\index{\texttt{InputLayerWavelet}}
\noindent \texttt{InputLayerWavelet} contains input for \texttt{LALGetLayerWavelet} function.
\begin{description}
\item[\texttt{Wavelet *wavelet}] Wavelet.
\item[\texttt{UINT4 index}] Index of a layer to extract.
%\item[\texttt{}] .
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputLayerWavelet
{
  Wavelet *wavelet;
  UINT4 index;
}
InputLayerWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputLayerWavelet}}
\index{\texttt{OutputLayerWavelet}}
\noindent \texttt{OutputLayerWavelet} contains output for \texttt{LALGetLayerWavelet} function.
\begin{description}
\item[\texttt{REAL4TimeSeries *layer}] The extracted layer.
\item[\texttt{INT4 status}] Status. Not used at the moment.
%\item[\texttt{}] .
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputLayerWavelet
{
  INT4 status;
  REAL4TimeSeries *layer;
}
OutputLayerWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputGetMaxLayerWavelet}}
\index{\texttt{InputGetMaxLayerWavelet}}
\noindent \texttt{InputGetMaxLayerWavelet} contains input for \texttt{LALGetMaxLayerWavelet} function.
\begin{description}
\item[\texttt{Wavelet *wavelet}] Wavelet.
%\item[\texttt{}] .
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputGetMaxLayerWavelet
{
  Wavelet *wavelet;
}
InputGetMaxLayerWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputGetMaxLayerWavelet}}
\index{\texttt{OutputGetMaxLayerWavelet}}
\noindent \texttt{OutputGetMaxLayerWavelet} contains output for \texttt{LALGetMaxLayerWavelet} function.
\begin{description}
\item[\texttt{UINT4 maxLayer}] Maximum layer.
%\item[\texttt{}] .
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputGetMaxLayerWavelet
{
  UINT4 maxLayer;
}
OutputGetMaxLayerWavelet;

/* typedef struct */
/* tagInputFractionWavelet */
/* { */
/*   Wavelet *in; */
/*   INT4 nF; */
/*   INT4 nL; */
/*   INT4 dim; */
/*   UINT4 seed; */
/* } */
/* InputFractionWavelet; */

/* typedef struct */
/* tagOutputFractionWavelet */
/* { */
/*   ClusterWavelet *out; */
/*   REAL4 nonZeroFraction; */
/* } */
/* OutputFractionWavelet; */


/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputPercentileWavelet}}
\index{\texttt{InputPercentileWavelet}}
\noindent \texttt{InputPercentileWavelet} contains input for \texttt{LALPercentileWavelet} function.
\begin{description}
\item[\texttt{Wavelet *in}] Original wavelet.
\item[\texttt{REAL4 nonZeroFraction}] Fraction of pixels that survive percentil transform.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputPercentileWavelet
{
  Wavelet *in;
  REAL4 nonZeroFraction;
}
InputPercentileWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputPercentileWavelet}}
\index{\texttt{OutputPercentileWavelet}}
\noindent \texttt{OutputPercentileWavelet} contains output for \texttt{LALPercentileWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *out}] Transformed wavelet embedded into \texttt{ClusterWavelet} structure.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputPercentileWavelet
{
  ClusterWavelet *out;
}
OutputPercentileWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputPixelSwapWavelet}}
\index{\texttt{InputPixelSwapWavelet}}
\noindent \texttt{InputPixelSwapWavelet} contains input for \texttt{LALPixelSwapWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *in}] Original wavelet.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputPixelSwapWavelet
{
  ClusterWavelet *in;
}
InputPixelSwapWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputPixelSwapWavelet}}
\index{\texttt{OutputPixelSwapWavelet}}
\noindent \texttt{OutputPixelSwapWavelet} contains output for \texttt{LALPixelSwapWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *out}] Wavelet in which each layer is swapped in time. 
It is used to estimate the accidental coincidence of glitches in two channels.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputPixelSwapWavelet
{
  ClusterWavelet *out;
}
OutputPixelSwapWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputPixelMixerWavelet}}
\index{\texttt{InputPixelMixerWavelet}}
\noindent \texttt{InputPixelMixerWavelet} contains input for \texttt{LALPixelMixerWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *in}] Original wavelet.
\item[\texttt{UINT4 seed}] Seed for random number generator.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputPixelMixerWavelet
{
  ClusterWavelet *in;
  RandomParams *rparams;
}
InputPixelMixerWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputPixelMixerWavelet}}
\index{\texttt{OutputPixelMixerWavelet}}
\noindent \texttt{OutputPixelMixerWavelet} contains output for \texttt{LALPixelMixerWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *out}] Wavelet in which pixels are randomly mixed. 
It is used to estimate the accidental coincidence of glitches in two channels.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputPixelMixerWavelet
{
  ClusterWavelet *out;
}
OutputPixelMixerWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputCoincidenceWavelet}}
\index{\texttt{InputCoincidenceWavelet}}
\noindent \texttt{InputCoincidenceWavelet} contains input for \texttt{LALCoincidenceWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *one}] Wavelet from channel one.
\item[\texttt{ClusterWavelet *two}] Wavelet from channel two.
\item[\texttt{UINT4 timeWindowNanoSec}] Maximum time shift in nanoseconds between two channels.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputCoincidenceWavelet
{
  ClusterWavelet *one;
  ClusterWavelet *two;
  INT4 timeWindowPixels;
  INT4 freqWindowPixels;
  INT4 coincidenceLevel;
  REAL4 minAmp4SinglePixels;
  REAL4 minAmp4ClusterExtension;
}
InputCoincidenceWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputCoincidenceWavelet}}
\index{\texttt{OutputCoincidenceWavelet}}
\noindent \texttt{OutputCoincidenceWavelet} contains output for \texttt{LALCoincidenceWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *one}] Wavelet from channel one filtered by coincidence with channel two.
\item[\texttt{ClusterWavelet *two}] Wavelet from channel two filtered by coincidence with channel one.
\item[\texttt{REAL4 occupancyOne}] The resulting fraction of non-zero pixels in channel one.
\item[\texttt{REAL4 occupancyTwo}] The resulting fraction of non-zero pixels in channel two.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputCoincidenceWavelet
{
  ClusterWavelet *one;
  ClusterWavelet *two;
  REAL4 occupancyOne;
  REAL4 occupancyTwo;
}
OutputCoincidenceWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{InputClusterWavelet}}
\index{\texttt{InputClusterWavelet}}
\noindent \texttt{InputClusterWavelet} contains input for \texttt{LALClusterWavelet} function.
\begin{description}
\item[\texttt{BOOLEAN aura}] This flag defines if we should use halo when clustering.
\item[\texttt{UINT4 minClusterSize}] Minimum cluster size.
\item[\texttt{UINT4 maxClusterSize}] Maximum cluster size.
\item[\texttt{ClusterWavelet *w}] Wavelet to cluster (with percentile amplitudes, after coincidence).
\item[\texttt{Wavelet *original}] Original wavelet.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagInputClusterWavelet
{
  BOOLEAN aura;
  UINT4 minClusterSize;
  UINT4 maxClusterSize;
  ClusterWavelet *w;
  /*Wavelet *original;*/
}
InputClusterWavelet;

/*************************************<lalLaTeX file="WaveburstStructs">
\subsubsection*{struct \texttt{OutputClusterWavelet}}
\index{\texttt{OutputClusterWavelet}}
\noindent \texttt{OutputClusterWavelet} contains output for \texttt{LALClusterWavelet} function.
\begin{description}
\item[\texttt{ClusterWavelet *w}] The resulting set of clusters.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct
tagOutputClusterWavelet
{
  ClusterWavelet *w;
}
OutputClusterWavelet;


typedef struct
tagInputReuseClusterWavelet
{
  BOOLEAN aura;
  UINT4 minClusterSize;
  UINT4 maxClusterSize;
  ClusterWavelet *w;
  ClusterWavelet *another;
  /*  Wavelet *original;
      Wavelet *originalAnother;*/
}
InputReuseClusterWavelet;


typedef struct
tagInputGetClusterParameters
{
  ClusterWavelet *w;
}
InputGetClusterParameters;

typedef struct
tagOutputGetClusterParameters
{
  int nothing;
}
OutputGetClusterParameters;

typedef struct 
tagInputForwardWavelet
{
  Wavelet *w;
  INT4 level;
  INT4 layer;
}
InputForwardWavelet;

typedef struct 
tagInputInverseWavelet
{
  Wavelet *w;
  INT4 level;
  INT4 layer;
}
InputInverseWavelet;

typedef struct tagInputt2wWavelet
{
  Wavelet *w;
  INT4 ldeep;
}
Inputt2wWavelet;

typedef struct tagOutputt2wWavelet
{
  Wavelet *w;
}
Outputt2wWavelet;

typedef struct tagInputw2tWavelet
{
  Wavelet *w;
  INT4 ldeep;
} 
Inputw2tWavelet;

typedef struct tagOutputw2tWavelet
{
  Wavelet *w;
}
Outputw2tWavelet;

void
LALGetLayerWavelet(LALStatus *status,
		   OutputLayerWavelet **output,
		   InputLayerWavelet *input);


void
LALGetMaxLayerWavelet(LALStatus *status,
		      OutputGetMaxLayerWavelet **output,
		      InputGetMaxLayerWavelet *input);


/* void */
/* LALFractionWavelet( LALStatus *status, */
/* 		    OutputFractionWavelet **output, */
/* 		    InputFractionWavelet  *input); */


void
LALPercentileWavelet( LALStatus *status,
		      OutputPercentileWavelet **output,
		      InputPercentileWavelet  *input);


void
LALPixelSwapWavelet(LALStatus *status,
		    OutputPixelSwapWavelet **output,
		    InputPixelSwapWavelet *input);

void
LALPixelMixerWavelet(LALStatus *status,
		     OutputPixelMixerWavelet **output,
		     InputPixelMixerWavelet *input);

void
LALCoincidenceWavelet(LALStatus *status,
		      OutputCoincidenceWavelet **output,
		      InputCoincidenceWavelet *input);

void
LALClusterWavelet(LALStatus *status,
		  OutputClusterWavelet **output,
		  InputClusterWavelet *input);

void
LALReuseClusterWavelet(LALStatus *status,
                       OutputClusterWavelet **output,
                       InputReuseClusterWavelet *input);

void
LALSetAmplitudesWavelet(LALStatus *status,
			ClusterWavelet *w);

void
LALAssignREAL4TimeSeries(LALStatus *status,
			 REAL4TimeSeries **left,
			 REAL4TimeSeries *right);

void
LALAllocateWavelet(LALStatus *status,
		   Wavelet **wavelet);
void
LALFreeWavelet(LALStatus *status,
	       Wavelet **wavelet);

void
LALFreeREAL4TimeSeries(LALStatus *status,
		       REAL4TimeSeries **t);

void
LALFreeClusterWavelet(LALStatus *status,
		      ClusterWavelet **w);

void
LALFreeOutPercentile(LALStatus *status,
		     OutputPercentileWavelet **p);

void
LALFreeOutCoincidence(LALStatus *status,
		      OutputCoincidenceWavelet **co);

void
LALFreeOutCluster(LALStatus *status,
		  OutputClusterWavelet **cl);

void LALForwardWavelet(LALStatus *status,
		       InputForwardWavelet *input);

void LALInverseWavelet(LALStatus *status,
		       InputInverseWavelet *input);

void LALt2wWavelet(LALStatus *status,
		   Inputt2wWavelet *input,
		   Outputt2wWavelet **output);

void LALw2tWavelet(LALStatus *status,
		   Inputw2tWavelet *input,
		   Outputw2tWavelet **output);

#endif
