/********************************* <lalVerbatim file="LALWaveletHV">
Author: Klimenko, Sergey and Yakushin, Igor
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{LALWavelet.h}}

[One sentence briefly defining scope of the header]

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALWavelet.h>
\end{verbatim}

[Generic documentation on the header; this is the main place to
document any stuff not specific to the module]

\subsection*{Error conditions}
\input{LALWaveletHE}

\subsection*{Structures}

[Document here any structures defined in the header.  
Also include any of them in the index; e.g.:]
% \index{\texttt{LALWaveletOutput}}
% \index{\texttt{LALWaveletInput}}
% \index{\texttt{LALWaveletParams}}

\vfill{\footnotesize\input{LALWaveletHV}}
\newpage\input{LALWaveletC}
\newpage\input{clusterTestC.tex}
\newpage\input{LALGetLayerWaveletTestC.tex}
\newpage\input{setMaskTestC.tex}
********************************** </lalLaTeX> */

#ifndef _LALWAVELET_H
#define _LALWAVELET_H

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>

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

enum COINCIDENCETYPE {GG,GV};


typedef struct
tagSlice
{
  /*  enum SLICETYPE type;*/
  UINT4 start;
  UINT4 size;
  UINT4 step;
}
Slice;

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
}
Wavelet;

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

typedef struct 
tagClusterWavelet
{
  Wavelet *wavelet;

  REAL4 *medians;

  UINT4 pMaskCount;
  UINT4 clusterCount;
  UINT4 clusterCountFinal;
  INT4 clusterType;
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


typedef struct
tagInputLayerWavelet
{
  Wavelet *wavelet;
  UINT4 index;
}
InputLayerWavelet;

typedef struct
tagOutputLayerWavelet
{
  INT4 status;
  REAL4TimeSeries *layer;
}
OutputLayerWavelet;

typedef struct
tagInputGetMaxLayerWavelet
{
  Wavelet *wavelet;
}
InputGetMaxLayerWavelet;

typedef struct
tagOutputGetMaxLayerWavelet
{
  UINT4 maxLayer;
}
OutputGetMaxLayerWavelet;

typedef struct
tagInputFractionWavelet
{
  Wavelet *in;
  INT4 nF;
  INT4 nL;
  INT4 dim;
  UINT4 seed;
}
InputFractionWavelet;

typedef struct
tagOutputFractionWavelet
{
  ClusterWavelet *out;
  REAL4 nonZeroFraction;
}
OutputFractionWavelet;


typedef struct
tagInputPercentileWavelet
{
  Wavelet *in;
  REAL4 nonZeroFraction;
}
InputPercentileWavelet;

typedef struct
tagOutputPercentileWavelet
{
  ClusterWavelet *out;
}
OutputPercentileWavelet;

typedef struct
tagInputPixelSwapWavelet
{
  ClusterWavelet *in;
}
InputPixelSwapWavelet;

typedef struct
tagOutputPixelSwapWavelet
{
  ClusterWavelet *out;
}
OutputPixelSwapWavelet;

typedef struct
tagInputPixelMixerWavelet
{
  ClusterWavelet *in;
  UINT4 seed;
}
InputPixelMixerWavelet;

typedef struct
tagOutputPixelMixerWavelet
{
  ClusterWavelet *out;
}
OutputPixelMixerWavelet;

typedef struct
tagInputCoincidenceWavelet
{
  ClusterWavelet *one;
  ClusterWavelet *two;
  UINT4 timeWindowNanoSec;
}
InputCoincidenceWavelet;

typedef struct
tagOutputCoincidenceWavelet
{
  ClusterWavelet *one;
  ClusterWavelet *two;
  REAL4 occupancyOne;
  REAL4 occupancyTwo;
}
OutputCoincidenceWavelet;

typedef struct
tagInputClusterWavelet
{
  BOOLEAN aura;
  UINT4 minClusterSize;
  UINT4 maxClusterSize;
  ClusterWavelet *w;
  Wavelet *original;
}
InputClusterWavelet;

typedef struct
tagOutputClusterWavelet
{
  ClusterWavelet *w;
}
OutputClusterWavelet;

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


void
LALGetLayerWavelet(LALStatus *status,
		   OutputLayerWavelet **output,
		   InputLayerWavelet *input);


void
LALGetMaxLayerWavelet(LALStatus *status,
		      OutputGetMaxLayerWavelet **output,
		      InputGetMaxLayerWavelet *input);


void
LALFractionWavelet( LALStatus *status,
		    OutputFractionWavelet **output,
		    InputFractionWavelet  *input);


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


#endif
