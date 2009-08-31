/*
*  Copyright (C) 2009 Reinhard Prix, Oleg Korobkin, Bernd Machenschalk
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

#ifndef _COMPUTEFSTATREAL4_H
#define _COMPUTEFSTATREAL4_H

/* used in HierarchicalSearch.c */
#ifdef USE_CUDA
#define GPUREADY_DEFAULT 1
#define INITIALIZE_COPROCESSOR_DEVICE
#define UNINITIALIZE_COPROCESSOR_DEVICE
#define REARRANGE_SFT_DATA
#endif /*  USE_CUDA */

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/* ---------- includes ---------- */
#include <lal/LALDatatypes.h>
#include <lal/DetectorStates.h>
#include <lal/ComputeFstat.h>

/* ---------- exported defines and macros ---------- */
/* used in OpenCL */
#define NUM_IFOS 2           /**< number of detectors */
#define MAX_NUM_SFTS 64      /**< max number of SFTs  */

/* used in OpenCL */
#ifndef USE_OPENCL_KERNEL
#define USE_OPENCL_KERNEL 0
#endif
#ifndef USE_OPENCL_KERNEL_CPU
#define USE_OPENCL_KERNEL_CPU 0
#endif

#if (USE_OPENCL_KERNEL && USE_OPENCL_KERNEL_CPU)
#error Cannot use both USE_OPENCL_KERNEL and USE_OPENCL_KERNEL_CPU at the same time
#endif

/** ----- switch between different versions of XLALComputeFStatFreqBandVector() ----- */
#ifdef USE_CUDA
#define XLALComputeFStatFreqBandVector XLALComputeFStatFreqBandVectorCUDA
#else
#define XLALComputeFStatFreqBandVector XLALComputeFStatFreqBandVectorCPU
#endif

/* used in OpenCL */
#if (USE_OPENCL_KERNEL || USE_OPENCL_KERNEL_CPU)
#define INITIALIZE_COPROCESSOR_DEVICE   CLWorkspace clW = empty_CLWorkspace; XLALInitCLWorkspace (&clW, &stackMultiSFT); 
#define UNINITIALIZE_COPROCESSOR_DEVICE XLALDestroyCLWorkspace ( &clW, &stackMultiSFT );
#define REARRANGE_SFT_DATA              XLALRearrangeSFTData ( &clW, &fstatVector );
#endif

#if USE_OPENCL_KERNEL
#ifdef _WIN32
#include <CL/cl.h>
#endif
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif 
#else
#define cl_mem void*
#endif

/* ---------- exported types ---------- */
// FIXME: the following types might be useful to move into LAL at some point
/** sequence of MultiSFT vectors -- for each segment */
typedef struct tagMultiSFTVectorSequence {
  UINT4 length;     		/**< number of segments */
  MultiSFTVector **data; 	/**< the SFT vectors */
} MultiSFTVectorSequence;

/** sequence of Multi-noise weights vectors -- for each segment */
typedef struct tagMultiNoiseWeightsSequence {
  UINT4 length;     		/**< number of segments */
  MultiNoiseWeights **data; 	/**< the noise weights */
} MultiNoiseWeightsSequence;

/** sequence of Multi-detector state vectors -- for each segment */
typedef struct tagMultiDetectorStateSeriesSequence {
  UINT4 length;     		/**< number of segments */
  MultiDetectorStateSeries **data; /**< the detector state series */
} MultiDetectorStateSeriesSequence;


/** REAL4 version of pulsar spins fkdot[] */
typedef struct {
  REAL4 FreqMain;			/**< "main" part of frequency fkdot[0], normally just the integral part */
  REAL4 fkdot[PULSAR_MAX_SPINS];	/**< remaining spin-parameters, including *fractional* part of Freq = fkdot[0] */
  UINT4 spdnOrder;			/**< highest non-zero spindown order (ie maximal value of 'k' in fkdot = df^(k)/dt^k */
} PulsarSpinsREAL4;


/** Simple container for REAL4-vectors, holding the SSB-timings DeltaT_alpha  and Tdot_alpha,
 *  with one entry per SFT-timestamp. We also store the SSB reference-time tau0.
 * NOTE: this is a REAL4 version of SSBtimes, preserving the required precision by appropriate
 * 'splitting' of REAL8's into pairs of REAL4s.
 */
typedef struct {
  LIGOTimeGPS refTime;		/**< Reference time wrt to which the time-differences DeltaT are computed */
  REAL4Vector *DeltaT_int;	/**< Integral part of time-difference of SFT-alpha - tau0 in SSB-frame */
  REAL4Vector *DeltaT_rem;	/**< Remainder of time-difference of SFT-alpha - tau0 in SSB-frame */
  REAL4Vector *TdotM1;		/**< dT/dt - 1 : time-derivative of SSB-time wrt local time for SFT-alpha, of order O(1e-4) */
} SSBtimesREAL4;


/** Multi-IFO container for SSB timings in REAL4-representation */
typedef struct {
  UINT4 length;			/**< number of IFOs */
  SSBtimesREAL4 **data;		/**< array of SSBtimes (pointers) */
} MultiSSBtimesREAL4;


/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb (for ML-estimators).
 * NOTE: this is simply a REAL4 version of Fcomponents.
 */
typedef struct {
  REAL4 F;		/**< F-statistic value */
  COMPLEX8 Fa;		/**< complex amplitude Fa */
  COMPLEX8 Fb;		/**< complex amplitude Fb */
} FcomponentsREAL4;

/** Struct holding buffered XLALDriverFstatREAL4()-internal quantities
 * to avoid unnecessarily recomputing things that depend ONLY on the skyposition and detector-state series
 * (but not on the spins).
 */
typedef struct {
  REAL8 Alpha, Delta;					/**< target skyposition of previous search */
  const MultiDetectorStateSeries *multiDetStates;	/**< input detStates series used in previous search */
  MultiSSBtimesREAL4 *multiSSB;				/**< SSB timings computed in previous search */
  MultiAMCoeffs *multiAMcoef;				/**< antenna-pattern coeffs computed in previous search */
} ComputeFBufferREAL4;

/** Struct holding buffered XLALComputeFStatFreqBandVector()-internal quantities
 * to avoid unnecessarily recomputing things that depend ONLY on the skyposition and detector-state series
 * (but not on the spins).
 */
typedef struct {
  REAL8 Alpha, Delta;						/**< target skyposition of previous search */
  const MultiDetectorStateSeriesSequence *multiDetStatesV;	/**< input detStates series used in previous search */
  UINT4 numSegments;						/**< number of segments */
  MultiSSBtimesREAL4 **multiSSB4V;				/**< array[numSegments] of SSB timings computed in previous search */
  MultiAMCoeffs **multiAMcoefV;					/**< array[numSegments] of antenna-pattern coeffs computed in previous search */
} ComputeFBufferREAL4V;

#if (USE_OPENCL_KERNEL || USE_OPENCL_KERNEL_CPU)

typedef struct {
  REAL4 FreqMain;
  REAL4 fkdot0;
} REAL42;

typedef struct {
  REAL4 Ad;
  REAL4 Bd;
  REAL4 Cd;
  REAL4 InvDd;
} REAL44;

typedef struct {
  REAL4 s[16];   // HARDCODED VALUE
} PulsarSpins16;

/** Memory buffer structures to group logically the buffers on the host and
 * memory objects on the device */
typedef struct tagUINT4MemBuffer {
  UINT4 length;
  UINT4 *data;
  cl_mem memobj;
} UINT4MemBuffer;

typedef struct tagREAL4MemBuffer {
  UINT4 length;
  REAL4 *data;
  cl_mem memobj;
} REAL4MemBuffer;

typedef struct tagCOMPLEX8MemBuffer {
  UINT4 length;
  COMPLEX8 *data;
  cl_mem memobj;
} COMPLEX8MemBuffer;

typedef struct tagREAL42MemBuffer {
  UINT4 length;
  REAL42 *data;
  cl_mem memobj;
} REAL42MemBuffer;

typedef struct tagREAL44MemBuffer {
  UINT4 length;
  REAL44 *data;
  cl_mem memobj;
} REAL44MemBuffer;

/** Struct to store OpenCL context: platform, queue, kernel etc.
 * It can be declared in a main rouine of an application, i.e. HierarchicalSearch, 
 * and then filled in in top-level function XLALComputeFStatFreqBandVector()
 */
typedef struct {

#if USE_OPENCL_KERNEL
    cl_platform_id    *platform;
    cl_device_id      *device;
    cl_context        *context;
    cl_command_queue  *cmd_queue;
    cl_program        *program;
    cl_kernel         *kernel;
#endif

    UINT4              numBins;       // <200
    UINT4              numSegments;   // ~121
    UINT4              numIFOs;       // 2
    UINT4              maxNumSFTs;    // <64
    UINT4              sftLen;        // ~254 -- can vary between the runs and between the calls

    REAL4MemBuffer     Fstat;

    COMPLEX8MemBuffer  multiSFTsFlat; /* flattened array for all SFT COMPLEX8 data */
    UINT4MemBuffer     numSFTsV;      /* numSFTs for numSegments x numIFOs data slots */

    REAL42MemBuffer    fkdot4;
    REAL4MemBuffer     tSSB_DeltaT_int;
    REAL4MemBuffer     tSSB_DeltaT_rem;
    REAL4MemBuffer     tSSB_TdotM1;

    REAL4MemBuffer     amcoe_a;
    REAL4MemBuffer     amcoe_b;

    REAL44MemBuffer    ABCInvD;

} CLWorkspace;

#endif

/* ---------- exported global variables ---------- */
extern const ComputeFBufferREAL4 empty_ComputeFBufferREAL4;
extern const ComputeFBufferREAL4V empty_ComputeFBufferREAL4V;
#if (USE_OPENCL_KERNEL || USE_OPENCL_KERNEL_CPU)
extern const CLWorkspace empty_CLWorkspace;
#endif

/* ---------- exported API prototypes ---------- */
int
XLALComputeFStatFreqBandVector ( REAL4FrequencySeriesVector *fstatBandV,
				 const PulsarDopplerParams *doppler,
				 const MultiSFTVectorSequence *multiSFTsV,
				 const MultiNoiseWeightsSequence *multiWeightsV,
				 const MultiDetectorStateSeriesSequence *multiDetStatesV,
				 UINT4 Dterms,
				 ComputeFBufferREAL4V *cfvBuffer
				 );

#if (USE_OPENCL_KERNEL || USE_OPENCL_KERNEL_CPU)
int
XLALComputeFStatFreqBandVectorOpenCL ( REAL4FrequencySeriesVector *fstatBandV,
				       const PulsarDopplerParams *doppler,
				       const MultiSFTVectorSequence *multiSFTsV,
				       const MultiNoiseWeightsSequence *multiWeightsV,
				       const MultiDetectorStateSeriesSequence *multiDetStatesV,
				       UINT4 Dterms,
				       ComputeFBufferREAL4V *cfvBuffer,
				       CLWorkspace *clW
				       );
#endif

int
XLALDriverFstatREAL4 ( REAL4 *Fstat,
                     const PulsarDopplerParams *doppler,
                     const MultiSFTVector *multiSFTs,
                     const MultiNoiseWeights *multiWeights,
                     const MultiDetectorStateSeries *multiDetStates,
                     UINT4 Dterms,
                     ComputeFBufferREAL4 *cfBuffer
                     );

void
XLALCoreFstatREAL4 (REAL4 *Fstat,
                  PulsarSpinsREAL4 *fkdot4,
                  const MultiSFTVector *multiSFTs,
                  MultiSSBtimesREAL4 *multiSSB4,
                  MultiAMCoeffs *multiAMcoef,
                  UINT4 Dterms
                  );


void
XLALComputeFaFbREAL4 ( FcomponentsREAL4 *FaFb,
                       const SFTVector *sfts,
                       const PulsarSpinsREAL4 *fkdot4,
                       const SSBtimesREAL4 *tSSB,
                       const AMCoeffs *amcoe,
                       UINT4 Dterms);


MultiSSBtimesREAL4 *
XLALGetMultiSSBtimesREAL4 ( const MultiDetectorStateSeries *multiDetStates,
                            REAL8 Alpha, REAL8 Delta,
                            LIGOTimeGPS refTime
                            );

SSBtimesREAL4 *
XLALGetSSBtimesREAL4 ( const DetectorStateSeries *DetectorStates,
                       REAL8 Alpha, REAL8 Delta,
                       LIGOTimeGPS refTime
                       );

void XLALDestroySSBtimesREAL4 ( SSBtimesREAL4 *tSSB );
void XLALDestroyMultiSSBtimesREAL4 ( MultiSSBtimesREAL4 *multiSSB );

void XLALEmptyComputeFBufferREAL4 ( ComputeFBufferREAL4 *cfb);
void XLALEmptyComputeFBufferREAL4V ( ComputeFBufferREAL4V *cfbv);

#if (USE_OPENCL_KERNEL || USE_OPENCL_KERNEL_CPU)
int XLALInitCLWorkspace (CLWorkspace *clW, const MultiSFTVectorSequence *stackMultiSFT);
void XLALDestroyCLWorkspace (CLWorkspace *clW, const MultiSFTVectorSequence *stackMultiSFT);
#endif

#if USE_OPENCL_KERNEL
void freeCLMemoryObject (cl_mem *memobj);
#endif

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
