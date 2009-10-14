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

#ifndef _COMPUTEFSTATREAL4OPENCL_H
#define _COMPUTEFSTATREAL4OPENCL_H

/* ---------- exported defines and macros ---------- */
#define NUM_IFOS 2           /**< number of detectors */
#define MAX_NUM_SFTS 64      /**< max number of SFTs  */

#if (USE_OPENCL_KERNEL && USE_OPENCL_KERNEL_CPU)
#error Cannot use both USE_OPENCL_KERNEL and USE_OPENCL_KERNEL_CPU at the same time
#endif

#if USE_OPENCL_KERNEL
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif 
#else
#define cl_mem void*
#endif

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

    REAL42MemBuffer    Freq;
    REAL4MemBuffer     fkdot16;
    REAL4MemBuffer     tSSB_DeltaT_int;
    REAL4MemBuffer     tSSB_DeltaT_rem;
    REAL4MemBuffer     tSSB_TdotM1;

    REAL4MemBuffer     amcoe_a;
    REAL4MemBuffer     amcoe_b;

    REAL44MemBuffer    ABCInvD;

} CLWorkspace;

/* ---------- exported global variables ---------- */
extern const CLWorkspace empty_CLWorkspace;
extern CLWorkspace clW;

/* ---------- additional exported API prototypes ---------- */
extern int XLALInitCLWorkspace (CLWorkspace *clWp, const MultiSFTVectorSequence *stackMultiSFT);
extern void XLALDestroyCLWorkspace (CLWorkspace *clWp, const MultiSFTVectorSequence *stackMultiSFT);

#if USE_OPENCL_KERNEL
void freeCLMemoryObject (cl_mem *memobj);
#endif

#endif  /* Double-include protection. */
