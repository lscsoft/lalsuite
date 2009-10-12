/*
 * Copyright (C) 2009 Reinhard Prix
 * Copyright (C) 2007 Chris Messenger
 * Copyright (C) 2006 John T. Whelan, Badri Krishnan
 * Copyright (C) 2005, 2006, 2007 Reinhard Prix
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

/** \author R. Prix, J. T. Whelan
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space,
 * following the equations in \ref JKS98.
 *
 * This code is partly a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 * NOTE: this file contains specialized versions of the central CFSv2 routines, that are aimed at GPU optimization.
 * At the moment this only means they are internally using only single precision, but still aggree to within 1%
 * for Tobs ~ 1day and fmax ~ 1kHz.
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/ComputeFstat.h>
#include <lal/LogPrintf.h>

#include "ComputeFstatREAL4.h"

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define LD_SMALL4       ((REAL4)2.0e-4)		/**< "small" number for REAL4*/
#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

/*----- Macros ----- */
#define SQ(x) ( (x) * (x) )
#define REM(x) ( (x) - (INT4)(x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
static const REAL4 inv_fact[PULSAR_MAX_SPINS] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };
CLWorkspace clW;
CLWorkspace *clWp = &clW;
int gpu_device_id = 0; /* might be set on the command-line in hs_boinc_extras.c */
int gpu_platform_id = 0;

/* global sin-cos lookup table */
#define LUT_RES         	64      /* resolution of lookup-table */
#define OO_LUT_RES		(1.0f / LUT_RES )

static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

/* empty initializers  */
static const LALStatus empty_LALStatus;
static const AMCoeffs empty_AMCoeffs;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFBuffer empty_ComputeFBuffer;
const PulsarSpinsREAL4 empty_PulsarSpinsREAL4;
const ComputeFBufferREAL4 empty_ComputeFBufferREAL4;
const ComputeFBufferREAL4V empty_ComputeFBufferREAL4V;
const FcomponentsREAL4 empty_FcomponentsREAL4;
const CLWorkspace empty_CLWorkspace;

/*---------- internal prototypes ----------*/
int finite(double x);
void sin_cos_2PI_LUT_REAL4 (REAL4 *sin2pix, REAL4 *cos2pix, REAL4 x);
void sin_cos_LUT_REAL4 (REAL4 *sinx, REAL4 *cosx, REAL4 x);
void init_sin_cos_LUT_REAL4 (void);

#if USE_OPENCL_KERNEL_CPU
#include "FStatOpenCLKernel.cl"
#endif

/*==================== FUNCTION DEFINITIONS ====================*/

/** REAL4 and GPU-ready version of ComputeFStatFreqBand(), extended to loop over segments as well.
 *
 * Computes a vector of Fstatistic values for a number of frequency bins, for each segment
 */
int
XLALComputeFStatFreqBandVectorOpenCL (   REAL4FrequencySeriesVector *fstatBandV, 		/**< [out] Vector of Fstat frequency-bands */
					 const PulsarDopplerParams *doppler,			/**< parameter-space point to compute F for */
					 const MultiSFTVectorSequence *multiSFTsV, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
					 const MultiNoiseWeightsSequence *multiWeightsV,	/**< noise-weights of all SFTs */
					 const MultiDetectorStateSeriesSequence *multiDetStatesV,/**< 'trajectories' of the different IFOs */
					 UINT4 Dterms,						/**< number of Dirichlet kernel Dterms to use */
					 ComputeFBufferREAL4V *cfvBuffer			/**< buffer quantities that don't need to be recomputed */
					 )
{
  static const char *fn = "XLALComputeFStatFreqBandVector()";

  const UINT4 CALL_XLALCoreFstatREAL4 = 0;
#if USE_OPENCL_KERNEL
  cl_int err, err_total;
#endif
  static int call_count = 0;

  UINT4 numBins, k;
  UINT4 numSegments, n;
  REAL8 f0, deltaF;
  REAL4 Fstat;
  REAL8 Freq0, Freq;

  REAL8 constSftsDataDeltaF = multiSFTsV->data[0]->data[0]->data->deltaF;
  REAL8 constSftsDataF0 = multiSFTsV->data[0]->data[0]->data->f0;

  REAL4 constTsft = (REAL4)(1.0 / constSftsDataDeltaF);                       /* length of SFTs in seconds */
  REAL4 constDFreq = (REAL4)constSftsDataDeltaF;
  INT4 constFreqIndex0 = (UINT4)(constSftsDataF0 / constDFreq + 0.5);         /* index of first frequency-bin in SFTs */

  /* increment the function call counter */
  ++call_count;

  /* report which flavour of function is called */
  if (call_count == 1) {
    if (USE_OPENCL_KERNEL) {
      LogPrintf (LOG_DEBUG, "%s: using OpenCL call on GPUs. ", fn);
    } else if (USE_OPENCL_KERNEL_CPU) {
      LogPrintf (LOG_DEBUG, "%s: using OpenCL kernel as a regular function on CPU. ", fn);
    }
    if (CALL_XLALCoreFstatREAL4) {
      LogPrintf (LOG_DEBUG, "calling function XLALCoreFstatREAL4", fn);
    }
    LogPrintf (LOG_DEBUG, "\n");
  }

  /* check input consistency */
  if ( !doppler || !multiSFTsV || !multiWeightsV || !multiDetStatesV || !cfvBuffer ) {
    XLALPrintError ("%s: illegal NULL input pointer.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( !fstatBandV || !fstatBandV->length ) {
    XLALPrintError ("%s: illegal NULL or empty output pointer 'fstatBandV'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numSegments = fstatBandV->length;

  if ( (multiSFTsV->length != numSegments) || (multiDetStatesV->length != numSegments ) ) {
    XLALPrintError ("%s: inconsistent number of segments between fstatBandV (%d), multiSFTsV(%d) and multiDetStatesV (%d)\n",
                    fn, numSegments, multiSFTsV->length, multiDetStatesV->length );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( multiWeightsV->length != numSegments ) {
    XLALPrintError ("%s: inconsistent number of segments between fstatBandV (%d) and multiWeightsV (%d)\n",
                    fn, numSegments, multiWeightsV->length );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numBins = fstatBandV->data[0].data->length;
  f0      = fstatBandV->data[0].f0;
  deltaF  = fstatBandV->data[0].deltaF;

  /* a check that the f0 values from thisPoint and fstatVector are
   * at least close to each other -- this is only meant to catch
   * stupid errors but not subtle ones */
  if ( fabs(f0 - doppler->fkdot[0]) >= deltaF ) {
    XLALPrintError ("%s: fstatVector->f0 = %f differs from doppler->fkdot[0] = %f by more than deltaF = %g\n", fn, f0, doppler->fkdot[0], deltaF );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* ---------- prepare REAL4 version of PulsarSpins to be passed into core functions */
  PulsarSpinsREAL4 fkdot4 = empty_PulsarSpinsREAL4;
  UINT4 s;
  Freq = Freq0 = doppler->fkdot[0];
  fkdot4.FreqMain = (INT4)Freq0;
  /* fkdot[0] now only carries the *remainder* of Freq wrt FreqMain */
  fkdot4.fkdot[0] = (REAL4)( Freq0 - (REAL8)fkdot4.FreqMain );
  /* the remaining spins are simply down-cast to REAL4 */
  for ( s=1; s < PULSAR_MAX_SPINS; s ++ )
    fkdot4.fkdot[s] = (REAL4) doppler->fkdot[s];
  UINT4 maxs;
  /* find highest non-zero spindown-entry */
  for ( maxs = PULSAR_MAX_SPINS - 1;  maxs > 0 ; maxs --  )
    if ( fkdot4.fkdot[maxs] != 0 )
      break;
  fkdot4.spdnOrder = maxs;

  for (k=1; k<=fkdot4.spdnOrder; k++) clWp->fkdot16.data[k] = fkdot4.fkdot[k];
//  PulsarSpins16 fkdot16;
//  fkdot16.s[0] = fkdot4.spdnOrder;
//  for (k=1; k<=fkdot4.spdnOrder; k++) fkdot16.s[k] = fkdot4.fkdot[k];


  /* ---------- Buffering quantities that don't need to be recomputed ---------- */

  /* check if for this skyposition and data, the SSB+AMcoefs were already buffered */
  if ( cfvBuffer->Alpha != doppler->Alpha ||
       cfvBuffer->Delta != doppler->Delta ||
       cfvBuffer->multiDetStatesV != multiDetStatesV ||
       cfvBuffer->numSegments != numSegments
       )
    { /* no ==> compute and buffer */

      LogPrintf(LOG_DEBUG, "In function %s: buffering quantities that don't need to be recomputed...\n", fn);

      SkyPosition skypos;
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;

      /* prepare buffer */
      XLALEmptyComputeFBufferREAL4V ( cfvBuffer );

      /* store new precomputed quantities in buffer */
      cfvBuffer->Alpha = doppler->Alpha;
      cfvBuffer->Delta = doppler->Delta;
      cfvBuffer->multiDetStatesV = multiDetStatesV ;
      cfvBuffer->numSegments = numSegments;

      if ( (cfvBuffer->multiSSB4V = XLALCalloc ( numSegments, sizeof(*cfvBuffer->multiSSB4V) )) == NULL ) {
        XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", fn, numSegments, sizeof(*cfvBuffer->multiSSB4V) );
        XLAL_ERROR ( fn, XLAL_ENOMEM );
      }
      if ( (cfvBuffer->multiAMcoefV = XLALCalloc ( numSegments, sizeof(*cfvBuffer->multiAMcoefV) )) == NULL ) {
        XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
        XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", fn, numSegments, sizeof(*cfvBuffer->multiAMcoefV) );
        XLAL_ERROR ( fn, XLAL_ENOMEM );
      }

      for ( n=0; n < numSegments; n ++ )
        {
          /* compute new SSB timings over all segments */
          if ( (cfvBuffer->multiSSB4V[n] = XLALGetMultiSSBtimesREAL4 ( multiDetStatesV->data[n], doppler->Alpha, doppler->Delta, doppler->refTime)) == NULL ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError ( "%s: XLALGetMultiSSBtimesREAL4() failed. xlalErrno = %d.\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

          LALStatus status = empty_LALStatus;
          LALGetMultiAMCoeffs ( &status, &(cfvBuffer->multiAMcoefV[n]), multiDetStatesV->data[n], skypos );
          if ( status.statusCode ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError ("%s: LALGetMultiAMCoeffs() failed with statusCode=%d, '%s'\n", fn, status.statusCode, status.statusDescription );
            XLAL_ERROR ( fn, XLAL_EFAILED );
          }

          /* apply noise-weights to Antenna-patterns and compute A,B,C */
          if ( XLALWeighMultiAMCoeffs ( cfvBuffer->multiAMcoefV[n], multiWeightsV->data[n] ) != XLAL_SUCCESS ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError("%s: XLALWeighMultiAMCoeffs() failed with error = %d\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

          /*
           * OpenCL: copy the data to the flat 1D memory buffers
           */
          {
            MultiSSBtimesREAL4 *multiSSB4 = cfvBuffer->multiSSB4V[n];
            MultiAMCoeffs *multiAMcoeff = cfvBuffer->multiAMcoefV[n];
            SSBtimesREAL4 *tSSB;
            AMCoeffs *amcoe;
            UINT4 X, offset, len;

            for (X=0; X<clWp->numIFOs; X++) {
              offset = (n*clWp->numIFOs + X) * clWp->maxNumSFTs;
              tSSB = multiSSB4->data[X];
              amcoe = multiAMcoeff->data[X];
              len = tSSB->DeltaT_int->length * sizeof(REAL4);

              // 1D-index of an array element {n, X, s}:
              // ind = (n*clWp->numIFOs + X) * clWp->maxNumSFTs + s

              memcpy (clWp->tSSB_DeltaT_int.data + offset, tSSB->DeltaT_int->data, len);
              memcpy (clWp->tSSB_DeltaT_rem.data + offset, tSSB->DeltaT_rem->data, len);
              memcpy (clWp->tSSB_TdotM1.data + offset, tSSB->TdotM1->data, len);

              memcpy (clWp->amcoe_a.data + offset, amcoe->a->data, len);
              memcpy (clWp->amcoe_b.data + offset, amcoe->b->data, len);
            }

            clWp->ABCInvD.data[n].Ad = multiAMcoeff->Mmunu.Ad;
            clWp->ABCInvD.data[n].Bd = multiAMcoeff->Mmunu.Bd;
            clWp->ABCInvD.data[n].Cd = multiAMcoeff->Mmunu.Cd;
            clWp->ABCInvD.data[n].InvDd = 1.0f / multiAMcoeff->Mmunu.Dd;

          }

        } /* for n < numSegments */

        /*
         * OpenCL: initialize the array of REAL4-split frequencies
         */
        Freq = Freq0;

        for (k=0; k<clWp->numBins; k++) {
          clWp->Freq.data[k].FreqMain = (INT4)Freq;
          clWp->Freq.data[k].fkdot0 = (REAL4)( Freq - (REAL8)fkdot4.FreqMain );
          Freq += deltaF;
        }

#if USE_OPENCL_KERNEL
        /*
         * OpenCL: copy the buffers to the device memory
         */
        UINT4 l2 = clWp->numSegments * clWp->numIFOs;
        UINT4 l3 = l2 * clWp->maxNumSFTs;
        err_total = CL_SUCCESS;
        err = clEnqueueWriteBuffer ( *(clWp->cmd_queue),           // which queue
                                     clWp->tSSB_DeltaT_int.memobj, // destination device pointer
                                     CL_TRUE,                     // blocking write?
                                     0,                           // offset
                                     l3 * sizeof(REAL4),          // size in bytes
                                     clWp->tSSB_DeltaT_int.data,   // source pointer
                                     0,                           // cl_uint num_events_in_wait_list
                                     NULL,                        // const cl_event *event_wait_list
                                     NULL);                       // cl_event *event
        err_total += (err-CL_SUCCESS);

        err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->tSSB_DeltaT_rem.memobj,
                                     CL_TRUE, 0, l3 * sizeof(REAL4),
                                     clWp->tSSB_DeltaT_rem.data, 0, NULL, NULL);
        err_total += (err-CL_SUCCESS);

        err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->tSSB_TdotM1.memobj,
                                     CL_TRUE, 0, l3 * sizeof(REAL4),
                                     clWp->tSSB_TdotM1.data, 0, NULL, NULL);
        err_total += (err-CL_SUCCESS);

        err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->amcoe_a.memobj,
                                     CL_TRUE, 0, l3 * sizeof(REAL4),
                                     clWp->amcoe_a.data, 0, NULL, NULL);
        err_total += (err-CL_SUCCESS);

        err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->amcoe_b.memobj,
                                     CL_TRUE, 0, l3 * sizeof(REAL4),
                                     clWp->amcoe_b.data, 0, NULL, NULL);
        err_total += (err-CL_SUCCESS);

        err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->ABCInvD.memobj,
                                     CL_TRUE, 0, l2 * sizeof(REAL44),
                                     clWp->ABCInvD.data, 0, NULL, NULL);
        err_total += (err-CL_SUCCESS);

        if (err_total != CL_SUCCESS) {
          XLALPrintError ("%s: Error copying data to memory buffer, error code = %d\n", fn, err );
          XLALDestroyCLWorkspace (clWp, multiSFTsV);
          XLAL_ERROR ( fn, XLAL_EINVAL );
        }
#endif // #if USE_OPENCL_KERNEL

    } /* if we could NOT reuse previously buffered quantites */

#if USE_OPENCL_KERNEL
  /*
   * copy frequency arrays
   */
  err_total = CL_SUCCESS;
  err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->Freq.memobj,
                               CL_TRUE, 0, clWp->numBins * sizeof(REAL42),
                               clWp->Freq.data, 0, NULL, NULL);
  err_total += (err-CL_SUCCESS);

  err = clEnqueueWriteBuffer ( *(clWp->cmd_queue), clWp->fkdot16.memobj,
                               CL_TRUE, 0, 16 * sizeof(REAL4),
                               clWp->fkdot16.data, 0, NULL, NULL);
  err_total += (err-CL_SUCCESS);
  if (err_total != CL_SUCCESS) {
    XLALPrintError ("%s: Error copying frequency data to device memory, error code = %d\n", fn, err );
    XLALDestroyCLWorkspace (clWp, multiSFTsV);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }


  /*
   * OpenCL: set kernel arguments
   */
  err_total = CL_SUCCESS;
  err = clSetKernelArg(*(clWp->kernel),       // wchich kernel
                       0,                    // argument index
                       sizeof(cl_mem),       // argument data size
                       (void *)&(clWp->Fstat.memobj) ); // argument data
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 1, sizeof(cl_mem), (void *)&(clWp->multiSFTsFlat.memobj));
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 2, sizeof(cl_mem), (void *)&(clWp->numSFTsV.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 3, sizeof(UINT4), (void *)&(clWp->sftLen) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 4, sizeof(REAL4), (void *)&constTsft);
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 5, sizeof(REAL4), (void *)&constDFreq);
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 6, sizeof(INT4), (void *)&constFreqIndex0);
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 7, sizeof(cl_mem), (void *)&(clWp->Freq.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 8, sizeof(UINT4), (void *)&(fkdot4.spdnOrder) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 9, sizeof(cl_mem), (void *)&(clWp->fkdot16.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 10, sizeof(cl_mem), (void *)&(clWp->tSSB_DeltaT_int.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 11, sizeof(cl_mem), (void *)&(clWp->tSSB_DeltaT_rem.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 12, sizeof(cl_mem), (void *)&(clWp->tSSB_TdotM1.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 13, sizeof(cl_mem), (void *)&(clWp->amcoe_a.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 14, sizeof(cl_mem), (void *)&(clWp->amcoe_b.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 15, sizeof(cl_mem), (void *)&(clWp->ABCInvD.memobj) );
  err_total += (err-CL_SUCCESS);

  err = clSetKernelArg(*(clWp->kernel), 16, sizeof(FcomponentsREAL4)*clWp->numIFOs*clWp->maxNumSFTs, NULL );
  err_total += (err-CL_SUCCESS);

  if (err_total != CL_SUCCESS) {
    XLALPrintError ("%s: Error while setting the kernel arguments, error code = %d\n", fn, err );
    XLALDestroyCLWorkspace (clWp, multiSFTsV);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }


  /*
   * OpenCL: enqueue kernel for execution
   * block-thread geometry: (numSegments,numBins,1) x (maxNumSFTs, numIFOs(2), 1)
   */
  LogPrintf(LOG_DEBUG, "In function %s: launching the kernel...\n", fn);

  size_t local_work_size[2], global_work_size[2];
  local_work_size[0] = clWp->maxNumSFTs;
  local_work_size[1] = clWp->numIFOs;
  global_work_size[0] = local_work_size[0] * clWp->numSegments;
  global_work_size[1] = local_work_size[1] * numBins;

  err = clEnqueueNDRangeKernel(*(clWp->cmd_queue), *(clWp->kernel),
                               2, // Work dimensions
                               NULL, // must be NULL (work offset)
                               global_work_size,  // global work items grid dimension
                               local_work_size,   // local work group size
                               0, // no events to wait on
                               NULL, // event list
                               NULL); // event for this kernel
  if (err != CL_SUCCESS) {
    XLALPrintError ("%s: Error enqueueing the kernel, error code = %d\n", fn, err );
    XLALDestroyCLWorkspace (clWp, multiSFTsV);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /*
   * Read output memory buffer after the kernel call
   */
  err = clEnqueueReadBuffer( *(clWp->cmd_queue),  // which queue
                             clWp->Fstat.memobj,  // source device memory buffer
                             CL_TRUE,            // blocking read
                             0,                  // offset
                             sizeof(REAL4)*clWp->Fstat.length, // size
                             clWp->Fstat.data,    // pointer
                             0, NULL, NULL);     // events
  if (err != CL_SUCCESS) {
    XLALPrintError ("%s: Error reading output buffer, error code = %d\n", fn, err );
    XLALDestroyCLWorkspace (clWp, multiSFTsV);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /*
   * Store the results in fstatBandV
   */
  for ( n = 0; n < numSegments; n ++ ) {
    for ( k = 0; k < numBins; k++) {
      fstatBandV->data[n].data->data[k] = clWp->Fstat.data[k * numSegments + n];
    }
  }


#endif // #if USE_OPENCL_KERNEL

  /* loop over all segments and compute FstatVector over frequencies for each */

#if USE_OPENCL_KERNEL_CPU
  FcomponentsREAL4 *FaFb_components;
  if ( (FaFb_components = XLALCalloc ( clWp->numIFOs * clWp->maxNumSFTs, sizeof(FcomponentsREAL4) )) == NULL ) {
    XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
    XLALDestroyCLWorkspace (clWp, multiSFTsV);
    XLALPrintError ( "%s: memory allocation for FaFb_components failed\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
#endif

  for ( n = 0; n < numSegments; n ++ )
    {

      Freq = Freq0;	/* reset frequency to start value for next segment */

      /* loop over frequency values and fill up values in fstatVector */
      for ( k = 0; k < numBins; k++)
        {
          /* REAL4-split frequency value */
          fkdot4.FreqMain = (INT4)Freq;
          fkdot4.fkdot[0] = (REAL4)( Freq - (REAL8)fkdot4.FreqMain );

          /* call the core function to compute one multi-IFO F-statistic */

          if (CALL_XLALCoreFstatREAL4) {
	      /* >>>>> this function could run on the GPU device <<<<< */
	      XLALCoreFstatREAL4 ( &Fstat, &fkdot4, multiSFTsV->data[n], cfvBuffer->multiSSB4V[n], cfvBuffer->multiAMcoefV[n], Dterms );

	      if ( xlalErrno ) {
		XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
		XLALPrintError ("%s: XLALCoreFstatREAL4() failed with errno = %d in loop n=%d, k=%d.\n", fn, xlalErrno, n, k );
		XLAL_ERROR ( fn, XLAL_EFUNC );
	      }

	      fstatBandV->data[n].data->data[k] = Fstat;

	      /* increase frequency by deltaF */
	      Freq += deltaF;
          }

#if USE_OPENCL_KERNEL_CPU
          UINT4 X, alpha;
          for (X = 0; X < clWp->numIFOs; X++ ) {
            for (alpha = 0; alpha < clWp->maxNumSFTs; alpha++ ) {
              OpenCLComputeFstatFaFb(clWp->Fstat.data,
                                     n, // curSegment
                                     k, // curBin
                                     clWp->maxNumSFTs,
                                     alpha,
                                     X,
                                     numSegments,
                                     clWp->multiSFTsFlat.data,
                                     clWp->numSFTsV.data,
                                     clWp->sftLen,
                                     constTsft,
                                     constDFreq,
                                     constFreqIndex0,
                                     clWp->Freq.data,
                                     fkdot4.spdnOrder,
                                     clWp->fkdot16.data,
                                     clWp->tSSB_DeltaT_int.data,
                                     clWp->tSSB_DeltaT_rem.data,
                                     clWp->tSSB_TdotM1.data,
                                     clWp->amcoe_a.data,
                                     clWp->amcoe_b.data,
                                     clWp->ABCInvD.data,
                                     FaFb_components);

              if (alpha) {
                FaFb_components[X * clWp->maxNumSFTs].Fa.re += FaFb_components[X * clWp->maxNumSFTs + alpha].Fa.re;
                FaFb_components[X * clWp->maxNumSFTs].Fa.im += FaFb_components[X * clWp->maxNumSFTs + alpha].Fa.im;
                FaFb_components[X * clWp->maxNumSFTs].Fb.re += FaFb_components[X * clWp->maxNumSFTs + alpha].Fb.re;
                FaFb_components[X * clWp->maxNumSFTs].Fb.im += FaFb_components[X * clWp->maxNumSFTs + alpha].Fb.im;
              }

            }
          }


          REAL4 Fa_re = (FaFb_components[0].Fa.re + FaFb_components[clWp->maxNumSFTs].Fa.re) * OOTWOPI_FLOAT;
          REAL4 Fa_im = (FaFb_components[0].Fa.im + FaFb_components[clWp->maxNumSFTs].Fa.im) * OOTWOPI_FLOAT;
          REAL4 Fb_re = (FaFb_components[0].Fb.re + FaFb_components[clWp->maxNumSFTs].Fb.re) * OOTWOPI_FLOAT;
          REAL4 Fb_im = (FaFb_components[0].Fb.im + FaFb_components[clWp->maxNumSFTs].Fb.im) * OOTWOPI_FLOAT;

          REAL4 Ad =  clWp->ABCInvD.data[n].Ad;
          REAL4 Bd =  clWp->ABCInvD.data[n].Bd;
          REAL4 Cd =  clWp->ABCInvD.data[n].Cd;
          REAL4 Dd_inv =  clWp->ABCInvD.data[n].InvDd;


          clWp->Fstat.data[k * numSegments + n] = Dd_inv * ( Bd * (SQ(Fa_re) + SQ(Fa_im) )
                                               + Ad * ( SQ(Fb_re) + SQ(Fb_im) )
                                               - 2.0f * Cd *( Fa_re * Fb_re + Fa_im * Fb_im )
                                               );
          fstatBandV->data[n].data->data[k] = clWp->Fstat.data[k * numSegments + n];
#endif // #if USE_OPENCL_KERNEL_CPU

        } /* for k < numBins */

    } /* for n < numSegments */



#if USE_OPENCL_KERNEL_CPU
  LALFree (FaFb_components);
#endif

  return XLAL_SUCCESS;

} /* XLALComputeFStatFreqBandVector() */


/** Initialize OpenCL workspace
 * Create memory objects associated with OpenCL context
 * and memory buffers */
int
XLALInitCLWorkspace ( CLWorkspace *clW,
                      const MultiSFTVectorSequence *stackMultiSFT )
{
  static const char *fn = "XLALInitCLWorkspace()";
  static const char *cl_kernel_filepath = "/Users/oleg/lalsuite/lalapps/src/pulsar/FDS_isolated/kernel.cl"; //TODO: do something with hardcoded kernel path

#if USE_OPENCL_KERNEL
  cl_int err, err_total;
  const cl_uint max_num_platforms = 3;
  static cl_platform_id platforms[3];
  cl_uint num_platforms;
  char strInfo[100];
  static cl_context context;
  const cl_uint max_num_devices = 4;
  cl_device_id devices[4];
  cl_uint num_devices;
  static cl_command_queue cmd_queue;
  static cl_program program;
  static cl_kernel kernel;

  clW->platform  = NULL;
  clW->device    = NULL;
  clW->context   = NULL;
  clW->cmd_queue = NULL;
  clW->program   = NULL;
  clW->kernel    = NULL;
#endif

  clW->multiSFTsFlat.data = NULL;
  clW->numSFTsV.data = NULL;
  clW->tSSB_DeltaT_int.data = NULL;
  clW->tSSB_DeltaT_rem.data = NULL;
  clW->tSSB_TdotM1.data = NULL;
  clW->amcoe_a.data = NULL;
  clW->amcoe_b.data = NULL;
  clW->ABCInvD.data = NULL;
  clW->Fstat.data = NULL;
  clW->fkdot16.data = NULL;

  LogPrintf(LOG_DEBUG, "In function %s: initializing OpenCL workspace\n", fn);

#if USE_OPENCL_KERNEL
  // query the platform ID
  LogPrintf(LOG_DEBUG, "In function %s: query the platform ID\n", fn);
  clGetPlatformIDs(max_num_platforms, platforms, &num_platforms);
  clW->platform = &(platforms[gpu_platform_id]);

  LogPrintf(LOG_DEBUG, "In function %s: Found %d platforms, using platform id %d\n", fn, num_platforms, gpu_platform_id);

  // query OpenCL platform info
  LogPrintf(LOG_DEBUG, "In function %s: query the OpenCL platform info\n", fn);
  err = clGetPlatformInfo ( *(clW->platform), CL_PLATFORM_PROFILE, 100, strInfo, NULL );
  if (err != CL_SUCCESS) {
      XLALPrintError ("%s: Error calling clGetPlatformInfo.\n", fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  // create OpenCL GPU context
  LogPrintf(LOG_DEBUG, "In function %s: create the OpenCL GPU context\n", fn);
  context = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);
  if (context == (cl_context)0) {
      XLALPrintError ("%s: Failed to create context\n", fn );
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  clW->context = &context;

  // get the list of available GPU devices
  LogPrintf(LOG_DEBUG, "In function %s: get the list of all available GPU devices\n", fn);
  err = clGetDeviceIDs( *(clW->platform),
                    CL_DEVICE_TYPE_GPU,
                    max_num_devices,
                    devices,
                    &num_devices);
  if (err != CL_SUCCESS) {
      XLALPrintError ("%s: Error querying number of OpenCL devices\n", fn );
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
  } else if ( gpu_device_id >= num_devices ) {
      XLALPrintError ("%s: Invalid device id %d, found only %d devices\n", fn, num_devices, gpu_device_id);
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
  } else {
    LogPrintf(LOG_DEBUG, "In function %s: Found %d devices, using device id %d\n", fn, num_devices, gpu_device_id);
  }
  clW->device = &(devices[gpu_device_id]);

  // create a command-queue
  LogPrintf(LOG_DEBUG, "In function %s: create OpenCL command queue\n", fn);
  cmd_queue = clCreateCommandQueue(*(clW->context), *(clW->device),
                                   CL_QUEUE_PROFILING_ENABLE,
                                   &err);
  if (cmd_queue == (cl_command_queue)0) {
      XLALPrintError ("%s: Failed to create command queue\n", fn );
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  clW->cmd_queue = &cmd_queue;

#endif // #if USE_OPENCL_KERNEL

  // allocate flattened memory buffers on host

  // WARNING: HARDCODED VALUES FOR NOW (RUN S5R4)
  // TODO: include proper assertions here to ensure correctness
  clW->numSegments = stackMultiSFT->length;
  clW->numIFOs = NUM_IFOS;
  clW->maxNumSFTs = MAX_NUM_SFTS;
  clW->sftLen = stackMultiSFT->data[0]->data[0]->data[0].data->length;
  clW->numBins = 200;

  // allocate large buffer arrays on a host
  LogPrintf(LOG_DEBUG, "In function %s: allocate 1D buffer arrays\n", fn);
  UINT4 l2 = clW->numSegments * clW->numIFOs;
  UINT4 l3 = l2 * clW->maxNumSFTs;
  UINT4 l4 = l3 * clW->sftLen;

  clW->multiSFTsFlat.length = l4;
  clW->multiSFTsFlat.data = XLALMalloc( sizeof(COMPLEX8) * l4);

  clW->numSFTsV.length = l2;
  clW->numSFTsV.data = XLALMalloc( sizeof(UINT4) * l2 );

  clW->tSSB_DeltaT_int.length = l3;
  clW->tSSB_DeltaT_int.data = XLALMalloc( sizeof(REAL4) * l3 );

  clW->tSSB_DeltaT_rem.length = l3;
  clW->tSSB_DeltaT_rem.data = XLALMalloc( sizeof(REAL4) * l3 );

  clW->tSSB_TdotM1.length = l3;
  clW->tSSB_TdotM1.data = XLALMalloc( sizeof(REAL4) * l3 );

  clW->amcoe_a.length = l3;
  clW->amcoe_a.data = XLALMalloc( sizeof(REAL4) * l3 );

  clW->amcoe_b.length = l3;
  clW->amcoe_b.data = XLALMalloc( sizeof(REAL4) * l3 );

  clW->ABCInvD.length = l2;
  clW->ABCInvD.data = XLALMalloc( sizeof(REAL44) * l2 );

  clW->fkdot16.length = 16; //hardcoded value!
  clW->fkdot16.data = XLALMalloc (sizeof(REAL4) * clW->fkdot16.length);

  // clW->fkdot16: initialized in RearrangeSFTData
  // clW->Freq: initialized in RearrangeSFTData
  // clW->Fstat:  initialized in RearrangeSFTData

  if ( clW->multiSFTsFlat.data == NULL || clW->numSFTsV.data == NULL
       || clW->tSSB_DeltaT_int.data == NULL || clW->tSSB_DeltaT_rem.data == NULL || clW->tSSB_DeltaT_rem.data == NULL
       || clW->amcoe_a.data == NULL || clW->amcoe_b.data == NULL ) {
      XLALPrintError ("%s: XLALMalloc() failed.\n", fn );
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  { // SFT data rearrangement block
    UINT4 n, X, s;

    MultiSFTVector *multiSFT;
    SFTVector *sFT;
    COMPLEX8Vector *cV;
    COMPLEX8 *ptrData;
    //TODO: state your assertions here!!!
    /* Traverse stackMultiSFT and do the following:
     * 1. copy data segments sequentially to multiSFTsFlat;
     * 2. deallocate data segments in stackMultiSFT;
     * 3. repoint the data pointers to the middle of multiSFTsFlat.
     *
     * Index of a data element at position m in SFT s, detector X and segment n, is:
     * ind = m + cV->sftLen * (s + cV->maxNumSFTs * (X + cV->numIFOs * n))
     */
    LogPrintf(LOG_DEBUG, "In function %s: flatten the stackMultiSFT data structure\n", fn);
    ptrData = clW->multiSFTsFlat.data;

    if (clW->numSegments != stackMultiSFT->length) {
      XLALPrintError ("%s: internal error: inconsistent clW->numSegments\n", fn);
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }
    for (n=0; n<stackMultiSFT->length; n++) {
      multiSFT = stackMultiSFT->data[n];
      if (clW->numIFOs != multiSFT->length) {
        XLALPrintError ("%s: internal error: inconsistent clW->numIFOs for segment %d\n", fn, n);
        XLALDestroyCLWorkspace (clW, stackMultiSFT);
        XLAL_ERROR ( fn, XLAL_EINVAL );
      }
      for (X=0; X<multiSFT->length; X++) {
        sFT = multiSFT->data[X];
        if (clW->maxNumSFTs < sFT->length) {
          XLALPrintError ("%s: internal error: number of SFTs exceeds MAX_NUM_SFTS for segment %d, detector %d\n", fn, n, X);
          XLALDestroyCLWorkspace (clW, stackMultiSFT);
          XLAL_ERROR ( fn, XLAL_EINVAL );
        }
        clW->numSFTsV.data[n * multiSFT->length + X] = sFT->length;
        for (s=0; s<sFT->length; s++) {
          cV = sFT->data[s].data;
          if (clW->sftLen != cV->length) {
            XLALPrintError ("%s: internal error: inconsistent SFT length in segment=%d, detector=%d, SFT %d\n", fn, n, X, s);
            XLALDestroyCLWorkspace (clW, stackMultiSFT);
            XLAL_ERROR ( fn, XLAL_EINVAL );
          }
          memcpy (ptrData, cV->data, cV->length * sizeof(COMPLEX8));

          LALFree (cV->data);
          cV->data = ptrData;
          ptrData += clW->sftLen;
        }
        ptrData += (clW->maxNumSFTs - sFT->length) * clW->sftLen;
      }
    }
  } // end of SFT data rearrangement block

#if USE_OPENCL_KERNEL
  // allocate buffer arrays on the device
  // only SFT array is copied to the device; the rest of the arrays are created empty and will be filled later
  LogPrintf(LOG_DEBUG, "In function %s: allocate OpenCL device memory buffers\n", fn);
  err_total = CL_SUCCESS;
  clW->multiSFTsFlat.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(COMPLEX8) * clW->multiSFTsFlat.length, clW->multiSFTsFlat.data, &err);
  err_total += (err-CL_SUCCESS);

  clW->numSFTsV.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(UINT4) * clW->numSFTsV.length, clW->numSFTsV.data, &err);
  err_total += (err-CL_SUCCESS);

  clW->tSSB_DeltaT_int.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL4) * clW->tSSB_DeltaT_int.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  clW->tSSB_DeltaT_rem.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL4) * clW->tSSB_DeltaT_rem.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  clW->tSSB_TdotM1.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL4) * clW->tSSB_TdotM1.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  clW->amcoe_a.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL4) * clW->amcoe_a.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  clW->amcoe_b.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL4) * clW->amcoe_b.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  clW->ABCInvD.memobj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(REAL44) * clW->ABCInvD.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  clW->fkdot16.memobj = clCreateBuffer (*(clW->context), CL_MEM_READ_ONLY, sizeof(REAL4)*clW->fkdot16.length, NULL, &err);
  err_total += (err-CL_SUCCESS);

  if (err_total != CL_SUCCESS) {
      XLALPrintError ("%s: Error creating OpenCL memory buffer, error code = %d\n", fn, err );
      XLALDestroyCLWorkspace (clW, stackMultiSFT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }


#ifdef OPENCL_KERNEL_TEXT

  extern char*opencl_kernel_text;
#define cl_kernel_strings opencl_kernel_text

#else

  // read kernel source into memory
  LogPrintf(LOG_DEBUG, "In function %s: read kernel source into memory\n", fn);

  FILE *fd;
  char *cl_kernel_strings;
  const size_t max_kernel_bytes = 50000;
  size_t bytes_read;

  if ( ( cl_kernel_strings  = XLALMalloc( max_kernel_bytes * sizeof(char))) == NULL ) {
      XLALPrintError ("%s: XLALMalloc() failed.\n", fn );
      XLALDestroyCLWorkspace (clW, stackMultiSFT );
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if((fd=fopen(cl_kernel_filepath, "rb"))==NULL) {
    XLALPrintError ("%s: ERROR: Cannot open OpenCL kernel file at location \"%s\".\n", fn, cl_kernel_filepath );
    XLALDestroyCLWorkspace (clW, stackMultiSFT );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  fread (cl_kernel_strings, max_kernel_bytes, 1, fd);
  if (! feof(fd) || ferror(fd)) {
    XLALPrintError ("%s: ERROR: Cannot read OpenCL kernel file at location \"%s\".\n", fn, cl_kernel_filepath );
    XLALDestroyCLWorkspace (clW, stackMultiSFT );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  bytes_read = ftell(fd);
  fclose (fd);
  cl_kernel_strings[bytes_read] = '\0'; // null-terminated string

#endif

  // create the program
  LogPrintf(LOG_DEBUG, "In function %s: create OpenCL program\n", fn);
  program = clCreateProgramWithSource( *(clW->context),
                                       1, // string count
                                       (const char **) &cl_kernel_strings, // program strings
                                       NULL, // string lengths
                                       &err); // error code

#ifndef OPENCL_KERNEL_TEXT
  LALFree(cl_kernel_strings);
#endif

  if (program == (cl_program)0) {
    XLALPrintError( "%s: ERROR: failed to create OpenCL program\n", fn);
    XLALDestroyCLWorkspace (clW, stackMultiSFT );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  clW->program = &program;

  // build the program
  LogPrintf(LOG_DEBUG, "In function %s: build OpenCL program...\n", fn);
  err = clBuildProgram(*(clW->program),
                       0,     // num devices in device list
                       NULL,  // device list
                       NULL,  // options
                       NULL,  // notifier callback function ptr
                       NULL); // error code
  if (err != CL_SUCCESS) {
    size_t len;
    char debug_buffer[2048];
    XLALPrintError( "%s: ERROR: failed to compile OpenCL program\n", fn);
    clGetProgramBuildInfo(*(clW->program), *(clW->device), CL_PROGRAM_BUILD_LOG,
                          sizeof(debug_buffer), debug_buffer, &len);
    XLALPrintError("%s\n", debug_buffer);
    XLALDestroyCLWorkspace (clW, stackMultiSFT );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  // finally, create the kernel
  LogPrintf(LOG_DEBUG, "In function %s: create kernel...\n", fn);
  kernel = clCreateKernel(*(clW->program), "OpenCLComputeFstatFaFb", NULL);
  if (kernel == (cl_kernel)0) {
    XLALPrintError( "%s: ERROR: failed to create kernel\n", fn);
    XLALDestroyCLWorkspace (clW, stackMultiSFT );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  clW->kernel = &kernel;
#endif // #if USE_OPENCL_KERNEL

  return 0;
} /* XLALInitCLWorkspace() */



/** Rearrange SFT data structures
 * Flatten the SFT data: combine small chunks of memory into a single
 * contiguous array, accessable via 4d-index */
void
XLALRearrangeSFTData ( CLWorkspace *clW,
                       const REAL4FrequencySeriesVector *fstatBandV )
{
  static const char *fn = "XLALRearrangeSFTData()";
  static int call_count = 0;
  ++call_count;

  LogPrintf(LOG_DEBUG, "In function %s: rearrange SFT data structures\n", fn);

  clW->numBins = fstatBandV->data[0].data->length;

  // deallocate previously allocated memory
  if (clW->Freq.data) LALFree(clW->Freq.data);
  if (clW->Fstat.data)  LALFree(clW->Fstat.data);

  // allocate memory for new arrays
  clW->Freq.length = clW->numBins;
  clW->Freq.data = XLALMalloc( sizeof(REAL42) * clW->numBins );

  clW->Fstat.length = clW->numSegments * clW->numBins;
  clW->Fstat.data = XLALMalloc( sizeof(REAL4) * clW->Fstat.length );

  if ( clW->Fstat.data == NULL || clW->Freq.data == NULL ) {
      XLALPrintError ("%s: XLALMalloc() failed.\n", fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
  }

#if USE_OPENCL_KERNEL
  { // create memory buffers on the device
    cl_int err, err_total = CL_SUCCESS;

    if (call_count > 1) {
      freeCLMemoryObject(&(clW->Freq.memobj));
      freeCLMemoryObject(&(clW->Fstat.memobj));
    }

    clW->Fstat.memobj = clCreateBuffer (*(clW->context), CL_MEM_READ_WRITE, sizeof(REAL4)*clW->Fstat.length, NULL, &err);
    err_total += (err-CL_SUCCESS);

    clW->Freq.memobj = clCreateBuffer(*(clW->context), CL_MEM_READ_ONLY, sizeof(REAL42)*clW->Freq.length, NULL, &err);
    err_total += (err-CL_SUCCESS);

    if (err_total != CL_SUCCESS) {
        XLALPrintError ("%s: Error creating OpenCL memory buffer, error code = %d\n", fn, err_total );
        XLAL_ERROR ( fn, XLAL_EINVAL );
    }
  }
#endif // #if USE_OPENCL_KERNEL

} /* XLALRearrangeSFTData */



/** Close OpenCL workspace
 * Free all objects and memory associated with the OpenCL Workspace */
void
XLALDestroyCLWorkspace ( CLWorkspace *clW,
                         const MultiSFTVectorSequence *stackMultiSFT )
{
  static const char *fn = "XLALDestroyCLWorkspace()";

  LogPrintf(LOG_DEBUG, "In function %s: deallocate memory, release OpenCL context\n", fn);

#if USE_OPENCL_KERNEL
  freeCLMemoryObject(&(clW->multiSFTsFlat.memobj));
  freeCLMemoryObject(&(clW->numSFTsV.memobj));
  freeCLMemoryObject(&(clW->Freq.memobj));
  freeCLMemoryObject(&(clW->tSSB_DeltaT_int.memobj));
  freeCLMemoryObject(&(clW->tSSB_DeltaT_rem.memobj));
  freeCLMemoryObject(&(clW->tSSB_TdotM1.memobj));
  freeCLMemoryObject(&(clW->amcoe_a.memobj));
  freeCLMemoryObject(&(clW->amcoe_b.memobj));
  freeCLMemoryObject(&(clW->ABCInvD.memobj));
  freeCLMemoryObject(&(clW->Fstat.memobj));
  freeCLMemoryObject(&(clW->fkdot16.memobj));
#endif

  if (clW->multiSFTsFlat.data) {
    // deallocate the array which contains flattened data for stackMultiSFTs
    LALFree (clW->multiSFTsFlat.data);

    UINT4 n, X, s;

    MultiSFTVector *multiSFT;
    SFTVector *sFT;
    COMPLEX8Vector *cV;

    // set all data pointers on the lowest level to NULL, since their memory has been
    // already deallocated
    for (n=0; n<stackMultiSFT->length; n++) {
      multiSFT = stackMultiSFT->data[n];
      for (X=0; X<multiSFT->length; X++) {
        sFT = multiSFT->data[X];
        for (s=0; s<sFT->length; s++) {
          cV = sFT->data[s].data;
          cV->data = NULL;
        }
      }
    }
  } // if clW->multiSFTsFlat.data != NULL

  if (clW->tSSB_DeltaT_int.data) LALFree(clW->tSSB_DeltaT_int.data);
  if (clW->tSSB_DeltaT_rem.data) LALFree(clW->tSSB_DeltaT_rem.data);
  if (clW->tSSB_TdotM1.data)     LALFree(clW->tSSB_TdotM1.data);
  if (clW->amcoe_a.data)         LALFree(clW->amcoe_a.data);
  if (clW->amcoe_b.data)         LALFree(clW->amcoe_b.data);
  if (clW->ABCInvD.data)         LALFree(clW->ABCInvD.data);
  if (clW->Freq.data)            LALFree(clW->Freq.data);
  if (clW->numSFTsV.data)        LALFree(clW->numSFTsV.data);
  if (clW->Fstat.data)           LALFree(clW->Fstat.data);
  if (clW->fkdot16.data)         LALFree(clW->fkdot16.data);

#if USE_OPENCL_KERNEL
  if (clW->kernel)    clReleaseKernel(*(clW->kernel));
  if (clW->program)   clReleaseProgram(*(clW->program));
  if (clW->cmd_queue) clReleaseCommandQueue (*(clW->cmd_queue));
  if (clW->context)   clReleaseContext (*(clW->context));
#endif

  return;
} /* XLALDestroyCLWorkspace() */


#if USE_OPENCL_KERNEL
/** A helper function to release OpenCL memory objects */
void freeCLMemoryObject (cl_mem *memobj) {
  cl_uint ref_count, i;

  // get an object's reference count
  clGetMemObjectInfo (*memobj, CL_MEM_REFERENCE_COUNT, sizeof(ref_count), &ref_count, NULL);

  // decrement reference count of a memory object unless its destroyed
  for (i=0;i<ref_count;i++) clReleaseMemObject(*memobj);
}
#endif
