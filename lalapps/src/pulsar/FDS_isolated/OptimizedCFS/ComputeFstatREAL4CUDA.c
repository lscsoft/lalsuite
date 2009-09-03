/*
 * Copyright (C) 2009 Anton Obukhov, Bernd Machenschalk
 * Copyright (C) 2009 NVIDIA Corporation ("NVIDIA")
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

/** \author Anton Obukhov
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * CUDA specific implementation of F-stat computation
 *
 */

/*---------- INCLUDES ----------*/
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/ComputeFstat.h>
#include <lal/LogPrintf.h>

#include "ComputeFstatREAL4.h"

#if defined(__GNUC__) && defined(WIN32)
#include "cuda_win32_gcc_host_defines.h"
#endif

#include <cuda_runtime_api.h>

/*---------- local DEFINES ----------*/

/*----- Macros ----- */

/*---------- internal types ----------*/
typedef struct {
  REAL4 fkdot16[PULSAR_MAX_SPINS-1];     /**< remaining spin-parameters, excluding *fractional* part of Freq = fkdot[0] */
} PulsarSpinsExREAL4;

/*---------- Global variables ----------*/
int cuda_device_id = -1;
static const LALStatus empty_LALStatus;


/*---------- internal prototypes ----------*/
extern void HostWrapperCUDAComputeFstatFaFb    (REAL4 *Fstat,
                                                UINT4 Fstat_pitch,

                                                UINT4 sfts_data_data_length,
                                                COMPLEX8 *ifo0_sfts_data_data_data,
                                                COMPLEX8 *ifo1_sfts_data_data_data,
                                                UINT4 sfts_data_data_data_pitch,
                                                REAL4 Tsft,
                                                REAL4 dFreq,
                                                INT4 freqIndex0,

                                                REAL4 *fkdot4_FreqMain,
                                                REAL4 *fkdot4_fkdot0,
                                                PulsarSpinsExREAL4 fkdot4ex,

                                                REAL4 *ifo0_tSSB_DeltaT_int,
                                                REAL4 *ifo1_tSSB_DeltaT_int,
                                                REAL4 *ifo0_tSSB_DeltaT_rem,
                                                REAL4 *ifo1_tSSB_DeltaT_rem,
                                                REAL4 *ifo0_tSSB_TdotM1,
                                                REAL4 *ifo1_tSSB_TdotM1,

                                                REAL4 *ifo0_amcoe_a,
                                                REAL4 *ifo1_amcoe_a,
                                                REAL4 *ifo0_amcoe_b,
                                                REAL4 *ifo1_amcoe_b,

                                                UINT4 *ifo0_sfts_length,
                                                UINT4 *ifo1_sfts_length,

                                                REAL4 *Ad,
                                                REAL4 *Bd,
                                                REAL4 *Cd,
                                                REAL4 *InvDd,

                                                UINT4 Dterms,

                                                UINT4 numBins,
                                                UINT4 numSegments);



int
XLALComputeFStatFreqBandVectorCUDA (REAL4FrequencySeriesVector *fstatBandV,         /**< [out] Vector of Fstat frequency-bands */
                                    const PulsarDopplerParams *doppler,             /**< parameter-space point to compute F for */
                                    const MultiSFTVectorSequence *multiSFTsV,       /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                                    const MultiNoiseWeightsSequence *multiWeightsV, /**< noise-weights of all SFTs */
                                    const MultiDetectorStateSeriesSequence *multiDetStatesV, /**< 'trajectories' of the different IFOs */
                                    UINT4 Dterms,                                   /**< number of Dirichlet kernel Dterms to use */
                                    ComputeFBufferREAL4V *cfvBuffer                 /**< buffer quantities that don't need to be recomputed */
                                    )
{
    static const char *fn = "XLALComputeFStatFreqBandVector()";

    UINT4 numBins, k;
    UINT4 numSegments, n;
    REAL8 f0, deltaF;
    REAL8 Freq0, Freq;
    PulsarSpinsExREAL4 fkdot4ex;
    UINT4 s;
    UINT4 m;

    //////////////////////////////////////////////////////////////////////////
    //CUDA related

    cudaError_t res;

    UINT4 constIFOs = multiSFTsV->data[0]->length;
    UINT4 constSftsDataDataLength = multiSFTsV->data[0]->data[0]->data->data->length;
    UINT4 maxNumSFTs = 0;

    REAL8 constSftsDataDeltaF = multiSFTsV->data[0]->data[0]->data->deltaF;
    REAL8 constSftsDataF0 = multiSFTsV->data[0]->data[0]->data->f0;

    REAL4 constTsft = (REAL4)(1.0 / constSftsDataDeltaF);                       /* length of SFTs in seconds */
    REAL4 constDFreq = (REAL4)constSftsDataDeltaF;
    INT4 constFreqIndex0 = (UINT4)(constSftsDataF0 / constDFreq + 0.5);         /* index of first frequency-bin in SFTs */

    struct cudaDeviceProp curDevProps;


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

    /* ---------- prepare REAL4 version of PulsarSpins (without FreqMain and fkdot[0] to be passed into core functions */
    Freq = Freq0 = doppler->fkdot[0];
    /* the remaining spins are simply down-cast to REAL4 */
    for (s=0; s < PULSAR_MAX_SPINS - 1; s++)
    {
        fkdot4ex.fkdot16[s] = (REAL4)doppler->fkdot[s+1];
    }

    /* ---------- Buffering quantities that don't need to be recomputed ---------- */

    /* check if for this skyposition and data, the SSB+AMcoefs were already buffered */
    if ( cfvBuffer->Alpha != doppler->Alpha ||
        cfvBuffer->Delta != doppler->Delta ||
        cfvBuffer->multiDetStatesV != multiDetStatesV ||
        cfvBuffer->numSegments != numSegments
        )
    { /* no ==> compute and buffer */

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
            LALStatus status = empty_LALStatus;

            /* compute new SSB timings over all segments */
            if ( (cfvBuffer->multiSSB4V[n] = XLALGetMultiSSBtimesREAL4 ( multiDetStatesV->data[n], doppler->Alpha, doppler->Delta, doppler->refTime)) == NULL ) {
                XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
                XLALPrintError ( "%s: XLALGetMultiSSBtimesREAL4() failed. xlalErrno = %d.\n", fn, xlalErrno );
                XLAL_ERROR ( fn, XLAL_EFUNC );
            }

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

        } /* for n < numSegments */

    } /* if we could NOT reuse previously buffered quantities */


    //////////////////////////////////////////////////////////////////////////
    // Initializations and Assertions
    //

    // Check numIFOs is equal to hardcoded value
    if (constIFOs != 2)
    {
        XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
        XLALPrintError ("%s: numIFOs is not equal to 2. Consider upgrading the client\n", fn );
        XLAL_ERROR ( fn, XLAL_EINVAL );
    }

    for (n = 0; n < numSegments; n++)
    {
        // Check numIFOs is constant
        if (constIFOs != multiSFTsV->data[n]->length)
        {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError ("%s: numIFOs is not constant across data set\n", fn );
            XLAL_ERROR ( fn, XLAL_EINVAL );
        }

        for (k = 0; k < constIFOs; k++)
        {
            // Determine max sfts length
            if (maxNumSFTs < multiSFTsV->data[n]->data[k]->length)
            {
                maxNumSFTs = multiSFTsV->data[n]->data[k]->length;
            }

            // Assert constant fields consistency
            for (m=0; m<multiSFTsV->data[n]->data[k]->length; m++)
            {
                // Check deltaF is constant
                if (constSftsDataDeltaF != multiSFTsV->data[n]->data[k]->data[m].deltaF)
                {
                    XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
                    XLALPrintError ("%s: deltaF is not constant across data set\n", fn );
                    XLAL_ERROR ( fn, XLAL_EINVAL );
                }

                if (constSftsDataF0 != multiSFTsV->data[n]->data[k]->data[m].f0)
                {
                    XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
                    XLALPrintError ("%s: f0 is not constant across data set\n", fn );
                    XLAL_ERROR ( fn, XLAL_EINVAL );
                }

                if (constSftsDataDataLength != multiSFTsV->data[n]->data[k]->data[m].data->length)
                {
                    XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
                    XLALPrintError ("%s: sftsDataDataDataLength is not constant across data set\n", fn );
                    XLAL_ERROR ( fn, XLAL_EINVAL );
                }
            }
        }
    }

    // If no cuda_device_id specified yet, find device with max SM
    if ( cuda_device_id < 0 ) {
      int nodevices, deviceid, maxsmdevice=0;
      size_t maxsm=0;
      
      if (cudaSuccess != cudaGetDeviceCount(&nodevices) ) {
        XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
        XLALPrintError ("%s: cudaGetDeviceCount failed, no CUDA devices found\n", fn );
        XLAL_ERROR ( fn, XLAL_EINVAL );
      }

      for(deviceid=0; deviceid<nodevices; deviceid++) {
	if (cudaSuccess != cudaGetDeviceProperties (&curDevProps, deviceid) ) {
	  XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
	  XLALPrintError ("%s: cudaGetDeviceProperties failed\n", fn );
	  XLAL_ERROR ( fn, XLAL_EINVAL );
	}
	if (maxsm < curDevProps.sharedMemPerBlock) {
	  maxsm = curDevProps.sharedMemPerBlock;
	  maxsmdevice = deviceid;
	}
      }

      cuda_device_id = maxsmdevice;
    }

    fprintf (stderr, "Using CUDA device %d\n", cuda_device_id);

    if (cudaSuccess != cudaGetDeviceProperties (&curDevProps, cuda_device_id) )
    {
        XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
        XLALPrintError ("%s: cudaGetDeviceProperties failed\n", fn );
        XLAL_ERROR ( fn, XLAL_EINVAL );
    }

    // align to warp size
    maxNumSFTs = (maxNumSFTs + curDevProps.warpSize - 1) & (~(curDevProps.warpSize - 1));

    // Check that (sfts_length rounded to warp size) x numIFOs <= MAX_THREADS_IN_BLOCK
    if (maxNumSFTs * constIFOs > (UINT4)curDevProps.maxThreadsPerBlock)
    {
        XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
        XLALPrintError ("%s: cudaGetDeviceProperties failed\n", fn );
        XLAL_ERROR ( fn, XLAL_EINVAL );
    }

    //////////////////////////////////////////////////////////////////////////
    // Flatten structures
    //
    {
    #define UNIPITCH256(type, numElems)     ((numElems * sizeof(type) + 255) & (~255))

    // calculate pitch for Fstat output (256-bytes alignment)
    UINT4 unipitch_Fstat = UNIPITCH256(REAL4, numBins);

    // calculate pitch for sfts_data_data_data (256-bytes alignment)
    UINT4 unipitch_sfts_data_data_data = UNIPITCH256(COMPLEX8, constSftsDataDataLength);

    //////////////////////////////////////////////////////////////////////////
    // HOST Declaration

    REAL4 *Fstat;

    COMPLEX8 *ifo0_sfts_data_data_data;
    COMPLEX8 *ifo1_sfts_data_data_data;

    REAL4 *fkdot4_FreqMain;
    REAL4 *fkdot4_fkdot_0;

    REAL4 *ifo0_tSSB_DeltaT_int;
    REAL4 *ifo1_tSSB_DeltaT_int;
    REAL4 *ifo0_tSSB_DeltaT_rem;
    REAL4 *ifo1_tSSB_DeltaT_rem;
    REAL4 *ifo0_tSSB_TdotM1;
    REAL4 *ifo1_tSSB_TdotM1;

    REAL4 *ifo0_amcoe_a;
    REAL4 *ifo1_amcoe_a;
    REAL4 *ifo0_amcoe_b;
    REAL4 *ifo1_amcoe_b;

    UINT4 *ifo0_sfts_length;
    UINT4 *ifo1_sfts_length;
    REAL4 *antenna_pattern_Ad;
    REAL4 *antenna_pattern_Bd;
    REAL4 *antenna_pattern_Cd;
    REAL4 *antenna_pattern_DdInv;

    //REAL4 *debugREAL4;
    //UINT4 *debugUINT4;

    //////////////////////////////////////////////////////////////////////////
    // DEVICE Declaration

    //1
    REAL4 *dev_Fstat;

    //2
    COMPLEX8 *dev_ifo0_sfts_data_data_data;
    COMPLEX8 *dev_ifo1_sfts_data_data_data;

    //3
    REAL4 *dev_fkdot4_FreqMain;
    REAL4 *dev_fkdot4_fkdot0;

    //4
    REAL4 *dev_ifo0_tSSB_DeltaT_int;
    REAL4 *dev_ifo1_tSSB_DeltaT_int;
    REAL4 *dev_ifo0_tSSB_DeltaT_rem;
    REAL4 *dev_ifo1_tSSB_DeltaT_rem;
    REAL4 *dev_ifo0_tSSB_TdotM1;
    REAL4 *dev_ifo1_tSSB_TdotM1;

    //5
    REAL4 *dev_ifo0_amcoe_a;
    REAL4 *dev_ifo1_amcoe_a;
    REAL4 *dev_ifo0_amcoe_b;
    REAL4 *dev_ifo1_amcoe_b;

    UINT4 *dev_ifo0_sfts_length;
    UINT4 *dev_ifo1_sfts_length;
    REAL4 *dev_antenna_pattern_Ad;
    REAL4 *dev_antenna_pattern_Bd;
    REAL4 *dev_antenna_pattern_Cd;
    REAL4 *dev_antenna_pattern_DdInv;

    //REAL4 *dev_debugREAL4;
    //UINT4 *dev_debugUINT4;


    //////////////////////////////////////////////////////////////////////////
    // HOST Allocation

    Fstat                    = malloc(unipitch_Fstat * numSegments);

    ifo0_sfts_data_data_data = malloc(unipitch_sfts_data_data_data * maxNumSFTs * numSegments);
    ifo1_sfts_data_data_data = malloc(unipitch_sfts_data_data_data * maxNumSFTs * numSegments);

    fkdot4_FreqMain          = malloc(numBins * sizeof(REAL4));
    fkdot4_fkdot_0           = malloc(numBins * sizeof(REAL4));

    ifo0_tSSB_DeltaT_int     = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo1_tSSB_DeltaT_int     = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo0_tSSB_DeltaT_rem     = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo1_tSSB_DeltaT_rem     = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo0_tSSB_TdotM1         = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo1_tSSB_TdotM1         = malloc(numSegments * maxNumSFTs * sizeof(REAL4));

    ifo0_amcoe_a             = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo1_amcoe_a             = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo0_amcoe_b             = malloc(numSegments * maxNumSFTs * sizeof(REAL4));
    ifo1_amcoe_b             = malloc(numSegments * maxNumSFTs * sizeof(REAL4));

    ifo0_sfts_length         = malloc(numSegments * sizeof(UINT4));
    ifo1_sfts_length         = malloc(numSegments * sizeof(UINT4));
    antenna_pattern_Ad       = malloc(numSegments * sizeof(REAL4));
    antenna_pattern_Bd       = malloc(numSegments * sizeof(REAL4));
    antenna_pattern_Cd       = malloc(numSegments * sizeof(REAL4));
    antenna_pattern_DdInv    = malloc(numSegments * sizeof(REAL4));
//#define NUM_INTERMEDIATES_PER_THREAD    4
//    debugREAL4 = malloc((64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(REAL4));
//    debugUINT4 = malloc((64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(UINT4));

    //////////////////////////////////////////////////////////////////////////
    // HOST meminit

    memset(ifo0_sfts_data_data_data, 0, unipitch_sfts_data_data_data * maxNumSFTs * numSegments);
    memset(ifo1_sfts_data_data_data, 0, unipitch_sfts_data_data_data * maxNumSFTs * numSegments);

    memset(ifo0_tSSB_DeltaT_int, 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo1_tSSB_DeltaT_int, 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo0_tSSB_DeltaT_rem, 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo1_tSSB_DeltaT_rem, 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo0_tSSB_TdotM1    , 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo1_tSSB_TdotM1    , 0, numSegments * maxNumSFTs * sizeof(REAL4));

    memset(ifo0_amcoe_a        , 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo1_amcoe_a        , 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo0_amcoe_b        , 0, numSegments * maxNumSFTs * sizeof(REAL4));
    memset(ifo1_amcoe_b        , 0, numSegments * maxNumSFTs * sizeof(REAL4));

    //memset(debugUINT4        , 0xCC, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(UINT4));
    //memset(debugREAL4        , 0xCC, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(REAL4));

    /* REAL4-split frequency value */
    for (k = 0, Freq = Freq0; k < numBins; k++)
    {
        // TODO: Are we cool with this line?
        fkdot4_FreqMain[k] = (REAL4)((INT4)Freq);
        fkdot4_fkdot_0[k] = (REAL4)(Freq - (REAL8)fkdot4_FreqMain[k]);

        /* increase frequency by deltaF */
        Freq += deltaF;
    }

    for (n=0; n<numSegments; n++)
    {
        UINT4 i;

        ifo0_sfts_length[n] = multiSFTsV->data[n]->data[0]->length;
        ifo1_sfts_length[n] = multiSFTsV->data[n]->data[1]->length;

        for (i=0; i<ifo0_sfts_length[n]; i++)
        {
            memcpy(((char *)ifo0_sfts_data_data_data) + (n * maxNumSFTs + i) * unipitch_sfts_data_data_data,
                   multiSFTsV->data[n]->data[0]->data[i].data->data,
                   constSftsDataDataLength * sizeof(COMPLEX8));
        }
        for (i=0; i<ifo1_sfts_length[n]; i++)
        {
            memcpy(((char *)ifo1_sfts_data_data_data) + (n * maxNumSFTs + i) * unipitch_sfts_data_data_data,
                   multiSFTsV->data[n]->data[1]->data[i].data->data,
                   constSftsDataDataLength * sizeof(COMPLEX8));
        }

        memcpy(ifo0_tSSB_DeltaT_int + n * maxNumSFTs, cfvBuffer->multiSSB4V[n]->data[0]->DeltaT_int->data, multiSFTsV->data[n]->data[0]->length * sizeof(REAL4));
        memcpy(ifo1_tSSB_DeltaT_int + n * maxNumSFTs, cfvBuffer->multiSSB4V[n]->data[1]->DeltaT_int->data, multiSFTsV->data[n]->data[1]->length * sizeof(REAL4));
        memcpy(ifo0_tSSB_DeltaT_rem + n * maxNumSFTs, cfvBuffer->multiSSB4V[n]->data[0]->DeltaT_rem->data, multiSFTsV->data[n]->data[0]->length * sizeof(REAL4));
        memcpy(ifo1_tSSB_DeltaT_rem + n * maxNumSFTs, cfvBuffer->multiSSB4V[n]->data[1]->DeltaT_rem->data, multiSFTsV->data[n]->data[1]->length * sizeof(REAL4));
        memcpy(ifo0_tSSB_TdotM1     + n * maxNumSFTs, cfvBuffer->multiSSB4V[n]->data[0]->TdotM1->data,     multiSFTsV->data[n]->data[0]->length * sizeof(REAL4));
        memcpy(ifo1_tSSB_TdotM1     + n * maxNumSFTs, cfvBuffer->multiSSB4V[n]->data[1]->TdotM1->data,     multiSFTsV->data[n]->data[1]->length * sizeof(REAL4));

        memcpy(ifo0_amcoe_a + n * maxNumSFTs, cfvBuffer->multiAMcoefV[n]->data[0]->a->data, multiSFTsV->data[n]->data[0]->length * sizeof(REAL4));
        memcpy(ifo1_amcoe_a + n * maxNumSFTs, cfvBuffer->multiAMcoefV[n]->data[1]->a->data, multiSFTsV->data[n]->data[1]->length * sizeof(REAL4));
        memcpy(ifo0_amcoe_b + n * maxNumSFTs, cfvBuffer->multiAMcoefV[n]->data[0]->b->data, multiSFTsV->data[n]->data[0]->length * sizeof(REAL4));
        memcpy(ifo1_amcoe_b + n * maxNumSFTs, cfvBuffer->multiAMcoefV[n]->data[1]->b->data, multiSFTsV->data[n]->data[1]->length * sizeof(REAL4));

        antenna_pattern_Ad[n] = (REAL4)cfvBuffer->multiAMcoefV[n]->Mmunu.Ad;
        antenna_pattern_Bd[n] = (REAL4)cfvBuffer->multiAMcoefV[n]->Mmunu.Bd;
        antenna_pattern_Cd[n] = (REAL4)cfvBuffer->multiAMcoefV[n]->Mmunu.Cd;
        antenna_pattern_DdInv[n] = 1.0f / (REAL4)cfvBuffer->multiAMcoefV[n]->Mmunu.Dd;
    }

    //////////////////////////////////////////////////////////////////////////
    // DEVICE Allocation (gridDim.x = numBins, gridDim.x = numSegments)

    //1
    res = cudaMalloc((void **)&dev_Fstat, unipitch_Fstat * numSegments);

    //2
    res = cudaMalloc((void **)&dev_ifo0_sfts_data_data_data, unipitch_sfts_data_data_data * maxNumSFTs * numSegments);
    res = cudaMalloc((void **)&dev_ifo1_sfts_data_data_data, unipitch_sfts_data_data_data * maxNumSFTs * numSegments);

    //3
    res = cudaMalloc((void **)&dev_fkdot4_FreqMain, numBins * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_fkdot4_fkdot0, numBins * sizeof(REAL4));

    //4
    res = cudaMalloc((void **)&dev_ifo0_tSSB_DeltaT_int, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo1_tSSB_DeltaT_int, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo0_tSSB_DeltaT_rem, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo1_tSSB_DeltaT_rem, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo0_tSSB_TdotM1, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo1_tSSB_TdotM1, numSegments * maxNumSFTs * sizeof(REAL4));

    //5
    res = cudaMalloc((void **)&dev_ifo0_amcoe_a, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo1_amcoe_a, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo0_amcoe_b, numSegments * maxNumSFTs * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_ifo1_amcoe_b, numSegments * maxNumSFTs * sizeof(REAL4));

    //extra
    res = cudaMalloc((void **)&dev_ifo0_sfts_length, numSegments * sizeof(UINT4));
    res = cudaMalloc((void **)&dev_ifo1_sfts_length, numSegments * sizeof(UINT4));
    res = cudaMalloc((void **)&dev_antenna_pattern_Ad, numSegments * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_antenna_pattern_Bd, numSegments * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_antenna_pattern_Cd, numSegments * sizeof(REAL4));
    res = cudaMalloc((void **)&dev_antenna_pattern_DdInv, numSegments * sizeof(REAL4));

    //res = cudaMalloc((void **)&dev_debugREAL4, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(REAL4));
    //res = cudaMalloc((void **)&dev_debugUINT4, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(UINT4));

    //////////////////////////////////////////////////////////////////////////
    // DEVICE to HOST memcpy

    //2
    res = cudaMemcpy(dev_ifo0_sfts_data_data_data, ifo0_sfts_data_data_data, maxNumSFTs * unipitch_sfts_data_data_data * numSegments, cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_sfts_data_data_data, ifo1_sfts_data_data_data, maxNumSFTs * unipitch_sfts_data_data_data * numSegments, cudaMemcpyHostToDevice);

    //3
    res = cudaMemcpy(dev_fkdot4_FreqMain, fkdot4_FreqMain, numBins * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_fkdot4_fkdot0,   fkdot4_fkdot_0,  numBins * sizeof(REAL4), cudaMemcpyHostToDevice);

    //4
    res = cudaMemcpy(dev_ifo0_tSSB_DeltaT_int, ifo0_tSSB_DeltaT_int, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_tSSB_DeltaT_int, ifo1_tSSB_DeltaT_int, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo0_tSSB_DeltaT_rem, ifo0_tSSB_DeltaT_rem, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_tSSB_DeltaT_rem, ifo1_tSSB_DeltaT_rem, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo0_tSSB_TdotM1,     ifo0_tSSB_TdotM1,     numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_tSSB_TdotM1,     ifo1_tSSB_TdotM1,     numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);

    //5
    res = cudaMemcpy(dev_ifo0_amcoe_a, ifo0_amcoe_a, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_amcoe_a, ifo1_amcoe_a, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo0_amcoe_b, ifo0_amcoe_b, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_amcoe_b, ifo1_amcoe_b, numSegments * maxNumSFTs * sizeof(REAL4), cudaMemcpyHostToDevice);

    //extra
    res = cudaMemcpy(dev_ifo0_sfts_length, ifo0_sfts_length, numSegments * sizeof(UINT4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_ifo1_sfts_length, ifo1_sfts_length, numSegments * sizeof(UINT4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_antenna_pattern_Ad, antenna_pattern_Ad, numSegments * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_antenna_pattern_Bd, antenna_pattern_Bd, numSegments * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_antenna_pattern_Cd, antenna_pattern_Cd, numSegments * sizeof(REAL4), cudaMemcpyHostToDevice);
    res = cudaMemcpy(dev_antenna_pattern_DdInv, antenna_pattern_DdInv, numSegments * sizeof(REAL4), cudaMemcpyHostToDevice);

    res = cudaMemset(dev_Fstat, 0, unipitch_Fstat * numSegments);

    //res = cudaMemcpy(dev_debugREAL4, debugREAL4, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(REAL4), cudaMemcpyHostToDevice);
    //res = cudaMemcpy(dev_debugUINT4, debugUINT4, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(UINT4), cudaMemcpyHostToDevice);

    //////////////////////////////////////////////////////////////////////////
    // Kernel Launch

    HostWrapperCUDAComputeFstatFaFb/*_debug*/(dev_Fstat,
                                    unipitch_Fstat,

                                    constSftsDataDataLength,
                                    dev_ifo0_sfts_data_data_data,
                                    dev_ifo1_sfts_data_data_data,
                                    unipitch_sfts_data_data_data,
                                    constTsft,
                                    constDFreq,
                                    constFreqIndex0,

                                    dev_fkdot4_FreqMain,
                                    dev_fkdot4_fkdot0,
                                    fkdot4ex,

                                    dev_ifo0_tSSB_DeltaT_int,
                                    dev_ifo1_tSSB_DeltaT_int,
                                    dev_ifo0_tSSB_DeltaT_rem,
                                    dev_ifo1_tSSB_DeltaT_rem,
                                    dev_ifo0_tSSB_TdotM1,
                                    dev_ifo1_tSSB_TdotM1,

                                    dev_ifo0_amcoe_a,
                                    dev_ifo1_amcoe_a,
                                    dev_ifo0_amcoe_b,
                                    dev_ifo1_amcoe_b,

                                    dev_ifo0_sfts_length,
                                    dev_ifo1_sfts_length,

                                    dev_antenna_pattern_Ad,
                                    dev_antenna_pattern_Bd,
                                    dev_antenna_pattern_Cd,
                                    dev_antenna_pattern_DdInv,

                                    Dterms,
                                    //dev_debugREAL4,
                                    //dev_debugUINT4,

                                    numBins,
                                    numSegments);
    res = cudaGetLastError();
    res = cudaThreadSynchronize();

    //////////////////////////////////////////////////////////////////////////
    // HOST to DEVICE memcpy

    res = cudaMemcpy(Fstat, dev_Fstat, unipitch_Fstat * numSegments, cudaMemcpyDeviceToHost);
    for (n=0; n<numSegments; n++)
    {
        for (k=0; k<numBins; k++)
        {
            fstatBandV->data[n].data->data[k] = ((REAL4 *)(((char *)Fstat) + n * unipitch_Fstat))[k];
        }
    }

    //res = cudaMemcpy(debugREAL4, dev_debugREAL4, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(REAL4), cudaMemcpyDeviceToHost);
    //res = cudaMemcpy(debugUINT4, dev_debugUINT4, (64*2* NUM_INTERMEDIATES_PER_THREAD + 100) * sizeof(REAL4), cudaMemcpyDeviceToHost);

    //////////////////////////////////////////////////////////////////////////
    // DEVICE Deallocation

    //1
    cudaFree(dev_Fstat);

    //2
    cudaFree(dev_ifo0_sfts_data_data_data);
    cudaFree(dev_ifo1_sfts_data_data_data);

    //3
    cudaFree(dev_fkdot4_FreqMain);
    cudaFree(dev_fkdot4_fkdot0);

    //4
    cudaFree(dev_ifo0_tSSB_DeltaT_int);
    cudaFree(dev_ifo1_tSSB_DeltaT_int);
    cudaFree(dev_ifo0_tSSB_DeltaT_rem);
    cudaFree(dev_ifo1_tSSB_DeltaT_rem);
    cudaFree(dev_ifo0_tSSB_TdotM1);
    cudaFree(dev_ifo1_tSSB_TdotM1);

    //5
    cudaFree(dev_ifo0_amcoe_a);
    cudaFree(dev_ifo1_amcoe_a);
    cudaFree(dev_ifo0_amcoe_b);
    cudaFree(dev_ifo1_amcoe_b);

    //extra
    cudaFree(dev_ifo0_sfts_length);
    cudaFree(dev_ifo1_sfts_length);
    cudaFree(dev_antenna_pattern_Ad);
    cudaFree(dev_antenna_pattern_Bd);
    cudaFree(dev_antenna_pattern_Cd);
    cudaFree(dev_antenna_pattern_DdInv);

    //////////////////////////////////////////////////////////////////////////
    // HOST Deallocation

    free(Fstat);

    free(ifo0_sfts_data_data_data);
    free(ifo1_sfts_data_data_data);

    free(fkdot4_FreqMain);
    free(fkdot4_fkdot_0);

    free(ifo0_tSSB_DeltaT_int);
    free(ifo1_tSSB_DeltaT_int);
    free(ifo0_tSSB_DeltaT_rem);
    free(ifo1_tSSB_DeltaT_rem);
    free(ifo0_tSSB_TdotM1);
    free(ifo1_tSSB_TdotM1);

    free(ifo0_amcoe_a);
    free(ifo1_amcoe_a);
    free(ifo0_amcoe_b);
    free(ifo1_amcoe_b);

    free(ifo0_sfts_length);
    free(ifo1_sfts_length);
    free(antenna_pattern_Ad);
    free(antenna_pattern_Bd);
    free(antenna_pattern_Cd);
    free(antenna_pattern_DdInv);
    }

    return XLAL_SUCCESS;

} /* XLALComputeFStatFreqBandVector_cuda() */


