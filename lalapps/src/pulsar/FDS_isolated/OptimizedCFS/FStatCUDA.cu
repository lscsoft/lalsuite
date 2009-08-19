/*
 * Copyright (c) 2009, NVIDIA Corporation ("NVIDIA") All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright notice, 
 *       this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in the 
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 */

//////////////////////////////////////////////////////////////////////////
// Magic values
//////////////////////////////////////////////////////////////////////////

#define PULSAR_MAX_SPINS    7


//////////////////////////////////////////////////////////////////////////
// CUDA type redefinitions
//////////////////////////////////////////////////////////////////////////

typedef float           REAL4;
typedef int             INT4;
typedef unsigned int    UINT4;


/** Single-precision floating-point complex number (8 bytes total) */
typedef struct tagCOMPLEX8
{
    REAL4 re;                                                                   /**< The real part. */
    REAL4 im;                                                                   /**< The imaginary part. */
} COMPLEX8;


/**< REAL4 version of pulsar spins fkdot[] */
typedef struct
{
    REAL4 FreqMain;                                                             /**< "main" part of frequency fkdot[0], normally just the integral part */
    REAL4 fkdot[PULSAR_MAX_SPINS];                                              /**< remaining spin-parameters, including *fractional* part of Freq = fkdot[0] */
} PulsarSpinsREAL4;


typedef struct {
    REAL4 fkdot16[PULSAR_MAX_SPINS-1];                                          /**< remaining spin-parameters, excluding *fractional* part of Freq = fkdot[0] */
} PulsarSpinsExREAL4;


/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb 
 * (for ML-estimators). NOTE: this is simply a REAL4 version of Fcomponents.
 */
typedef struct {
    REAL4 F;                                                                    /**< F-statistic value */
    COMPLEX8 Fa;                                                                /**< complex amplitude Fa */
    COMPLEX8 Fb;                                                                /**< complex amplitude Fb */
} FcomponentsREAL4;


//////////////////////////////////////////////////////////////////////////
// Macros
//////////////////////////////////////////////////////////////////////////

#define NAN_UINT4                   0x7fc00000
#define LD_SMALL4                   ((REAL4)2.0e-4)                             /**< "small" number for REAL4*/
#define PI_F                        3.1415926535897932384626f
#define TWOPI_F                     (PI_F * 2)
#define INV_TWOPI_F                 (1.0f / TWOPI_F)


#define SINCOS_TRIM_X(y,x)          y = (x) - floorf(x)
//#define SINCOS_TRIM_X(y,x)          y = (x)

#define SINCOS_2PI_TRIMMED(s,c,x)   sincosf(TWOPI_F*(x), const_cast<float *>(s), const_cast<float *>(c))
//#define SINCOS_2PI_TRIMMED(s,c,x)   {*s = sinf(TWOPI_F*(x)); *c = cosf(TWOPI_F*(x));}

#define REMAINDER(x)                fmodf(x, 1.0f)
//#define REMAINDER(x)                ( (x) - (INT4)(x) )

//#define RECIPROCAL(x)               ( 1.0f / (x) )
#define RECIPROCAL(x)               ( __frcp_rn (x) )

//#define REAL4TOINT4(x)              ( (INT4) (x) )
#define REAL4TOINT4(x)              ( __float2int_rz(x) )

#define FMOD(x,y)                   fmodf(x,y)

#define SQUARE(x)                   ( (x) * (x) )


//////////////////////////////////////////////////////////////////////////
// Constants
//////////////////////////////////////////////////////////////////////////

__constant__ REAL4 inv_fact[PULSAR_MAX_SPINS] = { 1.0f, 1.0f, (1.0/2.0), (1.0f/6.0f), (1.0f/24.0f), (1.0f/120.0f), (1.0f/720.0f) };


//////////////////////////////////////////////////////////////////////////
// Kernels
//////////////////////////////////////////////////////////////////////////

// TODO: move Dterms to template
// TODO: eliminate spdnOrder variable, unroll loop
// TODO: check how global mem is read, consider reading with texture fetches
__global__ void CUDAComputeFstatFaFb(REAL4 *Fstat,
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

                                     UINT4 Dterms)
{
    __shared__ FcomponentsREAL4 FaFb_components[2][64];

    COMPLEX8 Fa, Fb;
    Fa.re = 0.0f;
    Fa.im = 0.0f;
    Fb.re = 0.0f;
    Fb.im = 0.0f;

    int curSegment = blockIdx.y;
    int curBin = blockIdx.x;
    int maxSfts = blockDim.x;
    int curSfts = threadIdx.x;
    int curIFO = threadIdx.y;

    if (curIFO == 0 && curSfts < ifo0_sfts_length[curSegment] ||
        curIFO == 1 && curSfts < ifo1_sfts_length[curSegment])
    {
        //2D arrays pack
        COMPLEX8 *sfts_data_data_data;
        REAL4 *tSSB_DeltaT_int;
        REAL4 *tSSB_DeltaT_rem;
        REAL4 *tSSB_TdotM1;
        REAL4 *amcoe_a;
        REAL4 *amcoe_b;

        if (threadIdx.y)
        {
            sfts_data_data_data = (COMPLEX8 *)(((char *)ifo1_sfts_data_data_data) + (curSegment * maxSfts + curSfts) * sfts_data_data_data_pitch);
            tSSB_DeltaT_int = ifo1_tSSB_DeltaT_int + curSegment * maxSfts;
            tSSB_DeltaT_rem = ifo1_tSSB_DeltaT_rem + curSegment * maxSfts + curSfts;
            tSSB_TdotM1 = ifo1_tSSB_TdotM1 + curSegment * maxSfts + curSfts;
            amcoe_a = ifo1_amcoe_a + curSegment * maxSfts + curSfts;
            amcoe_b = ifo1_amcoe_b + curSegment * maxSfts + curSfts;
        }
        else
        {
            sfts_data_data_data = (COMPLEX8 *)(((char *)ifo0_sfts_data_data_data) + (curSegment * maxSfts + curSfts) * sfts_data_data_data_pitch);
            tSSB_DeltaT_int = ifo0_tSSB_DeltaT_int + curSegment * maxSfts;
            tSSB_DeltaT_rem = ifo0_tSSB_DeltaT_rem + curSegment * maxSfts + curSfts;
            tSSB_TdotM1 = ifo0_tSSB_TdotM1 + curSegment * maxSfts + curSfts;
            amcoe_a = ifo0_amcoe_a + curSegment * maxSfts + curSfts;
            amcoe_b = ifo0_amcoe_b + curSegment * maxSfts + curSfts;
        }

        INT4 freqIndex1 = freqIndex0 + sfts_data_data_length;                   /* index of last frequency-bin in SFTs */

        REAL4 f0 = fkdot4_FreqMain[curBin];
        REAL4 df = fkdot4_fkdot0[curBin];
        REAL4 tau = RECIPROCAL(df);
        REAL4 Freq = f0 + df;

        //TODO: coalesce
        REAL4 T0offs = tSSB_DeltaT_int[0];
        tSSB_DeltaT_int += curSfts;
        REAL4 phi_alpha_offs = df * FMOD(T0offs, tau);

        UINT4 spdnOrder;

        /* find highest non-zero spindown-entry */
        //TODO: this can be actually avoided by direct checking fkdot4->fkdot[i] in a loop of interest
        for (spdnOrder = 0; spdnOrder < PULSAR_MAX_SPINS - 1; spdnOrder++)
        {
            if (!fkdot4ex.fkdot16[spdnOrder])
            {
                break;
            }
        }

        /* Process one SFT per thread  */
        {
            COMPLEX8 *Xalpha = sfts_data_data_data;                             /* pointer to current SFT-data */
            COMPLEX8 *Xalpha_l;                                                 /* pointer to frequency-bin k in current SFT */
            REAL4 a_alpha, b_alpha;
            REAL4 s_alpha, c_alpha;                                             /* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
            REAL4 realQ, imagQ;                                                 /* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
            REAL4 realXP, imagXP;                                               /* Re/Im of sum_k X_ak * P_ak */
            REAL4 realQXP, imagQXP;                                             /* Re/Im of Q_alpha R_alpha */
            REAL4 lambda_alpha;
            INT4 kstar;                                                         /* central frequency-bin k* = round(xhat_alpha) */
            REAL4 kappa_max, kappa_star;
            INT4 k0, k1;

            /* ----- calculate kappa_max and lambda_alpha */
            {
                /* 1st oder: s = 0 */
                REAL4 TdotM1 = *tSSB_TdotM1;                                    /* defined as Tdot_al - 1 */
                REAL4 dT = *tSSB_DeltaT_rem;

                REAL4 T0 = *tSSB_DeltaT_int - T0offs;
                REAL4 deltaT = *tSSB_DeltaT_int + dT;

                /* phi_alpha = f * Tas; */
                REAL4 T0rem = FMOD( T0, tau );

                REAL4 phi_alpha_rem = phi_alpha_offs + f0 * dT + T0 * df + df * dT;

                REAL4 Dphi_alpha_int = f0;
                REAL4 Dphi_alpha_rem = df + Freq * TdotM1;

                /* higher-order spindowns */
                REAL4 Tas = deltaT;
                for (UINT4 s=0; s < spdnOrder; s++)
                {
                    REAL4 fsdot = fkdot4ex.fkdot16[s];
                    Dphi_alpha_rem += fsdot * Tas * inv_fact[s+1];              /* here: Tas = DT^s */
                    Tas *= deltaT;                                              /* now:  Tas = DT^(s+1) */
                    phi_alpha_rem += fsdot * Tas * inv_fact[s+2];
                } /* for s <= spdnOrder */

                /* Step 3: apply global factor of Tsft to complete Dphi_alpha */
                Dphi_alpha_int *= Tsft;
                Dphi_alpha_rem *= Tsft;

                REAL4 tmp = REMAINDER( 0.5f * Dphi_alpha_int ) + REMAINDER ( 0.5f * Dphi_alpha_rem );
                lambda_alpha = phi_alpha_rem - tmp;

                /* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
                kstar = REAL4TOINT4(Dphi_alpha_int) + REAL4TOINT4(Dphi_alpha_rem);
                kappa_star = REMAINDER(Dphi_alpha_int) + REMAINDER(Dphi_alpha_rem);
                kappa_max = kappa_star + 1.0f * Dterms - 1.0f;

                /* ----- check that required frequency-bins are found in the SFTs ----- */
                k0 = kstar - Dterms + 1;
                k1 = k0 + 2 * Dterms - 1;
                if ((k0 < freqIndex0) || (k1 > freqIndex1))
                {
                    ((int*)&(FaFb_components[threadIdx.y][threadIdx.x].Fa.re))[0] = NAN_UINT4;
                }

            } /* compute kappa_star, lambda_alpha */

            /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
            * the trig-functions need to be calculated only once!
            * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
            * closest to zero and will pose no numerical difficulties !
            */

            /* ---------- calculate the (truncated to Dterms) sum over k ---------- */
            Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */

            realXP = 0;
            imagXP = 0;

            /* if no danger of denominator -> 0 */
            if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
            {
                /* improved hotloop algorithm by Fekete Akos:
                * take out repeated divisions into a single common denominator,
                * plus use extra cleverness to compute the nominator efficiently...
                */
                REAL4 Sn = (*Xalpha_l).re;
                REAL4 Tn = (*Xalpha_l).im;
                REAL4 pn = kappa_max;
                REAL4 qn = pn;
                REAL4 U_alpha, V_alpha;

                /* recursion with 2*Dterms steps */
                for (UINT4 l = 1; l < 2*Dterms; l ++ )
                {
                    Xalpha_l ++;
                    pn = pn - 1.0f;                                             /* p_(n+1) */
                    Sn = pn * Sn + qn * (*Xalpha_l).re;                         /* S_(n+1) */
                    Tn = pn * Tn + qn * (*Xalpha_l).im;                         /* T_(n+1) */
                    qn *= pn;                                                   /* q_(n+1) */
                } /* for l <= 2*Dterms */

                REAL4 qn_inv = RECIPROCAL(qn);
                U_alpha = Sn * qn_inv;
                V_alpha = Tn * qn_inv;

                SINCOS_2PI_TRIMMED(&s_alpha, &c_alpha, kappa_star);
                c_alpha -= 1.0f;

                realXP = s_alpha * U_alpha - c_alpha * V_alpha;
                imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

                {
                    REAL4 _lambda_alpha = -lambda_alpha;

                    //TODO: consider removing
                    SINCOS_TRIM_X (_lambda_alpha, _lambda_alpha);

                    SINCOS_2PI_TRIMMED(&imagQ, &realQ, _lambda_alpha);
                }
            }                                                                   /* if |remainder| > LD_SMALL4 */
            else
            {                                                                   /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
                UINT4 ind0;

                /* realQ/imagQ are calculated in the hotloop in the other case; in this case we have to do it too */
                REAL4 _lambda_alpha = -lambda_alpha;

                //TODO: consider removing
                SINCOS_TRIM_X(_lambda_alpha, _lambda_alpha);

                SINCOS_2PI_TRIMMED(&imagQ, &realQ, _lambda_alpha);

                if (kappa_star <= LD_SMALL4)
                {
                    ind0 = Dterms - 1;
                }
                else
                {
                    ind0 = Dterms;
                }

                realXP = TWOPI_F * Xalpha_l[ind0].re;
                imagXP = TWOPI_F * Xalpha_l[ind0].im;
            }                                                                   /* if |remainder| <= LD_SMALL4 */

            realQXP = realQ * realXP - imagQ * imagXP;
            imagQXP = realQ * imagXP + imagQ * realXP;

            /* we're done: ==> combine these into Fa and Fb */
            a_alpha = (*amcoe_a);
            b_alpha = (*amcoe_b);

            Fa.re += a_alpha * realQXP;
            Fa.im += a_alpha * imagQXP;

            Fb.re += b_alpha * realQXP;
            Fb.im += b_alpha * imagQXP;
        }

    } //endif threadIdx

    __syncthreads();

    /* store intermediate result in shared memory */
    FaFb_components[threadIdx.y][threadIdx.x].Fa.re = Fa.re;
    FaFb_components[threadIdx.y][threadIdx.x].Fa.im = Fa.im;
    FaFb_components[threadIdx.y][threadIdx.x].Fb.re = Fb.re;
    FaFb_components[threadIdx.y][threadIdx.x].Fb.im = Fb.im;

    // Perform reduction. Assuming blockDim = [64,2,1]
    __syncthreads();
    if (threadIdx.x < 32)
    {
        FaFb_components[threadIdx.y][threadIdx.x].Fa.re += FaFb_components[threadIdx.y][threadIdx.x + 32].Fa.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fa.im += FaFb_components[threadIdx.y][threadIdx.x + 32].Fa.im;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.re += FaFb_components[threadIdx.y][threadIdx.x + 32].Fb.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.im += FaFb_components[threadIdx.y][threadIdx.x + 32].Fb.im;
    }
    __syncthreads();
    if (threadIdx.x < 16)
    {
        FaFb_components[threadIdx.y][threadIdx.x].Fa.re += FaFb_components[threadIdx.y][threadIdx.x + 16].Fa.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fa.im += FaFb_components[threadIdx.y][threadIdx.x + 16].Fa.im;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.re += FaFb_components[threadIdx.y][threadIdx.x + 16].Fb.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.im += FaFb_components[threadIdx.y][threadIdx.x + 16].Fb.im;
    }
    __syncthreads();
    if (threadIdx.x < 8)
    {
        FaFb_components[threadIdx.y][threadIdx.x].Fa.re += FaFb_components[threadIdx.y][threadIdx.x + 8].Fa.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fa.im += FaFb_components[threadIdx.y][threadIdx.x + 8].Fa.im;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.re += FaFb_components[threadIdx.y][threadIdx.x + 8].Fb.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.im += FaFb_components[threadIdx.y][threadIdx.x + 8].Fb.im;
    }
    __syncthreads();
    if (threadIdx.x < 4)
    {
        FaFb_components[threadIdx.y][threadIdx.x].Fa.re += FaFb_components[threadIdx.y][threadIdx.x + 4].Fa.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fa.im += FaFb_components[threadIdx.y][threadIdx.x + 4].Fa.im;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.re += FaFb_components[threadIdx.y][threadIdx.x + 4].Fb.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.im += FaFb_components[threadIdx.y][threadIdx.x + 4].Fb.im;
    }
    __syncthreads();
    if (threadIdx.x < 2)
    {
        FaFb_components[threadIdx.y][threadIdx.x].Fa.re += FaFb_components[threadIdx.y][threadIdx.x + 2].Fa.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fa.im += FaFb_components[threadIdx.y][threadIdx.x + 2].Fa.im;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.re += FaFb_components[threadIdx.y][threadIdx.x + 2].Fb.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.im += FaFb_components[threadIdx.y][threadIdx.x + 2].Fb.im;
    }
    __syncthreads();
    if (threadIdx.x < 1)
    {
        FaFb_components[threadIdx.y][threadIdx.x].Fa.re += FaFb_components[threadIdx.y][threadIdx.x + 1].Fa.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fa.im += FaFb_components[threadIdx.y][threadIdx.x + 1].Fa.im;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.re += FaFb_components[threadIdx.y][threadIdx.x + 1].Fb.re;
        FaFb_components[threadIdx.y][threadIdx.x].Fb.im += FaFb_components[threadIdx.y][threadIdx.x + 1].Fb.im;
    }
    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0)
    {
        FaFb_components[0][0].Fa.re += FaFb_components[1][0].Fa.re;
        FaFb_components[0][0].Fa.im += FaFb_components[1][0].Fa.im;
        FaFb_components[0][0].Fb.re += FaFb_components[1][0].Fb.re;
        FaFb_components[0][0].Fb.im += FaFb_components[1][0].Fb.im;

        // Compute final FaFb value
        FaFb_components[0][0].Fa.re *= INV_TWOPI_F;
        FaFb_components[0][0].Fa.im *= INV_TWOPI_F;
        FaFb_components[0][0].Fb.re *= INV_TWOPI_F;
        FaFb_components[0][0].Fb.im *= INV_TWOPI_F;

        // compute final F-statistic value

        // NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
        // therefore there is a factor of 2 difference with respect to the equations in JKS, which
        // where based on the single-sided PSD.

        ((REAL4 *)(((char *)Fstat) + curSegment * Fstat_pitch))[curBin] =
                 InvDd[curSegment] * ( Bd[curSegment] * ( SQUARE(FaFb_components[0][0].Fa.re) + SQUARE(FaFb_components[0][0].Fa.im) )
                                     + Ad[curSegment] * ( SQUARE(FaFb_components[0][0].Fb.re) + SQUARE(FaFb_components[0][0].Fb.im) )
            - 2.0f * Cd[curSegment] *( FaFb_components[0][0].Fa.re * FaFb_components[0][0].Fb.re + FaFb_components[0][0].Fa.im * FaFb_components[0][0].Fb.im )
              );
    }
}


//////////////////////////////////////////////////////////////////////////
// Wrappers
//////////////////////////////////////////////////////////////////////////

extern "C" void HostWrapperCUDAComputeFstatFaFb(REAL4 *Fstat,
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
                                                UINT4 numSegments)
{
    dim3 blockSize;
    dim3 gridSize;
    blockSize.x = 64; //max(numSFTs)
    blockSize.y = 2;  //IFOs
    blockSize.z = 1;
    gridSize.x = numBins;
    gridSize.y = numSegments;
    gridSize.z = 1;

    CUDAComputeFstatFaFb<<<gridSize,blockSize>>>(Fstat,
                                                 Fstat_pitch,

                                                 sfts_data_data_length,
                                                 ifo0_sfts_data_data_data,
                                                 ifo1_sfts_data_data_data,
                                                 sfts_data_data_data_pitch,
                                                 Tsft,
                                                 dFreq,
                                                 freqIndex0,

                                                 fkdot4_FreqMain,
                                                 fkdot4_fkdot0,
                                                 fkdot4ex,

                                                 ifo0_tSSB_DeltaT_int,
                                                 ifo1_tSSB_DeltaT_int,
                                                 ifo0_tSSB_DeltaT_rem,
                                                 ifo1_tSSB_DeltaT_rem,
                                                 ifo0_tSSB_TdotM1,
                                                 ifo1_tSSB_TdotM1,

                                                 ifo0_amcoe_a,
                                                 ifo1_amcoe_a,
                                                 ifo0_amcoe_b,
                                                 ifo1_amcoe_b,

                                                 ifo0_sfts_length,
                                                 ifo1_sfts_length,

                                                 Ad,
                                                 Bd,
                                                 Cd,
                                                 InvDd,

                                                 Dterms);
}
