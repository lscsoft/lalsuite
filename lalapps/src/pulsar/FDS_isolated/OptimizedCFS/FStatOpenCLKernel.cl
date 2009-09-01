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

#define PULSAR_MAX_SPINS   7
#define DTERMS             8

#if USE_OPENCL_KERNEL_CPU

#  define __kernel
#  define __global
#  define __local
#  define __constant const

#else // if use OpenCL kernel on GPU

typedef int           INT4;
typedef unsigned int  UINT4;
typedef float         REAL4;

// TODO rename to more descriptive
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

//typedef float2        REAL42;    // Needed to represent array of frequencies, split into REAL4 pairs
//typedef float4        REAL44;    // Vector with 4 float numbers: {A, B, C, D}
//typedef struct {REAL4 s[16];} PulsarSpins16; // Pulsar spindown values (up to 15)

/** Single-precision floating-point complex number (8 bytes total) */
typedef struct tagCOMPLEX8
{
    REAL4 re;                                                                   /**< The real part. */
    REAL4 im;                                                                   /**< The imaginary part. */
} COMPLEX8;


/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb 
 * (for ML-estimators). NOTE: this is simply a REAL4 version of Fcomponents.
 */
typedef struct {
    REAL4 F;                                                                    /**< F-statistic value */
    COMPLEX8 Fa;                                                                /**< complex amplitude Fa */
    COMPLEX8 Fb;                                                                /**< complex amplitude Fb */
} FcomponentsREAL4;

#endif // if USE_OPENCL_KERNEL_CPU

/*
 * Macros
 */

// #define numSegs 121
#define NUM_IFOS 2
// #define maxSfts 64
// #define sftLen 254

#define NAN_UINT4                   0x7fc00000
#define LD_SMALL4                   ((REAL4)2.0e-4)                             /**< "small" number for REAL4*/
#define PI_F                        3.1415926535897932384626f
#define TWOPI_F                     (PI_F * 2)
#define INV_TWOPI_F                 (1.0f / TWOPI_F)


#if USE_OPENCL_KERNEL_CPU
#  define SINCOS_TRIM_X(y,x)          y = (x) - floorf(x)
#else
#  define SINCOS_TRIM_X(y,x)          y = (x) - floor(x)
#endif
//#define SINCOS_TRIM_X(y,x)          y = (x)

#if USE_OPENCL_KERNEL_CPU
#  define SINCOS_2PI_TRIMMED(s,c,x)   sin_cos_2PI_LUT_REAL4(s,c,x);
//#  define SINCOS_2PI_TRIMMED(s,c,x)   {*s = sin(TWOPI_F*(x)); *c = cos(TWOPI_F*(x));}
#else
#  define SINCOS_2PI_TRIMMED(s,c,x)   *s = sincos(TWOPI_F*(x), (c))
//#  define SINCOS_2PI_TRIMMED(s,c,x)   {*s = sinf(TWOPI_F*(x)); *c = cosf(TWOPI_F*(x));}
#endif
//#define SINCOS_2PI_TRIMMED(s,c,x)   sincosf(TWOPI_F*(x), const_cast<float *>(s), const_cast<float *>(c))

//#define REMAINDER(x)                fmod(x, 1.0f)
#define REMAINDER(x)                ( (x) - (INT4)(x) )

#define RECIPROCAL(x)               ( 1.0f / (x) )
//#define RECIPROCAL(x)               ( __frcp_rn (x) )

#define REAL4TOINT4(x)              ( (INT4) (x) )
//#define REAL4TOINT4(x)              ( __float2int_rz(x) )

#if USE_OPENCL_KERNEL_CPU
#  define FMOD(x,y)                   fmodf(x,y)
#else
#  define FMOD(x,y)                   fmod(x,y)
#endif

#define SQ(x)                   ( (x) * (x) )


/*
 * Constants
 */

__constant REAL4 invfact[PULSAR_MAX_SPINS] = { 1.0f, 1.0f, (1.0/2.0), (1.0f/6.0f), (1.0f/24.0f), (1.0f/120.0f), (1.0f/720.0f) };

/*
 * OpenCL compute kernel 
 */
__kernel void OpenCLComputeFstatFaFb(__global REAL4 *Fstat,
#if USE_OPENCL_KERNEL_CPU
                                     UINT4 curSegment,
                                     UINT4 curBin,
                                     UINT4 maxSfts,
                                     UINT4 curSfts,
                                     UINT4 curIFO,
                                     UINT4 numSegs,
#endif

                                     __global COMPLEX8 *arg_sfts_data,
                                     __global UINT4 *arg_sfts_length,
                                     UINT4 sftLen,
                                     REAL4 Tsft,
                                     REAL4 dFreq,
                                     INT4 freqIndex0,

                                     __global REAL42 *fkdot4,   //TODO: simplify here
                                     PulsarSpins16 fkdot4ex,

                                     __global REAL4 *arg_tSSB_DeltaT_int,
                                     __global REAL4 *arg_tSSB_DeltaT_rem,
                                     __global REAL4 *arg_tSSB_TdotM1,

                                     __global REAL4 *arg_amcoe_a,
                                     __global REAL4 *arg_amcoe_b,

                                     __global REAL44 *ABCInvDd,
                                     __local FcomponentsREAL4 *FaFb_components)
{
    COMPLEX8 Fa, Fb;
    Fa.re = 0.0f;
    Fa.im = 0.0f;
    Fb.re = 0.0f;
    Fb.im = 0.0f;

#if !(USE_OPENCL_KERNEL_CPU)
    int curSegment = get_group_id(0);
    int curBin     = get_group_id(1);
    int maxSfts    = get_global_size(0) / get_num_groups(0);
    int curSfts    = get_local_id(0);
    int curIFO     = get_local_id(1);
    int numSegs    = get_num_groups(0);
#endif    

    int offset = curIFO + NUM_IFOS * curSegment;
    if (curSfts < arg_sfts_length[offset])
    {
        //2D arrays pack
        __global COMPLEX8 *sfts_data;           // m + sftLen * (s + NAM_NUM_SFTS * (X + NUM_IFOS * n))
        __global REAL4 *tSSB_DeltaT_int;        // s + maxSfts * (X + NUM_IFOS * n)
        __global REAL4 *tSSB_DeltaT_rem;
        __global REAL4 *tSSB_TdotM1;
        __global REAL4 *amcoe_a;
        __global REAL4 *amcoe_b;

        offset *= maxSfts;
        offset += curSfts;
        sfts_data = arg_sfts_data + sftLen * offset;
        tSSB_DeltaT_int = arg_tSSB_DeltaT_int + offset;
        tSSB_DeltaT_rem = arg_tSSB_DeltaT_rem + offset;
        tSSB_TdotM1 = arg_tSSB_TdotM1 + offset;
        amcoe_a = arg_amcoe_a + offset;
        amcoe_b = arg_amcoe_b + offset;

        INT4 freqIndex1 = freqIndex0 + sftLen;                   // index of last frequency-bin in SFTs

        REAL4 f0 = fkdot4[curBin].FreqMain;
        REAL4 df = fkdot4[curBin].fkdot0;
        REAL4 tau = RECIPROCAL(df);
        REAL4 Freq = f0 + df;

//        //TODO: coalesce
//        REAL4 T0offs = tSSB_DeltaT_int[0];
//        tSSB_DeltaT_int += curSfts;
//        REAL4 phi_alpha_offs = df * FMOD(T0offs, tau);

        UINT4 spdnOrder = (UINT4)fkdot4ex.s[0];

        // Process one SFT per thread 
        {
            __global COMPLEX8 *Xalpha = sfts_data;                              // pointer to current SFT-data
            __global COMPLEX8 *Xalpha_l;                                        // pointer to frequency-bin k in current SFT
            REAL4 a_alpha, b_alpha;
            REAL4 s_alpha, c_alpha;                                             // sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1)
            REAL4 realQ, imagQ;                                                 // Re and Im of Q = e^{-i 2 pi lambda_alpha} 
            REAL4 realXP, imagXP;                                               // Re/Im of sum_k X_ak * P_ak
            REAL4 realQXP, imagQXP;                                             // Re/Im of Q_alpha R_alpha
            REAL4 lambda_alpha;
            INT4 kstar;                                                         // central frequency-bin k* = round(xhat_alpha)
            REAL4 kappa_max, kappa_star;
            INT4 k0, k1;

            REAL4 _lambda_alpha_;
            // ----- calculate kappa_max and lambda_alpha 
            {
                // 1st oder: s = 0 
                REAL4 TdotM1 = *tSSB_TdotM1;                                    // defined as Tdot_al - 1 
                REAL4 dT = *tSSB_DeltaT_rem;

                REAL4 T0 = *tSSB_DeltaT_int;
                REAL4 deltaT = T0 + dT;

                // phi_alpha = f * Tas; 
                REAL4 T0rem = FMOD( T0, tau );

//                REAL4 phi_alpha_rem = phi_alpha_offs + f0 * dT + T0 * df + df * dT;  // shouldn't it be "phi_alpha_offs + f0 * dT + T0rem * df + df * dT" ?
                REAL4 phi_alpha_rem = f0 * dT;
                phi_alpha_rem += T0rem * df;
                phi_alpha_rem += df * dT;

                REAL4 Dphi_alpha_int = f0;
                REAL4 Dphi_alpha_rem = df + Freq * TdotM1;

                // higher-order spindowns 
                REAL4 Tas = deltaT;
                UINT4 s;
                for (s=1; s <= spdnOrder; s++)
                {
                  REAL4 fsdot = fkdot4ex.s[s];
                  Dphi_alpha_rem += fsdot * Tas * invfact[s]; 	/* here: Tas = DT^s */
                  Tas *= deltaT;					/* now:  Tas = DT^(s+1) */
                  phi_alpha_rem += fsdot * Tas * invfact[s+1];
                } /* for s <= spdnOrder */

                // Step 3: apply global factor of Tsft to complete Dphi_alpha 
                Dphi_alpha_int *= Tsft;
                Dphi_alpha_rem *= Tsft;

                REAL4 tmp = REMAINDER( 0.5f * Dphi_alpha_int ) + REMAINDER ( 0.5f * Dphi_alpha_rem );
                lambda_alpha = phi_alpha_rem - tmp;

                // real- and imaginary part of e^{-i 2 pi lambda_alpha } 
                SINCOS_2PI_TRIMMED(&imagQ, &realQ, - lambda_alpha);


                kstar = REAL4TOINT4(Dphi_alpha_int) + REAL4TOINT4(Dphi_alpha_rem);
                kappa_star = REMAINDER(Dphi_alpha_int) + REMAINDER(Dphi_alpha_rem);
                kappa_max = kappa_star + 1.0f * DTERMS - 1.0f;

                // ----- check that required frequency-bins are found in the SFTs ----- 
                k0 = kstar - DTERMS + 1;
                k1 = k0 + 2 * DTERMS - 1;
                if ((k0 < freqIndex0) || (k1 > freqIndex1))
                {
                    REAL4 NaN;
                    *((int*)&NaN) = NAN_UINT4;
                    FaFb_components[maxSfts * curIFO + curSfts].Fa.re = NaN;
                }

            } // compute kappa_star, lambda_alpha 

            // NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
            // the trig-functions need to be calculated only once!
            SINCOS_2PI_TRIMMED(&s_alpha, &c_alpha, kappa_star);
            c_alpha -= 1.0f;

            // We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
            // closest to zero and will pose no numerical difficulties !

            // ---------- calculate the (truncated to DTERMS) sum over k ---------- 
            Xalpha_l = Xalpha + k0 - freqIndex0;  // first frequency-bin in sum 

            realXP = 0;
            imagXP = 0;

            // if no danger of denominator -> 0 
            if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
            {
                // improved hotloop algorithm by Fekete Akos:
                // take out repeated divisions into a single common denominator,
                // plus use extra cleverness to compute the nominator efficiently...
                
                REAL4 Sn = (*Xalpha_l).re;
                REAL4 Tn = (*Xalpha_l).im;
                REAL4 pn = kappa_max;
                REAL4 qn = pn;
                REAL4 U_alpha, V_alpha;

                // recursion with 2*DTERMS steps 
                UINT4 l;
                for (l = 1; l < 2*DTERMS; l ++ )
                {
                    Xalpha_l ++;
                    pn = pn - 1.0f;                                            // p_(n+1) 
                    Sn = pn * Sn + qn * (*Xalpha_l).re;                         // S_(n+1) 
                    Tn = pn * Tn + qn * (*Xalpha_l).im;                         // T_(n+1) 
                    qn *= pn;                                                  // q_(n+1) 
                } // for l <= 2*DTERMS 

                REAL4 qn_inv = RECIPROCAL(qn);
                U_alpha = Sn * qn_inv;
                V_alpha = Tn * qn_inv;

                realXP = s_alpha * U_alpha - c_alpha * V_alpha;
                imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

            } // if |remainder| > LD_SMALL4 
            else
            { // otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar}
                UINT4 ind0;
                if ( kappa_star <= LD_SMALL4 ) ind0 = DTERMS - 1;
                else ind0 = DTERMS;
                realXP = TWOPI_F * Xalpha_l[ind0].re;
                imagXP = TWOPI_F * Xalpha_l[ind0].im;
            } // if |remainder| <= LD_SMALL4

            realQXP = realQ * realXP - imagQ * imagXP;
            imagQXP = realQ * imagXP + imagQ * realXP;

            // we're done: ==> combine these into Fa and Fb 
            a_alpha = (*amcoe_a);
            b_alpha = (*amcoe_b);

            Fa.re = a_alpha * realQXP;
            Fa.im = a_alpha * imagQXP;

            Fb.re = b_alpha * realQXP;
            Fb.im = b_alpha * imagQXP;

        }

    } //endif if (curSfts < sfts_length)

    offset = curIFO * maxSfts + curSfts;

    // store intermediate result in shared memory
    FaFb_components[offset].Fa.re = Fa.re;
    FaFb_components[offset].Fa.im = Fa.im;
    FaFb_components[offset].Fb.re = Fb.re;
    FaFb_components[offset].Fb.im = Fb.im;

#if !(USE_OPENCL_KERNEL_CPU)
    // Perform reduction. Assuming blockDim = [64,2,1]
    int m, o2;
    m = maxSfts;
    o2 = m/2;
    for (; m; m>>=1, o2>>=1) {
       barrier(CLK_LOCAL_MEM_FENCE);
       if (curSfts < o2)
       {
           FaFb_components[offset].Fa.re += FaFb_components[offset + o2].Fa.re;
           FaFb_components[offset].Fa.im += FaFb_components[offset + o2].Fa.im;
           FaFb_components[offset].Fb.re += FaFb_components[offset + o2].Fb.re;
           FaFb_components[offset].Fb.im += FaFb_components[offset + o2].Fb.im;
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    if (curSfts == 0 && curIFO == 0)
    {
        // Compute final FaFb value (implicit reduction over two detectors)
        Fa.re = (FaFb_components[0].Fa.re + FaFb_components[maxSfts].Fa.re) * INV_TWOPI_F;
        Fa.im = (FaFb_components[0].Fa.im + FaFb_components[maxSfts].Fa.im) * INV_TWOPI_F;
        Fb.re = (FaFb_components[0].Fb.re + FaFb_components[maxSfts].Fb.re) * INV_TWOPI_F;
        Fb.im = (FaFb_components[0].Fb.im + FaFb_components[maxSfts].Fb.im) * INV_TWOPI_F;

        // compute final F-statistic value

        // convenient shortcuts
        REAL4 Ad = ABCInvDd[curSegment].Ad;
        REAL4 Bd = ABCInvDd[curSegment].Bd;
        REAL4 Cd = ABCInvDd[curSegment].Cd;
        REAL4 Dd_inv = ABCInvDd[curSegment].InvDd;

        // NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
        // therefore there is a factor of 2 difference with respect to the equations in JKS, which
        // where based on the single-sided PSD.
        //
        Fstat[ curBin * numSegs + curSegment ] = Dd_inv * ( Bd * (SQ(Fa.re) + SQ(Fa.im) )
                                                + Ad * ( SQ(Fb.re) + SQ(Fb.im) )
                                                - 2.0f * Cd *( Fa.re * Fb.re + Fa.im * Fb.im )
                                                );
    }
    
#endif // if !(USE_OPENCL_KERNEL_CPU)

}

