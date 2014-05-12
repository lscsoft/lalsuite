// this function definition 'template' requires 2 macros to be set:
// FUNC: the function name
// HOTLOOP_SOURCE: the filename to be included containing the hotloop source

/* Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
static int
FUNC ( Fcomponents *FaFb,                    /* [out] Fa,Fb (and possibly atoms) returned */
       const SFTVector *sfts,                /* [in] input SFTs */
       const PulsarSpins fkdot,              /* [in] frequency and derivatives fkdot = d^kf/dt^k */
       const SSBtimes *tSSB,                 /* [in] SSB timing series for particular sky-direction */
       const AMCoeffs *amcoe,                /* [in] antenna-pattern coefficients for this sky-direction */
       const ComputeFParams *params )        /* addition computational params */
{
  /* ----- check validity of input */
  const UINT4 UNUSED Dterms = params->Dterms;
  RUNTIME_CHECK

  if ( !FaFb ) {
    XLALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      XLALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > LAL_FACT_MAX )
    {
      XLALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
                     LAL_FACT_MAX, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( XLAL_EINVAL);
    }

  UINT4 alpha;                  /* loop index over SFTs */
  UINT4 spdnOrder;              /* maximal spindown-orders */
  UINT4 numSFTs;                /* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 Tsft;                   /* length of SFTs in seconds */
  INT4 freqIndex0;              /* index of first frequency-bin in SFTs */
  INT4 freqIndex1;              /* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;           /* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;   /* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;              /* SFT alpha  */

  const REAL8 norm = OOTWOPI;

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  // locally initialize sin/cos lookuptable, as some hotloops use that directly
  static int firstcall = 1;
  if ( firstcall ) {
    XLALSinCosLUTInit();
    firstcall = 0;
  }

  /* ----- prepare return of 'FstatAtoms' if requested */
  if ( params->returnAtoms )
    {
      if ( (FaFb->multiFstatAtoms = LALMalloc ( sizeof(*FaFb->multiFstatAtoms) )) == NULL ){
        XLAL_ERROR ( XLAL_ENOMEM );
      }
      FaFb->multiFstatAtoms->length = 1;        /* in this function: single-detector only */
      if ( (FaFb->multiFstatAtoms->data = LALMalloc ( 1 * sizeof( *FaFb->multiFstatAtoms->data) )) == NULL ){
        LALFree (FaFb->multiFstatAtoms);
        XLAL_ERROR ( XLAL_ENOMEM );
      }
      if ( (FaFb->multiFstatAtoms->data[0] = XLALCreateFstatAtomVector ( numSFTs )) == NULL ) {
        LALFree ( FaFb->multiFstatAtoms->data );
        LALFree ( FaFb->multiFstatAtoms );
        XLAL_ERROR( XLAL_ENOMEM );
      }

      FaFb->multiFstatAtoms->data[0]->TAtom = Tsft;     /* time-baseline of returned atoms is Tsft */

    } /* if returnAtoms */

  // ----- find highest non-zero spindown-entry ----------
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] != 0.0 )
      break;

  Fa = 0.0f;
  Fb = 0.0f;

  a_al = amcoe->a->data;        /* point to beginning of alpha-arrays */
  b_al = amcoe->b->data;
  DeltaT_al = tSSB->DeltaT->data;
  Tdot_al = tSSB->Tdot->data;
  SFT_al = sfts->data;

  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      REAL4 a_alpha, b_alpha;

      INT4 kstar;               /* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      REAL4 realQ, imagQ;       /* Re and Im of Q = e^{i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;     /* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;   /* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_star;

      /* ----- calculate lambda_alpha */
      {
        UINT4 s;                /* loop-index over spindown-order */
        REAL8 phi_alpha, Dphi_alpha, DT_al;
        REAL8 Tas;      /* temporary variable to calculate (DeltaT_alpha)^s */
        REAL8 TAS_invfact_s;

        /* init for s=0 */
        phi_alpha = 0.0;
        Dphi_alpha = 0.0;
        DT_al = (*DeltaT_al);
        Tas = 1.0;              /* DeltaT_alpha ^ 0 */
        TAS_invfact_s=1.0;    /* TAS / s! */

        for (s=0; s <= spdnOrder; s++)
          {
            REAL8 fsdot = fkdot[s];
            Dphi_alpha += fsdot * TAS_invfact_s;  /* here: DT^s/s! */
            Tas *= DT_al;                         /* now: DT^(s+1) */
            TAS_invfact_s= Tas * LAL_FACT_INV[s+1];
            phi_alpha += fsdot * TAS_invfact_s;
          } /* for s <= spdnOrder */

        /* Step 3: apply global factors to complete Dphi_alpha */
        Dphi_alpha *= Tsft * (*Tdot_al);                /* guaranteed > 0 ! */
        lambda_alpha = 0.5 * Dphi_alpha - phi_alpha;

        kstar = (INT4) (Dphi_alpha);    /* k* = floor(Dphi_alpha) for positive Dphi */
        kappa_star = Dphi_alpha - 1.0 * kstar;  /* remainder of Dphi_alpha: >= 0 ! */

        /* ----- check that required frequency-bins are found in the SFTs ----- */
        k0 = kstar - Dterms + 1;
        k1 = kstar + Dterms;

        if ( (k0 < freqIndex0) || (k1 > freqIndex1) ) {
          XLAL_ERROR ( XLAL_EDOM, "Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n"
                       "\t\t[Parameters: alpha:%d, Dphi_alpha:%e, Tsft:%e, *Tdot_al:%e]\n",
                       k0, k1, freqIndex0, freqIndex1, alpha, Dphi_alpha, Tsft, *Tdot_al );
        }

      } /* compute kappa_star, lambda_alpha */

      /* ---------- calculate the (truncated to Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be
       * executed many millions of times, so anything in here
       * has a HUGE impact on the whole performance of the code.
       *
       * DON'T touch *anything* in here unless you really know
       * what you're doing !!
       *------------------------------------------------------------
       */
      { // ----- start: hotloop ----------
        COMPLEX8 *Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */

        /* somehow the branch prediction of gcc-4.1.2 terribly failes
           with the current case distinction in the hot-loop,
           having a severe impact on runtime of the E@H Linux App.
           So let's allow to give gcc a hint which path has a higher probablility */
        if ( likely( (kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4) ) )
          { /* if |remainder| > LD_SMALL4, ie no danger of denominator -> 0 */

#include HOTLOOP_SOURCE

          } /* if  */
        else
          { /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
            UINT4 ind0;

            /* real- and imaginary part of e^{i 2 pi lambda_alpha } */
            XLALSinCos2PiLUT ( &imagQ, &realQ, lambda_alpha );

            if ( kappa_star <= LD_SMALL4 ) {
              ind0 = Dterms - 1;
            }
            else {
              ind0 = Dterms;
            }

            realXP = TWOPI_FLOAT * crealf(Xalpha_l[ind0]);
            imagXP = TWOPI_FLOAT * cimagf(Xalpha_l[ind0]);

          } /* if |remainder| <= LD_SMALL4 */

      } // ----- end: hotloop ----------

      /* real- and imaginary part of e^{i 2 pi lambda_alpha } */
      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      COMPLEX16 Fa_alpha = crect( a_alpha * realQXP, a_alpha * imagQXP );
      Fa += Fa_alpha;

      COMPLEX16 Fb_alpha = crect( b_alpha * realQXP, b_alpha * imagQXP );
      Fb += Fb_alpha;

      /* store per-SFT F-stat 'atoms' for transient-CW search */
      if ( params->returnAtoms )
        {
          FaFb->multiFstatAtoms->data[0]->data[alpha].timestamp = (UINT4)XLALGPSGetREAL8( &SFT_al->epoch );
          FaFb->multiFstatAtoms->data[0]->data[alpha].a2_alpha   = a_alpha * a_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].b2_alpha   = b_alpha * b_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].ab_alpha   = a_alpha * b_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].Fa_alpha   = norm * Fa_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].Fb_alpha   = norm * Fb_alpha;
        }

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa = norm * Fa;
  FaFb->Fb = norm * Fb;

  return XLAL_SUCCESS;

} // XLALComputeFaFb<VARIANT>()

#undef FUNC
#undef HOTLOOP_SOURCE
#undef RUNTIME_CHECK
