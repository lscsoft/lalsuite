/* Modified version of ComputeFaFb() based on Xavies trick:
 * need sufficiently oversampled SFTs and uses ZERO Dterms.
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
static int
XLALComputeFaFbXavie ( Fcomponents *FaFb,               /* [out] Fa,Fb (and possibly atoms) returned */
                       const SFTVector *sfts,           /* [in] input SFTs */
                       const PulsarSpins fkdot,         /* [in] frequency and derivatives fkdot = d^kf/dt^k */
                       const SSBtimes *tSSB,            /* [in] SSB timing series for particular sky-direction */
                       const AMCoeffs *amcoe,           /* [in] antenna-pattern coefficients for this sky-direction */
                       const ComputeFParams *params     /* additional computational params */
                       )
{
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

  REAL4 Upsampling;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
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
  if ( params->returnAtoms )
    {
      XLALPrintError ("%s: using the option 'returnAtoms' is not supported in this function!\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  Upsampling = (REAL4) params->upsampling;

  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
    freqIndex0 *= Upsampling;
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
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

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 Xalpha_l;        /* frequency-bin k in current SFT */
      REAL4 realQ, imagQ;       /* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realQXP, imagQXP;   /* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha;       /* !NOTE!: this MUST be REAL8!!! otherwise you lose the signal! */

      /* ----- calculate kappa_max and lambda_alpha */
      {
        UINT4 s;                /* loop-index over spindown-order */
        REAL8 phi_alpha, Dphi_alpha, DT_al;
        REAL8 Tas;      /* temporary variable to calculate (DeltaT_alpha)^s */

        /* init for s=0 */
        phi_alpha = 0.0;
        Dphi_alpha = 0.0;
        DT_al = (*DeltaT_al);
        Tas = 1.0;              /* DeltaT_alpha ^ 0 */

        for (s=0; s <= spdnOrder; s++)
          {
            REAL8 fsdot = fkdot[s];
            Dphi_alpha += fsdot * Tas * LAL_FACT_INV[s];        /* here: DT^s/s! */
            Tas *= DT_al;                               /* now: DT^(s+1) */
            phi_alpha += fsdot * Tas * LAL_FACT_INV[s+1];
          } /* for s <= spdnOrder */

        /* Step 3: apply global factors to complete Dphi_alpha */
        Dphi_alpha *= Tsft * (*Tdot_al);                /* guaranteed > 0 ! */

        lambda_alpha = phi_alpha - 0.5 * Dphi_alpha;

        /* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
        if ( XLALSinCos2PiLUT ( &imagQ, &realQ, - lambda_alpha ) ) {
          XLAL_ERROR (XLAL_EFUNC);
        }

        kstar = (INT4) (Dphi_alpha * Upsampling + 0.5f - freqIndex0);   /* k* = round(Dphi_alpha*chi) for positive Dphi */

        /* ----- check that required frequency-bins are found in the SFTs ----- */
        if ( (kstar < 0) || (kstar > freqIndex1 - freqIndex0) )
          {
            XLALPrintError ("Required frequency-bin [%d] not covered by SFT-interval [%d, %d]\n\n",
                           freqIndex0 + kstar, freqIndex0, freqIndex1 );
            XLAL_ERROR(XLAL_EDOM);
          }

      } /* compute kstar, lambda_alpha */

      /* ---------- calculate the (truncated to ZERO Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be
       * executed many millions of times, so anything in here
       * has a HUGE impact on the whole performance of the code.
       *
       * DON'T touch *anything* in here unless you really know
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_l = Xalpha[kstar];  /* frequency-bin to use */

      /* lim_{kappa_star->0}P_alpha,k  = 2pi delta_{k,kstar} */

      /* combine with e^-i 2pi lambda_alpha */
      realQXP = realQ * crealf(Xalpha_l) - imagQ * cimagf(Xalpha_l);
      imagQXP = realQ * cimagf(Xalpha_l) + imagQ * crealf(Xalpha_l);

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa += crect( a_alpha * realQXP, a_alpha * imagQXP );
      Fb += crect( b_alpha * realQXP, b_alpha * imagQXP );

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa = Fa;
  FaFb->Fb = Fb;

  return XLAL_SUCCESS;

} /* XLALComputeFaFbXavie() */
