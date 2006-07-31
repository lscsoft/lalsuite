          {
	    NRCSID (CFSLOOPVECTTAG, "$Id$");

	    /* THIS IS DANGEROUS!! It relies on current implementation of COMPLEX8 type!! */
	    REAL4 *Xalpha_kR4 = &(Xalpha[sftIndex].re);

	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
	    REAL4 XRes, XIms;
	    REAL8 combAF, combAF1, combAF2, combAF3;

	    /* The main idea of the vectorization is that

	       [0] of a vector holds the real part of an even-index element
	       [1] of a vector holds the imaginary part of an even-index element
	       [2] of a vector holds the real part of an odd-index element
	       [3] of a vector holds the imaginary part of an odd-index element

	       The calculations for the four vector elements are performed independently
	       possibly using vector-operations, and the sumsfor even- and odd-indexed
	       element are combined at the very end.

	       aFreq[0] = aFreq[1] and aFreq[2] = aFreq[3], as well as
	       tFreq[0] = tFreq[1] and tFreq[2] = tFreq[3] throughout the whole loop.

	       Note that compared to the "old" loop this one is only ran klim/2=DTerms
	       times.
 	    */

	    REAL4 tFreq[4]; /* tempFreq1 as a vector */
	    REAL4 aFreq[4]; /* accFreq   as a vector */
	    REAL4  Xsum[4]; /* vector holding partial sums */

	    /* the same in double precision, two elements at once */
	    REAL8 tFreqD[2]; /* tempFreq1 as a vector */
	    REAL8 aFreqD[2]; /* accFreq   as a vector */
	    REAL8  XsumD[2]; /* vector holding partial sums */

	    /* init vectors */
	    aFreq[0] = 1.0;
	    aFreq[1] = 1.0;
	    aFreq[2] = 1.0;
	    aFreq[3] = 1.0;

	    tFreq[0] = tempFreq0 + klim/2 - 1;
	    tFreq[1] = tempFreq0 + klim/2 - 1;
	    tFreq[2] = tempFreq0 + klim/2 - 2;
	    tFreq[3] = tempFreq0 + klim/2 - 2;

	    Xsum[0] = 0.0;
	    Xsum[1] = 0.0;
	    Xsum[2] = 0.0;
	    Xsum[3] = 0.0;

	    aFreqD[0] = 1.0;
	    aFreqD[1] = 1.0;
	    XsumD[0] = 0.0;
	    XsumD[1] = 0.0;

	    /* Vectorized version of the "hot loop" */
	    /* This loop is now unrolled manually */
            /* for(k=0; k < klim / 2; k++) { */
	    {
	      UINT4 ve;
#define VEC_LOOP(n)\
	      for(ve=0;ve<4;ve++) {\
                Xsum[ve] = Xsum[ve] * tFreq[ve] + (Xalpha_kR4+n)[ve] * aFreq[ve];\
		aFreq[ve] *= tFreq[ve];\
                tFreq[ve] -= 2.0;\
	      }

#define VEC_LOOP_D(n)\
	      for(ve=0;ve<2;ve++) {\
                XsumD[ve] = XsumD[ve] * tFreqD[ve] + (Xalpha_kR4+n)[ve] * aFreqD[ve];\
		aFreqD[ve] *= tFreqD[ve];\
                tFreqD[ve] -= 1.0;\
	      }

	      VEC_LOOP(00+0); VEC_LOOP(00+4); VEC_LOOP(00+8); VEC_LOOP(00+12); 
	      VEC_LOOP(16+0); VEC_LOOP(16+4); VEC_LOOP(16+8);

              /* calculating the inner elements
		 VEC_LOOP(16+12); VEC_LOOP(32+0);
		 in double precision */

	      /* skip the values in single precision calculation */

	      for(ve=0;ve<4;ve++)
		tFreq[ve] -= 4.0;

	      // VEC_LOOP(16+12); VEC_LOOP(32+0);

	      tFreqD[0] = tempFreq0 + klim/2 - 15;
	      tFreqD[1] = tFreqD[0];

	      /* double precision vectorization */
	      VEC_LOOP_D(16+12);
	      VEC_LOOP_D(16+14);
	      VEC_LOOP_D(32+0);
	      VEC_LOOP_D(32+2);

                              VEC_LOOP(32+4); VEC_LOOP(32+8); VEC_LOOP(32+12); 
	      VEC_LOOP(48+0); VEC_LOOP(48+4); VEC_LOOP(48+8); VEC_LOOP(48+12); 
	    }

	    

	    /* conbination:
	       Xs / aF = ( Xs_even * aF_odd + Xs_odd * aF_even ) / ( aF_even * aF_odd )
	    */
	    combAF1 =            aFreq[2] * aFreqD[0];
	    combAF2 = aFreq[0]            * aFreqD[0];
	    combAF3 = aFreq[0] * aFreq[2];
	    combAF  =             combAF3 * aFreqD[0];

	    XRes = (Xsum[0] * combAF1 + Xsum[2] * combAF2 + XsumD[0] * combAF3) / combAF;
	    XIms = (Xsum[1] * combAF1 + Xsum[3] * combAF2 + XsumD[1] * combAF3) / combAF;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
