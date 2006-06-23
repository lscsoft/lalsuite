          {
	    NRCSID (CFSLOOPVECTTAG, "$Id$");

	    /* THIS IS DANGEROUS!! It relies on current implementation of COMPLEX8 type!! */
	    REAL4 *Xalpha_kR4 = &(Xalpha[sftIndex].re);

	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
	    REAL4 XRes, XIms;

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

	    REAL4 tFreq[4]; /* tempFreq0 as a vector */
	    REAL4 tFint[4]; /* the integer part of the former tempFreq1 */
	    REAL4 aFreq[4]; /* accFreq   as a vector */
	    REAL4  Xsum[4]; /* vector holding partial sums */

	    /* init vectors */
	    aFreq[0] = 1.0;
	    aFreq[1] = 1.0;
	    aFreq[2] = 1.0;
	    aFreq[3] = 1.0;

	    tFreq[0] = tempFreq0;
	    tFreq[1] = tempFreq0;
	    tFreq[2] = tempFreq0;
	    tFreq[3] = tempFreq0;

	    tFint[0] = klim/2 - 1;
	    tFint[1] = klim/2 - 1;
	    tFint[2] = klim/2 - 2;
	    tFint[3] = klim/2 - 2;

	    Xsum[0] = 0.0;
	    Xsum[1] = 0.0;
	    Xsum[2] = 0.0;
	    Xsum[3] = 0.0;

	    /* Vectorized version of the "hot loop" */
	    /* This loop is now unrolled manually */
            /* for(k=0; k < klim / 2; k++) { */
	    {
	      UINT4 ve;
#define VEC_LOOP(n)\
	      for(ve=0;ve<4;ve++) {\
                Xsum[ve] = Xsum[ve] * (tFreq[ve] + tFint[ve]) + (Xalpha_kR4+n)[ve] * aFreq[ve];\
		aFreq[ve] *= (tFreq[ve] + tFint[ve]);\
                tFint[ve] -= 2.0;\
	      }

	      VEC_LOOP(00+0); VEC_LOOP(00+4); VEC_LOOP(00+8); VEC_LOOP(00+12); 
	      VEC_LOOP(16+0); VEC_LOOP(16+4); VEC_LOOP(16+8); VEC_LOOP(16+12); 
	      VEC_LOOP(32+0); VEC_LOOP(32+4); VEC_LOOP(32+8); VEC_LOOP(32+12); 
	      VEC_LOOP(48+0); VEC_LOOP(48+4); VEC_LOOP(48+8); VEC_LOOP(48+12); 
	    }

	    /* conbination:
	       Xs / aF = ( Xs_even * aF_odd + Xs_odd * aF_even ) / ( aF_even * aF_odd )
	    */
	    aFreq[1] *= aFreq[2];
	    XRes = (Xsum[0] * aFreq[2] + Xsum[2] * aFreq[0]) / aFreq[1];
	    XIms = (Xsum[1] * aFreq[2] + Xsum[3] * aFreq[0]) / aFreq[1];

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
