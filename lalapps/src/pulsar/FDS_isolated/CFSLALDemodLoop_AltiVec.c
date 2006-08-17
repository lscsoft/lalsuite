          {
	    NRCSID (CFSLOOPVECTTAG, "$Id$");

	    /* THIS IS DANGEROUS!! It relies on current implementation of COMPLEX8 type!! */
	    REAL4 *Xalpha_kR4 = &(Xalpha[sftIndex].re);

	    /* temporary variables to prevent double calculations */
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

	       There was added a double-precision part that calculates the "inner"
	       elements of the loop (that contribute the most to the result) in
	       double precision. If vectorized, a complex calculation (i.e. two
	       real ones) should be performed on a 128Bit vector unit capable of double
	       precision, such as SSE or AltiVec
 	    */

	    float XsumS[4]  __attribute__ ((aligned (16))); /* aligned output */
	    vector unsigned char perm; /* = vec_lvsl(0,(float*)Xalpha_k); */
	    vector float fdval, reTFreq;
	    vector float load0, load1, load2, load3;
	    vector float load4, load5, load6, load7;
	    vector float Xsum  = {0,0,0,0};
	    vector float four2 = {2,2,2,2};
	    vector float four4 = {4,4,4,4};
	    vector float tFreq = { ((float)(tempFreq0 + klim/2 - 1)),
				   ((float)(tempFreq0 + klim/2 - 1)),
				   ((float)(tempFreq0 + klim/2 - 2)),
				   ((float)(tempFreq0 + klim/2 - 2)) };

	    REAL8 tFreqD, aFreqD = 1;
	    REAL8 XsumD[2] = {0,0}; /* vector holding partial sums */

	    /* Vectorized version of the Kernel loop */
	    /* This loop has now been unrolled manually */
            /* for(k=0; k < klim / 2; k++) { */
	    {
	      /* single precision vector "loop" (isn't actually a loop anymore) */
#define VEC_LOOP(n,a,b)\
              load##b = vec_ld(0,(Xalpha_kR4+(n)+4));\
	      fdval   = vec_perm(load##a,load##b,perm);\
	      reTFreq = vec_re(tFreq);\
	      tFreq   = vec_sub(tFreq,four2);\
	      Xsum    = vec_madd(fdval, reTFreq, Xsum);

	      /* non-vectorizing double-precision "loop" */
#define VEC_LOOP_D(n)\
              XsumD[0] = XsumD[0] * tFreqD + aFreqD * Xalpha_kR4[n];\
              XsumD[1] = XsumD[1] * tFreqD + aFreqD * Xalpha_kR4[n+1];\
	      aFreqD *= tFreqD;\
              tFreqD -= 1.0;

              load0 = vec_ld(0,(Xalpha_kR4));
	      perm  = vec_lvsl(0,(Xalpha_kR4));

	      /* seven single-precision calculations first */
	      VEC_LOOP(00+0,0,1); VEC_LOOP(00+4,1,2); VEC_LOOP(00+8,2,3); VEC_LOOP(00+12,3,4); 
	      VEC_LOOP(16+0,4,5); VEC_LOOP(16+4,5,6); VEC_LOOP(16+8,6,7);

              /* calculating the inner elements
		 VEC_LOOP(16+12); VEC_LOOP(32+0);
		 in double precision */

	      /* skip the values in single precision calculation */
	      tFreq   = vec_sub(tFreq,four4);

	      tFreqD = tempFreq0 + klim/2 - 15; /* start at the 14th element */

	      /* double precision calculations */
	      VEC_LOOP_D(16+12); VEC_LOOP_D(16+14);
	      VEC_LOOP_D(32+0);  VEC_LOOP_D(32+2);

	      /* the rest is done in single precision again */
              load0 = vec_ld(0,(Xalpha_kR4+36));

	      VEC_LOOP(32+4,0,1); VEC_LOOP(32+8,1,2); VEC_LOOP(32+12,2,3); 
	      VEC_LOOP(48+0,3,4); VEC_LOOP(48+4,4,5); VEC_LOOP(48+8,5,6); VEC_LOOP(48+12,6,7); 
	    }

	    /* output the vector */
	    vec_st(Xsum,0,XsumS);

	    /* conbination of the three partial sums: */
	    combAF  = 1.0 / aFreqD;
	    XRes = XsumD[0] * combAF + XsumS[0] + XsumS[2];
	    XIms = XsumD[1] * combAF + XsumS[1] + XsumS[3];

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
