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

	    float XsumS[4]  __attribute__ ((aligned (16))); /* for output */
	    float align[4]  __attribute__ ((aligned (16))); /* for output */
	    /* vector unsigned char permute_v = vec_lvsl( 0, (float *) Xalpha_k ); */
	    vector float fdval, reTFreq;
	    vector float Xsum = {0,0,0,0};
	    vector float four2 = {2,2,2,2};
	    vector float four4 = {4,4,4,4};
	    vector float tFreq = { ((float)(tempFreq0 + klim/2 - 1)),
				   ((float)(tempFreq0 + klim/2 - 1)),
				   ((float)(tempFreq0 + klim/2 - 2)),
				   ((float)(tempFreq0 + klim/2 - 2)) };

	    /* Unfortunately AltiVec / Velocity Engine doesn't do double precision
	       vectors, so we have to distinguish here again for the type of
	       double precision calculation (not really gaining readability...)
	    */
	    REAL8 tFreqD, aFreqD = 1;
	    REAL8 XsumD[2] = {0,0}; /* vector holding partial sums */

	    /* 
	       tFreq[0] = tempFreq0 + klim/2 - 1; tFreq[1] = tFreq[0];
	       tFreq[2] = tempFreq0 + klim/2 - 2; tFreq[3] = tFreq[2];

	       aFreqD = 1.0;
	       XsumD[0]  = 0.0;  XsumD[1] = 0.0;
	    */

	    /* Vectorized version of the Kernel loop */
	    /* This loop has now been unrolled manually */

            /* for(k=0; k < klim / 2; k++) { */
	    {
	      UINT4 ve; /* this var should be optimized away... */

	      /* single precision vector loop */
#define VEC_LOOP(n)\
              for(ve=0;ve<4;ve++)\
	        align[ve]=(Xalpha_kR4+(n))[ve];\
	      fdval   = vec_ld(0,align);\
	      reTFreq = vec_re(tFreq);\
	      tFreq   = vec_sub(tFreq,four2);\
	      Xsum    = vec_madd(fdval, reTFreq, Xsum);

#if 0
#define ADD4SSEA(a,b) \
                  "movlps "#a"(%[Xalpha_kX]),%%xmm4 \n\t"\
		  "movhps "#b"(%[Xalpha_kX]),%%xmm4 \n\t"\
		  "rcpps %%xmm0,%%xmm1       \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */\
                  "subps %%xmm5,%%xmm0       \n\t" /* XMM0: f-3 f-3 f-2 f-2 */\
                  "mulps %%xmm4,%%xmm1       \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
                  "addps %%xmm1,%%xmm2       \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */

	      mm0 = tFreq;
	      mm1 = reTF;
	      mm2 = Xsum;
	      mm4 = tmp;
	      mm5 = (2 2 2 2);
#endif

	      /* non-vectorizing double-precision "loop" */
#define VEC_LOOP_D(n)\
              XsumD[0] = XsumD[0] * tFreqD + aFreqD * Xalpha_kR4[n];\
              XsumD[1] = XsumD[1] * tFreqD + aFreqD * Xalpha_kR4[n+1];\
	      aFreqD *= tFreqD;\
              tFreqD -= 1.0;

	      VEC_LOOP(00+0); VEC_LOOP(00+4); VEC_LOOP(00+8); VEC_LOOP(00+12); 
	      VEC_LOOP(16+0); VEC_LOOP(16+4); VEC_LOOP(16+8);

              /* calculating the inner elements
		 VEC_LOOP(16+12); VEC_LOOP(32+0);
		 in double precision */

	      /* skip the values in single precision calculation */
	      tFreq   = vec_sub(tFreq,four4);

	      tFreqD = tempFreq0 + klim/2 - 15; /* start at the 14th element */

	      /* double precision vectorization */
	      VEC_LOOP_D(16+12);
	      VEC_LOOP_D(16+14);
	      VEC_LOOP_D(32+0);
	      VEC_LOOP_D(32+2);

                              VEC_LOOP(32+4); VEC_LOOP(32+8); VEC_LOOP(32+12); 
	      VEC_LOOP(48+0); VEC_LOOP(48+4); VEC_LOOP(48+8); VEC_LOOP(48+12); 
	    }

	    /* output the vector */
	    vec_st(Xsum,0,XsumS);

	    /* conbination of three partial sums: */
	    combAF  = 1.0 / aFreqD;
	    XRes = XsumD[0] * combAF + XsumS[0] + XsumS[2];
	    XIms = XsumD[1] * combAF + XsumS[1] + XsumS[3];

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
