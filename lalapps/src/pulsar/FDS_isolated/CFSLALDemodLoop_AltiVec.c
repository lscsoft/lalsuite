/*
*  Copyright (C) 2007 Bernd Machenschalk
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

          {

	    /* THIS IS DANGEROUS!! It relies on current implementation of COMPLEX8 type!! */
	    REAL4 *Xalpha_kR4 = &(Xalpha[sftIndex].re);

	    /* temporary variables to prevent double calculations */
	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
	    REAL4 XRes, XIms;
	    REAL8 combAF;

	    /* The main idea of the vectorization is that

	       [0] of a vector holds the real part of an even-index element
	       [1] of a vector holds the imaginary part of an even-index element
	       [2] of a vector holds the real part of an odd-index element
	       [3] of a vector holds the imaginary part of an odd-index element

	       The calculations for the four vector elements are performed independently
	       possibly using vector-operations, and the sumsfor even- and odd-indexed
	       element are combined at the very end.

	       There was added a double-precision part that calculates the "inner"
	       elements of the loop (that contribute the most to the result) in
	       double precision. If vectorized, a complex calculation (i.e. two
	       real ones) should be performed on a 128Bit vector unit capable of double
	       precision, such as SSE or AltiVec
 	    */

	    float XsumS[4]  __attribute__ ((aligned (16))); /* aligned for vector output */
	    /* the following vectors actually become registers in the AVUnit */
	    vector unsigned char perm;     /* holds permutation pattern for unaligned memory access */
	    vector float load0, load1, load2, load3;  /* temp registers for unaligned memory access */
	    vector float load4, load5, load6, load7;
	    vector float fdval, reTFreq;              /* temporary variables */
	    vector float Xsum  = {0,0,0,0};           /* collects the sums */
	    vector float four2 = {2,2,2,2};           /* vector constants */
	    vector float four4 = {4,4,4,4};
	    vector float tFreq = { ((float)(tempFreq0 + klim/2 - 1)), /* tempFreq as vector */
				   ((float)(tempFreq0 + klim/2 - 1)),
				   ((float)(tempFreq0 + klim/2 - 2)),
				   ((float)(tempFreq0 + klim/2 - 2)) };

	    REAL8 tFreqD, aFreqD = 1;      /* tempFreq and accFreq for double precision */
	    REAL8 XsumD[2] = {0,0};        /* partial sums */

	    /* Vectorized version of the Kernel loop */
	    /* This loop has now been unrolled manually */
	    {
	      /* single precision vector "loop" (isn't actually a loop anymore) */
#define VEC_LOOP(n,a,b)\
	      perm    = vec_lvsl(0,(Xalpha_kR4+(n)));\
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

	      /* init the memory access */
              load0 = vec_ld(0,(Xalpha_kR4));

	      /* seven single-precision calculations first */
	      VEC_LOOP(00+0,0,1); VEC_LOOP(00+4,1,2); VEC_LOOP(00+8,2,3); VEC_LOOP(00+12,3,4); 
	      VEC_LOOP(16+0,4,5); VEC_LOOP(16+4,5,6); VEC_LOOP(16+8,6,7);

              /* calculating the inner elements
		 VEC_LOOP(16+12); VEC_LOOP(32+0);
		 in double precision */

	      /* skip these values in single precision calculation */
	      tFreq   = vec_sub(tFreq,four4);

	      tFreqD = tempFreq0 + klim/2 - 15; /* start at the 14th element */

	      /* double precision calculations */
	      VEC_LOOP_D(16+12); VEC_LOOP_D(16+14);
	      VEC_LOOP_D(32+0);  VEC_LOOP_D(32+2);

	      /* the rest is done in single precision again */
	      /* init the memory access as above */
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
