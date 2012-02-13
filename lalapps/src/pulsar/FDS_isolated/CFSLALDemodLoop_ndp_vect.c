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

	    REAL4 tFreq[4]; /* tempFreq1 as a vector */
	    REAL4 aFreq[4]; /* accFreq   as a vector */
	    REAL4  Xsum[4]; /* vector holding partial sums */

	    /* Unfortunately AltiVec / Velocity Engine doesn't do double precision
	       vectors, so we have to distinguish here again for the type of
	       double precision calculation (not really gaining readability...)
	    */
#ifdef USE_DOUBLE_VECT
	    /* the same in double precision, two elements (one complex) at once */
	    REAL8  trans[2]; /* for type conversion REAL4->REAL8 */
	    REAL8 tFreqD[2]; /* tempFreq1 as a vector */
	    REAL8 aFreqD[2]; /* accFreq   as a vector */
#else
	    REAL8 tFreqD, aFreqD;
#endif
	    REAL8  XsumD[2]; /* vector holding partial sums */

	    /* init vectors */
	    aFreq[0] = 1.0; aFreq[1] = 1.0;
	    aFreq[2] = 1.0; aFreq[3] = 1.0;

	    tFreq[0] = tempFreq0 + klim/2 - 1; tFreq[1] = tFreq[0];
	    tFreq[2] = tempFreq0 + klim/2 - 2; tFreq[3] = tFreq[2];

	    Xsum[0] = 0.0; Xsum[1] = 0.0;
	    Xsum[2] = 0.0; Xsum[3] = 0.0;

#ifdef USE_DOUBLE_VECT
	    aFreqD[0] = 1.0; aFreqD[1] = 1.0;
#else
	    aFreqD = 1.0;
#endif
	    XsumD[0]  = 0.0;  XsumD[1] = 0.0;

	    /* Vectorized version of the "hot loop" */
	    /* This loop has now been unrolled manually */

            /* for(k=0; k < klim / 2; k++) { */
	    {
	      UINT4 ve; /* this var should be optimized away... */

	      /* single precision vector loop */
#define VEC_LOOP(n)\
	      for(ve=0;ve<4;ve++) {\
                Xsum[ve] = Xsum[ve] * tFreq[ve] + (Xalpha_kR4+n)[ve] * aFreq[ve];\
		aFreq[ve] *= tFreq[ve];\
                tFreq[ve] -= 2.0;\
	      }

	      /* double precision vector loop */
#ifdef USE_DOUBLE_VECT
#define VEC_LOOP_D(n)\
              trans[0] = (Xalpha_kR4+n)[0];\
              trans[1] = (Xalpha_kR4+n)[1];\
	      for(ve=0;ve<2;ve++) {\
                XsumD[ve] = XsumD[ve] * tFreqD[ve] + trans[ve] * aFreqD[ve];\
		aFreqD[ve] *= tFreqD[ve];\
                tFreqD[ve] -= 1.0;\
	      }
#else
	      /* non-vectorizing double-precision loop */
#define VEC_LOOP_D(n)\
              XsumD[0] = XsumD[0] * tFreqD + aFreqD * Xalpha_kR4[n];\
              XsumD[1] = XsumD[1] * tFreqD + aFreqD * Xalpha_kR4[n+1];\
	      aFreqD *= tFreqD;\
              tFreqD -= 1.0;
#endif

	      VEC_LOOP(00+0); VEC_LOOP(00+4); VEC_LOOP(00+8); VEC_LOOP(00+12); 
	      VEC_LOOP(16+0); VEC_LOOP(16+4); VEC_LOOP(16+8);

              /* calculating the inner elements
		 VEC_LOOP(16+12); VEC_LOOP(32+0);
		 in double precision */

	      /* skip the values in single precision calculation */

	      for(ve=0;ve<4;ve++)
		tFreq[ve] -= 4.0; /* skip 4 elements */

	      /* VEC_LOOP(16+12); VEC_LOOP(32+0); */

#ifdef USE_DOUBLE_VECT
	      tFreqD[0] = tempFreq0 + klim/2 - 15; /* start at the 14th element */
	      tFreqD[1] = tFreqD[0];
#else
	      tFreqD = tempFreq0 + klim/2 - 15; /* start at the 14th element */
#endif
	      /* double precision vectorization */
	      VEC_LOOP_D(16+12);
	      VEC_LOOP_D(16+14);
	      VEC_LOOP_D(32+0);
	      VEC_LOOP_D(32+2);

                              VEC_LOOP(32+4); VEC_LOOP(32+8); VEC_LOOP(32+12); 
	      VEC_LOOP(48+0); VEC_LOOP(48+4); VEC_LOOP(48+8); VEC_LOOP(48+12); 
	    }

	    
	    /* conbination of three partial sums:
	       Xsum = (X1*a2*a3 + a1*X2*a3 + a1*a2*X3) / (a1*a2*a3)
	    */
#ifdef USE_DOUBLE_VECT
	    combAF1 =            aFreq[2] * aFreqD[0];
	    combAF2 = aFreq[0]            * aFreqD[0];
#else
	    combAF1 =            aFreq[2] * aFreqD;
	    combAF2 = aFreq[0]            * aFreqD;
#endif
	    combAF3 = aFreq[0] * aFreq[2];
	    combAF  = 1.0 / (aFreq[0] * combAF1);

	    XRes = (Xsum[0] * combAF1 + Xsum[2] * combAF2 + XsumD[0] * combAF3) * combAF;
	    XIms = (Xsum[1] * combAF1 + Xsum[3] * combAF2 + XsumD[1] * combAF3) * combAF;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
