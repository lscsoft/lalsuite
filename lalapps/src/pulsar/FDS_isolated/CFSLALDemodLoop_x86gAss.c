          {
            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;

	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;

	    REAL4 XRes=0.0, XIms=0.0;   /* sums of Xa.re and Xa.im */
	    REAL4 XResX=0.0, XImsX=0.0;   /* sums of Xa.re and Xa.im */

	    REAL4 accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */

	    /* prepare values for SSE */
	    REAL4 tempFreqX = tempFreq1; /* REAL4 because of SSE */
	    COMPLEX8 *Xalpha_kX = Xalpha_k; /* -> SSE values */

	    tempFreq1 = tempFreq1 - 14;
	    Xalpha_k = Xalpha_k + 14; /* -> FPU values */
 
	    /* let the compiler code the x87 part */
	    for(k=0; k < 4 ; k++)
	      {
		XRes = tempFreq1 * XRes + (*Xalpha_k).re * accFreq;
		XIms = tempFreq1 * XIms + (*Xalpha_k).im * accFreq;
		
		accFreq *= (REAL4)tempFreq1;
		tempFreq1 --;
		Xalpha_k ++;
	      } /* for k < klim */
	    
	    XRes /= accFreq;
	    XIms /= accFreq;

	    /* The SSE part is coded in Assembler */
	    __asm __volatile
	      (
	       "mov %2,%%edx             \n\t"
	       "movss %3,%%xmm0          \n\t"

	       "jmp contcode             \n"
	       ".string \"$Id$\"\n"

	       "V1100:                   \n\t"
	       ".float 0.0               \n\t"
	       ".float 0.0               \n\t"
	       ".float 1.0               \n\t"
	       ".float 1.0               \n"
	       "V2222:                   \n\t"
	       ".float 2.0               \n\t"
	       ".float 2.0               \n\t"
	       ".float 2.0               \n\t"
	       ".float 2.0               \n"
	       "V4444:                   \n\t"
	       ".float 4.0               \n\t"
	       ".float 4.0               \n\t"
	       ".float 4.0               \n\t"
	       ".float 4.0               \n"
	       "contcode:                \n\t"
	       "movups V1100,%%xmm5      \n\t"
	       "movups V2222,%%xmm6      \n\t"
	       "movups V4444,%%xmm7      \n\t"

	       "shufps $0,%%xmm0,%%xmm0  \n\t" /* XMM0: f   f   f   f */
	       "subps %%xmm5,%%xmm0      \n\t" /* XMM0: f-1 f-1 f   f */
	       "xorps %%xmm2,%%xmm2      \n\t" /* XMM2 will collect the low-precision values */

	       /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
	       /* one SSE register will consist 4 REAL4 values */
	       /* 4 REAL4 vaules = 2 ReIm pairs */
 
	       ".rept 7                  \n\t" /* for (i=0; i<7; i++) { */
	       "movlps (%%edx),%%xmm4    \n\t"
	       "movhps +8(%%edx),%%xmm4  \n\t"
	       "rcpps %%xmm0,%%xmm1      \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
	       "subps %%xmm6,%%xmm0      \n\t" /* XMM2: f-3 f-3 f-2 f-2 */
	       "mulps %%xmm4,%%xmm1      \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
	       "addl  $16,%%edx          \n\t" /* Xalpha_kX = Xalpha_kX + 2; */
	       "addps %%xmm1,%%xmm2      \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */
	       ".endr                    \n\t" /* } for */
 
	       "subps %%xmm7,%%xmm0      \n\t" /* JUMP OVER FPU CALCULATED VALUES */
	       "addl  $32,%%edx          \n\t" /* Xalpha_kX = Xalpha_kX + 2; */

	       ".rept 7                  \n\t" /* for (i=0; i<7; i++) { */
	       "movlps (%%edx),%%xmm4    \n\t"
	       "movhps +8(%%edx),%%xmm4  \n\t"
	       "rcpps %%xmm0,%%xmm1      \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
	       "subps %%xmm6,%%xmm0      \n\t" /* XMM2: f-3 f-3 f-2 f-2 */
	       "mulps %%xmm4,%%xmm1      \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
	       "addl  $16,%%edx          \n\t" /* Xalpha_kX = Xalpha_kX + 2; */
	       "addps %%xmm1,%%xmm2      \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */
	       ".endr                    \n\t" /* } for */
 
 	       /* XMM2 consists the low precision values */
	       /* XRes and Xims consist the high precision values */
 
	       "movhlps %%xmm2,%%xmm3    \n\t" /* XMM3: ? ? C_ImH C_ReH */
	       "addps %%xmm3,%%xmm2      \n\t" /* XMM2: - - C_Im C_Re */

	       "movss %%xmm2,%0          \n\t" /* SAVE Re part */
	       "shufps $1,%%xmm2,%%xmm2  \n\t" /* XMM0: f   f   f   f */
	       "movss %%xmm2,%1          \n\t" /* SAVE Im part */
 
	       : "=m" (XResX), "=m" (XImsX)       /* output  (here: to memory)*/
	       : "m" (Xalpha_kX), "m" (tempFreqX) /* input (here: from memory)*/
	       : "xmm0","xmm1","xmm2","xmm3","xmm4",
	         "xmm5","xmm6","xmm7","edx"       /* clobbered registers */
	       );

	    /* And last, we add the single and double precision values */
	    
	    XRes = XRes + XResX;
	    XIms = XIms + XImsX;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
