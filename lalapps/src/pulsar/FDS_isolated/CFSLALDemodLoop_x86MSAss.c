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
	    
	    accFreq = 1.0 / accFreq;
	    XRes *= accFreq;
	    XIms *= accFreq;

#define AD_FLOAT ".float "
#define AD_ASCII ".string "

	    /* The SSE part is coded in Assembler */
	    __asm __volatile
	      (
	       "mov %2,%%edx\n\t"
	       "movss %3,%%xmm0\n\t"

	       "jmp contcode\n"
	       AD_ASCII "\"$Id$\"\n"

	       "V0011:\n\t"
	       AD_FLOAT "0.0\n\t"
	       AD_FLOAT "0.0\n\t"
	       AD_FLOAT "1.0\n\t"
	       AD_FLOAT "1.0\n"
	       "V2222:\n\t"
	       AD_FLOAT "2.0\n\t"
	       AD_FLOAT "2.0\n\t"
	       AD_FLOAT "2.0\n\t"
	       AD_FLOAT "2.0\n"
	       "contcode:\n\t"

	       "MOVUPS XMM5,V0011\n\t"
	       "SHUFPS XMM0,XMM0,0\n\t" /* XMM0: f   f   f   f */
	       "SUBPS  XMM0,XMM5\n\t" /* XMM0: f-1 f-1 f   f */
	       "XORPS  XMM2,XMM2\n\t" /* XMM2 will collect the low-precision values */
	       "MOVUPS XMM5,V2222\n\t"

	       /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
	       /* one SSE register will consist 4 REAL4 values */
	       /* 4 REAL4 vaules = 2 ReIm pairs */

	       /* we want to do this 14 times */
#define ADD4SSE \
	       "MOVLPS XMM4,[EDX]\n\t"\
	       "MOVHPS XMM4,[EDX+8]\n\t"\
	       "RCPPS  XMM1,XMM0\n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */\
	       "SUBPS  XMM0,XMM5\n\t" /* XMM2: f-3 f-3 f-2 f-2 */\
	       "MULPS  XMM1,XMM4\n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
	       "ADD    EDX,16\n\t" /* Xalpha_kX = Xalpha_kX + 2; */\
	       "ADDPS  XMM2,XMM1\n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */

	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
 
	       "SUBPS  XMM0,XMM5\n\t" /* JUMP OVER FPU CALCULATED VALUES */
	       "ADD    EDX,32\n\t" /* Xalpha_kX = Xalpha_kX + 2; */
	       "SUBPS  XMM0,XMM5\n\t" /* JUMP OVER FPU CALCULATED VALUES */

	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
	       ADD4SSE
 
 	       /* XMM2 consists the low precision values */
	       /* XRes and Xims consist the high precision values */
 
	       "MOVHLPS XMM3,XMM2\n\t" /* XMM3: ? ? C_ImH C_ReH */
	       "ADDPS   XMM2,XMM3\n\t" /* XMM2: - - C_Im C_Re */

	       "MOVSS   %0,XMM2\n\t" /* SAVE Re part */
	       "SHUFPS  XMM2,XMM2,1\n\t" /* XMM0: f   f   f   f */
	       "MOVSS   %1,XMM2\n\t" /* SAVE Im part */
 
	       : "=m" (XResX), "=m" (XImsX)       /* output  (here: to memory)*/
	       : "m" (Xalpha_kX), "m" (tempFreqX) /* input (here: from memory)*/
	       : "edx"                            /* clobbered registers */
	       );

	    /* And last, we add the single and double precision values */
	    
	    XRes = XRes + XResX;
	    XIms = XIms + XImsX;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
