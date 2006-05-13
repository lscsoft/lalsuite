          {
	    /* Options */
#define USEFloatFreQ1 0
#define SetValuesAtEnd 0
#define AddXResInCCode 0
	    /* when Xres + XResX shall be added in the C-Code */
	    /* there is a small difference in values */
	    /* C-Code:   -398873.94062822102569043636 */
	    /*                        ............... */
	    /* Asm-Code: -398873.94062570674577727914 */

	    
#ifdef __APPLE__ /* specials for Apples assembler */
#define AD_FLOAT ".single "
#define AD_ASCII ".ascii "
#else /* standard x86 gas */
#define AD_FLOAT ".float "
#define AD_ASCII ".string "
#endif
	    
#define reOffset "0"
#define imOffset "4"

            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;

	    /* prepare values for SSE */
	    REAL4 tempFreqX = tempFreq1; /* REAL4 because of SSE */
	    /* if this statement in the loop "accFreq *= (REAL4)tempFreq1;" */
	    /* shall use floats, the define USEFloatFreQ1 with 1 */
#if USEFloatFreQ1
	    REAL4 floatFreq1; /* REAL4 */
#endif
	    COMPLEX8 *Xalpha_kX = Xalpha_k; /* -> SSE values */

	    REAL4 XRes, XIms;   /* sums of Xa.re and Xa.im */
#if AddXResInCCode == 1
	    REAL4 XResX, XImsX; /* sums of Xa.re and Xa.im for SSE */
#endif

	    Xalpha_k = Xalpha_k + 14; /* -> FPU values */


	    /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
	    /* one SSE register will consist 4 REAL4 values */
	    /* 4 REAL4 vaules = 2 ReIm pairs */

	    /* we want to do this 14 times */
#define ADD4SSE(a,b) \
	       "movlps "#a"(%%edx),%%xmm4 \n\t"\
	       "movhps "#b"(%%edx),%%xmm4 \n\t"\
	       "rcpps %%xmm0,%%xmm1       \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */\
	       "subps %%xmm5,%%xmm0       \n\t" /* XMM2: f-3 f-3 f-2 f-2 */\
	       "mulps %%xmm4,%%xmm1       \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
	       "addps %%xmm1,%%xmm2       \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */

	    /* x87 FPU part */
	    __asm __volatile
	      (
	       // "addl $14,%[Xalpha_k]     \n\t"
	       /* REAL4 XRes=0.0, XIms=0.0;   /* sums of Xa.re and Xa.im */
	       "fld1                     \n\t" /* 1 */
	       "fldz                     \n\t" /* XRes, 1 */
	       "mov %[Xalpha_kX],%%edx   \n\t"
	       "fld %%st(0)              \n\t" /* XIms, XRes, 1 */

	       /* REAL4 accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */
	       "fld %%st(2)              \n\t" /* accFreq, XIms, XRes, 1 */

	       /* tempFreq1 = tempFreq1 - 14; */
	       ".section .rodata         \n"
	       ".align 8                 \n\t"
	       "C14:                     \n\t"
	       AD_FLOAT "14.0            \n\t"
	       ".text                    \n\t"

	       "fldl %[tempFreq1]        \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "fsub C14                 \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
 
	       /* for(k=0; k < 4 ; k++) */
	       /* { */
	       /*   XRes = tempFreq1 * XRes + (*Xalpha_k).re * accFreq; */
	       /*   XIms = tempFreq1 * XIms + (*Xalpha_k).im * accFreq; */
	       "fmul %%st, %%st(3)       \n\t" /* tempFreq1, accFreq, XIms, XRes*tempFreq1, 1 */
	       "fmul %%st, %%st(2)       \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "reOffset" (%[Xalpha_k])     \n\t" /* (*Xalpha_k).re, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).re*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "faddp %%st, %%st(4)      \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "imOffset" (%[Xalpha_k])     \n\t" /* (*Xalpha_k).im, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).im*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "faddp %%st, %%st(3)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */

	       /*   accFreq *= (REAL4)tempFreq1; */
#if USEFloatFreQ1
	       "fsts %[floatFreq1]       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "flds %[floatFreq1]       \n\t" /* (REAL4)tempFreq1, tempFreq1, accFreq, XIms, XRes, 1 */
	       "fmulp %%st, %%st(2)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#else
	       "fmul %%st, %%st(1)       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#endif

	       /*   tempFreq1 --; */
	       "fsub %%st(4)             \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       /*   Xalpha_k ++; */
	       /* "add $8, %[Xalpha_k]      \n\t" done as offset */
	       /* } /* for k < klim */
	       "fmul %%st, %%st(3)       \n\t" /* tempFreq1, accFreq, XIms, XRes*tempFreq1, 1 */
	       "fmul %%st, %%st(2)       \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "reOffset"+8 (%[Xalpha_k])     \n\t" /* (*Xalpha_k).re, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).re*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "faddp %%st, %%st(4)      \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "imOffset"+8 (%[Xalpha_k])     \n\t" /* (*Xalpha_k).im, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).im*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "faddp %%st, %%st(3)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#if USEFloatFreQ1
	       "fsts %[floatFreq1]       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "flds %[floatFreq1]       \n\t" /* (REAL4)tempFreq1, tempFreq1, accFreq, XIms, XRes, 1 */
	       "fmulp %%st, %%st(2)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#else
	       "fmul %%st, %%st(1)       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#endif
	       "fsub %%st(4)             \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	    
	       "fmul %%st, %%st(3)       \n\t" /* tempFreq1, accFreq, XIms, XRes*tempFreq1, 1 */
	       "fmul %%st, %%st(2)       \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "reOffset"+16 (%[Xalpha_k])     \n\t" /* (*Xalpha_k).re, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).re*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "faddp %%st, %%st(4)      \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "imOffset"+16 (%[Xalpha_k])     \n\t" /* (*Xalpha_k).im, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).im*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "faddp %%st, %%st(3)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#if USEFloatFreQ1
	       "fsts %[floatFreq1]       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "flds %[floatFreq1]       \n\t" /* (REAL4)tempFreq1, tempFreq1, accFreq, XIms, XRes, 1 */
	       "fmulp %%st, %%st(2)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#else
	       "fmul %%st, %%st(1)       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#endif
	       "fsub %%st(4)             \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */

	       "fmul %%st, %%st(3)       \n\t" /* tempFreq1, accFreq, XIms, XRes*tempFreq1, 1 */
	       "fmul %%st, %%st(2)       \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "reOffset"+24 (%[Xalpha_k])     \n\t" /* (*Xalpha_k).re, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).re*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes*tempFreq1, 1 */
	       "faddp %%st, %%st(4)      \n\t" /* tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "flds "imOffset"+24 (%[Xalpha_k])     \n\t" /* (*Xalpha_k).im, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "fmul %%st(2), %%st       \n\t" /* (*Xalpha_k).im*accFreq, tempFreq1, accFreq, XIms*tempFreq1, XRes, 1 */
	       "faddp %%st, %%st(3)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#if SetValuesAtEnd
#if USEFloatFreQ1
	       "fsts %[floatFreq1]       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "flds %[floatFreq1]       \n\t" /* (REAL4)tempFreq1, tempFreq1, accFreq, XIms, XRes, 1 */
	       "fmulp %%st, %%st(2)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#else
	       "fmul %%st, %%st(1)       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#endif
	       "fsub %%st(4)             \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */

	       /* accFreq = 1.0 / accFreq; */
	       "fld %%st(4)              \n\t" /* 1, tempFreq1, accFreq, XIms, XRes, 1 */
	       "fdivp %%st(2)            \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "fstpl %[tempFreq1]       \n\t" /* accFreq, XIms, XRes, 1 */
#else
#if USEFloatFreQ1
	       "fsts %[floatFreq1]       \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
	       "flds %[floatFreq1]       \n\t" /* (REAL4)tempFreq1, tempFreq1, accFreq, XIms, XRes, 1 */
	       "fmulp %%st, %%st(2)      \n\t" /* tempFreq1, accFreq, XIms, XRes, 1 */
#else
	       "fmulp %%st, %%st(1)      \n\t" /* accFreq, XIms, XRes, 1 */
#endif
	       /* accFreq = 1.0 / accFreq; */
	       "fld %%st(3)              \n\t" /* 1, accFreq, XIms, XRes, 1 */
	       "fdivp %%st(1)            \n\t" /* accFreq, XIms, XRes, 1 */
#endif
	       /* XRes *= accFreq; */
	       "fmul %%st, %%st(2)       \n\t" /* accFreq, XIms, XRes, 1 */
	       /* XIms *= accFreq; */
	       "fmulp %%st, %%st(1)      \n\t" /* XIms, XRes, 1 */

	       "movss %[tempFreqX],%%xmm0\n\t"

	       "jmp cntcode              \n"
	       ".align 64                \n\t"
	       "V0011:                   \n\t"
	       AD_FLOAT "0.0             \n\t"
	       AD_FLOAT "0.0             \n\t"
	       AD_FLOAT "1.0             \n\t"
	       AD_FLOAT "1.0             \n"
	       "V2222:                   \n\t"
	       AD_FLOAT "2.0             \n\t"
	       AD_FLOAT "2.0             \n\t"
	       AD_FLOAT "2.0             \n\t"
	       AD_FLOAT "2.0             \n"
	       ".align 64                \n\t"
	       "cntcode:                 \n\t"

	       "shufps $0,%%xmm0,%%xmm0  \n\t" /* XMM0: f   f   f   f */
	       "subps  V0011,%%xmm0      \n\t" /* XMM0: f-1 f-1 f   f */
	       "xorps  %%xmm2,%%xmm2     \n\t" /* XMM2 will collect the low-precision values */
	       "movups V2222,%%xmm5      \n\t"

	       ADD4SSE(0,8)
	       ADD4SSE(16,24)
	       ADD4SSE(32,40)
	       ADD4SSE(48,56)
	       ADD4SSE(64,72)
	       ADD4SSE(80,88)
	       ADD4SSE(96,104) 

	       "subps %%xmm5,%%xmm0      \n\t" /* JUMP OVER FPU CALCULATED VALUES */
	       "addl  $144,%%edx         \n\t" /* Xalpha_kX = Xalpha_kX + 2; */
	       "subps %%xmm5,%%xmm0      \n\t" /* JUMP OVER FPU CALCULATED VALUES */

	       ADD4SSE(0,8)
	       ADD4SSE(16,24)
	       ADD4SSE(32,40)
	       ADD4SSE(48,56)
	       ADD4SSE(64,72)
	       ADD4SSE(80,88)
	       ADD4SSE(96,104) 

 
 	       /* XMM2 consists the low precision values */
	       /* XRes and Xims consist the high precision values */
 
	       "movhlps %%xmm2,%%xmm3    \n\t" /* XMM3: ? ? C_ImH C_ReH */
	       "addps %%xmm3,%%xmm2      \n\t" /* XMM2: - - C_Im C_Re */

	       "ffree %%st(2)            \n\t" /* XIms, XRes */
#if AddXResInCCode == 0
	       "movss %%xmm2,%[XRes]     \n\t" /* SAVE Re part */
	       "shufps $1,%%xmm2,%%xmm2  \n\t" /* XMM0: f   f   f   f */
	       "movss %%xmm2,%[XIms]     \n\t" /* SAVE Im part */

	       "fadds %[XIms]            \n\t" /* XIms, XRes */
#else
	       "movss %%xmm2,%[XResX]    \n\t" /* SAVE Re part */
	       "shufps $1,%%xmm2,%%xmm2  \n\t" /* XMM0: f   f   f   f */
	       "movss %%xmm2,%[XImsX]    \n\t" /* SAVE Im part */

#endif
	       "fstps %[XIms]            \n\t" /* XRes */
#if AddXResInCCode == 0
	       "fadds %[XRes]            \n\t" /* XRes */
#endif
	       "fstps %[XRes]            \n\t" /* FPU Stack is empty */
 
	       : /* output  (here: to memory)*/
#if AddXResInCCode == 1
	       [XResX] "=m" (XResX),
	       [XImsX] "=m" (XImsX),
#endif
	       [XRes] "=m" (XRes),
	       [XIms] "=m" (XIms)
	       : /* input */
#if USEFloatFreQ1
	       [floatFreq1] "m" (floatFreq1),  
#endif
	       [tempFreq1] "m" (tempFreq1),
	       [Xalpha_k]  "r" (Xalpha_k),
	       [Xalpha_kX] "m" (Xalpha_kX),
	       [tempFreqX] "m" (tempFreqX)
	       : /* clobbered registers */
	       "edx",
	       "st","st(1)","st(2)","st(3)","st(4)","st(5)","st(6)","st(7)",
	       "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"
	       );

#if AddXResInCCode
	    /* And last, we add the single and double precision values */
	    XRes = XRes + XResX;
	    XIms = XIms + XImsX;
#endif
            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          }
