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

/* LALDemod variant with large assembler part
 * Authors see ComputeFStatistic.c

   This version was derived from CFSLALDemod_SSEgas.c rev. 1.6
   Here the sin/cos taylor expansion was parallized / vectorized for SSE2
   For this to work, some things had to be changed:
   - The LUT are now built from complex elements, carrying cos calculation
     values in the real, sin calculation values in the imaginary part
   - To make the calculations identical for sin and cos, the signs of the
     terms are now in the values in the tables, not in the calculation
   - The "special treatment" of the results (tcos =- 1; imagQ *= -1) have
     now gone into where these values are used, rather than in the expansion

                                                         Bernd Machenschalk */

/* gcc version and Apples assembler specials aren't needed currently, but kept for reference */
#ifdef __GNUC__
#ifndef __GNUC_PATCHLEVEL__
#define __GNUC_PATCHLEVEL__ 0
#endif
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

/* specials for Apples assembler */
#ifdef __APPLE__
#define AD_FLOAT   ".single "
#define AD_ASCII   ".ascii "
#define AD_ALIGN16 ".align 4"
#define AD_ALIGN64 ".align 6"
#else /* x86 gas */
#define AD_FLOAT   ".float "
#define AD_ASCII   ".string "
#define AD_ALIGN16 ".align 16"
#define AD_ALIGN64 ".align 64"
#endif

/* define a CPU_TYPE if there isn't one */
#ifndef CPU_TYPE_S
#define CPU_TYPE_S "0"
#endif

/* Lookup tables for fast sin/cos calculation */
/* we'll combine sin/cos values in complex items, so we'll define
   alialses for better readability */
#define cs re
#define sn im
static COMPLEX16 csVal[LUT_RES+1];
static COMPLEX16 csVal2PI[LUT_RES+1];
static COMPLEX16 csVal2PIPI[LUT_RES+1];
static REAL8 diVal[LUT_RES+1];
static BOOLEAN sincos_initialized = 1; /* reset after initializing the sin/cos tables */

void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 
{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8 *xSum=NULL, *ySum=NULL; /* temp variables for computation of fs*as and fs*bs */
  INT4 s;                       /* local variable for spinDwn calcs. */
  REAL8 xTemp;                  /* temp variable for phase model */
  REAL4 xTInt;                  /* integer part of xTemp */
  REAL8 deltaF;                 /* width of SFT band */
  INT4 k=0;                     /* loop counter */
  REAL8 *skyConst;              /* vector of sky constants data */
  REAL8 *spinDwn;               /* vector of spinDwn parameters (maybe a structure? */
  INT4  spOrder;                /* maximum spinDwn order */
  REAL8 realXP, imagXP;         /* temp variables used in computation of */
  INT4  nDeltaF;                /* number of frequency bins per SFT band */

  INT4  sftIndex;               /* more temp variables */
  COMPLEX16 Q __attribute__ ((aligned (16))); /* was REAL8 realQ, imagQ; */
  INT4 *tempInt1;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
  REAL8 f;

  UINT4 klim = 2*params->Dterms;
  /* #define klim 32               this actually is only a hint. The assembler Kernel loop
				   doesn't use this value at all, but relies on it being 32 */ 

  REAL4 tsin, tcos;
  COMPLEX16 tsincos __attribute__ ((aligned (16)));       /* tsin and tcos as complex pair */
  REAL4 realP, imagP;

  REAL8 A=params->amcoe->A;
  REAL8 B=params->amcoe->B;
  REAL8 C=params->amcoe->C;
  REAL8 D=params->amcoe->D;

  UINT4 M=params->SFTno;

  INITSTATUS(status);

  /* catch some obvious programming errors */
  ASSERT ( (Fs != NULL)&&(Fs->F != NULL), status, COMPUTEFSTAT_ENULL, COMPUTEFSTAT_MSGENULL );
  if (params->returnFaFb)
    {
      ASSERT ( (Fs->Fa != NULL)&&(Fs->Fb != NULL), status, COMPUTEFSTAT_ENULL, COMPUTEFSTAT_MSGENULL );
    }

  /* variable redefinitions for code readability */
  spOrder=params->spinDwnOrder;
  spinDwn=params->spinDwn;
  skyConst=params->skyConst;
  deltaF=(*input)->fft->deltaF;
  nDeltaF=(*input)->fft->data->length;

  /* res=10*(params->mCohSFT); */
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
  /* using three tables with values including PI is simply faster than doing
     the multiplications in the taylor expansion*/
  /* we include the signs of the expansion in the tables here to make the
     calculations for both sin and cos equal so we can do them in parallel */
  if ( sincos_initialized )
    {
      for (k=0; k <= LUT_RES; k++) {
	/* sin */
        csVal[k].sn      =   sin((LAL_TWOPI*k)/(LUT_RES));
        csVal2PI[k].cs   = - csVal[k].sn    * LAL_TWOPI;
        csVal2PIPI[k].sn =   csVal2PI[k].cs * LAL_PI;
	/* cos */
        csVal[k].cs      =   cos((LAL_TWOPI*k)/(LUT_RES));
        csVal2PI[k].sn   =   csVal[k].cs    * LAL_TWOPI;
        csVal2PIPI[k].cs = - csVal2PI[k].sn * LAL_PI;

	/* this additional table saves a "costly" division in sin/cos calculation */
	diVal[k] = (REAL8)k/(REAL8)(LUT_RES);
      }

      sincos_initialized = 0;
    }

  /* this loop computes the values of the phase model */
  xSum=(REAL8 *)LALMalloc(params->SFTno*sizeof(REAL8));
  ySum=(REAL8 *)LALMalloc(params->SFTno*sizeof(REAL8));
  tempInt1=(INT4 *)LALMalloc(params->SFTno*sizeof(INT4));
  for(alpha=0;alpha<params->SFTno;alpha++){
    tempInt1[alpha]=2*alpha*(spOrder+1)+1;
    xSum[alpha]=0.0;
    ySum[alpha]=0.0;
    for(s=0; s<spOrder;s++) {
      xSum[alpha] += spinDwn[s] * skyConst[tempInt1[alpha]+2+2*s];      
      ySum[alpha] += spinDwn[s] * skyConst[tempInt1[alpha]+1+2*s];
    }
  }


  /* Loop over frequencies to be demodulated */
  for(i=0 ; i< params->imax  ; i++ )
  {
    Fa.re =0.0;
    Fa.im =0.0;
    Fb.re =0.0;
    Fb.im =0.0;

    f = params->f0 + i*params->df;

    /* Loop over SFTs that contribute to F-stat for a given frequency */
    for(alpha=0;alpha<params->SFTno;alpha++)
      {
        REAL8 tempFreq0, tempFreq1;
        COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
        REAL4 a = params->amcoe->a->data[alpha];
        REAL4 b = params->amcoe->b->data[alpha];
        REAL8 x;
        REAL8 yTemp;

        /* NOTE: sky-constants are always positive!!
         * this can be seen from their definition (-> documentation)
         * we will use this fact in the following! 
         */
        xTemp = f * skyConst[ tempInt1[ alpha ] ] + xSum[ alpha ];       /* >= 0 !! */
        
        /* this will now be assumed positive, but we double-check this to be sure */
        if  (!finite(xTemp)) {
            fprintf (stderr, "xTemp is not finite\n");
            fprintf (stderr, "DEBUG: loop=%d, xTemp=%f, f=%f, alpha=%d, tempInt1[alpha]=%d\n", 
                     i, xTemp, f, alpha, tempInt1[alpha]);
            fprintf (stderr, "DEBUG: skyConst[ tempInt1[ alpha ] ] = %f, xSum[ alpha ]=%f\n",
                     skyConst[ tempInt1[ alpha ] ], xSum[ alpha ]);
            exit (COMPUTEFSTAT_EXIT_DEMOD);
        }
        if (xTemp < 0) {
            fprintf (stderr, "xTemp >= 0 failed\n");
            fprintf (stderr, "DEBUG: loop=%d, xTemp=%f, f=%f, alpha=%d, tempInt1[alpha]=%d\n", 
                     i, xTemp, f, alpha, tempInt1[alpha]);
            fprintf (stderr, "DEBUG: skyConst[ tempInt1[ alpha ] ] = %f, xSum[ alpha ]=%f\n",
                     skyConst[ tempInt1[ alpha ] ], xSum[ alpha ]);
            exit (COMPUTEFSTAT_EXIT_DEMOD);
        }

	{
#ifdef USE_ASS_XTEMP
	  INT4 xTInt;
	  __asm __volatile
#if (defined(USE_SSE2) && !(defined(USE_SSE3)))
	    (
	     "movsd     %[xTemp],%%xmm0 \n\t"
	     "cvttsd2si %%xmm0,%[xTInt] \n\t"
	     "cvtsi2sd  %[xTInt],%%xmm1 \n\t"
	     "subsd     %%xmm1,%%xmm0   \n\t"
	     "movsd     %%xmm0,%[tF0]   \n\t"

	     :
	     [xTInt] "=m" (xTInt),
	     [tF0]   "=m" (tempFreq0)
	     :
	     [xTemp] "m"  (xTemp)
#ifndef IGNORE_XMM_REGISTERS
	     : "xmm0","xmm1"
#endif
	     );
#else /* USE_SSE2 */
	    (
#ifdef USE_SSE3
	     "fld     %%st(0)          \n\t" /* xT xT */ /* clone xTemp on the stack */
	     "fisttpl %[xTInt]         \n\t" /* write integer w. truncation reagrdless of rounding */
	     "fisubl  %[xTInt]         \n\t" /* tempFreq0 = xTemp - xTInt */
	     "fstpl   %[tF0]           \n\t" /* store tempFreq0 */
#else
	     "fistl   %[xTInt]         \n\t" /* xT */           /* save the rounded value, keep the original in FPU */
	     "fildl   %[xTInt]         \n\t" /* xTi xT */       /* load it back on the stack as float */
	     "fxch                     \n\t" /* xT xTi */       /* flip to keep xTi after substraction */
	     "fsub    %%st(1),%%st(0)  \n\t" /* (xT-xTi) xTi */ /* substract */
	     "fsts    %[xTInt]         \n\t" /* (xT-xTi) xTi */ /* save the rounded value */
	     "sub     %%eax,%%eax      \n\t"                    /* EAX=0 */
	     "orl     %[xTInt],%%eax   \n\t"                    /* sets the S (sign) flag from xT-xTi */
	     "jns     tf0" CPU_TYPE_S "\n\t"                    /* jump if xTi<xT, i.e. rounding down */

	     "fld1                     \n\t" /* 1 (xT-xTi) xTi */     /* correct values if rounded up */
	     "fadd    %%st(0),%%st(1)  \n\t" /* 1 (xT-xTi+1) xTi */
	     "fsubrp  %%st(0),%%st(2)  \n\t" /* (xT-xTi+1) (xTi-1) */

	     "tf0" CPU_TYPE_S ":       \n\t"
	     "fstpl   %[tF0]           \n\t" /* xTi */          /* store tempFreq0 */
	     "fistpl  %[xTInt]         \n\t" /*% */             /* store xTInt */
#endif
	     :
	     [xTInt] "=m" (xTInt),
	     [tF0]   "=m" (tempFreq0)
	     :       "t"  (xTemp)
	     : "st(1)"
#ifndef USE_SSE3
	     ,"st(2)","eax"
#endif
#endif /* USE_SSE2 */
	     );
#else /* USE_ASS_XTEMP */

#ifdef USE_FLOOR
	  REAL4 xTInt =  floor(xTemp);
#else
	  REAL4 xTInt =  (UINT4)xTemp;
#endif
	  tempFreq0 = xTemp - xTInt;   /* lies in [0, +1) by definition */
#endif
	  sftIndex = xTInt - params->Dterms + 1 - params->ifmin;
	}

        if(sftIndex < 0){
              fprintf(stderr,"ERROR! sftIndex = %d < 0 in TestLALDemod run %d\n", sftIndex, cfsRunNo);
              fprintf(stderr," alpha=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                      alpha, xTemp, params->Dterms, params->ifmin);
              ABORT(status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
        }

        tempFreq1 = tempFreq0 + params->Dterms - 1;     /* positive if Dterms > 1 (trivial) */

        x = LAL_TWOPI * tempFreq1;      /* positive! */

        yTemp = f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ];

        /* we branch now (instead of inside the central loop)
         * depending on wether x can ever become SMALL in the loop or not, 
         * because it requires special treatment in the Dirichlet kernel
         */
        if ( tempFreq0 >= LD_SMALL ) {

          /* Variables for assembler Kernel loop */
          COMPLEX8 *Xalpha_kX = Xalpha + sftIndex;
          REAL4    tempFreqX  = tempFreq1;        /* REAL4 because of SSE */
          REAL4    XRes, XIms;                    /* sums of Xa.re and Xa.im */
          COMPLEX8 XSumsX;                        /* sums of Xa.re and Xa.im for SSE */
          COMPLEX8 *Xalpha_kF = Xalpha_kX + 16;   /* shift values for FPU calculation */
          REAL8    tempFreqF  = tempFreq1 - 15.0;
          
          /* constants in memory */
          static REAL8 lutr  = LUT_RES; /* LUT_RES in memory */
	  static REAL8 big   = 1<<30;   /* real(ly) big number */
	  static REAL4 half  = 0.5;
	  static REAL4 one   = 1.0;
          static REAL4 two   = 2.0;

          /* vector constants */
          /* having these not aligned will crash the assembler code */
          static REAL4 V0011[4] __attribute__ ((aligned (16))) = { 0,0,1,1 };
          static REAL4 V2222[4] __attribute__ ((aligned (16))) = { 2,2,2,2 };
          static REAL4 V4444[4] __attribute__ ((aligned (16))) = { 4,4,4,4 };
 
          static REAL8 yRem;

          __asm __volatile
      (
       /* calculation of tsin and tcos */

       "movsd    %[x],%%xmm0                    \n\t" /* xmm0L := x */
       "movsd    %[lutr],%%xmm1                 \n\t" /* xmm1L := LUT_RES */
       "mulsd    %%xmm0,%%xmm1                  \n\t" /* xmm1L := (x*LUT_RES) */
       "cvtsd2si %%xmm1,%%eax                   \n\t" /* EAX := (int)(x*LUT_RES) */
       "subsd    %[diVal](,%%eax,8),%%xmm0      \n\t" /* xmm0L := d = x - diVal[i] */
       "unpcklpd %%xmm0,%%xmm0                  \n\t" /* xmm0 := (d,d) */
       "add      %%eax,%%eax                    \n\t" /* EAX*2 as scale is 2,4,8 (need 16) */
       "movapd   %%xmm0,%%xmm1                  \n\t" /* xmm1 := (d,d) */
       "mulpd    %%xmm1,%%xmm1                  \n\t" /* xmm1 := (d*d,d*d) */
       "mulpd    %[csVal2PI]  (,%%eax,8),%%xmm0 \n\t"
       "mulpd    %[csVal2PIPI](,%%eax,8),%%xmm1 \n\t"
       "addpd    %[csVal]     (,%%eax,8),%%xmm0 \n\t"
       "addpd    %%xmm1,%%xmm0                  \n\t"
       "movapd   %%xmm0,%[cosv]                 \n\t"

       /* calculation of yRem */

       "movsd     %[yTemp],%%xmm0               \n\t" /* xmm0L := x */
       "addsd     %[big],%%xmm0                 \n\t" /* make sure xmm0L is >0 */
       "cvttsd2si %%xmm0,%%eax                  \n\t" /* EAX := (int)x */
       "cvtsi2sd  %%eax,%%xmm1                  \n\t" /* xmm1L := (int)x */
       "subsd     %%xmm1,%%xmm0                 \n\t" /* xmm0L := x-(int)x */

       /* calculation of realQ and imagQ */

       "movsd    %[lutr],%%xmm1                 \n\t" /* xmm1L := LUT_RES */
       "mulsd    %%xmm0,%%xmm1                  \n\t" /* xmm1L := (x*LUT_RES) */
       "cvtsd2si %%xmm1,%%eax                   \n\t" /* EAX := (int)(x*LUT_RES) */
       "subsd    %[diVal](,%%eax,8),%%xmm0      \n\t" /* xmm0L := d = x - diVal[i] */
       "unpcklpd %%xmm0,%%xmm0                  \n\t" /* xmm0 := (d,d) */
       "add      %%eax,%%eax                    \n\t" /* EAX*2 as scale is 2,4,8 (need 16) */
       "movapd   %%xmm0,%%xmm1                  \n\t" /* xmm1 := (d,d) */
       "mulpd    %%xmm1,%%xmm1                  \n\t" /* xmm1 := (d*d,d*d) */
       "mulpd    %[csVal2PI]  (,%%eax,8),%%xmm0 \n\t"
       "mulpd    %[csVal2PIPI](,%%eax,8),%%xmm1 \n\t"
       "addpd    %[csVal]     (,%%eax,8),%%xmm0 \n\t"
       "addpd    %%xmm1,%%xmm0                  \n\t"
       "movapd   %%xmm0,%[cosq]                 \n\t"

       /* Kernel loop */

       /* the "inner" four complex values that contribute the most to the sum
          will be calculated by the FPU with 80 Bit precision,
          the SSE taking care of the others ("wings") */

       /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
       /* one SSE register will consist 4 REAL4 values */
       /* 4 REAL4 vaules = 2 ReIm pairs */

       /* we want to do this SSE code 14 times */
#define ADD4SSEA(a,b) \
       "movlps "#a"(%[Xalpha_kX]),%%xmm4 \n\t"\
       "movhps "#b"(%[Xalpha_kX]),%%xmm4 \n\t"\
       "rcpps %%xmm0,%%xmm1       \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */\
       "subps %%xmm5,%%xmm0       \n\t" /* XMM0: f-3 f-3 f-2 f-2 */\
       "mulps %%xmm4,%%xmm1       \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
       "addps %%xmm1,%%xmm2       \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */

#define ADD4SSEB(a,b) \
       "movlps "#a"(%[Xalpha_kX]),%%xmm7 \n\t"\
       "movhps "#b"(%[Xalpha_kX]),%%xmm7 \n\t"\
       "rcpps %%xmm0,%%xmm6       \n\t" /* XMM6: 1/(f-1) 1/(f-1) 1/f 1/f */\
       "subps %%xmm5,%%xmm0       \n\t" /* XMM0: f-3 f-3 f-2 f-2 */\
       "mulps %%xmm7,%%xmm6       \n\t" /* XMM6: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
       "addps %%xmm6,%%xmm2       \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */

       /* SSE prelude */
       "movss   %[tempFreqX],%%xmm0   \n\t"
       "shufps  $0,%%xmm0,%%xmm0      \n\t" /* XMM0: f   f   f   f */
       "subps   %[V0011],%%xmm0       \n\t" /* XMM0: f-1 f-1 f   f */
       "xorps   %%xmm2,%%xmm2         \n\t" /* XMM2 will collect the low-precision values */
       "movaps  %[V2222],%%xmm5       \n\t"
       
       ADD4SSEA(0,8)

       "FLDL    %[tempFreq1]          \n\t" /* FLD D [TempFreqMinus15]   ;A15 */
       "FLDS    -16(%[Xalpha_k])      \n\t" /* FLD D [EBP-10h]   ;X14 A15 */
       "FMUL    %%ST(1),%%ST          \n\t" /* FMUL ST,ST1       ;X14A15 A15 */
       
       ADD4SSEB(16,24)
       
       "FLDS    -12(%[Xalpha_k])      \n\t" /* FLD D [EBP-0Ch]   ;Y14 X14A15 A15 */
       "FMUL    %%ST(2),%%ST          \n\t" /* FMUL ST,ST2       ;Y14A15 X14A15 A15 */
       "FXCH    %%ST(2)               \n\t" /* FXCH ST2          ;A15 X14A15 Y14A15 */

       ADD4SSEA(32,40)

       "FLDS    %[one]                \n\t" /* FLD D [_ONE_]     ;1 A15 X14A15 Y14A15 */
       "FADD    %%ST(1),%%ST          \n\t" /* FADD ST,ST1       ;A14 A15 X14A15 Y14A15 */
       
       ADD4SSEB(48,56)
       
       "FLDS    -8(%[Xalpha_k])       \n\t" /* FLD D [EBP-08h]   ;X15 A14 A15 X14A15 Y14A15 */
       "FMUL    %%ST(1),%%ST          \n\t" /* FMUL ST,ST1       ;X15A14 A14 A15 X14A15 Y14A15 */

       ADD4SSEA(64,72)

       "FADDP   %%ST,%%ST(3)          \n\t" /* FADDP ST3,ST      ;A14 A15 X' Y14A15 */
       "FMUL    %%ST,%%ST(1)          \n\t" /* FMUL ST1,ST       ;A14 Q145 X' Y14A15 */

       ADD4SSEB(80,88)

       "FLDS    -4(%[Xalpha_k])       \n\t" /* FLD D [EBP-04h]   ;Y15 A14 Q145 X' Y14A15 */
       "FMUL    %%ST(1),%%ST          \n\t" /* FMUL ST,ST1       ;Y15A14 A14 Q145 X' Y14A15 */

       ADD4SSEA(96,104) 

       "FADDP   %%ST,%%ST(4)          \n\t" /* FADDP ST4,ST      ;A14 Q145 X' Y' */
       "FSUBS   %[two]                \n\t" /* FSUB D [_TWO_]    ;A16 Q145 X' Y' */
       
       /* SSE: skip FPU calculated values */
       "addl    $144,%[Xalpha_kX]     \n\t" /* Xalpha_kX = Xalpha_kX + 4; */
       "subps   %[V4444],%%xmm0       \n\t"
       /* "subps   %%xmm5,%%xmm0         \n\t" */

       "FMUL    %%ST,%%ST(2)          \n\t" /* FMUL ST2,ST       ;A16 Q145 X'A16 Y' */
       "FMUL    %%ST,%%ST(3)          \n\t" /* FMUL ST3,ST       ;A16 Q145 X'A16 Y'A16 */
       
       ADD4SSEB(0,8)
       
       "FLDS    0(%[Xalpha_k])        \n\t" /* FLD D [EBP+00h]   ;X16 A16 Q145 X'A16 Y'A16 */
       "FMUL    %%ST(2),%%ST          \n\t" /* FMUL ST,ST2       ;X16Q145 A16 Q145 X'A16 Y'A16 */

       ADD4SSEA(16,24)
       
       "FADDP   %%ST,%%ST(3)          \n\t" /* FADDP ST3,ST      ;A16 Q145 X" Y'A16 */

       ADD4SSEB(32,40)

       "FLDS    4(%[Xalpha_k])        \n\t" /* FLD D [EBP+04h]   ;Y16 A16 Q145 X" Y'A16 */
       "FMUL    %%ST(2),%%ST          \n\t" /* FMUL ST,ST2       ;Y16Q145 A16 Q145 X" Y'A16 */
       
       ADD4SSEA(48,56)
       
       "FADDP   %%ST,%%ST(4)          \n\t" /* FADDP ST4,ST      ;A16 Q145 X" Y" */
       "FMUL    %%ST,%%ST(1)          \n\t" /* FMUL ST1,ST       ;A16 Q146 X" Y" */

       ADD4SSEB(64,72)

       "FSUBS   %[one]                \n\t" /* FSUB D [_ONE_]    ;A17 Q146 X" Y" */

       ADD4SSEA(80,88)

       "FMUL    %%ST,%%ST(2)          \n\t" /* FMUL ST2,ST       ;A17 Q146 X"A17 Y" */
       "FMUL    %%ST,%%ST(3)          \n\t" /* FMUL ST3,ST       ;A17 Q146 X"A17 Y"A17 */

       ADD4SSEB(96,104) 
       
       "FLDS    8(%[Xalpha_k])        \n\t" /* FLD D [EBP+08h]   ;X17 A17 Q146 X"A17 Y"A17 */
       "FMUL    %%ST(2),%%ST          \n\t" /* FMUL ST,ST2       ;X17Q146 A17 Q146 X"A17 Y"A17 */
       "FADDP   %%ST,%%ST(3)          \n\t" /* FADDP ST3,ST      ;A17 Q146 X! Y"A17 */
       
       /* add the two calculated parts of the real and imaginary part */
       "movhlps %%xmm2,%%xmm3         \n\t" /* XMM3: ? ? C_ImH C_ReH */
       "addps   %%xmm3,%%xmm2         \n\t" /* XMM2: - - C_Im C_Re */
       
       "FLDS    12(%[Xalpha_k])       \n\t" /* FLD D [EBP+0Ch]   ;Y17 A17 Q146 X! Y"A17 */
       "FMUL    %%ST(2),%%ST          \n\t" /* FMUL ST,ST2       ;Y17Q146 A17 Q146 X! Y"A17 */
       "FADDP   %%ST,%%ST(4)          \n\t" /* FADDP ST4,ST      ;A17 Q146 X! Y! */
       "FMULP   %%ST,%%ST(1)          \n\t" /* FMULP ST1,ST      ;Q147 X! Y! */
       
       /* store the result */
       "movlps  %%xmm2, %[XSumsX]     \n\t"
         
       "FDIVRS  %[one]                \n\t" /* FDIVR D [_ONE_]   ;1/Q x y */
       "FMUL    %%ST,%%ST(1)          \n\t" /* FMUL ST1,ST       ;1/Q xq y */
       "FMULP   %%ST,%%ST(2)          \n\t" /* FMULP ST2,ST      ;xq yq */
       
       "FSTPS   %[XRes]               \n\t"
       "FSTPS   %[XIms]               \n\t"
       
       /* interface */
       :
       /* output  (here: to memory)*/
       [cosq]        "=m" (Q.re),
       [cosv]        "=m" (tsincos.cs),
       [yRem]        "=m" (yRem),
       [XSumsX]      "=m" (XSumsX),
       [XRes]        "=m" (XRes),
       [XIms]        "=m" (XIms),

       /* input */
       [Xalpha_kX]   "+r" (Xalpha_kX) /* is changed by the code, so put into output section */
       :
       [Xalpha_k]    "r"  (Xalpha_kF),
       [tempFreq1]   "m"  (tempFreqF),
       [tempFreqX]   "m"  (tempFreqX),
       [x]           "m"  (tempFreq0),
       [yTemp]       "m"  (yTemp),

       /* constants */
       [lutr]        "m"  (lutr),
       [half]        "m"  (half),
       [one]         "m"  (one),
       [big]         "m"  (big),
       [two]         "m"  (two),
       [V0011]       "m"  (V0011[0]),
       [V2222]       "m"  (V2222[0]),
       [V4444]       "m"  (V4444[0]),

       /* tables */
       [csVal]       "m"  (csVal[0]),
       [csVal2PI]    "m"  (csVal2PI[0]),
       [csVal2PIPI]  "m"  (csVal2PIPI[0]),
       [diVal]       "m"  (diVal[0])

       :
       /* clobbered registers */
#ifndef IGNORE_XMM_REGISTERS
       "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7",
#endif
       "eax","ebx",
       "st","st(1)","st(2)","st(3)","st(4)"
       );           
      
          /* And last, we add the single and double precision values */
          XRes = XRes + XSumsX.re;
          XIms = XIms + XSumsX.im;
    
          {
            REAL4 tsin2pi =  tsincos.sn    * OOTWOPI;
            REAL4 tcos2pi = (tsincos.cs-1) * OOTWOPI;
            
            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          }
        }
        
        else /* if ( tempFreq0 >= LD_SMALL ) */
    
          {
/*	    fprintf(stderr,"small x\n"); */
	    
	    /* C version of the sin/cos calculations */

            /* calculation of tsin and tcos */
            {
              UINT4 idx  = tempFreq0 * LUT_RES +.5;
              REAL8 d    = tempFreq0 - diVal[idx];
              REAL8 d2   = d*d;
        
	      /* this originally was
		 tsin = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
		 tcos = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
		 However the sign is now incorporated in the tables, so we get:
	      */
              tsin = csVal[idx].sn + d * csVal2PI[idx].sn + d2 * csVal2PIPI[idx].sn;
              tcos = csVal[idx].cs + d * csVal2PI[idx].cs + d2 * csVal2PIPI[idx].cs;
              tcos -= 1.0;
            }
      
            /* calculation of yRem */
            {
#ifdef USE_FLOOR
              REAL8 yRem;
              if (yTemp >= 0) {
                yRem = yTemp - floor(yTemp);
              } else {
                /* yRem = yTemp - ceil(yTemp) + 1.0; */
                yRem = yTemp + floor(- yTemp) + 1.0;
              }
#else
              REAL8 yRem = yTemp - (INT4)(yTemp);
              if (yRem < 0) { yRem += 1.0f; } /* make sure this is in [0..1) */
#endif
        
              /* calculation of Q.re and Q.im */
              {
                UINT4 idx  = yRem * LUT_RES + .5;
                REAL8 d    = yRem - diVal[idx];
                REAL8 d2   = d*d;
                Q.im = csVal[idx].sn + d * csVal2PI[idx].sn + d2 * csVal2PIPI[idx].sn;
                Q.re = csVal[idx].cs + d * csVal2PI[idx].cs + d2 * csVal2PIPI[idx].cs;
              }
            }


            /* special version of the Kernel loop*/

            realXP=0.0;
            imagXP=0.0;
      
            /* Loop over terms in Dirichlet Kernel */
            for(k=0; k < klim ; k++)
              {
                COMPLEX8 Xalpha_k = Xalpha[sftIndex];
                sftIndex ++;
                /* If x is small we need correct x->0 limit of Dirichlet kernel */
                if( fabs(x) <  SMALL) 
                  {
                    realXP += Xalpha_k.re;
                    imagXP += Xalpha_k.im;
                  }      
                else
                  {
                    realP = tsin / x;
                    imagP = tcos / x;
                    /* these four lines compute P*xtilde */
                    realXP += Xalpha_k.re * realP;
                    realXP -= Xalpha_k.im * imagP;
                    imagXP += Xalpha_k.re * imagP;
                    imagXP += Xalpha_k.im * realP;
                  }
                
                tempFreq1 --;
                x = LAL_TWOPI * tempFreq1;
                
              } /* for k < klim */
            
          } /* if x could become close to 0 */
        
        /* implementation of amplitude demodulation */
	/* pulled in sign change of Q.im */
        {
          REAL8 realQXP = realXP*Q.re+imagXP*Q.im;
          REAL8 imagQXP = imagXP*Q.re-realXP*Q.im;
          Fa.re += a*realQXP;
          Fa.im += a*imagQXP;
          Fb.re += b*realQXP;
          Fb.im += b*imagQXP;
        }
      }      

    FaSq = Fa.re*Fa.re+Fa.im*Fa.im;
    FbSq = Fb.re*Fb.re+Fb.im*Fb.im;
    FaFb = Fa.re*Fb.re+Fa.im*Fb.im;
                        
    Fs->F[i] = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb);
    if (params->returnFaFb)
      {
        Fs->Fa[i] = Fa;
        Fs->Fb[i] = Fb;
      }


  }
  /* Clean up */
  LALFree(tempInt1);
  LALFree(xSum);
  LALFree(ySum);
  
  RETURN( status );
}
