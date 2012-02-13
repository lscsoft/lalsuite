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

/* LALDemod variant with sin/cos calculation coded in (GNU) assembler
 * Authors see ComputeFStatistic.c
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
#else /* x86 GNU as */
#define AD_FLOAT   ".float "
#define AD_ASCII   ".string "
#define AD_ALIGN16 ".align 16"
#define AD_ALIGN64 ".align 64"
#endif

#ifndef CPU_TYPE_S
#define CPU_TYPE_S "0"
#endif

void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 
{ 

  INT4  alpha,i,k=0;            /* loop indices */
  REAL8 *xSum=NULL, *ySum=NULL; /* temp variables for computation of fs*as and fs*bs */
  INT4  s;                      /* local variable for spinDwn calcs. */
  REAL8 xTemp;                  /* temp variable for phase model */
  REAL4 xTInt;                  /* integer part of xTemp */
  REAL8 deltaF;                 /* width of SFT band */
  REAL8 *skyConst;              /* vector of sky constants data */
  REAL8 *spinDwn;               /* vector of spinDwn parameters (maybe a structure? */
  INT4  spOrder;                /* maximum spinDwn order */
  INT4  nDeltaF;                /* number of frequency bins per SFT band */

  INT4  sftIndex;               /* more temp variables */
  REAL8 realXP, imagXP;
  REAL4 realP, imagP;
  static REAL8 realQ, imagQ;
  INT4  *tempInt1;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  REAL8 f;
  COMPLEX16 Fa, Fb;


#ifdef USE_BOINC                /* this actually is only a hint. The assembler Kernel loop */
#define klim 32                 /* doesn't use this value at all, but relies on it being 32 */ 
#else
  UINT4 klim = 2*params->Dterms;
#endif

  static REAL4 tsin, tcos;

  /* Lookup tables for fast sin/cos calculation */
  static REAL8 sinVal[LUT_RES+1];
  static REAL8 sinVal2PI[LUT_RES+1];
  static REAL8 sinVal2PIPI[LUT_RES+1];
  static REAL8 cosVal[LUT_RES+1];
  static REAL8 cosVal2PI[LUT_RES+1];
  static REAL8 cosVal2PIPI[LUT_RES+1];
  static REAL8 diVal[LUT_RES+1];
  static BOOLEAN sincos_initialized = 1; /* reset after initializing the sin/cos tables */

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
  /* using three tables with values including PI is simply faster */
  if ( sincos_initialized )
    {
      for (k=0; k <= LUT_RES; k++) {
        sinVal[k]      = sin((LAL_TWOPI*k)/(LUT_RES));
        sinVal2PI[k]   = sinVal[k]    * LAL_TWOPI;
        sinVal2PIPI[k] = sinVal2PI[k] * LAL_PI;
        cosVal[k]      = cos((LAL_TWOPI*k)/(LUT_RES));
        cosVal2PI[k]   = cosVal[k]    * LAL_TWOPI;
        cosVal2PIPI[k] = cosVal2PI[k] * LAL_PI;
      }

      /* this additional table saves another "costly" division in sin/cos calculation */
      for (k=0; k <= LUT_RES; k++)
        diVal[k] = (REAL8)k/(REAL8)(LUT_RES);

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
#ifndef USE_BOINC
            fprintf (stderr, "\n*** PLEASE report this bug to pulgroup@gravity.phys.uwm.edu *** \n\n");
#endif
            exit (COMPUTEFSTAT_EXIT_DEMOD);
        }
        if (xTemp < 0) {
            fprintf (stderr, "xTemp >= 0 failed\n");
            fprintf (stderr, "DEBUG: loop=%d, xTemp=%f, f=%f, alpha=%d, tempInt1[alpha]=%d\n", 
                     i, xTemp, f, alpha, tempInt1[alpha]);
            fprintf (stderr, "DEBUG: skyConst[ tempInt1[ alpha ] ] = %f, xSum[ alpha ]=%f\n",
                     skyConst[ tempInt1[ alpha ] ], xSum[ alpha ]);
#ifndef USE_BOINC
            fprintf (stderr, "\n*** PLEASE report this bug to pulgroup@gravity.phys.uwm.edu *** \n\n");
#endif
            exit (COMPUTEFSTAT_EXIT_DEMOD);
        }

#ifdef USE_FLOOR
        xTInt =  floor(xTemp);
#else
        xTInt =  (UINT4)xTemp;
#endif
        tempFreq0 = xTemp - xTInt;   /* lies in [0, +1) by definition */

        sftIndex = xTInt - params->Dterms + 1 - params->ifmin;

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

          /* Variables for Kernel loop */
	  COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
          REAL4    tempFreqX  = tempFreq1;        /* REAL4 because of SSE */
	  REAL4    XRes=0.0, XIms=0.0;            /* sums of Xa.re and Xa.im */
	  REAL4    accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */

          /* constants in memory */
          static REAL4 lutr = LUT_RES; /* LUT_RES in memory */
          static REAL4 half = .5;
          static REAL4 one  = 1.0;
          static REAL4 two  = 2.0;
#ifdef USE_ASS_2PI_MULT
	  static REAL4 pi2  = OOTWOPI; 
#endif
          /* vector constants */
          /* having these not aligned will crash the assembler code */
          static REAL4 V0011[4] __attribute__ ((aligned (16))) = { 0,0,1,1 };
          static REAL4 V2222[4] __attribute__ ((aligned (16))) = { 2,2,2,2 };
          static REAL4 V4444[4] __attribute__ ((aligned (16))) = { 4,4,4,4 };
 
          static REAL8 yRem;

          __asm __volatile
      (
       /* calculation of tsin and tcos */
       /* ---------------------------- */

       /* calculate index and put into EAX */
       /*                                    vvvvv-- these comments keep track of the FPU stack */
       "fldl   %[x]                   \n\t" /* x */
       "flds   %[lutr]                \n\t" /* LUT_RES x */
       "fmul   %%st(1),%%st(0)        \n\t" /* (x*LUT_RES) x */
       
#ifdef USE_DEFAULT_ROUNDING_MODE
       /* The default rounding mode is round-to-nearest, which is exactly what we want here.
          It should be restored by the compiler after every operation that involves floating-point-to-integer
          conversion (rounding). Switching the rounding mode is slow, so the following code is the fastest
          possible, as it simply expects the default rounding mode being active. However relying on this is
          kind of dangerous, so we don't do this by default. */
#ifndef USE_DEPENDENCIES
       "fistpl %[sinv]                \n\t" /* x */
       "movl   %[sinv],%%eax          \n\t" /* x */
#else
       /* build dependency chain to avoid memory re-read
	  Based on idea from Agner Fog
	  He thinks that should be faster on P4 */
       "fistl  %[sinv]                \n\t" /* store without popping */
       "fcomip %%st,%%st              \n\t" /* compare and pop, make flags depend on ST */
       "setc   %%al                   \n\t" /* make AL depend on flags */
       "andl   $0,%%eax               \n\t" /* set to 0 */
       "movl   %[sinv](%%eax),%%eax   \n\t" /* make dependent on EAX */
#endif
#else
       "fadds  %[half]                \n\t" /* (x*LUT_RES+.5) x */
       /* Implementation of floor() that doesn't rely on a rounding mode being set
          The current way works for positive values only! */
       /* This code temporary stores integer values in memory locations of float variables,
          but they're overwritten later anyway */
       "fistl  %[sinv]                \n\t" /* (x*LUT_RES+.5) x */ /* saving the rounded value, the original in FPU */
       "fildl  %[sinv]                \n\t" /* ((int)(x*LUT_RES+.5)) (x*LUT_RES+.5) x */
       "fsubrp                        \n\t" /* ((x*LUT_RES+.5)-((int)(x*LUT_RES+.5))) x */
       "movl   %[sinv],%%eax          \n\t" /* ((x*LUT_RES+.5)-((int)(x*LUT_RES+.5))) x */ /* get the result in EAX */
       "fstps  %[cosv]                \n\t" /* x */                /* value - round(value) will be negative if was rounding up */
       "shll   $1,%[cosv]             \n\t" /* x */                /* shift the sign bit into the carry bit */
       "sbb    $0,%%eax               \n\t" /* x */                /* substract the carry in EAX */
#endif
       
       /* calculate d = x - diVal[idx] in st(0) */
       "fsubl  %[diVal](,%%eax,8)     \n\t" /* (d = x - diVal[i]) */
       
       /* this calculates d*d on _top_ of the stack */ 
       /* copy d on the stack to prepare for calculating d*d */
       "fld    %%st(0)                \n\t" /* d d */
       "fmul   %%st(0),%%st(0)        \n\t" /* (d*d) d */
       
       /* three-term Taylor expansion for sin value, leaving d and d*d on stack, idx kept in eax */
       /* some reordering by Akos */
       "fldl  %[cosVal2PI](,%%eax,8)  \n\t" /* cosVal2PI[i] (d*d) d */
       "fmul  %%st(2),%%st(0)         \n\t" /* (d*cosVal2PI[i]) (d*d) d */
       "faddl %[sinVal](,%%eax,8)     \n\t" /* (sinVal[i]+d*cosVal2PI[i]) (d*d) d */
       "fldl  %[sinVal2PIPI](,%%eax,8)\n\t" /* sinVal2PIPI[i] (sinVal[i]+d*cosVal2PI[i]) (d*d) d */
       "fmul  %%st(2),%%st(0)         \n\t" /* (d*d*sinVal2PIPI[i]) (sinVal[i]+d*cosVal2PI[i]) (d*d) d */
       "fsubrp                        \n\t" /* ((sinVal[i]+d*cosVal2PI[i])-d*d*sinVal2PIPI[i]) (d*d) d */

       /* special: * TWOPI */
#ifdef USE_ASS_2PI_MULT
       "fmuls %[pi2]                  \n\t" /* (((sinVal[i]+d*cosVal2PI[i])-d*d*sinVal2PIPI[i])*2*pi) (d*d) d */
#endif
       "fstps %[sinv]                 \n\t" /* (d*d) d */
       
       /* calculation for cos value, this time popping the stack */
       /* this ordering mimics the calculation of cos as done by gcc (4.1.0)  */
       
       "fmull %[cosVal2PIPI](,%%eax,8)\n\t" /* (d*d*cosVal2PIPI[i]) d */
       "fsubrl %[cosVal](,%%eax,8)    \n\t" /* (cosVal[i]-d*d*cosVal2PIPI[i]) d */
       "fxch                          \n\t" /* d (cosVal[i]-d*d*cosVal2PIPI[i]) d */
       "fmull %[sinVal2PI](,%%eax,8)  \n\t" /* (d*sinVal2PI[i]) (cosVal[i]-d*d*cosVal2PIPI[i]) d */
       "fsubrp                        \n\t" /* ((cosVal[i]-d*d*cosVal2PIPI[i])-d*sinVal2PI[i]) */
       
#ifndef SKIP_COS_ROUNDING
       /* this stores the cos value in the single precision output variable before substracting 1.0,
          again to mimic the compilers code to get the same result. It is faster and might be better
          to simply leave this value on the stack, so we add a switch for skipping this */
       "fstps %[cosv]                 \n\t" /* % */
       "flds  %[cosv]                 \n\t" /* % */
#endif
       
       /* special here: tcos = (tcos - 1.0) * TWOPI */
       "fsubs %[one]                  \n\t" /* (cosVal[i]-d*(sinVal2PI[i]-d*cosVal2PIPI[i])-1) */
#ifdef USE_ASS_2PI_MULT
       "fmuls %[pi2]                  \n\t" /* (((sinVal[i]+d*cosVal2PI[i])-d*d*sinVal2PIPI[i])*2*pi) (d*d) d */
#endif
       "fstps %[cosv]                 \n\t" /* % */
       

       /* calculating yRem */
       /* ---------------- */

       "fldl   %[yTemp]               \n\t" /* yT */
       /* it requires some mind-bending to come from the algorithm used in the USE_FLOOR case
          of the yRem calculation in the C code below (in case yTemp < 0) to the case distinction
          implemented here.
          Hints:  yTemp<0 => yRem = yTemp - ceil(yTemp) + 1.0;
                  y<0 => ceil(y) = ((y-(int)y)<-.5) ? ((int)y-1.0) : ((int)y)
                  -1.0 + 1.0 = 0; 0 + 1.0 = +1.0 */
       /* The two cases (yTemp < 0) and (yTemp >= 0) then come out identical,
          so a case distinction is unnecessary. Furthermore this case turns out to be the
          same code independent of relying on the default rounding mode or not */
       "fistl  %[yRem]                \n\t" /* yT */
       "fisubl %[yRem]                \n\t" /* yT-(int)yT */ /* value - round(value) will be positive if was rounding up */
       "fsts   %[yRem]                \n\t" /* yT-(int)yT */ /* we will check the sign in integer registers */
       "sub    %%eax,%%eax            \n\t" /* yT-(int)yT */ /* EAX=0 */
       "orl    %[yRem],%%eax          \n\t" /* yT-(int)yT */ /* sets the S (sign) flag */
       "jns    yrem" CPU_TYPE_S      "\n\t" /* yT-(int)yT */ /* append the CPU_TYPE to label to allow more instances */
       "fadds  %[one]                 \n"   /* (yT-(int)yT+1) */
       "yrem" CPU_TYPE_S ":           \n\t" /* yR */
       
       /* calculation of realQ and imagQ */
       /* ------------------------------ */

       /* yRem still in st(0), following named x */
       /* see also comments in calculation of tsin and tcos */
       
       "flds   %[lutr]                \n\t" /* LUT_RES x */
       "fmul   %%st(1),%%st(0)        \n\t" /* (x*LUT_RES) x */
       
#ifdef USE_DEFAULT_ROUNDING_MODE
#ifndef USE_DEPENDENCIES
       "fistpl %[sinq]                \n\t" /* x */
       "movl   %[sinq],%%eax          \n\t" /* x */
#else
       "fistl  %[sinq]                \n\t" /* store without popping */
       "fcomip %%st,%%st              \n\t" /* compare and pop, make flags depend on ST */
       "setc   %%al                   \n\t" /* make AL depend on flags */
       "andl   $0,%%eax               \n\t" /* set to 0 */
       "movl   %[sinq](%%eax),%%eax   \n\t" /* make dependent on EAX */
#endif
#else
       "fadds  %[half]                \n\t" /* (x*LUT_RES+.5) x */
       "fistl  %[sinq]                \n\t" /* (x*LUT_RES+.5) x */ /* saving the rounded value, the original in FPU */
       "fildl  %[sinq]                \n\t" /* ((int)(x*LUT_RES+.5)) (x*LUT_RES+.5) x */
       "fsubrp                        \n\t" /* ((x*LUT_RES+.5)-((int)(x*LUT_RES+.5))) x */
       "movl   %[sinq],%%eax          \n\t" /* ((x*LUT_RES+.5)-((int)(x*LUT_RES+.5))) x */ /* get the result in EAX */
       "fstps  %[cosq]                \n\t" /* x */                /* value - round(value) will be negative if was rounding up */
       "shll   $1,%[cosq]             \n\t" /* x */                /* shift the sign bit into the carry bit */
       "sbb    $0,%%eax               \n\t" /* x */                /* substract the carry in EAX */
#endif
       "fsubl  %[diVal](,%%eax,8)     \n\t" /* (d = x - diVal[i]) */
       "fld    %%st(0)                \n\t" /* d d */
       "fmul   %%st(0),%%st(0)        \n\t" /* (d*d) d */

       "fldl  %[cosVal2PI](,%%eax,8)  \n\t" /* cosVal2PI[i] (d*d) d */
       "fmul  %%st(2),%%st(0)         \n\t" /* (d*cosVal2PI[i]) (d*d) d */
       "faddl %[sinVal](,%%eax,8)     \n\t" /* (sinVal[i]+d*cosVal2PI[i]) (d*d) d */
       "fldl  %[sinVal2PIPI](,%%eax,8)\n\t" /* sinVal2PIPI[i] (sinVal[i]+d*cosVal2PI[i]) (d*d) d */
       "fmul  %%st(2),%%st(0)         \n\t" /* (d*d*sinVal2PIPI[i]) (sinVal[i]+d*cosVal2PI[i]) (d*d) d */
       "fsubrp                        \n\t" /* ((sinVal[i]+d*cosVal2PI[i])-d*d*sinVal2PIPI[i]) (d*d) d */

       /* special here: sign change of the sin() result */
       "fchs                          \n\t" /* -(sinVal[i]+d*cosVal2PI[i]-d*d*sinVal2PIPI[i]) (d*d) d */
       "fstpl %[sinq]                 \n\t" /* (d*d) d */
       
       "fmull %[cosVal2PIPI](,%%eax,8)\n\t" /* (d*d*cosVal2PIPI[i]) d */
       "fsubrl %[cosVal](,%%eax,8)    \n\t" /* (cosVal[i]-d*d*cosVal2PIPI[i]) d */
       "fxch                          \n\t" /* d (cosVal[i]-d*d*cosVal2PIPI[i]) d */
       "fmull %[sinVal2PI](,%%eax,8)  \n\t" /* (d*sinVal2PI[i]) (cosVal[i]-d*d*cosVal2PIPI[i]) d */
       "fsubrp                        \n\t" /* ((cosVal[i]-d*d*cosVal2PIPI[i])-d*sinVal2PI[i]) */
       "fstpl %[cosq]                 \n\t" /* % */
       
       /* interface */
       :
       /* output  (here: to memory)*/
       [sinq]        "=m" (imagQ),
       [cosq]        "=m" (realQ),
       [sinv]        "=m" (tsin),
       [cosv]        "=m" (tcos),
       [yRem]        "=m" (yRem)
       :
       /* input */
       [x]           "m"  (tempFreq0),
       [yTemp]       "m"  (yTemp),

       /* constants */
       [lutr]        "m"  (lutr),
       [half]        "m"  (half),
       [one]         "m"  (one),
       [two]         "m"  (two),
#ifdef USE_ASS_2PI_MULT
       [pi2]         "m"  (pi2),
#endif

       /* tables */
       [sinVal]      "m"  (sinVal[0]),
       [cosVal]      "m"  (cosVal[0]),
       [sinVal2PI]   "m"  (sinVal2PI[0]),
       [cosVal2PI]   "m"  (cosVal2PI[0]),
       [sinVal2PIPI] "m"  (sinVal2PIPI[0]),
       [cosVal2PIPI] "m"  (cosVal2PIPI[0]),
       [diVal]       "m"  (diVal[0])

       :
       /* clobbered registers */
       "eax",
       "st","st(1)","st(2)","st(3)","st(4)"
       );           
      
	  for(k=0; k < klim ; k++)
	    {
	      XRes = tempFreq1 * XRes + (*Xalpha_k).re * accFreq;
	      XIms = tempFreq1 * XIms + (*Xalpha_k).im * accFreq;
	      
	      accFreq *= (REAL4)tempFreq1;
	      tempFreq1 --;
	      Xalpha_k ++;
	    } /* for k < klim */
	  
	  XRes /= accFreq;
	  XIms /= accFreq;
	  
#ifndef USE_ASS_2PI_MULT
          {
            REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
            REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
            
            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
	  }
#else
	  realXP = tsin * XRes - tcos * XIms;
	  imagXP = tcos * XRes + tsin * XIms;
#endif
        }
        
        else /* if ( tempFreq0 >= LD_SMALL ) */
    
          {
/*            fprintf(stderr,"small x\n"); */
            
            /* C version of the same calculations */

            /* calculation of tsin and tcos */
            {
              UINT4 idx  = tempFreq0 * LUT_RES +.5;
              REAL8 d    = tempFreq0 - diVal[idx];
              REAL8 d2   = d*d;
        
              tsin = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
              tcos = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
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
        
              /* calculation of realQ and imagQ */
              {
                UINT4 idx  = yRem * LUT_RES + .5;
                REAL8 d    = yRem - diVal[idx];
                REAL8 d2   = d*d;
                imagQ = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
                realQ = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
          
                imagQ = -imagQ;
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
        {
          REAL8 realQXP = realXP*realQ-imagXP*imagQ;
          REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
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
