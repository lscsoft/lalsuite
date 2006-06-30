/* LALDemod variant with structs and signatures for modifications by Feket Akos
 * Authors see ComputeFStatistic.c
                                                         Bernd Machenschalk */
RCSID( "$Id$");

static REAL8 sinTab[LUT_RES+1]; /* Lookup tables for fast sin/cos calculation */
static REAL8 sinTab2PI[LUT_RES+1];
static REAL8 sinTab2PIPI[LUT_RES+1];
static REAL8 cosTab[LUT_RES+1];
static REAL8 cosTab2PI[LUT_RES+1];
static REAL8 cosTab2PIPI[LUT_RES+1];
static REAL8 diviTab[LUT_RES+1];
static BOOLEAN firstCall = 1;


#define klim 32

/* __attribute__ ((always_inline)) */
void LALDemodSub(COMPLEX8* Xalpha, INT4 sftIndex,
		 REAL8 tempFreq0, REAL8 tempFreq1, REAL8 x, REAL8 yTemp,
		 REAL8* realXPo, REAL8* imagXPo, REAL8* realQo, REAL8* imagQo)
{
  REAL4 tsin, tcos;
  REAL4 realP, imagP;
  INT4 k;                             /* loop counter */
  REAL8 realXP, imagXP, realQ, imagQ; /* output parameters */
#ifdef DEBUG
  REAL4 tsin2, tcos2;
  REAL8 t1,t2,t3,t4;
#endif

  /* calculate tsin, tcos, realQ and imagQ from sin/cos LUT */

#ifdef USE_SINCOS_GAS

  static REAL4 lutr = LUT_RES; /* LUT_RES in memory */
  static REAL4 half = .5;

  // fprintf(stderr,"LALDemodSub called\n");

  __asm
    (
     /* calculate index and put into EAX */
     /*                                    vvvvv-- these comments keep track of the FPU stack */
     "fldl   %[x]                   \n\t" /* x */
     "fmul   %[lutr]                \n\t" /* LUT_RES x */

#ifdef USE_DEFAULT_ROUNDING_MODE
     "fistpl %[sinv]                \n\t" /* x */
     "movl   %[sinv],%%ebx          \n\t" /* x */
#else
     "fadds  %[half]                \n\t" /* (x*LUT_RES+.5) x */
                                          /* NOTE: this temporary stores an integer in a memory location of a
					           float variable, but it's overwritten later anyway */
     "fistl  %[sinv]                \n\t" /* (x*LUT_RES+.5) x */ /* saving the rounded value, the original in FPU */
     "fisubl %[sinv]                \n\t" /* (x*LUT_RES+.5) x */ /* value - round(value) will be negative if was rounding up */
     "fstps  %[cosv]                \n\t" /* x */                /* we will check the sign in integer registers */
     "movl   %[sinv],%%ebx          \n\t" /* x */
     "sub    %%eax,%%eax            \n\t" /* x */                /* EAX=0 */
     "orl    %[cosv],%%eax          \n\t" /* x */                /* it will set the S (sign) flag */
     "jns    sincos1                \n\t" /* x */                /* the result is ok, rounding = truncation */
     "dec    %%ebx                  \n\t" /* x */                /* sinv=sinv-1.0 (it was rounded up) */
     "sincos1:                      \n\t" /* x */
#endif

     "sal    $3,%%ebx               \n\t" /* x */ /* <<3 = *sizeof(REAL8) */

     /* calculate d = x - diviTab[idx] in st(0) */
     "fsubr  %[diviTab](%%ebx)      \n\t" /* (d = x - diviTab[i]) */

     /* add d*d on the stack */
     "fld    %%st(0)                 \n\t" /* d d */
     "fmul   %%st(0),%%st(0)         \n\t" /* (d*d) d */

     /* three-term Taylor expansion for sin value, starting with the last term,
	leaving d and d*d on stack, idx kept in EAX */
     "fldl   %[sinTab2PIPI](%%ebx)  \n\t" /* sinTab2PIPI[i] (d*d) d */
     "fmul   %%st(1),%%st(0)        \n\t" /* (d*d*sinTab2PIPI[i]) (d*d) d */

     "fldl   %[cosTab2PI](%%ebx)    \n\t" /* cosTab2PI[i] (d*d*sinTab2PIPI[i]) (d*d) d */
     "fmul   %%st(3),%%st(0)        \n\t" /* (d*cosTab2PI[i]) (d*d*sinTab2PIPI[i]) (d*d) d */
     "fsubp                         \n\t" /* (d*cosTab2PI[i]-d*d*sinTab2PIPI[i]) (d*d) d */

     "fadd   %[sinTab](%%ebx)       \n\t" /* (sinTab[i]+d*cosTab2PI[i]-d*d*sinTab2PIPI[i]) (d*d) d */
     "fstps  %[sinv]                \n\t" /* (d*d) d */

     /* similar calculation for cos value, this time popping the stack */
     "fmul   %[cosTab2PIPI](%%ebx)  \n\t" /* (d*d*cosTab2PIPI[i]) d */

     "fxch                          \n\t" /* d (d*d*cosTab2PIPI[i]) */
     "fmul   %[sinTab2PI](%%ebx)    \n\t" /* (d*sinTab2PI[i]) (d*d*cosTab2PIPI[i]) */
     "fsubp                         \n\t" /* (d*sinTab2PI[i]-d*d*cosTab2PIPI[i]) */
     
     "fsubr  %[cosTab](%%ebx)       \n\t" /* (cosTab[i]-d*sinTab2PI[i]-d*d*cosTab2PIPI[i])) */
#ifdef DEBUG
     "fstl  %[t4]\n\t"
#endif

#ifndef SKIP_COS_ROUNDING
     "fstps %[cosv]                 \n\t" /* % */

     /* special here: tcos -= 1.0 */
     "flds  %[cosv]                 \n\t" /* % */
#endif
     "fld1                          \n\t" /* 1 (cosTab[i]-d*(sinTab2PI[i]-d*cosTab2PIPI[i])) */
     "fsubrp                        \n\t" /* (cosTab[i]-d*(sinTab2PI[i]-d*cosTab2PIPI[i])-1) */
     "fstps %[cosv]                 \n\t" /* % */

     : /* output */
#ifdef DEBUG
     [t1]    "=m" (t1),
     [t2]    "=m" (t2),
     [t3]    "=m" (t3),
     [t4]    "=m" (t4),
#endif
     [sinv]  "=m" (tsin),
     [cosv]  "=m" (tcos)

     : /* input */
     [x]           "m" (tempFreq0),
     [lutr]        "m" (lutr),
     [half]        "m" (half),
     [sinTab]      "m" (sinTab[0]),
     [cosTab]      "m" (cosTab[0]),
     [sinTab2PI]   "m" (sinTab2PI[0]),
     [cosTab2PI]   "m" (cosTab2PI[0]),
     [sinTab2PIPI] "m" (sinTab2PIPI[0]),
     [cosTab2PIPI] "m" (cosTab2PIPI[0]),
     [diviTab]     "m" (diviTab[0])

     : /* clobbered registers */
     "eax", "ebx", "st","st(1)","st(2)","st(3)"
     );	

#ifdef DEBUG
  {
    UINT4 idx  = tempFreq0 * LUT_RES +.5;
    REAL8 d    = tempFreq0 - diviTab[idx];
    REAL8 d2   = d*d;
    REAL8 t11,t12,t13,t14;
    
    t11 = d2;

    tsin2 = sinTab[idx] + d * cosTab2PI[idx] - d2 * sinTab2PIPI[idx];
    tcos2 = (t14= (t13= cosTab[idx]) - (t12= d * sinTab2PI[idx] - d2 * cosTab2PIPI[idx]));

    tcos2 -= 1.0;
    
    if (fabs(tsin - tsin2) > 1e-18)
      fprintf(stderr,"\n%f: %e\n",tempFreq0,tsin-tsin2);

    if (fabs(tcos - tcos2) > 1e-18){
      fprintf(stderr,"\n%f: %e\n",tempFreq0,tcos-tcos2);
      // fprintf(stderr,"%.20f\n%.20f\n",t1,t11);
      fprintf(stderr,"%.20f\n%.20f\n",t2,t12);
      fprintf(stderr,"%.20f\n%.20f\n",t3,t13);
      fprintf(stderr,"%.20f\n%.20f\n",t4,t14);
      fprintf(stderr,"%.20f\n%.20f\n",tcos,tcos2);
    }
  }
#endif

#else

  {
    UINT4 idx  = tempFreq0 * LUT_RES +.5;
    REAL8 d    = tempFreq0 - diviTab[idx];
    REAL8 d2   = d*d;

    tsin = sinTab[idx] + d * cosTab2PI[idx] - d2 * sinTab2PIPI[idx];
    tcos = cosTab[idx] - d * sinTab2PI[idx] - d2 * cosTab2PIPI[idx];

    tcos -= 1.0;
    
  }
#endif
  
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
    {
      UINT4 idx  = yRem * LUT_RES + .5;
      REAL8 d    = yRem - diviTab[idx];
      REAL8 d2   = d*d;
      imagQ = sinTab[idx] + d * cosTab2PI[idx] - d2 * sinTab2PIPI[idx];
      realQ = cosTab[idx] - d * sinTab2PI[idx] - d2 * cosTab2PIPI[idx];
      
      imagQ = -imagQ;
    }
  }
  
  /* we branch now (instead of inside the central loop)
   * depending on wether x can ever become SMALL in the loop or not, 
   * because it requires special treatment in the Dirichlet kernel
   */
  if ( tempFreq0 < LD_SMALL ) 
    {
      
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
  else
    /* when optimizing for a specific architecture, we usually don't want to rewrite
       the whole file, but only the following block of C code, so we insert another
       hook here, in the hope that we still maintain readability. Also too we'd like
       to avoid the necessarity to keep changes in other parts of the file in sync
       between these versions.                                                  BM */
#if defined(USE_X86_GAS)
#include "CFSLALDemodLoop_x86gAss.c"
#elif defined(USE_X86_MAS)
#include "CFSLALDemodLoop_x86MSAss.c"
#elif defined(USE_ALTIVEC)
#include "CFSLALDemodLoop_AltiVec.c"
#elif defined(USE_NDP_VECT)
#include "CFSLALDemodLoop_ndp_vect.c"
#elif defined(USE_NEW_DIV_PART)
#include "CFSLALDemodLoop_div_part.c"
#else
    {
      COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
      
      realXP=0.0;
      imagXP=0.0;
      
      /* Loop over terms in Dirichlet Kernel */
      
      
      for(k=0; k < klim ; k++)
	{
	  REAL4 xinv = (REAL4)OOTWOPI / (REAL4)tempFreq1;
	  COMPLEX8 Xa = *Xalpha_k;
	  Xalpha_k ++;
	  tempFreq1 --;
	  
	  realP = tsin * xinv;
	  imagP = tcos * xinv;
	  /* these lines compute P*xtilde */
	  realXP += Xa.re * realP - Xa.im * imagP;
	  imagXP += Xa.re * imagP + Xa.im * realP;
	  
	} /* for k < klim */
      
    } /* if x cannot be close to 0 */
#endif

  *realXPo = realXP;
  *imagXPo = imagXP;
  *realQo  = realQ;
  *imagQo  = imagQ;
}


void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 
{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8 *xSum=NULL, *ySum=NULL; /* temp variables for computation of fs*as and fs*bs */
  INT4 s;                       /* local variable for spinDwn calcs. */
  REAL8 xTemp;                  /* temp variable for phase model */
  REAL4 xTInt;                  /* integer part of xTemp */
  REAL8 deltaF;                 /* width of SFT band */
  UINT4 k=0;
  REAL8 *skyConst;              /* vector of sky constants data */
  REAL8 *spinDwn;               /* vector of spinDwn parameters (maybe a structure? */
  INT4  spOrder;                /* maximum spinDwn order */
  REAL8 realXP, imagXP;         /* temp variables used in computation of */
  INT4  nDeltaF;                /* number of frequency bins per SFT band */
  INT4  sftIndex;               /* more temp variables */
  REAL8 realQ, imagQ;
  INT4 *tempInt1;
  /* UINT4 index; */
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
  REAL8 f;

  REAL8 A=params->amcoe->A;
  REAL8 B=params->amcoe->B;
  REAL8 C=params->amcoe->C;
  REAL8 D=params->amcoe->D;

  UINT4 M=params->SFTno;

  INITSTATUS( status, "TestLALDemod", rcsid );

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
  if ( firstCall )
    {
      for (k=0; k <= LUT_RES; k++) {
	sinTab[k]      = sin((LAL_TWOPI*k)/(LUT_RES));
	sinTab2PI[k]   = sinTab[k]    * LAL_TWOPI;
	sinTab2PIPI[k] = sinTab2PI[k] * LAL_PI;
	cosTab[k]      = cos((LAL_TWOPI*k)/(LUT_RES));
	cosTab2PI[k]   = cosTab[k]    * LAL_TWOPI;
	cosTab2PIPI[k] = cosTab2PI[k] * LAL_PI;
      }

      for (k=0; k <= LUT_RES; k++)
	diviTab[k] = (REAL8)k/(REAL8)(LUT_RES);
      firstCall = 0;
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

    f=params->f0+i*params->df;

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
         * this can be seen from there definition (-> documentation)
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

	LALDemodSub(Xalpha, sftIndex,
		    tempFreq0, tempFreq1, x, yTemp,
		    &realXP, &imagXP, &realQ, &imagQ);

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
