/* LALDemod variant with structs and signatures for modifications by Feket Akos
 * Authors see ComputeFStatistic.c
                                                         Bernd Machenschalk */
RCSID( "$Id$");

static REAL8 *sinVal, *sinVal2PI, *sinVal2PIPI;
static REAL8 *cosVal, *cosVal2PI, *cosVal2PIPI;
static REAL8 *divTab;

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

  /* calculate tsin, tcos, realQ and imagQ from sin/cos LUT */

#ifdef USE_SINCOS_GAS

  REAL4 lutr = LUT_RES; /* LUT_RES in memory */

  __asm __volatile
    (

     /* calculate idx and put into EAX */
     "FLDL   %[x]                   \n\t" /* LUT_RES */
     "FLDS   %[lutr]                \n\t" /* LUT_RES x */
     "FMUL   %%st(1),%%st(0)        \n\t" /* (x*LUT_RES) x */
     "FISTPL %[sinv]                \n\t" /* x */
     "MOVL   %[sinv],%%EAX          \n\t" /* x */
     "SAL    $3,%%EAX               \n\t" /* x */ /* <<3 */

     /* calculate d = divLUTtab[idx] - x in st(0) */
     "MOVL  %[divLUTtab],%%EDX      \n\t" /* x */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* divLUTtab[i] x */
     "FSUBRP                        \n\t" /* (d = tempFreq0 - divLUTtab[i]) */

     /* three-term Taylor expansion for sin value, leaving d on stack, idx kept in EAX */
     "MOV   %[sinVal2PIPI],%%EDX    \n\t" /* d */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* sinVal2PIPI[i] d */
     "FMUL  %%st(1),%%st(0)         \n\t" /* (d*sinVal2PIPI[i]) d */

     "MOV   %[cosVal2PI],%%EDX      \n\t" /* (d*sinVal2PIPI[i]) d */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* cosVal2PI[i] (d*sinVal2PIPI[i]) d */
     "FSUBP                         \n\t" /* (cosVal2PI[i]-d*sinVal2PIPI[i]) d */
     "FMUL  %%st(1),%%st(0)         \n\t" /* (d*(cosVal2PI[i]-d*sinVal2PIPI[i])) d */

     "MOV   %[sinVal],%%EDX         \n\t" /* (d*(cosVal2PI[i]-d*sinVal2PIPI[i])) d */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* sinVal[i] (d*(cosVal2PI[i]-d*sinVal2PIPI[i])) d */
     "FADDP                         \n\t" /* (sinVal[i]+d*(cosVal2PI[i]-d*sinVal2PIPI[i])) d */
     "FSTPS %[sinv]                 \n\t" /* d */

     /* the same for cos value, this time popping d */
     "MOV   %[cosVal2PIPI],%%EDX    \n\t" /* cosVal2PIPI[i] d */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* cosVal2PIPI[i] d */
     "FMUL  %%st(1),%%st(0)         \n\t" /* (d*cosVal2PIPI[i]) d */

     "MOV   %[sinVal2PI],%%EDX      \n\t" /* sinVal2PI[i] (d*cosVal2PIPI[i]) d */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* sinVal2PI[i] (d*cosVal2PIPI[i]) d */
     "FSUBP                         \n\t" /* (sinVal2PI[i]-d*cosVal2PIPI[i]) d */
     "FMULP                         \n\t" /* (d*(sinVal2PI[i]-d*cosVal2PIPI[i])) */

     "MOV   %[cosVal],%%EDX         \n\t" /* cosVal[i] (d*(sinVal2PI[i]-d*cosVal2PIPI[i])) */
     "FLDL  (%%EDX,%%EAX)           \n\t" /* cosVal[i] (d*(sinVal2PI[i]-d*cosVal2PIPI[i])) */
     "FSUBP                         \n\t" /* (cosVal[i]-d*(sinVal2PI[i]-d*cosVal2PIPI[i])) */

     /* special here: tcos -= 1.0 */
     "FLD1                          \n\t" /* 1 (cosVal[i]-d*(sinVal2PI[i]-d*cosVal2PIPI[i])) */
     "FSUBRP                        \n\t" /* (cosVal[i]-d*(sinVal2PI[i]-d*cosVal2PIPI[i])-1) */

     "FSTPS %[cosv]                 \n\t" /* % */

     : /* output */
     [sinv]  "=m" (tsin), /* used for passing a constant, too */
     [cosv]  "=m" (tcos)

     : /* input */
     [x]           "m" (tempFreq0),
     [lutr]        "m" (lutr),
     [sinVal]      "m" (sinVal),
     [cosVal]      "m" (cosVal),
     [sinVal2PI]   "m" (sinVal2PI),
     [cosVal2PI]   "m" (cosVal2PI),
     [sinVal2PIPI] "m" (sinVal2PIPI),
     [cosVal2PIPI] "m" (cosVal2PIPI),
     [divLUTtab]   "m" (divTab)

     : /* clobbered registers */
     "eax", "edx", "st","st(1)","st(2)"
     );	

#else
  
  {
    UINT4 idx  = tempFreq0 * LUT_RES +.5;
    REAL8 d    = tempFreq0 - divTab[idx];
    REAL8 d2   = d*d;

    tsin = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
    tcos = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
    
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
      REAL8 d    = yRem - divTab[idx];
      REAL8 d2   = d*d;
      imagQ = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
      realQ = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
      
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

  static REAL8 sinTab[LUT_RES+(LUT_RES/4)+1]; /* Lookup tables for fast sin/cos calculation */
  static REAL8 sinTab2PI[LUT_RES+(LUT_RES/4)+1];
  static REAL8 sinTab2PIPI[LUT_RES+(LUT_RES/4)+1];
  static REAL8 divLUTtab[LUT_RES+1];
  static BOOLEAN firstCall = 1;


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
      for (k=0; k <= (LUT_RES/4)*5; k++) {
	sinTab[k] = sin((LAL_TWOPI*k)/(LUT_RES));
	sinTab2PI[k] = sinTab[k]  *  LAL_TWOPI;
	sinTab2PIPI[k] = sinTab2PI[k] * LAL_PI;
      }
      sinVal = sinTab;
      sinVal2PI = sinTab2PI;
      sinVal2PIPI = sinTab2PIPI;
      cosVal = sinVal+(LUT_RES/4);
      cosVal2PI = sinVal2PI+(LUT_RES/4);
      cosVal2PIPI = sinVal2PIPI+(LUT_RES/4);

      for (k=0; k <= LUT_RES; k++)
	divLUTtab[k] = (REAL8)k/(REAL8)(LUT_RES);
      firstCall = 0;
      divTab = divLUTtab;
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
