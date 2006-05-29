/* LALDemod variant with structs and signatures for modifications by Feket Akos
 * Authors see ComputeFStatistic.c
                                                         Bernd Machenschalk */
RCSID( "$Id$");

__declspec(align(64)) static char tabspace[10*1024];

__declspec(align(64)) static struct marked_struct {
  REAL8 xTemp,yTemp;
  REAL8 realQ,imagQ;  
  REAL8 tempFreq0;
  REAL4 tsin,tcos;
  REAL4 XRes, XIms;
  COMPLEX8* Xalpha_k;
} marked;  

#include "2500nops.h"

/* <lalVerbatim file="LALDemodCP"> */
void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 
/* </lalVerbatim> */
{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8 *xSum=NULL, *ySum=NULL; /* temp variables for computation of fs*as and fs*bs */
  INT4 s;                       /* local variable for spinDwn calcs. */
  REAL8 deltaF;                 /* width of SFT band */
  INT4  k1;                     /* defining the sum over which is calculated */
  UINT4 k=0;
  REAL8 *skyConst;              /* vector of sky constants data */
  REAL8 *spinDwn;               /* vector of spinDwn parameters (maybe a structure? */
  INT4  spOrder;                /* maximum spinDwn order */
  REAL8 realXP, imagXP;         /* temp variables used in computation of */
  INT4  nDeltaF;                /* number of frequency bins per SFT band */
  INT4  sftIndex;               /* more temp variables */
  INT4 *tempInt1;
  /* UINT4 index; */
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
#ifdef USE_BOINC
#define klim 32
#else
  UINT4 klim = 2*params->Dterms;
#endif
  REAL8 f;
  static REAL8 sinVal[LUT_RES+(LUT_RES/4)+1]; /* Lookup tables for fast sin/cos calculation */
  static REAL8 *cosVal;
#ifdef USE_4_LUT
  static REAL8 sinVal2PI[LUT_RES+(LUT_RES/4)+1];
  static REAL8 sinVal2PIPI[LUT_RES+(LUT_RES/4)+1];
  static REAL8 *cosVal2PI, *cosVal2PIPI;
  static REAL8 divLUTtab[LUT_RES+1];
#endif
  static BOOLEAN firstCall = 1;

  REAL8 A=params->amcoe->A;
  REAL8 B=params->amcoe->B;
  REAL8 C=params->amcoe->C;
  REAL8 D=params->amcoe->D;

  UINT4 M=params->SFTno;

  INITSTATUS( status, "TestLALDemod", rcsid );

  /* catch some obvious programming errors */
  ASSERT ( (Fs != NULL)&&(Fs->F != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  if (params->returnFaFb)
    {
      ASSERT ( (Fs->Fa != NULL)&&(Fs->Fb != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
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
        sinVal[k] = sin((LAL_TWOPI*k)/(LUT_RES));
#ifdef USE_4_LUT
        sinVal2PI[k] = sinVal[k]  *  LAL_TWOPI;
        sinVal2PIPI[k] = sinVal2PI[k] * LAL_PI;
#endif
      }
      cosVal = sinVal+(LUT_RES/4);
#ifdef USE_4_LUT
      cosVal2PI = sinVal2PI+(LUT_RES/4);
      cosVal2PIPI = sinVal2PIPI+(LUT_RES/4);

      for (k=0; k <= LUT_RES; k++)
        divLUTtab[k] = (REAL8)k/(REAL8)(LUT_RES);
#endif      
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
        REAL8 tempFreq1;
        COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
        REAL4 a = params->amcoe->a->data[alpha];
        REAL4 b = params->amcoe->b->data[alpha];
        REAL8 x;
#ifndef USE_LUT_Y
        REAL8 y;
#endif
        REAL4 realP, imagP;             /* real and imaginary parts of P, see CVS */

        /* NOTE: sky-constants are always positive!!
         * this can be seen from there definition (-> documentation)
         * we will use this fact in the following! 
         */
        marked.xTemp= f * skyConst[ tempInt1[ alpha ] ] + xSum[ alpha ];       /* >= 0 !! */
        
        /* this will now be assumed positive, but we double-check this to be sure */
        if  (!finite(marked.xTemp)) {
            fprintf (stderr, "xTemp is not finite\n");
            fprintf (stderr, "DEBUG: loop=%d, xTemp=%f, f=%f, alpha=%d, tempInt1[alpha]=%d\n", 
                     i, marked.xTemp, f, alpha, tempInt1[alpha]);
            fprintf (stderr, "DEBUG: skyConst[ tempInt1[ alpha ] ] = %f, xSum[ alpha ]=%f\n",
                     skyConst[ tempInt1[ alpha ] ], xSum[ alpha ]);
#ifndef USE_BOINC
            fprintf (stderr, "\n*** PLEASE report this bug to pulgroup@gravity.phys.uwm.edu *** \n\n");
#endif
            exit (COMPUTEFSTAT_EXIT_DEMOD);
        }
        if (marked.xTemp < 0) {
            fprintf (stderr, "xTemp >= 0 failed\n");
            fprintf (stderr, "DEBUG: loop=%d, xTemp=%f, f=%f, alpha=%d, tempInt1[alpha]=%d\n", 
                     i, marked.xTemp, f, alpha, tempInt1[alpha]);
            fprintf (stderr, "DEBUG: skyConst[ tempInt1[ alpha ] ] = %f, xSum[ alpha ]=%f\n",
                     skyConst[ tempInt1[ alpha ] ], xSum[ alpha ]);
#ifndef USE_BOINC
            fprintf (stderr, "\n*** PLEASE report this bug to pulgroup@gravity.phys.uwm.edu *** \n\n");
#endif
            exit (COMPUTEFSTAT_EXIT_DEMOD);
        }

        /* find correct index into LUT -- pick closest point */
        marked.tempFreq0 = marked.xTemp - (UINT4)(marked.xTemp);   /* lies in [0, +1) by definition */
#ifdef USE_LUT_Y
        /* use LUT here, too */
        marked.yTemp = f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ];
#endif

        MARK(1);

        {
          UINT4  idx  = marked.tempFreq0*(REAL8)(LUT_RES)+.5;
#ifdef USE_4_LUT
          REAL8 d    = marked.tempFreq0 - divLUTtab[idx];
          REAL8 d2   = d*d;
                
          marked.tsin = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
          marked.tcos = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
#else
          REAL8 d=LAL_TWOPI*(marked.tempFreq0-(REAL8)idx/(REAL8)LUT_RES);
          REAL8 d2=0.5*d*d;
          REAL8 ts=sinVal[idx];
          REAL8 tc=cosVal[idx];
                
          marked.tsin = ts+d*tc-d2*ts;
          marked.tcos = tc-d*ts-d2*tc;
#endif
          marked.tcos -= 1.0;
        }

#ifdef USE_LUT_Y
        {
          REAL8 yRem = marked.yTemp - (INT4)(marked.yTemp);
          if (yRem < 0) { yRem += 1.0f; } /* make sure this is in [0..1) */
          {
#ifdef USE_4_LUT
            UINT4 idx  = yRem*(REAL8)(LUT_RES)+.5;
            REAL8 d    = yRem-divLUTtab[idx];
            REAL8 d2   = d*d;
          
            marked.imagQ = sinVal[idx] + d * cosVal2PI[idx] - d2 * sinVal2PIPI[idx];
            marked.realQ = cosVal[idx] - d * sinVal2PI[idx] - d2 * cosVal2PIPI[idx];
#else
            UINT4 idx = (UINT4)( yRem * LUT_RES + 0.5 );
            REAL8 d = LAL_TWOPI*(yRem - (REAL8)idx/(REAL8)LUT_RES);
            REAL8 d2=0.5*d*d;
            REAL8 ts = sinVal[idx];
            REAL8 tc = cosVal[idx];
            
            marked.imagQ = ts + d * tc - d2 * ts;
            marked.realQ = tc - d * ts - d2 * tc;
#endif
            marked.imagQ = -marked.imagQ;
          }
        }
#else
        marked.yTemp = - LAL_TWOPI * ( f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ] );
        marked.realQ = cos(marked.yTemp);
        marked.imagQ = sin(marked.yTemp);
#endif


        MARK(2);

        k1 = (UINT4)marked.xTemp - params->Dterms + 1;

        sftIndex = k1 - params->ifmin;

        if(sftIndex < 0){
              fprintf(stderr,"ERROR! sftIndex = %d < 0 in TestLALDemod run %d\n", sftIndex, cfsRunNo);
              fprintf(stderr," alpha=%d, k1=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                      alpha, k1, marked.xTemp, params->Dterms, params->ifmin);
              ABORT(status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
        }

        tempFreq1 = marked.tempFreq0 + params->Dterms - 1;     /* positive if Dterms > 1 (trivial) */

        x = LAL_TWOPI * tempFreq1;      /* positive! */

        /* we branch now (instead of inside the central loop)
         * depending on wether x can ever become SMALL in the loop or not, 
         * because it requires special treatment in the Dirichlet kernel
         */
        if ( marked.tempFreq0 < LD_SMALL ) 
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
                    realP = marked.tsin / x;
                    imagP = marked.tcos / x;
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

          {
            typedef struct R44Vtag {
              REAL4 a,b,c,d;
            } R44V;

            NRCSID (CFSLOOPx86MAS, "$Id$");

            REAL4 XResX=0.0, XImsX=0.0;   /* sums of Xa.re and Xa.im */
            REAL4 accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */

            /* prepare values for SSE */
            REAL4 tempFreqX;
            COMPLEX8 *Xalpha_kX;
            R44V V0011,V2222;

            MARK(3);

            V0011.a = 0.0;
            V0011.b = 0.0;
            V0011.c = 1.0;
            V0011.d = 1.0;

            V2222.a = 2.0;
            V2222.b = 2.0;
            V2222.c = 2.0;
            V2222.d = 2.0;

            marked.XRes=0.0; marked.XIms=0.0;   /* sums of Xa.re and Xa.im */
            tempFreqX = tempFreq1; /* REAL4 because of SSE */
            tempFreq1 = tempFreq1 - 14;
            marked.Xalpha_k = Xalpha + sftIndex;
            Xalpha_kX = marked.Xalpha_k; /* -> SSE values */
            marked.Xalpha_k = marked.Xalpha_k + 14; /* -> FPU values */
 
            /* let the compiler code the x87 part */
            for(k=0; k < 4 ; k++)
              {
                marked.XRes = tempFreq1 * marked.XRes + (*marked.Xalpha_k).re * accFreq;
                marked.XIms = tempFreq1 * marked.XIms + (*marked.Xalpha_k).im * accFreq;
                
                accFreq *= (REAL4)tempFreq1;
                tempFreq1 --;
                marked.Xalpha_k ++;
              } /* for k < klim */
            
            accFreq = 1.0 / accFreq;
            marked.XRes *= accFreq;
            marked.XIms *= accFreq;

            /* The SSE part is coded in Assembler */
            __asm {
               MOV    EBX,Xalpha_kX
               MOVSS  XMM0,tempFreqX

               MOVUPS XMM5,V0011
               SHUFPS XMM0,XMM0,0 /* XMM0: f   f   f   f */
               SUBPS  XMM0,XMM5 /* XMM0: f-1 f-1 f   f */
               XORPS  XMM2,XMM2 /* XMM2 will collect the low-precision values */
               MOVUPS XMM5,V2222

               /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
               /* one SSE register will consist 4 REAL4 values */
               /* 4 REAL4 vaules = 2 ReIm pairs */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */
 
               SUBPS  XMM0,XMM5 /* JUMP OVER FPU CALCULATED VALUES */
               ADD    EBX,32 /* Xalpha_kX = Xalpha_kX + 2; */
               SUBPS  XMM0,XMM5 /* JUMP OVER FPU CALCULATED VALUES */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               MOVLPS XMM4,[EBX]
               MOVHPS XMM4,[EBX+8]
               RCPPS  XMM1,XMM0 /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */
               SUBPS  XMM0,XMM5 /* XMM2: f-3 f-3 f-2 f-2 */
               MULPS  XMM1,XMM4 /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */
               ADD    EBX,16 /* Xalpha_kX = Xalpha_kX + 2; */
               ADDPS  XMM2,XMM1 /* XMM2: C_ImH C_ReH C_ImL C_ReL */

               /* XMM2 consists the low precision values */
               /* XRes and Xims consist the high precision values */
 
               MOVHLPS XMM3,XMM2 /* XMM3: ? ? C_ImH C_ReH */
               ADDPS   XMM2,XMM3 /* XMM2: - - C_Im C_Re */

               MOVSS   XResX,XMM2 /* SAVE Re part */
               SHUFPS  XMM2,XMM2,1 /* XMM0: f   f   f   f */
               MOVSS   XImsX,XMM2 /* SAVE Im part */
                 }

            /* And last, we add the single and double precision values */
            
            marked.XRes = marked.XRes + XResX;
            marked.XIms = marked.XIms + XImsX;

            MARK(4);

            {
              REAL4 tsin2pi = marked.tsin * (REAL4)OOTWOPI;
              REAL4 tcos2pi = marked.tcos * (REAL4)OOTWOPI;

              realXP = tsin2pi * marked.XRes - tcos2pi * marked.XIms;
              imagXP = tcos2pi * marked.XRes + tsin2pi * marked.XIms;
            }
          } /* if x cannot be close to 0 */


        if(sftIndex-1 > maxSFTindex) {
          fprintf(stderr,"ERROR! sftIndex = %d > %d in TestLALDemod\nalpha=%d,"
                 "k1=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                 sftIndex-1, maxSFTindex, alpha, k1, marked.xTemp, params->Dterms, params->ifmin);
          ABORT(status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
        }


        /* implementation of amplitude demodulation */
        {
          REAL8 realQXP = realXP*marked.realQ-imagXP*marked.imagQ;
          REAL8 imagQXP = realXP*marked.imagQ+imagXP*marked.realQ;
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

