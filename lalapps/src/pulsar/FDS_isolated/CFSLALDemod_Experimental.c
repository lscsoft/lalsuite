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

#if defined(EXTERNAL_INNER_LOOP)
extern void cfs_inner_loop(UINT4 klim, REAL8 tempFreq1, REAL4 tsin, REAL4 tcos, COMPLEX8*Xalpha_k, REAL8*realXP, REAL8*imagXP);
#endif

#if defined(LOCAL_INNER_LOOP)
void cfs_inner_loop(UINT4 klim, REAL8 tempFreq1, REAL4 tsin, REAL4 tcos, COMPLEX8*Xalpha_k, REAL8*realXP, REAL8*imagXP) {
  UINT4 k;
  for(k=0; k < klim ; k++) {
    REAL8 realP, imagP;
    REAL4 xinv = (REAL4)OOTWOPI / (REAL4)tempFreq1;
    COMPLEX8 Xa = *Xalpha_k;

    Xalpha_k ++;
    tempFreq1 --;
    
    realP = tsin * xinv;
    imagP = tcos * xinv;

    /* these lines compute P*xtilde */
    *realXP += Xa.re * realP - Xa.im * imagP;
    *imagXP += Xa.re * imagP + Xa.im * realP;
    
  } /* for k < klim */
}
#endif

#ifdef USE_R4LALDEMOD

void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 

{ 

/* ---------------------------------------------------------------------- */
/* new test-version of LAL-Demod using single precision  */

  INT4 alpha,i;                 /* loop indices */
  REAL8 *xSum=NULL, *ySum=NULL; /* temp variables for computation of fs*as and fs*bs */
  INT4 s;                       /* local variable for spinDwn calcs. */
  REAL8 xTemp;                  /* temp variable for phase model */
  REAL8 deltaF;                 /* width of SFT band */
  INT4  k1;                     /* defining the sum over which is calculated */
  UINT4 k=0;
  REAL8 *skyConst;              /* vector of sky constants data */
  REAL8 *spinDwn;               /* vector of spinDwn parameters (maybe a structure? */
  INT4  spOrder;                /* maximum spinDwn order */
  REAL8 realXP, imagXP;         /* temp variables used in computation of */
  INT4  nDeltaF;                /* number of frequency bins per SFT band */
  INT4  sftIndex;               /* more temp variables */
  REAL8 realQ, imagQ;
  INT4 *tempInt1;
  UINT4 index;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
  UINT4 klim = 2*params->Dterms;
  REAL8 f;
  static REAL8 sinVal[LUT_RES+1], cosVal[LUT_RES+1];        /*LUT values computed by the routine do_trig_lut*/
  static BOOLEAN firstCall = 1;

  REAL8 A=params->amcoe->A;
  REAL8 B=params->amcoe->B;
  REAL8 C=params->amcoe->C;
  REAL8 D=params->amcoe->D;

  UINT4 M=params->SFTno;

  INITSTATUS( status, "TestLALDemod2", rcsid );

  /* catch some obvious programming errors */
  ASSERT ( (Fs != NULL)&&(Fs->F != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  if (params->returnFaFb)
    {
      ASSERT ( (Fs->Fa != NULL)&&(Fs->Fb != NULL), status, 
	       COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
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
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( (LAL_TWOPI*k)/LUT_RES );
          cosVal[k] = cos( (LAL_TWOPI*k)/LUT_RES );
        }
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
        REAL4 tsin, tcos;
        COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
        REAL4 a = params->amcoe->a->data[alpha];
        REAL4 b = params->amcoe->b->data[alpha];
        REAL8 x,y;
        REAL4 realP, imagP;             /* real and imaginary parts of P, see CVS */

        /* NOTE: sky-constants are always positive!!
         * this can be seen from there definition (-> documentation)
         * we will use this fact in the following! 
         */
        xTemp= f * skyConst[ tempInt1[ alpha ] ] + xSum[ alpha ];       /* >= 0 !! */
        
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

        /* find correct index into LUT -- pick closest point */
        tempFreq0 = xTemp - (UINT4)xTemp;  /* lies in [0, +1) by definition */

        index = (UINT4)( tempFreq0 * LUT_RES + 0.5 );   /* positive! */
        {
          REAL8 d=LAL_TWOPI*(tempFreq0 - (REAL8)index/(REAL8)LUT_RES);
          REAL8 d2=0.5*d*d;
          REAL8 ts=sinVal[index];
          REAL8 tc=cosVal[index];
                
          tsin = ts+d*tc-d2*ts;
          tcos = tc-d*ts-d2*tc-1.0;
        }

        y = - LAL_TWOPI * ( f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ] );
        realQ = cos(y);
        imagQ = sin(y);

        /*
        REAL8 yTemp = f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ];
        REAL8 yRem = yTemp - (UINT4)yTemp;

        index = (UINT4)( yRem * LUT_RES + 0.5 );
        {
          REAL8 d = LAL_TWOPI*(yRem - (REAL8)index/(REAL8)LUT_RES);
          REAL8 d2=0.5*d*d;
          REAL8 ts = sinVal[index];
          REAL8 tc = cosVal[index];
                
          imagQ = ts + d * tc - d2 * ts;
          imagQ = -imagQ;
          realQ = tc - d * ts - d2 * tc;
        }
        */

        k1 = (UINT4)xTemp - params->Dterms + 1;

        sftIndex = k1 - params->ifmin;

	if(sftIndex < 0){
              fprintf(stderr,"ERROR! sftIndex = %d < 0 in TestLALDemod run %d\n", sftIndex, cfsRunNo);
              fprintf(stderr," alpha=%d, k1=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                      alpha, k1, xTemp, params->Dterms, params->ifmin);
	      ABORT(status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
	}

        tempFreq1 = tempFreq0 + params->Dterms - 1;     /* positive if Dterms > 1 (trivial) */

        x = LAL_TWOPI * tempFreq1;      /* positive! */

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
          {
            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;

            realXP=0.0;
            imagXP=0.0;

            /* Loop over terms in dirichlet Kernel */
            for(k=0; k < klim ; k++)
              {
                REAL4 xinv = (REAL4)OOTWOPI / (REAL4)tempFreq1;
                COMPLEX8 Xa = *Xalpha_k;
                Xalpha_k ++;
                tempFreq1 --;
                
                realP = tsin * xinv;
                imagP = tcos * xinv;
                /* these four lines compute P*xtilde */
                realXP += Xa.re * realP - Xa.im * imagP;
                imagXP += Xa.re * imagP + Xa.im * realP;

              } /* for k < klim */

          } /* if x cannot be close to 0 */
        

        if(sftIndex-1 > maxSFTindex) {
          fprintf(stderr,"ERROR! sftIndex = %d > %d in TestLALDemod\nalpha=%d,"
                 "k1=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                 sftIndex-1, maxSFTindex, alpha, k1, xTemp, params->Dterms, params->ifmin);
	  ABORT(status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
	}

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
#else /* USE_R4LALDEMOD, USE_ALTIVEC */


void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 

{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8 *xSum=NULL, *ySum=NULL; /* temp variables for computation of fs*as and fs*bs */
  INT4 s;                       /* local variable for spinDwn calcs. */
  REAL8 xTemp;                  /* temp variable for phase model */
  REAL8 deltaF;                 /* width of SFT band */
  INT4  k1;                     /* defining the sum over which is calculated */
  UINT4 k=0;
  REAL8 *skyConst;              /* vector of sky constants data */
  REAL8 *spinDwn;               /* vector of spinDwn parameters (maybe a structure? */
  INT4  spOrder;                /* maximum spinDwn order */
  REAL8 realXP, imagXP;         /* temp variables used in computation of */
  INT4  nDeltaF;                /* number of frequency bins per SFT band */
  INT4  sftIndex;               /* more temp variables */
  REAL8 realQ, imagQ;
  INT4 *tempInt1;
  UINT4 index;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
  UINT4 klim = 2*params->Dterms;
  REAL8 f;
  static REAL8 sinVal[LUT_RES+1], cosVal[LUT_RES+1];        /*LUT values computed by the routine do_trig_lut*/
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
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( (LAL_TWOPI*k)/LUT_RES );
          cosVal[k] = cos( (LAL_TWOPI*k)/LUT_RES );
        }
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
        REAL4 tsin, tcos;
        COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
        REAL4 a = params->amcoe->a->data[alpha];
        REAL4 b = params->amcoe->b->data[alpha];
        REAL8 x,y;
        REAL4 realP, imagP;             /* real and imaginary parts of P, see CVS */

        /* NOTE: sky-constants are always positive!!
         * this can be seen from there definition (-> documentation)
         * we will use this fact in the following! 
         */
        xTemp= f * skyConst[ tempInt1[ alpha ] ] + xSum[ alpha ];       /* >= 0 !! */
        
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

        /* find correct index into LUT -- pick closest point */
        tempFreq0 = xTemp - (UINT4)xTemp;  /* lies in [0, +1) by definition */

        index = (UINT4)( tempFreq0 * LUT_RES + 0.5 );   /* positive! */
        {
          REAL8 d=LAL_TWOPI*(tempFreq0 - (REAL8)index/(REAL8)LUT_RES);
          REAL8 d2=0.5*d*d;
          REAL8 ts=sinVal[index];
          REAL8 tc=cosVal[index];
                
          tsin = ts+d*tc-d2*ts;
          tcos = tc-d*ts-d2*tc-1.0;
        }

        y = - LAL_TWOPI * ( f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ] );
        realQ = cos(y);
        imagQ = sin(y);

        /*
        REAL8 yTemp = f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ];
        REAL8 yRem = yTemp - (UINT4)yTemp;

        index = (UINT4)( yRem * LUT_RES + 0.5 );
        {
          REAL8 d = LAL_TWOPI*(yRem - (REAL8)index/(REAL8)LUT_RES);
          REAL8 d2=0.5*d*d;
          REAL8 ts = sinVal[index];
          REAL8 tc = cosVal[index];
                
          imagQ = ts + d * tc - d2 * ts;
          imagQ = -imagQ;
          realQ = tc - d * ts - d2 * tc;
        }
        */

        k1 = (UINT4)xTemp - params->Dterms + 1;

        sftIndex = k1 - params->ifmin;

	if(sftIndex < 0){
              fprintf(stderr,"ERROR! sftIndex = %d < 0 in TestLALDemod run %d\n", sftIndex, cfsRunNo);
              fprintf(stderr," alpha=%d, k1=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                      alpha, k1, xTemp, params->Dterms, params->ifmin);
	      ABORT(status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
	}

        tempFreq1 = tempFreq0 + params->Dterms - 1;     /* positive if Dterms > 1 (trivial) */

        x = LAL_TWOPI * tempFreq1;      /* positive! */

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
#if defined(EXTERNAL_INNER_LOOP)
	    cfs_inner_loop(klim, tempFreq1, tsin, tcos, Xalpha+sftIndex, &realXP, &imagXP);
#elif defined(USE_SSE2)
          {
            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
	    
	    /* scalar variables to hold partial sums */
	    float realSums[4];
	    float imagSums[4];
	
	    __m128 realXP_v, imagXP_v;
	
	    /* init partial sums in SSE registers to zero */
	    realXP_v = _mm_setzero_ps();
	    imagXP_v = _mm_setzero_ps();
	    
            realXP=0.0;
            imagXP=0.0;

            /* Loop over terms in dirichlet Kernel */
            for(k=0; k+3 < klim ; k+=4)
              {
		/* subtract index for next 4 terms */
		__m128 tempFreq_v = _mm_set_ps((REAL4)(tempFreq1 - 3.0), (REAL4)(tempFreq1 - 2.0), (REAL4)(tempFreq1 - 1.0), (REAL4)(tempFreq1 - 0.0));
		__m128 xinv_v = _mm_div_ps(_mm_set1_ps(OOTWOPI), tempFreq_v);
		/* load 4 pairs of re/im */
		__m128 xAlpha0 = _mm_loadu_ps((float *)&Xalpha_k[0]);
		__m128 xAlpha1 = _mm_loadu_ps((float *)&Xalpha_k[2]);

                tempFreq1 = tempFreq1 - 4.0;
		Xalpha_k += 4;

		__m128 realP_v = _mm_mul_ps(_mm_set1_ps(tsin), xinv_v);
		__m128 imagP_v = _mm_mul_ps(_mm_set1_ps(tcos), xinv_v);
		
		/* de-interleave real and imag data -- this will create a vector of reals and a vector of imag */
		__m128 Xa_re = _mm_shuffle_ps(xAlpha0, xAlpha1, _MM_SHUFFLE(2,0,2,0));
		__m128 Xa_im = _mm_shuffle_ps(xAlpha0, xAlpha1, _MM_SHUFFLE(3,1,3,1));
		
		/* multiply and add to sums */
		realXP_v = _mm_add_ps( realXP_v, _mm_sub_ps( _mm_mul_ps(Xa_re, realP_v), _mm_mul_ps(Xa_im, imagP_v)));
		imagXP_v = _mm_add_ps( imagXP_v, _mm_add_ps( _mm_mul_ps(Xa_re, imagP_v), _mm_mul_ps(Xa_im, realP_v)));
              } /* for k < klim */

	    _mm_storeu_ps(realSums, realXP_v);
	    _mm_storeu_ps(imagSums, imagXP_v);
	    /* add up partial sums from SSE registers using scalar math */
	    realXP = realSums[0] + realSums[1] + realSums[2] + realSums[3];
	    imagXP = imagSums[0] + imagSums[1] + imagSums[2] + imagSums[3];
	    
	    /* continue on to catch any leftover terms */
	    for( ; k<klim; k++) {
                REAL4 xinv = (REAL4)OOTWOPI / (REAL4)tempFreq1;
                COMPLEX8 Xa = *Xalpha_k;
                Xalpha_k ++;
                tempFreq1 --;
                
                realP = tsin * xinv;
                imagP = tcos * xinv;
                /* these four lines compute P*xtilde */
                realXP += Xa.re * realP - Xa.im * imagP;
                imagXP += Xa.re * imagP + Xa.im * realP;

	    }

          }
#else
          {
            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;

            realXP=0.0;
            imagXP=0.0;

            /* Loop over terms in Dirichlet Kernel */

#ifdef SPLITVUNITS
/* The following are variants of the main loop split up for CPU SIMD units such as SSE.
   They are switched "on" by defining the macro SPLITVUNITS. They all have a special
   treatment for groups of VUNITS (4 by default) items in the Dirichlet Kernel, and then
   process the remaining ones with the usual generic loop.
*/

#ifndef VUNITS
#define VUNITS 4
#endif
            for(k=0; k+VUNITS < klim; k+=VUNITS) {
		int ki;
		float realXP_v[VUNITS];
		float imagXP_v[VUNITS];

#ifdef VUNITS_TWOLOOPS
/* This version splits up the main loop into two loops (counting the reduction below)
   It has shown to be the fastest version when compiling without auto-vectorization,
   but with SSE(2) support, e.g. with MSC 7.1
*/
		float xinv;
		float Xar,Xai;

		/* inner loop */
		for(ki=0;ki<VUNITS;ki++) {
		    COMPLEX8 Xa = *Xalpha_k;
		    Xar = Xa.re;
		    Xai = Xa.im;
		    xinv = (float)OOTWOPI / (float)(tempFreq1 - ki);
		    realP = tsin * xinv;
		    imagP = tcos * xinv;
		    realXP_v[ki] = Xar * realP - Xai * imagP;
		    imagXP_v[ki] = Xar * imagP + Xai * realP;
		    Xalpha_k ++;
		}

#else /* VUNITS_TWOLOOPS */
/* This version splits up the main loop into three loops (counting the reduction below)
   It's the only version (so far) that convinces the icc 8.1 to actually vectorize the
   "inner" loop.
*/
		float xinv_v[VUNITS];
		float Xar_v[VUNITS];
		float Xai_v[VUNITS];

		/* prepare Xa */
		for(ki=0;ki<VUNITS;ki++) {
		    COMPLEX8 Xa = *Xalpha_k;
		    Xar_v[ki] = Xa.re;
		    Xai_v[ki] = Xa.im;
		    Xalpha_k ++;
		}

		/* inner loop */
		for(ki=0;ki<VUNITS;ki++) {
		    xinv_v[ki] = (float)OOTWOPI / (float)(tempFreq1 - ki);
		    realP = tsin * xinv_v[ki];
		    imagP = tcos * xinv_v[ki];
		    realXP_v[ki] = Xar_v[ki] * realP - Xai_v[ki] * imagP;
		    imagXP_v[ki] = Xar_v[ki] * imagP + Xai_v[ki] * realP;
		}

#endif /* VUNITS_TWOLOOPS */

		/* reduction */
		for(ki=0;ki<VUNITS;ki++) {
		    realXP += realXP_v[ki];
		    imagXP += imagXP_v[ki];
		}

		/* increment */
		tempFreq1 -= VUNITS;
	    }
            for(; k < klim; k++)
#else /* SPLITVUNITS */

            for(k=0; k < klim ; k++)
#endif /* SPLITVUNITS */
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

        if(sftIndex-1 > maxSFTindex) {
          fprintf(stderr,"ERROR! sftIndex = %d > %d in TestLALDemod\nalpha=%d,"
                 "k1=%d, xTemp=%20.17f, Dterms=%d, ifmin=%d\n",
                 sftIndex-1, maxSFTindex, alpha, k1, xTemp, params->Dterms, params->ifmin);
	  ABORT(status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
	}


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
#endif /* USE_R4LALDEMOD */

