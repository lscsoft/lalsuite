static inline void gc_hotloop (REAL4 * fgrid2F, REAL4 * cgrid2F, UCHAR * fgridnc, REAL4 TwoFthreshold, UINT4 length  ) __attribute__ ((hot));
static inline void gc_hotloop_no_nc (REAL4 * fgrid2F, REAL4 * cgrid2F, UINT4 length  ) __attribute__ ((hot));
static inline void gc_hotloop_2Fmax_tracking (REAL4 * fgrid2F, REAL4 * fgrid2Fmax, UINT4 * fgrid2FmaxIdx, REAL4 * cgrid2F, UINT4 k, UINT4 length  ) __attribute__ ((hot));




#ifdef __APPLE__

/* Apple's gcc aligns ok */
#define ALRealloc LALRealloc
#define ALFree LALFree

#elif defined (__MINGW32__)

extern void *__mingw_aligned_realloc(void *ptr, size_t size, size_t align);
#define ALRealloc(p,s) __mingw_aligned_realloc(p,s,16)
#define ALFree __mingw_aligned_free

#else // neither APPLE nor MinGW

#include <stdlib.h>

#define ALFree free

static void *ALRealloc(void *ptr, size_t size);

/* in our case there is no need to keep the data,
   so we can simply do a free() and malloc() */
void *ALRealloc(void *ptr, size_t size) {
  if(ptr)
    free(ptr);
  if(posix_memalign(&ptr, 16, size))
    return(NULL);
  return(ptr);
}

#endif

void gc_hotloop_2Fmax_tracking (REAL4 * fgrid2F, REAL4 * fgrid2Fmax, UINT4 * fgrid2FmaxIdx, REAL4 * cgrid2F, UINT4 k, UINT4 length  ) {

  UINT4 ifreq_fg;
  int newMax;

  UINT4 VIIII[4] __attribute__ ((aligned (16))) = { k,k,k,k };
  UINT4 V1111[4] __attribute__ ((aligned (16))) = { 0xffffffff,0xffffffff,0xffffffff,0xffffffff};



  /* if this is the first segment (seg 0), then all we need to do is copy the cg vector into the fg 
   * vector and initialize the loudest segment data to also the first segment. This assumes 
   * we are calling this function in ascending order of segments */

  if(k==0) {
    memcpy(fgrid2F,cgrid2F,sizeof(REAL4)*length);
    memcpy(fgrid2Fmax,cgrid2F,sizeof(REAL4)*length);
    memset(fgrid2FmaxIdx,0,sizeof(UINT4)*length);

    return;
  }



  for(ifreq_fg=0 ; ifreq_fg +16 < length; ifreq_fg+=16 ) {
#ifdef EXP_NO_ASM

#pragma ivdep
    for(int j=0 ; j < 16; j++ ) {
            fgrid2F[0] += cgrid2F[0] ;

            newMax=(cgrid2F[0] >= fgrid2Fmax[0]);
            fgrid2Fmax[0]=fmaxf(fgrid2Fmax[0],cgrid2F[0]);
            fgrid2FmaxIdx[0]=fgrid2FmaxIdx[0]*(1-newMax)+k*newMax;
            fgrid2F++;
            cgrid2F++;
            fgrid2Fmax++;
            fgrid2FmaxIdx++;
    }

#else


    __asm __volatile (
         "MOVAPS %[Vk],%%xmm0 \n\t"
         "MOVAPS %[V1],%%xmm1 \n\t"
/* iteration 0,1,2,3 */
         "MOVUPS  (%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  (%[fg2F]),%%xmm3 \n\t"
         "MOVAPS  (%[fg2Fmax]),%%xmm4 \n\t"
         "MOVAPS  (%[fg2FmaxIdx]),%%xmm5 \n\t"
/* create mask with comparison result of former max 2F */
         "MOVAPS %%xmm4,%%xmm7 \n\t"
	 "CMPLEPS %%xmm2,%%xmm7 \n\t"     /* -1 if previous 2Fmax is <= coarse grid value */
/* summing */
         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */
         "MOVAPS  %%xmm3,(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */
	 "MOVAPS %%xmm0,%%xmm6 \n\t"
         "ANDPS   %%xmm7,%%xmm2  \n\t"
	 "ANDPS   %%xmm7,%%xmm6  \n\t"

/* negate the bitmask */
         "XORPS   %%xmm1,%%xmm7 \n\t"
         "ANDPS   %%xmm7,%%xmm4  \n\t"
         "ANDPS   %%xmm7,%%xmm5  \n\t"
/*get the new entries for the max 2F values by ORing the masked values */
         "ORPS   %%xmm2,%%xmm4  \n\t"
         "ORPS   %%xmm6,%%xmm5  \n\t"
/* write back */
	 "MOVAPS  %%xmm4,(%[fg2Fmax]) \n\t"
         "MOVAPS  %%xmm5,(%[fg2FmaxIdx]) \n\t"



/* iteration 4,5,6,7 */ 


         "MOVUPS  0x10(%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  0x10(%[fg2F]),%%xmm3 \n\t"
         "MOVAPS  0x10(%[fg2Fmax]),%%xmm4 \n\t"
         "MOVAPS  0x10(%[fg2FmaxIdx]),%%xmm5 \n\t"
/* create mask with comparison result of former max 2F */
         "MOVAPS %%xmm4,%%xmm7 \n\t"
         "CMPLEPS %%xmm2,%%xmm7 \n\t"     /* -1 if previous 2Fmax is <= coarse grid value */
/* summing */
         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */
         "MOVAPS  %%xmm3,0x10(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */
         "MOVAPS %%xmm0,%%xmm6 \n\t"
         "ANDPS   %%xmm7,%%xmm2  \n\t"
         "ANDPS   %%xmm7,%%xmm6  \n\t"

/* negate the bitmask */
         "XORPS   %%xmm1,%%xmm7 \n\t"
         "ANDPS   %%xmm7,%%xmm4  \n\t"
         "ANDPS   %%xmm7,%%xmm5  \n\t"
/*get the new entries for the max 2F values by ORing the masked values */
         "ORPS   %%xmm2,%%xmm4  \n\t"
         "ORPS   %%xmm6,%%xmm5  \n\t"
/* write back */
         "MOVAPS  %%xmm4,0x10(%[fg2Fmax]) \n\t"
         "MOVAPS  %%xmm5,0x10(%[fg2FmaxIdx]) \n\t"

/* iteration 8,9,10,11 */ 


         "MOVUPS  0x20(%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  0x20(%[fg2F]),%%xmm3 \n\t"
         "MOVAPS  0x20(%[fg2Fmax]),%%xmm4 \n\t"
         "MOVAPS  0x20(%[fg2FmaxIdx]),%%xmm5 \n\t"
/* create mask with comparison result of former max 2F */
         "MOVAPS %%xmm4,%%xmm7 \n\t"
         "CMPLEPS %%xmm2,%%xmm7 \n\t"     /* -1 if previous 2Fmax is <= coarse grid value */
/* summing */
         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */
         "MOVAPS  %%xmm3,0x20(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */
         "MOVAPS %%xmm0,%%xmm6 \n\t"
         "ANDPS   %%xmm7,%%xmm2  \n\t"
         "ANDPS   %%xmm7,%%xmm6  \n\t"

/* negate the bitmask */
         "XORPS   %%xmm1,%%xmm7 \n\t"
         "ANDPS   %%xmm7,%%xmm4  \n\t"
         "ANDPS   %%xmm7,%%xmm5  \n\t"
/*get the new entries for the max 2F values by ORing the masked values */
         "ORPS   %%xmm2,%%xmm4  \n\t"
         "ORPS   %%xmm6,%%xmm5  \n\t"
/* write back */
         "MOVAPS  %%xmm4,0x20(%[fg2Fmax]) \n\t"
         "MOVAPS  %%xmm5,0x20(%[fg2FmaxIdx]) \n\t"


/* iteration 12,13,14,15 */ 


         "MOVUPS  0x30(%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  0x30(%[fg2F]),%%xmm3 \n\t"
         "MOVAPS  0x30(%[fg2Fmax]),%%xmm4 \n\t"
         "MOVAPS  0x30(%[fg2FmaxIdx]),%%xmm5 \n\t"
/* create mask with comparison result of former max 2F */
         "MOVAPS %%xmm4,%%xmm7 \n\t"
         "CMPLEPS %%xmm2,%%xmm7 \n\t"     /* -1 if previous 2Fmax is <= coarse grid value */
/* summing */
         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */
         "MOVAPS  %%xmm3,0x30(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */
         "MOVAPS %%xmm0,%%xmm6 \n\t"
         "ANDPS   %%xmm7,%%xmm2  \n\t"
         "ANDPS   %%xmm7,%%xmm6  \n\t"

/* negate the bitmask */
         "XORPS   %%xmm1,%%xmm7 \n\t"
         "ANDPS   %%xmm7,%%xmm4  \n\t"
         "ANDPS   %%xmm7,%%xmm5  \n\t"
/*get the new entries for the max 2F values by ORing the masked values */
         "ORPS   %%xmm2,%%xmm4  \n\t"
         "ORPS   %%xmm6,%%xmm5  \n\t"
/* write back */
         "MOVAPS  %%xmm4,0x30(%[fg2Fmax]) \n\t"
         "MOVAPS  %%xmm5,0x30(%[fg2FmaxIdx]) \n\t"



/*  ---------------------------------------------------*/
:
      /* output */

      :
      /* input */
      [cg2F]       "r"  (cgrid2F),
      [fg2F]       "r"  (fgrid2F)

      ,

      [fg2Fmax]     "r"  (fgrid2Fmax),
      [fg2FmaxIdx]  "r"  (fgrid2FmaxIdx),


      [Vk]         "m"  (VIIII[0]),
      [V1]         "m"  (V1111[0])


      : /* clobbered */
      "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"

    ) ;


    fgrid2F+=16;
    cgrid2F+=16;
    fgrid2Fmax+=16;
    fgrid2FmaxIdx+=16;

#endif // EXP_No_ASM


  }

  /* take care of remaining iterations, length  modulo 16 */
  for( ; ifreq_fg < length; ifreq_fg++ ) {

            fgrid2F[0] += cgrid2F[0] ;

            newMax=(cgrid2F[0] >= fgrid2Fmax[0]);
            fgrid2Fmax[0]=fmaxf(fgrid2Fmax[0],cgrid2F[0]);
            fgrid2FmaxIdx[0]=fgrid2FmaxIdx[0]*(1-newMax)+k*newMax;
            fgrid2F++;
            cgrid2F++;
            fgrid2Fmax++;
            fgrid2FmaxIdx++;

  } /* for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) { */


}




void gc_hotloop(REAL4 * fgrid2F, REAL4 * cgrid2F, UCHAR * fgridnc, REAL4 TwoFthreshold, UINT4 length  )  {
  UINT4 ifreq_fg;

  REAL4 VTTTT[4] __attribute__ ((aligned (16))) = { TwoFthreshold,TwoFthreshold,TwoFthreshold,TwoFthreshold };


for(ifreq_fg=0 ; ifreq_fg +16 < length; ifreq_fg+=16 ) {
    /* unrolled loop (16 iterations of original loop) */
#ifdef EXP_NO_ASM

#pragma ivdep
  for(int j=0 ; j < 16; j++ ) {
	    fgrid2F[0] += cgrid2F[0] ;

	    fgridnc[0] += (TwoFthreshold < cgrid2F[0]);
	    fgridnc++;

	    fgrid2F++;
	    cgrid2F++;
  }

#else
    __asm __volatile (
         "MOVUPS  (%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  (%[fg2F]),%%xmm3 \n\t"

         "MOVAPS %[Vthresh2F],%%xmm7 \n\t"

         "MOVUPS  0x10(%[cg2F]),%%xmm4 \n\t"
         "MOVUPS  0x20(%[cg2F]),%%xmm5 \n\t"
         "MOVUPS  0x30(%[cg2F]),%%xmm6 \n\t"

         /* Loop iterations 1...4 */

         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */

         "MOVAPS  (%[fgnc]),%%xmm1 \n\t"  /* vector of 16 (!) number count values (unsigned bytes) */

         "MOVAPS  %%xmm3,(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */

         "MOVAPS %%xmm7,%%xmm3 \n\t"
         "CMPLEPS %%xmm2,%%xmm3   \n\t"   /* compare the four coarse grid 2F values to four */
                                          /* copies of threshold value in parallel          */
                                          /* result is a vector of 4 integer values:        */
                                          /*        -1 if TwoFthreshold < cgrid2F[i]        */
                                          /*         0 otherwise                            */
	 				  /* (saved in xmm3 for later processing)           */

         /* Loop iterations 5...8  (same as above) */

         "MOVAPS  0x10(%[fg2F]),%%xmm2 \n\t"
         "ADDPS   %%xmm4,%%xmm2 \n\t"
         "MOVAPS  %%xmm2,0x10(%[fg2F]) \n\t"

         "MOVAPS %%xmm7,%%xmm0 \n\t"
         "CMPLEPS %%xmm4,%%xmm0   \n\t"

         "PACKSSDW %%xmm0,%%xmm3 \n\t" /* combine two vectors of 4 double words (0/-1) */
				       /* to a vector of 8 words of 0/-1 in %%xmm3 */

         /* Loop iterations 9...12  (same as above) */

         "MOVAPS  0x20(%[fg2F]),%%xmm4 \n\t"
         "ADDPS   %%xmm5,%%xmm4 \n\t"
         "MOVAPS  %%xmm4,0x20(%[fg2F]) \n\t"

         "MOVAPS %%xmm7,%%xmm4 \n\t"
         "CMPLEPS %%xmm5,%%xmm4   \n\t"

         /* Loop iterations 13...16  (same as above) */

         "MOVAPS  0x30(%[fg2F]),%%xmm2 \n\t"
         "ADDPS   %%xmm6,%%xmm2 \n\t"
         "MOVAPS  %%xmm2,0x30(%[fg2F]) \n\t"

         "MOVAPS %%xmm7,%%xmm0 \n\t"
         "CMPLEPS %%xmm6,%%xmm0   \n\t"

         "PACKSSDW %%xmm0,%%xmm4 \n\t" /* 8 words of 0/-1 in %%xmm4 */

         "PACKSSWB %%xmm4,%%xmm3 \n\t" /* 16 unsigned bytes of 0/-1 in %%xmm3 */

         "PSUBB %%xmm3, %%xmm1   \n\t"   /* subtracting vector from number count vector */
         "MOVAPS  %%xmm1,(%[fgnc]) \n\t" /* to increment number count if threshold reached */

/*  ---------------------------------------------------*/
:
      /* output */

      :
      /* input */
      [cg2F]       "r"  (cgrid2F),
      [fg2F]       "r"  (fgrid2F)

      ,

      [fgnc]       "r"  (fgridnc),

      [Vthresh2F]   "m"  (VTTTT[0])

      : /* clobbered */
      "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"

    ) ;

#endif // EXP_No_ASM
    fgrid2F+=16;
    cgrid2F+=16;
    fgridnc+=16;

  }
  /* take care of remaining iterations, length  modulo 16 */
  for( ; ifreq_fg < length; ifreq_fg++ ) {
	    fgrid2F[0] += cgrid2F[0] ;
	    fgridnc[0] += (TwoFthreshold < cgrid2F[0]);
	    fgridnc++;
	    fgrid2F++;
	    cgrid2F++;
  } /* for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) { */

}

void gc_hotloop_no_nc(REAL4 * fgrid2F, REAL4 * cgrid2F, UINT4 length  )  {
  UINT4 ifreq_fg;



for( ifreq_fg=0 ; ifreq_fg +16 < length; ifreq_fg+=16 ) {
    /* unrolled loop (16 iterations of original loop) */
#ifdef EXP_NO_ASM

#pragma ivdep
  for(int j=0 ; j < 16; j++ ) {
	    fgrid2F[0] += cgrid2F[0] ;
	    fgrid2F++;
	    cgrid2F++;
  }

#else
    __asm __volatile (
         "MOVUPS  (%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  (%[fg2F]),%%xmm3 \n\t"
         "MOVUPS  0x10(%[cg2F]),%%xmm4 \n\t"
         "MOVUPS  0x20(%[cg2F]),%%xmm5 \n\t"
         "MOVUPS  0x30(%[cg2F]),%%xmm6 \n\t"

         /* Loop iterations 1...4 */

         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */
         "MOVAPS  %%xmm3,(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */

         /* Loop iterations 5...8  (same as above) */

         "MOVAPS  0x10(%[fg2F]),%%xmm2 \n\t"
         "ADDPS   %%xmm4,%%xmm2 \n\t"
         "MOVAPS  %%xmm2,0x10(%[fg2F]) \n\t"

         /* Loop iterations 9...12  (same as above) */

         "MOVAPS  0x20(%[fg2F]),%%xmm4 \n\t"
         "ADDPS   %%xmm5,%%xmm4 \n\t"
         "MOVAPS  %%xmm4,0x20(%[fg2F]) \n\t"

         /* Loop iterations 13...16  (same as above) */

         "MOVAPS  0x30(%[fg2F]),%%xmm2 \n\t"
         "ADDPS   %%xmm6,%%xmm2 \n\t"
         "MOVAPS  %%xmm2,0x30(%[fg2F]) \n\t"

/*  ---------------------------------------------------*/
:
      /* output */

      :
      /* input */
      [cg2F]       "r"  (cgrid2F),
      [fg2F]       "r"  (fgrid2F)

      : /* clobbered */
      "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"

    ) ;

#endif // EXP_No_ASM
    fgrid2F+=16;
    cgrid2F+=16;

  }
  /* take care of remaining iterations, length  modulo 16 */
  for( ; ifreq_fg < length; ifreq_fg++ ) {
	    fgrid2F[0] += cgrid2F[0] ;
	    fgrid2F++;
	    cgrid2F++;
  } /* for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) { */

}
