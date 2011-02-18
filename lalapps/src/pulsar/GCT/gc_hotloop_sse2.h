static inline void gc_hotloop (REAL4 * fgrid2F, REAL4 * cgrid2F, UCHAR * fgridnc, REAL4 TwoFthreshold, UINT4 length  ) __attribute__ ((hot));

#ifdef __APPLE__

/* Apple's gcc aligns ok */
#define ALRealloc LALRealloc
#define ALFree LALFree

#elif defined (__MINGW32__)

#define ALRealloc(p,s) __mingw_aligned_realloc(p,s,16)
#define ALFree __mingw_aligned_free

#else // neither APPLE nor MinGW

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


void gc_hotloop(REAL4 * fgrid2F, REAL4 * cgrid2F, UCHAR * fgridnc, REAL4 TwoFthreshold, UINT4 length  )  {
  UINT4 ifreq_fg;

  REAL4 VTTTT[4] __attribute__ ((aligned (16))) = { TwoFthreshold,TwoFthreshold,TwoFthreshold,TwoFthreshold };

  /* ensure alignment on fine grid. note fgridnc is UCHAR*
     while fgrid2F is UINT4*  */

  int offset = ((UINT4)fgridnc & 0xf) ;

  if(offset != 0) offset = 16-offset;

  for(ifreq_fg=0; offset > 0 && ifreq_fg < length; ifreq_fg++, offset-- ) {

	    fgrid2F[0] += cgrid2F[0] ;

#ifndef EXP_NO_NUM_COUNT	    
	    fgridnc[0] += (TwoFthreshold < cgrid2F[0]);
	    fgridnc++;
#endif
	    fgrid2F++;
	    cgrid2F++;

  } /* for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) { */


for( ; ifreq_fg +16 < length; ifreq_fg+=16 ) {
    /* unrolled loop (16 iterations of original loop) */
#ifdef EXP_NO_ASM

#pragma ivdep
  for(int j=0 ; j < 16; j++ ) {
	    fgrid2F[0] += cgrid2F[0] ;
#ifndef EXP_NO_NUM_COUNT	    
	    fgridnc[0] += (TwoFthreshold < cgrid2F[0]);
	    fgridnc++;
#endif
	    fgrid2F++;
	    cgrid2F++;
  }	    

#else
    __asm __volatile (
         "MOVUPS  (%[cg2F]),%%xmm2 \n\t"  /* load coarse grid values, possibly unaligned */
         "MOVAPS  (%[fg2F]),%%xmm3 \n\t"
#ifndef EXP_NO_NUM_COUNT	             
         "MOVAPS %[Vthresh2F],%%xmm7 \n\t"
#endif         
         "MOVUPS  0x10(%[cg2F]),%%xmm4 \n\t"
         "MOVUPS  0x20(%[cg2F]),%%xmm5 \n\t"
         "MOVUPS  0x30(%[cg2F]),%%xmm6 \n\t"

         /* Loop iterations 1...4 */

         "ADDPS   %%xmm2,%%xmm3 \n\t"     /* Add four coarse grid 2F values to fine grid sums */
#ifndef EXP_NO_NUM_COUNT	             
         "MOVAPS  (%[fgnc]),%%xmm1 \n\t"  /* vector of 16 (!) number count values (unsigned bytes) */
#endif         
         "MOVAPS  %%xmm3,(%[fg2F]) \n\t"  /* store 4 values in fine grid 2F sum array */

#ifndef EXP_NO_NUM_COUNT	             
         "MOVAPS %%xmm7,%%xmm3 \n\t"
         "CMPLEPS %%xmm2,%%xmm3   \n\t"   /* compare the four coarse grid 2F values to four */
                                          /* copies of threshold value in parallel          */
                                          /* result is a vector of 4 integer values:        */
                                          /*        -1 if TwoFthreshold < cgrid2F[i]        */
                                          /*         0 otherwise                            */
	 				  /* (saved in xmm3 for later processing)           */
#endif

         /* Loop iterations 5...8  (same as above) */

         "MOVAPS  0x10(%[fg2F]),%%xmm2 \n\t"
         "ADDPS   %%xmm4,%%xmm2 \n\t"
         "MOVAPS  %%xmm2,0x10(%[fg2F]) \n\t"
#ifndef EXP_NO_NUM_COUNT	                      
         "MOVAPS %%xmm7,%%xmm0 \n\t"
         "CMPLEPS %%xmm4,%%xmm0   \n\t"

         "PACKSSDW %%xmm0,%%xmm3 \n\t" /* combine two vectors of 4 double words (0/-1) */
				       /* to a vector of 8 words of 0/-1 in %%xmm3 */
#endif
         /* Loop iterations 9...12  (same as above) */

         "MOVAPS  0x20(%[fg2F]),%%xmm4 \n\t"
         "ADDPS   %%xmm5,%%xmm4 \n\t"
         "MOVAPS  %%xmm4,0x20(%[fg2F]) \n\t"

#ifndef EXP_NO_NUM_COUNT	             
         "MOVAPS %%xmm7,%%xmm4 \n\t"
         "CMPLEPS %%xmm5,%%xmm4   \n\t"
#endif

         /* Loop iterations 13...16  (same as above) */

         "MOVAPS  0x30(%[fg2F]),%%xmm2 \n\t"
         "ADDPS   %%xmm6,%%xmm2 \n\t"
         "MOVAPS  %%xmm2,0x30(%[fg2F]) \n\t"

#ifndef EXP_NO_NUM_COUNT	             
         "MOVAPS %%xmm7,%%xmm0 \n\t"
         "CMPLEPS %%xmm6,%%xmm0   \n\t"

         "PACKSSDW %%xmm0,%%xmm4 \n\t" /* 8 words of 0/-1 in %%xmm4 */

         "PACKSSWB %%xmm4,%%xmm3 \n\t" /* 16 unsigned bytes of 0/-1 in %%xmm3 */

         "PSUBB %%xmm3, %%xmm1   \n\t"   /* subtracting vector from number count vector */
         "MOVAPS  %%xmm1,(%[fgnc]) \n\t" /* to increment number count if threshold reached */
#endif

/*  ---------------------------------------------------*/
:
      /* output */

      :
      /* input */
      [cg2F]       "r"  (cgrid2F),
      [fg2F]       "r"  (fgrid2F)
      
#ifndef EXP_NO_NUM_COUNT	            
      ,

      [fgnc]       "r"  (fgridnc),

      [Vthresh2F]   "m"  (VTTTT[0])
#endif
      : /* clobbered */
      "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","memory"

    ) ;
    
#endif // EXP_No_ASM    
    fgrid2F+=16;
    cgrid2F+=16;
#ifndef EXP_NO_NUM_COUNT	        
    fgridnc+=16;
#endif

  }
  /* take care of remaining iterations, length  modulo 16 */
  for( ; ifreq_fg < length; ifreq_fg++ ) {
	    fgrid2F[0] += cgrid2F[0] ;
#ifndef EXP_NO_NUM_COUNT	    
	    fgridnc[0] += (TwoFthreshold < cgrid2F[0]);
	    fgridnc++;
#endif
	    fgrid2F++;
	    cgrid2F++;
  } /* for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) { */

}
