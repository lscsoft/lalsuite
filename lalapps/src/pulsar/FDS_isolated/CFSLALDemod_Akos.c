/* LALDemod variant with structs and signatures for modifications by Feket Akos
 * Authors see ComputeFStatistic.c
                                                         Bernd Machenschalk */
RCSID( "$Id$");

#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#ifdef __APPLE__ /* specials for Apples assembler */
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
            
static REAL8 sinVal[LUT_RES+1]; /* Lookup tables for fast sin/cos calculation */
static REAL8 sinVal2PI[LUT_RES+1];
static REAL8 sinVal2PIPI[LUT_RES+1];
static REAL8 cosVal[LUT_RES+1];
static REAL8 cosVal2PI[LUT_RES+1];
static REAL8 cosVal2PIPI[LUT_RES+1];
static REAL8 diVal[LUT_RES+1];

#define klim 32

/* __attribute__ ((always_inline)) */
void LALDemodSub(COMPLEX8* Xalpha, INT4 sftIndex,
                 REAL8 tempFreq0, REAL8 tempFreq1, REAL8 x, REAL8 yTemp,
                 REAL8* realXPo, REAL8* imagXPo, REAL8* realQo, REAL8* imagQo)
{
  REAL4 tsin, tcos;
  REAL4 realP, imagP;
  INT4 k;                                 /* loop counter */
  REAL8 realXP, imagXP, realQ, imagQ;     /* output parameters */


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

    /* calculate tsin, tcos, realQ and imagQ from sin/cos LUT */

    static REAL4 lutr = LUT_RES; /* LUT_RES in memory */
    static REAL4 one  = 1.0;
    static REAL4 half = .5;
    static REAL8 yRem;

    __asm __volatile
      (
       /* calculation of tsin and tcos */
         
       /* calculate index and put into EAX */
       /*                                    vvvvv-- these comments keep track of the FPU stack */
       "fldl   %[x]                   \n\t" /* x */
       "flds   %[lutr]                \n\t" /* LUT_RES x */
       "fmul   %%st(1),%%st(0)        \n\t" /* (x*LUT_RES) x */
       
#ifdef USE_DEFAULT_ROUNDING_MODE
       /* The default rounding mode is truncation (round-to-zero), which is exactly what we want here.
          It should be restored by the compiler after every operation that involves floating-point-to-integer
          conversion (rounding). Switching the rounding mode is slow, so the following code is the fastest
          possible, as it simply expects the default rounding mode being active. However relying on this is
          kind of dangerous, so we don't do this by default. */
       "fistpl %[sinv]                \n\t" /* x */
       "movl   %[sinv],%%ebx          \n\t" /* x */
#else
       "fadds  %[half]                \n\t" /* (x*LUT_RES+.5) x */
       /* Implementation of floor() that doesn't rely on a rounding mode being set
          The current way works for positive values only! */
       /* This code temporary stores integer values in memory locations of float variables,
          but they're overwritten later anyway */
       "fistl  %[sinv]                \n\t" /* (x*LUT_RES+.5) x */ /* saving the rounded value, the original in FPU */
       "fisubl %[sinv]                \n\t" /* (x*LUT_RES+.5) x */ /* value - round(value) will be negative if was rounding up */
       "fstps  %[cosv]                \n\t" /* x */                /* we will check the sign in integer registers */
       "movl   %[sinv],%%ebx          \n\t" /* x */
       "sub    %%eax,%%eax            \n\t" /* x */                /* EAX=0 */
       "orl    %[cosv],%%eax          \n\t" /* x */                /* it will set the S (sign) flag */
       "jns    sincos1                \n\t" /* x */                /* the result is ok, rounding = truncation */
       "dec    %%ebx                  \n"   /* x */                /* sinv=sinv-1.0 (it was rounded up) */
       "sincos1:                      \n\t" /* x */
#endif
       
       /* calculate d = x - diVal[idx] in st(0) */
       "fsubl  %[diVal](,%%ebx,8)     \n\t" /* (d = x - diVal[i]) */
       
       /* copy d on the stack to prepare for calculating d*d */
       "fld    %%st(0)                \n\t" /* d d */
       
       /* this calculates d*d on _top_ of the stack */ 
       "fmul   %%st(0),%%st(0)        \n\t" /* (d*d) d */
       
       /* three-term Taylor expansion for sin value, starting with the last term,
          leaving d and d*d on stack, idx kept in ebx */
       "fldl  %[sinVal2PIPI](,%%ebx,8)\n\t" /* sinVal2PIPI[i] (d*d) d */
       "fmul  %%st(1),%%st(0)         \n\t" /* (d*d*sinVal2PIPI[i]) (d*d) d */
       
       "fldl  %[cosVal2PI](,%%ebx,8)  \n\t" /* cosVal2PI[i] (d*d*sinVal2PIPI[i]) (d*d) d */
       "fmul  %%st(3),%%st(0)         \n\t" /* (d*cosVal2PI[i]) (d*d*sinVal2PIPI[i]) (d*d) d */
       "fsubp                         \n\t" /* (d*cosVal2PI[i]-d*d*sinVal2PIPI[i]) (d*d) d */
       
       "faddl %[sinVal](,%%ebx,8)     \n\t" /* (sinVal[i]+d*cosVal2PI[i]-d*d*sinVal2PIPI[i]) (d*d) d */
       "fstps %[sinv]                 \n\t" /* (d*d) d */
       
       /* calculation for cos value, this time popping the stack */
       /* this ordering mimics the calculation of cos as done by gcc (4.1.0)  */
       /* WARNING: changing the order of the substractions gives an error up to 1e-3 in the result!! */
       
       "fmull %[cosVal2PIPI](,%%ebx,8)\n\t" /* (d*d*cosVal2PIPI[i]) d */
       "fxch                          \n\t" /* d (d*d*cosVal2PIPI[i]) */
       "fmull %[sinVal2PI](,%%ebx,8)  \n\t" /* (d*sinVal2PI[i]) (d*d*cosVal2PIPI[i]) d */
       "fsubrl %[cosVal](,%%ebx,8)    \n\t" /* (cosVal[i]-d*sinVal2PI[i]) (d*d*cosVal2PIPI[i]) */
       "fsubp                         \n\t" /* (cosVal[i]-d*sinVal2PI[i]-d*d*cosVal2PIPI[i]) */
       
#ifndef SKIP_COS_ROUNDING
       /* this stores the cos value in the single precision output variable before substracting 1.0,
          again to mimic the compilers code to get the same result. It is faster and might be better
          to simply leave this value on the stack, so we add a switch for skipping this */
       "fstps %[cosv]                 \n\t" /* % */
       "flds  %[cosv]                 \n\t" /* % */
#endif
       
       /* special here: tcos -= 1.0 */
       "fsub  %[one]                  \n\t" /* (cosVal[i]-d*(sinVal2PI[i]-d*cosVal2PIPI[i])-1) */
       "fstps %[cosv]                 \n\t" /* % */
       

       /* calculating yRem */
       
       "fldl   %[yTemp]               \n\t" /* yT */
#ifdef USE_DEFAULT_ROUNDING_MODE
       /* see comment about USE_DEFAULT_ROUNDING_MODE above */
       "fistl  %[yRem]                \n\t" /* yT */
       "fisubl %[yRem]                \n\t" /* (yT-(int)yT) */
       "fsts   %[yRem]                \n\t" /* % */
       "sub    %%eax,%%eax            \n\t" /*   EAX=0 */
       "orl    %[yRem],%%eax          \n\t" /*   sets the S (sign) flag */
       "jns    sincos3                \n\t" /*   jump if not negative */
       "fadds  %[one]                 \n\t" /* (yT-(int)yT+1) */
       "sincos3:                      \n\t"
       "fstpl  %[yRem]                \n\t" /* % */
#else
       "fsts   %[yRem]                \n\t" /* yT */
       "sub    %%eax,%%eax            \n\t" /*   EAX=0 */
       "orl    %[yRem],%%eax          \n\t" /*   sets the S (sign) flag */
       "js     sincos3                \n\t" /*   jump if negative */
       
       /* yTemp >= 0 */
       "fistl  %[yRem]                \n\t" /* yT */         /* saving the rounded value, the original in FPU */
       "fisubl %[yRem]                \n\t" /* yT-(int)yT */ /* value - round(value) will be negative if was rounding up */
       "fsts   %[yRem]                \n\t" /* yT-(int)yT */ /* we will check the sign in integer registers */
       "sub    %%eax,%%eax            \n\t" /* yT-(int)yT */ /* EAX=0 */
       "orl    %[yRem],%%eax          \n\t" /* yT-(int)yT */ /* sets the S (sign) flag */
       "jns    sincos4                \n\t" /* yT-(int)yT */ /* the result is ok, rounding = truncation */
       "fadds  %[one]                 \n\t" /* (yT-(int)yT+1) */
       "jmp    sincos4                \n"   /* (yT-(int)yT+1) */
       
       /* yTemp < 0 */
       /* it requires some mind-bending to come from the algorithm used in the USE_FLOOR case
          of the yRem calculation in the C code below to the case distinction implemented here.
          Hints:  yTemp<0 => yRem = yTemp - ceil(yTemp) + 1.0;
                  y<0 => ceil(y) = ((y-(int)y)<-.5) ? ((int)y-1.0) : ((int)y)
                  -1.0 + 1.0 = 0; 0 + 1.0 = +1.0
       */
       "sincos3:                      \n\t"
       "fistl  %[yRem]                \n\t" /* yT */
       "fisubl %[yRem]                \n\t" /* yT-(int)yT */ /* value - round(value) will be positive if was rounding up */
       "fsts   %[yRem]                \n\t" /* yT-(int)yT */ /* we will check the sign in integer registers */
       "sub    %%eax,%%eax            \n\t" /* yT-(int)yT */ /* EAX=0 */
       "orl    %[yRem],%%eax          \n\t" /* yT-(int)yT */ /* sets the S (sign) flag */
       "jns    sincos4                \n\t" /* yT-(int)yT */ /* */
       "fadds  %[one]                 \n"   /* (yT-(int)yT+1) */
       
       "sincos4:                      \n\t" /* yR */
       "fstpl  %[yRem]                \n"   /* % */
#endif
       
       /* calculation of realQ and imagQ */
       /* (see also comments in calculation of tsin and tcos */
       
       "fldl   %[yRem]                \n\t" /* x */
       "flds   %[lutr]                \n\t" /* LUT_RES x */
       "fmul   %%st(1),%%st(0)        \n\t" /* (x*LUT_RES) x */
       
#ifdef USE_DEFAULT_ROUNDING_MODE
       "fistpl %[sinq]                \n\t" /* x */
       "movl   %[sinq],%%ebx          \n\t" /* x */
#else
       "fadds  %[half]                \n\t" /* (x*LUT_RES+.5) x */
       "fistl  %[sinq]                \n\t" /* (x*LUT_RES+.5) x */ /* saving the rounded value, the original in FPU */
       "fisubl %[sinq]                \n\t" /* (x*LUT_RES+.5) x */ /* value - round(value) will be negative if was rounding up */
       "fstps  %[cosq]                \n\t" /* x */                /* we will check the sign in integer registers */
       "movl   %[sinq],%%ebx          \n\t" /* x */
       "sub    %%eax,%%eax            \n\t" /* x */                /* EAX=0 */
       "orl    %[cosq],%%eax          \n\t" /* x */                /* it will set the S (sign) flag */
       "jns    sincos2                \n\t" /* x */                /* the result is ok, rounding = truncation */
       "dec    %%ebx                  \n"   /* x */                /* sinv=sinv-1.0 (it was rounded up) */
       "sincos2:                      \n\t" /* x */
#endif
       "fsubl  %[diVal](,%%ebx,8)     \n\t" /* (d = x - diVal[i]) */
       "fld    %%st(0)                \n\t" /* d d */
       "fmul   %%st(0),%%st(0)        \n\t" /* (d*d) d */
       "fldl  %[sinVal2PIPI](,%%ebx,8)\n\t" /* sinVal2PIPI[i] (d*d) d */
       "fmul  %%st(1),%%st(0)         \n\t" /* (d*d*sinVal2PIPI[i]) (d*d) d */
       
       "fldl  %[cosVal2PI](,%%ebx,8)  \n\t" /* cosVal2PI[i] (d*d*sinVal2PIPI[i]) (d*d) d */
       "fmul  %%st(3),%%st(0)         \n\t" /* (d*cosVal2PI[i]) (d*d*sinVal2PIPI[i]) (d*d) d */
       "fsubp                         \n\t" /* (d*cosVal2PI[i]-d*d*sinVal2PIPI[i]) (d*d) d */
       "faddl %[sinVal](,%%ebx,8)     \n\t" /* (sinVal[i]+d*cosVal2PI[i]-d*d*sinVal2PIPI[i]) (d*d) d */

       /* the following sign change of the sin() result is special here,
	  otherwise the calculation is the same as the previous one */
       "fchs                          \n\t" /* -(sinVal[i]+d*cosVal2PI[i]-d*d*sinVal2PIPI[i]) (d*d) d */
       "fstpl %[sinq]                 \n\t" /* (d*d) d */
       
       "fmull %[cosVal2PIPI](,%%ebx,8)\n\t" /* (d*d*cosVal2PIPI[i]) d */
       "fxch                          \n\t" /* d (d*d*cosVal2PIPI[i]) */
       "fmull %[sinVal2PI](,%%ebx,8)  \n\t" /* (d*sinVal2PI[i]) (d*d*cosVal2PIPI[i]) d */
       "fsubrl %[cosVal](,%%ebx,8)    \n\t" /* (cosVal[i]-d*sinVal2PI[i]) (d*d*cosVal2PIPI[i]) */
       "fsubp                         \n\t" /* (cosVal[i]-d*sinVal2PI[i]-d*d*cosVal2PIPI[i]) */
       "fstpl %[cosq]                 \n\t" /* % */


       /* Kernel loop */

       /* the "inner" four complex values that contribute the most to the sum
	  will be calculated by the FPU with 80 Bit precision,
	  the SSE taking care of the others ("wings") */

       /* calculation (for 28 REAL4 values (7x(2 ReIm pairs))) */
       /* one SSE register will consist 4 REAL4 values */
       /* 4 REAL4 vaules = 2 ReIm pairs */

       /* we want to do this SSE code 14 times */
#define ADD4SSE(a,b) \
       "movlps "#a"(%[Xalpha_kX]),%%xmm4 \n\t"\
       "movhps "#b"(%[Xalpha_kX]),%%xmm4 \n\t"\
       "rcpps %%xmm0,%%xmm1       \n\t" /* XMM1: 1/(f-1) 1/(f-1) 1/f 1/f */\
       "subps %%xmm5,%%xmm0       \n\t" /* XMM2: f-3 f-3 f-2 f-2 */\
       "mulps %%xmm4,%%xmm1       \n\t" /* XMM1: ImH/(f-1) ReH/(f-1) ImL/f ReL/f */\
       "addps %%xmm1,%%xmm2       \n\t" /* XMM2: C_ImH C_ReH C_ImL C_ReL */
    

       /* constants */
       "jmp cntcode              \n"
       AD_ALIGN16 "              \n\t"
       "C_V0011:                 \n\t"
       AD_FLOAT "0.0             \n\t"
       AD_FLOAT "0.0             \n\t"
       "C_F1:                    \n\t"
       AD_FLOAT "1.0             \n\t"
       AD_FLOAT "1.0             \n"
       "C_F2:                    \n\t"
       "C_V2222:                 \n\t"
       AD_FLOAT "2.0             \n\t"
       AD_FLOAT "2.0             \n\t"
       AD_FLOAT "2.0             \n\t"
       AD_FLOAT "2.0             \n"
       
       AD_ASCII "\"$Id$\"\n"
       AD_ALIGN16 "              \n\t"
       "cntcode:                 \n\t"
       
       /* SSE prelude */
       "movss %[tempFreqX],%%xmm0\n\t"
       "shufps $0,%%xmm0,%%xmm0  \n\t" /* XMM0: f   f   f   f */
       "subps  C_V0011,%%xmm0    \n\t" /* XMM0: f-1 f-1 f   f */
       "xorps  %%xmm2,%%xmm2     \n\t" /* XMM2 will collect the low-precision values */
       "movups C_V2222,%%xmm5    \n\t"
       /* add two complex elements at a time */
       ADD4SSE(0,8)
       ADD4SSE(16,24)
       ADD4SSE(32,40)
       
       "FLDL %[tempFreq1]        \n\t" /* FLD D [TempFreqMinus15]   ;A15 */
       "FLDS -16(%[Xalpha_k])    \n\t" /* FLD D [EBP-10h]   ;X14 A15 */
       "FMUL %%ST(1),%%ST        \n\t" /* FMUL ST,ST1       ;X14A15 A15 */
       
       ADD4SSE(48,56)
       
       "FLDS -12(%[Xalpha_k])    \n\t" /* FLD D [EBP-0Ch]   ;Y14 X14A15 A15 */
       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;Y14A15 X14A15 A15 */
       "FXCH %%ST(2)             \n\t" /* FXCH ST2          ;A15 X14A15 Y14A15 */
       "FLDS C_F1                \n\t" /* FLD D [_ONE_]     ;1 A15 X14A15 Y14A15 */
       "FADD %%ST(1),%%ST        \n\t" /* FADD ST,ST1       ;A14 A15 X14A15 Y14A15 */
       
       ADD4SSE(64,72)
       
       "FLDS -8(%[Xalpha_k])     \n\t" /* FLD D [EBP-08h]   ;X15 A14 A15 X14A15 Y14A15 */
       "FMUL %%ST(1),%%ST        \n\t" /* FMUL ST,ST1       ;X15A14 A14 A15 X14A15 Y14A15 */
       "FADDP %%ST,%%ST(3)       \n\t" /* FADDP ST3,ST      ;A14 A15 X' Y14A15 */
       "FMUL %%ST,%%ST(1)        \n\t" /* FMUL ST1,ST       ;A14 Q145 X' Y14A15 */
       "FLDS -4(%[Xalpha_k])     \n\t" /* FLD D [EBP-04h]   ;Y15 A14 Q145 X' Y14A15 */
       
       ADD4SSE(80,88)
       
       "FMUL %%ST(1),%%ST        \n\t" /* FMUL ST,ST1       ;Y15A14 A14 Q145 X' Y14A15 */
       "FADDP %%ST,%%ST(4)       \n\t" /* FADDP ST4,ST      ;A14 Q145 X' Y' */
       "FSUBS C_F2               \n\t" /* FSUB D [_TWO_]    ;A16 Q145 X' Y' */
       "FMUL %%ST,%%ST(2)        \n\t" /* FMUL ST2,ST       ;A16 Q145 X'A16 Y' */
       "FMUL %%ST,%%ST(3)        \n\t" /* FMUL ST3,ST       ;A16 Q145 X'A16 Y'A16 */
       
       ADD4SSE(96,104) 
       
       /* SSE: skip FPU calculated values */
       "subps %%xmm5,%%xmm0      \n\t"
       "addl  $144,%[Xalpha_kX]  \n\t" /* Xalpha_kX = Xalpha_kX + 4; */
       "subps %%xmm5,%%xmm0      \n\t"
       
       "FLDS 0(%[Xalpha_k])      \n\t" /* FLD D [EBP+00h]   ;X16 A16 Q145 X'A16 Y'A16 */
       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;X16Q145 A16 Q145 X'A16 Y'A16 */
       "FADDP %%ST,%%ST(3)       \n\t" /* FADDP ST3,ST      ;A16 Q145 X" Y'A16 */
       "FLDS 4(%[Xalpha_k])      \n\t" /* FLD D [EBP+04h]   ;Y16 A16 Q145 X" Y'A16 */
       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;Y16Q145 A16 Q145 X" Y'A16 */
       
       ADD4SSE(0,8)
       
       "FADDP %%ST,%%ST(4)       \n\t" /* FADDP ST4,ST      ;A16 Q145 X" Y" */
       "FMUL %%ST,%%ST(1)        \n\t" /* FMUL ST1,ST       ;A16 Q146 X" Y" */
       "FSUBS C_F1               \n\t" /* FSUB D [_ONE_]    ;A17 Q146 X" Y" */
       "FMUL %%ST,%%ST(2)        \n\t" /* FMUL ST2,ST       ;A17 Q146 X"A17 Y" */
       "FMUL %%ST,%%ST(3)        \n\t" /* FMUL ST3,ST       ;A17 Q146 X"A17 Y"A17 */
       
       ADD4SSE(16,24)
       
       "FLDS 8(%[Xalpha_k])      \n\t" /* FLD D [EBP+08h]   ;X17 A17 Q146 X"A17 Y"A17 */
       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;X17Q146 A17 Q146 X"A17 Y"A17 */
       "FADDP %%ST,%%ST(3)       \n\t" /* FADDP ST3,ST      ;A17 Q146 X! Y"A17 */
       
       ADD4SSE(32,40)
       
       "FLDS 12(%[Xalpha_k])     \n\t" /* FLD D [EBP+0Ch]   ;Y17 A17 Q146 X! Y"A17 */
       "FMUL %%ST(2),%%ST        \n\t" /* FMUL ST,ST2       ;Y17Q146 A17 Q146 X! Y"A17 */
       "FADDP %%ST,%%ST(4)       \n\t" /* FADDP ST4,ST      ;A17 Q146 X! Y! */
       "FMULP %%ST,%%ST(1)       \n\t" /* FMULP ST1,ST      ;Q147 X! Y! */
       
       ADD4SSE(48,56)
       
       "FDIVRS C_F1              \n\t" /* FDIVR D [_ONE_]   ;1/Q x y */
       "FMUL %%ST,%%ST(1)        \n\t" /* FMUL ST1,ST       ;1/Q xq y */
       "FMULP %%ST,%%ST(2)       \n\t" /* FMULP ST2,ST      ;xq yq */
       
       ADD4SSE(64,72)
       ADD4SSE(80,88)
       
       "FSTPS %[XRes]            \n\t"
       "FSTPS %[XIms]            \n\t"
       
       ADD4SSE(96,104) 
       
       /* add two complex elements at a time */
       
       /* add the two calculated parts of the real and imaginary part */
       "movhlps %%xmm2,%%xmm3    \n\t" /* XMM3: ? ? C_ImH C_ReH */
       "addps %%xmm3,%%xmm2      \n\t" /* XMM2: - - C_Im C_Re */
       
       /* store the result */
       "movlps %%xmm2, %[XSumsX] \n\t"
         

       /* interface */
       : /* output  (here: to memory)*/
       [sinq]        "=m" (imagQ),
       [cosq]        "=m" (realQ),
       [sinv]        "=m" (tsin),
       [cosv]        "=m" (tcos),
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
       [lutr]        "m"  (lutr),
       [one]         "m"  (one),
       [half]        "m"  (half),
       [sinVal]      "m"  (sinVal[0]),
       [cosVal]      "m"  (cosVal[0]),
       [sinVal2PI]   "m"  (sinVal2PI[0]),
       [cosVal2PI]   "m"  (cosVal2PI[0]),
       [sinVal2PIPI] "m"  (sinVal2PIPI[0]),
       [cosVal2PIPI] "m"  (cosVal2PIPI[0]),
       [diVal]       "m"  (diVal[0])

       : /* clobbered registers */
#ifndef IGNORE_XMM_REGISTERS
       "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5",
#endif
       "eax", "ebx",
       "st","st(1)","st(2)","st(3)","st(4)"
       );           
      
    /* And last, we add the single and double precision values */
    XRes = XRes + XSumsX.re;
    XIms = XIms + XSumsX.im;
    
    {
      REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
      REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
      
      realXP = tsin2pi * XRes - tcos2pi * XIms;
      imagXP = tcos2pi * XRes + tsin2pi * XIms;
    }
  }
  
  else /* if ( tempFreq0 >= LD_SMALL ) */
    
    {
      fprintf(stderr,"small x\n");

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
      for (k=0; k <= LUT_RES; k++) {
        sinVal[k] = sin((LAL_TWOPI*k)/(LUT_RES));
        sinVal2PI[k] = sinVal[k]  *  LAL_TWOPI;
        sinVal2PIPI[k] = sinVal2PI[k] * LAL_PI;
        cosVal[k] = cos((LAL_TWOPI*k)/(LUT_RES));
        cosVal2PI[k] = cosVal[k]  *  LAL_TWOPI;
        cosVal2PIPI[k] = cosVal2PI[k] * LAL_PI;
      }

      for (k=0; k <= LUT_RES; k++)
        diVal[k] = (REAL8)k/(REAL8)(LUT_RES);
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
