          {
	    typedef struct R44Vtag {
	      REAL4 a,b,c,d;
	    } R44V;

	    NRCSID (CFSLOOPx86MAS, "$Id$");

	    COMPLEX8 *Xalpha_k = Xalpha + sftIndex;

	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;

	    REAL4 XRes=0.0, XIms=0.0;   /* sums of Xa.re and Xa.im */
	    REAL4 XResX=0.0, XImsX=0.0;   /* sums of Xa.re and Xa.im */

	    REAL4 accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */

	    /* prepare values for SSE */
	    REAL4 tempFreqX = tempFreq1; /* REAL4 because of SSE */
	    COMPLEX8 *Xalpha_kX = Xalpha_k; /* -> SSE values */
	    R44V V0011,V2222;
	    
	    V0011.a = 0.0;
	    V0011.b = 0.0;
	    V0011.c = 1.0;
	    V0011.d = 1.0;

	    V2222.a = 2.0;
	    V2222.b = 2.0;
	    V2222.c = 2.0;
	    V2222.d = 2.0;

	    tempFreq1 = tempFreq1 - 14;
	    Xalpha_k = Xalpha_k + 14; /* -> FPU values */
 
	    /* let the compiler code the x87 part */
	    for(k=0; k < 4 ; k++)
	      {
		XRes = tempFreq1 * XRes + (*Xalpha_k).re * accFreq;
		XIms = tempFreq1 * XIms + (*Xalpha_k).im * accFreq;
		
		accFreq *= (REAL4)tempFreq1;
		tempFreq1 --;
		Xalpha_k ++;
	      } /* for k < klim */
	    
	    accFreq = 1.0 / accFreq;
	    XRes *= accFreq;
	    XIms *= accFreq;

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
	    
	    XRes = XRes + XResX;
	    XIms = XIms + XImsX;

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
