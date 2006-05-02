          {
	    NRCSID (CFSLOOPX86GASS, "$Id$");

            /* COMPLEX8 *Xalpha_k = Xalpha + sftIndex; */
	    /* THIS IS DANGEROUS!! It relies on current implementation of COMPLEX8 type!! */
	    REAL4 *Xalpha_kR4 = &(Xalpha[sftIndex].re);

	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;

	    REAL4 accFreq = 1.0;        /* accumulating frequency factor, becomes common denominator */
	    REAL4 accFreq1;             /* accumulating frequency factor, becomes common denominator */
	    REAL4 XRes=0.0, XIms=0.0;   /* sums of Xa.re and Xa.im */
	    REAL4 XRes1=0.0, XIms1=0.0; /* sums of Xa.re and Xa.im */

	    REAL4 tFreq = tempFreq0;
	    INT4  tFint = params->Dterms - 1;

	    
	    /* For now this is just for testing inline assembly:
	       just do accFreq1=accFreq in assembly using SSE
	    */

	    __asm __volatile
	      (
	       "movss %1,%%xmm0 \n\t"
	       "movss %%xmm0,%0 \n\t"

	       ".intel_syntax   \n\t"
	       /*
		 one can write intel syntax in here, but the external references
		 (input and output variables) have to be dealt with outside
		 these directives in AT&T syntax :-(
	       */
	       ".att_syntax     \n\t"

	       : "=m" (accFreq1) /* output (here: memory)*/
	       : "m"  (accFreq)  /* input  (here: memory)*/
	       : "%xmm0"         /* clobbered */
	       );


            for(k=0; k < klim / 2; k++)
	      {
                XRes = tFreq * XRes + tFint * XRes + *Xalpha_kR4 * accFreq;
                Xalpha_kR4 ++;
                XIms = tFreq * XIms + tFint * XIms + *Xalpha_kR4 * accFreq;
                Xalpha_kR4 ++;

		accFreq *= (tFreq + tFint);
                tFint--;

                XRes1 = tFreq * XRes1 + tFint * XRes1 + *Xalpha_kR4 * accFreq1;
                Xalpha_kR4 ++;
                XIms1 = tFreq * XIms1 + tFint * XIms1 + *Xalpha_kR4 * accFreq1;
                Xalpha_kR4 ++;

		accFreq1 *= (tFreq + tFint);
                tFint--;
              }

	    XRes = (XRes * accFreq1 + XRes1 * accFreq) / (accFreq * accFreq1);
	    XIms = (XIms * accFreq1 + XIms1 * accFreq) / (accFreq * accFreq1);

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
