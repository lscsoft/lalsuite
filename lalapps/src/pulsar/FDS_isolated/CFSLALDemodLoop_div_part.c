          {
	    static const char* CFSLoopRCSID = "$Id$";

            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;
	    REAL4 accFreq = 1.0; /* accumulating frequency factor, becomes common denominator */
	    REAL4 tsin2pi = tsin * (REAL4)OOTWOPI;
	    REAL4 tcos2pi = tcos * (REAL4)OOTWOPI;
	    REAL4 XRes=0.0, XIms=0.0; /* sums of Xa.re and Xa.im */

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

            realXP = tsin2pi * XRes - tcos2pi * XIms;
            imagXP = tcos2pi * XRes + tsin2pi * XIms;
          } /* if x cannot be close to 0 */
