/********************************** <lalVerbatim file="LALBarycenterCV">
Author: Cutler, C.
$Id$
*********************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALBarycenter.c}}
\label{ss:LALBarycenter.c}

Converts from detector arrival time (recorded by GPS clock) to
pulse emission time, in TDB.

\subsubsection*{Prototypes}
\input{LALBarycenterCP}
\idx{LALBarycenter()}
\idx{LALBarycenter()}

\subsubsection*{Description}

\verb@LALBarycenterEarth()@ computes the position and orientation
of the Earth, at some arrival time $t_a$ , specified
\verb@LIGOTimeGPS@ input structure. The Einstein delay is also 
calculated. The results are stored in the \verb@EarthState@ output 
structure, which can then be fed as input to \verb@LALBarycenter()@.
The idea is that \verb@LALBarycenterEarth()@ calculates quantities that
are independent of the source location and detector position
on Earth. Thus this function is called ONCE for every desired
arrival time; the results are then re-used as one steps around the
sky (and/or changes detector) at that time.

\verb@LALBarycenter()@ transforms from detector arrival time $t_a$
in GPS (as specified in the \verb@LIGOTimeGPS@ structure) to pulse 
emission time $t_e$, in TDB. (Actually, the returned $t_e$ is
the emission time plus the constant light-travel-time from
source to SSB.) The inputs to \verb@LALBarycenter()@, through
the \verb@BarycenterInput@ structure, are the source location,
detector site identifier, and GPS arrival time.
The function returns the  emission time $t_e(t_a)$, the 
derivative $d t_e/d t_a$, and the difference 
$t_e(t_a) - t_a $ through the \verb@EmissionTime@ structure.
The emission time $t_e(t_a)$ is returned in the \verb@LIGOTimeGPS@ format,
while the other two quantities are \verb@REAL8@s.

\subsubsection*{Algorithm}
The function ``corrects'' the pulse arrival time
by removing the Roemer delay, 
(including effects of Moon, planets, and the 
Earth's time-varying spin axis and spin rate),
Einstein delay, and Shapiro delay. Accuracy is
better than 5 $\mu$s. Full details will be in monograph
by Cutler in Los Alamos preprint archive.


\subsubsection*{Uses}

\subsubsection*{Notes}
\vfill{\footnotesize\input{LALBarycenterCV}}

</lalLaTeX> */

#include <lal/LALBarycenter.h>
		 
NRCSID(LALBARYCENTERC,"$Id$");

/* <lalVerbatim file="LALBarycenterCP"> */
void 
LALBarycenterEarth(LALStatus *stat, EarthState *earth, LIGOTimeGPS *tGPS, 
           EphemerisData *edat) 
{ /* </lalVerbatim> */

  REAL8 tgps[2];   /*I convert from two-integer representation to 
                      two REAL8s (just because I initially wrote my code for
                      REAL8s). GPS time(sec) = tgps[0] + 1.e-9*tgps[1] */ 

/*int ihour;                 most recent hour to current time, tGPS,
                               in look-up table */
/*REAL8 t0,t0hr; */

  REAL8 dtEtable;  /*time between consecutive entries in Earth ephem. table */
  REAL8 t0e;        /*time since first entry in Earth ephem. table */  
  INT4 ientryE;      /*entry in look-up table closest to current time, tGPS */

  REAL8 tinitE;           /*time (GPS) of first entry in Earth ephem table*/
  REAL8 tdiffE;           /*current time tGPS minus time of nearest entry
                             in Earth ephem look-up table */
  REAL8 tdiff2E;          /*tdiff2=tdiffE*tdiffE */
  /*next entries are same as above, but for Sun table */
  REAL8 dtStable;
  REAL8 t0s;      
  INT4 ientryS;   
  REAL8 tinitS;  
  REAL8 tdiffS;  
  REAL8 tdiff2S;

  INT4 j; /*dummy index */

  INITSTATUS(stat,"LALBarycenter",LALBARYCENTERC); 

    tgps[0] = (REAL8)tGPS->gpsSeconds; /*convert from INT4 to REAL8 */
    tgps[1] = (REAL8)tGPS->gpsNanoSeconds;

    tinitE = edat->gpsE[0];
    tinitS = edat->gpsS[0];
    
    dtEtable = edat->gpsE[1] - tinitE;
    dtStable = edat->gpsS[1] - tinitS;
 
    t0e=tgps[0]-tinitE;
    ientryE = floor((t0e/dtEtable) +0.5e0);  /*finding Earth table entry 
                                                closest to arrival time*/
    /*Curt: make sure dtEtable and dtStable aren't zero */
    t0s=tgps[0]-tinitS;
    ientryS = floor((t0s/dtStable) +0.5e0);  /*finding Sun table entry 
                                                closest to arrival time*/

    tdiffE = t0e -dtEtable*ientryE + tgps[1]*1.e-9; /*tdiff is arrival 
            time minus closest Earth table entry; tdiff can be pos. or neg. */
    tdiff2E = tdiffE*tdiffE;

    tdiffS=t0s -dtStable*ientryS + tgps[1]*1.e-9; /*same for Sun*/
    tdiff2S=tdiffS*tdiffS;


    /********************************************************************
     *Calucate position and vel. of center of Earth.
     *We extrapolate from a table produced using JPL DE405 ephemeris.
     *---------------------------------------------------------------------
     */
  {  

    REAL8* pos=edat->position[ientryE]; /*Cartesian coords of center of Earth
                                      from DE405 ephem, in sec. 0=x,1=y,2=z */
    REAL8* vel=edat->velocity[ientryE];
    REAL8* acc=edat->acceleration[ientryE]; 

        for (j=0;j<3;j++){
       earth->posNow[j]=pos[j] + vel[j]*tdiffE + 0.5*acc[j]*tdiff2E;
       earth->velNow[j]=vel[j] + acc[j]*tdiffE;
        }
  }
       
       /********************************************************************
	* Now calculating Earth's rotational state.
	*---------------------------------------------------------------------
	*/
    {
   INT2 leapsSince2000; /*number of leap secs added to UTC calendar since
                            Jan 1, 2000 00:00:00 UTC; leap is a NEGATIVE number
                            for dates before Jan 1, 2000. */
   
   REAL8 tu0JC;   /*time elapsed on UT1 clock (in Julian centuries)
                     between Jan 1.5 2000 (= J2000 epoch)
                     and midnight before gps[0,j] */
   
   REAL8 tuJC;    /*time elapsed on UT1 clock (in Julian centuries)
                     between Jan 1.5 2000 (= J2000 epoch)
                     and the present time (tGPS) */
   
   INT4 tuInt;     /*number of full seconds (integer) on elapsed 
                     on UT1 clock from Jan 1, 2000 00:00:00 UTC to present*/
   REAL8 dtu; /*see comments below*/
   
   REAL8 gmst;  /*Greenwich mean sidereal time at present (tGPS), in sec */
   REAL8 gmst0; /*Greenwich mean sidereal time at prev. midnight, in sec */
   REAL8 gmstRad; /*gmst, but in radians*/

   INT4  ut1secSince1Jan2000, fullUt1days;
   REAL8 daysSinceJ2000;
   REAL8 eps0= 0.40909280422232891e0;/*obliquity of ecliptic at JD 245145.0*/  
 
   leapsSince2000 = 0; /*right # for Jan. 1, 2000 */
   leapsSince2000 = edat->leap - 13;   
   tuInt=tGPS->gpsSeconds -630720013; /*first subtract off value 
                                 of gps clock at Jan 1, 2000 00:00:00 UTC */
   ut1secSince1Jan2000 = tuInt-leapsSince2000;  /*next subtract 
      off # leap secs added to UTC calendar since Jan 1, 2000 00:00:00 UTC */     
   tuJC = (ut1secSince1Jan2000 + tgps[1]*1.e-9 - 43200)/(8.64e4*36525);
           /*UT1 time elapsed, in Julian centuries, since  Jan 1.5 2000 
	     (= J2000 epoch) and pulse arrival time*/

   fullUt1days=floor(ut1secSince1Jan2000/8.64e4);/*full days on ut1 clock 
                                       since Jan 1, 2000 00:00:00 UTC */


   tu0JC = (fullUt1days-0.5e0)/36525.0;  /* UT1 time elapsed, 
              in Julian centuries between Jan 1.5 2000 (= J2000 epoch) 
              and midnight before gps[0,j] */

   dtu= tuJC - tu0JC; /*time, in Jul. centuries, between pulse arrival time                            and previous midnight (UT1) */          
   
   daysSinceJ2000 = (tuInt -43200)/8.64e4; /*days(not int) on GPS 
					    (or TAI,etc.) clock since J2000 */

   /*---------------------------------------------------------------------- 
    *in below formula for gmst0 & gmst, on rhs we SHOULD use time elapsed 
    * in UT1, but instead we insert elapsed in UTC. This will lead to 
    * errors in tdb of order 1 microsec. 
    *----------------------------------------------------------------------
    */

   gmst0=24110.54841e0 + tu0JC*(8640184.812866e0 + tu0JC*(0.093104e0
							    -tu0JC*6.2e-6));      /*formula for Greenwich mean sidereal time
      on p.50 of Explanatory Supp. to Ast Alm. ...
      gmst unit is sec, where 2pi rad = 24*60*60sec */
   
   /*-----------------------------------------------------------------------
    *Explan of tu0JC*8640184.812866e0 term: when tu0JC0=0.01 Julian centuries
    *(= 1 year), Earth has spun an extra 86401.84 sec (approx one revolution)
    *with respect to the stars.  
    *
    *Note: this way of organizing calc of gmst was stolen from TEMPO;
    *see subroutine lmst (= local mean sidereal time) .
    *------------------------------------------------------------------------
    */
   
   
   gmst=gmst0+ dtu*(8.64e4*36525 + 8640184.812866e0 
		    +0.093104e0*(tuJC + tu0JC) 
		    -6.2e-6*(tuJC*tuJC + tuJC*tu0JC + tu0JC* tu0JC));
   
   gmstRad=gmst*LAL_PI/43200;
   
   /*------------------------------------------------------------------------
    *Now adding calculating M and N, describing  lunisolar precession;  
    *see Eq. 3.213-1 in Exlpanatory Supp 
    *-------------------------------------------------------------------------
    */
   
   
   earth->M=(1.2812323e0*LAL_PI/1.8e2)*daysSinceJ2000/36525;
   earth->N=(0.5567530e0*LAL_PI/1.8e2)*daysSinceJ2000/36525;
   
   /*--------------------------------------------------------------------------
    * Now adding approx nutation (= short-period,forced motion, by definition).
    * These two dominant terms, with periods 18.6 yrs (big term) and 
    * 0.500 yrs (small term),respp., give nutation to around 1 arc sec; see
    * p. 120 of Explan. Supp. The forced nutation amplitude 
    *  is around 17 arcsec.
    *
    * Note the unforced motion or Chandler wobble (called ``polar motion'' 
    * in Explanatory Supp) is not included here. However its amplitude is
    * order of (and a somewhat less than) 1 arcsec; see plot on p. 270 of
    * Explanatory Supplement to Ast. Alm. 
    *
    * Note obliquity epsilon also changes with nutation, as 
    * given in 3.222-1, p.114 of Explan. Supp. 
    * So far, this has been ignored.
    *--------------------------------------------------------------------------
    */
   /* define variables for storing sin/cos of deltaE, eps0 etc */

   earth->delpsi = (-0.0048e0*LAL_PI/180.e0)*
     sin( (125.e0 - 0.05295e0*daysSinceJ2000)*LAL_PI/180.e0 )
     - (4.e-4*LAL_PI/180.e0)*
     sin( (200.9e0 + 1.97129e0*daysSinceJ2000)*LAL_PI/180.e0 );
   

   earth->deleps = (0.0026e0*LAL_PI/180.e0)*
     cos( (125.e0 - 0.05295e0*daysSinceJ2000)*LAL_PI/180.e0 )
     + (2.e-4*LAL_PI/180.e0)*
     cos( (200.9e0 + 1.97129e0*daysSinceJ2000)*LAL_PI/180.e0 );
   
     earth->gastRad = gmstRad + earth->delpsi*cos(eps0); 
                             /* ``equation of the equinoxes''*/
}
    /**********************************************************************
     *Now calculating Einstein delay. This is just difference between TDB and TDT.
     *We steal from TEMPO the approx 20 biggest terms in an expansion.
     *----------------------------------------------------------------------- 
     * jedtdt is Julian date (TDT) MINUS 2451545.0 (TDT);
     *  -7300.5 = 2444244.5 (Julian date when gps clock started) 
     *        - 2451545.0 TDT (approx. J2000.0); 
     * Note the 51.184 s difference between TDT and GPS 
     *-----------------------------------------------------------------------
     */
    {
    REAL8 jedtdt = -7300.5e0 + (tgps[0] + 51.184e0 + tgps[1]*1.e-9)/8.64e4; 
    REAL8 jt=jedtdt/3.6525e5; /*converting to TEMPO expansion param 
                               = Julian millenium, NOT Julian century */
    earth->einstein = 1.e-6*( 
       1656.674564e0 * sin(6283.075849991e0*jt + 6.240054195e0 )+
       22.417471e0 * sin(5753.384884897e0*jt + 4.296977442e0 )  +
       13.839792e0 * sin(12566.151699983e0*jt + 6.196904410e0 )  +
       4.770086e0 * sin(529.690965095e0*jt + 0.444401603e0 )   +
       4.676740e0 * sin(6069.776754553e0 *jt + 4.021195093e0 )   +
       2.256707e0 * sin(213.299095438e0 *jt + 5.543113262e0 )   +
       1.694205e0 * sin(-3.523118349e0 *jt + 5.025132748e0 )   +
       1.554905e0 * sin(77713.771467920e0 *jt + 5.198467090e0 )   +
       1.276839e0 * sin(7860.419392439e0 *jt + 5.988822341e0 )   +
       1.193379e0 * sin(5223.693919802e0 *jt + 3.649823730e0 )   +
       1.115322e0 * sin(3930.209696220e0 *jt + 1.422745069e0 )   +
       0.794185e0 * sin(11506.769769794e0 *jt + 2.322313077e0 )   +
       0.447061e0 * sin(26.298319800e0 *jt + 3.615796498e0 )   +
       0.435206e0 * sin(-398.149003408e0 *jt + 4.349338347e0 )   +
       0.600309e0 * sin(1577.343542448e0 *jt + 2.678271909e0 )   +
       0.496817e0 * sin(6208.294251424e0 *jt + 5.696701824e0 )   +
       0.486306e0 * sin(5884.926846583e0 *jt + 0.520007179e0 )   +
       0.432392e0 * sin(74.781598567e0 *jt + 2.435898309e0 )   +
       0.468597e0 * sin(6244.942814354e0 *jt + 5.866398759e0 )   +
       0.375510e0 * sin(5507.553238667e0 *jt + 4.103476804e0 )
       );

/* below, I've just taken derivative of above expression for einstein,
   then commented out terms that contribute less than around 10^{-12}
   to tDotBary */

    earth->deinstein = 1.e-6*( 
       1656.674564e0*6283.075849991e0*
                 cos(6283.075849991e0*jt + 6.240054195e0 )+

       22.417471e0*5753.384884897e0*
               cos(5753.384884897e0*jt + 4.296977442e0 )  +

       13.839792e0*12566.151699983e0*
               cos(12566.151699983e0*jt + 6.196904410e0 ) +

/* neglect next term
       4.770086e0*529.690965095e0*
               cos(529.690965095e0*jt + 0.444401603e0 )   +
*/
       4.676740e0*6069.776754553e0* 
              cos(6069.776754553e0*jt + 4.021195093e0 )   +

/* neglect next 2 terms
       2.256707e0*213.299095438e0* 
              cos(213.299095438e0*jt + 5.543113262e0 )    +

       1.694205e0*-3.523118349e0* 
              cos(-3.523118349e0*jt + 5.025132748e0 )     +
*/

       1.554905e0*77713.771467920e0*
              cos(77713.771467920e0*jt + 5.198467090e0 )  

/* neglect all the rest
       +   1.276839e0*7860.419392439e0* 
              cos(7860.419392439e0*jt + 5.988822341e0 )   +

       1.193379e0*5223.693919802e0*
              cos(5223.693919802e0*jt + 3.649823730e0 )   +

       1.115322e0*3930.209696220e0* 
              cos(3930.209696220e0*jt + 1.422745069e0 )   +

       0.794185e0*11506.769769794e0*
              cos(11506.769769794e0 *jt + 2.322313077e0 )  +

       0.447061e0*26.298319800e0*
              cos(26.298319800e0*jt + 3.615796498e0 )     +

       0.435206e0*-398.149003408e0* 
              cos(-398.149003408e0*jt + 4.349338347e0 )   +

       0.600309e0*1577.343542448e0*
              cos(1577.343542448e0*jt + 2.678271909e0 )   +

       0.496817e0*6208.294251424e0*
              cos(6208.294251424e0*jt + 5.696701824e0 )   +

       0.486306e0*5884.926846583e0*
              cos(5884.926846583e0*jt + 0.520007179e0 )   +

       0.432392e0*74.781598567e0*
              cos(74.781598567e0*jt + 2.435898309e0 )     +

       0.468597e0*6244.942814354e0* 
              cos(6244.942814354e0*jt + 5.866398759e0 )   +

       0.375510e0*5507.553238667e0*
              cos(5507.553238667e0*jt + 4.103476804e0 )
*/

       )/(8.64e4*3.6525e5);
    }
    /********************************************************************
     *Now calculating Earth-Sun separation vector, as needed
     *for Shapiro delay calculation. 
     *--------------------------------------------------------------------
     */
    {
      REAL8 rse2;
/* Curt: lately trying to re-write this with loop, to make more compact */
       REAL8* sunPos=edat->sunPos[ientryS];
       REAL8* sunVel=edat->sunVel[ientryS];
       REAL8* sunAcc=edat->sunAccel[ientryS];
       REAL8 sunPosNow[3], sunVelNow[3];
/*       REAL8 se[3],dse[3];  vector pointing from Sun to Earth (in sec),
                              and its time deriv (sec/sec)  */
       
       rse2=earth->drse=0.0;
       for (j=0;j<3;j++){
       sunPosNow[j]=sunPos[j] + sunVel[j]*tdiffS + 0.5*sunAcc[j]*tdiff2S;
       sunVelNow[j]=sunVel[j] + sunAcc[j]*tdiffS; 
      
       earth->se[j]=earth->posNow[j]-sunPosNow[j];
       earth->dse[j]=earth->velNow[j]-sunVelNow[j];
       rse2 += earth->se[j]*earth->se[j];
       earth->drse += earth->se[j]*earth->dse[j];
       }
       earth->rse=sqrt(rse2);
       earth->drse = earth->drse/earth->rse;
       /* Curt: make sure rse, b, (rse +seDotN) aren't zero */
     }
 RETURN(stat);
  }


/* <lalVerbatim file="LALBarycenterCP"> */
void 
LALBarycenter(LALStatus *stat, EmissionTime *emit, BarycenterInput *baryinput, 
           EarthState *earth) 
{ /* </lalVerbatim> */

  REAL8 longitude,latitude,rd;  /*geocentric (not geodetic!!) longitude
                                  and latitude of detector vertex, and
                                  dist. rd from center of Earth (sec). */
  REAL8 alpha,delta;  /* RA and DEC (radians) in 
                           ICRS realization of J2000 coords.*/

  REAL8 sinTheta;  /* sin(theta) = sin(pi/2 - delta) */   

  REAL8 tgps[2];   /*I convert from two-integer representation to 
                      two REAL8s (just because I initially wrote my code for
                      REAL8s). GPS time(sec) = tgps[0] + 1.e-9*tgps[1] */ 

  REAL8 roemer,droemer;  /*Roemer delay and its time derivative*/
  REAL8 erot,derot;      /*delay from center-of-Earth to detector (sec),
			   and its time deriv */

/*  REAL8 einstein;  /TDB-TDT; largest terms in TEMPO expansion
  REAL8 deinstein;   d(TDB-TDT)/d(TDT) */
 
  REAL8 shapiro, dshapiro; /*Shapiro delay due to Sun, and its time deriv. */
  REAL8 s[3]; /*unit vector pointing at source, in J2000 Cartesian coords */
  INT4 j; /*dummy index */

  INITSTATUS(stat,"LALBarycenter",LALBARYCENTERC); 

  ASSERT (baryinput!=NULL, stat, LALBARYCENTERH_ENULL, LALBARYCENTERH_MSGENULL);

  ASSERT (emit!= NULL, stat, LALBARYCENTERH_ENULL, LALBARYCENTERH_MSGENULL);

  tgps[0] = (REAL8)baryinput->tgps.gpsSeconds; /*convert from INT4 to REAL8 */
  tgps[1] = (REAL8)baryinput->tgps.gpsNanoSeconds; 
    
    alpha = baryinput->alpha;
    delta = baryinput->delta;
    sinTheta=sin(LAL_PI/2.0-delta);
    s[2]=cos(LAL_PI/2.0-delta);   /* s is vector that points towards source */
    s[1]=sinTheta*sin(alpha);  /* in Cartesian coords based on J2000 */
    s[0]=sinTheta*cos(alpha);  /* 0=x,1=y,2=z */
    
    rd = sqrt( baryinput->site.location[0]*baryinput->site.location[0]
              +baryinput->site.location[1]*baryinput->site.location[1]
              +baryinput->site.location[2]*baryinput->site.location[2]);
    latitude = LAL_PI/2.e0 - acos(baryinput->site.location[2]/rd);
    longitude= atan2(baryinput->site.location[1],baryinput->site.location[0]);
    rd /= LAL_C_SI;

    /********************************************************************
     *Calucate Roemer delay for detector at center of Earth.
     *We extrapolate from a table produced using JPL DE405 ephemeris.
     *---------------------------------------------------------------------
     */
  {  
       roemer=droemer=0.0;
       for (j=0;j<3;j++){
       roemer += s[j]*earth->posNow[j];
       droemer+= s[j]*earth->velNow[j]; 
       }
  }
       
       /********************************************************************
	* Now including Earth's rotation 
	*---------------------------------------------------------------------
	*/
  {
    REAL8 OMEGA = 7.29211510e-5;  /*ang. vel. of Earth (rad/sec)*/
    
    REAL8 eps0= 0.40909280422232891e0;/*obliquity of ecliptic at JD 245145.0, 
					in radians. NOT! to be confused with 
					LAL_EPSILON0_SI = permittivity of free
					space; value from Explan. Supp. 
					to Astronom. Almanac */
    /*  REAL8 eps0= (23.e0 + 26.e0/60 + 21.448e0/3.6e3)*LAL_PI/180;  */
    
   REAL8 NdotD;
   REAL8 alphaE,deltaE;
   
   alphaE = alpha + earth->M + earth->N*sin(alpha)*tan(delta) ;
   deltaE = delta + earth->N*cos((alphaE+ alpha)/2.e0);
   /*--------------------------------------------------------------------------
    * Now adding approx nutation (= short-period,forced motion, by definition).
    * These two dominant terms, with periods 18.6 yrs (big term) and 
    * 0.500 yrs (small term),respp., give nutation to around 1 arc sec; see
    * p. 120 of Explan. Supp. The forced nutation amplitude 
    *  is around 17 arcsec.
    *
    * Note the unforced motion or Chandler wobble (called ``polar motion'' 
    * in Explanatory Supp) is not included here. However its amplitude is
    * order of (and a somewhat less than) 1 arcsec; see plot on p. 270 of
    * Explanatory Supplement to Ast. Alm. 
    *
    * Note obliquity epsilon also changes with nutation, as 
    * given in 3.222-1, p.114 of Explan. Supp. 
    * So far, this has been ignored.
    *--------------------------------------------------------------------------
    */
 
   alphaE = alphaE + earth->delpsi*(cos(eps0) 
			       + sin(eps0)*sin(alphaE)*tan(deltaE))
     -earth->deleps*cos(alphaE)*tan(deltaE);
   
   deltaE = deltaE + earth->delpsi*sin(eps0)*cos(alphaE)
     + earth->deleps*sin(alphaE);
 
/*-----------------------------------------------------------------------
 * now taking dot product of detector location on earth and N 
 *-----------------------------------------------------------------------
 */
   NdotD = sin(latitude)*sin(deltaE)
                + cos(latitude)*cos(deltaE)*cos(earth->gastRad+longitude-alphaE);
   erot= rd*NdotD;
   derot= -OMEGA*rd*cos(latitude)*cos(deltaE)*sin(earth->gastRad+longitude-alphaE);
}

    /********************************************************************
     *Now adding Shapiro delay. Note according to J. Taylor review article
     *on pulsar timing, max value of Shapiro delay (when rays just graze sun)
     *is 120 microsec.
     * 
     *Here we calculate Shapiro delay 
     *for a detector at the center of the Earth.
     *Causes errors of order 10^{-4}sec * 4 * 10^{-5} = 4*10^{-9} sec
     *--------------------------------------------------------------------
     */
    {
      REAL8 seDotN,dseDotN,b,db;
       /* Curt: make sure rse, b, (rse +seDotN) aren't zero */
       seDotN = earth->se[2]*sin(delta)+ (earth->se[0]*cos(alpha) 
				     + earth->se[1]*sin(alpha))*cos(delta);
       dseDotN=earth->dse[2]*sin(delta)+(earth->dse[0]*cos(alpha)
				     +earth->dse[1]*sin(alpha))*cos(delta);

       b=sqrt(earth->rse*earth->rse-seDotN*seDotN);
       db=(earth->rse*earth->drse-seDotN*dseDotN)/b;

       shapiro  = 9.852e-6*log( (LAL_AU_SI/LAL_C_SI)/(earth->rse +seDotN));
       dshapiro = -9.852e-6*(earth->drse+ dseDotN)/(earth->rse +seDotN);
    }
    /*-----------------------------------------------------------------------
     *Now adding it all up. 
     *tBary (what we want to know) is pulse emission time in TDB coords.
     *The lab measures GPS time, but we set tGPS =0 at first instant of 1998
     *-----------------------------------------------------------------------
     */
    {  

  INT4 deltaTint; /* floor of deltaT */  

  emit->deltaT = roemer + erot + earth->einstein - shapiro;

  emit->tDot = 1.e0 + droemer + derot + earth->deinstein - dshapiro;
  
  deltaTint = (INT4)floor(emit->deltaT);

  if ( (1.e-9*tgps[1] + emit->deltaT - deltaTint) >= 1.e0 )
    {emit->te.gpsSeconds = baryinput->tgps.gpsSeconds + deltaTint + 1;
     emit->te.gpsNanoSeconds = (INT4)floor(1.e9*(tgps[1]*1.e-9 
                        + emit->deltaT - deltaTint -1.e0));
       }
   else
    {emit->te.gpsSeconds = baryinput->tgps.gpsSeconds + deltaTint;
     emit->te.gpsNanoSeconds = (INT4)floor(1.e9*(tgps[1]*1.e-9 
                        + emit->deltaT - deltaTint));
       }

  {  
       for (j=0;j<3;j++){
       emit->rDetector[j] = earth->posNow[j];
       emit->vDetector[j] = earth->velNow[j]; 
       }
  }
       /*Curt: rem to add Cartesian components of Earth rotation to above--
	 see your old notes*/
    }              
RETURN(stat);
  }










