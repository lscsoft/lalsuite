/*
 * Copyright (C) 2012 Miroslav Shaltev, R Prix
*  Copyright (C) 2007 Curt Cutler, Jolien Creighton, Reinhard Prix, Teviet Creighton
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

#include <lal/Date.h>
#include <lal/LALBarycenter.h>

#define OBLQ 0.40909280422232891e0; /* obliquity of ecliptic at JD 245145.0* in radians */;

/// ---------- internal buffer type for optimized Barycentering function ----------
typedef struct tagfixed_sky
{
  REAL8 sinAlpha;	/// sin(alpha)
  REAL8 cosAlpha;	/// cos(alpha)
  REAL8 sinDelta;	/// sin(delta)
  REAL8 cosDelta;	/// cos(delta)
  REAL8 n[3];		/// unit vector pointing from SSB to the source, in J2000 Cartesian coords, 0=x,1=y,2=z
} fixed_sky_t;

typedef struct tagfixed_site
{
  REAL8 rd;		/// distance 'rd' from center of Earth, in light seconds
  REAL8 longitude;	/// geocentric (not geodetic!!) longitude of detector vertex
  REAL8 latitude;	/// geocentric latitude of detector vertex
  REAL8 sinLat;		/// sin(latitude)
  REAL8 cosLat;		/// cos(latitude);
  REAL8 rd_sinLat;	/// rd * sin(latitude)
  REAL8 rd_cosLat;	/// rd * cos(latitude)
} fixed_site_t;

struct tagBarycenterBuffer
{
  REAL8 alpha;			/// buffered sky-location: right-ascension in rad
  REAL8 delta;			/// buffered sky-location: declination in rad
  fixed_sky_t fixed_sky;	/// fixed-sky buffered quantities

  LALDetector site;		/// buffered detector site
  fixed_site_t fixed_site;	/// fixed-site buffered quantities

  BOOLEAN active;		/// switch set on TRUE of buffer has been filled
}; // struct tagBarycenterBuffer


/** \author Curt Cutler
 * \brief Computes the position and orientation of the Earth, at some arrival time
 * \f$t_a\f$, specified <tt>LIGOTimeGPS</tt> input structure.
 *
 * The Einstein delay is also calculated. The results are stored in the
 * <tt>EarthState</tt> output structure, which can then be fed as input to
 * <tt>LALBarycenter()</tt>.
 * The idea is that <tt>LALBarycenterEarth()</tt> calculates quantities that
 * are independent of the source location and detector position
 * on Earth. Thus this function is called ONCE for every desired
 * arrival time; the results are then re-used as one steps around the
 * sky (and/or changes detector) at that time.
 */
int
XLALBarycenterEarth ( EarthState *earth, 		/**< [out] the earth's state at time tGPS */
                      const LIGOTimeGPS *tGPS, 		/**< [in] GPS time tgps */
                      const EphemerisData *edat) 	/**< [in] ephemeris-files */
{
  REAL8 tgps[2];   /*I convert from two-integer representation to
                      two REAL8s (just because I initially wrote my code for
                      REAL8s). GPS time(sec) = tgps[0] + 1.e-9*tgps[1] */

/*int ihour;                 most recent hour to current time, tGPS,
                               in look-up table */
/*REAL8 t0,t0hr; */

  REAL8 t0e;        /*time since first entry in Earth ephem. table */
  INT4 ientryE;      /*entry in look-up table closest to current time, tGPS */

  REAL8 tinitE;           /*time (GPS) of first entry in Earth ephem table*/
  REAL8 tdiffE;           /*current time tGPS minus time of nearest entry
                             in Earth ephem look-up table */
  REAL8 tdiff2E;          /*tdiff2=tdiffE*tdiffE */
  /*next entries are same as above, but for Sun table */
  REAL8 t0s;
  INT4 ientryS;
  REAL8 tinitS;
  REAL8 tdiffS;
  REAL8 tdiff2S;

  INT4 j; /*dummy index */
  
  earth->ttype = TIMECORRECTION_ORIGINAL; /* using the original XLALBarycenterEarth function */
  
  /* check input */
  if ( !earth || !tGPS || !edat || !edat->ephemE || !edat->ephemS ) {
    XLALPrintError ("%s: invalid NULL input 'earth', 'tGPS', 'edat', 'edat->ephemE' or 'edat->ephemS'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

    tgps[0] = (REAL8)tGPS->gpsSeconds; /*convert from INT4 to REAL8 */
    tgps[1] = (REAL8)tGPS->gpsNanoSeconds;

    tinitE = edat->ephemE[0].gps;
    tinitS = edat->ephemS[0].gps;

    t0e=tgps[0]-tinitE;
    ientryE = floor((t0e/edat->dtEtable) +0.5e0);  /*finding Earth table entry
                                                closest to arrival time*/
    t0s=tgps[0]-tinitS;
    ientryS = floor((t0s/edat->dtStable) +0.5e0);  /*finding Sun table entry
                                                closest to arrival time*/

    /*Making sure tgps is within earth and sun ephemeris arrays*/

    if ( ( ientryE < 0 ) || ( ientryE >=  edat->nentriesE )) {
      XLALPrintError ("%s: input GPS time %f outside of Earth ephem range [%f, %f]\n", __func__, tgps[0], tinitE, tinitE + edat->nentriesE * edat->dtEtable );
      XLAL_ERROR ( XLAL_EDOM );
    }
    if ( ( ientryS < 0 ) || ( ientryS >=  edat->nentriesS ) ){
      XLALPrintError ("%s: input GPS time %f outside of Sun ephem range [%f, %f]\n", __func__, tgps[0], tinitS, tinitS + edat->nentriesS * edat->dtStable );
      XLAL_ERROR ( XLAL_EDOM );
    }

    tdiffE = t0e -edat->dtEtable*ientryE + tgps[1]*1.e-9; /*tdiff is arrival
	      time minus closest Earth table entry; tdiff can be pos. or neg. */
    tdiff2E = tdiffE*tdiffE;

    tdiffS=t0s -edat->dtStable*ientryS + tgps[1]*1.e-9; /*same for Sun*/
    tdiff2S=tdiffS*tdiffS;


    /********************************************************************
     *Calucate position and vel. of center of Earth.
     *We extrapolate from a table produced using JPL DE405 ephemeris.
     *---------------------------------------------------------------------
     */
    {

      REAL8* pos=edat->ephemE[ientryE].pos; /*Cartesian coords of center of Earth
					     from DE405 ephem, in sec. 0=x,1=y,2=z */
      REAL8* vel=edat->ephemE[ientryE].vel;
      REAL8* acc=edat->ephemE[ientryE].acc;

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
   REAL8 eps0= OBLQ;/*obliquity of ecliptic at JD 245145.0*/

   leapsSince2000 = 0; /*right # for Jan. 1, 2000 */
   leapsSince2000 = XLALGPSLeapSeconds( tGPS->gpsSeconds ) - 13;
   {
     INT4 err = xlalErrno;
     if ( err != XLAL_SUCCESS ) {
       XLALPrintError ("%s: XLALGPSLeapSeconds (%d) failed.\n", __func__, tGPS->gpsSeconds );
       XLAL_ERROR ( XLAL_EINVAL );
     }
   }

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

   earth->gmstRad = gmstRad;

   /*------------------------------------------------------------------------
    *Now calculating 3 terms, describing  lunisolar precession;
    *see Eqs. 3.212 on pp. 104-5 in Explanatory Supp.
    *-------------------------------------------------------------------------
    */

   earth->tzeA = tuJC*(2306.2181e0 + (0.30188e0 + 0.017998e0*tuJC)*tuJC )
             *LAL_PI/6.48e5;

   earth->zA = tuJC*(2306.2181e0 + (1.09468e0 + 0.018203e0*tuJC)*tuJC )
             *LAL_PI/6.48e5;

   earth->thetaA = tuJC*(2004.3109e0 - (0.42665e0 + 0.041833*tuJC)*tuJC )
             *LAL_PI/6.48e5;

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
     *Now calculating Einstein delay. This is just difference between
     *TDB and TDT.
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

/* now adding NEXT biggest (2nd-tier) terms from Tempo */
    earth->einstein = earth->einstein + 1.e-6*(
      0.243085 * sin(-775.522611324 *jt + 3.651837925 )   +
      0.173435 * sin(18849.227549974 *jt + 6.153743485 )   +
      0.230685 * sin(5856.477659115 *jt + 4.773852582 )   +
      0.203747 * sin(12036.460734888 *jt + 4.333987818 )   +
      0.143935 * sin(-796.298006816 *jt + 5.957517795 )   +
      0.159080 * sin(10977.078804699 *jt + 1.890075226 )   +
      0.119979 * sin(38.133035638 *jt + 4.551585768 )   +
      0.118971 * sin(5486.777843175 *jt + 1.914547226 )   +
      0.116120 * sin(1059.381930189 *jt + 0.873504123 )   +
      0.137927 * sin(11790.629088659 *jt + 1.135934669 )   +
      0.098358 * sin(2544.314419883 *jt + 0.092793886 )   +
      0.101868 * sin(-5573.142801634 *jt + 5.984503847 )   +
      0.080164 * sin(206.185548437 *jt + 2.095377709 )   +
      0.079645 * sin(4694.002954708 *jt + 2.949233637 )   +
      0.062617 * sin(20.775395492 *jt + 2.654394814 )   +
      0.075019 * sin(2942.463423292 *jt + 4.980931759 )   +
      0.064397 * sin(5746.271337896 *jt + 1.280308748 )   +
      0.063814 * sin(5760.498431898 *jt + 4.167901731 )   +
      0.048042 * sin(2146.165416475 *jt + 1.495846011 )   +
      0.048373 * sin(155.420399434 *jt + 2.251573730 )
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

/*Note I don't bother adding 2nd-tier terms to deinstein_tempo either */
    }
    
    /********************************************************************
     *Now calculating Earth-Sun separation vector, as needed
     *for Shapiro delay calculation.
     *--------------------------------------------------------------------
     */
  {
    REAL8 rse2;
    REAL8* sunPos=edat->ephemS[ientryS].pos;
    REAL8* sunVel=edat->ephemS[ientryS].vel;
    REAL8* sunAcc=edat->ephemS[ientryS].acc;
    REAL8 sunPosNow[3], sunVelNow[3];

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

  return XLAL_SUCCESS;

} /* XLALBarycenterEarth() */


int
XLALBarycenterEarthNew ( EarthState *earth,                /**< [out] the earth's state at time tGPS */
                         const LIGOTimeGPS *tGPS,          /**< [in] GPS time tgps */
                         const EphemerisData *edat,        /**< [in] ephemeris-files */
                         const TimeCorrectionData *tdat,   /**< [in] time correction file data */
                         TimeCorrectionType ttype          /**< [in] time correction type */
                       ){
  earth->ttype = ttype; /* set TimeCorrectionType */

  /* if wanted ORIGNAL version of the code just use XLALBarycenterEarth */
  if ( ttype == TIMECORRECTION_ORIGINAL )
    return XLALBarycenterEarth( earth, tGPS, edat );

  /* start off doing the same as XLALBarycenterEarth() */
    REAL8 tgps[2];   /*I convert from two-integer representation to
                      two REAL8s (just because I initially wrote my code for
                      REAL8s). GPS time(sec) = tgps[0] + 1.e-9*tgps[1] */

    /*int ihour;                 most recent hour to current time, tGPS,
                               in look-up table */
    /*REAL8 t0,t0hr; */

    REAL8 t0e;        /*time since first entry in Earth ephem. table */
    INT4 ientryE;      /*entry in look-up table closest to current time, tGPS */

    REAL8 tinitE;           /*time (GPS) of first entry in Earth ephem table*/
    REAL8 tdiffE;           /*current time tGPS minus time of nearest entry
                             in Earth ephem look-up table */
    REAL8 tdiff2E;          /*tdiff2=tdiffE*tdiffE */
    /*next entries are same as above, but for Sun table */
    REAL8 t0s;
    INT4 ientryS;
    REAL8 tinitS;
    REAL8 tdiffS;
    REAL8 tdiff2S;

    static REAL8 scorr; /* SI second/metre correction factor */
    
    INT4 j; /*dummy index */
  
    /* check input */
    if ( !earth || !tGPS || !edat || !edat->ephemE || !edat->ephemS || !tdat || !tdat->timeCorrs ) {
      XLALPrintError ("%s: invalid NULL input 'earth', 'tGPS', 'edat', 'edat->ephemE', 'edat->ephemS', 'tdat' or 'tdat->timeCorrs'\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL );
    }
    tgps[0] = (REAL8)tGPS->gpsSeconds; /*convert from INT4 to REAL8 */
    tgps[1] = (REAL8)tGPS->gpsNanoSeconds;

    tinitE = edat->ephemE[0].gps;
    tinitS = edat->ephemS[0].gps;

    t0e=tgps[0]-tinitE;
    ientryE = floor((t0e/edat->dtEtable) +0.5e0);  /*finding Earth table entry
                                                closest to arrival time*/
    t0s=tgps[0]-tinitS;
    ientryS = floor((t0s/edat->dtStable) +0.5e0);  /*finding Sun table entry
                                                closest to arrival time*/
    /*Making sure tgps is within earth and sun ephemeris arrays*/
    if ( ( ientryE < 0 ) || ( ientryE >=  edat->nentriesE )) {
      XLALPrintError ("%s: input GPS time %f outside of Earth ephem range [%f, %f]\n", __func__, tgps[0], tinitE, tinitE + edat->nentriesE * edat->dtEtable );
      XLAL_ERROR ( XLAL_EDOM );
    }
    if ( ( ientryS < 0 ) || ( ientryS >=  edat->nentriesS ) ){
      XLALPrintError ("%s: input GPS time %f outside of Sun ephem range [%f, %f]\n", __func__, tgps[0], tinitS, tinitS + edat->nentriesS * edat->dtStable );
      XLAL_ERROR ( XLAL_EDOM );
    }

    tdiffE = t0e -edat->dtEtable*ientryE + tgps[1]*1.e-9; /*tdiff is arrival
              time minus closest Earth table entry; tdiff can be pos. or neg. */
    tdiff2E = tdiffE*tdiffE;

    tdiffS=t0s -edat->dtStable*ientryS + tgps[1]*1.e-9; /*same for Sun*/
    tdiff2S=tdiffS*tdiffS;

    /* in the TCB system the gravitational potential well of the solar system is
     * removed, so clocks run slightly faster than the SI second by a small
     * correction factor */
    if( ttype == TIMECORRECTION_TEMPO2 || ttype == TIMECORRECTION_TCB ) scorr = IFTE_K;
    else scorr = 1.;
    /********************************************************************
     *Calculate position and vel. of center of Earth.
     *We extrapolate from a table produced using a JPL ephemeris.
     *---------------------------------------------------------------------
     */
    {
      REAL8* pos=edat->ephemE[ientryE].pos; /*Cartesian coords of center of Earth
                                              from ephem, in sec. 0=x,1=y,2=z */
      REAL8* vel=edat->ephemE[ientryE].vel;
      REAL8* acc=edat->ephemE[ientryE].acc;
      
      for (j=0;j<3;j++){
        earth->posNow[j] = scorr * (pos[j] + vel[j]*tdiffE + 0.5*acc[j]*tdiff2E);
        earth->velNow[j] = scorr * (vel[j] + acc[j]*tdiffE);
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
      REAL8 eps0= OBLQ;/*obliquity of ecliptic at JD 245145.0*/

      leapsSince2000 = 0; /*right # for Jan. 1, 2000 */
      leapsSince2000 = XLALGPSLeapSeconds( tGPS->gpsSeconds ) - 13;
      {
        INT4 err = xlalErrno;
        if ( err != XLAL_SUCCESS ) {
          XLALPrintError ("%s: XLALGPSLeapSeconds (%d) failed.\n", __func__, tGPS->gpsSeconds );
          XLAL_ERROR ( XLAL_EINVAL );
        }
      }

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

      dtu= tuJC - tu0JC; /*time, in Jul. centuries, between pulse arrival time and previous midnight (UT1) */

      daysSinceJ2000 = (tuInt -43200.)/8.64e4; /*days(not int) on GPS
                                            (or TAI,etc.) clock since J2000 */

      /*----------------------------------------------------------------------
       * in below formula for gmst0 & gmst, on rhs we SHOULD use time elapsed
       * in UT1, but instead we insert elapsed in UTC. This will lead to
       * errors in tdb of order 1 microsec.
       *----------------------------------------------------------------------
       */

      gmst0=24110.54841e0 + tu0JC*(8640184.812866e0 + tu0JC*(0.093104e0 -tu0JC*6.2e-6));
      /*formula for Greenwich mean sidereal time
      on p.50 of Explanatory Supp. to Ast Alm. ...
      gmst unit is sec, where 2pi rad = 24*60*60sec */

      /*-----------------------------------------------------------------------
      * Explan of tu0JC*8640184.812866e0 term: when tu0JC0=0.01 Julian centuries
      * (= 1 year), Earth has spun an extra 86401.84 sec (approx one revolution)
      * with respect to the stars.
      *
      * Note: this way of organizing calc of gmst was stolen from TEMPO;
      * see subroutine lmst (= local mean sidereal time) .
      *------------------------------------------------------------------------
      */

      gmst=gmst0+ dtu*(8.64e4*36525. + 8640184.812866e0
                      +0.093104e0*(tuJC + tu0JC)
                      -6.2e-6*(tuJC*tuJC + tuJC*tu0JC + tu0JC* tu0JC));

      gmstRad=gmst*LAL_PI/43200.;

      earth->gmstRad = gmstRad;

      /*------------------------------------------------------------------------
       * Now calculating 3 terms, describing  lunisolar precession;
       * see Eqs. 3.212 on pp. 104-5 in Explanatory Supp.
       *-------------------------------------------------------------------------
       */

      earth->tzeA = tuJC*(2306.2181e0 + (0.30188e0 + 0.017998e0*tuJC)*tuJC )
                *LAL_PI/6.48e5;

      earth->zA = tuJC*(2306.2181e0 + (1.09468e0 + 0.018203e0*tuJC)*tuJC )
                *LAL_PI/6.48e5;

      earth->thetaA = tuJC*(2004.3109e0 - (0.42665e0 + 0.041833*tuJC)*tuJC )
                *LAL_PI/6.48e5;

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
     * Now calculating Einstein delay using the look-up table and linearly 
     * interpolating between points. Depending on the ephemeris type this is 
     * either the difference between TDB and TDT, or TCB and TDT. The 
     * TCB calculation includes a linear drift term correcting the SSB time
     * to be outside the solar systems gravitational potential. These 
     * calculations are taken from the TEMPO2 tt2tdb.C file.
     * 
     * Note the 51.184 s difference between TDT and GPS
     *-----------------------------------------------------------------------
     */
    {
      /* get deltaT from the look-up table */
      INT4 cidx = floor( ((tgps[0] + tgps[1]*1.e-9) - tdat->timeCorrStart)/tdat->dtTtable );
      
      if( cidx < 0 || cidx > (INT4)tdat->nentriesT-2 ){
        XLALPrintError ("%s: input GPS time %f outside of time ephem range\n", __func__, tgps[0]);
        XLAL_ERROR ( XLAL_EDOM );
      }

      REAL8 dtidx = (tgps[0] + tgps[1]*1.e-9) - (tdat->timeCorrStart + (REAL8)cidx*tdat->dtTtable);
      REAL8 grad = (tdat->timeCorrs[cidx+1] - tdat->timeCorrs[cidx]) / tdat->dtTtable;
    
      REAL8 deltaT = tdat->timeCorrs[cidx] + grad*dtidx;
      
      REAL8 correctionTT_Teph; 
      
      correctionTT_Teph = IFTE_TEPH0 + deltaT / (1.0-IFTE_LC);
      /* the observatory term (obsTerm) correction from TEMPO2's tt2tdb.C code
         is added in XLALBarycenter */
      if( ttype == TIMECORRECTION_TEMPO || ttype == TIMECORRECTION_TDB ){ /* conversion for TDT to TDB using lookup table */
        /* add correction offset as in TEMPO2 tt2tdb.C to 
         * account for its readdition in correctionTT_Teph */
        correctionTT_Teph -= IFTE_TEPH0 / ( 1.0 - IFTE_LC );
        earth->einstein = correctionTT_Teph;
        earth->deinstein = 0.;
      }
      /* conversion for TDT to TCB using lookup table */
      else if( ttype == TIMECORRECTION_TEMPO2 || ttype == TIMECORRECTION_TCB ){
        /* GPS time in MJD */
        REAL8 mjdtt = 44244. + ((tgps[0] + tgps[1]*1.e-9) + 51.184)/86400.;
        
        earth->einstein = IFTE_KM1 * (mjdtt-IFTE_MJD0)*86400.0 /* linear drift term */
                + IFTE_K * (correctionTT_Teph - (long double)IFTE_TEPH0);
        earth->deinstein = IFTE_KM1; /* account for the drift in the derivative */
      }
      
      /* now get the time derivative using Curt's method: NOTE: this really should be
       checked as the differences may not be negligble */
      REAL8 jedtdt = -7300.5e0 + (tgps[0] + 51.184e0 + tgps[1]*1.e-9)/8.64e4;
      REAL8 jt=jedtdt/3.6525e5; /*converting to TEMPO expansion param
                               = Julian millenium, NOT Julian century */

      /* below, I've just taken derivative of above expression for einstein,
         then commented out terms that contribute less than around 10^{-12}
         to tDotBary */
      earth->deinstein += 1.e-6*(
        1656.674564e0*6283.075849991e0 * cos(6283.075849991e0*jt + 6.240054195e0 )+
        22.417471e0*5753.384884897e0 * cos(5753.384884897e0*jt + 4.296977442e0 )  +
        13.839792e0*12566.151699983e0 * cos(12566.151699983e0*jt + 6.196904410e0 ) +
        4.676740e0*6069.776754553e0 * cos(6069.776754553e0*jt + 4.021195093e0 )   +
        1.554905e0*77713.771467920e0 * cos(77713.771467920e0*jt + 5.198467090e0 )
        )/(8.64e4*3.6525e5);
    }
    
    /********************************************************************
     *Now calculating Earth-Sun separation vector, as needed
     *for Shapiro delay calculation.
     *--------------------------------------------------------------------
     */
    {
      REAL8 rse2;
      REAL8* sunPos=edat->ephemS[ientryS].pos;
      REAL8* sunVel=edat->ephemS[ientryS].vel;
      REAL8* sunAcc=edat->ephemS[ientryS].acc;
      REAL8 sunPosNow[3], sunVelNow[3];

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

  return XLAL_SUCCESS;

} /* XLALBarycenterEarthNew() */


/** \author Curt Cutler
 * \brief Transforms from detector arrival time \f$t_a\f$ in GPS (as specified in the
 * LIGOTimeGPS structure) to pulse emission time \f$t_e\f$, in TDB.
 *
 * Actually, the returned \f$t_e\f$ is
 * the emission time plus the constant light-travel-time from
 * source to SSB.) The inputs to LALBarycenter(), through
 * the BarycenterInput structure, are the source location,
 * detector site identifier, and GPS arrival time.
 * The function returns the  emission time \f$t_e(t_a)\f$, the
 * derivative \f$d t_e/d t_a\f$, and the difference
 * \f$t_e(t_a) - t_a \f$ through the EmissionTime structure.
 * The emission time \f$t_e(t_a)\f$ is returned in the LIGOTimeGPS format,
 * while the other two quantities are REAL8's.
 */
int
XLALBarycenter ( EmissionTime *emit, 			/**< [out] emission-time information */
                 const BarycenterInput *baryinput, 	/**< [in] info about detector and source-location */
                 const EarthState *earth) 		/**< [in] earth-state (from LALBarycenterEarth()) */
{
  REAL8 longitude,latitude,rd;  /*geocentric (not geodetic!!) longitude
                                  and latitude of detector vertex, and
                                  dist. rd from center of Earth (sec). */
  REAL8 OMEGA = 7.29211510e-5;  /*ang. vel. of Earth (rad/sec)*/

  REAL8 alpha,delta;  /* RA and DEC (radians) in
                           ICRS realization of J2000 coords.*/

  REAL8 sinTheta;  /* sin(theta) = sin(pi/2 - delta) */

  REAL8 tgps[2];   /*I convert from two-integer representation to
                      two REAL8s (just because I initially wrote my code for
                      REAL8s). GPS time(sec) = tgps[0] + 1.e-9*tgps[1] */

  REAL8 roemer,droemer;  /*Roemer delay and its time derivative*/
  REAL8 erot,derot;      /*delay from center-of-Earth to detector (sec),
			   and its time deriv */
  REAL8 shapiro, dshapiro; /*Shapiro delay due to Sun, and its time deriv. */
  REAL8 finiteDistCorr, dfiniteDistCorr; /*correction to Roemer delay due
		                           to finite dist D to source;
                                           important for D < 100pc */
  REAL8 s[3]; /*unit vector pointing at source, in J2000 Cartesian coords */
  INT4 j; /*dummy index */

  /* check input */
  if ( !emit || !baryinput || !earth ) {
    XLALPrintError ("%s: invalid NULL input 'baryinput, 'emit' or 'earth'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  tgps[0] = (REAL8)baryinput->tgps.gpsSeconds; /*convert from INT4 to REAL8 */
  tgps[1] = (REAL8)baryinput->tgps.gpsNanoSeconds;

    alpha = baryinput->alpha;
    delta = baryinput->delta;

/* check that alpha and delta are in reasonable range */
    if ( ( fabs(alpha) > LAL_TWOPI) || ( fabs(delta) > LAL_PI_2 ) ){
      XLALPrintError ("%s: alpha = %f outside of [-2pi,2pi] or delta = %f outside of [-pi/2,pi/2]\n", __func__, alpha, delta );
      XLAL_ERROR ( XLAL_EDOM );
    }

    sinTheta=sin(LAL_PI/2.0-delta);
    s[2]=cos(LAL_PI/2.0-delta);   /* s is vector that points towards source */
    s[1]=sinTheta*sin(alpha);  /* in Cartesian coords based on J2000 */
    s[0]=sinTheta*cos(alpha);  /* 0=x,1=y,2=z */

    rd = sqrt( baryinput->site.location[0]*baryinput->site.location[0]
	      +baryinput->site.location[1]*baryinput->site.location[1]
	      +baryinput->site.location[2]*baryinput->site.location[2]);
    if ( rd == 0.0 )
      latitude = LAL_PI_2;
    else
      latitude = LAL_PI_2 - acos(baryinput->site.location[2]/rd);
    longitude= atan2(baryinput->site.location[1],baryinput->site.location[0]);

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
       emit->roemer = roemer;
       emit->droemer = droemer;
  }

  /* get the observatory term (if in TDB) */
  REAL8 obsTerm = 0;    /* observatory term correction from TEMPO2 tt2tb.C */
  if ( earth->ttype != TIMECORRECTION_ORIGINAL )
    {
      REAL8 obsEarth[3];
      observatoryEarth( obsEarth, baryinput->site, &baryinput->tgps, earth->gmstRad, earth->delpsi, earth->deleps );

      for ( j = 0; j<3; j++ ) obsTerm += obsEarth[j]*earth->velNow[j];

      obsTerm /= (1.0-IFTE_LC)*(REAL8)IFTE_K;
    }

       /********************************************************************
	* Now including Earth's rotation
	*---------------------------------------------------------------------
	*/
  {

    REAL8 eps0= OBLQ;/*obliquity of ecliptic at JD 245145.0,
					in radians. NOT! to be confused with
					LAL_EPSILON0_SI = permittivity of free
					space; value from Explan. Supp.
					to Astronom. Almanac */
    /*  REAL8 eps0= (23.e0 + 26.e0/60 + 21.448e0/3.6e3)*LAL_PI/180;  */

   REAL8 NdotD;
   REAL8 cosDeltaSinAlphaMinusZA, cosDeltaCosAlphaMinusZA, sinDelta;
            /* above used in calculating effect of luni-solar precession*/

   REAL8 delXNut,delYNut,delZNut,NdotDNut; /* used in calculating effect
            of forced nutation */

   cosDeltaSinAlphaMinusZA = sin(alpha + earth->tzeA)*cos(delta);

   cosDeltaCosAlphaMinusZA = cos(alpha + earth->tzeA)*cos(earth->thetaA)
        *cos(delta) - sin(earth->thetaA)*sin(delta);

   sinDelta = cos(alpha + earth->tzeA)*sin(earth->thetaA)*cos(delta)
             + cos(earth->thetaA)*sin(delta);

/*now taking NdotD, including lunisolar precession, using
  Eqs. 3.212-2 of Explan. Supp.
  Basic idea for incorporating luni-solar precession
  is to change the (alpha,delta) of source to compensate for
  Earth's time-changing spin axis.
*/

   NdotD = sin(latitude)*sinDelta +cos(latitude)*(
            cos(earth->gastRad+longitude-earth->zA)*cosDeltaCosAlphaMinusZA
           +sin(earth->gastRad+longitude-earth->zA)*cosDeltaSinAlphaMinusZA  );

   erot = rd*NdotD;

   derot = OMEGA*rd*cos(latitude)*(
     - sin(earth->gastRad+longitude-earth->zA)*cosDeltaCosAlphaMinusZA
     + cos(earth->gastRad+longitude-earth->zA)*cosDeltaSinAlphaMinusZA  );


   /*--------------------------------------------------------------------------
    * Now adding approx nutation (= short-period,forced motion, by definition).
    * These two dominant terms, with periods 18.6 yrs (big term) and
    * 0.500 yrs (small term),resp., give nutation to around 1 arc sec; see
    * p. 120 of Explan. Supp. The forced nutation amplitude
    *  is around 17 arcsec.
    *
    * Note the unforced motion or Chandler wobble (called ``polar motion''
    * in Explanatory Supp) is not included here. However its amplitude is
    * order of (and a somewhat less than) 1 arcsec; see plot on p. 270 of
    * Explanatory Supplement to Ast. Alm.
    *
    *Below correction for nutation from Eq.3.225-2 of Explan. Supp.
    *Basic idea is to change the (alpha,delta) of source to
    *compensate for Earth's time-changing spin axis.
    *--------------------------------------------------------------------------
    */

   delXNut = -(earth->delpsi)*(cos(delta)*sin(alpha)*cos(eps0)
                      +sin(delta)*sin(eps0));

   delYNut = cos(delta)*cos(alpha)*cos(eps0)*(earth->delpsi)
               -sin(delta)*(earth->deleps);

   delZNut = cos(delta)*cos(alpha)*sin(eps0)*(earth->delpsi)
               +cos(delta)*sin(alpha)*(earth->deleps);


   NdotDNut= sin(latitude)*delZNut
               + cos(latitude)*cos(earth->gastRad+longitude)*delXNut
               + cos(latitude)*sin(earth->gastRad+longitude)*delYNut;

   erot = erot + rd*NdotDNut;

   derot = derot +OMEGA*rd*(
               -cos(latitude)*sin(earth->gastRad+longitude)*delXNut
               +cos(latitude)*cos(earth->gastRad+longitude)*delYNut );

   emit->erot = erot;
   emit->derot = derot;
/* Note erot has a periodic piece (P=one day) AND a constant piece,
   since z-component (parallel to North pole) of vector from
   Earth-center to detector is constant */

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
      REAL8 rsun = 2.322e0; /*radius of sun in sec */
       seDotN = earth->se[2]*sin(delta)+ (earth->se[0]*cos(alpha)
				     + earth->se[1]*sin(alpha))*cos(delta);
       dseDotN=earth->dse[2]*sin(delta)+(earth->dse[0]*cos(alpha)
				     +earth->dse[1]*sin(alpha))*cos(delta);

       b=sqrt(earth->rse*earth->rse-seDotN*seDotN);
       db=(earth->rse*earth->drse-seDotN*dseDotN)/b;


       if ((b < rsun) && (seDotN < 0.e0))/*if gw travels thru interior of Sun*/       {
        shapiro  = 9.852e-6*log( (LAL_AU_SI/LAL_C_SI)/
	                  (seDotN + sqrt(rsun*rsun + seDotN*seDotN)))
                    +19.704e-6*(1.e0 - b/rsun);
	dshapiro = - 19.704e-6*db/rsun;
       }
       else /*else the usual expression*/
       {
        shapiro  = 9.852e-6*log( (LAL_AU_SI/LAL_C_SI)/(earth->rse +seDotN));
        dshapiro = -9.852e-6*(earth->drse+ dseDotN)/(earth->rse +seDotN);
       }
       
       emit->shapiro = shapiro;
       emit->dshapiro = dshapiro;
    }
    /********************************************************************
     *Now correcting Roemer delay for finite distance to source.
     *Timing corrections are order 10 microsec
     *for sources closer than about 100 pc = 10^10 sec.
     *--------------------------------------------------------------------
     */
    {
      REAL8 r2 = 0.e0; /*squared dist from SSB to center of earth, in sec^2 */
      REAL8 dr2 = 0.e0; /*time deriv of r2 */
       if (baryinput->dInv > 1.0e-11){ /*implement if corr.  > 1 microsec*/
       for (j=0;j<3;j++){
       r2 += earth->posNow[j]*earth->posNow[j];
       dr2 += 2.e0*earth->posNow[j]*earth->velNow[j];
       }
       finiteDistCorr = -0.5e0*(r2 - roemer*roemer)*baryinput->dInv;
       dfiniteDistCorr = -(0.5e0*dr2 - roemer*droemer)*baryinput->dInv;
       }
       else{
	 finiteDistCorr = 0.e0;
         dfiniteDistCorr =0.e0;
       }
    }
    /*-----------------------------------------------------------------------
     *Now adding it all up.
     *emit.te is pulse emission time in TDB coords
     *(up to a constant roughly equal to ligh travel time from source to SSB).
     *emit->deltaT = emit.te - tgps.
     *-----------------------------------------------------------------------
     */
    {

  INT4 deltaTint; /* floor of deltaT */

  emit->deltaT = roemer + erot + obsTerm + earth->einstein - shapiro + finiteDistCorr;

  emit->tDot = 1.e0 + droemer + derot + earth->deinstein - dshapiro + dfiniteDistCorr;

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
       /*Now adding Earth's rotation to rDetector and
         vDetector, but NOT yet including effects of
         lunisolar prec. or nutation (though DO use gast instead of gmst)
         For next 10 years, will give rDetector to 10 microsec and
         v/c to 10^-9 */
       emit->rDetector[0] += rd*cos(latitude)*cos(longitude+earth->gastRad);
       emit->vDetector[0] += -OMEGA*rd*cos(latitude)*
                               sin(longitude+earth->gastRad);
       emit->rDetector[1] += rd*cos(latitude)*sin(longitude+earth->gastRad);
       emit->vDetector[1] += OMEGA*rd*cos(latitude)*
                               cos(longitude+earth->gastRad);
       emit->rDetector[2] += rd*sin(latitude);
       /*no change to emit->vDetector[2] = component along spin axis,
         if ignore prec. and nutation */
   }

    }

    return XLAL_SUCCESS;

} /* XLALBarycenter() */

/** \author Curt Cutler, Miroslav Shaltev, R Prix
 * \brief Speed optimized version of XLALBarycenter(),
 * should be fully equivalent except for the additional buffer argument.
 * The 'buffer' is used to keep sky-specific and detector-specific values that can can potentially be re-used
 * across subsequent calls to this function.
 *
 * NOTE: The 'buffer' argument is a pointer to a buffer-pointer, which is allowed to be buffer==NULL (=> no buffering),
 * otherwise a given non-NULL (*buffer) is re-used, or if (*buffer==NULL) we create a new one and return in (*buffer).
 *
 */
int
XLALBarycenterOpt ( EmissionTime *emit, 		/**< [out] emission-time information */
                    const BarycenterInput *baryinput, 	/**< [in] info about detector and source-location */
                    const EarthState *earth, 		/**< [in] earth-state (from LALBarycenterEarth()) */
                    BarycenterBuffer **buffer		/**< [in/out] internal buffer for speed optimization */
                    )
{
  /* ---------- check input sanity ---------- */
  XLAL_CHECK ( emit != NULL, XLAL_EINVAL, "Invalid input: emit == NULL");
  XLAL_CHECK ( baryinput != NULL, XLAL_EINVAL, "Invalid input: baryinput == NULL");
  XLAL_CHECK ( earth != NULL, XLAL_EINVAL, "Invalid input: earth == NULL");

  // physical constants used by Curt (slightly different from LAL's Constants, but kept for binary-equivalence with XLALBarycenter()
  const REAL8 OMEGA = 7.29211510e-5;  /* ang. vel. of Earth (rad/sec)*/
  // REAL8 eps0 = 0.40909280422232891;	/* obliquity of ecliptic at JD 245145.0, in radians. Value from Explan. Supp. to Astronom. Almanac */
  const REAL8 sinEps0 = 0.397777155931914; 	// sin ( eps0 );
  const REAL8 cosEps0 = 0.917482062069182;	// cos ( eps0 );

  REAL8 tgps[2];
  tgps[0] = baryinput->tgps.gpsSeconds;
  tgps[1] = baryinput->tgps.gpsNanoSeconds;

  // ---------- sky-position dependent quantities: compute or re-use from buffer is applicable
  REAL8 alpha,delta;  /* RA and DEC (radians) in ICRS realization of J2000 coords.*/
  alpha = baryinput->alpha;
  delta = baryinput->delta;

  REAL8 sinAlpha, cosAlpha, sinDelta, cosDelta;
  REAL8 n[3]; /*unit vector pointing from SSB to the source, in J2000 Cartesian coords, 0=x,1=y,2=z */

  // ---------- handle buffering of recurring computed quantities
  // 3 possible scenarios:
  //   i)   caller gave as an initialized buffer: we use that one
  //   ii)  caller gave as a pointer to a NULL-initialized buffer-pointer: allocate buffer, initialize and pass back to caller
  //   iii) caller gave as a NULL pointer: use buffer internally, but destroy at the end of this function
  BarycenterBuffer *myBuffer = NULL;
  if ( (buffer) && (*buffer) ) // caller gave as an allocated buffer
    myBuffer = (*buffer);
  else	// we need to create a buffer
    XLAL_CHECK ( (myBuffer = XLALCalloc(1,sizeof(*myBuffer))) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,sizeof(*myBuffer))\n" );

  // now we're (basically) guaranteed to have a valid buffer in 'myBuffer'

  // use buffered sky-quantities if same sky-position as stored in buffer
  if ( myBuffer->active && (alpha == myBuffer->alpha) && (delta == myBuffer->delta ) )
    {
      sinDelta = myBuffer->fixed_sky.sinDelta;
      cosDelta = myBuffer->fixed_sky.cosDelta;
      sinAlpha = myBuffer->fixed_sky.sinAlpha;
      cosAlpha = myBuffer->fixed_sky.cosAlpha;
      n[0] = myBuffer->fixed_sky.n[0];
      n[1] = myBuffer->fixed_sky.n[1];
      n[2] = myBuffer->fixed_sky.n[2];
    }
  else // not buffered, recompute sky-dependent quantities
    {
      /* check that alpha and delta are in reasonable range */
      XLAL_CHECK ( fabs(alpha) <= LAL_TWOPI, XLAL_EDOM, "alpha = %f outside of allowed range [-2pi,2pi]\n", alpha );
      XLAL_CHECK ( fabs(delta) <= LAL_PI_2,  XLAL_EDOM, "delta = %f outside of allowed range [-pi/2,pi/2]\n", delta );

      sinDelta = cos ( LAL_PI/2.0 - delta );	// this weird way of computing it is required to stay binary identical to Curt's function
      cosDelta = sin ( LAL_PI/2.0 - delta );
      sinAlpha = sin ( alpha );
      cosAlpha = cos ( alpha );

      n[0] = cosDelta * cosAlpha;
      n[1] = cosDelta * sinAlpha;
      n[2] = sinDelta;

      // ... and store them in the myBuffer
      myBuffer->alpha 		 = alpha;
      myBuffer->delta 		 = delta;
      myBuffer->fixed_sky.sinDelta = sinDelta;
      myBuffer->fixed_sky.cosDelta = cosDelta;
      myBuffer->fixed_sky.sinAlpha = sinAlpha;
      myBuffer->fixed_sky.cosAlpha = cosAlpha;
      myBuffer->fixed_sky.n[0] 	 = n[0];
      myBuffer->fixed_sky.n[1] 	 = n[1];
      myBuffer->fixed_sky.n[2] 	 = n[2];
      myBuffer->active = 1;
    } // if not re-using sky-buffered quantities

  // ---------- detector site-position dependent quantities: compute or re-use from buffer is applicable
  REAL8 rd;   /* distance 'rd' of detector from center of Earth, in light seconds */
  REAL8 longitude, latitude; 	/* geocentric (not geodetic!!) longitude and latitude of detector vertex */
  REAL8 sinLat, cosLat;
  REAL8 rd_sinLat, rd_cosLat;	// shortcuts for 'rd * sin(latitude)' and 'rd * cos(latitude)' respectively

  // use buffered site-quantities if same site as stored in myBuffer
  if ( myBuffer->active && (memcmp ( &baryinput->site, &myBuffer->site, sizeof(myBuffer->site) ) == 0) )
    {
      rd	= myBuffer->fixed_site.rd;
      longitude	= myBuffer->fixed_site.longitude;
      latitude	= myBuffer->fixed_site.latitude;
      sinLat	= myBuffer->fixed_site.sinLat;
      cosLat	= myBuffer->fixed_site.cosLat;
      rd_sinLat	= myBuffer->fixed_site.rd_sinLat;
      rd_cosLat	= myBuffer->fixed_site.rd_cosLat;
    }
  else
    {
      rd = sqrt( + baryinput->site.location[0]*baryinput->site.location[0]
                 + baryinput->site.location[1]*baryinput->site.location[1]
                 + baryinput->site.location[2]*baryinput->site.location[2] );

      longitude = atan2 ( baryinput->site.location[1], baryinput->site.location[0] );
      if ( rd == 0.0 )
        latitude = LAL_PI_2;	// avoid division by 0, for detector at center of earth
      else
        latitude = LAL_PI_2 - acos ( baryinput->site.location[2] / rd );

      sinLat = sin ( latitude );
      cosLat = cos ( latitude );
      rd_sinLat = rd * sinLat;
      rd_cosLat = rd * cosLat;

      // ... and store them in the myBuffer
      memcpy ( &myBuffer->site, &baryinput->site , sizeof(myBuffer->site) );
      myBuffer->fixed_site.rd 		= rd;
      myBuffer->fixed_site.longitude 	= longitude;
      myBuffer->fixed_site.latitude 	= latitude;
      myBuffer->fixed_site.sinLat 	= sinLat;
      myBuffer->fixed_site.cosLat 	= cosLat;
      myBuffer->fixed_site.rd_sinLat 	= rd_sinLat;
      myBuffer->fixed_site.rd_cosLat 	= rd_cosLat;
      myBuffer->active = 1;
    } // if not re-using site-buffered quantities

  /*---------------------------------------------------------------------
   * Calucate Roemer delay for detector at center of Earth.
   * We extrapolate from a table produced using JPL DE405 ephemeris.
   *---------------------------------------------------------------------
   */
  /*Roemer delay and its time derivative*/
  emit->roemer = 0;
  emit->droemer = 0;
  for ( UINT4 j = 0; j < 3; j++ )
    {
      emit->roemer  += n[j] * earth->posNow[j];
      emit->droemer += n[j] * earth->velNow[j];
    }

  /* get the observatory term (if in TDB) */
  REAL8 obsTerm = 0;    /* observatory term correction from TEMPO2 tt2tb.C */
  if ( earth->ttype != TIMECORRECTION_ORIGINAL )
    {
      REAL8 obsEarth[3];
      observatoryEarth( obsEarth, baryinput->site, &baryinput->tgps, earth->gmstRad, earth->delpsi, earth->deleps );

      for ( UINT4 j = 0; j < 3; j++ )
        obsTerm += obsEarth[j] * earth->velNow[j];

      obsTerm /= (1.0-IFTE_LC)*(REAL8)IFTE_K;
    }

  /*---------------------------------------------------------------------
   * Now including Earth's rotation
   *---------------------------------------------------------------------
   */
  /* calculating effect of luni-solar precession */
  REAL8 sinAlphaMinusZA = sin ( alpha + earth->tzeA );
  REAL8 cosAlphaMinusZA = cos ( alpha + earth->tzeA );
  REAL8 cosThetaA = cos ( earth->thetaA );
  REAL8 sinThetaA = sin ( earth->thetaA );

  REAL8 cosDeltaSinAlphaMinusZA = sinAlphaMinusZA * cosDelta;

  REAL8 cosDeltaCosAlphaMinusZA = cosAlphaMinusZA * cosThetaA * cosDelta - sinThetaA * sinDelta;

  REAL8 sinDeltaCurt = cosAlphaMinusZA * sinThetaA * cosDelta + cosThetaA * sinDelta;

  /* now taking NdotD, including lunisolar precession, using
     Eqs. 3.212-2 of Explan. Supp.
     Basic idea for incorporating luni-solar precession
     is to change the (alpha,delta) of source to compensate for
     Earth's time-changing spin axis.
  */
  REAL8 cosGastZA = cos ( earth->gastRad + longitude-earth->zA );
  REAL8 sinGastZA = sin ( earth->gastRad + longitude-earth->zA );

  REAL8 rd_NdotD = rd_sinLat * sinDeltaCurt + rd_cosLat * ( cosGastZA * cosDeltaCosAlphaMinusZA + sinGastZA * cosDeltaSinAlphaMinusZA );

  /* delay from center-of-Earth to detector (sec), and its time deriv */
  emit->erot = rd_NdotD;
  emit->derot = OMEGA * rd_cosLat * ( - sinGastZA * cosDeltaCosAlphaMinusZA + cosGastZA * cosDeltaSinAlphaMinusZA );

  /*--------------------------------------------------------------------------
   * Now adding approx nutation (= short-period,forced motion, by definition).
   * These two dominant terms, with periods 18.6 yrs (big term) and
   * 0.500 yrs (small term),resp., give nutation to around 1 arc sec; see
   * p. 120 of Explan. Supp. The forced nutation amplitude
   *  is around 17 arcsec.
   *
   * Note the unforced motion or Chandler wobble (called ``polar motion''
   * in Explanatory Supp) is not included here. However its amplitude is
   * order of (and a somewhat less than) 1 arcsec; see plot on p. 270 of
   * Explanatory Supplement to Ast. Alm.
   *
   * Below correction for nutation from Eq.3.225-2 of Explan. Supp.
   * Basic idea is to change the (alpha,delta) of source to
   * compensate for Earth's time-changing spin axis.
   *--------------------------------------------------------------------------
   */
  REAL8 delXNut = - earth->delpsi * ( cosDelta * sinAlpha * cosEps0 + sinDelta * sinEps0 );

  REAL8 delYNut = cosDelta * cosAlpha * cosEps0 * earth->delpsi - sinDelta * earth->deleps;

  REAL8 delZNut = cosDelta * cosAlpha * sinEps0 * earth->delpsi + cosDelta * sinAlpha * earth->deleps;

  REAL8 cosGastLong = cos ( earth->gastRad + longitude );
  REAL8 sinGastLong = sin ( earth->gastRad + longitude );

  REAL8 rd_NdotDNut = rd_sinLat * delZNut + rd_cosLat * cosGastLong * delXNut + rd_cosLat * sinGastLong * delYNut;

  emit->erot += rd_NdotDNut;
  emit->derot += OMEGA * ( - rd_cosLat * sinGastLong * delXNut + rd_cosLat * cosGastLong * delYNut );

  /* Note erot has a periodic piece (P=one day) AND a constant piece,
     since z-component (parallel to North pole) of vector from
     Earth-center to detector is constant
  */

  /*--------------------------------------------------------------------
   * Now adding Shapiro delay. Note according to J. Taylor review article
   * on pulsar timing, max value of Shapiro delay (when rays just graze sun)
   * is 120 microsec.
   *
   * Here we calculate Shapiro delay
   * for a detector at the center of the Earth.
   * Causes errors of order 10^{-4}sec * 4 * 10^{-5} = 4*10^{-9} sec
   *--------------------------------------------------------------------
   */
  REAL8 rsun = 2.322; /*radius of sun in sec */
  REAL8 seDotN  = earth->se[2] * sinDelta + ( earth->se[0]  * cosAlpha + earth->se[1] * sinAlpha ) * cosDelta;
  REAL8 dseDotN = earth->dse[2]* sinDelta + ( earth->dse[0] * cosAlpha + earth->dse[1] * sinAlpha ) * cosDelta;

  REAL8 b = sqrt ( earth->rse * earth->rse - seDotN * seDotN );
  REAL8 db = ( earth->rse * earth->drse - seDotN * dseDotN ) / b;

  /* if gw travels thru interior of Sun*/
  if ( ( b < rsun ) && ( seDotN < 0 ) )
    {
      emit->shapiro  = 9.852e-6 * log ( (LAL_AU_SI/LAL_C_SI) / ( seDotN + sqrt ( rsun*rsun + seDotN*seDotN ) ) ) + 19.704e-6 * ( 1.0 - b / rsun );
      emit->dshapiro = - 19.704e-6 * db / rsun;
    }
  else /* else the usual expression*/
    {
      emit->shapiro  =  9.852e-6 * log( (LAL_AU_SI/LAL_C_SI) / ( earth->rse + seDotN ) );
      emit->dshapiro = -9.852e-6 * ( earth->drse + dseDotN ) / ( earth->rse + seDotN );
    }

  /*--------------------------------------------------------------------
   * Now correcting Roemer delay for finite distance to source.
   * Timing corrections are order 10 microsec
   * for sources closer than about 100 pc = 10^10 sec.
   *--------------------------------------------------------------------
   */
  REAL8 r2 = 0; 	/* squared dist from SSB to center of earth, in sec^2 */
  REAL8 dr2 = 0; 	/* time deriv of r2 */
  REAL8 finiteDistCorr, dfiniteDistCorr; /*correction to Roemer delay due to finite dist D to source; important for D < 100pc */

  if (baryinput->dInv > 1.0e-11)	/* implement if corr.  > 1 microsec */
    {
      for ( UINT4 j=0; j<3; j++ )
        {
          r2  += earth->posNow[j] * earth->posNow[j];
          dr2 += 2.0 * earth->posNow[j] * earth->velNow[j];
        }
      finiteDistCorr  = - 0.5 * ( r2 - emit->roemer * emit->roemer ) * baryinput->dInv;
      dfiniteDistCorr = - ( 0.5 * dr2 - emit->roemer * emit->droemer ) * baryinput->dInv;
    }
  else
    {
      finiteDistCorr = 0;
      dfiniteDistCorr = 0;
    }

  /* -----------------------------------------------------------------------
   * Now adding it all up.
   * emit.te is pulse emission time in TDB coords
   * (up to a constant roughly equal to ligh travel time from source to SSB).
   * emit->deltaT = emit.te - tgps.
   * -----------------------------------------------------------------------
   */
  emit->deltaT = emit->roemer + emit->erot + earth->einstein - emit->shapiro + finiteDistCorr + obsTerm;
  emit->tDot = 1.0 + emit->droemer + emit->derot + earth->deinstein - emit->dshapiro + dfiniteDistCorr;

  INT4 deltaTint = floor ( emit->deltaT );

  if ( ( 1e-9 * tgps[1] + emit->deltaT - deltaTint ) >= 1.e0 )
    {
      emit->te.gpsSeconds     = baryinput->tgps.gpsSeconds + deltaTint + 1;
      emit->te.gpsNanoSeconds = floor ( 1e9 * ( tgps[1] * 1e-9 + emit->deltaT - deltaTint - 1.0 ) );
    }
  else
    {
      emit->te.gpsSeconds     = baryinput->tgps.gpsSeconds + deltaTint;
      emit->te.gpsNanoSeconds = floor ( 1e9 * ( tgps[1] * 1e-9 + emit->deltaT - deltaTint ) );
    }

  for ( UINT4 j=0; j<3; j++)
    {
      emit->rDetector[j] = earth->posNow[j];
      emit->vDetector[j] = earth->velNow[j];
    }

  /* Now adding Earth's rotation to rDetector and
     vDetector, but NOT yet including effects of
     lunisolar prec. or nutation (though DO use gast instead of gmst)
     For next 10 years, will give rDetector to 10 microsec and
     v/c to 10^-9
  */
  emit->rDetector[0] += rd_cosLat * cosGastLong;
  emit->vDetector[0] += - OMEGA * rd_cosLat * sinGastLong;
  emit->rDetector[1] += rd_cosLat * sinGastLong;
  emit->vDetector[1] += OMEGA * rd_cosLat * cosGastLong;
  emit->rDetector[2] += rd_sinLat;
  /*no change to emit->vDetector[2] = component along spin axis, if ignore prec. and nutation */

  // finish buffer handling
  // *) if user passed a pointer to NULL, return our internal buffer
  if ( buffer && (*buffer == NULL ) )
    (*buffer) = myBuffer;
  // *) if passed NULL as 'buffer', then we destroy the internal buffer now
  if ( buffer == NULL )
    XLALFree ( myBuffer );

  return XLAL_SUCCESS;

} /* XLALBarycenterOpt() */



/* ==================== deprecated LAL interface (only wrappers to XLAL-fcts now) ==================== */

/** Deprecated LAL wrapper to XLALBarycenterEarth()
 * \deprecated use XLALBarycenterEarth() instead.
 */
void
LALBarycenterEarth(LALStatus *status,		/**< [in/out] LAL status structure pointer */
		   EarthState *earth, 		/**< [out] the earth's state at time tGPS */
		   const LIGOTimeGPS *tGPS, 	/**< [in] GPS time tgps */
		   const EphemerisData *edat) 	/**< [in] ephemeris-files */
{
  INITSTATUS(status);

  if ( XLALBarycenterEarth ( earth, tGPS, edat) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALBarycenterEarth() failed with xlalErrno = %d\n", __func__, xlalErrno );
    ABORT ( status, LALBARYCENTERH_EXLAL, LALBARYCENTERH_MSGEXLAL);
  }

  RETURN(status);

} /* LALBarycenterEarth() */

/** Deprecated LAL wrapper to XLALBarycenter()
 * \deprecated use XLALBarycenter() instead.
 */
void
LALBarycenter(LALStatus *status,		/**< [in/out] LAL status structure pointer */
	      EmissionTime *emit, 		/**< [out] emission-time information */
	      const BarycenterInput *baryinput, /**< [in] info about detector and source-location */
	      const EarthState *earth) 		/**< [in] earth-state (from LALBarycenterEarth()) */
{
  INITSTATUS(status);

  if ( XLALBarycenter ( emit, baryinput, earth ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALBarycenter() failed with xlalErrno = %d\n", __func__, xlalErrno );
    ABORT ( status, LALBARYCENTERH_EXLAL, LALBARYCENTERH_MSGEXLAL);
  }

  RETURN(status);

} /* LALBarycenter() */

/** Function to calculate the precession matrix give Earth nutation values
  * depsilon and dpsi for a given MJD time.
  *
  * This function is a slightly modified version of the TEMPO2 function 
  * get_precessionMatrix in get_ObsCoord.C (itself translated from the 
  * TEMPO FORTRAN functions for precession and nutation). */
void precessionMatrix( REAL8 prn[3][3], /**< [out] The precession matrix */
                       REAL8 mjd,       /**< [in] The UT1 local mean siderial time (in MJD format) */ 
                       REAL8 dpsi,      /**< [in] The dpsi value for the Earth nutation */
                       REAL8 deps )     /**< [in] The deps value for the Earth nutation */
{
  INT4 i, j;

  /* For precession */
  REAL8 t, trad, zeta, z, theta;
  static const REAL8 par_zeta[3] = {2306.2181, 0.30188, 0.017998};
  static const REAL8 par_z[3] = {2306.2181, 1.09468, 0.018203};
  static const REAL8 par_theta[3] = {2004.3109, -0.42665, -0.041833};
  static const REAL8 seconds_per_rad = 3600.0/LAL_PI_180;
  
  REAL8 czeta, szeta, cz, sz, ctheta, stheta;

  REAL8 nut[3][3], prc[3][3];

  /* cosine and sine of the obliquity of the ecliptic */
  const REAL8 ceps = 0.917482062069182;
  const REAL8 seps = 0.397777155931914;

  t = (mjd - 51544.5)/36525.0;
  trad = t / seconds_per_rad;
  
  zeta = trad * (par_zeta[0]+t*(par_zeta[1]+t*par_zeta[2]));
  z = trad * (par_z[0]+t*(par_z[1]+t*par_z[2]));
  theta = trad * (par_theta[0]+t*(par_theta[1]+t*par_theta[2]));

  czeta = cos(zeta);
  szeta = sin(zeta);
  cz = cos(z);
  sz = sin(z);
  ctheta = cos(theta);
  stheta = sin(theta);

  prc[0][0] = czeta*ctheta*cz - szeta*sz;
  prc[1][0] = czeta*ctheta*sz + szeta*cz;
  prc[2][0] = czeta*stheta;
  prc[0][1] = -szeta*ctheta*cz - czeta*sz;
  prc[1][1] = -szeta*ctheta*sz + czeta*cz;
  prc[2][1] = -szeta*stheta;
  prc[0][2] = -stheta*cz;
  prc[1][2] = -stheta*sz;
  prc[2][2] = ctheta;

  nut[0][0] = 1.0;
  nut[0][1] = -dpsi*ceps;
  nut[0][2] = -dpsi*seps;
  nut[1][0] = -nut[0][1];
  nut[1][1] = 1.0;
  nut[1][2] = -deps;
  nut[2][0] = -nut[0][2];
  nut[2][1] = -nut[1][2];
  nut[2][2] = 1.0;

  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      prn[j][i] = nut[i][0]*prc[0][j] + nut[i][1]*prc[1][j] + nut[i][2]*prc[2][j];
}

/** Function to get the observatory site location with respect to the
  * centre of the Earth, taking into account precession and nutation.
  * This is based on the routines in the TEMPO2 get_ObsCoord.C code */ 
void observatoryEarth( REAL8 obsEarth[3],     /**< [out] The x, y, z coordinates of the observatory from the centre of the Earth */
                       const LALDetector det, /**< [in] The LAL detector site */
                       const LIGOTimeGPS *tgps,     /**< [in] The GPS time at the detector */
                       REAL8 gmst,            /**< [in] The Greenwich Mean Sidereal Time */
                       REAL8 dpsi,            /**< [in] dpsi for Earth nutation */
                       REAL8 deps             /**< [in] deps for Earth nutation */
                      ){        
  static REAL8 erad; /* observatory distance from Earth centre */
  static REAL8 hlt;  /* observatory latitude */
  static REAL8 alng; /* observatory longitude */
  REAL8 tmjd = 44244. + ( XLALGPSGetREAL8( tgps ) + 51.184 )/86400.;
  
  INT4 j = 0;
  
  /* note that det.location is in light-seconds */
  erad = sqrt( det.location[0]*det.location[0]
             + det.location[1]*det.location[1]
             + det.location[2]*det.location[2] );
  if ( erad == 0.0 ) hlt = LAL_PI_2;
  else hlt = asin(det.location[2] / erad);
  
  alng = atan2(-det.location[1], det.location[0]);

  static REAL8 siteCoord[3];
  REAL8 eeq[3], prn[3][3];

  siteCoord[0] = erad * cos(hlt);
  siteCoord[1] = siteCoord[0]*tan(hlt);
  siteCoord[2] = alng;

  /* PC,PS equatorial and meridional components of nutations of longitude */      
  REAL8 toblq = (tmjd - 5.15445e4)/36525.0;
  REAL8 oblq = (((1.813e-3*toblq-5.9e-4)*toblq-4.6815e1)*toblq +84381.448)*LAL_PI_180/3600.0;

  REAL8 pc = cos(oblq+deps)*dpsi;

  /* Compute the local, true sidereal time */
  REAL8 ph = gmst + pc - siteCoord[2];  
  
  /* Get X-Y-Z coordinates taking out earth rotation */
  eeq[0] = siteCoord[0]*cos(ph); 
  eeq[1] = siteCoord[0]*sin(ph);
  eeq[2] = siteCoord[1];

  /* Now obtain PRN -- the precession matrix */
  precessionMatrix(prn, tmjd, dpsi, deps);
  /* Calculate the position after precession/nutation*/
  for ( j = 0; j < 3 ; j++ )
    obsEarth[j] = prn[j][0]*eeq[0] + prn[j][1]*eeq[1] + prn[j][2]*eeq[2];
}
