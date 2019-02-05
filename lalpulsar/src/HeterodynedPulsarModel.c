/*
 * Copyright (C) 2018 Matthew Pitkin
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

#include <lal/HeterodynedPulsarModel.h>

/** Macro to round a value to the nearest integer. */
#define ROUND(x) (floor(x+0.5))

/** Macro to square a value. */
#define SQUARE(x) ( (x) * (x) )

 /**
 * \brief The phase evolution difference compared to a heterodyned phase (for a pulsar) 
 *
 * This function will calculate the difference in the phase evolution of a
 * source at a particular sky location as observed in a detector on Earth
 * compared with that used to heterodyne (or SpectralInterpolate) the data. If
 * the phase evolution is described by a Taylor expansion:
 * \f[
 * \phi(T) = \sum_{k=1}^n \frac{f^{(k-1)}}{k!} T^k,
 * \f]
 * where \f$f^{(x)}\f$ is the xth time derivative of the gravitational-wave
 * frequency, and \f$T\f$ is the pulsar proper time, then the phase difference
 * is given by
 * \f[
 * \Delta\phi(t) = \sum_{k=1}^n \left( \frac{\Delta f^{(k-1)}}{k!}(t+\delta t_1)^k + \frac{f^{(k-1)}_2}{k!} \sum_{i=0}^{i<k} \left(\begin{array}{c}k \\ i\end{array}\right) (\Delta t)^{k-i} (t+\delta t_1)^i \right),
 * \f]
 * where \f$t\f$ is the signal arrival time at the detector minus the given
 * pulsar period epoch, \f$\delta t_1\f$ is the barycentring time delay (from
 * both solar system and binary orbital effects) calculated at the heterodyned
 * values, \f$\Delta f^{(x)} = f_2^{(x)}-f1^{(x)}\f$ is the diffence in
 * frequency (derivative) between the current value (\f$f_2^{(x)}\f$) and the
 * heterodyne value (\f$f_1^{(x)}\f$), and
 * \f$\Delta t = \delta t_2 - \delta t_1\f$ is the difference between the
 * barycentring time delay calculated at the current values (\f$\delta t_1\f$)
 * and the heterodyned values. Frequency time derivatives up to any order are
 * allowed. The pulsar proper time is calculated by correcting the time of
 * arrival at Earth, \f$t\f$ to the solar system barycentre and if necessary
 * the binary system barycenter, so
 * \f$T = t + \delta{}t_{\rm SSB} + \delta{}t_{\rm BSB}\f$.
 *
 * In this function the time delay needed to correct to the solar system
 * barycenter is only calculated if required, i.e., if an update is required
 * due to a change in the sky position. The same is true for the binary system
 * time delay, which is only calculated if it needs updating due to a change in
 * the binary system parameters.
 *
 * This function can just return the phase evolution (and not phase evolution
 * difference) if \c ssbdts and \c bsbdts are NULL, and no \c DELTAF values are
 * set in the \c params structure (or if the \c DELTAF values are all set to
 * zero).
 *
 * \param params [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times at which to calculate the phase difference
 * \param freqfactor [in] The multiplicative factor on the pulsar frequency for a particular model
 * \param ssbdts [in] The vector of SSB time delays used for the original heterodyne. If this is
 * \c NULL then the SSB delay for the given position will be calculated.
 * \param updateSSBDelay [in] Set to a non-zero value if the SSB delay needs to be recalculated at
 * an updated sky position compared to that used for the heterodyne.
 * \param bsbdts [in] The vector of BSB time delays used for the original heterodyne. If this is
 * \c NULL then the BSB delay for the given binary parameters will be calculated.
 * \param updateBSBDelay [in] Set to a non-zero value if the BSB delay needs to be recalulated at
 * a set of updated binary system parameters.
 * \param detector [in] A pointer to a \c LALDetector structure for a particular detector
 * \param ephem [in] A pointer to an \c EphemerisData structure containing solar system ephemeris information
 * \param tdat [in] A pointer to a \c TimeCorrectionData structure containing time system correction information
 * \param ttype [in] The \c TimeCorrectionType value
 *
 * \return A vector of rotational phase difference values (in cycles NOT radians)
 *
 * \sa XLALHeterodynedPulsarGetSSBDelay
 * \sa XLALHeterodynedPulsarGetBSBDelay
 */
REAL8Vector *XLALHeterodynedPulsarPhaseDifference( PulsarParameters *params,
                                                   const LIGOTimeGPSVector *datatimes,
                                                   REAL8 freqfactor,
                                                   REAL8Vector *ssbdts,
                                                   UINT4 updateSSBDelay,
                                                   REAL8Vector *bsbdts,
                                                   UINT4 updateBSBDelay,
                                                   const LALDetector *detector,
                                                   const EphemerisData *ephem,
                                                   const TimeCorrectionData *tdat,
                                                   TimeCorrectionType ttype ){
  /* check inputs */
  XLAL_CHECK_NULL( params != NULL, XLAL_EFUNC, "PulsarParameters must not be NULL" );
  XLAL_CHECK_NULL( datatimes != NULL, XLAL_EFUNC, "datatimes must not be NULL" );
  XLAL_CHECK_NULL( freqfactor > 0., XLAL_EFUNC, "freqfactor must be greater than zero" );
  XLAL_CHECK_NULL( detector != NULL, XLAL_EFUNC, "LALDetector must not be NULL" );
  XLAL_CHECK_NULL( ephem != NULL, XLAL_EFUNC, "EphemerisData must not be NULL" );

  UINT4 i = 0, j = 0, k = 0, length = 0, isbinary = 0;

  REAL8 DT = 0., deltat = 0., deltatpow = 0., deltatpowinner = 1., taylorcoeff = 1., Ddelay = 0., Ddelaypow = 0.;

  REAL8Vector *phis = NULL, *dts = NULL, *fixdts = NULL, *bdts = NULL, *fixbdts = NULL;

  REAL8 pepoch = PulsarGetREAL8ParamOrZero(params, "PEPOCH"); /* time of ephem info */
  REAL8 cgw = PulsarGetREAL8ParamOrZero(params, "CGW");
  REAL8 T0 = pepoch;

  /* glitch parameters */
  REAL8 *glep = NULL, *glph = NULL, *glf0 = NULL, *glf1 = NULL, *glf2 = NULL, *glf0d = NULL, *gltd = NULL, *deltafs = NULL;
  UINT4 glnum = 0;

  length = datatimes->length;

  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );

  /* get solar system barycentring time delays */
  if ( ssbdts == NULL ){
    /* calculate SSB delay at the given parameters */
    fixdts = XLALHeterodynedPulsarGetSSBDelay( params, datatimes, detector, ephem, tdat, ttype );
  }
  else{
    fixdts = ssbdts; /* use SSB delays passed to function as calculated from the heterodyne parameters */
    if( updateSSBDelay ){ /* get SSB delay for updated sky position parameters */
      dts = XLALHeterodynedPulsarGetSSBDelay( params, datatimes, detector, ephem, tdat, ttype );
    }
  }
  XLAL_CHECK_NULL( length == fixdts->length, XLAL_EFUNC, "Lengths of time stamp vector and SSB delay vector are not the same" );

  if ( PulsarCheckParam( params, "BINARY" ) ){ 
    isbinary = 1; /* see if pulsar is in binary */

    if ( bsbdts == NULL ){
      /* calculate BSB delay at the given parameters */
      if ( dts != NULL ){ fixbdts = XLALHeterodynedPulsarGetBSBDelay( params, datatimes, dts, ephem ); }
      else { fixbdts = XLALHeterodynedPulsarGetBSBDelay( params, datatimes, fixdts, ephem ); }
    }
    else{
      fixbdts = bsbdts; /* use BSB delays passed to function as calculated from the heterodyne parameters */
      if ( updateBSBDelay || updateSSBDelay ){
        if ( dts != NULL ){ bdts = XLALHeterodynedPulsarGetBSBDelay( params, datatimes, dts, ephem ); }
        else { bdts = XLALHeterodynedPulsarGetBSBDelay( params, datatimes, fixdts, ephem ); }
      }
    }
    XLAL_CHECK_NULL( length == fixbdts->length, XLAL_EFUNC, "Lengths of time stamp vector and BSB delay vector are not the same" );
  }

  /* get vector of frequencies and frequency differences */
  const REAL8Vector *freqs = PulsarGetREAL8VectorParam( params, "F" );
  deltafs = XLALCalloc( freqs->length, sizeof(REAL8) );  // initialise deltafs to zero
  if ( PulsarCheckParam( params, "DELTAF" ) ){
    const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "DELTAF" );

    XLAL_CHECK_NULL( tmpvec->length == freqs->length, XLAL_EFUNC, "Number of frequencies is different from number of delta fs" );

    for ( i=0; i<tmpvec->length; i++ ){ deltafs[i] = tmpvec->data[i]; }
  }
  else{
    /* set deltafs to (negative) frequencies */
    for ( i=0; i<freqs->length; i++ ){ deltafs[i] = -freqs->data[i]; }
  }

  if ( PulsarCheckParam( params, "GLEP" ) ){ /* see if pulsar has glitch parameters */
    const REAL8Vector *glpars = PulsarGetREAL8VectorParam( params, "GLEP" );
    glnum = glpars->length;

    /* get epochs */
    glep = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    for ( i=0; i<glpars->length; i++ ){ glep[i] = glpars->data[i]; }

    /* get phase offsets */
    glph = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    if ( PulsarCheckParam( params, "GLPH" ) ){
      const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "GLPH" );
      for ( i=0; i<tmpvec->length; i++ ){ glph[i] = tmpvec->data[i]; }
    }

    /* get frequencies offsets */
    glf0 = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    if ( PulsarCheckParam( params, "GLF0" ) ){
      const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "GLF0" );
      for ( i=0; i<tmpvec->length; i++ ){ glf0[i] = tmpvec->data[i]; }
    }

    /* get frequency derivative offsets */
    glf1 = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    if ( PulsarCheckParam( params, "GLF1" ) ){
      const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "GLF1" );
      for ( i=0; i<tmpvec->length; i++ ){ glf1[i] = tmpvec->data[i]; }
    }

    /* get second frequency derivative offsets */
    glf2 = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    if ( PulsarCheckParam( params, "GLF2" ) ){
      const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "GLF2" );
      for ( i=0; i<tmpvec->length; i++ ){ glf2[i] = tmpvec->data[i]; }
    }

    /* get decaying frequency component offset derivative */
    glf0d = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    if ( PulsarCheckParam( params, "GLF0D" ) ){
      const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "GLF0D" );
      for ( i=0; i<tmpvec->length; i++ ){ glf0d[i] = tmpvec->data[i]; }
    }

    /* get decaying frequency component decay time constant */
    gltd = XLALCalloc(glnum, sizeof(REAL8)); /* initialise to zeros */
    if ( PulsarCheckParam( params, "GLTD" ) ){
      const REAL8Vector *tmpvec = PulsarGetREAL8VectorParam( params, "GLTD" );
      for ( i=0; i<tmpvec->length; i++ ){ gltd[i] = tmpvec->data[i]; }
    }
  }

  for( i=0; i<length; i++){
    REAL8 deltaphi = 0., innerphi = 0.; /* change in phase */
    Ddelay = 0.;                        /* change in SSB/BSB delay */

    REAL8 realT = XLALGPSGetREAL8( &datatimes->data[i] ); /* time of data */
    DT = realT - T0; /* time diff between data and start of data */

    /* get difference in solar system barycentring time delays */
    if ( dts != NULL ){ Ddelay += ( dts->data[i] - fixdts->data[i] ); }
    deltat = DT + fixdts->data[i];

    if ( isbinary ){
      /* get difference in binary system barycentring time delays */
      if ( bdts != NULL && fixbdts != NULL ) { Ddelay += ( bdts->data[i] - fixbdts->data[i] ); }
      deltat += fixbdts->data[i];
    }

    /* correct for speed of GW compared to speed of light */
    if ( cgw > 0.0 && cgw < 1. ) {
      deltat /= cgw;
      Ddelay /= cgw;
    }

    /* get the change in phase (compared to the heterodyned phase) */
    deltatpow = deltat;
    for ( j=0; j<freqs->length; j++ ){
      taylorcoeff = gsl_sf_fact(j+1);
      deltaphi += deltafs[j]*deltatpow/taylorcoeff;
      if ( Ddelay != 0. ){
        innerphi = 0.;
        deltatpowinner = 1.; /* this starts as one as it is first raised to the power of zero */
        Ddelaypow = pow(Ddelay, j+1);
        for ( k=0; k<j+1; k++ ){
          innerphi += gsl_sf_choose(j+1, k) * Ddelaypow * deltatpowinner;
          deltatpowinner *= deltat; /* raise power */
          Ddelaypow /= Ddelay;      /* reduce power */
        }
        deltaphi += innerphi*freqs->data[j]/taylorcoeff;
      }
      deltatpow *= deltat;
    }

    /* check for glitches */
    if ( glnum > 0 ){
      /* get glitch phase - based on equations in formResiduals.C of TEMPO2 from Eqn. 1 of Yu et al (2013) http://ukads.nottingham.ac.uk/abs/2013MNRAS.429..688Y */
      for ( j=0; j<glnum; j++ ){
        if ( deltat >= (glep[j]-T0) ){
          REAL8 dtg = 0, expd = 1.;
          dtg = deltat - (glep[j]-T0); /* time since glitch */
          if ( gltd[j] != 0. ) { expd = exp(-dtg/gltd[j]); } /* decaying part of glitch */
          deltaphi += glph[j] + glf0[j]*dtg + 0.5*glf1[j]*dtg*dtg + (1./6.)*glf2[j]*dtg*dtg*dtg + glf0d[j]*gltd[j]*(1.-expd);
        }
      }
    }

    deltaphi *= freqfactor; /* multiply by frequency factor */
    phis->data[i] = deltaphi - floor(deltaphi); /* only need to keep the fractional part of the phase */
  }

  /* free memory */
  if ( dts != NULL ){ XLALDestroyREAL8Vector( dts ); }
  if ( bdts != NULL ){ XLALDestroyREAL8Vector( bdts ); }
  if ( ssbdts == NULL ){ XLALDestroyREAL8Vector( fixdts ); }
  if ( bsbdts == NULL ){ XLALDestroyREAL8Vector( fixbdts ); }
  XLALFree(deltafs);

  if ( glnum > 0 ){
    XLALFree( glep );
    XLALFree( glph );
    XLALFree( glf0 );
    XLALFree( glf1 );
    XLALFree( glf2 );
    XLALFree( glf0d );
    XLALFree( gltd );
  }

  return phis;
}


/**
 * \brief Computes the delay between a GPS time at Earth and the solar system barycentre
 *
 * This function calculates the time delay between a GPS time at a specific
 * location (e.g., a gravitational-wave detector) on Earth and the solar system
 * barycentre. The delay consists of three components: the geometric time delay
 * (Roemer delay) \f$t_R = \mathbf{r}(t)\hat{n}/c\f$ (where \f$\mathbf{r}(t)\f$
 * is the detector's position vector at time \f$t\f$), the special relativistic
 * Einstein delay \f$t_E\f$, and the general relativistic Shapiro delay
 * \f$t_S\f$.
 *
 * \param pars [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times at Earth
 * \param detector [in] Information on the detector position on the Earth
 * \param ephem [in] Information on the solar system ephemeris
 * \param tdat [in] Informaion on the time system corrections
 * \param ttype [in] The type of time system corrections to perform
 *
 * \return A vector of time delays in seconds
 *
 * \sa XLALBarycenter
 * \sa XLALBarycenterEarthNew
 */
REAL8Vector *XLALHeterodynedPulsarGetSSBDelay( PulsarParameters *pars,
                                               const LIGOTimeGPSVector *datatimes,
                                               const LALDetector *detector,
                                               const EphemerisData *ephem,
                                               const TimeCorrectionData *tdat,
                                               TimeCorrectionType ttype ){
  /* check inputs */
  XLAL_CHECK_NULL( pars != NULL, XLAL_EFUNC, "PulsarParameters must not be NULL" );
  XLAL_CHECK_NULL( datatimes != NULL, XLAL_EFUNC, "datatimes must not be NULL" );
  XLAL_CHECK_NULL( detector != NULL, XLAL_EFUNC, "LALDetector must not be NULL" );
  XLAL_CHECK_NULL( ephem != NULL, XLAL_EFUNC, "EphemerisData must not be NULL" );

  INT4 i = 0, length = 0;

  BarycenterInput bary;

  REAL8Vector *dts = NULL;

  /* copy barycenter and ephemeris data */
  bary.site.location[0] = detector->location[0]/LAL_C_SI;
  bary.site.location[1] = detector->location[1]/LAL_C_SI;
  bary.site.location[2] = detector->location[2]/LAL_C_SI;

  REAL8 ra = 0.;
  if ( PulsarCheckParam( pars, "RA" ) ) { ra = PulsarGetREAL8Param( pars, "RA" ); }
  else if ( PulsarCheckParam( pars, "RAJ" ) ) { ra = PulsarGetREAL8Param( pars, "RAJ" ); }
  else {
    XLAL_ERROR_NULL( XLAL_EINVAL, "No source right ascension specified!" );
  }
  REAL8 dec = 0.;
  if ( PulsarCheckParam( pars, "DEC" ) ) { dec = PulsarGetREAL8Param( pars, "DEC" ); }
  else if ( PulsarCheckParam( pars, "DECJ" ) ) { dec = PulsarGetREAL8Param( pars, "DECJ" ); }
  else {
    XLAL_ERROR_NULL( XLAL_EINVAL, "No source declination specified!" );
  }
  REAL8 pmra = PulsarGetREAL8ParamOrZero( pars, "PMRA" );
  REAL8 pmdec = PulsarGetREAL8ParamOrZero( pars, "PMDEC" );
  REAL8 pepoch = PulsarGetREAL8ParamOrZero( pars, "PEPOCH" );
  REAL8 posepoch = PulsarGetREAL8ParamOrZero( pars, "POSEPOCH" );
  REAL8 px = PulsarGetREAL8ParamOrZero( pars, "PX" );     /* parallax */

   /* set the position and frequency epochs if not already set */
  if( pepoch == 0. && posepoch != 0.) { pepoch = posepoch; }
  else if( posepoch == 0. && pepoch != 0. ) { posepoch = pepoch; }

  length = datatimes->length;

  /* allocate memory for times delays */
  dts = XLALCreateREAL8Vector( length );

  /* set 1/distance if parallax value is given (1/sec) */
  if( px != 0. ) { bary.dInv = px*(LAL_C_SI/LAL_AU_SI); }
  else { bary.dInv = 0.; }

  /* make sure ra and dec are wrapped within 0--2pi and -pi.2--pi/2 respectively */
  ra = fmod(ra, LAL_TWOPI);
  REAL8 absdec = fabs(dec);
  if ( absdec > LAL_PI_2 ){
    UINT4 nwrap = floor((absdec+LAL_PI_2)/LAL_PI);
    dec = (dec > 0 ? 1. : -1.)*(nwrap%2 == 1 ? -1. : 1.)*(fmod(absdec + LAL_PI_2, LAL_PI) - LAL_PI_2);
    ra = fmod(ra + (REAL8)nwrap*LAL_PI, LAL_TWOPI); /* move RA by pi */
  }

  EarthState earth;
  EmissionTime emit;
  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &datatimes->data[i] );

    bary.tgps = datatimes->data[i];
    bary.delta = dec + ( realT - posepoch ) * pmdec;
    bary.alpha = ra + ( realT - posepoch ) * pmra / cos( bary.delta );

    /* call barycentring routines */
    XLAL_CHECK_NULL( XLALBarycenterEarthNew( &earth, &bary.tgps, ephem, tdat, ttype ) == XLAL_SUCCESS, XLAL_EFUNC, "Barycentring routine failed" );
    XLAL_CHECK_NULL( XLALBarycenter( &emit, &bary, &earth ) == XLAL_SUCCESS, XLAL_EFUNC, "Barycentring routine failed" );

    dts->data[i] = emit.deltaT;
  }

  return dts;
}


/**
 * \brief Computes the delay between a pulsar in a binary system and the barycentre of the system
 *
 * This function uses \c XLALBinaryPulsarDeltaTNew to calculate the time delay
 * between for a pulsar in a binary system between the time at the pulsar and
 * the time at the barycentre of the system. This includes Roemer delays and
 * relativistic delays. The orbit may be described by different models and can
 * be purely Keplarian or include various relativistic corrections.
 *
 * \param pars [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times
 * \param dts [in] A vector of solar system barycentre time delays
 * \param edat [in] Solar system ephemeris information
 *
 * \return A vector of time delays in seconds
 *
 * \sa XLALBinaryPulsarDeltaTNew
 */
REAL8Vector *XLALHeterodynedPulsarGetBSBDelay( PulsarParameters *pars,
                                               const LIGOTimeGPSVector *datatimes,
                                               const REAL8Vector *dts,
                                               const EphemerisData *edat ){
  /* check inputs */
  XLAL_CHECK_NULL( pars != NULL, XLAL_EFUNC, "PulsarParameters must not be NULL" );
  XLAL_CHECK_NULL( datatimes != NULL, XLAL_EFUNC, "datatimes must not be NULL" );
  XLAL_CHECK_NULL( dts != NULL, XLAL_EFUNC, "SSB delay times must not be NULL" );
  XLAL_CHECK_NULL( edat != NULL, XLAL_EFUNC, "EphemerisData must not be NULL" );

  REAL8Vector *bdts = NULL;
  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;
  EarthState earth;

  INT4 i = 0, length = datatimes->length;

  /* check whether there's a binary model */
  if ( PulsarCheckParam( pars, "BINARY" ) ){
    bdts = XLALCreateREAL8Vector( length );

    for ( i = 0; i < length; i++ ){
      binput.tb = XLALGPSGetREAL8( &datatimes->data[i] ) + dts->data[i];

      XLALGetEarthPosVel( &earth, edat, &datatimes->data[i] );

      binput.earth = earth; /* current Earth state */
      XLALBinaryPulsarDeltaTNew( &boutput, &binput, pars );
      bdts->data[i] = boutput.deltaT;
    }
  }
  return bdts;
}


/**
 * \brief Get the position and velocity of the Earth at a given time
 *
 * This function will get the position and velocity of the Earth from the
 * ephemeris data at the time t. It will be returned in an EarthState
 * structure. This is based on the start of the \c XLALBarycenterEarth
 * function.
 *
 * \param earth [out] The returned \c EarthState structure containing the Earth's position and velocity
 * \param edat [in] The solar system ephemeris data
 * \param tGPS [in] A \c LIGOTimeGPS structure containing a GPS time
 *
 * \sa XLALBarycenterEarth
 */
void XLALGetEarthPosVel( EarthState *earth,
                         const EphemerisData *edat,
                         const LIGOTimeGPS *tGPS ){
  REAL8 tgps[2];

  REAL8 t0e;        /* time since first entry in Earth ephem. table */
  INT4 ientryE;     /* entry in look-up table closest to current time, tGPS */

  REAL8 tinitE;     /* time (GPS) of first entry in Earth ephem table */
  REAL8 tdiffE;     /* current time tGPS minus time of nearest entry in Earth ephem look-up table */
  REAL8 tdiff2E;    /* tdiff2 = tdiffE * tdiffE */

  INT4 j;

  /* check input */
  if ( !earth || !tGPS || !edat || !edat->ephemE || !edat->ephemS ) {
    XLAL_ERROR_VOID( XLAL_EINVAL, "Invalid NULL input 'earth', 'tGPS', 'edat','edat->ephemE' or 'edat->ephemS'" );
  }

  tgps[0] = (REAL8)tGPS->gpsSeconds; /* convert from INT4 to REAL8 */
  tgps[1] = (REAL8)tGPS->gpsNanoSeconds;

  tinitE = edat->ephemE[0].gps;

  t0e = tgps[0] - tinitE;
  ientryE = ROUND(t0e/edat->dtEtable);  /* finding Earth table entry */

  if ( ( ientryE < 0 ) || ( ientryE >=  edat->nentriesE )) {
    XLAL_ERROR_VOID( XLAL_EDOM, "Input GPS time %f outside of Earth ephem range [%f, %f]\n", tgps[0], tinitE, tinitE +
edat->nentriesE * edat->dtEtable );
  }

  /* tdiff is arrival time minus closest Earth table entry; tdiff can be pos. or neg. */
  tdiffE = t0e -edat->dtEtable*ientryE + tgps[1]*1.e-9;
  tdiff2E = tdiffE*tdiffE;

  REAL8* pos = edat->ephemE[ientryE].pos;
  REAL8* vel = edat->ephemE[ientryE].vel;
  REAL8* acc = edat->ephemE[ientryE].acc;

  for (j=0;j<3;j++){
    earth->posNow[j]=pos[j] + vel[j]*tdiffE + 0.5*acc[j]*tdiff2E;
    earth->velNow[j]=vel[j] + acc[j]*tdiffE;
  }
}


/**
 * \brief The amplitude model of a complex heterodyned signal from the
 * \f$l=2, m=1,2\f$ harmonics of a rotating neutron star.
 *
 * This function calculates the complex heterodyned time series model for a
 * rotating neutron star. It will currently calculate the model for emission
 * from the \f$l=m=2\f$ harmonic (which gives emission at twice the rotation
 * frequency) and/or the \f$l=2\f$ and \f$m=1\f$ harmonic (which gives emission
 * at the rotation frequency). See LIGO T1200265-v3. Further harmonics can be
 * added and are defined by the \c freqFactor value, which is the multiple of
 * the spin-frequency at which emission is produced.
 *
 * The antenna pattern functions are contained in a 1D lookup table, so within
 * this function the correct value for the given time is interpolated from this
 * lookup table using linear interpolation.
 *
 * \param pars [in] A set of pulsar parameters
 * \param freqfactor [in] The harmonic frequency of the signal in units of the
 * pulsar rotation frequency
 * \param varyphase [in] Set to a non-zero value is the signal phase is
 * different to the heterodyne phase (or if wanting the signal output at all
 * time stamps).
 * \param useroq [in] Set to a non-zero value if a reduced order quadrature
 * likelihood is being used
 * \param nonGR [in] Set to a non-zero value to indicate a non-GR polarisation
 * is used
 * \param timestamps [in] A vector of GPS times at which to calculate the signal
 * \param resp [in] A detector response function look-up table
 * 
 * \return A \c COMPLEX16TimeSeries the amplitude time series
 */
COMPLEX16TimeSeries* XLALHeterodynedPulsarGetAmplitudeModel( PulsarParameters *pars,
                                                             REAL8 freqfactor,
                                                             UINT4 varyphase,
                                                             UINT4 useroq,
                                                             UINT4 nonGR,
                                                             const LIGOTimeGPSVector *timestamps,
                                                             const DetResponseTimeLookupTable *resp ){
  COMPLEX16TimeSeries *csignal = NULL;

  /* create signal if not already allocated */
  LIGOTimeGPS gpstime; /* dummy time stamp */
  gpstime.gpsSeconds = 0;
  gpstime.gpsNanoSeconds = 0;

  if ( varyphase || useroq ){
    XLAL_CHECK_NULL( timestamps->length > 0, XLAL_EFUNC, "length must be a positive integer" );

    /* in these cases use signal length to specify the length of the signal */
    csignal = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, timestamps->length );
  }
  else{
    if ( nonGR ){
      /* for non-GR signals there are six values that need to be held */
      csignal = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 6 );
    }
    else{
      /* for GR signals just two values are required */
      csignal = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 2 );
    }
  }

  /* switch to waveform parameters if required */
  XLALPulsarSourceToWaveformParams( pars );

  /* check inputs */
  XLAL_CHECK_NULL( pars != NULL, XLAL_EFUNC, "PulsarParameters is NULL" );
  XLAL_CHECK_NULL( freqfactor > 0., XLAL_EFUNC, "freqfactor must be a positive number" );
  XLAL_CHECK_NULL( resp != NULL, XLAL_EFUNC, "Response look-up table is NULL");
  XLAL_CHECK_NULL( freqfactor != 1. || nonGR != 1, XLAL_EFUNC, "A non-GR model is not defined for a frequency harmonic of 1" );

  UINT4 i = 0;

  REAL8 T, twopsi, t0;
  REAL8 cosiota = PulsarGetREAL8ParamOrZero( pars, "COSIOTA" );
  REAL8 siniota = sin(acos(cosiota));
  REAL8 s2psi = 0., c2psi = 0., spsi = 0., cpsi = 0.;

  twopsi = 2.*PulsarGetREAL8ParamOrZero( pars, "PSI" );
  s2psi = sin(twopsi);
  c2psi = cos(twopsi);

  t0 = XLALGPSGetREAL8( &resp->t0 );  /* epoch of detector response look-up table */

  if ( nonGR == 1 ){
    spsi = sin(PulsarGetREAL8ParamOrZero( pars, "PSI" ));
    cpsi = cos(PulsarGetREAL8ParamOrZero( pars, "PSI" ));
  }

  COMPLEX16 expPhi;
  COMPLEX16 Cplus = 0., Ccross = 0., Cx = 0., Cy = 0., Cl = 0., Cb = 0.;

  /* get the amplitude and phase factors */
  if( freqfactor == 1. ){
    /* the l=2, m=1 harmonic at the rotation frequency. NOTE: there is currently no
       amplitude model defined for a non-GR signal at this harmonic */
    expPhi = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PHI21" ));
    Cplus = -0.25 * PulsarGetREAL8ParamOrZero( pars, "C21" ) * siniota * cosiota * expPhi;
    Ccross = 0.25 * I * PulsarGetREAL8ParamOrZero( pars, "C21" ) * siniota * expPhi;
  }
  else if( freqfactor == 2. ){
    /* the l=2, m=2 harmonic at twice the rotation frequency */
    if ( nonGR ){ /* amplitude if nonGR is specifiec */
      COMPLEX16 expPhiTensor, expPsiTensor, expPhiScalar, expPsiScalar, expPhiVector, expPsiVector;

      expPhiTensor = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PHI0TENSOR" ) );
      expPsiTensor = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PSITENSOR" ) );
      expPhiScalar = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PHI0SCALAR" ) );
      expPsiScalar = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PSISCALAR" ) );
      expPhiVector = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PHI0VECTOR" ) );
      expPsiVector = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PSIVECTOR" ) );

      Cplus = 0.5 * expPhiTensor * PulsarGetREAL8ParamOrZero( pars, "HPLUS" );
      Ccross = 0.5 * expPhiTensor * PulsarGetREAL8ParamOrZero( pars, "HCROSS" ) * expPsiTensor;
      Cx = 0.5 * expPhiVector * PulsarGetREAL8ParamOrZero( pars, "HVECTORX" );
      Cy = 0.5 * expPhiVector * PulsarGetREAL8ParamOrZero( pars, "HVECTORY" ) * expPsiVector;
      Cb = 0.5 * expPhiScalar * PulsarGetREAL8ParamOrZero( pars, "HSCALARB" );
      Cl = 0.5 * expPhiScalar * PulsarGetREAL8ParamOrZero( pars, "HSCALARL" ) * expPsiScalar;
    }
    else{ /* just GR tensor mode amplitudes */
      expPhi = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PHI22" ) );
      Cplus = -0.5 * PulsarGetREAL8ParamOrZero( pars, "C22" ) * ( 1. + cosiota * cosiota ) * expPhi;
      Ccross = I * PulsarGetREAL8ParamOrZero( pars, "C22" ) * cosiota * expPhi;
    }
  }
  else{
    XLAL_ERROR_NULL( XLAL_EINVAL, "Error... currently unknown frequency factor (%.2lf) for models.", freqfactor );
  }

  if ( varyphase || useroq ){ /* have to compute the full time domain signal */
    REAL8 tsv;

    /* set lookup table parameters */
    tsv = LAL_DAYSID_SI / (REAL8)resp->ntimebins;

    for( i=0; i<csignal->data->length; i++ ){
      REAL8 plus00, plus01, cross00, cross01, plus = 0., cross = 0.;
      REAL8 x00, x01, y00, y01, b00, b01, l00, l01;
      REAL8 timeScaled, timeMin, timeMax;
      REAL8 plusT = 0., crossT = 0., x = 0., y = 0., xT = 0., yT = 0., b = 0., l = 0.;
      INT4 timebinMin, timebinMax;

      /* set the time bin for the lookup table */
      /* sidereal day in secs*/
      T = fmod( XLALGPSGetREAL8( &timestamps->data[i] ) - t0, LAL_DAYSID_SI );  /* convert GPS time to sidereal day */
      timebinMin = (INT4)fmod( floor(T / tsv), resp->ntimebins );
      timeMin = timebinMin*tsv;
      timebinMax = (INT4)fmod( timebinMin + 1, resp->ntimebins );
      timeMax = timeMin + tsv;

      /* get values of matrix for linear interpolation */
      plus00 = resp->fplus->data[timebinMin];
      plus01 = resp->fplus->data[timebinMax];

      cross00 = resp->fcross->data[timebinMin];
      cross01 = resp->fcross->data[timebinMax];

      /* rescale time for linear interpolation on a unit square */
      timeScaled = (T - timeMin)/(timeMax - timeMin);

      plus = plus00 + (plus01-plus00)*timeScaled;
      cross = cross00 + (cross01-cross00)*timeScaled;

      plusT = plus*c2psi + cross*s2psi;
      crossT = cross*c2psi - plus*s2psi;

      if ( nonGR ){
        x00 = resp->fx->data[timebinMin];
        x01 = resp->fx->data[timebinMax];
        y00 = resp->fy->data[timebinMin];
        y01 = resp->fy->data[timebinMax];
        b00 = resp->fb->data[timebinMin];
        b01 = resp->fb->data[timebinMax];
        l00 = resp->fl->data[timebinMin];
        l01 = resp->fl->data[timebinMax];

        x = x00 + (x01-x00)*timeScaled;
        y = y00 + (y01-y00)*timeScaled;
        b = b00 + (b01-b00)*timeScaled;
        l = l00 + (l01-l00)*timeScaled;

        xT = x*cpsi + y*spsi;
        yT = y*cpsi - x*spsi;
      }

      /* create the complex signal amplitude model appropriate for the harmonic */
      csignal->data->data[i] = ( Cplus * plusT ) + ( Ccross * crossT );

      /* add non-GR components if required */
      if ( nonGR ){ csignal->data->data[i] += ( Cx*xT ) + ( Cy*yT ) + Cb*b + Cl*l; }
    }
  }
  else{ /* just have to calculate the values to multiply the pre-summed data */
    /* for tensor-only models (e.g. the default of GR) calculate the two components of
     * the single model value - things multiplied by a(t) and things multiplied by b(t)
     * (both these will have real and imaginary components). NOTE: the values input into
     * csignal->data->data are not supposed to be identical to the above
     * relationships between the amplitudes and polarisation angles, as these are the
     * multiplicative coefficients of the a(t) and b(t) summations.
     */

    /* put multiples of a(t) in first value and b(t) in second */
    if ( !nonGR ){
      csignal->data->data[0] = (Cplus*c2psi - Ccross*s2psi);
      csignal->data->data[1] = (Cplus*s2psi + Ccross*c2psi);
    }
    else{
      csignal->data->data[0] = (Cplus*c2psi - Ccross*s2psi);
      csignal->data->data[1] = (Cplus*s2psi + Ccross*c2psi);
      csignal->data->data[2] = (Cx*cpsi - Cy*spsi);
      csignal->data->data[3] = (Cx*spsi + Cy*cpsi);
      csignal->data->data[4] = Cb;
      csignal->data->data[5] = Cl;
    }
  }

  return csignal;
}


/**
 * \brief Generate the model of the neutron star signal
 *
 * Firstly the time varying amplitude of the signal will be calculated based
 * on the antenna pattern and amplitude parameters. Then, if searching over
 * phase parameters, the phase evolution of the signal will be calculated. The
 * difference between the new phase model, \f$\phi(t)_n\f$, and that used to
 * heterodyne the data, \f$\phi(t)_h\f$, will be calculated and the complex
 * signal model, \f$M\f$, modified accordingly:
 * \f[
 * M'(t) = M(t)\exp{i((\phi(t)_n - \phi(t)_h))}.
 * \f]
 * This does not try to undo the signal modulation in the data, but instead
 * replicates the modulation in the model, hence the positive phase difference
 * rather than a negative phase in the exponential function.
 *
 * The model is described in Equations 7 and 8 of \cite Pitkin2017 .
 *
 * \param pars [in] A \c PulsarParameters structure containing the model parameters
 * \param freqfactor [in] The harmonic frequency of the signal in units of the
 * pulsar rotation frequency
 * \param usephase [in] Set to a non-zero value is the signal phase is
 * different to the heterodyne phase (or if wanting the signal output at all
 * time stamps).
 * \param useroq [in] Set to a non-zero value if a reduced order quadrature
 * likelihood is being used
 * \param nonGR [in] Set to a non-zero value to indicate a non-GR polarisation
 * is used
 * \param timestamps [in] A vector of GPS times at which to calculate the signal
 * \param hetssbdelays [in] The vector of SSB time delays used for the original heterodyne. If this is
 * \c NULL then the SSB delay for the given position will be calculated.
 * \param updateSSBDelay [in] Set to a non-zero value if the SSB delay needs to be recalculated at
 * an updated sky position compared to that used for the heterodyne.
 * \param hetbsbdelays [in] The vector of BSB time delays used for the original heterodyne. If this is
 * \c NULL then the BSB delay for the given binary parameters will be calculated.
 * \param updateBSBDelay [in] Set to a non-zero value if the BSB delay needs to be recalulated at
 * a set of updated binary system parameters.
 * \param resp [in] A loop-up table for the detector response function.
 * \param ephem [in] A pointer to an \c EphemerisData structure containing solar system ephemeris information
 * \param tdat [in] A pointer to a \c TimeCorrectionData structure containing time system correction information
 * \param ttype [in] The \c TimeCorrectionType value
 * 
 * \sa XLALHeterodynedPulsarGetAmplitudeModel
 * \sa XLALHeterodynedPulsarGetPhaseModel
 */
COMPLEX16TimeSeries* XLALHeterodynedPulsarGetModel( PulsarParameters *pars,
                                                    REAL8 freqfactor,
                                                    UINT4 usephase,
                                                    UINT4 useroq,
                                                    UINT4 nonGR,
                                                    const LIGOTimeGPSVector *timestamps,
                                                    REAL8Vector *hetssbdelays,
                                                    UINT4 updateSSBDelay,
                                                    REAL8Vector *hetbsbdelays,
                                                    UINT4 updateBSBDelay,
                                                    const DetResponseTimeLookupTable *resp,
                                                    const EphemerisData *ephem,
                                                    const TimeCorrectionData *tdat,
                                                    TimeCorrectionType ttype ){
  UINT4 i = 0;
  COMPLEX16TimeSeries *csignal = NULL;

  /* get amplitude model */
  csignal = XLALHeterodynedPulsarGetAmplitudeModel( pars,
                                                    freqfactor,
                                                    usephase,
                                                    useroq,
                                                    nonGR,
                                                    timestamps,
                                                    resp );

  // include phase change if required
  if ( usephase ){
    REAL8Vector *dphi = NULL;
    dphi = XLALHeterodynedPulsarPhaseDifference( pars,
                                                 timestamps,
                                                 freqfactor,
                                                 hetssbdelays,
                                                 updateSSBDelay,
                                                 hetbsbdelays,
                                                 updateBSBDelay,
                                                 resp->det,
                                                 ephem,
                                                 tdat,
                                                 ttype );

    COMPLEX16 M = 0., expp = 0.;

    for ( i = 0; i < dphi->length; i++ ){
      /* phase factor by which to multiply the (almost) DC signal model. NOTE: this does not try to undo
       * the signal modulation in the data, but instead replicates it in the model, hence the positive
       * phase rather than a negative phase in the cexp function. */
      expp = cexp( LAL_TWOPI * I * dphi->data[i] );

      M = csignal->data->data[i];

      /* heterodyne */
      csignal->data->data[i] = M * expp;
    }

    XLALDestroyREAL8Vector( dphi );
  }

  return csignal;
}


/**
 * \brief Creates a lookup table of the detector antenna pattern
 *
 * This function creates a lookup table of the tensor, vector and scalar
 * antenna patterns for a given detector orientation and source sky position.
 * For the tensor modes these are the functions given by equations 10-13 in
 * \cite JKS98 , whilst for the vector and scalar modes they are the \f$\psi\f$
 * independent parts of e.g. equations 5-8 of \cite Nishizawa2009 . We remove
 * the \f$\psi\f$ dependance by setting \f$\psi=0\f$.
 *
 * If \c avedt is a value over 60 seconds then the antenna pattern will
 * actually be the mean value from 60 second intervals within that timespan.
 * This accounts for the fact that each data point is actually an averaged
 * value over the given timespan.
 *
 * \param t0 [in] initial GPS time of the data
 * \param det [in] structure containing the detector
 * \param alpha [in] the source right ascension in radians
 * \param delta [in] the source declination in radians
 * \param timeSteps [in] the number of grid bins to use in time
 * \param avedt [in] average the antenna pattern over this timespan
 */
DetResponseTimeLookupTable* XLALDetResponseLookupTable( REAL8 t0,
                                                        const LALDetector *det,
                                                        REAL8 alpha,
                                                        REAL8 delta,
                                                        UINT4 timeSteps,
                                                        REAL8 avedt ){
  LIGOTimeGPS gps;

  // allocate structure
  DetResponseTimeLookupTable *resp;
  resp = XLALMalloc(sizeof(*resp));
  resp->det = XLALMalloc(sizeof(LALDetector*));
  memcpy(&resp->det, &det, sizeof(det));
  resp->alpha = alpha;
  resp->delta = delta;
  resp->psi = 0.;
  resp->ntimebins = timeSteps;
  resp->fplus = XLALCreateREAL8Vector(timeSteps);
  resp->fcross = XLALCreateREAL8Vector(timeSteps);
  resp->fx = XLALCreateREAL8Vector(timeSteps);
  resp->fy = XLALCreateREAL8Vector(timeSteps);
  resp->fb = XLALCreateREAL8Vector(timeSteps);
  resp->fl = XLALCreateREAL8Vector(timeSteps);
  XLALGPSSetREAL8(&resp->t0, t0);  // set the epoch

  REAL8 T = 0., Tstart = 0., Tav = 0., psi = 0.;

  REAL8 fplus = 0., fcross = 0., fx = 0., fy = 0., fb = 0., fl = 0.;
  REAL8 tsteps = (REAL8)timeSteps;

  UINT4 j = 0, k = 0, nav = 0;

  /* number of points to average */
  if ( avedt == 60. ) { nav = 1; }
  else{ nav = floor(avedt/60.) + 1; }

  for( j = 0 ; j < timeSteps ; j++ ){
    resp->fplus->data[j] = 0.;
    resp->fcross->data[j] = 0.;
    resp->fx->data[j] = 0.;
    resp->fy->data[j] = 0.;
    resp->fb->data[j] = 0.;
    resp->fl->data[j] = 0.;

    /* central time of lookup table point */
    T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

    if ( nav % 2 ){ /* is odd */
      Tstart = T - 0.5*((REAL8)nav-1.)*60.;
    }
    else{ /* is even */
      Tstart = T - (0.5*(REAL8)nav - 1.)*60. - 30.;
    }

    for ( k = 0; k < nav; k++ ){
      Tav = Tstart + 60.*(REAL8)k;

      XLALGPSSetREAL8(&gps, Tav);
      XLALComputeDetAMResponseExtraModes( &fplus, &fcross, &fb, &fl, &fx, &fy, det->response,
                                          alpha, delta, psi,
                                          XLALGreenwichMeanSiderealTime( &gps ) );

      resp->fplus->data[j] += fplus;
      resp->fcross->data[j] += fcross;
      resp->fx->data[j] += fx;
      resp->fy->data[j] += fy;
      resp->fb->data[j] += fb;
      resp->fl->data[j] += fl;
    }

    resp->fplus->data[j] /= (REAL8)nav;
    resp->fcross->data[j] /= (REAL8)nav;
    resp->fx->data[j] /= (REAL8)nav;
    resp->fy->data[j] /= (REAL8)nav;
    resp->fb->data[j] /= (REAL8)nav;
    resp->fl->data[j] /= (REAL8)nav;
  }

  return resp;
}


/** \brief Free memory for antenna response look-up table structure
 *
 * \param resp [in] the look-up table structure to be freed
 */ 
void XLALDestroyDetResponseTimeLookupTable(DetResponseTimeLookupTable* resp){
  if ( resp == NULL ){ return; }

  if ( resp->fplus ){
    XLALDestroyREAL8Vector( resp->fplus );
  }
  if ( resp->fcross ){
    XLALDestroyREAL8Vector( resp->fcross );
  }
  if ( resp->fx ){
    XLALDestroyREAL8Vector( resp->fx );
  }
  if ( resp->fy ){
    XLALDestroyREAL8Vector( resp->fy );
  }
  if ( resp->fb ){
    XLALDestroyREAL8Vector( resp->fb );
  }
  if ( resp->fl ){
    XLALDestroyREAL8Vector( resp->fl );
  }

  XLALFree(resp);
}


/**
 * \brief Convert sources parameters into amplitude and phase notation parameters
 *
 * Convert the physical source parameters into the amplitude and phase notation
 * given in Eqns 62-65 of \cite Jones2015 .
 *
 * Note that \c phi0 is essentially the rotational phase of the pulsar. Also,
 * note that if using \f$h_0\f$, and therefore the convention for a signal as
 * defined in \cite JKS98 , the sign of the waveform model is the opposite of
 * that in \cite Jones2015 , and therefore a sign flip is required in the
 * amplitudes.
 */
void XLALPulsarSourceToWaveformParams( PulsarParameters *params ){
  REAL8 sinlambda, coslambda, sinlambda2, coslambda2, sin2lambda;
  REAL8 theta, sintheta, costheta2, sintheta2, sin2theta;
  REAL8 phi0 = PulsarGetREAL8ParamOrZero( params, "PHI0" );
  REAL8 h0 = PulsarGetREAL8ParamOrZero( params, "H0" );
  REAL8 I21 = PulsarGetREAL8ParamOrZero( params, "I21" );
  REAL8 I31 = PulsarGetREAL8ParamOrZero( params, "I31" );
  REAL8 C21 = PulsarGetREAL8ParamOrZero( params, "C21" );
  REAL8 C22 = PulsarGetREAL8ParamOrZero( params, "C22" );
  REAL8 phi21 = PulsarGetREAL8ParamOrZero( params, "PHI21" );
  REAL8 phi22 = PulsarGetREAL8ParamOrZero( params, "PHI22" );
  REAL8 lambda = PulsarGetREAL8ParamOrZero( params, "LAMBDA" );
  REAL8 costheta = PulsarGetREAL8ParamOrZero( params, "COSTHETA" );

  if ( h0 != 0.){
    phi22 = 2.*phi0;
    phi22 = phi22 - LAL_TWOPI*floor(phi22/LAL_TWOPI);
    PulsarAddParam( params, "PHI22", &phi22, PULSARTYPE_REAL8_t );

    C22 = -0.5*h0; /* note the change in sign so that the triaxial model conforms to the convertion in JKS98 */
    PulsarAddParam( params, "C22", &C22, PULSARTYPE_REAL8_t );
  }
  else if ( ( I21 != 0. || I31 != 0. ) && ( C22 == 0. && C21 == 0. ) ) {
    sinlambda = sin( lambda );
    coslambda = cos( lambda );
    sin2lambda = sin( 2. * lambda );
    sinlambda2 = SQUARE( sinlambda );
    coslambda2 = SQUARE( coslambda );

    theta = acos( costheta );
    sintheta = sin( theta );
    sin2theta = sin( 2. * theta );
    sintheta2 = SQUARE( sintheta );
    costheta2 = SQUARE( costheta );

    REAL8 A22 = I21 * ( sinlambda2 - coslambda2 * costheta2 ) - I31 * sintheta2;
    REAL8 B22 = I21 * sin2lambda * costheta;
    REAL8 A222 = SQUARE( A22 );
    REAL8 B222 = SQUARE( B22 );

    REAL8 A21 = I21 * sin2lambda * sintheta;
    REAL8 B21 = sin2theta * ( I21 * coslambda2 - I31 );
    REAL8 A212 = SQUARE( A21 );
    REAL8 B212 = SQUARE( B21 );

    C22 = 2.*sqrt( A222 + B222 );
    C21 = 2.*sqrt( A212 + B212 );

    PulsarAddParam( params, "C22", &C22, PULSARTYPE_REAL8_t );
    PulsarAddParam( params, "C21", &C21, PULSARTYPE_REAL8_t );

    phi22 = fmod( 2.*phi0 - atan2( B22, A22 ), LAL_TWOPI );
    phi21 = fmod( phi0 - atan2( B21, A21 ), LAL_TWOPI );

    PulsarAddParam( params, "PHI22", &phi22, PULSARTYPE_REAL8_t );
    PulsarAddParam( params, "PHI21", &phi21, PULSARTYPE_REAL8_t );
  }
}
