/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, Colin Gill, John Veitch
 *
 * \brief Pulsar model functions for use in parameter estimation codes for targeted pulsar searches.
 */

#include "ppe_models.h"
#include <lal/CWFastMath.h>

#define SQUARE(x) ( (x) * (x) );

static BinaryPulsarParams empty_BinaryPulsarParams;

/******************************************************************************/
/*                            MODEL FUNCTIONS                                 */
/******************************************************************************/

/**
 * \brief Defines the pulsar model/template to use
 *
 * This function is the wrapper for functions defining the pulsar model template to be used in the analysis. It also
 * uses \c rescale_parameter to scale any parameters back to their true values for use in the model and places them into
 * a \c BinaryPulsarParams structure.
 *
 * Note: Any additional models should be added into this function.
 *
 * \param data [in] The data structure hold data and current parameter info
 *
 * \sa rescale_parameter
 * \sa pulsar_model
 */
void get_pulsar_model( LALInferenceIFOData *data ){
  BinaryPulsarParams pars = empty_BinaryPulsarParams; /* initialise as empty */

  /* set model parameters (including rescaling) */
  //pars.h0 = rescale_parameter( data, "h0" );
  pars.cosiota = rescale_parameter( data, "cosiota" );

  /* check whether new psi-phi0 coordinates are used */
  //if ( LALInferenceCheckVariable( data->modelParams, "psiprime" ) &&
  //     LALInferenceCheckVariable( data->modelParams, "phi0prime" ) ){
  //  REAL8 phi0prime = rescale_parameter( data, "phi0prime" );
  //  REAL8 psiprime = rescale_parameter( data, "psiprime" );

  //  /* convert phi0' and psi' into phi0 and psi */
  //  inverse_phi0_psi_transform( phi0prime, psiprime, &pars.phi0, &pars.psi );
  //}
  //else{
  pars.psi = rescale_parameter( data, "psi" );
  //  pars.phi0 = rescale_parameter( data, "phi0" );
  //}

  if( LALInferenceCheckVariable( data->dataParams, "jones-model" ) ){
    /* use parameterisation from Ian Jones's original model */
    pars.I21 = rescale_parameter( data, "I21" );
    pars.I31 = rescale_parameter( data, "I31" );
    pars.lambda = rescale_parameter( data, "lambda" );
    pars.costheta = rescale_parameter( data, "costheta" );
    pars.phi0 = rescale_parameter( data, "phi0" );

    invert_source_params( &pars );
  }
  else{
    pars.C21 = rescale_parameter( data, "C21" );
    pars.C22 = rescale_parameter( data, "C22" );
    pars.phi21 = rescale_parameter( data, "phi21" );

    if( LALInferenceCheckVariable( data->dataParams, "biaxial" ) ){
      /* use complex amplitude parameterisation, but set up for a biaxial star */
      pars.phi22 = 2.*pars.phi21;
    }
    else{
      pars.phi22 = rescale_parameter( data, "phi22" );
    }
  }

  /* set the potentially variable parameters */
  pars.pepoch = rescale_parameter( data, "pepoch" );
  pars.posepoch = rescale_parameter( data, "posepoch" );

  pars.ra = rescale_parameter( data, "ra" );
  pars.pmra = rescale_parameter( data, "pmra" );
  pars.dec = rescale_parameter( data, "dec" );
  pars.pmdec = rescale_parameter( data, "pmdec" );

  pars.f0 = rescale_parameter( data, "f0" );
  pars.f1 = rescale_parameter( data, "f1" );
  pars.f2 = rescale_parameter( data, "f2" );
  pars.f3 = rescale_parameter( data, "f3" );
  pars.f4 = rescale_parameter( data, "f4" );
  pars.f5 = rescale_parameter( data, "f5" );

  /* speed of GWs as a fraction of speed of light LAL_C_SI */
  pars.cgw = rescale_parameter( data, "cgw" );

  /* check if there are binary parameters */
  if( LALInferenceCheckVariable(data->modelParams, "model") ){
    /* binary system model - NOT pulsar model */
    pars.model = *(CHAR**)LALInferenceGetVariable( data->modelParams, "model" );

    pars.e = rescale_parameter( data, "e" );
    pars.w0 = rescale_parameter( data, "w0" );
    pars.Pb = rescale_parameter( data, "Pb" );
    pars.x = rescale_parameter( data, "x" );
    pars.T0 = rescale_parameter( data, "T0" );

    pars.e2 = rescale_parameter( data, "e2" );
    pars.w02 = rescale_parameter( data, "w02" );
    pars.Pb2 = rescale_parameter( data, "Pb2" );
    pars.x2 = rescale_parameter( data, "x2" );
    pars.T02 = rescale_parameter( data, "T02" );

    pars.e3 = rescale_parameter( data, "e3" );
    pars.w03 = rescale_parameter( data, "w03" );
    pars.Pb3 = rescale_parameter( data, "Pb3" );
    pars.x3 = rescale_parameter( data, "x3" );
    pars.T03 = rescale_parameter( data, "T03" );

    pars.xpbdot = rescale_parameter( data, "xpbdot" );
    pars.eps1 = rescale_parameter( data, "eps1" );
    pars.eps2 = rescale_parameter( data, "eps2" );
    pars.eps1dot = rescale_parameter( data, "eps1dot" );
    pars.eps2dot = rescale_parameter( data, "eps2dot" );
    pars.Tasc = rescale_parameter( data, "Tasc" );

    pars.wdot = rescale_parameter( data, "wdot" );
    pars.gamma = rescale_parameter( data, "gamma" );
    pars.Pbdot = rescale_parameter( data, "Pbdot" );
    pars.xdot = rescale_parameter( data, "xdot" );
    pars.edot = rescale_parameter( data, "edot" );

    pars.s = rescale_parameter( data, "s" );
    pars.dr = rescale_parameter( data, "dr" );
    pars.dth = rescale_parameter( data, "dth" );
    pars.a0 = rescale_parameter( data, "a0" );
    pars.b0 = rescale_parameter( data, "b0" );

    pars.M = rescale_parameter( data, "M" );
    pars.m2 = rescale_parameter( data, "m2" );
  }

  /* now get pulsar model */
  pulsar_model( pars, data );

}


/**
 * \brief Rescale parameter back to its true value
 *
 * This function will rescale a parameter to its true value using the scale factor and minimum scale value.
 *
 * \param data [in] data structure containing parameter information
 * \param parname [in] name of the parameter requiring rescaling
 *
 * \return Rescaled parameter value
 */
REAL8 rescale_parameter( LALInferenceIFOData *data, const CHAR *parname ){
  REAL8 par = 0., scale = 0., offset = 0.;
  CHAR scaleName[VARNAME_MAX] = "";
  CHAR offsetName[VARNAME_MAX] = "";

  sprintf(scaleName, "%s_scale", parname);
  sprintf(offsetName, "%s_scale_min", parname);

  scale = *(REAL8*)LALInferenceGetVariable( data->dataParams, scaleName );
  offset = *(REAL8*)LALInferenceGetVariable( data->dataParams, offsetName );

  par = *(REAL8*)LALInferenceGetVariable( data->modelParams, parname );

  par = par*scale + offset;

  return par;
}


/**
 * \brief Generate the model of the neutron star signal
 *
 * The function requires that the pulsar model is set using the \c model-type command line argument (this is set in \c
 * main, and if not specified defaults to a \c triaxial model). Currently the model can be \c triaxial for quadrupole
 * emission from a triaxial star at twice the rotation freqeuncy, or \c pinsf for a two component emission model with
 * emission at the rotation frequency <i>and</i> twice the rotation frequency. Depending on the specified model the
 * function calls the appropriate model function.
 *
 * Firstly the time varying amplitude of the signal will be calculated based on the antenna pattern and amplitude
 * parameters. Then, if searching over phase parameters, the phase evolution of the signal will be calculated. The
 * difference between the new phase model, \f$\phi(t)_n\f$, and that used to heterodyne the data, \f$\phi(t)_h\f$,
 * (stored in \c data->timeData->data) will be calculated and the complex signal model, \f$M\f$, modified accordingly:
 * \f[
 * M'(t) = M(t)\exp{i(-(\phi(t)_n - \phi(t)_h))}.
 * \f]
 *
 * \param params [in] A \c BinaryPulsarParams structure containing the model parameters
 * \param data [in] The data structure containing the detector data and additional info
 *
 * \sa get_amplitude_model
 * \sa get_amplitude_model
 * \sa get_phase_model
 */
void pulsar_model( BinaryPulsarParams params, LALInferenceIFOData *data ){
  INT4 i = 0, length = 0;
  UINT4 j = 0;
  REAL8 mm = 0.;

  /* check model type to get amplitude model */
  // modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );

  get_amplitude_model( params, data );

  /* if ( !strcmp( modeltype, "triaxial" ) ){
    get_triaxial_amplitude_model( params, data );
  }
  else if ( !strcmp( modeltype, "pinsf" ) ){
    get_pinsf_amplitude_model( params, data );
  } */
  /* ADD NEW MODELS HERE */
  /* else{
    fprintf(stderr, "Error... model '%s' is not defined!\n", modeltype);
    exit(0);
  } */

  /* get difference in phase for f component and perform extra heterodyne */
  REAL8Vector *freqFactors = NULL;
  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  if( LALInferenceCheckVariable( data->dataParams, "mismatch" ) ){
    mm = *(REAL8 *)LALInferenceGetVariable( data->dataParams, "mismatch" );
  }

  for( j = 0; j < freqFactors->length; j++ ){
    REAL8Vector *dphi = NULL;
    UINT4 nohet = 0; /* set if extra phase heterodyne is not required */

    /* move data pointer along one as one iteration of the model is held over j data structures, moved this from bottom
     * of loop so actioned for 2nd run through the loop.*/
    if ( j > 0 ) { data = data->next; }

    length = data->compModelData->data->length;
    /* the timeData vector within the LALIFOData structure contains the phase calculated using the initial (heterodyne)
     * values of the phase parameters */
    if ( varyphase ){
      /* check whether to recompute the full phase or not */
      if( LALInferenceCheckVariable( data->dataParams, "downsampled_times" ) ){
        REAL8Vector *dsdphi1 = NULL, *dsdphi2 = NULL;
        LIGOTimeGPSVector *downst =
          *(LIGOTimeGPSVector **)LALInferenceGetVariable( data->dataParams, "downsampled_times" );

        /* get the previous downsampled phase if it exists */
        if ( LALInferenceCheckVariable( data->dataParams, "ds_phase" ) ){
          dsdphi1 = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "ds_phase" );
        }
        else{
          XLALPrintError("Error, downsampled phase does not exist\n");
          XLAL_ERROR_VOID(XLAL_EFAILED);
        }

        /* get the downsampled phase for the current parameters */
        dsdphi2 = get_phase_model( params, data, freqFactors->data[j], 1 );

        /* work out phase mismatch (if any value in dsdphi1 is not zero it means
           ds_phase has been set) */
        if( dsdphi1->data[dsdphi1->length-1] != 0. && dsdphi2 ){
          REAL8 mmcalc = get_phase_mismatch( dsdphi1, dsdphi2, downst );

          /* if small mismatch then just use previous phase if available */
          if ( mmcalc < mm ) { nohet = 1; }
        }

        /* make sure the "previous" down sampled phase is the right one for comparison */
        if ( !nohet ) { memcpy(dsdphi1->data, dsdphi2->data, sizeof(REAL8)*dsdphi1->length ); }

        XLALDestroyREAL8Vector( dsdphi2 );
      }

      /* reheterodyne with the phase */
      if ( !nohet ){
        if ( (dphi = get_phase_model( params, data, freqFactors->data[j], 0 )) != NULL ){
          for( i=0; i<length; i++ ){
            COMPLEX16 M;
            REAL8 dphit;
            COMPLEX16 expp;

            dphit = -fmod(dphi->data[i] - data->timeData->data->data[i], 1.);

            expp = cexp( LAL_TWOPI * I * dphit );

            M = data->compModelData->data->data[i];

            /* heterodyne */
            data->compModelData->data->data[i] = M * expp;
          }

          XLALDestroyREAL8Vector( dphi );
        }
      }
    }
  }
}


/**
 * \brief The phase evolution of a source
 *
 * This function will calculate the phase evolution of a source at a particular sky location as observed at Earth. The
 * phase evolution is described by a Taylor expansion:
 * \f[
 * \phi(T) = \sum_{k=1}^n \frac{f^{(k-1)}}{k!} T^k,
 * \f]
 * where \f$f^x\f$ is the xth time derivative of the gravitational wave frequency, and \f$T\f$ is the pulsar proper
 * time. Frequency time derivatives are currently allowed up to the fifth derivative. The pulsar proper time is
 * calculated by correcting the time of arrival at Earth, \f$t\f$ to the solar system barycentre and if necessary the
 * binary system barycenter, so \f$T = t + \delta{}t_{\rm SSB} + \delta{}t_{\rm BSB}\f$.
 *
 * In this function the time delay caused needed to correct to the solar system barycenter is only calculated if
 * required i.e. if it's not been previously calculated and an update is required due to a change in the sky position.
 * The same is true for the binary system time delay, which is only calculated if it has not previously been obtained or
 * needs updating due to a change in the binary system parameters.
 *
 * The solar system barycentre delay does not have to be explicitly computed for every time stamp passed to it, but
 * instead will just use linear interpolation within a time range set by \c interptime.
 *
 * \param params [in] A set of pulsar parameters
 * \param data [in] The data structure containing the detector data and additional info
 * \param freqFactor [in] the multiplicative factor on the pulsar frequency for a particular model
 * \param downsampled *UNDOCUMENTED*
 *
 * \return A vector of rotational phase values
 *
 * \sa get_ssb_delay
 * \sa get_bsb_delay
 */
REAL8Vector *get_phase_model( BinaryPulsarParams params, LALInferenceIFOData *data, REAL8 freqFactor,
                              UINT4 downsampled ){
  INT4 i = 0, length = 0;

  REAL8 T0 = 0., DT = 0., deltat = 0., deltat2 = 0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */

  REAL8Vector *phis = NULL, *dts = NULL, *bdts = NULL;
  LIGOTimeGPSVector *datatimes = NULL;

  /* check if we want to calculate the phase at a the downsampled rate */
  if ( downsampled ){
    if( LALInferenceCheckVariable( data->dataParams, "downsampled_times" ) ){
      datatimes = *(LIGOTimeGPSVector **)LALInferenceGetVariable(
        data->dataParams, "downsampled_times" );
    }
    else{
      fprintf(stderr, "Error, no downsampled time series available\n");
      exit(1);
    }
  }
  else datatimes = data->dataTimes;

  /* if edat is NULL then return a NULL pointer */
  if( data->ephem == NULL ) return NULL;

  length = datatimes->length;

  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );

  /* get time delays */
  /*Why ==NULL, surely it will equal null if not set to get ssb delays?*/
  if( (dts = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "ssb_delays" )) == NULL || varyskypos == 1 ){
    /* get time delays with an interpolation of interptime (30 mins) */
    dts = get_ssb_delay( params, datatimes, data->ephem, data->tdat, data->ttype, data->detector, interptime );
  }

  if( (bdts = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "bsb_delays" )) == NULL || varybinary == 1 ){
    /* get binary system time delays */
    bdts = get_bsb_delay( params, datatimes, dts, data->ephem );
  }

  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &data->dataTimes->data[i] ); /*time of data*/

    T0 = params.pepoch; /*time of ephem info*/

    DT = realT - T0; /*time diff between data and ephem info*/

    if ( params.model != NULL ) { deltat = DT + dts->data[i] + bdts->data[i]; }
    else { deltat = DT + dts->data[i]; }

    deltat /= (1.-params.cgw); /* correct for speed of GW compared to speed of light */

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = freqFactor*deltat*(params.f0 +
      inv_fact[2]*params.f1*deltat +
      inv_fact[3]*params.f2*deltat2 +
      inv_fact[4]*params.f3*deltat*deltat2 +
      inv_fact[5]*params.f4*deltat2*deltat2 +
      inv_fact[6]*params.f5*deltat2*deltat2*deltat);
  }

  /* free memory */
  if ( !LALInferenceCheckVariable( data->dataParams, "ssb_delays") || varyskypos == 1 )
    XLALDestroyREAL8Vector( dts );

  if ( !LALInferenceCheckVariable( data->dataParams, "bsb_delays") || varybinary == 1 )
    XLALDestroyREAL8Vector( bdts );

  return phis;
}


/**
 * \brief Computes the delay between a GPS time at Earth and the solar system barycentre
 *
 * This function calculate the time delay between a GPS time at a specific location (e.g. a gravitational wave detector)
 * on Earth and the solar system barycentre. The delay consists of three components: the geometric time delay (Roemer
 * delay) \f$t_R = \mathbf{r}(t)\hat{n}/c\f$ (where \f$\mathbf{r}(t)\f$ is the detector's position vector at time
 * \f$t\f$), the special relativistic Einstein delay \f$t_E\f$, and the general relativistic Shapiro delay \f$t_S\f$.
 *
 * Rather than computing the time delay at every time stamp passed to the function it is instead (if requested) able to
 * perform linear interpolation to a point within a range given by \c interptime.
 *
 * \param pars [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times at Earth
 * \param ephem [in] Information on the solar system ephemeris
 * \param detector [in] Information on the detector position on the Earth
 * \param interptime [in] The time (in seconds) between explicit recalculations of the time delay
 * \param tdat *UNDOCUMENTED*
 * \param ttype *UNDOCUMENTED*
 * \return A vector of time delays in seconds
 *
 * \sa XLALBarycenter
 * \sa XLALBarycenterEarth
 */
REAL8Vector *get_ssb_delay( BinaryPulsarParams pars, LIGOTimeGPSVector *datatimes, EphemerisData *ephem,
                            TimeCorrectionData *tdat, TimeCorrectionType ttype, LALDetector *detector,
                            REAL8 interptime ){
  INT4 i = 0, length = 0;

  REAL8 T0 = 0., DT = 0., DTplus = 0.;

  EarthState earth, earth2;
  EmissionTime emit, emit2;

  BarycenterInput *bary = NULL;

  REAL8Vector *dts = NULL;

  /* if edat is NULL then return a NULL poniter */
  if( ephem == NULL ) { return NULL; }

  /* copy barycenter and ephemeris data */
  bary = (BarycenterInput*)XLALCalloc( 1, sizeof(BarycenterInput) );
  memcpy( &bary->site, detector, sizeof(LALDetector) );

  bary->alpha = pars.ra;
  bary->delta = pars.dec;

   /* set the position and frequency epochs if not already set */
  if( pars.pepoch == 0. && pars.posepoch != 0.) { pars.pepoch = pars.posepoch; }
  else if( pars.posepoch == 0. && pars.pepoch != 0. ) { pars.posepoch = pars.pepoch; }

  length = datatimes->length;

  /* allocate memory for times delays */
  dts = XLALCreateREAL8Vector( length );

  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( pars.px != 0. ) { bary->dInv = pars.px*1e-3*LAL_C_SI/LAL_PC_SI; }
  else if( pars.dist != 0. ) { bary->dInv = LAL_C_SI/(pars.dist*1e3*LAL_PC_SI); }
  else { bary->dInv = 0.; }

  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &datatimes->data[i] );

    T0 = pars.pepoch;

    DT = realT - T0;

    /* only do call to the barycentring routines once every interptime (unless
       interptime == 0), otherwise just linearly interpolate between them */
    if( i == 0 || DT > DTplus || interptime == 0 ){
      bary->tgps = datatimes->data[i];

      bary->delta = pars.dec + ( realT - pars.posepoch ) * pars.pmdec;
      bary->alpha = pars.ra + ( realT - pars.posepoch ) * pars.pmra / cos( bary->delta );

      /* call barycentring routines */
      XLAL_CHECK_NULL( XLALBarycenterEarthNew( &earth, &bary->tgps, ephem, tdat, ttype ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( XLALBarycenter( &emit, bary, &earth ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* add interptime to the time */
      if ( interptime > 0 ){
        DTplus = DT + interptime;
        XLALGPSAdd( &bary->tgps, interptime );

        /* No point in updating the positions as difference will be tiny */
        XLAL_CHECK_NULL( XLALBarycenterEarthNew( &earth2, &bary->tgps, ephem,tdat, ttype )
          == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK_NULL( XLALBarycenter( &emit2, bary, &earth2 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    /* linearly interpolate to get emitdt */
    if( interptime > 0. ){
      dts->data[i] = emit.deltaT + (DT - (DTplus - interptime)) *
        (emit2.deltaT - emit.deltaT)/interptime;
    }
    else { dts->data[i] = emit.deltaT; }
  }

  XLALFree( bary );

  return dts;
}


/**
 * \brief Computes the delay between a pulsar in a binary system and the barycentre of the system
 *
 * This function uses \c XLALBinaryPulsarDeltaT to calculate the time delay between for a pulsar in a binary system
 * between the time at the pulsar and the time at the barycentre of the system. This includes Roemer delays and
 * relativistic delays. The orbit may be described by different models and can be purely Keplarian or include various
 * relativistic corrections.
 *
 * \param pars [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times
 * \param dts [in] A vector of solar system barycentre time delays
 * \param edat *UNDOCUMENTED*
 * \return A vector of time delays in seconds
 *
 * \sa XLALBinaryPulsarDeltaT
 */
REAL8Vector *get_bsb_delay( BinaryPulsarParams pars, LIGOTimeGPSVector *datatimes, REAL8Vector *dts,
                            EphemerisData *edat ){
  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;
  REAL8Vector *bdts = NULL;

  INT4 i = 0, length = datatimes->length;

  bdts = XLALCreateREAL8Vector( length );

  for ( i = 0; i < length; i++ ){
    /* check whether there's a binary model */
    if ( pars.model ){
      EarthState earth;

      binput.tb = XLALGPSGetREAL8( &datatimes->data[i] ) + dts->data[i];

      get_earth_pos_vel( &earth, edat, &datatimes->data[i] );

      binput.earth = earth; /* current Earth state */
      XLALBinaryPulsarDeltaT( &boutput, &binput, &pars );
      bdts->data[i] = boutput.deltaT;
    }
    else { bdts->data[i] = 0.; }
  }

  return bdts;
}


/**
 * \brief The amplitude model of a complex heterodyned triaxial neutron star
 *
 * This function calculates the complex heterodyned time series model for a triaxial neutron star (see
 * [\cite DupuisWoan2005]). It is defined as:
 * \f{eqnarray}{
 * y(t) & = & \frac{h_0}{2} \left( \frac{1}{2}F_+(t,\psi)
 * (1+\cos^2\iota)\exp{i\phi_0} - iF_{\times}(t,\psi)\cos{\iota}\exp{i\phi_0}
 * \right),
 * \f}
 * where \f$F_+\f$ and \f$F_{\times}\f$ are the antenna response functions for the plus and cross polarisations.
 *
 * The antenna pattern functions are contained in a 2D lookup table, so within this function the correct value for the
 * given time and \f$\psi\f$ are interpolated from this lookup table using bilinear interpolation (e.g.):
 * \f{eqnarray}{
 * F_+(\psi, t) = F_+(\psi_i, t_j)(1-\psi)(1-t) + F_+(\psi_{i+1}, t_j)\psi(1-t)
 * + F_+(\psi_i, t_{j+1})(1-\psi)t + F_+(\psi_{i+1}, t_{j+1})\psi{}t,
 * \f}
 * where \f$\psi\f$ and \f$t\f$ have been scaled to be within a unit square, and \f$\psi_i\f$ and \f$t_j\f$ are the
 * closest points within the lookup table to the required values.
 *
 * \param pars [in] A set of pulsar parameters
 * \param data [in] The data parameters giving information on the data and detector
 *
 * DEPRECATED: use \c get_amplitude_model instead.
 */
void get_triaxial_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data ){
  INT4 i = 0, length;

  REAL8 psteps, tsteps, psv, tsv;
  INT4 psibinMin, psibinMax, timebinMin, timebinMax;
  REAL8 plus, cross;
  REAL8 plus00, plus01, plus10, plus11, cross00, cross01, cross10, cross11;
  REAL8 psiScaled, timeScaled;
  REAL8 psiMin, psiMax, timeMin, timeMax;
  REAL8 T;
  REAL8 Xplus, Xcross;
  COMPLEX16 expiphi, Xpexpphi, Xcexpphi;

  gsl_matrix *LU_Fplus, *LU_Fcross;
  REAL8Vector *sidDayFrac = NULL;

  length = data->dataTimes->length;

  /* set lookup table parameters */
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );

  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fplus" );
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fcross" );

  /* get the sidereal time since the initial data point % sidereal day */
  sidDayFrac = *(REAL8Vector**)LALInferenceGetVariable( data->dataParams, "siderealDay" );

  expiphi = cexp( I * pars.phi0 );

  /************************* CREATE MODEL *************************************/
  /* This model is a complex heterodyned time series for a triaxial neutron star emitting at twice its rotation
   * frequency (as defined in Dupuis and Woan, PRD, 2005):
   *    h(t) = (h0/2) * ((1/2)*F+(t)*(1+cos(iota)^2)*exp(i*phi0) - i*Fx(t)*cos(iota)*exp(i*phi0))
   ****************************************************************************/
  Xplus = 0.25*(1.+pars.cosiota*pars.cosiota)*pars.h0;
  Xcross = 0.5*pars.cosiota*pars.h0;
  Xpexpphi = Xplus*expiphi;
  Xcexpphi = Xcross*expiphi;

  /* set the psi bin for the lookup table - the lookup table runs from -pi/2 to pi/2, but for the triaxial case we only
   * require psi values from -pi/4 to pi/4 (the grid will be twice as coarse) */
  psv = LAL_PI / ( psteps - 1. );
  psibinMin = (INT4)floor( ( pars.psi + LAL_PI_2 )/psv );
  psiMin = -(LAL_PI_2) + psibinMin*psv;
  psibinMax = psibinMin + 1;
  psiMax = psiMin + psv;

  /* rescale psi for bilinear interpolation on a unit square */
  psiScaled = (pars.psi - psiMin)/(psiMax - psiMin);

  tsv = LAL_DAYSID_SI / tsteps;

  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = sidDayFrac->data[i];
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
   timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;

    /* get values of matrix for bilinear interpolation */
    plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
    plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
    plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
    plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );

    cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
    cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
    cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
    cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );

    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);

    plus = plus00*(1. - psiScaled)*(1. - timeScaled) +
      plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
      plus11*psiScaled*timeScaled;
    cross = cross00*(1. - psiScaled)*(1. - timeScaled) +
      cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
      + cross11*psiScaled*timeScaled;

    /* create the complex signal amplitude model */
    data->compModelData->data->data[i] = plus*Xpexpphi - I*cross*Xcexpphi;
  }
}


/**
 * \brief The amplitude model of a complex heterodyned signal from a NS rotating about the pinning axis of its pinned
 * superfluid component.
 *
 * This function calculates the complex heterodyned time series model for a triaxial neutron star rotating about the
 * pinning axis of its pinned superfluid component.
 *
 * Unlike the standard triaxial model, this model has emission at f and 2f, therefore this model function processes two
 * sets of data per detector. In this model the \f$\phi_0\f$ parameter is the initial rotational phase, rather than the
 * GW phase as in the triaxial model.
 *
 * As for the standard triaxial model, the antenna pattern functions are contained in a 2D lookup table, so within this
 * function the correct value for the given time and \f$\psi\f$ are interpolated from this lookup table using bilinear
 * interpolation (e.g.):
 * \f{eqnarray}{
 * F_+(\psi, t) = F_+(\psi_i, t_j)(1-\psi)(1-t) + F_+(\psi_{i+1}, t_j)\psi(1-t)
 * + F_+(\psi_i, t_{j+1})(1-\psi)t + F_+(\psi_{i+1}, t_{j+1})\psi{}t,
 * \f}
 * where \f$\psi\f$ and \f$t\f$ have been scaled to be within a unit square, and \f$\psi_i\f$ and \f$t_j\f$ are the
 * closest points within the lookup table to the required values.
 *
 * \param pars [in] A set of pulsar parameters
 * \param data [in] The data parameters giving information on the data and detector
 *
 * DEPRECATED: use \c get_amplitude_model instead.
 */
void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data ){
  INT4 i = 0, length;

  REAL8 psteps, tsteps, psv, tsv;
  INT4 psibinMin, psibinMax, timebinMin, timebinMax;
  REAL8 plus00, plus01, plus10, plus11, cross00, cross01, cross10, cross11;
  REAL8 psiScaled, timeScaled;
  REAL8 psiMin, psiMax, timeMin, timeMax;
  REAL8 plus, cross;
  REAL8 T;
  REAL8 Xplusf, Xcrossf, Xplus2f, Xcross2f;
  REAL8 A21, A22, B21, B22;
  COMPLEX16 ePhi, e2Phi;
  REAL8 iota = acos(pars.cosiota), theta = acos(pars.costheta);
  REAL8 siniota = sin(iota);
  REAL8 sintheta = sin(theta), sin2theta = sin( 2.*theta );
  REAL4 coslambda, sinlambda;
  REAL8 sin2lambda = sin( 2.*pars.lambda );
  REAL8 f2_r;

  gsl_matrix *LU_Fplus, *LU_Fcross;
  REAL8Vector *sidDayFrac1 = NULL;
  REAL8Vector *sidDayFrac2 = NULL;

  /* set lookup table parameters */
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );

  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fplus" );
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fcross");
  /* get the sidereal time since the initial data point % sidereal day */
  sidDayFrac1 = *(REAL8Vector**)LALInferenceGetVariable( data->dataParams, "siderealDay" );

  /* phi0 here is rotational phase not GW phase */
  ePhi = cexp( pars.phi0 * I );
  e2Phi = cexp( 2. * pars.phi0 * I );

  XLAL_CHECK_VOID( XLALSinCosLUT( &sinlambda, &coslambda, pars.lambda ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* f^2 / r */
  f2_r = pars.f0 * pars.f0 / pars.dist;

  /************************* CREATE MODEL *************************************
   * This model is a complex heterodyned time series for a pinned superfluid neutron star emitting at its roation
   * frequency and twice its rotation frequency (as originally defined in Jones 2009, but using Eqns 35-38 of Jones
   * 2012 LIGO DCC T1200265-v3).
   ****************************************************************************/
  Xplusf = -( f2_r / 2. ) * siniota * pars.cosiota;
  Xcrossf = ( f2_r / 2. ) * siniota;
  Xplus2f = -f2_r * ( 1. + pars.cosiota * pars.cosiota );
  Xcross2f = 2. * f2_r * pars.cosiota;

  A21 = pars.I21 * sin2lambda * sintheta;
  B21 = ( pars.I21 * coslambda * coslambda - pars.I31 ) * sin2theta;

  A22 = pars.I21 * ( sinlambda * sinlambda - coslambda * coslambda * pars.costheta * pars.costheta ) -
    pars.I31 * sintheta * sintheta;
  B22 = pars.I21 * sin2lambda * pars.costheta;

  /* set the psi bin for the lookup table (look-up table cover that fill -pi/2 to pi/2 range) */
  psv = LAL_PI / ( psteps - 1. );
  psibinMin = (INT4)floor( ( pars.psi + LAL_PI_2 )/psv );
  psiMin = -(LAL_PI_2) + psibinMin*psv;
  psibinMax = psibinMin + 1;
  psiMax = psiMin + psv;

  /* rescale psi for bilinear interpolation on a unit square */
  psiScaled = (pars.psi - psiMin)/(psiMax - psiMin);

  tsv = LAL_DAYSID_SI / tsteps;

  /*--------------------------------------------------------------------------*/
  /* set model for 1f component */
  length = data->dataTimes->length;

  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = sidDayFrac1->data[i];
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;

    /* get values of matrix for bilinear interpolation */
    plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
    plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
    plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
    plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );

    cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
    cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
    cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
    cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );

    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);

    plus = plus00*(1. - psiScaled)*(1. - timeScaled) +
      plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
      plus11*psiScaled*timeScaled;
    cross = cross00*(1. - psiScaled)*(1. - timeScaled) +
      cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
      + cross11*psiScaled*timeScaled;

    /* create the complex signal amplitude model at f */
    data->compModelData->data->data[i] = ( plus * Xplusf * ePhi * ( A21 - I * B21 ) ) +
                                         ( cross * Xcrossf * ePhi * ( B21 + I * A21 ) );
  }

  /*--------------------------------------------------------------------------*/
  /* set model for 2f component */
  length = data->next->dataTimes->length;

  sidDayFrac2 = *(REAL8Vector**)LALInferenceGetVariable( data->next->dataParams, "siderealDay" );

  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = sidDayFrac2->data[i];
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;

    /* get values of matrix for bilinear interpolation */
    plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
    plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
    plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
    plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );

    cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
    cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
    cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
    cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );

    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);

    plus = plus00*(1. - psiScaled)*(1. - timeScaled) +
      plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
      plus11*psiScaled*timeScaled;
    cross = cross00*(1. - psiScaled)*(1. - timeScaled) +
      cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
      + cross11*psiScaled*timeScaled;

    /* create the complex signal amplitude model at 2f */
    data->next->compModelData->data->data[i] = ( plus * Xplus2f * e2Phi * ( A22 - I * B22 ) ) +
                                               ( cross * Xcross2f * e2Phi * ( B22 + I * A22 ) );
  }
  /*--------------------------------------------------------------------------*/
}


/**
 * \brief The amplitude model of a complex heterodyned signal from the \f$l=2, m=1,2\f$ harmonics of a rotating neutron
 * star.
 *
 * This function calculates the complex heterodyned time series model for a rotating neutron star. It will currently
 * calculate the model for emission from the \f$l=m=2\f$ harmonic (which gives emission at twice the rotation frequency)
 * and/or the \f$l=2\f$ and \f$m=1\f$ harmonic (which gives emission at the rotation frequency). See LIGO T1200265-v3.
 * Further harmonics can be added and are defined by the \c freqFactor value, which is the multiple of the
 * spin-frequency at which emission is produced.
 *
 * The antenna pattern functions are contained in a 2D lookup table, so within this function the correct value for the
 * given time and \f$\psi\f$ is interpolated from this lookup table using bilinear interpolation (e.g.):
 * \f{eqnarray}{
 * F_+(\psi, t) = F_+(\psi_i, t_j)(1-\psi)(1-t) + F_+(\psi_{i+1}, t_j)\psi(1-t)
 * + F_+(\psi_i, t_{j+1})(1-\psi)t + F_+(\psi_{i+1}, t_{j+1})\psi{}t,
 * \f}
 * where \f$\psi\f$ and \f$t\f$ have been scaled to be within a unit square, and \f$\psi_i\f$ and \f$t_j\f$ are the
 * closest points within the lookup table to the required values.
 *
 * \param pars [in] A set of pulsar parameters
 * \param data [in] The data parameters giving information on the data and detector
 *
 */
void get_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data ){
  UINT4 i = 0, j = 0, length;

  REAL8 psteps, tsteps, psv, tsv;
  INT4 psibinMin, psibinMax, timebinMin, timebinMax;
  REAL8 plus00, plus01, plus10, plus11, cross00, cross01, cross10, cross11;
  REAL8 psiScaled, timeScaled;
  REAL8 psiMin, psiMax, timeMin, timeMax;
  REAL8 plus, cross;
  REAL8 T;
  REAL8 siniota = sin(acos(pars.cosiota));

  gsl_matrix *LU_Fplus, *LU_Fcross;
  REAL8Vector *freqFactors = NULL;

  /* set lookup table parameters */
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );

  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fplus" );
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fcross" );

  freqFactors = *(REAL8Vector**)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  psv = LAL_PI / ( psteps - 1. );
  psibinMin = (INT4)floor( ( pars.psi + LAL_PI_2 )/psv );
  psiMin = -(LAL_PI_2) + psibinMin*psv;
  psibinMax = psibinMin + 1;
  psiMax = psiMin + psv;

  /* rescale psi for bilinear interpolation on a unit square */
  psiScaled = (pars.psi - psiMin)/(psiMax - psiMin);

  tsv = LAL_DAYSID_SI / tsteps;

  /* loop over components in data as given by the frequency factors */
  for( j = 0; j < freqFactors->length; j++ ){
    REAL8Vector *sidDayFrac = NULL;
    COMPLEX16 expPhi;
    COMPLEX16 Cplus, Ccross;

    if ( !data ){
      XLALPrintError( "%s: Error... data not defined.\n", __func__ );
      XLAL_ERROR_VOID( XLAL_EINVAL );
    }

    /* get the sidereal time since the initial data point % sidereal day */
    sidDayFrac = *(REAL8Vector**)LALInferenceGetVariable( data->dataParams, "siderealDay" );

    length = data->dataTimes->length;

    /* get the amplitude and phase factors */
    if( freqFactors->data[j] == 1. ){
      /* the l=2, m=1 harmonic at the rotation frequency */
      expPhi = cexp( I * pars.phi21 );
      Cplus = -0.25 * pars.C21 * siniota * pars.cosiota * expPhi;
      Ccross = 0.25 * I * pars.C21 * siniota * expPhi;
    }
    else if( freqFactors->data[j] == 2. ){
      /* the l=2, m=2 harmonic at twice the rotation frequency */
      expPhi = cexp( I * pars.phi22 );
      Cplus = -0.5 * pars.C22 * ( 1. + pars.cosiota * pars.cosiota ) * expPhi;
      Ccross = I * pars.C22 * pars.cosiota * expPhi;
    }
    else{
      fprintf(stderr, "%s: Error... currently unknown frequency factor (%.2lf) for models.\n", __func__,
        freqFactors->data[j] );
      exit(1);
    }

    for( i=0; i<length; i++ ){
      /* set the time bin for the lookup table */
      /* sidereal day in secs*/
      T = sidDayFrac->data[i];
      timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
      timeMin = timebinMin*tsv;
      timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
      timeMax = timeMin + tsv;

      /* get values of matrix for bilinear interpolation */
      plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
      plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
      plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
      plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );

      cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
      cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
      cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
      cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );

      /* rescale time for bilinear interpolation on a unit square */
      timeScaled = (T - timeMin)/(timeMax - timeMin);

      plus = plus00*(1. - psiScaled)*(1. - timeScaled) +
        plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
        plus11*psiScaled*timeScaled;
      cross = cross00*(1. - psiScaled)*(1. - timeScaled) +
        cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
        + cross11*psiScaled*timeScaled;

      /* create the complex signal amplitude model appropriate for the harmonic */
      data->compModelData->data->data[i] = ( Cplus * plus ) + ( Ccross * cross );
    }

    data = data->next;
  }
}


/**
 * \brief Calculate the natural logarithm of the evidence that the data consists of only Gaussian noise
 * The function will calculate the natural logarithm of the evidence that the data (from one or more detectors) consists
 * of stationary segments/chunks described by a Gaussian with zero mean and unknown variance.
 *
 * The evidence is obtained from the joint likelihood given in \c pulsar_log_likelihood with the model term \f$y\f$ set
 * to zero.
 *
 * \param runState [in] The algorithm run state
 *
 * \return The natural logarithm of the noise only evidence
 */
REAL8 noise_only_model( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;

  REAL8 logL = 0.0;
  UINT4 i = 0;
  INT4 k = 0;

  REAL8Vector *freqFactors = NULL;
  FILE *fp = NULL;
  CHAR *Znoisefile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  /*-----------------------------*/
  /*get the outfile name*/
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );

  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  /* open the file to output noise evidence (or null signal evidence) for each individual data stream */
  /* set the Znoise filename to the outfile name with "_Znoise" appended */
  Znoisefile = XLALStringDuplicate( ppt->value );
  Znoisefile = XLALStringAppend( Znoisefile, "_Znoise" );

  if( (fp = fopen(Znoisefile, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open output Znoise file!\n");
    exit(0);
  }

  /*calculate the evidence */
  while ( data ){
    UINT4Vector *chunkLengths = NULL;
    REAL8Vector *sumDat = NULL;

    REAL8 chunkLength = 0.;
    REAL8 logLtmp = 0.;

    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams,  "chunkLength" );
    sumDat = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "sumData" );
    /*Sum the logL over the datachunks*/
    for (i=0; i<chunkLengths->length; i++){
      chunkLength = (REAL8)chunkLengths->data[i];

      logLtmp -= chunkLength * log(sumDat->data[i]) + LAL_LN2 * (chunkLength-1.) + gsl_sf_lnfact(chunkLength);
    }

    logL += logLtmp;

    /* output the noise evidence for each data stream for each detector  */
    fprintf(fp, "%s\t%.3lf\t%.16le\n", data->name, freqFactors->data[k], logLtmp);

    k += 1; /* advance counter now, as freqfactors array index starts at zero.*/

    /* reset k, freqfactor counter once all datastreamns for a detector are done */
    if( k >= (INT4)freqFactors->length ) { k = 0; }

    data = data->next;
  }

  fclose(fp);

  return logL;
}


/**
 * \brief Calculate the phase mismatch between two vectors of phases
 *
 * The function will calculate phase mismatch between two vectors of phases (with phases given in cycles rather than
 * radians).
 *
 * The mismatch is calculated as:
 * \f[
 * M = 1-\frac{1}{T}\int_0^T \cos{2\pi(\phi_1 - \phi_2)} dt.
 * \f]
 * In the function the integral is performed using the trapezium rule.
 *
 * PARAM phi1 [in] First phase vector
 * PARAM phi2 [in] Second phase vector
 * PARAM t [in] The time stamps of the phase points
 *
 * \return The mismatch
 */
REAL8 get_phase_mismatch( REAL8Vector *phi1, REAL8Vector *phi2, LIGOTimeGPSVector *t ){
  REAL8 mismatch = 0., dp1 = 0., dp2 = 0.;
  REAL4 sp, cp1, cp2;
  UINT4 i = 0;

  REAL8 T = 0., dt = 0.;

  /* data time span */
  T = XLALGPSGetREAL8(&t->data[t->length-1]) - XLALGPSGetREAL8(&t->data[0]);

  if ( phi1->length != phi2->length ){
    XLALPrintError("Phase lengths should be equal!\n");
    XLAL_ERROR_REAL8(XLAL_EFAILED);
  }

  /* calculate mismatch - integrate with trapezium rule */
  for( i = 0; i < phi1->length-1; i++ ){
    if ( i == 0 ){
      dp1 = fmod( phi1->data[i] - phi2->data[i], 1. );
      XLAL_CHECK_REAL8( XLALSinCos2PiLUT( &sp, &cp1, dp1 ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    else{
      dp1 = dp2;
      cp1 = cp2;
    }

    dp2 = fmod( phi1->data[i+1] - phi2->data[i+1], 1. );

    dt = XLALGPSGetREAL8(&t->data[i+1]) - XLALGPSGetREAL8(&t->data[i]);

    XLAL_CHECK_REAL8( XLALSinCos2PiLUT( &sp, &cp2, dp2 ) == XLAL_SUCCESS, XLAL_EFUNC );

    mismatch += (cp1 + cp2) * dt;
  }

  return (1. - fabs(mismatch)/(2.*T));
}


/**
 * \brief Get the position and velocity of the Earth at a given time
 *
 * This function will get the position and velocity of the Earth from the ephemeris data at the time t. It will be
 * returned in an EarthState structure. This is based on the start of the XLALBarycenterEarth function.
 */
void get_earth_pos_vel( EarthState *earth, EphemerisData *edat, LIGOTimeGPS *tGPS){
  REAL8 tgps[2];

  REAL8 t0e;        /* time since first entry in Earth ephem. table */
  INT4 ientryE;     /* entry in look-up table closest to current time, tGPS */

  REAL8 tinitE;     /* time (GPS) of first entry in Earth ephem table */
  REAL8 tdiffE;     /* current time tGPS minus time of nearest entry in Earth ephem look-up table */
  REAL8 tdiff2E;    /* tdiff2 = tdiffE * tdiffE */

  INT4 j;

  /* check input */
  if ( !earth || !tGPS || !edat || !edat->ephemE || !edat->ephemS ) {
    XLALPrintError ("%s: invalid NULL input 'earth', 'tGPS', 'edat','edat->ephemE' or 'edat->ephemS'\n", __func__ );
    XLAL_ERROR_VOID( XLAL_EINVAL );
  }

  tgps[0] = (REAL8)tGPS->gpsSeconds; /* convert from INT4 to REAL8 */
  tgps[1] = (REAL8)tGPS->gpsNanoSeconds;

  tinitE = edat->ephemE[0].gps;

  t0e = tgps[0] - tinitE;
  ientryE = ROUND(t0e/edat->dtEtable);  /* finding Earth table entry */

  if ( ( ientryE < 0 ) || ( ientryE >=  edat->nentriesE )) {
    XLALPrintError ("%s: input GPS time %f outside of Earth ephem range [%f, %f]\n", __func__, tgps[0], tinitE, tinitE +
edat->nentriesE * edat->dtEtable );
    XLAL_ERROR_VOID( XLAL_EDOM );
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
/*------------------------ END OF MODEL FUNCTIONS ----------------------------*/


/*----------------- FUNCTIONS TO CONVERT BETWEEN PARAMETERS ------------------*/
/**
 * \brief Convert \f$\phi_0\f$ and \f$\psi\f$ to a new coordinate system
 *
 * This function will convert the initial phase \f$\phi_0\f$ and polarisation angle \f$\psi\f$ into a new coordinate
 * system. As they are currently defined when \f$\psi\f$ wraps around at the limits of its range (\f$ \pm \pi/4 \f$
 * radians) it is equivalent to a \f$ \pi \f$ radians shift in \f$\phi_0\f$. This leads to a bimodal distribution in
 * \f$\phi_0\f$. A new coordinate system that is uni-modal and wraps around at the edges without introduction any phase
 * shift is given by:
 * \f[
 * \left( \begin{array}{c} {\phi'}_0 \\ {\psi}' \end{array} \right) =
 * \left( \begin{array}{cc} \sin{\theta} & \cos{\theta} \\ -\sin{\theta} &
 * \cos{\theta} \end{array} \right)
 * \left( \begin{array}{c} \phi_0 \\ \psi \end{array} \right),
 * \f]
 * where \f$\theta = \arctan{(1/2)}\f$.
 *
 * NOTE: This may want to be moved into LALInference at some point.
 *
 * \param phi0 [in] the initial phase parameter
 * \param psi [in] the polarisation angle parameter
 * \param phi0prime [in] the new coordinate axis
 * \param psiprime [in] the new coordinate axis
 */
void phi0_psi_transform( REAL8 phi0, REAL8 psi, REAL8 *phi0prime, REAL8 *psiprime ){
  REAL8 theta = atan2( 1., 2. );
  REAL8 st = sin(theta);
  REAL8 ct = cos(theta);

  /* check psi is in range: -pi/4 < psi < pi/4 */
  if( fabs(psi) > LAL_PI/4 ){
    XLALPrintError("Error... psi is not in range.\n");
    XLAL_ERROR_VOID(XLAL_EFUNC);
  }

  /* put phi0 in range -pi < phi0 < pi */
  if ( phi0 > 2.*LAL_PI ) { phi0 = fmod(phi0, LAL_TWOPI); }
  else { phi0 = LAL_TWOPI - fmod(LAL_TWOPI-phi0, LAL_TWOPI); }
  phi0 -= LAL_PI;

  *phi0prime = (st*phi0 + ct*psi);
  *psiprime = (-st*phi0 + ct*psi);
}


/**
 * \brief Convert new \f${\phi'}_0\f$ and \f$\psi'\f$ coordinate system back
 * to \f$\phi_0\f$ and \f$\psi\f$
 *
 * This function will convert the new parameters \f${\phi'}_0\f$ and \f$\psi'\f$, defined in \c phi0_psi_transform()
 * into the original \f$\phi_0\f$ and \f$\psi\f$ coordinates. This is done through the inverse transform:
 * \f{eqnarray}{
 * \left( \begin{array}{c} {\phi}_0 \\ {\psi} \end{array} \right) & = &
 * \left( \begin{array}{cc} \sin{\theta} & \cos{\theta} \\ -\sin(\theta) &
 * \cos{\theta} \end{array} \right)^{-1}
 * \left( \begin{array}{c} {\phi'}_0 \\ {\psi'} \end{array} \right), \\
 * & = & \left( \begin{array}{cc} \frac{1}{2\sin{\theta}} &
 * -\frac{1}{2\sin{\theta}} \\ \frac{1}{2\cos{\theta}} &
 * \frac{1}{2\cos{\theta}} \end{array} \right)
 * \left( \begin{array}{c} {\phi'}_0 \\ {\psi'} \end{array} \right),
 * \f}
 * where \f$\theta = \arctan{(1/2)}\f$.
 *
 * The \f${\phi'}_0\f$ and \f$\psi'\f$ should both be in the range \f$ \pm (\pi/2)\cos{\theta}\f$, which will return
 * \f$\psi\f$ in the range \f$ \pm \pi/2 \f$, and \f$\phi_0\f$ in the range \f$ \pm \pi \f$. These will then be
 * converted back into their original ranges.
 *
 * NOTE: This may want to be moved into LALInference at some point.
 *
 * \param phi0prime [in] the new coordinate axis
 * \param psiprime [in] the new coordinate axis
 * \param phi0 [in] the initial phase parameter
 * \param psi [in] the polarisation angle parameter
 */
void inverse_phi0_psi_transform( REAL8 phi0prime, REAL8 psiprime, REAL8 *phi0, REAL8 *psi ){
  REAL8 theta = atan2( 1., 2. );
  REAL8 ct = cos( theta );
  REAL8 o2st = 1. / ( 2. * sin( theta ) );
  REAL8 o2ct = 1. / ( 2. * ct );
  REAL8 phitmp = 0., psitmp = 0.;

  /* check psiprime and phi0prime is in range +/- (pi/2)cos(theta) */
  if ( fabs(phi0prime) > LAL_PI_2*ct || fabs(psiprime) > LAL_PI_2*ct ){
    fprintf(stderr, "phi0prime = %le, psiprime = %le\n", phi0prime, psiprime);
    XLALPrintError("Error... phi0prime or psiprime are not in range\n");
    XLAL_ERROR_VOID(XLAL_EFUNC);
  }

  phitmp = o2st*phi0prime - o2st*psiprime;
  psitmp = o2ct*phi0prime + o2ct*psiprime;

  /* get psi into +/- pi/4 range */
  if ( fabs(psitmp) > LAL_PI/4. ){
    phitmp += LAL_PI; /* rotate phase by pi */

    /* wrap around psi */
    if ( psitmp > LAL_PI/4. ) { psitmp = -(LAL_PI/4.) + fmod(psitmp+(LAL_PI/4.), LAL_PI_2); }
    else { psitmp = (LAL_PI/4.) - fmod((LAL_PI/4.)-psitmp, LAL_PI_2); }
  }

  *psi = psitmp;

  /* get phi0 into 0 -> 2pi range */
  if ( phitmp > LAL_TWOPI ) { phitmp = fmod(phitmp, LAL_TWOPI); }
  else { phitmp = LAL_TWOPI - fmod(LAL_TWOPI-phitmp, LAL_TWOPI); }

  *phi0 = phitmp;
}


/**
 * \brief Convert sources parameters into amplitude and phase notation parameters
 *
 * Convert the physical source parameters into the amplitude and phase notation given in Eqns
 * 76-79 of LIGO T1200265-v3.
 */
void invert_source_params( BinaryPulsarParams *params ){
  REAL8 sinlambda, coslambda, sinlambda2, coslambda2, sin2lambda;
  REAL8 theta, sintheta, costheta, costheta2, sintheta2, sin2theta;
  REAL8 phi0 = params->phi0;

  sinlambda = sin( params->lambda );
  coslambda = cos( params->lambda );
  sin2lambda = sin( 2. * params->lambda );
  sinlambda2 = SQUARE( sinlambda );
  coslambda2 = SQUARE( coslambda );

  theta = acos( params->costheta );
  sintheta = sin( theta );
  costheta = params->costheta;
  sin2theta = sin( 2. * theta );
  sintheta2 = SQUARE( sintheta );
  costheta2 = SQUARE( costheta );

  REAL8 A22 = params->I21 * ( sinlambda2 - coslambda2 * costheta2 ) - params->I31 * sintheta2;
  REAL8 B22 = params->I21 * sin2lambda * costheta;
  REAL8 A222 = SQUARE( A22 );
  REAL8 B222 = SQUARE( B22 );

  REAL8 A21 = params->I21 * sin2lambda * sintheta;
  REAL8 B21 = sin2theta * ( params->I21 * coslambda2 - params->I31 );
  REAL8 A212 = SQUARE( A21 );
  REAL8 B212 = SQUARE( B21 );

  params->C22 = 2.*sqrt( A222 + B222 );
  params->C21 = 2.*sqrt( A212 + B212 );

  params->phi22 = fmod( phi0 - atan2( B22, A22 ), LAL_TWOPI );
  params->phi21 = fmod( ( phi0/2. ) - atan2( B21, A21 ), LAL_TWOPI );
}

