/**
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin, Colin Gill, John Veitch
 *
 * \brief Pulsar model functions for use in parameter estimation codes for targeted pulsar searches.
 */

#include "ppe_models.h"
#include <lal/SinCosLUT.h>

#ifndef _OPENMP
#define omp ignore
#endif

#define SQUARE(x) ( (x) * (x) )

/******************************************************************************/
/*                            MODEL FUNCTIONS                                 */
/******************************************************************************/

/**
 * \brief Defines the pulsar model/template to use
 *
 * This function is the wrapper for functions defining the pulsar model template to be used in the analysis.
 * It places them into a \c PulsarParameters structure.
 *
 * Note: Any additional models should be added into this function.
 *
 * \param model [in] The model structure hold model information and current parameter info
 *
 * \sa add_pulsar_parameter
 * \sa pulsar_model
 */
void get_pulsar_model( LALInferenceModel *model ){
  PulsarParameters *pars = XLALCalloc(sizeof(*pars),1);

  /* set model parameters (including rescaling) */
  add_pulsar_parameter( model->params, pars, "PSI" );

  if( ( LALInferenceCheckVariableNonFixed( model->params, "H0" ) || LALInferenceCheckVariable( model->ifo->params, "jones-model" ) ) && !LALInferenceCheckVariable( model->ifo->params, "nonGR" ) ){
    add_pulsar_parameter( model->params, pars, "H0" );

    /* use parameterisation from Ian Jones's original model */
    add_pulsar_parameter( model->params, pars, "I21" );
    add_pulsar_parameter( model->params, pars, "I31" );
    add_pulsar_parameter( model->params, pars, "LAMBDA" );

    /* check whether using cos(theta) or theta as the variable (defaulting to cos(theta) being the value that is set */
    if ( LALInferenceCheckVariableNonFixed( model->params, "THETA" ) ){
      REAL8 costheta = cos(LALInferenceGetREAL8Variable( model->params, "THETA" ));
      PulsarAddParam( pars, "COSTHETA", &costheta, PULSARTYPE_REAL8_t );
    }
    else{ add_pulsar_parameter( model->params, pars, "COSTHETA" ); }
    add_pulsar_parameter( model->params, pars, "PHI0" ); /* note that this is the rotational phase */

    /* check whether using cos(theta) or theta as the variable (defaulting to cos(theta) being the value that is set */
    if ( LALInferenceCheckVariableNonFixed( model->params, "IOTA" ) ){
      REAL8 cosiota = cos(LALInferenceGetREAL8Variable( model->params, "IOTA" ));
      PulsarAddParam( pars, "COSIOTA", &cosiota, PULSARTYPE_REAL8_t );
    }
    else{ add_pulsar_parameter( model->params, pars, "COSIOTA" ); }

    invert_source_params( pars );
  }
  else if ( LALInferenceCheckVariable( model->ifo->params, "nonGR" ) ){
    /* speed of GWs as (1 - fraction of speed of light LAL_C_SI) */
    add_pulsar_parameter( model->params, pars, "CGW" );

    /* amplitudes for use with non-GR searches */
    /* tensor modes */
    add_pulsar_parameter( model->params, pars, "HPLUS" );
    add_pulsar_parameter( model->params, pars, "HCROSS" );
    /* scalar modes */
    add_pulsar_parameter( model->params, pars, "HSCALARB" );
    add_pulsar_parameter( model->params, pars, "HSCALARL" );
    /* vector modes */
    add_pulsar_parameter( model->params, pars, "HVECTORX" );
    add_pulsar_parameter( model->params, pars, "HVECTORY" );

    add_pulsar_parameter( model->params, pars, "PHI0SCALAR" );
    add_pulsar_parameter( model->params, pars, "PSISCALAR" );
    add_pulsar_parameter( model->params, pars, "PHI0VECTOR" );
    add_pulsar_parameter( model->params, pars, "PSIVECTOR" );
    add_pulsar_parameter( model->params, pars, "PHI0TENSOR" );
    add_pulsar_parameter( model->params, pars, "PSITENSOR" );

    /* parameters that might be needed for particular models */
    add_pulsar_parameter( model->params, pars, "H0" );
    add_pulsar_parameter( model->params, pars, "IOTA" );
    add_pulsar_parameter( model->params, pars, "COSIOTA" );

    /* check whether a specific nonGR model was requested */
    if ( LALInferenceCheckVariable( model->ifo->params, "nonGRmodel" ) ){
      char* nonGRmodel;
      nonGRmodel = *(char**)LALInferenceGetVariable( model->ifo->params, "nonGRmodel" );
      set_nonGR_model_parameters( pars, nonGRmodel );
    }
  }
  else{
    /* amplitude for usual GR search */
    add_pulsar_parameter( model->params, pars, "C21" );
    add_pulsar_parameter( model->params, pars, "C22" );
    add_pulsar_parameter( model->params, pars, "PHI21" );

    /* check whether using cos(theta) or theta as the variable (defaulting to cos(theta) being the value that is set */
    if ( LALInferenceCheckVariableNonFixed( model->params, "IOTA" ) ){
      REAL8 cosiota = cos(LALInferenceGetREAL8Variable( model->params, "IOTA" ));
      PulsarAddParam( pars, "COSIOTA", &cosiota, PULSARTYPE_REAL8_t );
    }
    else{ add_pulsar_parameter( model->params, pars, "COSIOTA" ); }

    if( LALInferenceCheckVariable( model->ifo->params, "biaxial" ) ){
      /* use complex amplitude parameterisation, but set up for a biaxial star */
      REAL8 phi22 = 2.*LALInferenceGetREAL8Variable( model->params, "PHI21" );
      PulsarAddParam( pars, "PHI22", &phi22, PULSARTYPE_REAL8_t );
    }
    else{
      add_pulsar_parameter( model->params, pars, "PHI22" );
    }
  }

  /* set the potentially variable parameters */
  add_pulsar_parameter( model->params, pars, "PEPOCH" );
  add_pulsar_parameter( model->params, pars, "POSEPOCH" );

  add_pulsar_parameter( model->params, pars, "RA" );
  add_pulsar_parameter( model->params, pars, "PMRA" );
  add_pulsar_parameter( model->params, pars, "DEC" );
  add_pulsar_parameter( model->params, pars, "PMDEC" );

  add_pulsar_parameter( model->params, pars, "F0" );
  add_pulsar_parameter( model->params, pars, "F1" );
  add_pulsar_parameter( model->params, pars, "F2" );
  add_pulsar_parameter( model->params, pars, "F3" );
  add_pulsar_parameter( model->params, pars, "F4" );
  add_pulsar_parameter( model->params, pars, "F5" );

  /* check if there are binary parameters */
  if( LALInferenceCheckVariable(model->ifo->params, "BINARY") ){
    /* binary system model - NOT pulsar model */
    CHAR *binary = NULL;
    binary = XLALStringDuplicate(*(CHAR**)LALInferenceGetVariable( model->ifo->params, "BINARY" ));
    PulsarAddParam( pars, "BINARY", &binary, PULSARTYPE_string_t );

    add_pulsar_parameter( model->params, pars, "ECC" );
    add_pulsar_parameter( model->params, pars, "OM" );
    add_pulsar_parameter( model->params, pars, "PB" );
    add_pulsar_parameter( model->params, pars, "A1" );
    add_pulsar_parameter( model->params, pars, "T0" );

    add_pulsar_parameter( model->params, pars, "ECC_2" );
    add_pulsar_parameter( model->params, pars, "OM_2" );
    add_pulsar_parameter( model->params, pars, "PB_2" );
    add_pulsar_parameter( model->params, pars, "A1_2" );
    add_pulsar_parameter( model->params, pars, "T0_2" );

    add_pulsar_parameter( model->params, pars, "ECC_3" );
    add_pulsar_parameter( model->params, pars, "OM_3" );
    add_pulsar_parameter( model->params, pars, "PB_3" );
    add_pulsar_parameter( model->params, pars, "A1_3" );
    add_pulsar_parameter( model->params, pars, "T0_3" );

    add_pulsar_parameter( model->params, pars, "XPBDOT" );
    add_pulsar_parameter( model->params, pars, "EPS1" );
    add_pulsar_parameter( model->params, pars, "EPS2" );
    add_pulsar_parameter( model->params, pars, "EPS1DOT" );
    add_pulsar_parameter( model->params, pars, "EPS2DOT" );
    add_pulsar_parameter( model->params, pars, "TASC" );

    add_pulsar_parameter( model->params, pars, "OMDOT" );
    add_pulsar_parameter( model->params, pars, "GAMMA" );
    add_pulsar_parameter( model->params, pars, "PBDOT" );
    add_pulsar_parameter( model->params, pars, "XDOT" );
    add_pulsar_parameter( model->params, pars, "EDOT" );

    add_pulsar_parameter( model->params, pars, "SINI" );
    add_pulsar_parameter( model->params, pars, "DR" );
    add_pulsar_parameter( model->params, pars, "DTHETA" );
    add_pulsar_parameter( model->params, pars, "A0" );
    add_pulsar_parameter( model->params, pars, "B0" );

    add_pulsar_parameter( model->params, pars, "MTOT" );
    add_pulsar_parameter( model->params, pars, "M2" );

    if ( LALInferenceCheckVariable( model->params, "FB" ) ){
      REAL8Vector *fb = NULL;
      fb = *(REAL8Vector **)LALInferenceGetVariable( model->params, "FB" );
      PulsarAddParam( pars, "FB", &fb, PULSARTYPE_REAL8Vector_t );
    }
  }

  pulsar_model( pars, model->ifo );

  PulsarFreeParams( pars ); /* free memory */
}

/**
 * \brief Set amplitude parameters for specific non-GR models.
 *
 * Turns physical parameters from a particular nonGR model into the corresponding antenna pattern amplitudes and phases. All nonGR models must be included here.
 * \param pars [in] parameter structure
 * \param nonGRmodel [in] name of model requested
 */
void set_nonGR_model_parameters( PulsarParameters *pars, char* nonGRmodel ){
  /* determine what model was requested */
  int isG4v = strcmp(nonGRmodel, "G4v") * strcmp(nonGRmodel, "g4v") * strcmp(nonGRmodel, "G4V");
  int isST = strcmp(nonGRmodel, "scalar-tensor") * strcmp(nonGRmodel, "ST");
  REAL8 h0 = PulsarGetREAL8ParamOrZero( pars, "H0" );

  if( isG4v == 0 ){
    /* \f$ h_{\rm x} = h_0 \sin \iota~,~\phi_{\rm x} = -\pi/2 \f$ */
    /* \f$ h_{\rm y} = h_0 \sin \iota \cos \iota~,~\phi_{\rm y} = 0 */
    REAL8 iota = PulsarGetREAL8ParamOrZero( pars, "IOTA" );
    REAL8 cosiota = cos(iota);
    REAL8 siniota = sin(iota);

    REAL8 hVectorX = h0 * siniota;
    REAL8 hVectorY = h0 * siniota * cosiota;
    REAL8 psiVector = LAL_PI_2;
    PulsarAddParam( pars, "HVECTORX", &hVectorX, PULSARTYPE_REAL8_t );
    PulsarAddParam( pars, "HVECTORY", &hVectorY, PULSARTYPE_REAL8_t );
    PulsarAddParam( pars, "PSIVECTOR", &psiVector, PULSARTYPE_REAL8_t );
  }
  else if( isST == 0 ) {
    /* GR plus an unconstrained breathing mode */
    /* TO DO: are there specific functional relations linking the HSCALARB to pulsar parameters like COSIOTA, etc.? */
    REAL8 cosiota = PulsarGetREAL8ParamOrZero( pars, "COSIOTA" );
    REAL8 hPlus = 0.5 * h0 * (1. + cosiota * cosiota);
    REAL8 hCross = h0 * cosiota;
    REAL8 psiTensor = -LAL_PI_2;
    PulsarAddParam( pars, "HPLUS", &hPlus, PULSARTYPE_REAL8_t );
    PulsarAddParam( pars, "HCROSS", &hCross, PULSARTYPE_REAL8_t );
    PulsarAddParam( pars, "PSITENSOR", &psiTensor, PULSARTYPE_REAL8_t );
  } else {
    XLALPrintError ("%s: unrecognized non-GR model. Currently supported: scalar-tensor (ST), G4v, or no argument for full search.\n", __func__ );
    XLAL_ERROR_VOID( XLAL_EINVAL );
  }
}


/**
 * \brief Add a \c REAL8 parameter from a \c LALInfernceVariable into a \c PulsarParameters variable
 *
 * \param var [in] \c LALInferenceVariable structure
 * \param params [in] \c PulsarParameters structure
 * \param parname [in] name of the parameter
 */
void add_pulsar_parameter( LALInferenceVariables *var, PulsarParameters *params, const CHAR *parname ){
  REAL8 par = LALInferenceGetREAL8Variable( var, parname );
  PulsarAddParam( params, parname, &par, PULSARTYPE_REAL8_t );
}


/**
 * \brief Add a \c REAL8 parameter from a \c PulsarParameters variable into a \c LALInfernceVariable
 *
 * This function will rescale a parameter to its true value using the scale factor and minimum scale value.
 *
 * \param params [in] \c PulsarParameters structure
 * \param var [in] \c LALInferenceVariable structure
 * \param parname [in] name of the parameter
 * \param vary [in] the parameter \c LALInferenceParamVaryType
 */
void add_variable_parameter( PulsarParameters *params, LALInferenceVariables *var, const CHAR *parname, LALInferenceParamVaryType vary ){
  REAL8 par = PulsarGetREAL8ParamOrZero( params, parname );
  LALInferenceAddVariable( var, parname, &par, LALINFERENCE_REAL8_t, vary );
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
 * (stored in \c ifo->timeData->data) will be calculated and the complex signal model, \f$M\f$, modified accordingly:
 * \f[
 * M'(t) = M(t)\exp{i(-(\phi(t)_n - \phi(t)_h))}.
 * \f]
 *
 * \param params [in] A \c BinaryPulsarParams structure containing the model parameters
 * \param ifo [in] The ifo model structure containing the detector paramters and buffers
 *
 * \sa get_amplitude_model
 * \sa get_phase_model
 */
void pulsar_model( PulsarParameters *params, LALInferenceIFOModel *ifo ){
  INT4 i = 0, length = 0;
  UINT4 j = 0;
  LALInferenceIFOModel *ifomodel1 = ifo, *ifomodel2 = ifo;

  /* get the amplitude model */
  get_amplitude_model( params, ifomodel1 );

  /* check whether to search over the phase parameters or not - this only needs to be set for the
   * first ifo linked list in at set for a given detector (i.e. it doesn't need to be set for
   * different frequency streams */
  if ( LALInferenceCheckVariable( ifomodel2->params, "varyphase" ) ) {
    /* get difference in phase for f component and perform extra heterodyne */
    REAL8Vector *freqFactors = NULL;
    freqFactors = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "freqfactors" );

    while ( ifomodel2 ){
      for( j = 0; j < freqFactors->length; j++ ){
        REAL8Vector *dphi = NULL;

        length = ifomodel2->compTimeSignal->data->length;

        /* the timeData vector within the LALIFOModel structure contains the phase calculated using the initial (heterodyne)
         * values of the phase parameters */

        /* reheterodyne with the phase */
        if ( (dphi = get_phase_model( params, ifomodel2, freqFactors->data[j] )) != NULL ){
          #pragma omp parallel for
          for( i=0; i<length; i++ ){
            COMPLEX16 M;
            REAL8 dphit;
            COMPLEX16 expp;

            dphit = fmod(dphi->data[i] - ifomodel2->timeData->data->data[i], 1.);
            expp = cexp( LAL_TWOPI * I * dphit );

            M = ifomodel2->compTimeSignal->data->data[i];

            /* heterodyne */
            ifomodel2->compTimeSignal->data->data[i] = M * expp;
          }

          XLALDestroyREAL8Vector( dphi );
        }

        ifomodel2 = ifomodel2->next;
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
 * \param ifo [in] The ifo model structure containing the detector parameters and buffers
 * \param freqFactor [in] the multiplicative factor on the pulsar frequency for a particular model
 *
 * \return A vector of rotational phase values
 *
 * \sa get_ssb_delay
 * \sa get_bsb_delay
 */
REAL8Vector *get_phase_model( PulsarParameters *params, LALInferenceIFOModel *ifo, REAL8 freqFactor ){
  INT4 i = 0, length = 0;

  REAL8 DT = 0., deltat = 0., deltat2 = 0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */

  REAL8Vector *phis = NULL, *dts = NULL, *bdts = NULL;
  LIGOTimeGPSVector *datatimes = NULL;

  REAL8 T0 = PulsarGetREAL8ParamOrZero(params, "PEPOCH"); /*time of ephem info*/
  REAL8 cgw = PulsarGetREAL8ParamOrZero(params, "CGW");

  /* check if we want to calculate the phase at a the downsampled rate */
  datatimes = ifo->times;

  /* if edat is NULL then return a NULL pointer */
  if( ifo->ephem == NULL ) return NULL;

  length = datatimes->length;

  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );

  /* get time delays */
  if( (dts = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "ssb_delays" )) == NULL ||
       LALInferenceCheckVariable( ifo->params, "varyskypos" ) ){
    /* get time delays with an interpolation of interptime (30 mins) */
    dts = get_ssb_delay( params, datatimes, ifo->ephem, ifo->tdat, ifo->ttype, ifo->detector, interptime );
  }

  if( (bdts = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "bsb_delays" )) == NULL ||
       LALInferenceCheckVariable( ifo->params, "varybinary" ) ){
    /* get binary system time delays */
    bdts = get_bsb_delay( params, datatimes, dts, ifo->ephem );
  }

  #pragma omp parallel for private(DT, deltat, deltat2)
  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &datatimes->data[i] ); /*time of data*/

    DT = realT - T0; /*time diff between data and ephem info*/

    if ( PulsarCheckParam( params, "BINARY" ) ) { deltat = DT + dts->data[i] + bdts->data[i]; }
    else { deltat = DT + dts->data[i]; }

    /* correct for speed of GW compared to speed of light */
    if ( cgw > 0.0 && cgw < 1. ) { deltat /= cgw; }

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = freqFactor*deltat*(PulsarGetREAL8ParamOrZero(params, "F0") +
      0.5*PulsarGetREAL8ParamOrZero(params, "F1")*deltat +
      (1./6.)*PulsarGetREAL8ParamOrZero(params, "F2")*deltat2 +
      (1./24.)*PulsarGetREAL8ParamOrZero(params, "F3")*deltat*deltat2 +
      (1./120.)*PulsarGetREAL8ParamOrZero(params, "F4")*deltat2*deltat2 +
      (1./720.)*PulsarGetREAL8ParamOrZero(params, "F5")*deltat2*deltat2*deltat);
  }

  /* free memory */
  if ( !LALInferenceCheckVariable( ifo->params, "ssb_delays") || LALInferenceCheckVariable( ifo->params, "varyskypos" ) )
    XLALDestroyREAL8Vector( dts );

  if ( !LALInferenceCheckVariable( ifo->params, "bsb_delays") || LALInferenceCheckVariable( ifo->params, "varybinary" ) )
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
REAL8Vector *get_ssb_delay( PulsarParameters *pars, LIGOTimeGPSVector *datatimes, EphemerisData *ephem,
                            TimeCorrectionData *tdat, TimeCorrectionType ttype, LALDetector *detector,
                            REAL8 interptime ){
  INT4 i = 0, length = 0;
  REAL8 T0 = 0.;

  BarycenterInput bary;
  UINT4 baryfail = 0;

  REAL8Vector *dts = NULL;

  /* if edat is NULL then return a NULL poniter */
  if( ephem == NULL ) { return NULL; }

  /* copy barycenter and ephemeris data */
  bary.site.location[0] = detector->location[0]/LAL_C_SI;
  bary.site.location[1] = detector->location[1]/LAL_C_SI;
  bary.site.location[2] = detector->location[2]/LAL_C_SI;

  REAL8 ra = 0.;
  if ( PulsarCheckParam( pars, "RA" ) ) { ra = PulsarGetREAL8Param( pars, "RA" ); }
  else if ( PulsarCheckParam( pars, "RAJ" ) ) { ra = PulsarGetREAL8Param( pars, "RAJ" ); }
  else {
    XLALPrintError ("%s: No source right ascension specified!", __func__ );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }
  REAL8 dec = 0.;
  if ( PulsarCheckParam( pars, "DEC" ) ) { dec = PulsarGetREAL8Param( pars, "DEC" ); }
  else if ( PulsarCheckParam( pars, "DECJ" ) ) { dec = PulsarGetREAL8Param( pars, "DECJ" ); }
  else {
    XLALPrintError ("%s: No source declination specified!", __func__ );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }
  REAL8 pmra = PulsarGetREAL8ParamOrZero( pars, "PMRA" );
  REAL8 pmdec = PulsarGetREAL8ParamOrZero( pars, "PMDEC" );
  REAL8 pepoch = PulsarGetREAL8ParamOrZero( pars, "PEPOCH" );
  REAL8 posepoch = PulsarGetREAL8ParamOrZero( pars, "POSEPOCH" );
  REAL8 px = PulsarGetREAL8ParamOrZero( pars, "PX" );     /* parallax */
  REAL8 dist = PulsarGetREAL8ParamOrZero( pars, "DIST" ); /* distance */

   /* set the position and frequency epochs if not already set */
  if( pepoch == 0. && posepoch != 0.) { pepoch = posepoch; }
  else if( posepoch == 0. && pepoch != 0. ) { posepoch = pepoch; }

  length = datatimes->length;

  /* allocate memory for times delays */
  dts = XLALCreateREAL8Vector( length );

  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( px != 0. ) { bary.dInv = px*1e-3*LAL_C_SI/LAL_PC_SI; }
  else if( dist != 0. ) { bary.dInv = LAL_C_SI/(dist*1e3*LAL_PC_SI); }
  else { bary.dInv = 0.; }

  T0 = pepoch;

  /* make sure ra and dec are wrapped withing 0--2pi and -pi.2--pi/2 respectively */
  ra = fmod(ra, LAL_TWOPI);
  REAL8 absdec = fabs(dec);
  if ( absdec > LAL_PI_2 ){
    UINT4 nwrap = floor((absdec+LAL_PI_2)/LAL_PI);
    dec = (dec > 0 ? 1. : -1.)*(nwrap%2 == 1 ? -1. : 1.)*(fmod(absdec + LAL_PI_2, LAL_PI) - LAL_PI_2);
    ra = fmod(ra + (REAL8)nwrap*LAL_PI, LAL_TWOPI); /* move RA by pi */
  }

  #pragma omp parallel for shared(baryfail) private(bary)
  for( i=0; i<length; i++){
    if ( !baryfail ){ /* need this to quickly exit loop while using openmp (as can't use break statements) */
      EarthState earth, earth2;
      EmissionTime emit, emit2;

      REAL8 DT = 0., DTplus = 0., realT = XLALGPSGetREAL8( &datatimes->data[i] );

      DT = realT - T0;

      /* only do call to the barycentring routines once every interptime (unless
         interptime == 0), otherwise just linearly interpolate between them */
      if( i == 0 || DT > DTplus || interptime == 0 ){
        bary.tgps = datatimes->data[i];

        bary.delta = dec + ( realT - posepoch ) * pmdec;
        bary.alpha = ra + ( realT - posepoch ) * pmra / cos( bary.delta );

        /* call barycentring routines */
        if ( XLALBarycenterEarthNew( &earth, &bary.tgps, ephem, tdat, ttype ) != XLAL_SUCCESS ){ baryfail = 1; continue; }
        if ( XLALBarycenter( &emit, &bary, &earth ) != XLAL_SUCCESS ){ baryfail = 1; continue; }

        /* add interptime to the time */
        if ( interptime > 0 ){
          DTplus = DT + interptime;
          XLALGPSAdd( &bary.tgps, interptime );

          /* No point in updating the positions as difference will be tiny */
          if( XLALBarycenterEarthNew( &earth2, &bary.tgps, ephem, tdat, ttype ) != XLAL_SUCCESS ){ baryfail = 1; continue; }
          if( XLALBarycenter( &emit2, &bary, &earth2 ) != XLAL_SUCCESS){ baryfail = 1; continue; }
        }
      }

      /* linearly interpolate to get emitdt */
      if( interptime > 0. ){
        dts->data[i] = emit.deltaT + (DT - (DTplus - interptime)) * (emit2.deltaT -  emit.deltaT)/interptime;
      }
      else { dts->data[i] = emit.deltaT; }
    }
  }

  /* check that for loop didn't fail early */
  XLAL_CHECK_NULL( !baryfail, XLAL_EFUNC );

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
REAL8Vector *get_bsb_delay( PulsarParameters *pars, LIGOTimeGPSVector *datatimes, REAL8Vector *dts,
                            EphemerisData *edat ){
  REAL8Vector *bdts = NULL;

  INT4 i = 0, length = datatimes->length;

  bdts = XLALCreateREAL8Vector( length );

  #pragma omp parallel for
  for ( i = 0; i < length; i++ ){
    /* check whether there's a binary model */
    if ( PulsarCheckParam( pars, "BINARY" ) ){
      BinaryPulsarInput binput;
      BinaryPulsarOutput boutput;
      EarthState earth;

      binput.tb = XLALGPSGetREAL8( &datatimes->data[i] ) + dts->data[i];

      get_earth_pos_vel( &earth, edat, &datatimes->data[i] );

      binput.earth = earth; /* current Earth state */
      XLALBinaryPulsarDeltaTNew( &boutput, &binput, pars );
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
 * The antenna pattern functions a(t) and b(t) are contained in a 1D lookup table, so within this function
 * the correct value for the given time is interpolated from this lookup table using linear interpolation.
 *
 * \param pars [in] A set of pulsar parameters
 * \param ifo [in] The ifo parameters giving information on the data and detector
 *
 * DEPRECATED: use \c get_amplitude_model instead.
 */
void get_triaxial_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo ){
  INT4 i = 0, length;

  REAL8 tsteps, tsv;
  INT4 timebinMin, timebinMax;
  REAL8 plus, cross, plusT, crossT;
  REAL8 plus00, plus01, cross00, cross01;
  REAL8 timeScaled;
  REAL8 timeMin, timeMax;
  REAL8 T;
  REAL8 Xplus, Xcross;
  COMPLEX16 expiphi, Xpexpphi, Xcexpphi;
  REAL8 twopsi = 0.;
  REAL4 c2psi = 0., s2psi = 0.;

  REAL8Vector *LU_Fplus = NULL, *LU_Fcross = NULL;
  REAL8Vector *sidDayFrac = NULL;

  length = ifo->times->length;

  /* set lookup table parameters */
  tsteps = *(INT4*)LALInferenceGetVariable( ifo->params, "timeSteps" );

  LU_Fplus = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "a_response" );
  LU_Fcross = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "b_response" );

  /* get the sidereal time since the initial data point % sidereal day */
  sidDayFrac = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "siderealDay" );

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
  twopsi = 2.*pars.psi;
  XLAL_CHECK_VOID( XLALSinCosLUT( &s2psi, &c2psi, twopsi ) == XLAL_SUCCESS, XLAL_EFUNC );

  tsv = LAL_DAYSID_SI / tsteps;

  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = sidDayFrac->data[i];
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;

    /* get values of matrix for linear interpolation */
    plus00 = LU_Fplus->data[timebinMin];
    plus01 = LU_Fplus->data[timebinMax];

    cross00 = LU_Fcross->data[timebinMin];
    cross01 = LU_Fcross->data[timebinMax];

    /* rescale time for linear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);

    plus = plus00 + (plus01-plus00)*timeScaled;
    cross = cross00 + (cross01-cross00)*timeScaled;

    plusT = plus*c2psi + cross*s2psi;
    crossT = cross*c2psi - plus*s2psi;

    /* create the complex signal amplitude model */
    ifo->compTimeSignal->data->data[i] = plusT*Xpexpphi - I*crossT*Xcexpphi;
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
 * sets of data per detector. In this model the \f$\phi_0\f$ parameter is the initial phase for the \f$l=m=2\f$ harmonic
 * (i.e. equivalent to the phase at 2f for the triaxial model).
 *
 * As for the standard triaxial model, the antenna pattern functions are contained in a 1D lookup table, so within this
 * function the correct value for the given time is interpolated from this lookup table using linear
 * interpolation.
 *
 * \param pars [in] A set of pulsar parameters
 * \param ifo  [in] The ifo model containing detector-specific parameters
 *
 * DEPRECATED: use \c get_amplitude_model instead.
 */
void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOModel *ifo ){
  INT4 i = 0, length;

  REAL8 tsteps, tsv;
  INT4 timebinMin, timebinMax;
  REAL8 plus00, plus01, cross00, cross01;
  REAL8 timeScaled;
  REAL8 timeMin, timeMax;
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

  REAL8Vector *LU_Fplus = NULL, *LU_Fcross = NULL;
  REAL8Vector *sidDayFrac1 = NULL;
  REAL8Vector *sidDayFrac2 = NULL;
  REAL8 twopsi = 0., plusT = 0., crossT = 0.;
  REAL4 c2psi = 0., s2psi = 0.;

  /* set lookup table parameters */
  tsteps = *(INT4*)LALInferenceGetVariable( ifo->params, "timeSteps" );

  twopsi = 2.*pars.psi;
  XLAL_CHECK_VOID( XLALSinCosLUT( &s2psi, &c2psi, twopsi ) == XLAL_SUCCESS, XLAL_EFUNC );

  LU_Fplus = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "a_response" );
  LU_Fcross = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "b_response");
  /* get the sidereal time since the initial data point % sidereal day */
  sidDayFrac1 = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "siderealDay" );

  /* phi0 here the GW phase for the l=m=2 mode */
  ePhi = cexp( 0.5 * pars.phi0 * I );
  e2Phi = cexp( pars.phi0 * I );

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

  tsv = LAL_DAYSID_SI / tsteps;

  /*--------------------------------------------------------------------------*/
  /* set model for 1f component */
  length = ifo->times->length;

  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = sidDayFrac1->data[i];
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;

    /* get values of matrix for linear interpolation */
    plus00 = LU_Fplus->data[timebinMin];
    plus01 = LU_Fplus->data[timebinMax];

    cross00 = LU_Fcross->data[timebinMin];
    cross01 = LU_Fcross->data[timebinMax];

    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);

    plus = plus00 + (plus01-plus00)*timeScaled;
    cross = cross00 + (cross01-cross00)*timeScaled;

    plusT = plus*c2psi + cross*s2psi;
    crossT = cross*c2psi - plus*s2psi;

    /* create the complex signal amplitude model at f */
    ifo->compTimeSignal->data->data[i] = ( plusT * Xplusf * ePhi * ( A21 - I * B21 ) ) +
                                         ( crossT * Xcrossf * ePhi * ( B21 + I * A21 ) );
  }

  /*--------------------------------------------------------------------------*/
  /* set model for 2f component */
  length = ifo->next->times->length;

  sidDayFrac2 = *(REAL8Vector**)LALInferenceGetVariable( ifo->next->params, "siderealDay" );

  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = sidDayFrac2->data[i];
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;

    /* get values of matrix for linear interpolation */
    plus00 = LU_Fplus->data[timebinMin];
    plus01 = LU_Fplus->data[timebinMax];

    cross00 = LU_Fcross->data[timebinMin];
    cross01 = LU_Fcross->data[timebinMax];

    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);

    plus = plus00 + (plus01-plus00)*timeScaled;
    cross = cross00 + (cross01-cross00)*timeScaled;

    plusT = plus*c2psi + cross*s2psi;
    crossT = cross*c2psi - plus*s2psi;

    /* create the complex signal amplitude model at 2f */
    ifo->next->compTimeSignal->data->data[i] = ( plusT * Xplus2f * e2Phi * ( A22 - I * B22 ) ) +
                                               ( crossT * Xcross2f * e2Phi * ( B22 + I * A22 ) );
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
 * The antenna pattern functions are contained in a 1D lookup table, so within this function the correct value for the
 * given time is interpolated from this lookup table using linear interpolation..
 *
 * \param pars [in] A set of pulsar parameters
 * \param ifo  [in] The ifo model containing detector-specific parameters
 *
 */
void get_amplitude_model( PulsarParameters *pars, LALInferenceIFOModel *ifo ){
  UINT4 i = 0, j = 0, length;

  REAL8 T, twopsi;
  REAL8 cosiota = PulsarGetREAL8ParamOrZero( pars, "COSIOTA" );
  REAL8 siniota = sin(acos(cosiota));
  REAL8 s2psi = 0., c2psi = 0., spsi = 0., cpsi = 0.;
  UINT4 nonGR = 0;

  REAL8Vector *freqFactors = NULL;
  INT4 varyphase = 0, roq = 0;

  freqFactors = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "freqfactors" );

  if( LALInferenceCheckVariable( ifo->params, "varyphase" ) ){ varyphase = 1; }
  if( LALInferenceCheckVariable( ifo->params, "roq" ) ){ roq = 1; }

  twopsi = 2.*PulsarGetREAL8ParamOrZero( pars, "PSI" );
  s2psi = sin(twopsi);
  c2psi = cos(twopsi);

  /* check for non-GR model */
  if ( LALInferenceCheckVariable( ifo->params, "nonGR" ) ){
    nonGR = *(UINT4*)LALInferenceGetVariable( ifo->params, "nonGR" );
  }

  if ( nonGR == 1 ){
    spsi = sin(PulsarGetREAL8ParamOrZero( pars, "PSI" ));
    cpsi = cos(PulsarGetREAL8ParamOrZero( pars, "PSI" ));
  }

  if ( nonGR == 1 && freqFactors->length > 1 ){
    XLALPrintError( "%s: Error... currently can only use non-GR parameters for l=m=2 harmonic.\n", __func__ );
    XLAL_ERROR_VOID( XLAL_EFAILED );
  }

  /* loop over all detectors */
  while( ifo ){
    /* loop over components in data as given by the frequency factors */
    for( j = 0; j < freqFactors->length; j++ ){
      COMPLEX16 expPhi;
      COMPLEX16 Cplus = 0., Ccross = 0., Cx = 0., Cy = 0., Cl = 0., Cb = 0.;

      if ( !ifo ){
        XLALPrintError( "%s: Error... ifo model not defined.\n", __func__ );
        XLAL_ERROR_VOID( XLAL_EINVAL );
      }

      /* get the amplitude and phase factors */
      if( freqFactors->data[j] == 1. ){
        /* the l=2, m=1 harmonic at the rotation frequency */
        expPhi = cexp( I * PulsarGetREAL8ParamOrZero( pars, "PHI21" ));
        Cplus = -0.25 * PulsarGetREAL8ParamOrZero( pars, "C21" ) * siniota * cosiota * expPhi;
        Ccross = 0.25 * I * PulsarGetREAL8ParamOrZero( pars, "C21" ) * siniota * expPhi;
      }
      else if( freqFactors->data[j] == 2. ){
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
        XLALPrintError("%s: Error... currently unknown frequency factor (%.2lf) for models.\n", __func__, freqFactors->data[j] );
        XLAL_ERROR_VOID( XLAL_EINVAL );
      }

      if ( varyphase || roq ){ /* have to compute the full time domain signal */
        REAL8Vector *sidDayFrac = NULL;
        REAL8 tsv, tsteps;

        REAL8Vector *LUfplus = NULL, *LUfcross = NULL, *LUfx = NULL, *LUfy = NULL, *LUfb = NULL, *LUfl = NULL;

        /* set lookup table parameters */
        tsteps = (REAL8)(*(INT4*)LALInferenceGetVariable( ifo->params, "timeSteps" ));
        tsv = LAL_DAYSID_SI / tsteps;

        LUfplus = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "a_response_tensor" );
        LUfcross = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "b_response_tensor" );

        if ( nonGR ){
          LUfx = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "a_response_vector" );
          LUfy = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "b_response_vector" );
          LUfb = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "a_response_scalar" );
          LUfl = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "b_response_scalar" );
        }

        /* get the sidereal time since the initial data point % sidereal day */
        sidDayFrac = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "siderealDay" );

        length = ifo->times->length;

        #pragma omp parallel for private(T)
        for( i=0; i<length; i++ ){
          REAL8 plus00, plus01, cross00, cross01, plus = 0., cross = 0.;
          REAL8 x00, x01, y00, y01, b00, b01, l00, l01;
          REAL8 timeScaled, timeMin, timeMax;
          REAL8 plusT = 0., crossT = 0., x = 0., y = 0., xT = 0., yT = 0., b = 0., l = 0.;
          INT4 timebinMin, timebinMax;

          /* set the time bin for the lookup table */
          /* sidereal day in secs*/
          T = sidDayFrac->data[i];
          timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
          timeMin = timebinMin*tsv;
          timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
          timeMax = timeMin + tsv;

          /* get values of matrix for linear interpolation */
          plus00 = LUfplus->data[timebinMin];
          plus01 = LUfplus->data[timebinMax];

          cross00 = LUfcross->data[timebinMin];
          cross01 = LUfcross->data[timebinMax];

          /* rescale time for linear interpolation on a unit square */
          timeScaled = (T - timeMin)/(timeMax - timeMin);

          plus = plus00 + (plus01-plus00)*timeScaled;
          cross = cross00 + (cross01-cross00)*timeScaled;

          plusT = plus*c2psi + cross*s2psi;
          crossT = cross*c2psi - plus*s2psi;

          if ( nonGR ){
            x00 = LUfx->data[timebinMin];
            x01 = LUfx->data[timebinMax];
            y00 = LUfy->data[timebinMin];
            y01 = LUfy->data[timebinMax];
            b00 = LUfb->data[timebinMin];
            b01 = LUfb->data[timebinMax];
            l00 = LUfl->data[timebinMin];
            l01 = LUfl->data[timebinMax];

            x = x00 + (x01-x00)*timeScaled;
            y = y00 + (y01-y00)*timeScaled;
            b = b00 + (b01-b00)*timeScaled;
            l = l00 + (l01-l00)*timeScaled;

            xT = x*cpsi + y*spsi;
            yT = y*cpsi - x*spsi;
          }

          /* create the complex signal amplitude model appropriate for the harmonic */
          ifo->compTimeSignal->data->data[i] = ( Cplus * plusT ) + ( Ccross * crossT );

          /* add non-GR components if required */
          if ( nonGR ){ ifo->compTimeSignal->data->data[i] += ( Cx*xT ) + ( Cy*yT ) + Cb*b + Cl*l; }
        }
      }
      else{ /* just have to calculate the values to multiply the pre-summed data */
        /* for tensor-only models (e.g. the default of GR) calculate the two components of
         * the single model value - things multplied by a(t) and things multiplied by b(t)
         * (both these will have real and imaginary components:
         * h(t) = a(t)*(Cplus*cos(2psi) - Ccross*sin(2psi)) +
         *        b(t)*(Cplus*sin(2psi) + Ccross*cos(2psi))
         *
         * For scalar and vector modes also calculate the four addiational components:
         * h(t)^V = aV(t)*(Cx*cos(psi) - Cy*sin(psi)) +
         *          bV(t)*(Cy*sin(psi) + Cy*cos(psi)),
         * h(t)^S = aS(t)*Cb + bS(t)*Cl.
         */

        /* put multiples of a(t) in first value and b(t) in second */
        if ( !nonGR ){
          /* first check that compTimeSignal has been reduced in size to just hold these two values */
          if ( ifo->compTimeSignal->data->length != 2 ){ /* otherwise resize it */
            ifo->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifo->compTimeSignal, 0, 2 );
          }

          ifo->compTimeSignal->data->data[0] = (Cplus*c2psi - Ccross*s2psi);
          ifo->compTimeSignal->data->data[1] = (Cplus*s2psi + Ccross*c2psi);
        }
        else{
          /* first check that compTimeSignal has been reduced in size to just hold these size values */
          if ( ifo->compTimeSignal->data->length != 6 ){ /* otherwise resize it */
            ifo->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifo->compTimeSignal, 0, 6 );
          }

          ifo->compTimeSignal->data->data[0] = (Cplus*c2psi - Ccross*s2psi);
          ifo->compTimeSignal->data->data[1] = (Cplus*s2psi + Ccross*c2psi);
          ifo->compTimeSignal->data->data[2] = (Cx*cpsi - Cy*spsi);
          ifo->compTimeSignal->data->data[3] = (Cx*spsi + Cy*cpsi);
          ifo->compTimeSignal->data->data[4] = Cb;
          ifo->compTimeSignal->data->data[5] = Cl;
        }
      }

      ifo = ifo->next;
    }
  }
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


/**
 * \brief Creates a lookup table of the detector antenna pattern
 *
 * This function creates a lookup table of the tensor, vector and scalar antenna patterns for a given
 * detector orientation and source sky position. For the tensor modes these are the functions given by
 * equations 10-13 in \cite JKS98 , whilst for the vector and scalar modes they are the \f$\psi\f$
 * independent parts of e.g. equations 5-8 of \cite Nishizawa2009 . We remove the \f$\psi\f$ dependent
 * by setting \f$\psi=0\f$.
 *
 * If \c avedt is a value over 60 seconds then the antenna pattern will actually be the mean value from
 * 60 second intervals within that timespan. This accounts for the fact that each data point is actually an
 * averaged value over the given timespan.
 *
 * \param t0 [in] initial GPS time of the data
 * \param detNSource [in] structure containing the detector and source orientations and locations
 * \param timeSteps [in] the number of grid bins to use in time
 * \param avedt [in] average the antenna pattern over this timespan
 * \param aT [in] a vector into which the a(t) Fplus tensor antenna pattern lookup table will be output
 * \param bT [in] a vector into which the b(t) Fcross tensor antenna pattern lookup table will be output
 * \param aV [in] a vector into which the a(t) Fx vector antenna pattern lookup table will be output
 * \param bV [in] a vector into which the b(t) Fy vector antenna pattern lookup table will be output
 * \param aS [in] a vector into which the a(t) Fb scalar antenna pattern lookup table will be output
 * \param bS [in] a vector into which the b(t) Fl scalar antenna pattern lookup table will be output
 */
void response_lookup_table( REAL8 t0, LALDetAndSource detNSource, INT4 timeSteps, REAL8 avedt, REAL8Vector *aT,
                            REAL8Vector *bT, REAL8Vector *aV, REAL8Vector *bV, REAL8Vector *aS,
                            REAL8Vector *bS ){
  LIGOTimeGPS gps;
  REAL8 T = 0, Tstart = 0., Tav = 0;

  REAL8 fplus = 0., fcross = 0., fx = 0., fy = 0., fb = 0., fl = 0.;
  REAL8 tsteps = (REAL8)timeSteps;

  INT4 j = 0, k = 0, nav = 0;

  /* number of points to average */
  if ( avedt == 60. ) { nav = 1; }
  else{ nav = floor(avedt/60.) + 1; }

  /* set the polarisation angle to zero to get the a(t) and b(t) antenna pattern functions */
  detNSource.pSource->orientation = 0.0;

  for( j = 0 ; j < timeSteps ; j++ ){
    aT->data[j] = 0.;
    bT->data[j] = 0.;
    aV->data[j] = 0.;
    bV->data[j] = 0.;
    aS->data[j] = 0.;
    bS->data[j] = 0.;

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

      XLALComputeDetAMResponseExtraModes( &fplus, &fcross, &fb, &fl, &fx, &fy, detNSource.pDetector->response,
                                          detNSource.pSource->equatorialCoords.longitude,
                                          detNSource.pSource->equatorialCoords.latitude,
                                          detNSource.pSource->orientation, XLALGreenwichMeanSiderealTime( &gps ) );

      aT->data[j] += fplus;
      bT->data[j] += fcross;
      aV->data[j] += fx;
      bV->data[j] += fy;
      aS->data[j] += fb;
      bS->data[j] += fl;
    }

    aT->data[j] /= (REAL8)nav;
    bT->data[j] /= (REAL8)nav;
    aV->data[j] /= (REAL8)nav;
    bV->data[j] /= (REAL8)nav;
    aS->data[j] /= (REAL8)nav;
    bS->data[j] /= (REAL8)nav;
  }
}


/*------------------------ END OF MODEL FUNCTIONS ----------------------------*/


/*----------------- FUNCTIONS TO CONVERT BETWEEN PARAMETERS ------------------*/

/**
 * \brief Convert sources parameters into amplitude and phase notation parameters
 *
 * Convert the physical source parameters into the amplitude and phase notation given in Eqns
 * 62-65 of \cite Jones:2015 .
 *
 * Note that \c phi0 is essentially the rotational phase of the pulsar.
 */
void invert_source_params( PulsarParameters *params ){
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
    C22 = 0.5*h0;
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
