#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiralTestGRParams.h>
#include <lal/Window.h>
#include <lal/LALSimBlackHoleRingdownTiger.h>
#define EPS LAL_REAL4_EPS
#define TINY LAL_REAL8_MIN

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// Top-level function to get combination of all modes
// Note: currently not used by LALInference for recovery
// Can be called from SWIG-wrapped lalsimulation for plotting
int XLALSimBlackHoleRingdownTigerAllModes(
    REAL8TimeSeries **hplus,      /**< plus-polarization waveform [returned] */
    REAL8TimeSeries **hcross,     /**< cross-polarization waveform [returned] */
    const LIGOTimeGPS *t0,                /**< start time of ringdown */
    REAL8 phi0,                   /**< initial phase of ringdown (rad) */
    REAL8 deltaT,                 /**< sampling interval (s) */
    REAL8 mass,                   /**< black hole mass (kg) */
    REAL8 a,      /**< black hole dimensionless spin parameter */
    REAL8 eta,         /**< symmetric mass ratio of progenitor */ 
    REAL8 spin1[3],    /**< initial spin for 1st component */
    REAL8 spin2[3],    /**< initial spin for 2nd component */
    REAL8 chiEff,      /**< effective spin parameter for initial spins */
    REAL8 distance,               /**< distance to source (m) */
    REAL8 inclination,            /**< inclination of source's spin axis (rad) */
    LALSimInspiralTestGRParam *nonGRparams  /**< testing GR parameters */
                                         )
{
    SphHarmTimeSeries *qnmodes=NULL;        /**< List containing empty Quasi-Normal Modes  */

    qnmodes = XLALSphHarmTimeSeriesAddMode(qnmodes, NULL, 2, 2);
    qnmodes = XLALSphHarmTimeSeriesAddMode(qnmodes, NULL, 2, 1);
    qnmodes = XLALSphHarmTimeSeriesAddMode(qnmodes, NULL, 3, 3);
    qnmodes = XLALSphHarmTimeSeriesAddMode(qnmodes, NULL, 4, 4);
   
    return XLALSimBlackHoleRingdownTiger(hplus, hcross, qnmodes, t0, phi0, deltaT, mass, a, eta, spin1, spin2, chiEff, distance, inclination, nonGRparams);
}

int XLALSimBlackHoleRingdownTiger(
				  REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
				  REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
				  SphHarmTimeSeries *hlms, /**< Head of linked list of waveform modes */
				  const LIGOTimeGPS *t0,		/**< start time of ringdown */
				  REAL8 phi0,			/**< initial phase of ringdown (rad) */
				  REAL8 deltaT,			/**< sampling interval (s) */
				  REAL8 mass,			/**< black hole mass (kg) */
				  REAL8 a,	/**< black hole dimensionless spin parameter */
				  //	REAL8 fractional_mass_loss,	/**< fraction of mass radiated in this mode */
				  REAL8 eta,         /**< symmetric mass ratio of progenitor */
				  REAL8 spin1[3],    /**< initial spin for 1st component */
				  REAL8 spin2[3],    /**< initial spin for 2nd component */
				  REAL8 chiEff,      /**< effective spin parameter for initial spins */
				  REAL8 distance,		/**< distance to source (m) */
				  REAL8 inclination,		/**< inclination of source's spin axis (rad) */
				  LALSimInspiralTestGRParam *nonGRparams  /**< testing GR parameters */
				  )
{	
  // const INT4 s = -2; /* spin weight for gravitational radiation */
  // COMPLEX16 UNUSED A, UNUSED omega, UNUSED Q, UNUSED tau;
  // COMPLEX16 Yplus, Ycross;
  size_t length;
  size_t j;
  size_t maxmodelength=0;
  size_t thismodelength=0;
  INT4 UNUSED errnum;
  SphHarmTimeSeries *modeList; 
  SphHarmTimeSeries *thisMode;
  //const LALUnit *hunits; // TODO: fill in units for hplus/hcross
  REAL8 dtau = 0.0, dfreq = 0.0 ;
  
  modeList = hlms;
  thisMode = hlms;
  
  char *nonGRParamName = malloc(512*sizeof(char)) ;
  
  /* Go through the list of quasi-normal modes and fill in the structure */
  
  while (thisMode){
    
    sprintf(nonGRParamName,"dtau%d%d",thisMode->l,thisMode->m) ;
    if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
      dtau = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
    sprintf(nonGRParamName,"dfreq%d%d",thisMode->l,thisMode->m) ;
    if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
      dfreq = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
    
    XLALSimBlackHoleRingdownModeTiger(&(thisMode->mode), t0, phi0, deltaT, mass, a, distance, inclination, eta, spin1, spin2, chiEff, thisMode->l, thisMode->m, dfreq, dtau);
    thismodelength = thisMode->mode->data->length;
    maxmodelength = thismodelength > maxmodelength ?  thismodelength : maxmodelength;
    thisMode = thisMode->next;
  }
  // thisMode->mode is a COMPLEX16TimeSeries
    
  /* Create the plus and cross polarization vectors */
  
  length = maxmodelength;
  *hplus = XLALCreateREAL8TimeSeries("hplus", t0, 0.0, deltaT, &lalStrainUnit, length);
  *hcross = XLALCreateREAL8TimeSeries("hcross", t0, 0.0, deltaT, &lalStrainUnit, length);
  
  /* Initialize vectors */
  memset( (*hplus)->data->data, 0, length*sizeof(REAL8) );
  memset( (*hcross)->data->data, 0, length*sizeof(REAL8) );
  
  /* prepare window */
  REAL8Window *window_rd;
  REAL8 start = 0.0; //TODO: starting time of window might be adjusted later on
  REAL8 rise_time = 0.001;
  window_rd = XLALCreatePlanckREAL8Window(length,start,length*deltaT-2.0,1.0/deltaT, rise_time);

  /* Fill in the plus and cross polarization vectors */
  thisMode = modeList;
  COMPLEX16 hpc = 0.0;
  while (thisMode){
    for (j=0; j<thisMode->mode->data->length; j++){
      hpc = (thisMode->mode->data->data[j]);
      (*hplus)->data->data[j] += creal(hpc)*(window_rd->data->data)[j]; // or * prefac(thisMode->l, thisMode->m) ?
      (*hcross)->data->data[j] -= cimag(hpc)*(window_rd->data->data)[j];
    }
    thisMode = thisMode->next;
  }
  
  
  free(nonGRParamName);
  free(window_rd); 
  return 0;
}

/**
 * Computes the waveform for the ringdown of a black hole
 * quasinormal mode (l,m).
 * The amplitudes and QNM frequen are calculated
 * according to Gossan et. al 2012 [arXiv: 1111.5819]
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdownModeTiger(
				      // REAL8TimeSeries **hpluslm,	/**< plus-polarization waveform [returned] */
				      // REAL8TimeSeries **hcrosslm,	/**< cross-polarization waveform [returned] */
				      COMPLEX16TimeSeries **hlmmode,  /**< complex waveform for lm mode */
				      const LIGOTimeGPS *t0,		/**< start time of ringdown */
				      REAL8 phi0,			/**< initial phase of ringdown (rad) */
				      REAL8 deltaT,			/**< sampling interval (s) */
				      REAL8 mass,			/**< black hole mass (kg) */
				      REAL8 a,	/**< black hole dimensionless spin parameter */
				      //	REAL8 fractional_mass_loss,	/**< fraction of mass radiated in this mode */
				      REAL8 distance,		/**< distance to source (m) */
				      REAL8 inclination,		/**< inclination of source's spin axis (rad) */
				      REAL8 eta,         /**< symmetric mass ratio of progenitor */
				      REAL8 UNUSED spin1[3],    /**< initial spin for 1st component */
				      REAL8 UNUSED spin2[3],    /**< initial spin for 2nd component */
				      REAL8 chiEff,      /**< effective spin parameter for initial spins */
				      UINT4 l,				/**< polar mode number */
				      INT4 m,				/**< azimuthal mode number */
				      REAL8 dfreq,         /**< relative shift in the real frequency parameter */
				      REAL8 dtau          /**< relative shift in the damping time parameter */
				      )
{ 
  // COMPLEX16TimeSeries *hlm; // IS THIS NECESSARY?
  // REAL8 spin1[3] = {0.0};
  // REAL8 spin2[3] = {0.0};
  const INT4 UNUSED s = -2; /* spin weight for gravitational radiation */
  // double cosi = cos(inclination);
  COMPLEX16 omega;
  REAL8 A, alphalm, tau, freq;
  COMPLEX16 Yplus, Ycross;
  size_t length;
  size_t j;
  INT4 UNUSED errnum;
  REAL8 mass_sec = mass*LAL_MTSUN_SI/LAL_MSUN_SI ;
  REAL8 dist_sec = distance/LAL_C_SI ;
  char name[256];     // TODO: replace with name[MAXNAMELENGTH]
  // const LALUnit *hunits; // TODO: fill in units
    
  sprintf(name, "h%u%d", l, m);

  alphalm = XLALSimRingdownQNMAmplitudes(l, m, eta, chiEff);
  A = alphalm*mass_sec/dist_sec;
    
  Yplus = XLALSimSphericalHarmonicPlus(l, m, inclination);
  Ycross = XLALSimSphericalHarmonicCross(l, m, inclination);

  omega = XLALSimRingdownFitOmega(l, m, 0, a);
  tau = XLALQNMTauOfOmega(omega, mass_sec)*(1 + dtau);
  freq = XLALQNMFreqOfOmega(omega, mass_sec)*(1 + dfreq);
  
  length = (size_t) ceil( -tau * log(1.0e-16) / deltaT ) ;

  /* allocate memory in hlm
   * freed in LALInferenceTemplateXLALSimBlackHoleRingdown using XLALDestroySphHarmTimeSeries(qnmodes)
   */
  *hlmmode = XLALCreateCOMPLEX16TimeSeries(name, t0, 0.0, deltaT, &lalStrainUnit, length);
  memset( (*hlmmode)->data->data, 0, length*sizeof(COMPLEX16) ) ;

  REAL8 t = 0.0 ;
	if (A > 0.0){
	  for (j=0; j<length; j++){
      t = ((REAL8) j)*deltaT ;
      (*hlmmode)->data->data[j] = A*Yplus*exp(-t/tau)*cos(t*freq - m*phi0);
      (*hlmmode)->data->data[j] -= I*A*Ycross*exp(-t/tau)*sin(t*freq - m*phi0);
	  }
  }
    
    // hlmmode = &hlm; 

/*
	if (XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	XLAL_TRY(sphwf1 = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(mu, a, l, m, s, A, omega), errnum);
	XLAL_TRY(sphwf2 = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(-mu, a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR(XLAL_EFUNC);
	omega *= 0.5; // convert from Leaver's convention 2M = 1 
	
	// compute length of waveform to compute 
	length = ceil(log(LAL_REAL8_EPS) * LAL_G_SI * mass / (pow(LAL_C_SI, 3.0) * cimag(omega) * deltaT));
	if (length < 1)
		XLAL_ERROR(XLAL_EBADLEN);

	// compute the amplitude factors for the +m and -m modes 
	A1 = A2 = -4.0 * (LAL_G_SI*mass/(pow(LAL_C_SI, 2.0)*distance))
		* sqrt(-0.5*cimag(omega)*fractional_mass_loss)/cabs(omega)
 		* cexp(I*m*phi0);
	A1 *= sphwf1;
	A2 = conj(A2*sphwf2);

	omega_dt = pow(LAL_C_SI, 3.0)*omega*deltaT/(LAL_G_SI*mass);

	*hplus = XLALCreateREAL8TimeSeries("H_PLUS", t0, 0.0, deltaT, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("H_CROSS", t0, 0.0, deltaT, &lalStrainUnit, length);
	if (hplus == NULL || hcross == NULL)
		XLAL_ERROR(XLAL_EFUNC);


	// compute the waveforms 
	for (j = 0; j < length; ++j) {
		COMPLEX16 h;
		h = A1*cexp(-I*omega_dt*j) + A2*cexp(I*conj(omega_dt)*j);
		(*hplus)->data->data[j] = creal(h);
		(*hcross)->data->data[j] = -cimag(h);
	}
*/
	return 0;
}

REAL8 XLALSimSphericalHarmonicPlus(UINT4 l, INT4 m, REAL8 iota){
	COMPLEX16 Yplus = 0.0;
	Yplus = XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, m) + (1.0 - 2.0*(l % 2))*XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, -m);
	return Yplus;
}

REAL8 XLALSimSphericalHarmonicCross(UINT4 l, INT4 m, REAL8 iota){
	COMPLEX16 Ycross = 0.0;
	Ycross = XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, m) - (1.0 - 2.0*(l % 2))*XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, -m);
	return Ycross;
}


/**
 * Calculates the amplitudes for the QNM l, m, n=0 for
 * a given symmetric mass-ratio eta of the initial binary.
 * Based on an interpolation for NON-SPINNING binary black
 * holes, derived in Gossan et al (2012) [arXiv: 1111.5819]
 * 
 * TO DO: rewrite this properly, add initial (effective) spin dependence
 **/
REAL8 XLALSimRingdownQNMAmplitudes(INT4 l, INT4 m, REAL8 eta, REAL8 chiEff){
	REAL8 A = 0.8639*eta;
	if (l==2 && m==2){ A *= 1.0; }
    // else if (l==2 && abs(m)==1){ A *= 0.52*pow(1.0 - 4.*eta, 0.71); }	
	else if (l==2 && abs(m)==1){ A *= 0.43*(sqrt(1.0 - 4.*eta) - chiEff); }	
	else if (l==3 && abs(m)==3){ A *= 0.44*pow(1.0 - 4.*eta, 0.45); }
	else if (l==3 && abs(m)==2){ A *= 3.69*(eta - 0.2)*(eta - 0.2) + 0.053; }
	else if (l==4 && abs(m)==4){ A *= 5.41*(eta - 0.22)*(eta - 0.22) + 0.04; }
	else A = 0.0;
	return A;
	}
	
/**
 * Computes the complex dimensionless eigen-frequency for the
 * quasinormal mode (l,m), given a dimensionless spin a \in [-1,1].
 * This is based on the interpolation formula from Berti et al 2006
 * [arXiv: gr-qc/0512160].
 * The dimensionful frequency is Re(omega)/M.
 * The dimensionful damping time is M/Im(omega).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
COMPLEX16 XLALSimRingdownFitOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a){
	COMPLEX16 omega=0.0;
	if (fabs(a) > 1.){
		fprintf(stderr, "ERROR: Dimensionless spin magnitude larger than 1! Aborting...\n");
		exit(-1);
	}
	switch (l){
	case 2:
		switch (abs(m)){
		case 2:
			omega = (1.5251 - 1.1568*pow(1.e0-a,0.1292));
			omega += I*omega/(2.*(0.7000 + 1.4187*pow(1.e0-a,-0.4990)));
			break;
		case 1:
			omega = (0.6000 - 0.2339*pow(1.e0-a,0.4175));
			omega += I*omega/(2.*(-0.3000 + 2.3561*pow(1.e0-a,-0.2277)));
			break;
		default:
			break;
		}
		break;
	case 3:
		switch (abs(m)){
		case 3:
			omega = (1.8956 - 1.3043*pow(1.e0-a,0.1818));
			omega += I*omega/(2.*(0.9000 + 2.3430*pow(1.e0-a,-0.4810)));
			break;
		case 2:
			omega = (1.1481 - 0.5552*pow(1.e0-a,0.3002));
			omega += I*omega/(2.*(0.8313 + 2.3773*pow(1.e0-a,-0.3655)));
			break;
		default:
			break;
		}
		break;
	case 4:
		switch (abs(m)){
		case 4:
			omega = (2.3000 - 1.5056*pow(1.e0-a,0.2244));
			omega += I*omega/(2.*(1.1929 + 3.1191*pow(1.e0-a,-0.4825)));
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
	if (omega == 0.0){
		fprintf(stderr, "ERROR: Invalid mode l = %u, m = %d, n = %u\n", l, m, n);
		exit(-1);
	}
	return omega;
}

/** 
 * Frequency in rad/sec given dimensionless complex frequency and M
 **/
REAL8 XLALQNMFreqOfOmega(COMPLEX16 omega, REAL8 mtot){
	return (creal(omega)/mtot);
}

/** 
 * Damping time in sec given dimensionless complex frequency and M
 **/
REAL8 XLALQNMTauOfOmega(COMPLEX16 omega, REAL8 mtot){
	REAL8 tau = 0;
	tau = mtot/cimag(omega);
	return tau;
}

