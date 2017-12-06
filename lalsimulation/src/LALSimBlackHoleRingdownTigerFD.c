#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h> 
#include <lal/LALConstants.h> 
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiral.h>
#include <lal/XLALError.h> 
#include <lal/Date.h>
#include <string.h>

#include <lal/LALSimInspiralTestGRParams.h> 
#include <lal/LALSimBlackHoleRingdownTiger.h>
#define EPS LAL_REAL4_EPS
#define TINY LAL_REAL8_MIN 
#ifdef __GNUC_
#define UNUSED __attribute__ ((unused)) 
#else 
#define UNUSED 
#endif
/*High level function to combine all the modes*/
int XLALSimBlackHoleRingdownTigerFD(
		COMPLEX16FrequencySeries **hptilde,    /**< FD plus polarization */
		COMPLEX16FrequencySeries **hctilde,    /**< FD cross polarization */
		REAL8 phi0,                            /**<initial phase of ringdown*/
		REAL8 deltaF,                         /**<sampling interval (Hz)*/
		REAL8 fEnd,
		REAL8 fStart,
		REAL8 mass,                           /**<blackhole mass*/
		REAL8 a,                              /**<blackhole dimensionless spin parameter*/
		REAL8 eta,                            /**<symmetric mass ratio*/
		REAL8 chiEff,
		REAL8 distance,                       /**distance(m)*/
		REAL8 inclination,                    /**< inclination of source's spin axis (rad)*/
		LALSimInspiralTestGRParam *nonGRparams  /**< testing GR parameters */
		)
{
	size_t j, n, jStart;
	REAL8 f_max;
	REAL8 f;
	REAL8 dtau22 = 0.0, dfreq22 = 0.0;
	REAL8 dtau21 = 0.0, dfreq21 = 0.0;
	REAL8 dtau33 = 0.0, dfreq33 = 0.0;
	REAL8 dtau32 = 0.0, dfreq32 = 0.0;
	REAL8 dtau44 = 0.0, dfreq44 = 0.0;
	COMPLEX16 *data_p = NULL;
	COMPLEX16 *data_c = NULL;
	
	char *nonGRParamName = malloc(512*sizeof(char)) ; 

	sprintf(nonGRParamName,"dtau22") ;
    	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dtau22 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq22") ;
	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dfreq22 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
	sprintf(nonGRParamName,"dtau21") ;
    	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dtau21 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq21") ;
	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dfreq21 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
	sprintf(nonGRParamName,"dtau33") ;
    	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dtau33 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq33") ;
	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dfreq33 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
	sprintf(nonGRParamName,"dtau32") ;
    	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dtau32 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq32") ;
	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dfreq32 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
	sprintf(nonGRParamName,"dtau44") ;
    	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dtau44 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq44") ;
	if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
	dfreq44 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;

	LIGOTimeGPS tC = {0, 0};
	 /* allocate htilde */
	if ( fEnd == 0. ) // End at 4096(default srate)/2
		f_max = 2048.0;
	else // End at user-specified freq.
		f_max = fEnd;
	n = (size_t) (f_max / deltaF + 1);
	jStart = (size_t) ceil(fStart / deltaF);
	XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */

	*hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
	if (!(*hptilde)) XLAL_ERROR(XLAL_EFUNC);
	memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));
	XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
	*hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
	if (!(*hctilde)) XLAL_ERROR(XLAL_EFUNC);
	memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
	XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);
	COMPLEX16FrequencySeries *hp22, *hc22, *hp21, *hc21, *hp33, *hc33, *hp32, *hc32, *hp44, *hc44; /* *hp2m2, *hc2m2, *hp2m1, *hc2m1, *hp3m3, *hc3m3, *hp3m2, *hc3m2, *hp4m4, *hc4m4;*/

	/*importing waveforms of different mode*/
	XLALSimBlackHoleRingdownModeTigerFD(&hp22, &hc22, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 2, 2, chiEff, dfreq22, dtau22);
	XLALSimBlackHoleRingdownModeTigerFD(&hp21, &hc21, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 2, 1, chiEff, dfreq21, dtau21);
	XLALSimBlackHoleRingdownModeTigerFD(&hp33, &hc33, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 3, 3, chiEff, dfreq33, dtau33);
	XLALSimBlackHoleRingdownModeTigerFD(&hp32, &hc32, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 3, 2, chiEff, dfreq32, dtau32);
	XLALSimBlackHoleRingdownModeTigerFD(&hp44, &hc44, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 4, 4, chiEff, dfreq44, dtau44);

/*	XLALSimBlackHoleRingdownModeTigerFD(&hp2m2, &hc2m2, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 2, -2, chiEff, dfreq22, dtau22);
	XLALSimBlackHoleRingdownModeTigerFD(&hp2m1, &hc2m1, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 2, -1, chiEff, dfreq21, dtau21);
	XLALSimBlackHoleRingdownModeTigerFD(&hp3m3, &hc3m3, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 3, -3, chiEff, dfreq33, dtau33);
	XLALSimBlackHoleRingdownModeTigerFD(&hp3m2, &hc3m2, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 3, -2, chiEff, dfreq32, dtau32);
	XLALSimBlackHoleRingdownModeTigerFD(&hp4m4, &hc4m4, deltaF, fStart, fEnd, phi0, mass, a, distance, inclination, eta, 4, -4, chiEff, dfreq44, dtau44);	*/

	data_p = (*hptilde)->data->data;
	data_c = (*hctilde)->data->data;	

	f = jStart*deltaF;
	for ( j=jStart;j<n;j++) { 

	COMPLEX16 hp22a = 0.0;
        COMPLEX16 hc22a = 0.0;
	COMPLEX16 hp21a = 0.0;
	COMPLEX16 hc21a = 0.0;
	COMPLEX16 hp33a = 0.0;
	COMPLEX16 hc33a = 0.0;
	COMPLEX16 hp32a = 0.0;
	COMPLEX16 hc32a = 0.0;  
	COMPLEX16 hp44a = 0.0;
	COMPLEX16 hc44a = 0.0;

/*	COMPLEX16 hp2m2a = 0.0;
        COMPLEX16 hc2m2a = 0.0;
	COMPLEX16 hp2m1a = 0.0;
	COMPLEX16 hc2m1a = 0.0;
	COMPLEX16 hp3m3a = 0.0;
	COMPLEX16 hc3m3a = 0.0;
	COMPLEX16 hp3m2a = 0.0;
	COMPLEX16 hc3m2a = 0.0;  
	COMPLEX16 hp4m4a = 0.0;
	COMPLEX16 hc4m4a = 0.0;*/

	hp32a = (hp32->data->data[j]);
	hc32a = (hc32->data->data[j]); 
	hp22a = (hp22->data->data[j]);
	hc22a = (hc22->data->data[j]); 
	hp21a = (hp21->data->data[j]); 
	hc21a = (hc21->data->data[j]);	       
	hp33a = (hp33->data->data[j]); 
	hc33a = (hc33->data->data[j]); 
	hp44a = (hp44->data->data[j]);
	hc44a = (hc44->data->data[j]);

/*	hp3m2a = (hp3m2->data->data[j]);
	hc3m2a = (hc3m2->data->data[j]); 
	hp2m2a = (hp2m2->data->data[j]);
	hc2m2a = (hc2m2->data->data[j]); 
	hp2m1a = (hp2m1->data->data[j]); 
	hc2m1a = (hc2m1->data->data[j]);	       
	hp3m3a = (hp3m3->data->data[j]); 
	hc3m3a = (hc3m3->data->data[j]); 
	hp4m4a = (hp4m4->data->data[j]);
	hc4m4a = (hc4m4->data->data[j]); */

/*	data_p[j] = hp22a+hp2m2a+hp33a+hp3m3a+hp21a+hp2m1a+hp32a+hp3m2a+hp44a+hp4m4a;
	data_c[j] = hc22a+hc2m2a+hc33a+hc3m3a+hc21a+hc2m1a+hc32a+hc3m2a+hc44a+hc4m4a;*/

	data_p[j] = hp22a+hp21a+hp32a+hp33a+hp44a;
	data_c[j] = hc22a+hc21a+hc32a+hc33a+hc44a;
	
	f+=deltaF;

	}
        
        if(hp32) XLALDestroyCOMPLEX16FrequencySeries(hp32);
        if(hc32) XLALDestroyCOMPLEX16FrequencySeries(hc32);
        if(hp22) XLALDestroyCOMPLEX16FrequencySeries(hp22);
        if(hc22) XLALDestroyCOMPLEX16FrequencySeries(hc22);
        if(hp21) XLALDestroyCOMPLEX16FrequencySeries(hp21);
        if(hc21) XLALDestroyCOMPLEX16FrequencySeries(hc21);
        if(hp33) XLALDestroyCOMPLEX16FrequencySeries(hp33);
        if(hc33) XLALDestroyCOMPLEX16FrequencySeries(hc33);
        if(hp44) XLALDestroyCOMPLEX16FrequencySeries(hp44);
        if(hc44) XLALDestroyCOMPLEX16FrequencySeries(hc44);

/*	if(hp3m2) XLALDestroyCOMPLEX16FrequencySeries(hp3m2);
        if(hc3m2) XLALDestroyCOMPLEX16FrequencySeries(hc3m2);
        if(hp2m2) XLALDestroyCOMPLEX16FrequencySeries(hp2m2);
        if(hc2m2) XLALDestroyCOMPLEX16FrequencySeries(hc2m2);
        if(hp2m1) XLALDestroyCOMPLEX16FrequencySeries(hp2m1);
        if(hc2m1) XLALDestroyCOMPLEX16FrequencySeries(hc2m1);
        if(hp3m3) XLALDestroyCOMPLEX16FrequencySeries(hp3m3);
        if(hc3m3) XLALDestroyCOMPLEX16FrequencySeries(hc3m3);
        if(hp4m4) XLALDestroyCOMPLEX16FrequencySeries(hp4m4);
        if(hc4m4) XLALDestroyCOMPLEX16FrequencySeries(hc4m4); */

	return XLAL_SUCCESS;
}
/*Create the waveform frequency series (+ and x polarization) with specific l and m. */
int XLALSimBlackHoleRingdownModeTigerFD(
		COMPLEX16FrequencySeries **hptilde_lm,    /**< FD plus polarization */
		COMPLEX16FrequencySeries **hctilde_lm,    /**< FD cross polarization */
		const REAL8 deltaF,                    /**< Frequency resolution */
		const REAL8 fStart,                    /**< Start GW frequency (Hz) */
		const REAL8 fEnd,                      /**< Highest GW frequency (Hz)*/
		REAL8 phi0,                       /**< initial phase of ringdown */
		REAL8 mass,                       /**< black hole mass (kg) */
		REAL8 a,  /**< black hole dimensionless spin parameter */
		REAL8 distance,           /**< distance to source (m) */
		REAL8 inclination,                /**< inclination of source's spin axis (rad) */ 
		REAL8 eta,         /**< symmetric mass ratio of progenitor */ 
		UINT4 l,
		INT4 m,
		REAL8 chiEff,      /**< effective spin parameter for initial spins */ 
		REAL8 dfreq,         /**< relative shift in thereal frequency parameter */
		REAL8 dtau          /**< relative shift in the damping time parameter */
		)
{
	REAL8 A_lm, alphalm, tau, freq;
	COMPLEX16 omegalm;
	REAL8 mass_sec = mass*LAL_MTSUN_SI/LAL_MSUN_SI;
	REAL8 dist_sec = distance/LAL_C_SI;
	REAL8 shift;
	COMPLEX16 Ypluslm, Ycrosslm;
	COMPLEX16 A_plus, A_cross;
	size_t j, n, jStart;
	REAL8 f_max;
	REAL8 f;
	COMPLEX16 *data_p = NULL;
	COMPLEX16 *data_c = NULL;
	
	 /* Perform some initial checks */ 
	if (mass <= 0) XLAL_ERROR(XLAL_EDOM); 
	if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
	if (distance <= 0) XLAL_ERROR(XLAL_EDOM);	

	/* The following functions are difeined in LALSimBlackHoleRingdownTiger.c*/

	alphalm = XLALSimRingdownQNMAmplitudes(l, m, eta, chiEff);
	
	A_lm = alphalm*mass_sec/dist_sec;
	
	Ypluslm = XLALSimSphericalHarmonicPlus(l, m, inclination); 
	Ycrosslm = XLALSimSphericalHarmonicCross(l, m, inclination);
	
	A_plus = A_lm*Ypluslm ;
	A_cross = A_lm*Ycrosslm ;

	omegalm = XLALSimRingdownFitOmega(l, m, 0, a);

	tau = XLALQNMTauOfOmega(omegalm, mass_sec)*(1 + dtau);
	freq = XLALQNMFreqOfOmega(omegalm, mass_sec)*(1 + dfreq);

	/* allocate htilde_p and htilde_c*/
	if ( fEnd == 0. )
	f_max = 2048.0;
	else
	f_max = fEnd;

	LIGOTimeGPS tC = {0,0};

	n = (size_t)(f_max/deltaF + 1);
	shift = LAL_TWOPI*(tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);
	XLALGPSAdd(&tC, -1 / deltaF);  /*coalesce at t=0*/
	
	*hptilde_lm = XLALCreateCOMPLEX16FrequencySeries("hptilde_lm: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
	if (!(*hptilde_lm)) XLAL_ERROR(XLAL_EFUNC);
	memset((*hptilde_lm)->data->data, 0, n * sizeof(COMPLEX16));
	XLALUnitMultiply(&(*hptilde_lm)->sampleUnits, &(*hptilde_lm)->sampleUnits, &lalSecondUnit);
	*hctilde_lm = XLALCreateCOMPLEX16FrequencySeries("hctilde_lm: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
	if (!(*hctilde_lm)) XLAL_ERROR(XLAL_EFUNC);
	memset((*hctilde_lm)->data->data, 0, n * sizeof(COMPLEX16)); 
	XLALUnitMultiply(&(*hctilde_lm)->sampleUnits, &(*hctilde_lm)->sampleUnits, &lalSecondUnit);
	jStart = (size_t) ceil(fStart / deltaF);
	f = jStart*deltaF;
	data_p = (*hptilde_lm)->data->data;
	data_c = (*hctilde_lm)->data->data; 

	for ( j=jStart;j<n;j++) {

	data_p[j] = (A_plus/2.)*(cos(f*shift) - I*sin(f*shift))*((tau-I*tau*tau*(freq+LAL_TWOPI*f))*(cos(m*phi0)+I*sin(m*phi0))/(1.+ tau*tau*(freq+LAL_TWOPI*f)*(freq+LAL_TWOPI*f)) + (tau-I*tau*tau*(-freq+LAL_TWOPI*f))*(cos(m*phi0)-I*sin(m*phi0))/(1.+ tau*tau*(-freq+LAL_TWOPI*f)*(-freq+LAL_TWOPI*f)));
	data_c[j] = (-I*A_cross/2.)*(cos(f*shift) - I*sin(f*shift))*((tau-I*tau*tau*(freq+LAL_TWOPI*f))*(cos(m*phi0)+I*sin(m*phi0))/(1.+ tau*tau*(freq+LAL_TWOPI*f)*(freq+LAL_TWOPI*f)) - (tau-I*tau*tau*(-freq+LAL_TWOPI*f))*(cos(m*phi0)-I*sin(m*phi0))/(1.+ tau*tau*(-freq+LAL_TWOPI*f)*(-freq+LAL_TWOPI*f)));


	f+=deltaF;

	}
	
	return 0;
} 
