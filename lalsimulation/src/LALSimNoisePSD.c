/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <math.h>
#include <stdio.h>

#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/FrequencySeries.h>
#include <lal/LALSimReadData.h>
#include <lal/LALSimNoise.h>

// Values for iLIGO

#define LAL_ILIGO_ARMLENGTH_SI 3995.0
#define LAL_ILIGO_LASER_POWER_BS_SI 250.0
#define LAL_ILIGO_LASER_WAVELENGTH_SI 1.064e-6
#define LAL_ILIGO_FINESSE 220.0
#define LAL_ILIGO_MIRROR_MASS_SI 11.0
#define LAL_ILIGO_TEMPERATURE_SI 290.0
#define LAL_ILIGO_THERMAL_STACK_FREQ_SI 10.0
#define LAL_ILIGO_THERMAL_SUSP_FREQ_SI 0.76
#define LAL_ILIGO_THERMAL_SUSP_QUAL 1e6
#define LAL_ILIGO_THERMAL_COAT_FREQ_SI 1e4
#define LAL_ILIGO_THERMAL_COAT_QUAL 1e6

// Values for aLIGO

// these ones seem vaguely well motivated
#define LAL_ALIGO_ARMLENGTH_SI 3995.0
#define LAL_ALIGO_LASER_POWER_LOW_SI 25.0
#define LAL_ALIGO_LASER_POWER_HIGH_SI 125.0
#define LAL_ALIGO_LASER_WAVELENGTH_SI 1.064e-6
#define LAL_ALIGO_MIRROR_MASS_SI 40.0
#define LAL_ALIGO_MIRROR_LOSS 37.5e-6
#define LAL_ALIGO_BS_LOSS 0.002
#define LAL_ALIGO_ITM_TRANSMITTANCE 0.014
#define LAL_ALIGO_PRM_TRANSMITTANCE 0.027
#define LAL_ALIGO_SRM_TRANSMITTANCE 0.2
// these ones are not physically motivated
// but seem to be phenomenologically rightish
#define LAL_ALIGO_TEMPERATURE_SI 290.0
#define LAL_ALIGO_THERMAL_SUSP_FREQ_SI 9.0
#define LAL_ALIGO_THERMAL_SUSP_QUAL 6e10
#define LAL_ALIGO_THERMAL_COAT_FREQ_SI 1e4
#define LAL_ALIGO_THERMAL_COAT_QUAL 6e6




/**
 * Provides a rather ad-hoc estimate of the seismic noise power spectral density
 * at a given frequency.
 *
 * This is a crude estimate based on characteristic frequencies for the
 * pendulum and the stack.  What is computed is
 * \f[
 * S_h(f) = L^{-2} S_g(f) (f_{\mathrm{pend}}/f)^4
 * (f_{\mathrm{stack}}/f)^{4n_{\mathrm{stack}}}
 * \f]
 * where
 * \f[
 * S_g(f) = 10^{-18}\,\mathrm{m}^2\,\mathrm{Hz}^{-1}\times(10\,\mathrm{Hz}/f)^4
 * \f]
 * is the displacement power spectrum of ground motion.
 *
 * Warning: the transfer function is only correct at frequencies above the
 * specified characteristic pendulum and stack frequencies.
 */
double XLALSimNoisePSDSeismic(
	double f,		/**< frequency (Hz) */
	double L,		/**< arm length (m) */
	double f_pend,		/**< characteristic frequency of pendulum suspension (Hz) */
	double f_stack,		/**< characteristic frequency of isolation stack (Hz) */
	double n_stack		/**< number of layers of stack */
)
{
	double A_pend   = pow(f_pend / f, 2);
	double A_stack  = pow(f_stack / f, 2 * n_stack);
	double S_ground = 1e-18; /* m^2 / Hz */

	if (f > 10.0)
		S_ground *= pow(10.0 / f, 4);

	return S_ground * pow(A_pend * A_stack / L, 2);
}


/**
 * Provides a rather ad-hoc estimate of the suspension thermal noise power
 * spectral density at a given frequency.
 *
 * This is a crude estimate based on the characteristic frequency of the
 * pendulum suspension and its quality factor (= 1 / loss angle).  What is
 * computed is
 * \f[
 * S_h(f) = L^{-2} \frac{2 k T}{\pi^3 f_0^3 M Q} \left( \frac{f_0}{f} \right)^5.
 * \f]
 *
 * Warning: this only describes the broadband noise at frequencies above the
 * pendulum frequency; it does not have the correct noise near the resonances.
 */
double XLALSimNoisePSDSuspTherm(
	double f,		/**< frequency (Hz) */
	double L,		/**< arm length (m) */
	double M,		/**< mirror mass (kg) */
	double T,		/**< temperature (K) */
	double f0,		/**< pendulum frequency */
	double Q		/**< pendulum quality */
)
{
	double fac = 2.0 * LAL_K_SI * T / (L * L * M * Q * pow(LAL_PI * f0, 3.0));
	return fac * pow(f0 / f, 5.0);
}


/**
 * Provides a rather ad-hoc estimate of the mirror thermal noise power spectral
 * density at a given frequency.
 *
 * This is a crude estimate based on the characteristic frequency of the
 * mirror/coating internal modes and their quality factor (= 1 / loss angle).
 * What is computed is
 * \f[
 * S_h(f) = L^{-2} \frac{2 k T}{\pi^3 f_0^3 M Q} \frac{f_0}{f}
 * \f]
 *
 * Warning: this only describes the broadband noise at frequencies below the
 * resonance frequency; it does not have the correct noise near the resonances.
 */
double XLALSimNoisePSDMirrorTherm(
	double f,		/**< frequency (Hz) */
	double L,		/**< arm length (m) */
	double M,		/**< mirror mass (kg) */
	double T,		/**< average per mirror power loss */
	double f0,		/**< average per mirror power loss */
	double Q		/**< average per mirror power loss */
)
{
	double fac = 2.0 * LAL_K_SI * T / (L * L * M * Q * pow(LAL_PI * f0, 3.0));
	return fac * (f0 / f);
}


/**
 * Computes the shot noise in strain-equivalent units using a conventional
 * model appropriate to initial interferometric detectors.
 *
 * Uses the formula for shot noise from
 *
 */
double XLALSimNoisePSDShot(
	double f,		/**< frequency (Hz) */
	double P_BS,		/**< laser power on beamsplitter (W) */
	double lambda,		/**< laser wavelength (m) */
	double L,		/**< arm length (m) */
	double finesse,		/**< arm cavity finesse */
	double eta		/**< effective quantum efficiency of photodiode */
)
{
	double tau_s   = L*finesse/(LAL_PI*LAL_C_SI); // cavity storage time
	double f_pole  = 1.0 / (4.0*LAL_PI*tau_s);
	double S_DC    = ((LAL_PI*LAL_HBAR_SI*lambda*f_pole*f_pole)/(LAL_C_SI*eta*P_BS)); // DC limit of shot noise
	COMPLEX16 C_FAC; // normalized sensing function
	C_FAC  = cexp(2.0*LAL_PI*I*f*L/LAL_C_SI) * sinh(2.0*LAL_PI*f_pole*L/LAL_C_SI) / csinh((2.0*LAL_PI*f_pole*L/LAL_C_SI)*(1.0 + I*(f/f_pole)));
	return S_DC / pow(cabs(C_FAC), 2);
}

/**
 * Computes the quantum noise (shot noise and radiation pressure noise)
 * according to Buonanno and Chen, Phys. Rev. D 64 0402006 (2001).
 *
 * This code is adapted from the GWINC matlab function shotrad.m which includes
 * updated losses by Kirk McKenzie.
 *
 * For simplicity, only losses from the mirrors are included.  Losses from
 * coupling and from the SRC are ignored.  (These could be included as
 * effective losses in A_BS if needed.) A fixed photdiode quantum efficiency of
 * eta = 0.9 is used.
 *
 * Note: this code is adapted from GWINC.
 */
double XLALSimNoisePSDQuantum(
	double f,		/**< frequency (Hz) */
	double I0,		/**< laser power (W) */
	double lambda,		/**< laser wavelength (m) */
	double L,		/**< arm length (m) */
	double M,		/**< mirror mass (kg) */
	double A,		/**< average per mirror power loss */
	double A_BS,		/**< power loss at beam splitter */
	double T_ITM,		/**< transmittance of ITM */
	double T_PRM,		/**< transmittance of PRM */
	double T_SRM,		/**< transmittance of SRM */
	double ds,		/**< detuning phase (rad) */
	double zeta,		/**< demod/detection/homodyne phase */
	double eta		/**< quantum efficiency of photodiode */
)
{
	/* This code is adapted from GWINC */
	double Omega   = 2.0*LAL_PI*f;		// [BC, table 1] Signal angular frequency [rad/s]
	double omega_0 = 2.0*LAL_PI*LAL_C_SI/lambda;	// [BC, table 1] Laser angular frequency [rad/s]
	double lambda_SR = A_BS;
	double lambda_PD = 1.0 - eta;
	double tau     = sqrt(T_SRM);		// SRM Transmittance [amplitude]
	double rho     = sqrt(1.0 - tau*tau);	// SRM Reflectivity [amplitude]
	double phi     = (LAL_PI-ds)/2.0;		// [BC, between 2.14 & 2.15] SR Detuning
	double lambda_arm = A*2.0;		// [BC, after 5.2] Round Trip loss in arm [Power]
	double gamma_ac = T_ITM*LAL_C_SI/(4.0*L);	// [KLMTV-PRD2001] Arm cavity half bandwidth [1/s]
	double epsilon = lambda_arm/(2.0*gamma_ac*L/LAL_C_SI);	// [BC, after 5.2] Loss coefficent for arm cavity
	double r1      = sqrt(1.0 - T_ITM);	// ITM Reflectivity
	double rarm    = r1 - T_ITM * sqrt(1.0 - 2.0*A) / (1.0 - r1*sqrt(1.0 - 2.0*A)); // Arm Cavity Reflectivity
	double G_PRC   = T_PRM/pow(1.0 + sqrt(1.0 - T_PRM)*rarm*sqrt(1.0 - A_BS), 2.0); // Power Recycling Gain
	double I_0     = G_PRC*I0;		// [BC, table 1] Power at BS [W]
	double I_SQL   = (M*L*L*pow(gamma_ac,4.0))/(4.0*omega_0);	// [BC, 2.14] Power to reach free mass SQL [W]
	double Kappa  = 2.0*((I_0/I_SQL)*pow(gamma_ac,4.0))/(Omega*Omega*(gamma_ac*gamma_ac+Omega*Omega));   // [BC, 2.13] Effective Radiation Pressure Coupling
	double beta    = atan(Omega/gamma_ac);	// [BnC, after 2.11] Phase shift of GW SB in arm
	double h_SQL   = sqrt(8.0*LAL_HBAR_SI/(M*pow(Omega*L,2.0)));	// [BnC, 2.12] SQL Strain

	// Coefficients [BC, Equations 5.8 to 5.12]
	COMPLEX16 C11_L  = sqrt(1.0-lambda_PD) * ( (1.0+rho*rho) * ( cos(2.0*phi) + Kappa/2.0 * sin(2.0*phi) ) - 2.0*rho*cos(2.0*beta) - 1.0/4.0*epsilon * ( -2.0 * (1.0+cexp(2.0*I*beta))*(1.0+cexp(2.0*I*beta)) * rho + 4.0 * (1.0+rho*rho) * pow(cos(beta),2.0)*cos(2.0*phi) + ( 3.0+cexp(I*2.0*beta) ) * Kappa * (1.0+rho*rho) * sin(2.0*phi) ) + lambda_SR * ( cexp(2.0*I*beta)*rho-1.0/2.0 * (1.0+rho*rho) * ( cos(2.0*phi)+Kappa/2.0 * sin(2.0*phi) ) ) );
	COMPLEX16 C22_L  = C11_L;
	COMPLEX16 C12_L  = sqrt(1.0-lambda_PD) * tau*tau * ( - ( sin(2.0*phi) + Kappa*pow(sin(phi),2.0) )+ 1.0/2.0*epsilon*sin(phi) * ( (3.0+cexp(2.0*I*beta)) * Kappa * sin(phi) + 4.0*pow(cos(beta),2.0) * cos(phi)) + 1.0/2.0*lambda_SR * ( sin(2.0*phi)+Kappa*pow(sin(phi),2.0)) );
	COMPLEX16 C21_L  = sqrt(1.0-lambda_PD) * tau*tau * ( (sin(2.0*phi)-Kappa*pow(cos(phi),2.0) ) + 1.0/2.0*epsilon*cos(phi) * ( (3.0+cexp(2.0*I*beta) )*Kappa*sin(phi) - 4.0*pow(cos(beta),2.0)*sin(phi) ) + 1.0/2.0*lambda_SR * ( -sin(2.0*phi) + Kappa*pow(cos(phi),2.0)) );

	COMPLEX16 D1_L   = sqrt(1.0-lambda_PD) * ( - (1.0+rho*cexp(2.0*I*beta) ) * sin(phi) + 1.0/4.0*epsilon * ( 3.0+rho+2.0*rho*cexp(4.0*I*beta) + cexp(2.0*I*beta)*(1.0+5.0*rho) ) * sin(phi)+ 1.0/2.0*lambda_SR * cexp(2.0*I*beta) * rho * sin(phi) );
	COMPLEX16 D2_L   = sqrt(1.0-lambda_PD) * ( - (-1.0+rho*cexp(2.0*I*beta) ) * cos(phi) + 1.0/4.0*epsilon * ( -3.0+rho+2.0*rho*cexp(4.0*I*beta) + cexp(2.0*I*beta) * (-1.0+5.0*rho) ) * cos(phi)+ 1.0/2.0*lambda_SR * cexp(2.0*I*beta) * rho * cos(phi) );

	COMPLEX16 P11    = 1.0/2.0*sqrt(1-lambda_PD) * sqrt(lambda_SR) * tau * ( -2.0*rho*cexp(2.0*I*beta)+2.0*cos(2.0*phi)+Kappa*sin(2.0*phi) );
	COMPLEX16 P22    = P11;
	COMPLEX16 P12    = -sqrt(1.0-lambda_PD) * sqrt(lambda_SR)*tau*sin(phi)*(2.0*cos(phi)+Kappa*sin(phi) );
	COMPLEX16 P21    =  sqrt(1.0-lambda_PD) * sqrt(lambda_SR)*tau*cos(phi)*(2.0*sin(phi)-Kappa*cos(phi) );

	COMPLEX16 Q11    = sqrt(lambda_PD) * ( cexp(-2.0*I*beta)+rho*rho*cexp(2.0*I*beta)-rho*(2.0*cos(2.0*phi)+Kappa*sin(2.0*phi)) + 1.0/2.0*epsilon*rho * (cexp(-2.0*I*beta)*cos(2.0*phi)+cexp(2.0*I*beta)* ( -2.0*rho-2.0*rho*cos(2.0*beta)+cos(2.0*phi)+Kappa*sin(2.0*phi) ) + 2.0*cos(2.0*phi)+3.0*Kappa*sin(2.0*phi))-1.0/2.0*lambda_SR*rho * ( 2.0*rho*cexp(2.0*I*beta)-2.0*cos(2.0*phi)-Kappa*sin(2.0*phi) ) );
	COMPLEX16 Q22    = Q11;
	COMPLEX16 Q12    = 0.0;
	COMPLEX16 Q21    = 0.0;

	COMPLEX16 N11    = sqrt(1.0-lambda_PD) * sqrt(epsilon/2.0)*tau *(Kappa*(1.0+rho*cexp(2.0*I*beta))*sin(phi)+2.0*cos(beta)*(cexp(-I*beta)*cos(phi)-rho*cexp(I*beta)*(cos(phi)+Kappa*sin(phi))));
	COMPLEX16 N22    = -sqrt(1.0-lambda_PD)*sqrt(2.0*epsilon)*tau*(-cexp(-I*beta)+rho*cexp(I*beta))*cos(beta)*cos(phi);
	COMPLEX16 N12    = -sqrt(1.0-lambda_PD)*sqrt(2.0*epsilon)*tau*(cexp(-I*beta)+rho*cexp(I*beta))*cos(beta)*sin(phi);
	COMPLEX16 N21    = sqrt(1.0-lambda_PD)*sqrt(2.0*epsilon)*tau*(-Kappa*(1.0+rho)*cos(phi)+2.0*cos(beta)*(cexp(-I*beta)+rho*cexp(I*beta))*cos(beta)*sin(phi));


	//>>>>>>>>    QUANTUM NOISE POWER SPECTRAL DENSITY [BC, 5.13]   <<<<<<<<<<<<<<<<<

	double n = h_SQL*h_SQL/(2.0*Kappa*tau*tau*pow(cabs(D1_L*sin(zeta)+D2_L*cos(zeta)),2.0))*(
    			pow(cabs(C11_L*sin(zeta)+C21_L*cos(zeta)),2.0)+
    			pow(cabs(C12_L*sin(zeta)+C22_L*cos(zeta)),2.0)+
    			pow(cabs(P11*sin(zeta)+P21*cos(zeta)),2.0)+
    			pow(cabs(P12*sin(zeta)+P22*cos(zeta)),2.0)+
    			pow(cabs(Q11*sin(zeta)+Q21*cos(zeta)),2.0)+
    			pow(cabs(Q12*sin(zeta)+Q22*cos(zeta)),2.0)+
    			pow(cabs(N11*sin(zeta)+N21*cos(zeta)),2.0)+
    			pow(cabs(N12*sin(zeta)+N22*cos(zeta)),2.0));

	return n;
}


/*
 *
 * FIRST GENERATION DETECTORS
 *
 */


/**
 * Provides the noise power spectrum based on a phenomenological fit
 * to the SRD curve for iLIGO.
 *
 * This is a fit to the data provided for the Science Requirements Document
 * (SRD) curve for initial LIGO given, which can be found at
 * http://www.ligo.caltech.edu/~jzweizig/distribution/LSC_Data/
 *
 * The Science Requirements Document is located at
 * http://www.ligo.caltech.edu/docs/E/E950018-02.pdf
 */
double XLALSimNoisePSDiLIGOSRD(double f /**< frequency (Hz) */)
{
	const double aseis = 1.57271;
	const double pseis = -14.0;
	const double athrm = 3.80591e-19;
	const double pthrm = -2.0;
	const double ashot = 1.12277e-23;
	const double fshot = 89.3676;
	double seis = aseis * aseis * pow(f, 2.0*pseis);
	double thrm = athrm * athrm * pow(f, 2.0*pthrm);
	double shot = ashot * ashot * (1.0 + pow(f / fshot, 2.0));
	return seis + thrm + shot;
}


/**
 * Provides the seismic noise power spectrum for iLIGO.
 *
 * Note: only valit for f > 10 Hz.
 * This is mostly a phenomenological fit.
 */
double XLALSimNoisePSDiLIGOSeismic(double f /**< frequency (Hz) */)
{
	double seismic;
	/* four layers of suspension */
	seismic = XLALSimNoisePSDSeismic(f,
		LAL_ILIGO_ARMLENGTH_SI,
		LAL_ILIGO_THERMAL_SUSP_FREQ_SI,
		LAL_ILIGO_THERMAL_STACK_FREQ_SI,
		4);
	return seismic;
}


/**
 * Provides the thermal noise (suspension + coating) power spectrum for iLIGO.
 *
 * Note: this is a phenomenological fit to the broadband component.
 */
double XLALSimNoisePSDiLIGOThermal(double f /**< frequency (Hz) */)
{
	double susp;
	double coat;
	susp = XLALSimNoisePSDSuspTherm(f,
			LAL_ILIGO_ARMLENGTH_SI,
			LAL_ILIGO_MIRROR_MASS_SI,
			LAL_ILIGO_TEMPERATURE_SI,
			LAL_ILIGO_THERMAL_SUSP_FREQ_SI,
			LAL_ILIGO_THERMAL_SUSP_QUAL);
	coat = XLALSimNoisePSDMirrorTherm(f,
			LAL_ILIGO_ARMLENGTH_SI,
			LAL_ILIGO_MIRROR_MASS_SI,
			LAL_ILIGO_TEMPERATURE_SI,
			LAL_ILIGO_THERMAL_COAT_FREQ_SI,
			LAL_ILIGO_THERMAL_COAT_QUAL);
	return susp + coat;
}


/**
 * Provides the shot noise power spectrum for iLIGO.
 *
 * Note: the effective quantum efficiency is one-third the actual quantum
 * efficiency owing to the RF readout scheme.  A fiducial value of 250 W
 * of power on the beamsplitter is used.
 */
double XLALSimNoisePSDiLIGOShot(double f /**< frequency (Hz) */)
{
	const double eta = 0.9/3.0; // effective quantum efficiency of photodiode for RF readout
	double shot;
	shot = XLALSimNoisePSDShot(f,
			LAL_ILIGO_LASER_POWER_BS_SI,
			LAL_ILIGO_LASER_WAVELENGTH_SI,
			LAL_ILIGO_ARMLENGTH_SI,
			LAL_ILIGO_FINESSE,
			eta);
	return shot;
}


/**
 * Provides the shot noise power spectrum for eLIGO.
 *
 * Note: A fiducial value of 250 W of power on the beamsplitter is used.
 */
double XLALSimNoisePSDeLIGOShot(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode for DC readout
	double shot;
	shot = XLALSimNoisePSDShot(f,
			LAL_ILIGO_LASER_POWER_BS_SI,
			LAL_ILIGO_LASER_WAVELENGTH_SI,
			LAL_ILIGO_ARMLENGTH_SI,
			LAL_ILIGO_FINESSE,
			eta);
	return shot;
}


/**
 * Provides the noise power spectrum for a model of the iLIGO detector.
 *
 * Warning: not all noise sources are correctly accounted for (in particular,
 * there is no actuation noise modelled) so this noise spectrum does not
 * correspond to the S5 spectrum.
 */
double XLALSimNoisePSDiLIGOModel(double f /**< frequency (Hz) */)
{
	double seismic;
	double thermal;
	double shot;

	seismic = XLALSimNoisePSDiLIGOSeismic(f);
	thermal = XLALSimNoisePSDiLIGOThermal(f);
	shot = XLALSimNoisePSDiLIGOShot(f);

	return shot + seismic + thermal;
}


/**
 * Provides the noise power spectrum for a model of the eLIGO detector.
 *
 * Warning: not all noise sources are correctly accounted for so this noise
 * spectrum does not correspond to the S6 spectrum.
 */
double XLALSimNoisePSDeLIGOModel(double f /**< frequency (Hz) */)
{
	double seismic;
	double thermal;
	double shot;

	seismic = XLALSimNoisePSDiLIGOSeismic(f);
	thermal = XLALSimNoisePSDiLIGOThermal(f);
	shot = XLALSimNoisePSDeLIGOShot(f);

	return shot + seismic + thermal;
}

/**
 * Provides the design noise power spectrum for Virgo based on a
 * phenomenological fit (from the Virgo webiste) that can be approximated by the
 * following:
 * \f{equation}{
 * S_h(f) =
 * s_0 \left ( \frac {7.87f}{f_0} \right )^{-4.8} + \frac{6}{17} \frac{f_0}{f}
 * + \left [1 + \left (\frac {f}{f_0} \right)^2 \right ],
 * \f}
 * where \f$s_0=10.2e-46\f$.
 *
 * Warning: This comes from the deprecated function LALVIRGOPsd in the lal
 * noisemodels package, which comes with no reference to the curve. An updated
 * version of this model, with a reference would be welcomed.
 */
double XLALSimNoisePSDVirgo(double f /**< frequency (Hz) */)
{
  REAL8 s0, x;

  x = f/500.;

  s0 = 10.2e-46;

  return s0*( pow(7.87*x,-4.8) + 6./17./x + 1. + x*x);
}

/**
 * Provides a GEO noise power spectrum based on that from Table IV of
 * \ref dis2001.
 *
 * The comes from the deprecated function LALGEOPsd in the lal noisemodels
 * package.
 */
double XLALSimNoisePSDGEO(double f /**< frequency (Hz) */)
{
  REAL8 x, seismic, thermal, shot;

  x = f/150.;
  seismic = pow(10.,-16.) *  pow(x,-30.);
  thermal = 34. / x;
  shot = 20. * (1 - pow(x,2.) + 0.5 * pow(x,4.)) / (1. + 0.5 * pow(x,2.));

  return 1e-46*(seismic + thermal + shot);
}


/**
 * Provides a GEO-HF noise power spectrum based on a fit to Figure 6
 * from \ref Grote2010.
 *
 * The fit is good between 50Hz to 8kHz and errors between the analytic
 * fit given and the <a href="https://intranet.aei.uni-hannover.de/geo600/geohflogbook.nsf/7e8722dffa24dea0c1256de900406c84/4837a612ac990060c12575ce004e70fd?OpenDocument">estimated curve</a> are less than 1%.
 */
double XLALSimNoisePSDGEOHF(double f /**< frequency (Hz) */)
{
  REAL8 f2 = f*f;

  return 7.18e-46*(1. + (f2/(1059.*1059.))) + (4.90e-41/f2) + (8.91e-43/f) + (1.6e-17/pow(f, 16.));
}


/**
 * Provides a TAMA300 noise power spectrum based on that from Table IV of
 * \ref dis2001.
 *
 * The comes from the deprecated function LALTAMAPsd in the lal noisemodels
 * package.
 */
double XLALSimNoisePSDTAMA(double f /**< frequency (Hz) */)
{
  REAL8 seismic, thermal, shot, x;

  x = f/400.;
  seismic = pow(x,-5);
  thermal = 13. / x;
  shot = 9. * (1. + x*x);

  return 75.e-46*(seismic + thermal + shot);
}

/*
 *
 * SECOND GENERATION DETECTORS
 *
 */


/**
 * Provides the thermal noise (suspension + coating) power spectrum for aLIGO.
 *
 * Note: this is a phenomenological fit to the broadband component.
 */
double XLALSimNoisePSDaLIGOThermal(double f /**< frequency (Hz) */)
{
	double susp;
	double coat;
	susp = XLALSimNoisePSDSuspTherm(f,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_TEMPERATURE_SI,
			LAL_ALIGO_THERMAL_SUSP_FREQ_SI,
			LAL_ALIGO_THERMAL_SUSP_QUAL);
	coat = XLALSimNoisePSDMirrorTherm(f,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_TEMPERATURE_SI,
			LAL_ALIGO_THERMAL_COAT_FREQ_SI,
			LAL_ALIGO_THERMAL_COAT_QUAL);
	return susp + coat;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the low-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled No SRM.
 */
double XLALSimNoisePSDaLIGOQuantumNoSRMLowPower(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double ds = 0.0; // detuning phase -- no SRM! (rad)
	const double zeta = 130.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f,
			LAL_ALIGO_LASER_POWER_LOW_SI,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			1.0 /* no SRM! */,
			ds, zeta, eta);

	return quantum;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the high-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is the same a No SRM but with 125 W laser power.
 */
double XLALSimNoisePSDaLIGOQuantumNoSRMHighPower(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double ds = 0.0; // detuning phase -- no SRM! (rad)
	const double zeta = 130.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f,
			LAL_ALIGO_LASER_POWER_HIGH_SI,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			1.0 /* no SRM! */,
			ds, zeta, eta);

	return quantum;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the low-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, Low Power.
 */
double XLALSimNoisePSDaLIGOQuantumZeroDetLowPower(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double ds = 0.0; // detuning phase (rad)
	const double zeta = 116.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f,
			LAL_ALIGO_LASER_POWER_LOW_SI,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			LAL_ALIGO_SRM_TRANSMITTANCE,
			ds, zeta, eta);

	return quantum;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the high-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, High Power.
 */
double XLALSimNoisePSDaLIGOQuantumZeroDetHighPower(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double ds = 0.0; // detuning phase (rad)
	const double zeta = 116.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f,
			LAL_ALIGO_LASER_POWER_HIGH_SI,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			LAL_ALIGO_SRM_TRANSMITTANCE,
			ds, zeta, eta);

	return quantum;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to NS-NS inspirals.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled NS-NS Opt.
 */
double XLALSimNoisePSDaLIGOQuantumNSNSOpt(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double ds = 11.0 * LAL_PI_180; // detuning phase (rad)
	const double zeta = 103.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f,
			LAL_ALIGO_LASER_POWER_HIGH_SI,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			LAL_ALIGO_SRM_TRANSMITTANCE,
			ds, zeta, eta);

	return quantum;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to 30+30 solar mass binary
 * black holes with fixed signal recycling cavity detuning of 20 degrees.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled BHBH 20-degree Detune.
 */
double XLALSimNoisePSDaLIGOQuantumBHBH20Deg(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double I0 = 20.0; // input power (W)
	const double ds = 20.0 * LAL_PI_180; // detuning phase (rad)
	const double zeta = 105.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f, I0,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			LAL_ALIGO_SRM_TRANSMITTANCE,
			ds, zeta, eta);

	return quantum;
}


/**
 * Provides the quantum noise power spectrum for aLIGO under the
 * configuration tuned to narrow-band high-frequency sensitivity around
 * 1 kHz.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled High Freq.
 */
double XLALSimNoisePSDaLIGOQuantumHighFrequency(double f /**< frequency (Hz) */)
{
	const double eta = 0.9; // quantum efficiency of photodiode
	const double T_SRM = 0.011; // SRM Transmittance
	const double ds = 4.7 * LAL_PI_180; // detuning phase (rad)
	const double zeta = 128.0 * LAL_PI_180; // homodyne detection phase (rad)
	double quantum;

	quantum = XLALSimNoisePSDQuantum(f,
			LAL_ALIGO_LASER_POWER_HIGH_SI,
			LAL_ALIGO_LASER_WAVELENGTH_SI,
			LAL_ALIGO_ARMLENGTH_SI,
			LAL_ALIGO_MIRROR_MASS_SI,
			LAL_ALIGO_MIRROR_LOSS,
			LAL_ALIGO_BS_LOSS,
			LAL_ALIGO_ITM_TRANSMITTANCE,
			LAL_ALIGO_PRM_TRANSMITTANCE,
			T_SRM, ds, zeta, eta);

	return quantum;
}


/**
 * Provides the noise power spectrum for aLIGO under the low-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled No SRM.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGONoSRMLowPower(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumNoSRMLowPower(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for aLIGO under the high-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is the same a No SRM but with 125 W laser power.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGONoSRMHighPower(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumNoSRMHighPower(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for aLIGO under the low-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, Low Power.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOZeroDetLowPower(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumZeroDetLowPower(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for aLIGO under the high-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, High Power.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOZeroDetHighPower(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumZeroDetHighPower(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to NS-NS inspirals.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled NS-NS Opt.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGONSNSOpt(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumNSNSOpt(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to 30+30 solar mass binary
 * black holes with fixed signal recycling cavity detuning of 20 degrees.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled BHBH 20-degree Detune.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOBHBH20Deg(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumBHBH20Deg(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for aLIGO under the
 * configuration tuned to narrow-band high-frequency sensitivity around
 * 1 kHz.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled High Freq.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOHighFrequency(double f /**< frequency (Hz) */)
{
	double quantum;
	double thermal;

	quantum = XLALSimNoisePSDaLIGOQuantumHighFrequency(f);
	thermal = XLALSimNoisePSDaLIGOThermal(f);

	return quantum + thermal;
}


/**
 * Provides the noise power spectrum for KAGRA based on that from Eqn 5 of
 * \ref md2012. This is a phenomenological fit to the KAGRA spectrum from
 * http://gwcenter.icrr.u-tokyo.ac.jp/en/researcher/parameter
 */
double XLALSimNoisePSDKAGRA(double f /**< frequency (Hz) */)
{
  REAL8 x = log(f / 100.);
  REAL8 x2 = x*x;
  REAL8 asd = 0.;

  /* calculate ASD from reference */
  asd = 6.499e-25 * ( 9.72e-9*exp(-1.43 - 9.88*x - 0.23*x2)
                     + 1.17*exp(0.14 - 3.10*x - 0.26*x2)
                     + 1.70*exp(0.14 + 1.09*x - 0.013*x2)
                     + 1.25*exp(0.071 + 2.83*x - 4.91*x2) );

  /* return PSD */
  return asd*asd;
}


/**
 * Provides the noise power spectrum for AdvVirgo based on that from Eqn 6 of
 * \ref md2012. This is a phenomenological fit to the AdvVirgo spectrum from
 * http://wwwcascina.virgo.infin.it/advirgo.
 */
double XLALSimNoisePSDAdvVirgo(double f /**< frequency (Hz) */)
{
  REAL8 x = log(f / 300.);
  REAL8 x2 = x*x;
  REAL8 asd = 0.;

  /* calculate ASD from reference */
  asd = 1.259e-24 * ( 0.07*exp(-0.142 - 1.437*x + 0.407*x2)
                     + 3.1*exp(-0.466 - 1.043*x - 0.548*x2)
                     + 0.4*exp(-0.304 + 2.896*x - 0.293*x2)
                     + 0.09*exp(1.466 + 3.722*x - 0.984*x2) );

  /* return PSD */
  return asd*asd;
}


/**
 * Evaluates a power spectral density function, psdfunc, at the frequencies required
 * to populate the frequency series psd, with a low frequency cutoff flow.
 */
int XLALSimNoisePSD(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow,			/**< low frequency cutoff (Hz) */
	double (*psdfunc)(double)	/**< function that provides the PSD at a specified frequency */
)
{
	size_t kmin;
	size_t k;

	/* set DC and Nyquist to zero */
	/* note: assumes last element is Nyquist */
	psd->data->data[0] = psd->data->data[psd->data->length - 1] = 0.0;

	/* determine low frequency cutoff */
	kmin = flow / psd->deltaF;

	psd->data->data[0] = 0.0; /* set DC to zero */
	for (k = 1; k < kmin; ++k) /* set low frequency components to zero */
		psd->data->data[k] = 0.0;
	for (; k < psd->data->length - 1; ++k) /* evaluate psdfn for frequencies in requested band */
		psd->data->data[k] = (*psdfunc)(k * psd->deltaF);
	psd->data->data[psd->data->length - 1] = 0.0; /* set Nyquist to zero (presume this is Nyquist!) */

	return 0;
}


/**
 * Reads file fname containing two-column amplitude spectral density data file
 * and interpolates at the frequencies required to populate the frequency
 * series psd, with a low frequency cutoff flow.
 */
int XLALSimNoisePSDFromFile(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow,			/**< low frequency cutoff (Hz) */
	const char *fname		/**< file containing amplitude spectral density data */
)
{
	double *f;
	double *h;
	size_t  n;
	size_t  i;
	size_t  kmin;
	size_t  k;
	LALFILE *fp;

	/* first, read the data form the datafile */
	fp = XLALSimReadDataFileOpen(fname);
	if (!fp)
		XLAL_ERROR(XLAL_EFUNC);
	n = XLALSimReadDataFile2Col(&f, &h, fp);
	XLALFileClose(fp);
	if (n == (size_t)(-1))
		XLAL_ERROR(XLAL_EFUNC);
	/* take the log of the amplitude spectral density data */
	for (i = 0; i < n; ++i)
		h[i] = log(h[i]);

	/* set DC and Nyquist to zero */
	/* note: assumes last element is Nyquist */
	psd->data->data[0] = psd->data->data[psd->data->length - 1] = 0.0;

	/* determine low frequency cutoff */
	kmin = flow / psd->deltaF;

	i = 1;
	psd->data->data[0] = 0.0; /* set DC to zero */
	for (k = 1; k < kmin; ++k) /* set low frequency components to zero */
		psd->data->data[k] = 0.0;
	for (; k < psd->data->length - 1; ++k) {
		double fk = k * psd->deltaF; /* target frequency */
		double hk;
		double x;
		/* interpolate data for this frequency value */
		while (f[i] < fk && i < n - 1)
			++i;
		x = (f[i] - fk) / (f[i] - f[i-1]);
		hk = x * h[i-1] + (1.0 - x) * h[i];
		/* power spectrum is exp( 2 * log(amplitude spectrum) ) */
		psd->data->data[k] = exp(2.0 * hk);
	}
	psd->data->data[psd->data->length - 1] = 0.0; /* set Nyquist to zero (presume this is Nyquist!) */

	XLALFree(h);
	XLALFree(f);
	return 0;
}

/* prefix for noise psd files provided by LIGO-T0900288 */
#define T0900288 "LIGO-T0900288-v3-"

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "NO_SRM.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGONoSRMLowPowerGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
)
{
	return XLALSimNoisePSDFromFile(psd, flow, T0900288 "NO_SRM.txt");
}

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "ZERO_DET_low_P.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOZeroDetLowPowerGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
)
{
	return XLALSimNoisePSDFromFile(psd, flow, T0900288 "ZERO_DET_low_P.txt");
}

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "ZERO_DET_high_P.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOZeroDetHighPowerGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
)
{
	return XLALSimNoisePSDFromFile(psd, flow, T0900288 "ZERO_DET_high_P.txt");
}

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "NSNS_Opt.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGONSNSOptGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
)
{
	return XLALSimNoisePSDFromFile(psd, flow, T0900288 "NSNS_Opt.txt");
}

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "BBH_20deg.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOBHBH20DegGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
)
{
	return XLALSimNoisePSDFromFile(psd, flow, T0900288 "BHBH_20deg.txt");
}

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "High_Freq.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOHighFrequencyGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
)
{
	return XLALSimNoisePSDFromFile(psd, flow, T0900288 "High_Freq.txt");
}



/*
 *
 * TEST CODE
 *
 */

#if 0

/*
 * This test produces plot data that should be similar to the theoretical noise
 * component curves for shot, seismic, suspension thermal, and mirror thermal
 * shown in Figure 7 of Rep. Prog. Phys. 72 (2009) 076901.  It also includes the
 * SRD curve for iLIGO.
 */
int test_iligo_psd(void)
{
	const double eta = 0.9/3.0;
	double f;
	FILE *fp = fopen("psd_iligo.dat", "w");
	for (f = 30.0; f < 8000.0; f *= 1.01) {
		double S_susp;
		double S_coat;
		double S_shot;
		double S_seis;

		S_shot = XLALSimNoisePSDShot(f,
				LAL_ILIGO_LASER_POWER_BS_SI,
				LAL_ILIGO_LASER_WAVELENGTH_SI,
				LAL_ILIGO_ARMLENGTH_SI,
				LAL_ILIGO_FINESSE,
				eta);

		S_seis = XLALSimNoisePSDSeismic(f,
				LAL_ILIGO_ARMLENGTH_SI,
				LAL_ILIGO_THERMAL_SUSP_FREQ_SI,
				LAL_ILIGO_THERMAL_STACK_FREQ_SI,
				4);

		S_susp = XLALSimNoisePSDSuspTherm(f,
				LAL_ILIGO_ARMLENGTH_SI,
				LAL_ILIGO_MIRROR_MASS_SI,
				LAL_ILIGO_TEMPERATURE_SI,
				LAL_ILIGO_THERMAL_SUSP_FREQ_SI,
				LAL_ILIGO_THERMAL_SUSP_QUAL);

		S_coat = XLALSimNoisePSDMirrorTherm(f,
				LAL_ILIGO_ARMLENGTH_SI,
				LAL_ILIGO_MIRROR_MASS_SI,
				LAL_ILIGO_TEMPERATURE_SI,
				LAL_ILIGO_THERMAL_COAT_FREQ_SI,
				LAL_ILIGO_THERMAL_COAT_QUAL);

		fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", f, sqrt(S_shot + S_seis + S_susp + S_coat), sqrt(S_shot), sqrt(S_seis), sqrt(S_susp), sqrt(S_coat), sqrt(XLALSimNoisePSDiLIGOSRD(f)));
	}
	fclose(fp);
	return 0;
}

/*
 * This test produces plot data that should be similar to the aLIGO noise
 * curves shown in LIGO-T0900288-v3.
 */
int test_aligo_psd(void)
{
	double f;
	FILE *fp = fopen("psd_aligo.dat", "w");
	for (f = 9.0; f < 3000.0; f *= 1.01)
		fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", f,
			sqrt(XLALSimNoisePSDaLIGONoSRMLowPower(f)),
			sqrt(XLALSimNoisePSDaLIGOZeroDetLowPower(f)),
			sqrt(XLALSimNoisePSDaLIGOZeroDetHighPower(f)),
			sqrt(XLALSimNoisePSDaLIGONSNSOpt(f)),
			sqrt(XLALSimNoisePSDaLIGOBHBH20Deg(f)),
			sqrt(XLALSimNoisePSDaLIGOHighFrequency(f)));
	fclose(fp);
	return 0;
}

int main(void)
{
	XLALSetErrorHandler(XLALAbortErrorHandler);
	test_iligo_psd();
	test_aligo_psd();
	LALCheckMemoryLeaks();
	return 0;
}

#endif
