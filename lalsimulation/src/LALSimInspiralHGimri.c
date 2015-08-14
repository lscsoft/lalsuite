
//====================================================================================================
// Code to generate a quasi-circular intermediate mass-ratio inspiral trajectory and graviational
// waveform following Huerta & Gair 2011 (1009.1985). The system's evolution is divided into four
// stages. See Huerta & Gair 2011 for discussion of each.
//
//	1. Adiabatic inspiral
//	2. Transition (see also Ori & Thorne 2000 -- gr-qc/0003032)
//	3. Geodesic Plunge
//	4. Quasi-Normal Ringdown
//
// System is evolved using 2PN angular momentum flux from Gair and Glampedakis 2006 (gr-qc/0510129)
// with higher order fits to Teukolsky results. Azimuthal phase is evolved using the 2PN conservative
// corrections discussed in Huerta & Gair 2008 (0812.4208), and the waveform is generated with the
// flat space quadrupole formula (see Babak et al 2007 -- gr-qc/0607007).
//
// Tom Callister
// thomas.a.callister@gmail.com
// 2014
//
//====================================================================================================

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>

#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#define c ((REAL8)(LAL_C_SI))
#define G ((REAL8)(LAL_G_SI))
#define pi ((REAL8)(LAL_PI))
#define pc ((REAL8)(LAL_PC_SI))
#define Msun ((REAL8)(LAL_MSUN_SI))
#define Msun_sec ((REAL8)(G*Msun/(c*c*c)))
#define GPC_sec ((REAL8)(pc*pow(10.,9.)/c))

//Function prototypes
static REAL8 XLALHGimri_AngMomFlux(REAL8 q, REAL8 r, REAL8 nu);
static INT4 XLALHGimri_PlusEquations(const gsl_vector * x, void *params, gsl_vector * f);
static INT4 XLALHGimri_CrossEquations(const gsl_vector * z, void *params, gsl_vector * g);
static REAL8 XLALHGimri_dFdr(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r);
static REAL8 XLALHGimri_initialP(REAL8 p0, void *params);
static REAL8 XLALHGimri_dLzdr(REAL8 q, REAL8 r);
static INT4 HGimri_start(REAL8 m, REAL8 M, REAL8 q, REAL8 D, REAL8 Sdotn, REAL8 phi0, REAL8 p0, REAL8Sequence *hplus, REAL8Sequence *hcross, REAL8 dt, UINT4 Npts);

INT4 XLALHGimriGenerator(REAL8TimeSeries **hplus,REAL8TimeSeries **hcross,REAL8 phi0,REAL8 dt,REAL8 m1,REAL8 m2,REAL8 f_min,REAL8 r,REAL8 inc,REAL8 s1z);

//Data type to hold current simulation regime
enum stage {INSPIRAL, TRANSITION, PLUNGE, FINALPLUNGE};

//Structure template used to match plunging waveform onto QNM ringdown at three points across the light-ring (LR)
struct rparams {
	REAL8 M;		//Total binary mass (units seconds)
	REAL8 Sdotn;		//Cosine of inclination angle (relative to detector)
	REAL8 final_mass;	//Final BH mass (units seconds)
	REAL8 final_q;		//Final BH reduced spin
	REAL8 w0, w1, w2;	//Frequencies of n=0,1,2 overtones of the (l,m)=(2,+/-2) QNMs, real parts
	REAL8 wi0, wi1, wi2;	//Frequencies of n=0,1,2 overtones of the (l,m)=(2,+/-2) QNMs, imaginary parts
	REAL8 ab_factor;	//Amplitude factor
	REAL8 hp_1, dhpdi_1;	//hplus and its derivative at point 1
	REAL8 hp_2, dhpdi_2;	//hplus and its derivative at point 2
	REAL8 hp_3, dhpdi_3;	//hplus and its derivative at point 3
	REAL8 hx_1, dhxdi_1;	//hcross and its derivative at point 1
	REAL8 hx_2, dhxdi_2;	//hcross and its derivative at point 2
	REAL8 hx_3, dhxdi_3;	//hcross and its derivative at point 3
	REAL8 dt;		//Dimensionless time step
	};

//Structure template used in computing in initial radius
struct pParams {
	REAL8 m;		//CO mass (seconds)
	REAL8 M;		//BH mass (seconds)
	REAL8 q;		//BH reduced spin
	REAL8 f_min;		//Desired initial GW frequency
	};

static REAL8 XLALHGimri_initialP(REAL8 p0, void *params) {

	//====================================================
	// Root-finding function used to find initial radius p
	// corresponding to a GW frequency f_min.
	//====================================================

	//Read params
	struct pParams *p = (struct pParams *) params;
	REAL8 m = p->m;
	REAL8 M = p->M;
	REAL8 q = p->q;
	REAL8 f_min = p->f_min;
	REAL8 nu = (m*M)/pow(m+M,2.);

	//Conservative corrections
	REAL8 d0	= 1./8.;
	REAL8 d1	= 1975./896.;
	REAL8 d1_5	= -27.*pi/10.;
	REAL8 l1_5	= -191./160.;
	REAL8 d2	= 1152343./451584.;

	//Match orbital frequency to twice the GW frequency f_min
	return 	(1./(pow(p0,3/2.)+q))*(1. + nu*(d0+d1/p0+(d1_5+q*l1_5)*pow(p0,-1.5)+d2*pow(p0,-2.0))) - pi*(f_min)*(m+M)*Msun_sec;

	}


static REAL8 XLALHGimri_AngMomFlux(REAL8 q, REAL8 r, REAL8 nu) {

	//===================================================================
	// Angular momentum flux from Gair & Glampedakis 2006 (gr-qc/0510129).
	// See Eq (45) for 2PN fluxes, and (57) for the higher order fits to
	// circular Teukolsky data.
	//===================================================================

        REAL8 pref = -6.4/(pow(r,7./2.));
        REAL8 rtr=sqrt(r);

	//Cosine of orbital inclination. 1 if prograde, -1 if retrograde
	REAL8 cosInc;
	REAL8 sign;

	if (q >= 0.0) {
		cosInc = 1.;
		sign = 1.;
		}
	else {
		cosInc = -1.;
		sign = -1.;
		}

	//Absolute value of spin
	REAL8 q_abs = fabs(q);

	//Higher order Teukolsky terms
	REAL8 c1a,c1b,c1c;
	REAL8 c2a,c2b,c2c;
	REAL8 c3a,c3b,c3c;
	REAL8 c4a,c4b,c4c;
	REAL8 c5a,c5b,c5c;
	REAL8 c6a,c6b,c6c;
	REAL8 c7a,c7b,c7c;
	REAL8 c8a,c8b,c8c;
	REAL8 c9a,c9b,c9c;
	REAL8 c10a,c10b,c10c;
	REAL8 c11a,c11b,c11c;
	REAL8 higherorderfit,c1,c2,correction,circbit;

	c1a=-10.741956;
	c1b=28.5942157;
	c1c=-9.077378144;
	c1=c1a+(c1b+c1c/rtr)/rtr;

	c2a=-1.428362761;
	c2b=10.70029768;
	c2c=-33.70903016;
	c2=(c2a+(c2b+c2c/rtr)/rtr);

	c3a=-28.15174147;
	c3b=60.9607071973;
	c3c=40.99984205;

	c4a=-0.348161211;
	c4b=2.37258476;
	c4c=-66.65840948;

	c5a=-0.715392387;
	c5b=3.21592568;
	c5c=5.28887649;

	c6a=-7.6103411;
	c6b=128.87778309;
	c6c=-475.4650442;

	c7a=12.290783385;
	c7b=-113.1250548;
	c7c=306.11883292;

	c8a=40.9258725;
	c8b=-347.2713496;
	c8c=886.50332051;

	c9a=-25.48313727;
	c9b=224.22721861;
	c9c=-490.98212316;

	c10a=-9.006337706;
	c10b=91.17666278;
	c10c=-297.001939215;

	c11a=-0.64500047;
	c11b=-5.13591989;
	c11c=47.19818628;

	correction = q_abs*(
			736.2086781-283.9553066*rtr+q_abs*(-1325.1852209+483.266206498*rtr)+pow(q_abs,2.0)*(634.49936445-219.223848944*rtr)
			+(82.07804475-25.82025864*rtr)+q_abs*(-904.16109275+301.477789146*rtr)+pow(q_abs,2.0)*(827.31891826-271.9659423*rtr)
			)/(r*rtr);

	higherorderfit = q_abs*c1 + pow(q_abs,3.)*c2
			+ cosInc*( c3a+(c3b+c3c/rtr)/rtr )
			+ pow(q_abs,2.)*cosInc*(c4a+(c4b+c4c/rtr)/rtr)
			+ pow(q_abs,4.)*cosInc*(c5a+(c5b+c5c/rtr)/rtr)
			+ q_abs*(c6a+(c6b+c6c/rtr)/rtr)
			+ pow(q_abs,3.)*(c7a+(c7b+c7c/rtr)/rtr)
			+ pow(q_abs,2.)*cosInc*(c8a+(c8b+c8c/rtr)/rtr)
			+ pow(q_abs,4.)*cosInc*(c9a+(c9b+c9c/rtr)/rtr)
			+ pow(q_abs,3.)*(c10a+(c10b+c10c/rtr)/rtr)
			+ pow(q_abs,4.)*cosInc*(c11a+(c11b+c11c/rtr)/rtr)
			+ correction;

	//Bracketed flux terms from Eq (45) with higher order corrections
	circbit = cosInc - (61./12.)*q_abs/(r*rtr) -1247.*cosInc/(336.*r) + 4.*pi*cosInc/(r*rtr) - 44711.*cosInc/(9072.*r*r)
			+ (33./16.)*q_abs*q_abs*cosInc/(r*r) + higherorderfit/(r*r*rtr);

	//Return full flux. (sign) factor needed to guarantee that Ldot is negative at leading order,
	//following our adopted spin convention.
	return(sign*nu*pref*circbit);

	}

static INT4 XLALHGimri_PlusEquations(const gsl_vector * x, void *params, gsl_vector * f) {

	//=============================================================================
	// Used by gsl_multiroots to match plus-polarized plunge waveform onto ringdown.
	// Ringdown waveform consists of (l,m)=(2,+/-2) QNMs, with n=0,1,2 overtones. The
	// RD waveform and its time derivative is matched to plunge waveform at three
	// points across the light ring, considering the n=0 fundamental tone at the first
	// time step, the n=0,1 tones at the second, and the n=0,1,2 tones at the third
	//=============================================================================

	//See definition rparams above for structure information
	REAL8 dt		= ((struct rparams *) params)->dt;
        REAL8 wi0 		= ((struct rparams *) params)->wi0;
        REAL8 w0 		= ((struct rparams *) params)->w0;
        REAL8 wi1 		= ((struct rparams *) params)->wi1;
        REAL8 w1		= ((struct rparams *) params)->w1;
        REAL8 wi2 		= ((struct rparams *) params)->wi2;
        REAL8 w2 		= ((struct rparams *) params)->w2;
        REAL8 final_mass 	= ((struct rparams *) params)->final_mass;
        REAL8 final_q 		= ((struct rparams *) params)->final_q;
        REAL8 Sdotn 		= ((struct rparams *) params)->Sdotn;
        REAL8 M 		= ((struct rparams *) params)->M;
	REAL8 ab_factor 	= ((struct rparams *) params)->ab_factor;

        REAL8 hp_1 		= ((struct rparams *) params)->hp_1;
        REAL8 dhpdi_1 		= ((struct rparams *) params)->dhpdi_1;
        REAL8 hp_2 		= ((struct rparams *) params)->hp_2;
        REAL8 dhpdi_2 		= ((struct rparams *) params)->dhpdi_2;
        REAL8 hp_3 		= ((struct rparams *) params)->hp_3;
        REAL8 dhpdi_3		= ((struct rparams *) params)->dhpdi_3;

	//Leading amplitude factors
        REAL8 aone		= 1.+ 2.*Sdotn + Sdotn*Sdotn;
        REAL8 atwo		= 2.+ Sdotn - 4.*Sdotn*Sdotn - 3.*Sdotn*Sdotn*Sdotn;

	//Unknown RD constants
        const REAL8 a0n 	= gsl_vector_get (x, 0);
        const REAL8 a0p 	= gsl_vector_get (x, 1);
	const REAL8 a1n 	= gsl_vector_get (x, 2);
        const REAL8 a1p 	= gsl_vector_get (x, 3);
	const REAL8 a2n 	= gsl_vector_get (x, 4);
        const REAL8 a2p 	= gsl_vector_get (x, 5);

	//Functions of the form: (Ringdown Quantity) - (Plunge Quantity). We seek the roots of these functions

	//h_plus at time step 1
	const REAL8 yy0 = ab_factor*a0n*(aone-(2./9.)*atwo*w0*final_q) - hp_1;

	//h_plus time derivative at time step 1
        const REAL8 yy1 = ab_factor*(aone-(2./9.)*atwo*w0*final_q)*(-a0p*w0 - a0n*wi0) - dhpdi_1*final_mass/(M*dt);

	//h_plus at time step 2
        const REAL8 yy2 = (((aone - (2./9.)*atwo*w0*final_q)*(a0n*cos((M*dt*w0)/final_mass) - 1.*a0p*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(a1n*cos((M*dt*w1)/final_mass) - 1.*a1p*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass))
			- hp_2/ab_factor;

	//h_plus time derivative at time step 2
        const REAL8 yy3 = (((-1.*(aone - (2./9.)*atwo*w0*final_q)*wi0*(a0n*cos((M*dt*w0)/final_mass) - 1.*a0p*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_q)*(-1.*a0p*w0*cos((M*dt*w0)/final_mass) - 1.*a0n*w0*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_q)*wi1*(a1n*cos((M*dt*w1)/final_mass) - 1.*a1p*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(-1.*a1p*w1*cos((M*dt*w1)/final_mass) - 1.*a1n*w1*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)))
			- dhpdi_2*final_mass/(M*dt*ab_factor);

	//h_plus at time step 3
        const REAL8 yy4 = (((aone - (2./9.)*atwo*w0*final_q)*(a0n*cos((2.*M*dt*w0)/final_mass) - 1.*a0p*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(a1n*cos((2.*M*dt*w1)/final_mass) - 1.*a1p*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_q)*(a2n*cos((2.*M*dt*w2)/final_mass) - 1.*a2p*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass))
			- hp_3/ab_factor;

	//h_plus time derivative at time step 3
        const REAL8 yy5 = (((-1.*(aone - (2./9.)*atwo*w0*final_q)*wi0*(a0n*cos((2.*M*dt*w0)/final_mass) - 1.*a0p*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_q)*(-1.*a0p*w0*cos((2.*M*dt*w0)/final_mass) - 1.*a0n*w0*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_q)*wi1*(a1n*cos((2.*M*dt*w1)/final_mass) - 1.*a1p*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(-1.*a1p*w1*cos((2.*M*dt*w1)/final_mass) - 1.*a1n*w1*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w2*final_q)*wi2*(a2n*cos((2.*M*dt*w2)/final_mass) - 1.*a2p*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_q)*(-1.*a2p*w2*cos((2.*M*dt*w2)/final_mass) - 1.*a2n*w2*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)))
			- dhpdi_3*final_mass/(M*dt*ab_factor);

        gsl_vector_set (f, 0, yy0);
        gsl_vector_set (f, 1, yy1);
	gsl_vector_set (f, 2, yy2);
        gsl_vector_set (f, 3, yy3);
        gsl_vector_set (f, 4, yy4);
        gsl_vector_set (f, 5, yy5);

        return GSL_SUCCESS;
	}


static INT4 XLALHGimri_CrossEquations(const gsl_vector * z, void *params, gsl_vector * g) {

	//=============================================================================
	// Used by gsl_multiroots to match cross-polarized plunge waveform onto ringdown.
	// Ringdown waveform consists of (l,m)=(2,+/-2) QNMs, with n=0,1,2 overtones. The
	// RD waveform and its time derivative is matched to plunge waveform at three
	// points across the light ring, considering the n=0 fundamental tone at the first
	// time step, the n=0,1 tones at the second, and the n=0,1,2 tones at the third
	//=============================================================================

	//See definition rparams above for structure information
	REAL8 dt		= ((struct rparams *) params)->dt;
        REAL8 wi0 		= ((struct rparams *) params)->wi0;
        REAL8 w0 		= ((struct rparams *) params)->w0;
        REAL8 wi1 		= ((struct rparams *) params)->wi1;
        REAL8 w1		= ((struct rparams *) params)->w1;
        REAL8 wi2 		= ((struct rparams *) params)->wi2;
        REAL8 w2 		= ((struct rparams *) params)->w2;
        REAL8 final_mass 	= ((struct rparams *) params)->final_mass;
        REAL8 final_q 		= ((struct rparams *) params)->final_q;
        REAL8 Sdotn 		= ((struct rparams *) params)->Sdotn;
        REAL8 M 		= ((struct rparams *) params)->M;
	REAL8 ab_factor 	= ((struct rparams *) params)->ab_factor;
        REAL8 hx_1 		= ((struct rparams *) params)->hx_1;
        REAL8 dhxdi_1 		= ((struct rparams *) params)->dhxdi_1;
        REAL8 hx_2 		= ((struct rparams *) params)->hx_2;
        REAL8 dhxdi_2 		= ((struct rparams *) params)->dhxdi_2;
        REAL8 hx_3 		= ((struct rparams *) params)->hx_3;
        REAL8 dhxdi_3 		= ((struct rparams *) params)->dhxdi_3;

	//Leading amplitude factors
        REAL8 aone		= 1.+ 2.*Sdotn + Sdotn*Sdotn;
        REAL8 atwo		= 2.+ Sdotn - 4.*Sdotn*Sdotn - 3.*Sdotn*Sdotn*Sdotn;

	//Unknown RD constants
        const REAL8 a0c 	= gsl_vector_get (z, 0);
        const REAL8 a0cp 	= gsl_vector_get (z, 1);
        const REAL8 a1c 	= gsl_vector_get (z, 2);
        const REAL8 a1cp 	= gsl_vector_get (z, 3);
        const REAL8 a2c 	= gsl_vector_get (z, 4);
        const REAL8 a2cp 	= gsl_vector_get (z, 5);

	//Functions of the form: (Ringdown Quantity) - (Plunge Quantity). We seek the roots of these functions

	//h_cross at time step 1
	const REAL8 yy0 = ab_factor*a0c*(aone-(2./9.)*atwo*w0*final_q) + hx_1;

	//h_cross time derivative at time step 1
        const REAL8 yy1 = ab_factor*(aone-(2./9.)*atwo*w0*final_q)*(a0cp*w0 - a0c*wi0) + dhxdi_1*final_mass/(M*dt);

	//h_cross at time step 2
        const REAL8 yy2 = (((aone - (2./9.)*atwo*w0*final_q)*(a0c*cos((M*dt*w0)/final_mass) + a0cp*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(a1c*cos((M*dt*w1)/final_mass) + a1cp*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass))
			+ hx_2/ab_factor;

	//h_cross time derivative at time step 2
        const REAL8 yy3 = (((-1.*(aone - (2./9.)*atwo*w0*final_q)*wi0*(a0c*cos((M*dt*w0)/final_mass) + a0cp*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_q)*(a0cp*w0*cos((M*dt*w0)/final_mass) - 1.*a0c*w0*sin((M*dt*w0)/final_mass)))/exp((M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_q)*wi1*(a1c*cos((M*dt*w1)/final_mass) + a1cp*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(a1cp*w1*cos((M*dt*w1)/final_mass) - 1.*a1c*w1*sin((M*dt*w1)/final_mass)))/exp((M*dt*wi1)/final_mass)))
			+ dhxdi_2*final_mass/(M*dt*ab_factor);

	//h_cross at time step 3
        const REAL8 yy4 = (((aone - (2./9.)*atwo*w0*final_q)*(a0c*cos((2.*M*dt*w0)/final_mass) + a0cp*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(a1c*cos((2.*M*dt*w1)/final_mass) + a1cp*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_q)*(a2c*cos((2.*M*dt*w2)/final_mass) + a2cp*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass))
			+ hx_3/ab_factor;

	//h_cross time derivative at time step 3
        const REAL8 yy5 = (((-1.*(aone - (2./9.)*atwo*w0*final_q)*wi0*(a0c*cos((2.*M*dt*w0)/final_mass) + a0cp*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			+ ((aone - (2./9.)*atwo*w0*final_q)*(a0cp*w0*cos((2.*M*dt*w0)/final_mass) - 1.*a0c*w0*sin((2.*M*dt*w0)/final_mass)))/exp((2.*M*dt*wi0)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w1*final_q)*wi1*(a1c*cos((2.*M*dt*w1)/final_mass) + a1cp*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			+ ((aone - (2./9.)*atwo*w1*final_q)*(a1cp*w1*cos((2.*M*dt*w1)/final_mass) - 1.*a1c*w1*sin((2.*M*dt*w1)/final_mass)))/exp((2.*M*dt*wi1)/final_mass)
			- (1.*(aone - (2./9.)*atwo*w2*final_q)*wi2*(a2c*cos((2.*M*dt*w2)/final_mass) + a2cp*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)
			+ ((aone - (2./9.)*atwo*w2*final_q)*(a2cp*w2*cos((2.*M*dt*w2)/final_mass) - 1.*a2c*w2*sin((2.*M*dt*w2)/final_mass)))/exp((2.*M*dt*wi2)/final_mass)))
			+ dhxdi_3*final_mass/(M*dt*ab_factor);

        gsl_vector_set (g, 0, yy0);
        gsl_vector_set (g, 1, yy1);
        gsl_vector_set (g, 2, yy2);
        gsl_vector_set (g, 3, yy3);
        gsl_vector_set (g, 4, yy4);
        gsl_vector_set (g, 5, yy5);

        return GSL_SUCCESS;
	}

static REAL8 XLALHGimri_dLzdr(REAL8 q, REAL8 r) {

	//==============================================
	//Partial derivative (\partial L_z)/(\partial r)
	//for circular geodesic orbits
	//==============================================

	//REAL8 sign;
	//if (q >= 0.) sign = 1.;
	//else sign = -1.;

	REAL8 denom1 = r*r*r-3.*r*r+2.*q*pow(r,3./2.);
	return ((2.*r-q/sqrt(r)-0.5*(r*r-2.*q*sqrt(r)+q*q)*(3.*r*r-6.*r+3.*q*sqrt(r))/denom1)/sqrt(denom1));

	}

static REAL8 XLALHGimri_dFdr(REAL8 E, REAL8 Lz, REAL8 a, REAL8 r) {

	//==================================================
	// Dimensionless (\partial F)/(\partial r), where F=R/(V_t)^2 and
	// R and V_t are the Kerr radial and time potentials.
	// Related to radial coordinate acceleration via
	// (d^2 r)/(dt^2) = (1/2)dFdr
	//==================================================

	REAL8 R = pow(E*(a*a+r*r)-a*Lz,2.) - (r*r-2.*r+a*a)*((Lz-a*E)*(Lz-a*E) + r*r);
	REAL8 Vt = a*(Lz-a*E) + ((r*r+a*a)/(r*r-2.*r+a*a))*(E*(r*r+a*a)-Lz*a);
	REAL8 dRdr = 4.*E*r*(E*(a*a+r*r)-a*Lz) - 2.*(r-1.)*(pow(Lz-a*E,2.)+r*r) - 2.*r*(r*r-2.*r+a*a);
	REAL8 dVtdr = 2.*r*(E*(r*r+a*a)-a*Lz)/(r*r-2.*r+a*a) + (r*r+a*a)*(2.*E*r)/(r*r-2.*r+a*a) - (r*r+a*a)*(E*(r*r+a*a)-a*Lz)*(2.*r-2.)/pow(r*r-2.*r+a*a,2.);

	return ( dRdr/pow(Vt,2.) - 2.*R*dVtdr/pow(Vt,3.) );
	}

static INT4 HGimri_start(REAL8 m, REAL8 M, REAL8 q, REAL8 D, REAL8 Sdotn, REAL8 phi0, REAL8 p0,
	REAL8Sequence *hplus, REAL8Sequence *hcross, REAL8 dt, UINT4 Npts) {

	//====================================================================================================
	// Function which evolves an inspiral trajectory and gravitational waveform for a circular intermediate
	// mass-ratio inspiral, following Huerta & Gair 2011 (see file header above for more info).
	//
	// INPUTS:
	//
	// m: 		Compact object mass (in units of seconds)
	// M: 		Central BH mass (seconds)
	// q: 		Dimensionless reduced spin of the central BH (-1<q<1)
	// D: 		Distance to source (in seconds)
	// Sdotn: 	Cosine of the inclination angle between the detector line of sight and the orbital
	// 		angular momentum vector of the IMRI
	// phi0: 	Initial azimuthal angle
	// p0:		Initial orbital radius
	// hplus:	REAL8TimeSequence object to store plus-polarized GW
	// hcross:	REAL8TimeSequence object to store cross-polarized GW
	// dt:		Dimensionless integration time step
	// Npts:	The maximum number of steps to evolve
	//
	// NOTE:
	//
	// Within this function, masses and source distance are in units of seconds. Central black hole
	// spin is expressed as the dimensionless reduced spin (-1<q<1). Orbital radius and time are dimensionless,
	// with dimensions removed by multiplying by the appropriate factor of total mass (m+M).
	//
	//====================================================================================================

	//Get intrinsic source parameters. Note: masses passed to HGimri_start() in units of seconds.

	REAL8 	mu;		//Reduced mass in seconds
	REAL8 	nu;		//Reduced mass ratio (mu/(m+M))
	REAL8	q_abs;		//Absolute value of spin q
	REAL8 	final_mass;	//Final post-merger mass in units of seconds, from Huerta & Gair (1009.1985) Eq.40
	REAL8	final_q;	//Final post-merger spin, from Huerta & Gair (1009.1985) Eq. 41

	mu		= ((m*M)/(m+M));
	nu		= mu/(m+M);
	q_abs		= fabs(q);
	final_mass 	= (M+m)*(1. + nu*(sqrt(8./9.)-1.) - 0.498*nu*nu);
	final_q		= q_abs - 0.129*q_abs*q_abs*nu - 0.384*q_abs*nu*nu - 2.686*q_abs*nu + 2.*sqrt(3.)*nu - 3.454*nu*nu + 2.353*nu*nu*nu;

	//Prograde (sign = 1) or retrograde (sign = -1)

	REAL8 sign;
	if (q < 0.0)
		sign = -1.0;
	else
		sign = 1.0;

	//Conservative corrections as defined in Huerta & Gair (1009.1985v3) Eq.9

	REAL8 d0	= 1./8.;
	REAL8 d1	= 1975./896.;
	REAL8 d1_5	= -27.*pi/10.;
	REAL8 l1_5	= -191./160.;
	REAL8 d2	= 1152343./451584.;

	//Find ISCO radius and constants

	REAL8	Z1,Z2;
	REAL8 	r_isco;		//Radius at ISCO
	REAL8 	E_isco;		//Energy at ISCO
	REAL8 	L_isco;		//Angular momentum at ISCO
	REAL8 	omega_isco;	//Orbital angular frequency at ISCO
	REAL8 	gamma_isco;	//Lorentz factor d\tau/dt at ISCO

	Z1 		= 1. + pow((1.-q*q),1./3.)*(pow(1.+q, 1./3.)+pow(1.-q, 1./3.));
	Z2 		= sqrt(3.*q*q + Z1*Z1);
	r_isco 		= 3. + Z2 - sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2));
	E_isco		= (1.+q/pow(r_isco,1.5)-2./r_isco)/sqrt(1.+(2.*q)/pow(r_isco,1.5)-3./r_isco);
	L_isco 		= sqrt(r_isco)*(1.+pow(q,2.)/pow(r_isco,2.)-(2.*q)/pow(r_isco,1.5))/sqrt(1.+(2.*q)/pow(r_isco,1.5)-3./r_isco);
	omega_isco 	= (1./(pow(r_isco,3/2.)+q))*(1. + nu*(d0+d1/r_isco+(d1_5+q*l1_5)*pow(r_isco,-1.5)+d2*pow(r_isco,-2.0)));
	gamma_isco 	= sqrt(1. - 3./r_isco + 2.*q/pow(r_isco,1.5))/(1. + q/pow(r_isco,1.5));

	//Find dimensionless constants governing transition behavior (Ori & Thorne /gr-qc/0003032 Equations 3.9,3.18,3.19)

	REAL8 	a_isco;
	REAL8 	b_isco;
	REAL8 	k_isco;
	REAL8 	eps_dot;

	eps_dot = (1.13197 + 0.0713235*q*q + (0.183783*q)/(-2.34484 + 2.15323*q));
	a_isco 	= (3./pow(r_isco,6.)) * (r_isco*r_isco + 2.*( q*q*(E_isco*E_isco-1.) - L_isco*L_isco )*r_isco + 10.*pow(L_isco-q*E_isco,2.));
	b_isco 	= (2./pow(r_isco,4.)) * ( (L_isco - q*q*E_isco*omega_isco)*r_isco - 3.*(L_isco - q*E_isco)*(1. - q*omega_isco) );
	k_isco 	= (32./5.)*pow(omega_isco,7./3.)*eps_dot/gamma_isco;

	//Calculate transition point, plunge point, and light ring radii

	REAL8 	r_match;			//Radius at which inspiral is matched onto transition (Corresponds to T=-1 in Ori & Thorne Eq. 3.23)
	REAL8	L_match;			//Angular momentum at match
	REAL8	r_plunge, E_plunge, L_plunge;	//Conditions at the start of plunge. Corresponds to T=2 (using Ori & Thorne Eq. 3.25)
	REAL8 	r_lightring;			//Radius of Light Ring

	r_match 	= r_isco + sqrt(1.0*pow(nu*b_isco*k_isco,4./5.)/pow(a_isco,6./5.));
	L_match 	= sqrt(r_match)*(1.+pow(q,2.)/pow(r_match,2.)-(2.*q)/pow(r_match,1.5))/sqrt(1.+(2.*q)/pow(r_match,1.5)-3./r_match);
	r_plunge	= r_isco + pow(nu*b_isco*k_isco,2./5.)*pow(a_isco,-3./5.)*(-6./pow(3.412-.75,2.));
	E_plunge	= E_isco - omega_isco*k_isco*pow(a_isco*b_isco*k_isco,-1./5.)*(3.412-0.75)*pow(nu,4./5.);
	L_plunge	= L_isco - k_isco*pow(a_isco*b_isco*k_isco,-1./5.)*(3.412-0.75)*pow(nu,4./5.);
	r_lightring 	= 2.*(1.+cos((2./3.)*acos(-q)));

	//Misc variable declaration

	REAL8 dr;			//Differential radial step corresponding to the time step dt
	REAL8 d2phidt2;			//Second coordinate time derivatives of phi
	REAL8 d2rdt2;			//Second coordinate time derivatives of radius
	INT4 i_match=0.;		//Index at start of transition.
	INT4 i_lightring=0;		//Index at which the CO reaches the light ring
	INT4 i_finalplunge = 0;		//Used to count the final 10 iterations inside the lightring. Plunge halted when i_finalplunge=10

	REAL8 r;					//Instantaneous radial coordinate
	REAL8 r_old;					//Radius at previous time step (used for sanity checks)
	REAL8 drdt=0.;					//Radial velocity
	REAL8 phi;					//Azimuthal position
	REAL8 dphidt;					//Orbital angular velocity
	REAL8 dphidt_old;				//Orbital angular velocity at previous time step (used to numerically compute d2phidt2)
	REAL8 r_K2,drdt_K2, dphidt_K2,d2rdt2_K2=0.;	//For Stage 2 of 4th-order Runge Kutta integration
	REAL8 r_K3,drdt_K3,dphidt_K3,d2rdt2_K3;		//For Stage 3 of RK4 integration
	REAL8 r_K4,drdt_K4,dphidt_K4,d2rdt2_K4;		//For stage 4 of RK4 integration

	//Initialize coordinates
	r 	= p0;
	r_old	= p0;
	phi 	= phi0;
	dphidt	= (1./(pow(r,3/2.)+q))*(1. + nu*(d0+d1/r+(d1_5+q*l1_5)*pow(r,-1.5)+d2*pow(r,-2.0)));
	enum stage currentStage = INSPIRAL;

	//SANITY CHECK: Verify that we're starting outside transition. If not, abort.
	if (r <= r_match) {
		XLALPrintError("XLAL Error: Beginning inside Transition (Stage 2). Specify larger initial frequency f_min\n");
		XLAL_ERROR(XLAL_EDOM);
		}

	//==========================================
	// Iterate through trajectory
	// Stage 1: Adiabatic Inspiral
	// Stage 2: Transition Phase
	// Stage 3: Plunge
	//
	// Evolve trajectory using an 4th order
	// Runge-Kutta integrator
	//==========================================

	for (UINT4 i=0; i<Npts; i++) {

		//Direct to the appropriate iterator
		switch (currentStage) {

			//=============================
			// Stage 1: Adiabatic Inspiral
			//=============================
			case INSPIRAL:

				//Dimensionless radial speed and angular velocity.
				drdt = XLALHGimri_AngMomFlux(q,r,nu)/XLALHGimri_dLzdr(q,r);
				dphidt = (1./(pow(r,3/2.)+q))*(1. + nu*(d0+d1/r+(d1_5+q*l1_5)*pow(r,-1.5)+d2*pow(r,-2.0)));

				//RK4 Step 2
				r_K2 = r + (dt/2.)*drdt;
				drdt_K2 = XLALHGimri_AngMomFlux(q,r_K2,nu)/XLALHGimri_dLzdr(q,r_K2);
				dphidt_K2 = (1./(pow(r_K2,3/2.)+q))*(1. + nu*(d0+d1/r_K2+(d1_5+q*l1_5)*pow(r_K2,-1.5)+d2*pow(r_K2,-2.0)));

				//RK4 Step 3
				r_K3 = r + (dt/2.)*drdt_K2;
				drdt_K3 = XLALHGimri_AngMomFlux(q,r_K3,nu)/XLALHGimri_dLzdr(q,r_K3);
				dphidt_K3 = (1./(pow(r_K3,3/2.)+q))*(1. + nu*(d0+d1/r_K3+(d1_5+q*l1_5)*pow(r_K3,-1.5)+d2*pow(r_K3,-2.0)));

				//RK4 Step 4
				r_K4 = r + (dt)*drdt_K3;
				drdt_K4 = XLALHGimri_AngMomFlux(q,r_K4,nu)/XLALHGimri_dLzdr(q,r_K4);
				dphidt_K4 = (1./(pow(r_K4,3/2.)+q))*(1. + nu*(d0+d1/r_K4+(d1_5+q*l1_5)*pow(r_K4,-1.5)+d2*pow(r_K4,-2.0)));

				//Gravitational wave amplitudes, following Huerta & Gair (1009.1985) Eqs.14,15
				hplus->data[i] = ((4.*mu*pow(dphidt*r,2.))/D) * ((1.+pow(Sdotn,2.))/2.) * cos(2.*phi);
				hcross->data[i] = ((4.*mu*pow(dphidt*r,2.))/D) * Sdotn * sin(2.*phi);

				//If we've hit the transition point, signal move into the transition regime.
				dr = (1./6.)*dt*(drdt+2.*drdt_K2+2.*drdt_K3+drdt_K4);
				if (r+dr < r_match) {

					//Record transition start
					i_match		= i;
					currentStage 	= TRANSITION;

					}

				//Update coordinates
				r_old	= r;
				r 	+= dr;
				phi 	+= (1./6.)*dt*(dphidt+2.*dphidt_K2+2.*dphidt_K3+dphidt_K4);

				break;

			//=============================
			// Stage 2: Transition
			//=============================
			case TRANSITION:

				dphidt_old = dphidt;

				//Update acceleration using H&G (1009.1985) Eq.29. The gamma^2 factor makes this (to very
				//good approximation) a coordinate t derivative. Directly  compute d^2phi/dt^2.
				d2rdt2 = gamma_isco*gamma_isco*(
						- a_isco*(r-r_isco)*(r-r_isco)
						+ b_isco*(L_match-L_isco-nu*k_isco*gamma_isco*(i-i_match)*dt));
				dphidt = (1./(pow(r,3/2.)+q))*(1. + nu*(d0+d1/r+(d1_5+q*l1_5)*pow(r,-1.5)+d2*pow(r,-2.0)));
				d2phidt2 = (dphidt-dphidt_old)/(dt);

				//Step 2
				r_K2 = r+(dt/2.)*drdt;
				drdt_K2 = drdt+(dt/2.)*d2rdt2;
				dphidt_K2 = (1./(pow(r_K2,3/2.)+q))*(1. + nu*(d0+d1/r_K2+(d1_5+q*l1_5)*pow(r_K2,-1.5)+d2*pow(r_K2,-2.0)));
				d2rdt2_K2 = gamma_isco*gamma_isco*(
						- a_isco*(r_K2-r_isco)*(r_K2-r_isco)
						+ b_isco*(L_match-L_isco-nu*k_isco*gamma_isco*((i-i_match)*dt+(dt/2.))));

				//Step 3
				r_K3 = r+(dt/2.)*drdt_K2;
				drdt_K3 = drdt+(dt/2.)*d2rdt2_K2;
				dphidt_K3 = (1./(pow(r_K3,3/2.)+q))*(1. + nu*(d0+d1/r_K3+(d1_5+q*l1_5)*pow(r_K3,-1.5)+d2*pow(r_K3,-2.0)));
				d2rdt2_K3 = gamma_isco*gamma_isco*(
						- a_isco*(r_K3-r_isco)*(r_K3-r_isco)
						+ b_isco*(L_match-L_isco-nu*k_isco*gamma_isco*((i-i_match)*dt+(dt/2.))));

				//Step 4
				r_K4 = r+(dt)*drdt_K3;
				drdt_K4 = drdt+(dt/2.)*d2rdt2_K3;
				dphidt_K4 = (1./(pow(r_K4,3/2.)+q))*(1. + nu*(d0+d1/r_K4+(d1_5+q*l1_5)*pow(r_K4,-1.5)+d2*pow(r_K4,-2.0)));
				d2rdt2_K4 = gamma_isco*gamma_isco*(
						- a_isco*(r_K4-r_isco)*(r_K4-r_isco)
						+ b_isco*(L_match-L_isco-nu*k_isco*gamma_isco*((i-i_match)*dt+(dt))));

				//GW amplitudes. Huerta & Gair (1009.1985) Eq. 37 and 38
				hplus->data[i] = (mu/(2.*D))*(
						(1. - 2.*(2.*Sdotn*Sdotn-1.)*pow(cos(phi),2.) - 3.*cos(2.*phi))*drdt*drdt
						+ (3. + (2.*Sdotn*Sdotn-1.))*(2.*cos(2.*phi)*dphidt*dphidt + sin(2.*phi)*d2phidt2)*r*r
						+ (4.*(3.+(2.*Sdotn*Sdotn-1.))*sin(2.*phi)*dphidt*drdt
							+ (1.-2.*(2.*Sdotn*Sdotn-1.)*pow(cos(phi),2.)-3.*cos(2.*phi))*d2rdt2)*r);

				hcross->data[i] = (-2.*mu/D)*Sdotn*( sin(2.*(phi))*drdt*drdt
						+ r*r*(-2.*sin(2.*(phi))*dphidt*dphidt + cos(2.*( phi))*d2phidt2)
						+ r*(4.*cos(2.*(phi))*dphidt*drdt + sin(2.*(phi))*d2rdt2));

				//Check if we've reached reached plunge
				dr = (1./6.)*dt*(drdt+2.*drdt_K2+2.*drdt_K3+drdt_K4);
				if (r+dr < r_plunge) {

					//Recompute dphidt using geodesic equation to avoid discontinuous d^2phi/dt^2
					dphidt = (L_plunge*(r-2.)+2.*E_plunge*q)/(E_plunge*(r*r*r+(2.+r)*q*q) - 2.*q*L_plunge);

					//Record plunge start
					currentStage = PLUNGE;

					}

				//Iterate
				r_old	= r;
				r	+= dr;
				phi 	+= (1./6.)*dt*(dphidt+2.*dphidt_K2+2.*dphidt_K3+dphidt_K4);
				drdt 	+= (1./6.)*dt*(d2rdt2+2.*d2rdt2_K2+2.*d2rdt2_K3+d2rdt2_K4);

				break;

			//====================================
			// Stage 3: Plunge
			//====================================

			//If inside LR, iterate i_finalplunge and move on to PLUNGE case below
			case FINALPLUNGE:

				//Update counter to track # of iterations inside LR
				i_finalplunge ++;

				//After 10 iterations inside light ring, stop integration.
				if (i_finalplunge == 10) {

					i = Npts;
					break;

					}

			//In all cases...
			case PLUNGE:

				//RK1
				dphidt_old = dphidt;
				dphidt = (L_plunge*(r-2.)+2.*E_plunge*q)/(E_plunge*(r*r*r+(2.+r)*q*q) - 2.*q*L_plunge);
				d2phidt2 = (dphidt-dphidt_old)/(dt);
				d2rdt2 	= (1./2.)*XLALHGimri_dFdr(E_plunge,L_plunge,q,r);

				//RK2
				r_K2 = r + (dt/2.)*drdt;
				drdt_K2 = drdt + (dt/2.)*d2rdt2;
				dphidt_K2 = (L_plunge*(r_K2-2.)+2.*E_plunge*q)/(E_plunge*(r_K2*r_K2*r_K2+(2.+r_K2)*q*q) - 2.*q*L_plunge);
				d2rdt2 = (1./2.)*XLALHGimri_dFdr(E_plunge,L_plunge,q,r_K2);

				//RK3
				r_K3 = r + (dt/2.)*drdt_K2;
				drdt_K3 = drdt + (dt/2.)*d2rdt2_K2;
				dphidt_K3 = (L_plunge*(r_K3-2.)+2.*E_plunge*q)/(E_plunge*(r_K3*r_K3*r_K3+(2.+r_K3)*q*q) - 2.*q*L_plunge);
				d2rdt2_K3 = (1./2.)*XLALHGimri_dFdr(E_plunge,L_plunge,q,r_K3);

				//RK4
				r_K4 = r + (dt)*drdt_K3;
				drdt_K4 = drdt + (dt)*d2rdt2_K3;
				dphidt_K4 = (L_plunge*(r_K4-2.)+2.*E_plunge*q)/(E_plunge*(r_K4*r_K4*r_K4+(2.+r_K4)*q*q) - 2.*q*L_plunge);
				d2rdt2_K4 = (1./2.)*XLALHGimri_dFdr(E_plunge,L_plunge,q,r_K4);

				//GW amplitudes
				hplus->data[i] = (mu/(2.*D))*((1. - 2.*(2.*Sdotn*Sdotn-1.)*pow(cos(phi),2.) - 3.*cos(2.*phi))*drdt*drdt
						+ (3. + (2.*Sdotn*Sdotn-1.))*(2.*cos(2.*phi)*dphidt*dphidt + sin(2.*phi)*d2phidt2)*r*r
						+ (4.*(3.+(2.*Sdotn*Sdotn-1.))*sin(2.*phi)*dphidt*drdt
							+ (1.-2.*(2.*Sdotn*Sdotn-1.)*pow(cos(phi),2.)-3.*cos(2.*phi))*d2rdt2)*r);

				hcross->data[i] = (-2.*mu/D)*Sdotn*( sin(2.*(phi))*drdt*drdt
						+ r*r*(-2.*sin(2.*(phi))*dphidt*dphidt + cos(2.*( phi))*d2phidt2)
						+ r*(4.*cos(2.*(phi))*dphidt*drdt + sin(2.*(phi))*d2rdt2));

				//Update coordinates
				r_old	= r;
				r	+= (1./6.)*dt*(drdt+2.*drdt_K2+2.*drdt_K3+drdt_K4);
				phi	+= (1./6.)*dt*(dphidt+2.*dphidt_K2+2.*dphidt_K3+dphidt_K4);
				drdt	+= (1./6.)*dt*(d2rdt2+2.*d2rdt2_K2+2.*d2rdt2_K3+d2rdt2_K4);

				//check to see if we've passed light ring
				if ((r < r_lightring) && (currentStage == PLUNGE)) {

					//If so, initiate end of plunge
					currentStage = FINALPLUNGE;
					i_lightring = i;

					}

				break;

			}

			//===============
			// SANITY CHECKS
			//===============

			//If r is negative or NaN, raise error.
			if (r != r || r < 0.0) {
				XLALPrintError("XLAL Error: Radius is negative or NaN. ");
				if (currentStage == INSPIRAL) XLALPrintError("Halting INSPIRAL\n.");
				if (currentStage == TRANSITION) XLALPrintError("Halting TRANSITION.\n");
				if (currentStage == PLUNGE) XLALPrintError("Halting PLUNGE.\n");
				if (currentStage == FINALPLUNGE) {
					if (i_finalplunge <= 10)
						XLALPrintError("Failed to complete 10 steps inside light ring. Use smaller time step dt. ");
					XLALPrintError("Halting FINALPLUNGE.\n");
					}
				XLAL_ERROR(XLAL_EFAILED);
				}

			//If r has increased since the last time step, raise error.
			if (r > r_old) {
				XLALPrintError("XLAL Error: Increasing radius. ");
				if (currentStage == INSPIRAL) XLALPrintError("Halting INSPIRAL\n.");
				if (currentStage == TRANSITION) XLALPrintError("Halting TRANSITION.\n");
				if (currentStage == PLUNGE) XLALPrintError("Halting PLUNGE.\n");
				if (currentStage == FINALPLUNGE) {
					if (i_finalplunge <= 10)
						XLALPrintError("Failed to complete 10 steps inside light ring. Use smaller time step dt. ");
					XLALPrintError("Halting FINALPLUNGE.\n");
					}
				XLAL_ERROR(XLAL_EFAILED);
				}

			//If we've reached the max array length and failed to plunge, either something has gone wrong numerically,
			//or Npts is too small. In the latter case, choose smaller time steps or increase Npts in XLALHGimri_generator() below
			if (i == Npts && currentStage != FINALPLUNGE) {
				XLALPrintError("XLAL Error: Reached Npts before finishing plunge. Binary either failed to plunge, or Npts is too small.\n");
				XLAL_ERROR(XLAL_EFAILED);
				}

		}


	//=====================================================================
	// Get ringdown frequencies by interpolating data presented in Table 2
	// of Berti et al 2006 (gr-qc/0512160). See Huerta & Gair for details.
	//=====================================================================

	INT4 N = 11;

	//Final q values
	REAL8 fq[11] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98};

	//Real and imaginary dimensionless frequencies for mode (l=2,m=2,n=0). From Berti et. al.
	REAL8 w0_arr[11] 	= {0.3737, 0.3870, 0.4021, 0.4195, 0.4398, 0.4641, 0.4940, 0.5326, 0.5860, 0.6716, 0.8254};
	REAL8 wi0_arr[11] 	= {0.0890, 0.0887, 0.0883, 0.0877, 0.0869, 0.0856, 0.0838, 0.0808, 0.0756, 0.0649, 0.0386};

	//Real and imaginary dimensionless frequencies for mode (l=2,m=2,n=1). From Berti et. al.
	REAL8 w1_arr[11] 	= {0.3467, 0.3619, 0.3790, 0.3984, 0.4208, 0.4474, 0.4798, 0.5212, 0.5779, 0.6677, 0.8249};
	REAL8 wi1_arr[11]	= {0.2739, 0.2725, 0.2705, 0.2680, 0.2647, 0.2602, 0.2538, 0.2442, 0.2281, 0.1953, 0.1159};

	//Real and imaginary dimensionless frequencies for mode (l=2,m=2,n=2). From Berti et. al.
	REAL8 w2_arr[11] 	= {0.3011, 0.3192, 0.3393, 0.3619, 0.3878, 0.4179, 0.4542, 0.4999, 0.5622, 0.6598, 0.8238};
	REAL8 wi2_arr[11]	= {0.4783, 0.4735, 0.4679, 0.4613, 0.4533, 0.4433, 0.4303, 0.4123, 0.3839, 0.3275, 0.1933};

	//Set up interpolator
	const gsl_interp_type *bn = gsl_interp_cspline;
	gsl_interp_accel *acc = gsl_interp_accel_alloc();

	//Define a bunch of interpolation splines
	gsl_spline *spline0 = gsl_spline_alloc(bn,N);
	gsl_spline *spline1 = gsl_spline_alloc(bn,N);
	gsl_spline *spline2 = gsl_spline_alloc(bn,N);
	gsl_spline *spline3 = gsl_spline_alloc(bn,N);
	gsl_spline *spline4 = gsl_spline_alloc(bn,N);
	gsl_spline *spline5 = gsl_spline_alloc(bn,N);

	//Initialize splines
	gsl_spline_init(spline0, fq, w0_arr, N);
	gsl_spline_init(spline1, fq, wi0_arr, N);
	gsl_spline_init(spline2, fq, w1_arr, N);
	gsl_spline_init(spline3, fq, wi1_arr, N);
	gsl_spline_init(spline4, fq, w2_arr, N);
	gsl_spline_init(spline5, fq, wi2_arr, N);

	//Get the real and imaginary frequencies for our final BH
	REAL8 w0, wi0;
	REAL8 w1, wi1;
	REAL8 w2, wi2;

	w0 	= gsl_spline_eval(spline0, final_q, acc);
	wi0 	= gsl_spline_eval(spline1, final_q, acc);
	w1	= gsl_spline_eval(spline2, final_q, acc);
	wi1	= gsl_spline_eval(spline3, final_q, acc);
	w2	= gsl_spline_eval(spline4, final_q, acc);
	wi2	= gsl_spline_eval(spline5, final_q, acc);

	//Free memory
	gsl_spline_free(spline0);
	gsl_spline_free(spline1);
	gsl_spline_free(spline2);
	gsl_spline_free(spline3);
	gsl_spline_free(spline4);
	gsl_spline_free(spline5);


	//==================================================================
	// Build interpolation functions for h_plus and h_cross at 20 time
	// steps centered at r_lightring.
	//==================================================================

	//Define arrays to hold the final 20 values of hplus, hcross, and their derivatives during plunge.
	//Note, derivatives are technically dh/di, not dh/dt. Get the latter using dh/dt = (dh/di)/dt
	REAL8 i_interp[20];
	REAL8 hp_interp[20];
	REAL8 hx_interp[20];

	//Unknown constant amplitudes in ringdown equation.
	REAL8 a0n, a0p;
	REAL8 a1n, a1p;
	REAL8 a2n, a2p;
	REAL8 a0c, a0cp;
	REAL8 a1c, a1cp;
	REAL8 a2c, a2cp;

	//Read out the final twenty values of hplus[] and hcross[].
	//The lightring bisects i_lightring (i_interp[10]) and i_lightring-1 (i_interp[9])
	for (INT4 i=0; i<20; i++) {
		i_interp[i]	= i+(i_lightring-10);
		hp_interp[i] 	= hplus->data[i+(i_lightring-10)]*D;
		hx_interp[i] 	= hcross->data[i+(i_lightring-10)]*D;
		}

	//Make some splines!
	gsl_spline *spline_hp = gsl_spline_alloc(gsl_interp_cspline, 20);	//To interpolate across hplus
	gsl_spline *spline_hx = gsl_spline_alloc(gsl_interp_cspline, 20);	//To interpolate across hcross
	gsl_spline_init(spline_hp, i_interp, hp_interp, 20);
	gsl_spline_init(spline_hx, i_interp, hx_interp, 20);

	//Define hp, hx, and their first, second, and third derivatives with respect to i at three steps across the light radius.
	REAL8 hp_1, dhpdi_1;
	REAL8 hx_1, dhxdi_1;
	REAL8 hp_2, dhpdi_2;
	REAL8 hx_2, dhxdi_2;
	REAL8 hp_3, dhpdi_3;
	REAL8 hx_3, dhxdi_3;

	//Evaluate first point at i_interp[9] = i_lightring-1
	hp_1		= gsl_spline_eval(spline_hp, i_interp[9], acc);
	dhpdi_1		= gsl_spline_eval_deriv(spline_hp, i_interp[9], acc);
	hx_1		= gsl_spline_eval(spline_hx, i_interp[9], acc);
	dhxdi_1		= gsl_spline_eval_deriv(spline_hx, i_interp[9], acc);

	//Evaluate second point at i_interp[10] = i_lightring
	hp_2		= gsl_spline_eval(spline_hp, i_interp[10], acc);
	dhpdi_2		= gsl_spline_eval_deriv(spline_hp, i_interp[10], acc);
	hx_2		= gsl_spline_eval(spline_hx, i_interp[10], acc);
	dhxdi_2		= gsl_spline_eval_deriv(spline_hx, i_interp[10], acc);

	//Evaluate third point at i_interp[11] = i_lightring+1
	hp_3		= gsl_spline_eval(spline_hp, i_interp[11], acc);
	dhpdi_3		= gsl_spline_eval_deriv(spline_hp, i_interp[11], acc);
	hx_3		= gsl_spline_eval(spline_hx, i_interp[11], acc);
	dhxdi_3		= gsl_spline_eval_deriv(spline_hx, i_interp[11], acc);

	//Free memory
	gsl_spline_free(spline_hp);
	gsl_spline_free(spline_hx);

	//==============================================================
	// Match plunging waveform onto QNM ringdown waveform by
	// numerically finding the QNM amplitudes which match the
	// waveform amplitudes and derivatives at the three points above.
	//
	// See comments in declarations of XLALHGimri_Plus/CrossEquations()
	// for more detail.
	// =============================================================

	//Create an instance of rparams, populate
	struct rparams p;
	p.dt		= dt;
	p.M		= (M+m);
	p.Sdotn		= Sdotn;
	p.final_mass	= final_mass;
	p.final_q	= final_q;
	p.w0		= w0;
	p.w1		= w1;
	p.w2		= w2;
	p.wi0		= wi0;
	p.wi1		= wi1;
	p.wi2		= wi2;
	p.ab_factor	= (1./8.)*sqrt(5./pi)*final_mass;
	p.hp_1		= hp_1;
	p.hp_2		= hp_2;
	p.hp_3		= hp_3;
	p.dhpdi_1	= dhpdi_1;
	p.dhpdi_2	= dhpdi_2;
	p.dhpdi_3	= dhpdi_3;
	p.hx_1		= hx_1;
	p.hx_2		= hx_2;
	p.hx_3		= hx_3;
	p.dhxdi_1	= dhxdi_1;
	p.dhxdi_2	= dhxdi_2;
	p.dhxdi_3	= dhxdi_3;

	//Initialize root solvers
	const gsl_multiroot_fsolver_type *T;
	const gsl_multiroot_fsolver_type *Q;
	gsl_multiroot_fsolver *plus_solver;
	gsl_multiroot_fsolver *cross_solver;
	INT4 statusp, statusc;
	size_t iter = 0;
	const size_t nmatch = 6;

	//Define target functions to be XLALHGimri_Plus/CrossEquations
	gsl_multiroot_function f_plus = {&XLALHGimri_PlusEquations, nmatch, &p};
	gsl_multiroot_function f_cross = {&XLALHGimri_CrossEquations, nmatch, &p};

	//Define two 6D vectors, and feed initial guesses
	REAL8 x_init[6] = {-0.392254,  4.575194, -0.870431,  5.673678,  0.979356, -3.637467};
	REAL8 z_init[6] = {-0.392254,  4.575194, -0.870431,  5.673678,  0.979356, -3.637467};
	gsl_vector *xmr = gsl_vector_alloc(nmatch);
	gsl_vector *z 	= gsl_vector_alloc(nmatch);

	for (size_t i=0; i<nmatch; i++) {
		gsl_vector_set(xmr, i, x_init[i]);
		gsl_vector_set(z, i, z_init[i]);
		}

	//Set up solvers
	T = gsl_multiroot_fsolver_hybrids;
	Q = gsl_multiroot_fsolver_hybrids;
	plus_solver = gsl_multiroot_fsolver_alloc(T, 6);
	cross_solver = gsl_multiroot_fsolver_alloc(Q, 6);
	gsl_multiroot_fsolver_set(plus_solver, &f_plus, xmr);
	gsl_multiroot_fsolver_set(cross_solver, &f_cross, z);

	REAL8 fminval = 1.e-5;
	size_t MaxITS = 200000;

	//Find root of XLALHGimri_PlusEquations
	do {
		iter++;
		statusp = gsl_multiroot_fsolver_iterate(plus_solver);
		if (statusp) break;
		statusp = gsl_multiroot_test_residual(plus_solver->f, fminval);
		}
	while (statusp == GSL_CONTINUE && iter < MaxITS);

	//Find root of XLALHGimri_CrossEquations
	iter = 0;
	do {
		iter++;
		statusc = gsl_multiroot_fsolver_iterate(cross_solver);
		if (statusc) break;
		statusc = gsl_multiroot_test_residual(cross_solver->f, fminval);
		}
	while (statusc == GSL_CONTINUE && iter < MaxITS);

	if (statusc != GSL_SUCCESS) {
		XLALPrintError("XLAL Error: Ringdown rootfinding failed.\n");
		XLAL_ERROR(XLAL_EFAILED);
		}

	//Read out fit parameters
	a0n	= gsl_vector_get(plus_solver->x,0);
	a0p	= gsl_vector_get(plus_solver->x,1);
	a1n	= gsl_vector_get(plus_solver->x,2);
	a1p	= gsl_vector_get(plus_solver->x,3);
	a2n	= gsl_vector_get(plus_solver->x,4);
	a2p	= gsl_vector_get(plus_solver->x,5);
	a0c	= gsl_vector_get(cross_solver->x,0);
	a0cp	= gsl_vector_get(cross_solver->x,1);
	a1c	= gsl_vector_get(cross_solver->x,2);
	a1cp	= gsl_vector_get(cross_solver->x,3);
	a2c	= gsl_vector_get(cross_solver->x,4);
	a2cp	= gsl_vector_get(cross_solver->x,5);

	//Free memory
	gsl_multiroot_fsolver_free(plus_solver);
	gsl_multiroot_fsolver_free(cross_solver);
	gsl_vector_free(xmr);
	gsl_vector_free(z);

	//=====================================================================================
	// Recompute wave amplitudes from (i_lightring-1), the first matching point, onwards
	//=====================================================================================

	REAL8 aone, atwo;
	REAL8 hp0, hp1, hp2;
	REAL8 hc0, hc1, hc2;
	REAL8 ringdown_amp;
	REAL8 hp_ringdown, hc_ringdown;

	aone = 1.+ 2.*Sdotn + Sdotn*Sdotn;
	atwo = 2.+ Sdotn - 4.*Sdotn*Sdotn - 3.*Sdotn*Sdotn*Sdotn;
	ringdown_amp = (1./8.)*sqrt(5./pi)*(final_mass/D);

	for (UINT4 i=i_lightring; i<Npts; i++) {

		//NOTE: final_q and the ringdown frequencies wi are nondimensionalized by final_mass. Their direct product is therefore already
		//dimensionless. When multiplying the wi's by dt, however, factors of (M/final_mass), however, are needed to properly rescale the wi's.
		hp0 = exp(-(i-(i_lightring-1.))*dt*(wi0*(m+M)/final_mass))*( aone-(2./9.)*final_q*w0*atwo )*( a0n*cos(w0*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)) - a0p*sin(w0*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)));
		hp1 = exp(-(i-(i_lightring-1.))*dt*(wi1*(m+M)/final_mass))*( aone-(2./9.)*final_q*w1*atwo )*( a1n*cos(w1*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)) - a1p*sin(w1*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)));
		hp2 = exp(-(i-(i_lightring-1.))*dt*(wi2*(m+M)/final_mass))*( aone-(2./9.)*final_q*w2*atwo )*( a2n*cos(w2*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)) - a2p*sin(w2*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)));
		hc0 = exp(-(i-(i_lightring-1.))*dt*(wi0*(m+M)/final_mass))*( aone-(2./9.)*final_q*w0*atwo )*( a0c*cos(w0*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)) + a0cp*sin(w0*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)));
		hc1 = exp(-(i-(i_lightring-1.))*dt*(wi1*(m+M)/final_mass))*( aone-(2./9.)*final_q*w1*atwo )*( a1c*cos(w1*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)) + a1cp*sin(w1*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)));
		hc2 = exp(-(i-(i_lightring-1.))*dt*(wi2*(m+M)/final_mass))*( aone-(2./9.)*final_q*w2*atwo )*( a2c*cos(w2*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)) + a2cp*sin(w2*(i-(i_lightring-1.))*(dt*(m+M)/final_mass)));

		//Get hplus and hcross
		hp_ringdown	= (hp0+hp1+hp2);
		hc_ringdown	= -(hc0+hc1+hc2);

		//Record wave amplitudes
		hplus->data[i] 	= ringdown_amp*hp_ringdown;
		hcross->data[i] = ringdown_amp*hc_ringdown;

		//Halt when wave amplitudes become sufficiently low
		if (fabs(hplus->data[i]) < 1.e-30) {

			hplus->length=i-1;	//Redefine REAL8TimeSeries length to be the current index.
			hcross->length=i-1;
			i = Npts;		//Halt.

			}

		}

	return 0;

	}


/**
 * @addtogroup LALSimInspiralHGimri_c
 * @brief Routines for generating the Huerta-Gair Intermediate-Mass-Ratio Inspiral model.
 * @{
 */

//===================================================
// Generator function for the Huerta-Gair IMRI model.
//===================================================

INT4 XLALHGimriGenerator(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 phi0,			//Initial phi
	REAL8 dt,			//Integration time step (in seconds)
	REAL8 m1,			//BH mass (passed in kg)
	REAL8 m2,			//CO mass (passed in kg)
	REAL8 f_min,			//Initial frequency (passed in Hz)
	REAL8 r,			//Distance to system (passed in meters)
	REAL8 inc,			//Inclination angle between detector line of sight n and central BH spin axis
	REAL8 s1z			//Central BH reduced spin
	) {

	//Parse input, adjusting dimensiones
	REAL8 q = s1z;				//Rename spin
	REAL8 Sdotn = cos(inc);			//Cosine(inclination)
	REAL8 M = m1/LAL_MSUN_SI;		//Units of solar masses
	REAL8 m = m2/LAL_MSUN_SI;		//Units of solar masses
	REAL8 dist = r/(LAL_PC_SI*1.e9);	//Units of Gpc

	//SANITY CHECKS
	if (f_min < 0.) {
		XLALPrintError("XLAL Error: Minimum frequency is %f. f_min must be greater than 0\n",f_min);
		XLAL_ERROR(XLAL_EDOM);
		}

	if (q >= 1 || q<=-1) {
		XLALPrintError("XLAL Error: Magnitude of BH spin %f must be less than 1.\n",q);
		XLAL_ERROR(XLAL_EDOM);
		}

	if (m/M > 1./10.) {
		XLALPrintWarning("XLAL Warning: Mass-ratio is %f. Code likely inaccurate for mass ratios m2/m1 > 1/10\n",m/M);
		}

	if (dist <= 0) {
		XLALPrintError("XLAL Error: Distance %f Mpc to source must be greater than 0 Mpc.\n",dist);
		XLAL_ERROR(XLAL_EDOM);
		}

	//==========================================================================
	//Numerically solve for initial radius corresponding to initial gravitational
	//wave frequency of f_min (or an initial orbital frequency of twice f_min).
	//==========================================================================

	//Construct a gsl rootfinder.
	INT4 status;
	size_t iter=0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	//Populate pParams object
	struct pParams p0Params;
	p0Params.m = m;
	p0Params.M = M;
	p0Params.q = q;
	p0Params.f_min = f_min;

	//Define target function
	gsl_function F;
	F.function = &XLALHGimri_initialP;
	F.params = &p0Params;

	//Initialize p0 and the search bounds.
	REAL8 p0 = 0.;
	REAL8 p_high = 500.;
	REAL8 p_low = 0.;

	//Set a reasonable lower limit by equating p_low to the ISCO radius.
	if (q==0.) p_low = 6.;
	else {
		REAL8 Z1 = 1. + pow((1.-q*q),1./3.)*(pow(1.+q, 1./3.)+pow(1.-q, 1./3.));
		REAL8 Z2 = sqrt(3.*q*q + Z1*Z1);
		p_low = 3. + Z2 - (q/fabs(q))*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2));
		}

	//Find p0
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, p_low, p_high);

	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		p0 = gsl_root_fsolver_root(s);
		p_low = gsl_root_fsolver_x_lower(s);
		p_high = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval (p_low, p_high, 0, 0.001);
		}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);

	if (status != GSL_SUCCESS) {
		XLALPrintError("XLAL Error: Rootfinding failed. Could not find inital radius ISCO < r < 500M satisfying specified initial frequency.\n");
		XLAL_ERROR(XLAL_EFAILED); //"generic failure"
		}

	if (p0 >= 20.) {
		XLALPrintWarning("XLAL Warning: Beginning at large initial radius p0 = %fM. Code will require long runtime to compute waveform.\n",p0);
		}

	//==============================================================
	//Create REAL8TimeSeries to hold waveform data, start simulation
	//==============================================================

	REAL8 fDyne = 0.0;
	size_t length = 1000/dt;
	static LIGOTimeGPS epoch;
	*hplus = XLALCreateREAL8TimeSeries("h_plus",&epoch,fDyne,dt,&lalStrainUnit,length);
	*hcross = XLALCreateREAL8TimeSeries("h_cross",&epoch,fDyne,dt,&lalStrainUnit,length);
	if (*hplus == NULL || *hcross == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	HGimri_start(m*Msun_sec,M*Msun_sec,q,dist*GPC_sec,Sdotn,phi0,p0,(*hplus)->data,(*hcross)->data,dt/((m+M)*Msun_sec),length);

	return 0;

	}

/** @} */
