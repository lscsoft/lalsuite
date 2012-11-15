#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lal/LALSimBlackHoleRingdown.h>

#define EPS LAL_REAL4_EPS
#define TINY LAL_REAL8_MIN
#define MAXITER 16384

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


struct LALSimBlackHoleRingdownModeLeaver { double a; int l; int m; int s; COMPLEX16 A; COMPLEX16 omega; };


/* Expansion coefficients for Schwarzschild modes. */
/* Equation 8 of Leaver (1985). */
static int XLALSimBlackHoleRingdownModeSphericalCoefficientsLeaver(COMPLEX16 *alp, COMPLEX16 *bet, COMPLEX16 *gam, double UNUSED a, int l, int UNUSED m, int n, int s, COMPLEX16 UNUSED A, COMPLEX16 omega)
{
	COMPLEX16 rho = -I*omega;
	int epsilon = s*s - 1;
	*alp = n*n + (2.0*rho + 2.0)*n + 2.0*rho + 1.0;
	*bet = -(2.0*n*n + (8.0*rho + 2.0)*n + 8.0*rho*rho + 4.0*rho + l*(l + 1.0) - epsilon);
	*gam = n*n + 4.0*rho*n + 4.0*rho*rho - epsilon - 1.0;
	return 0;
}


/* Angular expansion coefficients for Kerr modes. */
/* Equation 20 of Leaver (1985). */
static int XLALSimBlackHoleRingdownModeAngularCoefficientsLeaver(COMPLEX16 *alp, COMPLEX16 *bet, COMPLEX16 *gam, double a, int UNUSED l, int m, int n, int s, COMPLEX16 A, COMPLEX16 omega)
{
	double k1 = 0.5*abs(m - s);
	double k2 = 0.5*abs(m + s);
	*alp = -2.0*(n + 1.0)*(n + 2.0*k1 + 1.0);
	*bet = n*(n - 1.0) + 2.0*n*(k1 + k2 + 1.0 - 2.0*a*omega);
	*bet -= 2.0*a*omega*(2.0*k1 + s + 1.0) - (k1 + k2)*(k1 + k2 + 1.0);
	*bet -= a*a*omega*omega + s*(s + 1.0) + A;
	*gam = 2.0*a*omega*(n + k1 + k2 + s);
	return 0;
}


/* Radial expansion coefficients for Kerr modes. */
/* Equations 25 and 26 of Leaver (1985). */
static int XLALSimBlackHoleRingdownModeRadialCoefficientsLeaver(COMPLEX16 *alp, COMPLEX16 *bet, COMPLEX16 *gam, double a, int UNUSED l, int m, int n, int s, COMPLEX16 A, COMPLEX16 omega)
{
	COMPLEX16 c0, c1, c2, c3, c4;
	double b = sqrt(1.0 - 4.0*a*a);
	c0 = 1.0 - s - I*omega - 2.0*I*(0.5*omega - a*m)/b;
	c1 = -4.0 + 2.0*I*omega*(2.0 + b) + 4.0*I*(0.5*omega - a*m)/b;
	c2 = s + 3.0 - 3.0*I*omega - 2.0*I*(0.5*omega - a*m)/b;
	c3 = omega*omega*(4.0 + 2.0*b - a*a) - 2.0*a*m*omega - s - 1.0 + (2.0 + b)*I*omega - A
		+ (4.0*omega + 2.0*I)*(0.5*omega - a*m)/b;
	c4 = s + 1.0 - 2.0*omega*omega - (2.0*s + 3.0)*I*omega
		- (4.0*omega + 2.0*I)*(0.5*omega - a*m)/b;
	*alp = n*n + (c0 + 1.0)*n + c0;
	*bet = -2.0*n*n + (c1 + 2.0)*n + c3;
	*gam = n*n + (c2 - 3.0)*n + c4 - c2 + 2.0;
	return 0;
}


/* Evaluates the continued fractions in            */
/* Equations (13), (21), or (27) of Leaver (1985). */
/* TODO: if we ever want overtones, we will need   */
/* generalize this so that it can solve inversions */
/* such as Equation (14) of Leaver (1985)....      */
/* Uses the modified Lentz's method (see Numerical */
/* Recipes).                                       */
static COMPLEX16 XLALSimBlackHoleRingdownModeEigenvalueEvaluateContinuedFractionLeaver(double a, int l, int m, int s, COMPLEX16 A, COMPLEX16 omega, int (*coef)(COMPLEX16 *, COMPLEX16 *, COMPLEX16 *, double, int, int, int, int, COMPLEX16, COMPLEX16))
{
	int n = 0;
	COMPLEX16 alp, alpsv, bet, gam;
	COMPLEX16 afac, bfac;
	COMPLEX16 f, fsv, C, D, Delta;

	alpsv = 0;
	coef(&alp, &bet, &gam, a, l, m, n, s, A, omega);
	afac  = -alpsv*gam;
	alpsv = alp;
	bfac  = bet;
	f = bfac;
	if (cabs(f) < TINY)
		f = TINY;
	C = f;
	D = 0;

	while (n++ < MAXITER) {
		coef(&alp, &bet, &gam, a, l, m, n, s, A, omega);
		afac  = -alpsv*gam;
		alpsv = alp;
		bfac  = bet;
		D = bfac + afac*D;
		if (cabs(D) < TINY)
			D = TINY;
		C = bfac + afac/C;
		if (cabs(C) < TINY)
			C = TINY;
		D = 1.0/D;
		Delta = C*D;
		fsv = f;
		f *= Delta;
		if (cabs(f - fsv) < EPS)
			return f;
	}

	/* only get here if MAXITER is exceeded */
	XLAL_ERROR(XLAL_EMAXITER);
}


/* Residual function for solving the continued fraction equation */
/* Equation (13) of Leaver (1985). */
static int XLALSimBlackHoleRingdownModeSchwarzschildEigenvalueSolveResid(const gsl_vector *x, void *params, gsl_vector *f)
{
	struct LALSimBlackHoleRingdownModeLeaver *p = params;
	COMPLEX16 A = 0.0, omega;
	COMPLEX16 cf;
	int errnum;
	omega = gsl_vector_get(x, 0) + I*gsl_vector_get(x, 1);
	XLAL_TRY(cf = XLALSimBlackHoleRingdownModeEigenvalueEvaluateContinuedFractionLeaver(p->a, p->l, p->m, p->s, A, omega, XLALSimBlackHoleRingdownModeSphericalCoefficientsLeaver), errnum);
	if (errnum)
		XLAL_ERROR(XLAL_EFUNC);

	gsl_vector_set(f, 0, creal(cf));
	gsl_vector_set(f, 1, cimag(cf));
	return 0;
}


/* Residual function for solving the continued fraction equations */
/* Equation (21) and (27) of Leaver (1985). */
/* Note: simultaneously solves both of these equations together. */
static int XLALSimBlackHoleRingdownModeKerrEigenvalueSolveResid(const gsl_vector *x, void *params, gsl_vector *f)
{
	struct LALSimBlackHoleRingdownModeLeaver *p = params;
	COMPLEX16 A, omega;
	COMPLEX16 cf1, cf2;
	int errnum;
	A = gsl_vector_get(x,0) + I*gsl_vector_get(x,1);
	omega = gsl_vector_get(x,2) + I*gsl_vector_get(x,3);
	XLAL_TRY(cf1 = XLALSimBlackHoleRingdownModeEigenvalueEvaluateContinuedFractionLeaver(p->a, p->l, p->m, p->s, A, omega, XLALSimBlackHoleRingdownModeAngularCoefficientsLeaver), errnum);
	if (errnum)
		XLAL_ERROR(XLAL_EFUNC);
	XLAL_TRY(cf2 = XLALSimBlackHoleRingdownModeEigenvalueEvaluateContinuedFractionLeaver(p->a, p->l, p->m, p->s, A, omega, XLALSimBlackHoleRingdownModeRadialCoefficientsLeaver), errnum);
	if (errnum)
		XLAL_ERROR(XLAL_EFUNC);
	gsl_vector_set(f, 0, creal(cf1));
	gsl_vector_set(f, 1, cimag(cf1));
	gsl_vector_set(f, 2, creal(cf2));
	gsl_vector_set(f, 3, cimag(cf2));
	return 0;
}


/* Solves the continued fraction equation, */
/* Equation (13) of Leaver (1985), */
/* for the eigenfrequency omega of the */
/* quasinormal mode for Schwarzschild. */
static int XLALSimBlackHoleRingdownModeEigenvalueSolveSchwarzschild(COMPLEX16 *omega, int l, int m, int s)
{
	enum { ndim = 2 };
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *solver;
	int status;
	size_t iter = 0;
	int errnum;
	struct LALSimBlackHoleRingdownModeLeaver p;
	gsl_multiroot_function f = { &XLALSimBlackHoleRingdownModeSchwarzschildEigenvalueSolveResid, ndim, &p };
	gsl_vector *x = gsl_vector_alloc(ndim);

	if (!x)
		XLAL_ERROR(XLAL_ENOMEM);

	gsl_vector_set(x, 0, creal(*omega));
	gsl_vector_set(x, 1, cimag(*omega));
	p.a = 0.0;
	p.l = l;
	p.m = m;
	p.s = s;

	T = gsl_multiroot_fsolver_hybrids;
	solver = gsl_multiroot_fsolver_alloc(T, ndim);
	if (!solver) {
		gsl_vector_free(x);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	gsl_multiroot_fsolver_set(solver, &f, x);

	do {
		++iter;
		XLAL_TRY(status = gsl_multiroot_fsolver_iterate(solver), errnum);
		if (errnum) {
			gsl_multiroot_fsolver_free(solver);
			gsl_vector_free(x);
			XLAL_ERROR(XLAL_EFUNC);
		}
    		if (status)
      			break;
    		XLAL_TRY(status = gsl_multiroot_test_residual(solver->f, EPS), errnum);
		if (errnum) {
			gsl_multiroot_fsolver_free(solver);
			gsl_vector_free(x);
			XLAL_ERROR(XLAL_EFUNC);
		}
  	} while (status == GSL_CONTINUE && iter < MAXITER);
	if (iter >= MAXITER) {
		gsl_multiroot_fsolver_free(solver);
		gsl_vector_free(x);
		XLAL_ERROR(XLAL_EMAXITER);
	}

  	*omega = gsl_vector_get(solver->x, 0) + I*gsl_vector_get(solver->x, 1);

	gsl_multiroot_fsolver_free(solver);
	gsl_vector_free(x);
	return 0;
}


/* Solves the continued fraction equation, */
/* Equations (21) and (27) of Leaver (1985), */
/* for the eigenfrequency omega and angular */
/* separation constant A of the */
/* quasinormal mode for Kerr. */
static int XLALSimBlackHoleRingdownModeEigenvalueSolveKerr(COMPLEX16 *A, COMPLEX16 *omega, double a, int l, int m, int s)
{
	enum { ndim = 4 };
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *solver;
	int status;
	size_t iter = 0;
	struct LALSimBlackHoleRingdownModeLeaver p;
	int errnum;
	gsl_multiroot_function f = { &XLALSimBlackHoleRingdownModeKerrEigenvalueSolveResid, ndim, &p };
	gsl_vector *x = gsl_vector_alloc(ndim);

	if (!x)
		XLAL_ERROR(XLAL_ENOMEM);

	gsl_vector_set(x, 0, creal(*A));
	gsl_vector_set(x, 1, cimag(*A));
	gsl_vector_set(x, 2, creal(*omega));
	gsl_vector_set(x, 3, cimag(*omega));
	p.a = a;
	p.l = l;
	p.m = m;
	p.s = s;

	T = gsl_multiroot_fsolver_hybrids;
	solver = gsl_multiroot_fsolver_alloc(T, ndim);
	if (!solver) {
		gsl_vector_free(x);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	gsl_multiroot_fsolver_set(solver, &f, x);

	do {
		++iter;
		XLAL_TRY(status = gsl_multiroot_fsolver_iterate(solver), errnum);
		if (errnum) {
			gsl_multiroot_fsolver_free(solver);
			gsl_vector_free(x);
			XLAL_ERROR(XLAL_EFUNC);
		}
    		if (status)
      			break;
    		XLAL_TRY(status = gsl_multiroot_test_residual(solver->f, EPS), errnum);
		if (errnum) {
			gsl_multiroot_fsolver_free(solver);
			gsl_vector_free(x);
			XLAL_ERROR(XLAL_EFUNC);
		}
  	} while (status == GSL_CONTINUE && iter < MAXITER);
	if (iter >= MAXITER) {
		gsl_multiroot_fsolver_free(solver);
		gsl_vector_free(x);
		XLAL_ERROR(XLAL_EMAXITER);
	}

  	*A = gsl_vector_get(solver->x, 0) + I*gsl_vector_get(solver->x, 1);
  	*omega = gsl_vector_get(solver->x, 2) + I*gsl_vector_get(solver->x, 3);

  	/* problems occur if imaginary part of A is tiny but not zero */
  	if (fabs(cimag(*A)) < EPS)
    		*A = creal(*A);

	gsl_multiroot_fsolver_free(solver);
	gsl_vector_free(x);
	return 0;
}


/**
 * Low-level routine that computes the black hole quasinormal mode
 * eigenefrequency, omega, and angular separation constant A for a given
 * (l,m) mode and spin-weight s (s=-2 for gravitational perturbations).
 *
 * Implements Leaver's method by simultaneously
 * solving the continued fraction equations Eq. (21) and Eq. (27)
 * of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \warning The variables are represented in Leaver's conventions
 * in which G = c = 2M = 1.  In particular this means, |a| < 0.5.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdownModeEigenvaluesLeaver(
	COMPLEX16 *A,		/**< angular separation constant [returned] */
	COMPLEX16 *omega,		/**< eigenfrequency [returned] */
	double a,		/**< spin parameter (note: |a| < 0.5) */
	int l,			/**< mode value l */
	int m,			/**< mode value m */
	int s			/**< spin weight (s = -2 for gravitational perturbations) */
)
{
	double fac = 1.0/sqrt(27.0);
	double atry;
	int aneg = 0;

	/* negative values of a can be obtained from positive values by an identity */
	if (signbit(a)) {
		aneg = 1;
		a = fabs(a);
		m = -m;
	}

	if (a >= 0.5 || l < abs(s) || abs(m) > l || s > 0 || s < -2)
		XLAL_ERROR(XLAL_EINVAL);

	/* start at Schwarzschild values */
  	*A = l*(l+1) - s*(s+1);
	*omega = fac*(2*l+1-I); /* asymptotic value for large l */
	if (XLALSimBlackHoleRingdownModeEigenvalueSolveSchwarzschild(omega, l, m, s) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* step towards requested value of a */
  	for (atry = 0; atry < a; atry += 0.1 * (0.5 - atry))
		if (XLALSimBlackHoleRingdownModeEigenvalueSolveKerr(A, omega, atry, l, m, s) < 0)
			XLAL_ERROR(XLAL_EFUNC);

  	/* now use the current guess to get value at requested a */
  	if (XLALSimBlackHoleRingdownModeEigenvalueSolveKerr(A, omega, a, l, m, s) < 0)
		XLAL_ERROR(XLAL_EFUNC);

  	/* if a was negative, apply the identity */
  	if (aneg) {
    		*A = conj(*A);
    		*omega = I*conj(-I*(*omega));
  	}

	return 0;
}


/* Equations 18 and 19 of Leaver (1985) */
static COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver(double mu, double a, int l, int m, int s, COMPLEX16 A, COMPLEX16 omega)
{
	COMPLEX16 alp, bet, gam;
	double mup1 = mu + 1;
	double mum1 = mu - 1;
	COMPLEX16 prod = 1;
	COMPLEX16 sum;
	COMPLEX16 delta;
	COMPLEX16 a_n;
	COMPLEX16 a_nm1;
	COMPLEX16 a_np1;
	int n = 0;

	XLALSimBlackHoleRingdownModeAngularCoefficientsLeaver(&alp, &bet, &gam, a, l, m, n, s, A, omega);
	sum = 1.0;
	a_n = 1.0;
	a_np1 = -bet*a_n/alp; /* Eq. 19, first line */
	while (n++ < MAXITER) {
		a_nm1 = a_n;
		a_n   = a_np1;
		prod *= mup1;
		delta = a_n*prod;
		sum  += delta;
		if (cabs(delta)*cabs(delta) < /* FIXME EPS */ 1e-4)
			break;
		XLALSimBlackHoleRingdownModeAngularCoefficientsLeaver(&alp, &bet, &gam, a, l, m, n, s, A, omega);
		a_np1 = -(bet*a_n + gam*a_nm1)/alp; /* Eq. 19, second line */
	}
  	if (n >= MAXITER)
		XLAL_ERROR(XLAL_EMAXITER);
  	sum *= pow(mup1, 0.5*abs(m-s));
  	sum *= pow(-mum1, 0.5*abs(m+s));
  	sum *= cexp(a*omega*mu);
  	return sum;
}


/* Integrand function for gsl integrator used to compute
 * normalization factor for spheroidal wave functions. */
static double XLALSimBlackHoleRingdownSpheriodalWaveFunctionNormIntegrand(double mu, void *params)
{
	struct LALSimBlackHoleRingdownModeLeaver *p = params;
	COMPLEX16 sphwf;
	double r;
 	int errnum; 

	XLAL_TRY(sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver(mu, p->a, p->l, p->m, p->s, p->A, p->omega), errnum);
	if (errnum)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	r = cabs(sphwf);


	return r*r;
}


/* Computes the normalization factor for spheroidal wave functions. */
static COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunctionNormLeaver(double a, int l, int m, int s, COMPLEX16 A, COMPLEX16 omega)
{
	struct LALSimBlackHoleRingdownModeLeaver p;
	enum { WORKSPACESIZE = 1000 };
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(WORKSPACESIZE);
	double integral, error;
	COMPLEX16 sphwf;
	COMPLEX16 norm;
	int signneg;
	gsl_function f;
 	int errnum;
	int status;

	p.a = a;
	p.l = l;
	p.m = m;
	p.s = s;
	p.A = A;
	p.omega = omega;
	f.function = XLALSimBlackHoleRingdownSpheriodalWaveFunctionNormIntegrand;
	f.params = &p;
	XLAL_TRY(status = gsl_integration_qags(&f, -1.0, 1.0, 0.0, 1e-6, WORKSPACESIZE, w, &integral, &error), errnum);
	gsl_integration_workspace_free(w);
	if (status || errnum) /* there was an error during integration */
		XLAL_ERROR_REAL8(XLAL_EFUNC);

	/* get complex part so that sphwf is real at mu=0 */
	XLAL_TRY(sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver(0.0, a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	norm = cabs(sphwf)/sphwf;

	XLAL_TRY(sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver(-1.0 + EPS, a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	sphwf *= norm;
	signneg = signbit(creal(sphwf));
	signneg = (l - (m>s?m:s))%2 ? !signneg : signneg;

	norm /= sqrt(integral);
	if (signneg)
		norm *= -1;

	return norm;
}


/**
 * Low-level routine that evaluates the spheroidal wave function at a
 * specified value of mu = cos(theta) for a given (l,m) mode and
 * spin-weight s (s=-2 for gravitational perturbations).
 * Also requires the angular separation constant A and eigenvalue
 * omega for that mode, which are calculated by the routine
 * XLALSimBlackHoleRingdownModeEigenvaluesLeaver().
 *
 * Implements Leaver's method by simultaneously
 * solving the continued fraction equations Eq. (21) and Eq. (27)
 * of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \warning The variables are represented in Leaver's conventions
 * in which G = c = 2M = 1.  In particular this means, |a| < 0.5.
 *
 * \todo Extend so that overtones can be computed too.
 */
COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(
	double mu,	/**< cosine of polar angle */
	double a,	/**< spin parameter (note: |a| < 0.5) */
	int l,		/**< mode value l */
	int m,		/**< mode value m */
	int s,		/**< spin weight (s = -2 for gravitational perturbations) */
	COMPLEX16 A,	/**< angular separation constant */
	COMPLEX16 omega	/**< eigenfrequency */
)
{
	COMPLEX16 norm;
	COMPLEX16 sphwf;
	int errnum;

	if (fabs(mu) > 1.0 || fabs(a) >= 0.5 || l < abs(s) || abs(m) > l || s > 0 || s < -2)
		XLAL_ERROR(XLAL_EINVAL);

	XLAL_TRY(norm = XLALSimBlackHoleRingdownSpheroidalWaveFunctionNormLeaver(a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	XLAL_TRY(sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver(mu, a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	sphwf *= norm;

	return sphwf;
}


/**
 * Computes the frequency and quality factor of a specified quasinormal
 * mode (l,m) of spin weight s perturbations (s=-2 for gravitational
 * perturbations) of a black hole of a specified mass and spin. 
 *
 * Uses the method of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdownMode(
	double *frequency,		/**< mode frequency (Hz) [returned] */
	double *quality,		/**< mode quality factor [returned] */
	double mass,			/**< black hole mass (kg) */
	double dimensionless_spin,	/**< black hole dimensionless spin parameter (-1,+1) */
	int l,				/**< polar mode number */
	int m,				/**< azimuthal mode number */
	int s				/**< spin weight (s=-2 for gravitational radiation) */
)
{
	double a = 0.5*dimensionless_spin; /* convert to Leaver's convention 2M = 1 */
	COMPLEX16 A, omega;
	if (XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s) < 0)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	omega *= 0.5; /* convert from Leaver's convention 2M = 1 */
	*frequency = fabs(creal(omega)) / (LAL_TWOPI * mass);
	*quality   = fabs(creal(omega)) / (-2.0 * cimag(omega));
	return 0;
}


/**
 * Evaluates the value of spheroidal wave function at a given
 * polar angle theta for a specified mode (l,m) and spin weight s
 * (s=-2 for gravitational perturbations) and
 * dimensionless spin parameter.
 *
 * Uses the method of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunction(
	double theta,			/**< polar angle (radians) */
	double dimensionless_spin,	/**< black hole dimensionless spin parameter */
	int l,				/**< polar mode number */
	int m,				/**< azimuthal mode number */
	int s				/**< spin weight (s=-2 for gravitational radiation) */
)
{
	double a = 0.5*dimensionless_spin; /* convert to Leaver conventions 2M = 1 */
	double mu = cos(theta);
	COMPLEX16 A, omega;
	COMPLEX16 sphwf;
	
	if (XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s) < 0)
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(mu, a, l, m, s, A, omega);
	return sphwf;
}


/**
 * Computes the waveform for the ringdown of a black hole
 * quasinormal mode (l,m).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdown(
	REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
	REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
	const LIGOTimeGPS *t0,		/**< start time of ringdown */
	double phi0,			/**< initial phase of ringdown (rad) */
	double deltaT,			/**< sampling interval (s) */
	double mass,			/**< black hole mass (kg) */
	double dimensionless_spin,	/**< black hole dimensionless spin parameter */
	double fractional_mass_loss,	/**< fraction of mass radiated in this mode */
	double distance,		/**< distance to source (m) */
	double inclination,		/**< inclination of source's spin axis (rad) */
	int l,				/**< polar mode number */
	int m				/**< azimuthal mode number */
)
{
	const int s = -2; /* spin weight for gravitational radiation */
	double mu = cos(inclination);
	double a = 0.5*dimensionless_spin; /* convert to Leaver conventions 2M = 1 */
	COMPLEX16 A, omega;
	COMPLEX16 sphwf1, sphwf2;
	COMPLEX16 A1, A2;
	COMPLEX16 omega_dt;
	size_t length;
	size_t j;
	int errnum;

	if (XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	XLAL_TRY(sphwf1 = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(mu, a, l, m, s, A, omega), errnum);
	XLAL_TRY(sphwf2 = XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(-mu, a, l, m, s, A, omega), errnum);
	if (errnum)
		XLAL_ERROR(XLAL_EFUNC);
	omega *= 0.5; /* convert from Leaver's convention 2M = 1 */
	
	/* compute length of waveform to compute */
	length = ceil(log(LAL_REAL8_EPS) * LAL_G_SI * mass / (pow(LAL_C_SI, 3.0) * cimag(omega) * deltaT));
	if (length < 1)
		XLAL_ERROR(XLAL_EBADLEN);

	/* compute the amplitude factors for the +m and -m modes */
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

	/* compute the waveforms */
	for (j = 0; j < length; ++j) {
		COMPLEX16 h;
		h = A1*cexp(-I*omega_dt*j) + A2*cexp(I*conj(omega_dt)*j);
		(*hplus)->data->data[j] = creal(h);
		(*hcross)->data->data[j] = -cimag(h);
	}

	return 0;
}


/* TEST CODE */

#if 0 /* TEST CODE */

#include <stdio.h>
int ringdown_waveform(void)
{
	static LIGOTimeGPS epoch;
	const double M = 10.0*LAL_MSUN_SI;
	const double a = 0.98;
	const double r = 1e6*LAL_PC_SI;
	const double i = LAL_PI/6;
	const double dt = 1.0/16384.0;
	const double eps = 0.01;
	const int l = 2;
	const int m = 2;
	const char *fname = "ringdown_waveform.dat";
	REAL8TimeSeries *hplus = NULL;
	REAL8TimeSeries *hcross = NULL;
	size_t j;
	FILE *fp;

	fp = fopen(fname, "w");
	XLALSimBlackHoleRingdown(&hplus, &hcross, &epoch, 0.0, dt, M, a, eps, r, i, l, m);
	for (j = 0; j < hplus->data->length; ++j)
		fprintf(fp, "%f\t%e\t%e\n", j*dt, hplus->data->data[j], hcross->data->data[j]);
	fclose(fp);

	XLALDestroyREAL8TimeSeries(hcross);
	XLALDestroyREAL8TimeSeries(hplus);

	return 0;
}

int grasp_spherical_table(void)
{
	const double a = 0.0;
	const int l = 2;
	const int m = 2;
	const int s = -2;
	const char *fname = "grasp_spherical_table.txt";
	FILE *fp;
	COMPLEX16 A, omega;
	COMPLEX16 norm;
	COMPLEX16 sphwf;
	double fac = 1.0/sqrt(2.0*M_PI);
	double exactfac = sqrt(5.0/(64.0*M_PI));
	double muvec[] = {-0.99,-0.95,-0.75,-0.55,-0.35,-0.15,0.15,0.35,0.55,0.75,0.95,0.99};
	size_t i;

	fp = fopen(fname, "w");
	/* note: multiply a by 0.5 for Leaver conventions */
	XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, 0.5*a, l, m, s);
	norm = XLALSimBlackHoleRingdownSpheroidalWaveFunctionNormLeaver(0.5*a, l, m, s, A, omega);
	for (i = 0; i < sizeof(muvec)/sizeof(*muvec); ++i) {
		double mu = muvec[i];
		sphwf = norm * XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver( mu, a, l, m, s, A, omega );
		fprintf(fp, "%+e\t%e\t%e\n", mu, fac*creal(sphwf), exactfac*pow(1.0+mu,2));
	}
	fclose(fp);

	return 0;
}

int grasp_spheroid_figure(void) 
{
	const double a = 0.98;
	const int l = 2;
	const int m = 2;
	const int s = -2;
	const char *fname = "grasp_spheroid_figure.dat";
	FILE *fp;
	COMPLEX16 A, omega;
	COMPLEX16 norm;
	COMPLEX16 sphwf;
	double mu;

	fp = fopen(fname, "w");
	/* note: multiply a by 0.5 for Leaver conventions */
	XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, 0.5*a, l, m, s);
	norm = XLALSimBlackHoleRingdownSpheroidalWaveFunctionNormLeaver(0.5*a, l, m, s, A, omega);
	for (mu = -1.0; mu <= 1.0; mu += 0.01) {
		sphwf = norm * XLALSimBlackHoleRingdownSpheroidalWaveFunction1Leaver(mu, 0.5*a, l, m, s, A, omega);
		fprintf(fp, "%e\t%e\t%e\n", mu, creal(sphwf), cimag(sphwf));
	}
	fclose(fp);

	return 0;
}

int leaver_table_2(void)
{
	const int l = 2;
	const int m = 0;
	const int s = -2;
	const char *fname = "leaver_table_2.txt";
	FILE *fp;
	COMPLEX16 A, omega;
	double avec[] = {0.0,0.1,0.2,0.3,0.4,0.45,0.49,0.4999};
	size_t i;

	fp = fopen(fname, "w");

	/* positive values of a */
	for ( i = 0; i < sizeof(avec)/sizeof(*avec); ++i ) {
		double a = avec[i];
		XLALSimBlackHoleRingdownModeEigenvaluesLeaver( &A, &omega, a, l, m, s);
		fprintf(fp, "%+.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, creal(A), cimag(A), creal(omega), cimag(omega));
	}

	fprintf(fp, "\n");

	/* negative values of a */
	for (i = 0; i < sizeof(avec)/sizeof(*avec); ++i) {
		double a = -avec[i];
		XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s);
		fprintf(fp, "%+.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, creal(A), cimag(A), creal(omega), cimag(omega));
	}

	fclose(fp);
	return 0;
}

int leaver_table_3(void)
{
	const int l = 2;
	const int m = 1;
	const int s = -2;
	const char *fname = "leaver_table_3.txt";
	FILE *fp;
	COMPLEX16 A, omega;
	double avec[] = {0.0,0.1,0.2,0.3,0.4,0.45,0.49,0.4999};
	size_t i;

	fp = fopen(fname, "w");

	/* positive values of a */
	for (i = 0; i < sizeof(avec)/sizeof(*avec); ++i) {
		double a = avec[i];
		XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s);
		fprintf(fp, "%+.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, creal(A), cimag(A), creal(omega), cimag(omega));
	}

	fprintf(fp, "\n");

	/* negative values of a */
	for (i = 0; i < sizeof(avec)/sizeof(*avec); ++i) {
		double a = -avec[i];
		XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s);
		fprintf(fp, "%+.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, creal(A), cimag(A), creal(omega), cimag(omega));
	}

	fclose(fp);
	return 0;
}

int main(void)
{
	lalDebugLevel = 7;
	XLALSetErrorHandler(XLALAbortErrorHandler);
	leaver_table_2();
	leaver_table_3();
	grasp_spherical_table();
	grasp_spheroid_figure();
	ringdown_waveform();
	return 0;
}

#endif /* TEST CODE */


/**
 * Computes the final mass and spin of the black hole resulting from merger. 
 * They are given by fittings of NR simulations results. Specifically,
 * for EOBNR, Table I of Buonanno et al. PRD76, 104049;
 * for EOBNRv2 and EOBNRv2HM, Eqs. 29a and 29b of Pan et al. PRD84, 124052;
 * for SEOBNRv1, Eq. 8 of Tichy and Marronetti PRD78, 081501 and 
 *               Eqs. 1 and 3 of Barausse and Rezzolla ApJ704, L40.
 */
INT4 XLALSimIMREOBFinalMassSpin(
	REAL8    *finalMass,  /**<< OUTPUT, the final mass (scaled by original total mass) */
	REAL8    *finalSpin,  /**<< OUTPUT, the final spin (scaled by final mass) */
  const REAL8     mass1,      /**<< The mass of the 1st component of the system */
  const REAL8     mass2,      /**<< The mass of the 2nd component of the system */
  const REAL8     spin1[3],   /**<< The spin of the 1st object; only needed for spin waveforms */
  const REAL8     spin2[3],   /**<< The spin of the 2nd object; only needed for spin waveforms */
  Approximant     approximant /**<< The waveform approximant being used */
)
{
  static const REAL8 root9ovr8minus1 = -0.057190958417936644;
  static const REAL8 root12          = 3.4641016151377544;

  /* Constants used for the spinning EOB model */
  static const REAL8 p0 = 0.04826;
  static const REAL8 p1 = 0.01559;
  static const REAL8 p2 = 0.00485;
  static const REAL8 t0 = -2.8904;
  static const REAL8 t2 = -3.5171;
  static const REAL8 t3 = 2.5763;
  static const REAL8 s4 = -0.1229;
  static const REAL8 s5 = 0.4537;

  REAL8 totalMass;
  REAL8 eta, eta2, eta3;
  REAL8 a1, a2, chiS, q;
  REAL8 tmpVar;

  /* get a local copy of the intrinsic parameters */
  totalMass = mass1 + mass2;
  eta = mass1 * mass2 / (totalMass*totalMass);
  eta2 = eta * eta;
  eta3 = eta2 * eta;

  switch ( approximant )
  {
    case EOBNRv2:
    case EOBNRv2HM:
      eta3 = eta2 * eta;
      /* Final mass and spin given by Eqs. 29a and 29b of Pan et al. PRD84, 124052 */
      *finalMass = 1. + root9ovr8minus1 * eta - 0.4333 * eta2 - 0.4392 * eta3;
      *finalSpin = root12 * eta - 3.871 * eta2 + 4.028 * eta3;
      break;
    case EOBNR:
      /* Final mass and spin given by Table I of Buonanno et al. PRD76, 104049 */
      *finalMass = 1 - 0.057191 * eta - 0.498 * eta2;
      *finalSpin = 3.464102 * eta - 2.9 * eta2;
      break;
    case SEOBNRv1:
      /* Final mass/spin comes from Eq. 8 of Tichy and Marronetti PRD78, 081501 
         and from Eqs. 1 and 3 of Barausse and Rezzolla ApJ704, L40 */
      a1 = spin1[2];
      a2 = spin2[2];

      chiS = 0.5 * ( a1 + a2 );
      q    = mass1 / mass2;

      *finalMass = 1. - p0 + 2.*p1 * chiS  + 4.*p2*chiS*chiS;
      *finalMass = 1. + (0.9515 - 1.0)*4.*eta - 0.013*16.*eta2*(a1+a2);
      tmpVar     = ( a1 + a2 /q/q) / ( 1. + 1/q/q);
      *finalSpin = tmpVar + tmpVar * eta * (s4 * tmpVar + s5 * eta + t0 ) + eta * (2. * sqrt(3.) + t2*eta + t3*eta*eta );
      break;
    default:
      XLALPrintError( "XLAL Error %s - Unsupported approximant.\n", __func__ );
      XLAL_ERROR( XLAL_EINVAL );
  }

  /*printf( "Final mass = %e, Final spin = %e\n", *finalMass, *finalSpin );*/
  return XLAL_SUCCESS;
}


/**
 * This function generates the quasinormal mode frequencies for a black
 * hole ringdown. At present, this function works for the 22, 21, 33, 44
 * and 55 modes, and includes 8 overtones. The final frequencies are
 * computed by interpolating the data found on the webpage of 
 * Vitor Cardoso, http://centra.ist.utl.pt/~vitor/?page=ringdown
 * In this page, frequecy data are given for positive final spins only.
 * For a negative final spin chi<0 case, the (l,m) mode frequency is given by 
 * the (l,-m) mode frequency of the positive final spin -chi case.
 */
INT4 XLALSimIMREOBGenerateQNMFreqV2(
  COMPLEX16Vector *modefreqs, /**<< OUTPUT, complex freqs of overtones in unit of Hz */
  const REAL8      mass1,     /**<< The mass of the 1st component (in Solar masses) */
  const REAL8      mass2,     /**<< The mass of the 2nd component (in Solar masses) */
  const REAL8      spin1[3],  /**<< The spin of the 1st object; only needed for spin waveforms */
  const REAL8      spin2[3],  /**<< The spin of the 2nd object; only needed for spin waveforms */
  UINT4            l,         /**<< The l value of the mode in question */
  UINT4            m,         /**<< The m value of the mode in question */
  UINT4            nmodes,    /**<< The number of overtones that should be included (max 8) */
  Approximant      approximant/**<< The waveform approximant being used */
  )
{

  /* Data for interpolating the quasinormal mode frequencies is taken from 
   * The webpage of Vitor Cardoso, http://centra.ist.utl.pt/~vitor/?page=ringdown 
   * In the spin range of (-0.995, 0.999), 
   * interpolation error within 5e-5 for the first 3 overtones of the 22 mode 
   * interpolation error within 0.005 for the first 8 overtones of the 22 mode 
   * On the data points, interpolation error is within 1e-5, except for the 
   * highest overtone of the (4,4) and (5,5) modes.
   */

  static const double afinallist[107] = {-0.9996, -0.9995, -0.9994, -0.9992, -0.999, -0.9989, -0.9988, 
  -0.9987, -0.9986, -0.9985, -0.998, -0.9975, -0.997, -0.996, -0.995, -0.994, -0.992, -0.99, -0.988, 
  -0.986, -0.984, -0.982, -0.98, -0.975, -0.97, -0.96, -0.95, -0.94, -0.92, -0.9, -0.88, -0.86, -0.84, 
  -0.82, -0.8, -0.78, -0.76, -0.74, -0.72, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, 
  -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 
  0.65, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 
  0.975, 0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.995, 0.996, 0.997, 0.9975, 0.998, 
  0.9985, 0.9986, 0.9987, 0.9988, 0.9989, 0.999, 0.9992, 0.9994, 0.9995, 0.9996};

  /* 2, 2 mode */

  const double reomegaqnm22[8][107] = {{0.270228, 0.276562, 0.280636, 0.285234, 0.287548, 0.288282, 0.288845, 0.289287, 0.289639, 0.289924, 0.290781, 0.291189, 0.291418, 0.291658, 0.291785, 0.29187, 0.291998, 0.292111, 0.292221, 0.292331, 0.292441, 0.292552, 0.292664, 0.292943, 0.293223, 0.293787, 0.294354, 0.294925, 0.296077, 0.297244, 0.298426, 0.299624, 0.300837, 0.302067, 0.303313, 0.304577, 0.305857, 0.307156, 0.308473, 0.309808, 0.313232, 0.316784, 0.320473, 0.324307, 0.328299, 0.332458, 0.336798, 0.341333, 0.346079, 0.351053, 0.356275, 0.361768, 0.367557, 0.373672, 0.380146, 0.387018, 0.394333, 0.402145, 0.410518, 0.419527, 0.429264, 0.439842, 0.451402, 0.464123, 0.478235, 0.494045, 0.511969, 0.5326, 0.541794, 0.55163, 0.562201, 0.573616, 0.586017, 0.59958, 0.614539, 0.631206, 0.650018, 0.671614, 0.696995, 0.727875, 0.74632, 0.767674, 0.793208, 0.808235, 0.825429, 0.8331, 0.841343, 0.850272, 0.860046, 0.870893, 0.883162, 0.897446, 0.905664, 0.914902, 0.925581, 0.931689, 0.938524, 0.946385, 0.948123, 0.949929, 0.951813, 0.953784, 0.955854, 0.960358, 0.965514, 0.968438, 0.97169}, 
{0.243689, 0.247076, 0.248937, 0.250577, 0.251041, 0.251088, 0.251073, 0.251023, 0.250952, 0.250874, 0.25052, 0.250339, 0.250292, 0.250359, 0.250463, 0.250549, 0.250686, 0.250814, 0.250945, 0.251079, 0.251214, 0.251349, 0.251484, 0.25182, 0.252158, 0.252836, 0.253519, 0.254205, 0.255591, 0.256994, 0.258413, 0.259851, 0.261306, 0.26278, 0.264273, 0.265785, 0.267317, 0.268869, 0.270442, 0.272036, 0.276117, 0.280342, 0.284721, 0.289262, 0.293978, 0.29888, 0.303981, 0.309295, 0.314838, 0.320629, 0.326687, 0.333036, 0.339701, 0.346711, 0.354101, 0.36191, 0.370183, 0.378976, 0.388353, 0.39839, 0.409183, 0.420847, 0.433527, 0.447407, 0.462728, 0.479807, 0.499079, 0.521161, 0.53097, 0.541447, 0.552684, 0.564795, 0.577922, 0.592247, 0.608001, 0.625499, 0.645173, 0.667658, 0.693938, 0.725708, 0.744582, 0.766349, 0.792272, 0.807482, 0.824852, 0.832591, 0.840901, 0.849896, 0.859735, 0.870645, 0.882977, 0.897322, 0.905568, 0.914834, 0.925538, 0.931657, 0.938502, 0.946371, 0.94811, 0.949919, 0.951804, 0.953776, 0.955847, 0.960353, 0.96551, 0.968435, 0.971688}, 
{0.191626, 0.188311, 0.185222, 0.182107, 0.182532, 0.183173, 0.183784, 0.184279, 0.184634, 0.18486, 0.184811, 0.184353, 0.184263, 0.184581, 0.184713, 0.18472, 0.184862, 0.185055, 0.185205, 0.185357, 0.185521, 0.185687, 0.18585, 0.18625, 0.186658, 0.187473, 0.188294, 0.189121, 0.190789, 0.192478, 0.194189, 0.195923, 0.19768, 0.19946, 0.201263, 0.203092, 0.204945, 0.206823, 0.208728, 0.210659, 0.215608, 0.220737, 0.226056, 0.231577, 0.237312, 0.243272, 0.249472, 0.255926, 0.262651, 0.269665, 0.276988, 0.284642, 0.292654, 0.301053, 0.309873, 0.319153, 0.328939, 0.339285, 0.350255, 0.361927, 0.374396, 0.387779, 0.402225, 0.417925, 0.43513, 0.454179, 0.475545, 0.499906, 0.510701, 0.522217, 0.534557, 0.547847, 0.56224, 0.577926, 0.595148, 0.614222, 0.635576, 0.659827, 0.687923, 0.721496, 0.741239, 0.763831, 0.790516, 0.806079, 0.823779, 0.831643, 0.840075, 0.849189, 0.859144, 0.870168, 0.882612, 0.897067, 0.905368, 0.914687, 0.925444, 0.931587, 0.938454, 0.946343, 0.948085, 0.949897, 0.951785, 0.95376, 0.955833, 0.960343, 0.965503, 0.96843, 0.971683}, 
{0.127766, 0.134925, 0.137314, 0.135026, 0.132545, 0.132613, 0.133213, 0.133927, 0.134502, 0.134846, 0.134197, 
0.133839, 0.134383, 0.134311, 0.134305, 0.134548, 0.134538, 0.134764, 0.134855, 0.13499, 0.135153, 0.135276, 0.135401, 0.13576, 0.136094, 0.136791, 0.137494, 0.138199, 0.139634, 0.141097, 0.142587, 0.144106, 0.145656, 0.147236, 0.148849, 0.150493, 0.152171, 0.153884, 0.155632, 0.157416, 0.162043, 0.16692, 0.172065, 0.177495, 0.183228, 0.18928, 0.19567, 0.202416, 0.209537, 0.217053, 0.224984, 0.233352, 0.242183, 0.251505, 0.261348, 0.271749, 0.282749, 0.294397, 0.306753, 0.319886, 0.333884, 0.348856, 0.364938, 0.382309, 0.40121, 0.421977, 0.445099, 0.471336, 0.482957, 0.495375, 0.508729, 0.523187, 0.538956, 0.556285, 0.575464, 0.59683, 0.620788, 0.647869, 0.678894, 0.715333, 0.736426, 0.760279, 0.788101, 0.804173, 0.822338, 0.830374, 0.83897, 0.848242, 0.858347, 0.869515, 0.8821, 0.896696, 0.905068, 0.91446, 0.925289, 0.931468, 0.938371, 0.946292, 0.948041, 0.949858, 0.951752, 0.953732, 0.955811, 0.960329, 0.965496, 0.968424, 0.971679}, 
{0.111845, 0.106167, 0.107777, 0.111206, 0.108936, 0.108364, 0.108639, 0.109304, 0.109887, 0.110167, 0.109102, 0.109714, 0.109696, 0.109638, 0.109676, 0.109676, 0.109779, 0.109961, 0.109986, 0.110136, 0.11021, 0.110304, 0.110427, 0.110658, 0.110906, 0.111408, 0.111916, 0.112428, 0.113458, 0.114508, 0.115576, 0.116663, 0.117771, 0.118902, 0.120054, 0.12123, 0.122431, 0.123659, 0.124914, 0.126199, 0.12955, 0.133123, 0.136951, 0.141068, 0.145514, 0.150331, 0.155561, 0.161252, 0.167448, 0.174197, 0.18154, 0.189518, 0.198166, 0.207515, 0.217594, 0.228433, 0.24006, 0.252511, 0.265825, 0.280057, 0.295269, 0.311546, 0.328992, 0.347741, 0.367967, 0.389902, 0.413872, 0.440385, 0.451887, 0.464034, 0.476972, 0.490925, 0.506263, 0.523592, 0.543861, 0.568205, 0.597077, 0.629836, 0.666171, 0.707148, 0.730203, 0.755806, 0.785154, 0.801892, 0.82065, 0.828901, 0.837698, 0.847157, 0.857436, 0.868766, 0.881504, 0.896248, 0.904696, 0.914167, 0.925079, 0.931302, 0.938248, 0.946214, 0.947971, 0.949797, 0.951699, 0.953687, 0.955773, 0.960306, 0.965484, 0.968417, 0.971675}, 
{0.0973518, 0.0987518, 0.0952128, 0.097231, 0.0968051, 0.0959061, 0.0958383, 0.0963991, 0.0969597, 0.0971588, 0.0962907, 0.0968828, 0.0964446, 0.0967173, 0.0968133, 0.0967445, 0.0968961, 0.0970152, 0.0970355, 0.0971754, 0.0972259, 0.0973507, 0.097408, 0.0976397, 0.0978753, 0.0983115, 0.0987594, 0.099209, 0.100108, 0.101005, 0.101904, 0.102805, 0.103707, 0.104611, 0.10552, 0.106431, 0.107346, 0.108266, 0.109192, 0.110125, 0.112493, 0.114931, 0.117463, 0.120121, 0.122943, 0.125975, 0.129278, 0.132931, 0.137032, 0.14171, 0.147116, 0.153415, 0.160767, 0.169299, 0.179085, 0.19014, 0.202441, 0.215947, 0.230623, 0.246457, 0.263464, 0.281691, 0.301218, 0.322157, 0.344655, 0.368896, 0.39509, 0.42345, 0.43543, 0.447764, 0.460418, 0.473311, 0.486283, 0.499037,
0.511101, 0.522095, 0.532461, 0.536917, 0.527548, 0.521465, 0.518438, 0.515417, 0.51241, 0.510911, 0.509415, 0.508818, 0.508221, 0.507624, 0.507028, 0.506432, 0.505837, 0.505242, 0.504944, 0.504647, 0.50435, 0.504202, 0.504054, 0.503904, 0.503874, 0.503845, 0.503817, 0.503789, 0.50376, 0.503692, 0.503638, 0.503632, 0.503585},
{0.084098, 0.0868741, 0.0879148, 0.0858816, 0.0871769, 0.0862173, 0.0858779, 0.0863779, 0.0869438, 0.0870509, 0.0866021, 0.0864788, 0.0867758, 0.0866752, 0.0866613, 0.0868868, 0.086965, 0.0870574, 0.0871075, 0.0872251, 0.0873545, 0.0874293, 0.087563, 0.0878271, 0.088086, 0.0886037, 0.0891314, 0.0896576, 0.0906978, 0.0917299, 0.0927446, 0.093744, 0.0947274, 0.0956935, 0.0966414, 0.0975719, 0.0984848, 0.0993786, 0.100255, 0.101113, 0.103185, 0.105155, 0.107036, 0.108838, 0.110571, 0.112244, 0.113866, 0.115452, 0.117034, 0.118689, 0.120596, 0.123138, 0.127029, 0.133252, 0.142539, 0.154841, 0.169536, 0.185971, 0.203721, 0.222585, 0.242517, 0.263571, 0.285865, 0.309569, 0.334899, 0.362133, 0.391629, 0.423854, 0.437641, 0.45201, 0.467014, 0.482713, 0.499165, 0.516428, 0.534543, 0.553525, 0.573596, 0.602227, 0.649295, 0.697231, 0.722841, 0.750611, 0.781797, 0.799333,
0.818797, 0.827303, 0.836336, 0.846012, 0.856487, 0.867993, 0.880887, 0.895775, 0.904293, 0.913837, 0.92483, 0.931097, 0.938091, 0.946107, 0.947875, 0.949711, 0.951623, 0.953621, 0.955717, 0.960269, 0.965464, 0.968404, 0.971669}, 
{0.0783357, 0.0733796, 0.0766204, 0.0744319, 0.0764107, 0.0754793, 0.0749396, 0.0755031, 0.0761004, 0.0760294, 0.0759874, 0.075511, 0.0757765, 0.0758322, 0.0759714, 0.0759232, 0.0761145, 0.0762289, 0.0763205, 0.0764184, 0.0765972, 0.0767207, 0.0768409, 0.0771736, 0.0775187, 0.07821, 0.0789117, 0.0796137, 0.0810158, 0.0824117, 0.0837964, 0.0851528, 0.0864792, 0.0877698, 0.0890208, 0.0902295, 0.0913937, 0.0925132, 0.0935883, 0.0946171, 0.0969926, 0.0990908, 0.100921, 0.102489, 0.103795, 0.104821, 0.10553, 0.105854, 0.105679, 0.104817, 0.102968, 0.0997112, 0.094996, 0.0928224, 0.102996, 0.122468, 0.144233,  0.166117, 0.187881, 0.209753, 0.232042, 0.255048, 0.279057, 0.304354, 0.331244, 0.360084, 0.391319, 0.425537, 0.440234, 0.455604, 0.471729, 0.488705, 0.506652, 0.52572, 0.546107, 0.568081, 0.592057, 0.618799, 0.650217, 0.691189, 0.716737, 0.745623, 0.778341, 0.796657, 0.816858, 0.825638, 0.834931, 0.844846, 0.855539, 0.867237, 0.880294, 0.895318, 0.903898, 0.913504, 0.924565, 0.930871, 0.93791, 0.945976, 0.947755, 0.949602, 0.951525, 0.953535, 0.955643, 0.960217, 0.965434, 0.968355, 0.971639}};

  const double imomegaqnm22[8][107] = {{0.0784607, 0.0801007, 0.0816472, 0.0840037, 0.0854954, 0.0860129, 0.0864215, 0.0867455, 0.0870035, 0.08721, 0.0877802, 0.0879878, 0.088065, 0.0880965, 0.0880878, 0.0880757, 0.0880612, 0.0880574, 0.0880589, 0.0880627, 0.0880675, 0.0880727, 0.088078, 0.0880913, 0.0881045, 0.0881304, 0.088156, 0.0881813, 0.0882315, 0.0882807, 0.0883289, 0.0883763, 0.0884226, 0.0884679, 0.0885122, 0.0885555, 0.0885976, 0.0886386, 0.0886785, 0.0887172, 0.0888085, 0.0888917, 0.0889663, 0.0890315, 0.0890868, 0.0891313, 0.0891643, 0.0891846, 0.0891911, 0.0891825, 0.0891574, 0.0891138, 0.0890496, 0.0889623, 0.0888489, 0.0887057, 0.0885283, 0.0883112, 0.0880477, 0.0877293, 0.0873453, 0.086882, 0.0863212, 0.0856388, 0.0848021, 0.0837652, 0.0824618, 0.0807929, 0.0799908, 0.0790927, 0.0780817, 0.0769364, 0.0756296, 0.0741258, 0.072378, 0.0703215, 0.0678642, 0.0648692, 0.0611186, 0.0562313, 0.053149, 0.0494336, 0.0447904, 0.0419586, 0.0386302, 0.0371155, 0.0354676, 0.033659, 0.0316516, 0.0293904, 0.0267908, 0.0237095, 0.0219107, 0.0198661, 0.0174737, 0.0160919, 0.014534, 0.0127274, 0.0123259, 0.0119077, 0.0114708, 0.0110127, 0.0105306, 0.00947799, 0.00826689, 0.00757704, 0.0068074}, 
{0.284478, 0.282567, 0.281162, 0.279245, 0.278072, 0.27767, 0.277358, 0.277118, 0.276936, 0.2768, 0.276536, 0.276565, 0.276647, 0.276753, 0.276778, 0.276775, 0.276764, 0.276768, 0.276776, 0.276784, 0.276791, 0.276797, 0.276803, 0.276818, 0.276834, 0.276864, 0.276893, 0.276922, 0.276977, 0.277028, 0.277076, 0.27712, 0.27716, 0.277197, 0.27723, 0.277259, 0.277283, 0.277304, 0.277321, 0.277333, 0.277344, 0.277326, 0.277278, 0.277196, 0.27708, 0.276926, 0.276732, 0.276495, 0.276212, 0.275877, 0.275486, 0.275033, 0.274512, 0.273915, 0.273232, 0.272452, 0.271562, 0.270546, 0.269383, 0.268049, 0.266512, 0.264733, 0.262661, 0.260225, 0.257331, 0.253847, 0.24958, 0.244238, 0.241705, 0.238888, 0.235736, 0.232183, 0.228149, 0.223524, 0.218165, 0.211876, 0.204376, 0.195252, 0.183847, 0.169019, 0.159687, 0.148458, 0.134453, 0.125927, 0.115917, 0.111365, 0.106415, 0.100985, 0.0949599, 0.0881754, 0.080377, 0.0711342, 0.0657383, 0.0596046, 0.0524264, 0.04828, 0.043605, 0.0381837, 0.0369789, 0.0357239, 0.0344128, 0.0330384, 0.0315916, 0.0284334, 0.0247999, 0.0227304, 0.0204216}, 
{0.504534, 0.501693, 0.501238, 0.503711, 0.506007, 0.50649, 0.506633, 0.506559, 0.506365, 0.50612, 0.505173, 0.50512, 0.505381, 0.505519, 0.505371, 0.505318, 0.505354, 0.505311, 0.505262, 0.505236, 0.50521, 0.505177, 0.505141, 0.505059, 0.504978, 0.504812, 0.504645, 0.504476, 0.504135, 0.503786, 0.503432, 0.50307, 0.502702, 0.502326, 0.501944, 0.501554, 0.501156, 0.500751, 0.500338, 0.499917, 0.498828, 0.497686, 0.496488, 0.49523, 0.49391, 0.492523, 0.491065, 0.489531, 0.487913, 0.486206, 0.484399, 0.482484, 0.480448, 0.478277, 0.475955, 0.473463, 0.470778, 0.467872, 0.464713, 0.46126, 0.457463, 0.453259, 0.448569, 0.443287, 0.437269, 0.430315, 0.422131, 0.412262, 0.407689, 0.402664, 0.397103, 0.390894, 0.383895, 0.375919, 0.366715, 0.35594, 0.34311, 0.327518, 0.308058, 0.282827, 0.266998, 0.248008, 0.22441, 0.210086, 0.193307, 0.18569, 0.177413, 0.16834, 0.158283, 0.146966, 0.133965, 0.118563, 0.109572, 0.0993515, 0.0873891, 0.0804781, 0.0726848, 0.0636465, 0.0616378, 0.0595453, 0.0573593, 0.0550676, 0.0526555, 0.0473899, 0.0413327, 0.037883, 0.0340348}, 
{0.777682, 0.77817, 0.774907, 0.771079, 0.772456, 0.773571, 0.774249, 0.774423, 0.774236, 0.773865, 0.772678, 0.773376, 0.773489, 0.77299, 0.773229, 0.773121, 0.773005, 0.772942, 0.772818, 0.772776, 0.77268, 0.772579, 0.772506, 0.772294, 0.772084, 0.771654, 0.771224, 0.770787, 0.7699, 0.768995, 0.768072, 0.767129, 0.766168, 0.765187, 0.764186, 0.763165, 0.762123, 0.761061, 0.759977, 0.758872, 0.756015, 0.75302, 0.749885, 0.746606, 0.743182, 0.739611, 0.735891, 0.732017, 0.727985, 0.723788, 0.719419, 0.714865, 0.710114, 0.705148, 0.699946, 0.694481, 0.688721, 0.682628, 0.676155, 0.669243, 0.661823, 0.653808, 0.645091, 0.635536, 0.624969, 0.613154, 0.599764, 0.584301, 0.577364, 0.569885, 0.561756, 0.552827, 0.542888, 0.531648, 0.518699, 0.50348, 0.485224, 0.462864, 0.434823, 0.398473, 0.375741, 0.348572, 0.314974, 0.294663, 0.270945, 0.2602, 0.248542, 0.235779, 0.221649, 0.205769, 0.187547, 0.16598, 0.153396, 0.139094, 0.122354, 0.112682, 0.101773, 0.0891181, 0.0863053, 0.0833752, 0.0803139, 0.0771045, 0.0737264, 0.0663517, 0.0578682, 0.053037, 0.0476481}, 
{1.04431, 1.0466, 1.05065, 1.04828, 1.04719, 1.04809, 1.04894, 1.04925, 1.04906, 1.04863, 1.04812, 1.04862, 1.04806, 1.04843, 1.04808, 1.04824, 1.04803, 1.04799, 1.04788, 1.04781, 1.04768, 1.04762, 1.04752, 1.04729, 1.04704, 1.04657, 1.04609, 1.0456, 1.0446, 1.04357, 1.04251, 1.04141, 1.04028, 1.03911, 1.03791, 1.03667, 1.03539, 1.03406, 1.0327, 1.03129, 1.02756, 1.02353, 1.01918, 1.01448, 1.00943, 1.004, 0.9982, 0.992012, 0.985437, 0.978476, 0.971133, 0.963412, 0.955317, 0.946845, 0.937987, 0.928725, 0.919026, 0.908845, 0.898118, 0.886764, 0.874681, 0.861738, 0.847776, 0.832595, 0.815947, 0.797527, 0.776957, 0.7538, 0.74371, 0.733106, 0.721967, 0.710268, 0.697962, 0.684902, 0.670601, 0.65375, 0.632042, 0.603295, 0.56587, 0.517102, 0.486783, 0.450771, 0.406552, 0.379969, 0.349048, 0.33508, 0.319949, 0.303412, 0.285135, 0.264631, 0.241142, 0.213382, 0.197201, 0.178818, 0.157307, 0.144878, 0.130859, 0.114594, 0.110978, 0.107211, 0.103275, 0.099148, 0.0948041, 0.0853202, 0.0744092, 0.0681953, 0.0612642}, 
{1.32188, 1.31695, 1.31668, 1.31947, 1.31731, 1.3177, 1.31852, 1.31894, 1.31877, 1.31831, 1.31844, 1.31807, 1.31822, 1.318, 1.3182, 1.318, 1.31792, 1.31791, 1.31782, 1.3177, 1.31764, 1.31754, 1.31745, 1.31723, 1.31701, 1.31656, 1.31612, 1.31566, 1.31474, 1.31378, 1.3128, 1.31179, 1.31075, 1.30967, 1.30856, 1.3074, 1.3062, 1.30495, 1.30366, 1.30231, 1.2987, 1.2947, 1.29025, 1.2853, 1.2798, 1.27368, 1.26687, 1.25931, 1.25091, 1.24163, 1.23143, 1.22032, 1.20835, 1.19561, 1.18221, 1.16825, 1.15377, 1.13878, 1.1232, 1.10691, 1.08976, 1.07153, 1.05194, 1.03067, 1.0073, 0.981265, 0.951825, 0.91795, 0.902851, 0.886706, 0.869389, 0.850781, 0.830803, 0.809534, 0.787531, 0.766572, 0.750871, 0.748723, 0.741744, 0.732374, 0.728044, 0.723765, 0.719565, 0.717492, 0.715438, 0.714621, 0.713807, 0.712996, 0.712188, 0.711383, 0.710581, 0.709781, 0.709382, 0.708984, 0.708587, 0.708388, 0.70819, 0.707992, 0.707952, 0.707911, 0.707871, 0.707832, 0.707796, 0.70772, 0.70762, 0.707592, 0.707598}, 
{1.58204, 1.58635, 1.58318, 1.58511, 1.58355, 1.58343, 1.58418, 1.58468, 1.58448, 1.58398, 1.58441, 1.58383, 1.58417, 1.5841, 1.58387, 1.58384, 1.58373, 1.58367, 1.58357, 1.58342, 1.58336, 1.58323, 1.58315, 1.58289, 1.58264, 1.58213, 1.58164, 1.58114, 1.58016, 1.57919, 1.57821, 1.57723, 1.57624, 1.57524, 1.57422, 1.57317, 1.5721, 1.57099, 1.56984, 1.56865, 1.56545, 1.56189, 1.55789, 1.55337, 1.54825, 1.54244, 1.53582, 1.52823, 1.51949, 1.50932, 1.49738, 1.48328, 1.46673, 1.44791, 1.42768, 1.40711, 1.38682, 1.36683, 1.34689, 1.32662, 1.30561, 1.28347, 1.25976, 1.234, 1.20558, 1.17376, 1.13753, 1.09541, 1.07645, 1.05602, 1.03387, 1.0097, 0.983084, 0.953463, 0.919986, 0.881204, 0.833927, 0.770337, 0.707988, 0.641106, 0.601525, 0.555418, 0.499613, 0.466355, 0.427874, 0.41055, 0.391822, 0.371395, 0.348864, 0.323638, 0.2948, 0.260789, 0.240992, 0.218519, 0.192238, 0.177057, 0.159933, 0.140065, 0.135647, 0.131044, 0.126235, 0.121193, 0.115884, 0.104293, 0.0909554, 0.0833586, 0.0748845}, 
{1.84957, 1.84968, 1.85134, 1.85007, 1.84983, 1.84929, 1.85004, 1.8506, 1.85027, 1.84968, 1.85008, 1.85008, 1.84964, 1.84961, 1.84963, 1.84966, 1.84941, 1.84923, 1.84906, 1.84884, 1.84864, 1.84851, 1.84829, 1.84788, 1.84746, 1.84662, 1.84583, 1.84506, 1.84361, 1.84225, 1.84097, 1.83975, 1.83859, 1.83747, 1.83638, 1.8353, 1.83424, 1.83317, 1.8321, 1.831, 1.82815, 1.82502, 1.82153, 1.81757, 1.81306, 1.8079, 1.80199, 1.79518, 1.78728, 1.77801, 1.76686, 1.75278, 1.73323, 1.70384, 1.66882, 1.63893, 1.61423, 1.59206, 1.57064, 1.54892, 1.52619, 1.5019, 1.47555, 1.44657, 1.41428, 1.37781, 1.33598, 1.28712, 1.26509, 1.24133, 1.21558, 1.18751, 1.15671, 1.12263, 1.08453, 1.04141, 0.991742, 0.933178, 0.861906, 0.773444, 0.721837, 0.663505, 0.594647, 0.554166, 0.507666, 0.486825, 0.464348, 0.439888, 0.412969, 0.382898, 0.348601, 0.308246, 0.284798, 0.25821, 0.227147, 0.209213, 0.188989, 0.165524, 0.160306, 0.15487, 0.149189, 0.143232, 0.136961, 0.123266, 0.107504, 0.113694, 0.102135}};

  /* 2, 1 mode */

  const double reomegaqnm21[8][107] = {{0.336609, 0.339386, 0.340852, 0.342219, 0.342818, 0.343002, 
  0.343143, 0.343254, 0.343344, 0.343417, 0.343641, 0.343748, 
  0.343804, 0.343853, 0.343871, 0.343879, 0.343886, 0.343891, 
  0.343896, 0.343902, 0.343909, 0.343915, 0.343922, 0.343941, 0.34396,
   0.344002, 0.344049, 0.344101, 0.34422, 0.344359, 0.344517, 
  0.344696, 0.344896, 0.345115, 0.345356, 0.345617, 0.345899, 
  0.346201, 0.346525, 0.34687, 0.347824, 0.348911, 0.350132, 0.351491,
   0.35299, 0.354633, 0.356423, 0.358366, 0.360469, 0.362738, 
  0.365183, 0.367812, 0.370637, 0.373672, 0.376931, 0.380432, 
  0.384197, 0.388248, 0.392615, 0.39733, 0.402436, 0.407979, 0.41402, 
  0.420632, 0.427909, 0.435968, 0.444968, 0.455121, 0.459569, 
  0.464271, 0.469259, 0.474564, 0.480231, 0.486308, 0.492859, 
  0.499965, 0.507729, 0.516291, 0.525845, 0.536673, 0.542693, 
  0.549213, 0.556329, 0.560146, 0.564155, 0.565814, 0.567505, 
  0.569227, 0.570976, 0.572749, 0.574535, 0.576322, 0.577208, 
  0.578084, 0.578948, 0.579374, 0.579795, 0.580212, 0.580295, 
  0.580377, 0.58046, 0.580541, 0.580623, 0.580784, 0.580942, 0.581018,
   0.581093}, {0.313486, 0.314984, 0.315748, 0.31636, 0.316504, 
  0.316512, 0.316502, 0.316483, 0.31646, 0.316436, 0.316347, 0.316311,
   0.316305, 0.316315, 0.316325, 0.31633, 0.316338, 0.316345, 
  0.316352, 0.31636, 0.316368, 0.316375, 0.316383, 0.316401, 0.316419,
   0.316455, 0.316491, 0.316527, 0.316602, 0.316682, 0.316772, 
  0.316874, 0.316989, 0.317122, 0.317273, 0.317444, 0.317636, 
  0.317851, 0.31809, 0.318353, 0.319124, 0.320064, 0.321179, 0.322476,
   0.32396, 0.325636, 0.327508, 0.329582, 0.331865, 0.334364, 
  0.337088, 0.340045, 0.343249, 0.346711, 0.350448, 0.354478, 
  0.358823, 0.363506, 0.368558, 0.374014, 0.379916, 0.386313, 
  0.393268, 0.400856, 0.409172, 0.418337, 0.428509, 0.439897, 
  0.444855, 0.450076, 0.455587, 0.46142, 0.467612, 0.474208, 0.481261,
   0.488834, 0.497004, 0.505862, 0.515506, 0.526011, 0.531572, 
  0.537263, 0.542874, 0.545479, 0.547723, 0.548438, 0.548984, 0.54929,
   0.549241, 0.548642, 0.54714, 0.544185, 0.542606, 0.542325, 
  0.541626, 0.540962, 0.540617, 0.540156, 0.540066, 0.539986, 
  0.539906, 0.539818, 0.539726, 0.539553, 0.539377, 0.539288, 
  0.539198}, {0.270281, 0.269419, 0.268655, 0.267935, 0.267923, 
  0.268009, 0.268102, 0.268185, 0.268251, 0.268297, 0.268334, 
  0.268285, 0.268273, 0.268312, 0.268342, 0.268359, 0.268398, 
  0.268442, 0.268483, 0.268523, 0.268562, 0.268601, 0.26864, 0.268733,
   0.268822, 0.268989, 0.269141, 0.26928, 0.26952, 0.269718, 0.269881,
   0.270016, 0.270132, 0.270235, 0.270333, 0.270431, 0.270536, 
  0.270653, 0.270785, 0.270937, 0.271426, 0.2721, 0.272992, 0.274121, 
  0.275504, 0.277153, 0.279078, 0.281288, 0.283792, 0.286598, 
  0.289716, 0.293156, 0.296931, 0.301053, 0.30554, 0.31041, 0.315684, 
  0.321389, 0.327554, 0.334218, 0.341421, 0.349217, 0.357668, 
  0.366852, 0.376864, 0.387825, 0.39989, 0.413261, 0.419032, 0.425075,
   0.431414, 0.438075, 0.445088, 0.452488, 0.46031, 0.468593, 
  0.477369, 0.48665, 0.49639, 0.506348, 0.511184, 0.515607, 0.51903, 
  0.519992, 0.5201, 0.51983, 0.519361, 0.51871, 0.517943, 0.517204, 
  0.516728, 0.516696, 0.516098, 0.513376, 0.509813, 0.508083, 
  0.506412, 0.504789, 0.504468, 0.504148, 0.50383, 0.503511, 0.503194,
   0.502559, 0.501923, 0.501605, 0.501286}, {0.21269, 0.214089, 
  0.215208, 0.215364, 0.21474, 0.214584, 0.21456, 0.214625, 0.21473, 
  0.214837, 0.214992, 0.214885, 0.214922, 0.215045, 0.215071, 
  0.215139, 0.215267, 0.215383, 0.215507, 0.215624, 0.215741, 
  0.215858, 0.215974, 0.216256, 0.216533, 0.217064, 0.217566, 0.21804,
   0.218903, 0.219653, 0.220295, 0.220836, 0.221285, 0.22165, 
  0.221945, 0.22218, 0.222368, 0.22252, 0.222647, 0.222759, 0.223032, 
  0.223383, 0.2239, 0.224647, 0.225672, 0.22701, 0.228689, 0.23073, 
  0.233151, 0.235969, 0.239199, 0.242854, 0.246951, 0.251505, 
  0.256534, 0.262058, 0.268101, 0.274689, 0.281852, 0.289628, 
  0.298059, 0.307197, 0.317104, 0.327855, 0.339542, 0.352278, 
  0.366207, 0.38151, 0.388065, 0.394895, 0.402017, 0.409454, 0.417227,
   0.425361, 0.433878, 0.442797, 0.452127, 0.461855, 0.471918, 
  0.482195, 0.48739, 0.49271, 0.498495, 0.501813, 0.5056, 0.507251, 
  0.508939, 0.510585, 0.512031, 0.513028, 0.513289, 0.512758, 
  0.512367, 0.512091, 0.511921, 0.51135, 0.5096, 0.50715, 0.50665, 
  0.506154, 0.505662, 0.505175, 0.504694, 0.503744, 0.502809, 
  0.502344, 0.50188}, {0.171239, 0.169224, 0.168634, 0.170053, 
  0.169938, 0.169613, 0.169418, 0.169408, 0.169531, 0.169696, 
  0.169825, 0.169724, 0.16992, 0.169916, 0.170075, 0.17014, 0.170344, 
  0.170511, 0.17071, 0.170886, 0.171073, 0.171258, 0.171437, 0.171892,
   0.172337, 0.173213, 0.174064, 0.174889, 0.176453, 0.177895, 
  0.179203, 0.180372, 0.181395, 0.182273, 0.183007, 0.183606, 
  0.184081, 0.184444, 0.18471, 0.184896, 0.185105, 0.18512, 0.185113, 
  0.185222, 0.18555, 0.186176, 0.187164, 0.188562, 0.190413, 0.192751,
   0.195604, 0.199, 0.202962, 0.207515, 0.212679, 0.21848, 0.224939, 
  0.232083, 0.23994, 0.248543, 0.257928, 0.268139, 0.279228, 0.291258,
   0.304303, 0.318456, 0.33383, 0.350563, 0.357675, 0.365045, 
  0.372687, 0.380617, 0.388855, 0.397422, 0.406346, 0.41567, 0.425462,
   0.435857, 0.447149, 0.45999, 0.467349, 0.475546, 0.484594, 
  0.489367, 0.494295, 0.496339, 0.498459, 0.500692, 0.50308, 0.505635,
   0.508214, 0.510197, 0.510572, 0.510387, 0.509864, 0.509635, 
  0.509456, 0.50879, 0.508443, 0.507999, 0.50747, 0.506881, 0.506257, 
  0.504974, 0.503704, 0.503082, 0.502466}, {0.136159, 0.137914, 
  0.137016, 0.136626, 0.137329, 0.137055, 0.136755, 0.136667, 
  0.136793, 0.136997, 0.136954, 0.137117, 0.137154, 0.137279, 0.1373, 
  0.13747, 0.137649, 0.137842, 0.138069, 0.138258, 0.138476, 0.138671,
   0.138881, 0.139389, 0.139898, 0.14091, 0.14191, 0.142898, 0.144829,
   0.146684, 0.148447, 0.150099, 0.151622, 0.153001, 0.154221, 
  0.15527, 0.156142, 0.156839, 0.157368, 0.15774, 0.158094, 0.157845, 
  0.15725, 0.156526, 0.155847, 0.155354, 0.155162, 0.155369, 0.156061,
   0.157316, 0.159206, 0.161795, 0.165143, 0.169299, 0.174306, 
  0.180192, 0.186977, 0.194671, 0.203282, 0.212811, 0.223263, 
  0.234649, 0.246985, 0.260297, 0.274622, 0.290009, 0.306518, 
  0.324226, 0.331666, 0.339322, 0.347205, 0.355332, 0.363727, 
  0.372432, 0.381518, 0.391111, 0.401435, 0.412868, 0.425995, 
  0.441563, 0.45052, 0.460463, 0.471717, 0.477974, 0.484665, 0.487438,
   0.490248, 0.493085, 0.49596, 0.498923, 0.502084, 0.505527, 
  0.507222, 0.508526, 0.508875, 0.508627, 0.508251, 0.507926, 
  0.507863, 0.507785, 0.50767, 0.507484, 0.507181, 0.506113, 0.504624,
   0.503839, 0.50306}, {0.114342, 0.113963, 0.115141, 0.114074, 
  0.11486, 0.114789, 0.114485, 0.114338, 0.114452, 0.114665, 0.114512,
   0.114807, 0.11467, 0.114821, 0.114995, 0.115025, 0.115231, 
  0.115438, 0.115661, 0.115846, 0.11605, 0.116249, 0.116452, 0.116951,
   0.117454, 0.118454, 0.11946, 0.120464, 0.122464, 0.124443, 
  0.126383, 0.128266, 0.130073, 0.131781, 0.133368, 0.134813, 0.13609,
   0.137179, 0.138065, 0.138742, 0.139567, 0.139335, 0.138333, 
  0.136834, 0.135064, 0.133208, 0.131418, 0.129826, 0.128564, 
  0.127775, 0.127629, 0.128336, 0.130133, 0.133252, 0.137863, 
  0.144021, 0.151672, 0.160684, 0.170902, 0.182186, 0.194429, 
  0.207556, 0.221528, 0.236327, 0.251955, 0.268421, 0.285744, 0.30395,
   0.311492, 0.319195, 0.327076, 0.335168, 0.343522, 0.352232, 
  0.361443, 0.371382, 0.382361, 0.394771, 0.409089, 0.426092, 
  0.436065, 0.4474, 0.460415, 0.467675, 0.475563, 0.478922, 0.482402, 
  0.485993, 0.489665, 0.493377, 0.497113, 0.500995, 0.503086, 
  0.505273, 0.507186, 0.507677, 0.507628, 0.507177, 0.507077, 
  0.506981, 0.506892, 0.506806, 0.506714, 0.506383, 0.505391, 
  0.504576, 0.503665}, {0.0996362, 0.098809, 0.0988228, 0.0989612, 
  0.0991812, 0.0993225, 0.0990912, 0.0989133, 0.0990066, 0.0992056, 
  0.0990844, 0.0992013, 0.0993046, 0.0993698, 0.0994135, 0.0995406, 
  0.0997102, 0.0999057, 0.100109, 0.100292, 0.100464, 0.100668, 
  0.100843, 0.10132, 0.101793, 0.102742, 0.103705, 0.104671, 0.106617,
   0.108576, 0.11054, 0.112492, 0.114417, 0.116296, 0.118106, 
  0.119822, 0.121422, 0.122876, 0.12415, 0.125214, 0.12686, 0.127085, 
  0.126123, 0.124264, 0.121765, 0.118824, 0.115576, 0.112102, 0.10844,
   0.104598, 0.100602, 0.0966436, 0.0934548, 0.0928224, 0.0967608, 
  0.105267, 0.116748, 0.129817, 0.143719, 0.158124, 0.17292, 0.188094,
   0.203673, 0.219697, 0.236201, 0.253208, 0.270727, 0.288784, 
  0.296181, 0.303707, 0.311404, 0.319342, 0.327628, 0.336416, 
  0.345907, 0.35633, 0.367923, 0.380968, 0.395965, 0.413866, 0.424377,
   0.436366, 0.450442, 0.458506, 0.467375, 0.471171, 0.475126, 
  0.479254, 0.483568, 0.488052, 0.492639, 0.497232, 0.499557, 
  0.501993, 0.504598, 0.505838, 0.506697, 0.506742, 0.506649, 
  0.506536, 0.50641, 0.506281, 0.506156, 0.505924, 0.505541, 0.505047,
   0.504219}};

  const double imomegaqnm21[8][107] = {{0.0743116, 0.0772659, 0.0791824, 0.0812815, 0.0822652, 0.082556, 
  0.0827684, 0.082926, 0.0830445, 0.0831347, 0.0833596, 0.0834304, 
  0.0834555, 0.0834714, 0.0834806, 0.0834917, 0.0835195, 0.0835507, 
  0.0835831, 0.0836156, 0.0836481, 0.0836804, 0.0837126, 0.083792, 
  0.0838704, 0.0840241, 0.0841737, 0.0843193, 0.0845994, 0.084865, 
  0.085117, 0.0853561, 0.0855831, 0.0857987, 0.0860035, 0.086198, 
  0.0863828, 0.0865584, 0.0867253, 0.086884, 0.087247, 0.0875662, 
  0.0878462, 0.0880906, 0.0883025, 0.088484, 0.088637, 0.0887628, 
  0.0888621, 0.0889354, 0.0889828, 0.0890037, 0.0889973, 0.0889623, 
  0.0888968, 0.0887983, 0.0886636, 0.0884885, 0.0882679, 0.0879952, 
  0.0876618, 0.0872571, 0.086767, 0.086173, 0.0854501, 0.0845642, 
  0.0834665, 0.0820852, 0.0814304, 0.0807035, 0.0798926, 0.0789827, 
  0.077955, 0.0767847, 0.0754394, 0.0738744, 0.072027, 0.0698043, 
  0.067061, 0.0635492, 0.0613721, 0.058791, 0.0556436, 0.0537759, 
  0.0516427, 0.0506978, 0.0496908, 0.0486136, 0.0474565, 0.0462084, 
  0.0448564, 0.0433879, 0.0426068, 0.0417939, 0.0409501, 0.0405173, 
  0.0400777, 0.0396319, 0.039542, 0.0394519, 0.0393616, 0.039271, 
  0.0391802, 0.0389974, 0.0388111, 0.038714, 0.0386091}, {0.258785, 
  0.258526, 0.258215, 0.257703, 0.257385, 0.257284, 0.25721, 0.257158,
   0.257123, 0.257101, 0.257089, 0.257133, 0.257178, 0.257248, 
  0.257304, 0.257357, 0.257465, 0.257574, 0.257683, 0.257791, 
  0.257898, 0.258005, 0.258112, 0.258376, 0.258637, 0.259149, 
  0.259649, 0.260137, 0.261076, 0.261967, 0.262812, 0.263613, 
  0.264372, 0.265089, 0.265768, 0.26641, 0.267016, 0.267589, 0.26813, 
  0.26864, 0.269792, 0.270782, 0.271629, 0.272346, 0.272945, 0.273434,
   0.273821, 0.274112, 0.274309, 0.274415, 0.27443, 0.274354, 
  0.274183, 0.273915, 0.273543, 0.273058, 0.272452, 0.271711, 
  0.270819, 0.269754, 0.268489, 0.266992, 0.265218, 0.263108, 
  0.260587, 0.257547, 0.253838, 0.249238, 0.247078, 0.244692, 
  0.242045, 0.239089, 0.235766, 0.231999, 0.227686, 0.222686, 0.2168, 
  0.209731, 0.201006, 0.189803, 0.18282, 0.174484, 0.164205, 0.158029,
   0.150891, 0.1477, 0.14428, 0.140606, 0.136657, 0.132444, 0.128112, 
  0.124511, 0.123925, 0.123321, 0.12167, 0.121162, 0.120732, 0.12014, 
  0.120048, 0.119951, 0.119844, 0.119734, 0.119631, 0.119429, 
  0.119227, 0.119127, 0.119027}, {0.451082, 0.450126, 0.449946, 
  0.450393, 0.450836, 0.450954, 0.451016, 0.451037, 0.451032, 
  0.451015, 0.450927, 0.450953, 0.451028, 0.451151, 0.451242, 
  0.451336, 0.451538, 0.451736, 0.451932, 0.452129, 0.452326, 
  0.452521, 0.452716, 0.453201, 0.453683, 0.454634, 0.45557, 0.456489,
   0.458275, 0.459989, 0.461626, 0.463185, 0.464665, 0.466067, 
  0.46739, 0.468636, 0.469808, 0.470906, 0.471935, 0.472897, 0.475022,
   0.476783, 0.478216, 0.479352, 0.48022, 0.48084, 0.481231, 0.481407,
   0.481377, 0.481149, 0.480726, 0.480109, 0.479294, 0.478277, 
  0.477047, 0.475592, 0.473893, 0.471929, 0.469669, 0.467076, 
  0.464104, 0.460691, 0.456761, 0.45221, 0.446903, 0.440653, 0.433196,
   0.424144, 0.419951, 0.415358, 0.4103, 0.404693, 0.398434, 0.391383,
   0.383355, 0.374094, 0.363227, 0.350188, 0.334044, 0.313093, 
  0.299808, 0.283583, 0.262719, 0.249472, 0.233119, 0.225306, 
  0.216483, 0.206365, 0.194563, 0.180538, 0.16345, 0.141483, 0.12731, 
  0.111596, 0.0954876, 0.0868316, 0.0774193, 0.066857, 0.0645541, 
  0.0621712, 0.0596988, 0.057125, 0.0544355, 0.0486318, 0.0420655, 
  0.0383765, 0.0343031}, {0.670736, 0.672092, 0.671812, 0.670663, 
  0.670549, 0.67073, 0.670927, 0.671078, 0.671164, 0.671189, 0.671021,
   0.671075, 0.671219, 0.67131, 0.671427, 0.671572, 0.671807, 
  0.672062, 0.672309, 0.672557, 0.672808, 0.673057, 0.673305, 
  0.673929, 0.674552, 0.675799, 0.677043, 0.678283, 0.680743, 
  0.683166, 0.685536, 0.68784, 0.690067, 0.692205, 0.694247, 0.696187,
   0.69802, 0.699744, 0.701359, 0.702865, 0.706165, 0.708833, 
  0.710915, 0.71246, 0.71351, 0.714105, 0.714277, 0.714053, 0.713452, 
  0.71249, 0.711176, 0.709516, 0.707509, 0.705148, 0.702423, 0.699317,
   0.695803, 0.691851, 0.687416, 0.682445, 0.676867, 0.670592, 
  0.663504, 0.655451, 0.646227, 0.635554, 0.623036, 0.60809, 0.601244,
   0.593791, 0.58563, 0.576638, 0.566651, 0.555455, 0.54276, 0.528156,
   0.511038, 0.490465, 0.464847, 0.431189, 0.409563, 0.382936, 
  0.348824, 0.32773, 0.302947, 0.291759, 0.279694, 0.26659, 0.252176, 
  0.235931, 0.216835, 0.192946, 0.178198, 0.160681, 0.138951, 
  0.125667, 0.110911, 0.095128, 0.0917732, 0.0883219, 0.0847576, 
  0.0810612, 0.0772103, 0.0689279, 0.0595841, 0.0543422, 
  0.048559}, {0.910592, 0.910063, 0.91138, 0.911832, 0.911012, 
  0.911003, 0.911199, 0.911439, 0.911606, 0.91166, 0.911352, 0.911585,
   0.911648, 0.911704, 0.911854, 0.911933, 0.912192, 0.912425, 
  0.912663, 0.912897, 0.913141, 0.913376, 0.913616, 0.91422, 0.91483, 
  0.916061, 0.917309, 0.918574, 0.921148, 0.923769, 0.92642, 0.929083,
   0.931738, 0.934365, 0.936941, 0.939446, 0.941862, 0.944173, 
  0.946367, 0.948436, 0.953025, 0.956761, 0.959656, 0.961747, 0.96308,
   0.963701, 0.963651, 0.962967, 0.96168, 0.959815, 0.957389, 
  0.954417, 0.950902, 0.946845, 0.942236, 0.937057, 0.931283, 
  0.924873, 0.917775, 0.909918, 0.901208, 0.891524, 0.880707, 
  0.868545, 0.854753, 0.838939, 0.820547, 0.798756, 0.788822, 
  0.778034, 0.76625, 0.75329, 0.738922, 0.722836, 0.704612, 0.683652, 
  0.659084, 0.629562, 0.592913, 0.54542, 0.515651, 0.480063, 0.436051,
   0.409334, 0.377864, 0.363486, 0.347803, 0.330563, 0.311454, 
  0.290066, 0.265802, 0.237412, 0.220683, 0.201065, 0.17693, 0.162373,
   0.145371, 0.124789, 0.120176, 0.115422, 0.110541, 0.105527, 
  0.100358, 0.0893955, 0.0771778, 0.0703585, 0.0628494}, {1.16313, 
  1.16205, 1.16102, 1.16236, 1.1618, 1.16155, 1.16165, 1.16191, 
  1.16212, 1.16217, 1.16185, 1.16214, 1.162, 1.16222, 1.16228, 1.1624,
   1.16261, 1.16279, 1.163, 1.16321, 1.16342, 1.16362, 1.16384, 
  1.16436, 1.16489, 1.16598, 1.1671, 1.16824, 1.1706, 1.17306, 
  1.17562, 1.17826, 1.18098, 1.18376, 1.18658, 1.18941, 1.19222, 
  1.19499, 1.19769, 1.20029, 1.20624, 1.21127, 1.21528, 1.21824, 
  1.22017, 1.2211, 1.22105, 1.22006, 1.21814, 1.21533, 1.21164, 
  1.20711, 1.20176, 1.19561, 1.18868, 1.18098, 1.17251, 1.16324, 
  1.15314, 1.14213, 1.13009, 1.11686, 1.10224, 1.08592, 1.06753, 
  1.04652, 1.02214, 0.993243, 0.980058, 0.965724, 0.950047, 0.932779, 
  0.913601, 0.892095, 0.867699, 0.839644, 0.806856, 0.767797, 
  0.720146, 0.659944, 0.622848, 0.578732, 0.524237, 0.491317, 
  0.452927, 0.435542, 0.416661, 0.395943, 0.372912, 0.346898, 
  0.316949, 0.28167, 0.261342, 0.238436, 0.211307, 0.195088, 0.176156,
   0.153337, 0.148132, 0.142654, 0.136869, 0.130747, 0.124273, 
  0.110319, 0.0949385, 0.0864639, 0.0771851}, {1.41428, 1.41604, 
  1.4156, 1.41536, 1.4156, 1.41524, 1.41519, 1.41543, 1.41565, 
  1.41568, 1.41551, 1.41554, 1.41563, 1.41566, 1.41579, 1.41591, 
  1.41608, 1.41624, 1.41643, 1.41663, 1.4168, 1.417, 1.41717, 1.41765,
   1.41813, 1.41911, 1.42011, 1.42114, 1.42328, 1.42554, 1.42792, 
  1.43041, 1.43303, 1.43576, 1.4386, 1.44152, 1.44453, 1.44758, 
  1.45065, 1.45369, 1.46099, 1.46752, 1.473, 1.47733, 1.48042, 
  1.48226, 1.48282, 1.48207, 1.47997, 1.47646, 1.4715, 1.46505, 
  1.45715, 1.44791, 1.43752, 1.42619, 1.41412, 1.4014, 1.38801, 
  1.37384, 1.35872, 1.34237, 1.32445, 1.30456, 1.28213, 1.25643, 
  1.22645, 1.19065, 1.17422, 1.15628, 1.13659, 1.11481, 1.09053, 
  1.06323, 1.03225, 0.996712, 0.955487, 0.906959, 0.848541, 0.775494, 
  0.730792, 0.678016, 0.613349, 0.574379, 0.528885, 0.508284, 
  0.485946, 0.461517, 0.434493, 0.404108, 0.369135, 0.327547, 
  0.303278, 0.275875, 0.244209, 0.225953, 0.204982, 0.179727, 
  0.173965, 0.167908, 0.161519, 0.154752, 0.147552, 0.131555, 
  0.113035, 0.102763, 0.09161}, {1.6695, 1.6684, 1.66944, 1.66866, 
  1.66927, 1.66897, 1.66882, 1.66899, 1.6692, 1.66921, 1.66919, 
  1.66906, 1.66925, 1.66933, 1.66938, 1.66941, 1.6696, 1.66977, 
  1.66995, 1.67014, 1.67031, 1.67048, 1.67066, 1.6711, 1.67155, 
  1.67246, 1.6734, 1.67436, 1.67635, 1.67846, 1.68069, 1.68303, 
  1.68552, 1.68814, 1.69089, 1.69379, 1.69681, 1.69996, 1.70322, 
  1.70655, 1.71494, 1.72293, 1.73008, 1.73612, 1.7409, 1.74437, 
  1.74645, 1.74707, 1.74612, 1.74337, 1.73845, 1.73073, 1.71933, 
  1.70384, 1.68568, 1.6674, 1.65022, 1.63402, 1.61822, 1.60219, 
  1.58534, 1.56716, 1.54712, 1.52462, 1.49894, 1.46917, 1.43397, 
  1.39142, 1.37171, 1.35009, 1.32623, 1.29974, 1.27014, 1.23683, 
  1.19912, 1.15609, 1.10649, 1.04841, 0.97878, 0.892317, 0.839828, 
  0.778168, 0.702974, 0.657884, 0.605411, 0.581667, 0.555906, 
  0.527711, 0.496517, 0.46151, 0.421398, 0.373867, 0.346042, 0.314367,
   0.277467, 0.25637, 0.232725, 0.204918, 0.198597, 0.191947, 
  0.184925, 0.177485, 0.169568, 0.151985, 0.131247, 0.119312, 
  0.106197}};

  /* 3, 3 mode */
  const double reomegaqnm33[8][107] = {{0.445768, 0.452799, 0.456948, 0.460943, 0.462462, 0.462842, 
  0.463095, 0.463269, 0.463394, 0.463488, 0.463746, 0.463886, 
  0.463989, 0.464144, 0.464267, 0.464374, 0.464572, 0.464763, 
  0.464952, 0.46514, 0.465329, 0.465518, 0.465707, 0.466182, 0.466657,
   0.467612, 0.468573, 0.46954, 0.471491, 0.473465, 0.475464, 
  0.477487, 0.479535, 0.481609, 0.483709, 0.485837, 0.487991, 
  0.490174, 0.492386, 0.494627, 0.500363, 0.5063, 0.512449, 0.518826, 
  0.525445, 0.532323, 0.539479, 0.546934, 0.55471, 0.562834, 0.571335,
   0.580244, 0.5896, 0.599443, 0.609823, 0.620796, 0.632425, 0.644787,
   0.657972, 0.672086, 0.68726, 0.70365, 0.721455, 0.740921, 0.762369,
   0.786223, 0.813057, 0.843687, 0.857254, 0.871717, 0.887201, 
  0.90386, 0.921885, 0.941521, 0.963088, 0.987016, 1.01391, 1.04464, 
  1.08058, 1.1241, 1.14998, 1.17986, 1.21547, 1.23637, 1.26023, 
  1.27086, 1.28227, 1.29462, 1.30812, 1.32308, 1.33999, 1.35965, 
  1.37094, 1.38363, 1.39829, 1.40666, 1.41603, 1.42679, 1.42917, 
  1.43164, 1.43422, 1.43692, 1.43975, 1.44591, 1.45295, 1.45695, 
  1.46139}, {0.428357, 0.432238, 0.434351, 0.436362, 0.437184, 
  0.437404, 0.437552, 0.437654, 0.437725, 0.437773, 0.437868, 
  0.437898, 0.437929, 0.438015, 0.438115, 0.438217, 0.438421, 
  0.438623, 0.438826, 0.43903, 0.439233, 0.439437, 0.439641, 0.440152,
   0.440665, 0.441694, 0.442729, 0.443771, 0.445872, 0.447999, 
  0.45015, 0.452328, 0.454532, 0.456764, 0.459023, 0.461311, 0.463627,
   0.465974, 0.46835, 0.470758, 0.476918, 0.483289, 0.489884, 
  0.496717, 0.503805, 0.511164, 0.518814, 0.526776, 0.535074, 
  0.543734, 0.552785, 0.562261, 0.5722, 0.582644, 0.593642, 0.60525, 
  0.617535, 0.630573, 0.644453, 0.659285, 0.675198, 0.692352, 
  0.710943, 0.731221, 0.753507, 0.778225, 0.805952, 0.837504, 
  0.851449, 0.866295, 0.882167, 0.89922, 0.917645, 0.937687, 0.959667,
   0.984015, 1.01133, 1.0425, 1.07889, 1.12285, 1.14897, 1.17907, 
  1.21491, 1.23591, 1.25988, 1.27055, 1.28201, 1.2944, 1.30794, 
  1.32294, 1.33989, 1.35958, 1.37089, 1.38359, 1.39826, 1.40664, 
  1.41601, 1.42678, 1.42916, 1.43163, 1.43421, 1.43691, 1.43974, 
  1.4459, 1.45295, 1.45695, 1.46139}, {0.390162, 0.390933, 0.390855, 
  0.39013, 0.389593, 0.389452, 0.389377, 0.389349, 0.389352, 0.389374,
   0.389533, 0.389627, 0.38968, 0.38978, 0.389896, 0.390013, 0.390243,
   0.390473, 0.390704, 0.390935, 0.391167, 0.391398, 0.39163, 
  0.392211, 0.392793, 0.393963, 0.395139, 0.396321, 0.398706, 
  0.401118, 0.403557, 0.406025, 0.408522, 0.411048, 0.413604, 
  0.416191, 0.41881, 0.42146, 0.424144, 0.426861, 0.433807, 0.44098, 
  0.448397, 0.456071, 0.464019, 0.47226, 0.480815, 0.489704, 0.498954,
   0.50859, 0.518645, 0.529153, 0.540151, 0.551685, 0.563804, 
  0.576567, 0.590041, 0.604302, 0.619445, 0.635578, 0.652833, 
  0.671371, 0.691391, 0.713145, 0.736956, 0.763249, 0.792606, 
  0.825846, 0.840483, 0.856032, 0.872619, 0.890397, 0.90956, 0.930355,
   0.953102, 0.978234, 1.00636, 1.03835, 1.07558, 1.12042, 1.14697, 
  1.17752, 1.21379, 1.23502, 1.2592, 1.26995, 1.28149, 1.29396, 
  1.30758, 1.32266, 1.33967, 1.35943, 1.37078, 1.38351, 1.3982, 
  1.4066, 1.41598, 1.42676, 1.42914, 1.43162, 1.4342, 1.4369, 1.43973,
   1.44589, 1.45294, 1.45694, 1.46139}, {0.332321, 0.329612, 0.329447,
   0.330998, 0.331533, 0.331485, 0.331367, 0.331244, 0.331151, 
  0.331101, 0.331226, 0.331354, 0.331387, 0.331504, 0.331643, 
  0.331766, 0.332024, 0.33228, 0.332537, 0.332795, 0.333053, 0.333311,
   0.333569, 0.334216, 0.334865, 0.336168, 0.337478, 0.338796, 
  0.341454, 0.344143, 0.346863, 0.349614, 0.352398, 0.355215, 
  0.358066, 0.360951, 0.363872, 0.366829, 0.369823, 0.372855, 
  0.380605, 0.388612, 0.39689, 0.405455, 0.414327, 0.423525, 0.43307, 
  0.442986, 0.4533, 0.464039, 0.475237, 0.486928, 0.499155, 0.511962, 
  0.525402, 0.539534, 0.554428, 0.570163, 0.586836, 0.604556, 0.62346,
   0.643711, 0.665511, 0.689115, 0.71485, 0.743145, 0.774586, 
  0.809998, 0.82553, 0.841991, 0.859506, 0.878231, 0.898362, 0.920146,
   0.943908, 0.970086, 0.999286, 1.0324, 1.0708, 1.11686, 1.14404, 
  1.17523, 1.21215, 1.23369, 1.25818, 1.26906, 1.28072, 1.2933, 
  1.30704, 1.32223, 1.33936, 1.35922, 1.37061, 1.38339, 1.39812, 
  1.40653, 1.41593, 1.42673, 1.42911, 1.43159, 1.43418, 1.43688, 
  1.43971, 1.44588, 1.45294, 1.45694, 1.46138}, {0.279748, 0.281425, 
  0.280003, 0.278496, 0.279479, 0.279701, 0.279679, 0.279534, 0.27938,
   0.279284, 0.279511, 0.279556, 0.279569, 0.279745, 0.279852, 
  0.279994, 0.280246, 0.280507, 0.280766, 0.281026, 0.281285, 
  0.281546, 0.281806, 0.282459, 0.283115, 0.284431, 0.285757, 
  0.287091, 0.289787, 0.292519, 0.295288, 0.298095, 0.30094, 0.303825,
   0.306751, 0.309717, 0.312725, 0.315776, 0.318871, 0.322011, 
  0.330062, 0.338415, 0.347087, 0.356096, 0.365463, 0.37521, 0.385359,
   0.395936, 0.406968, 0.418488, 0.430528, 0.443126, 0.456326, 
  0.470174, 0.484725, 0.500042, 0.516195, 0.533269, 0.551359, 
  0.570583, 0.591077, 0.61301, 0.636588, 0.662069, 0.689786, 0.720173,
   0.753822, 0.791564, 0.808064, 0.825515, 0.844044, 0.863808, 
  0.885005, 0.907885, 0.932779, 0.960131, 0.990559, 1.02497, 1.06475, 
  1.11229, 1.14026, 1.17226, 1.21, 1.23195, 1.25686, 1.26789, 1.2797, 
  1.29245, 1.30633, 1.32167, 1.33893, 1.35893, 1.37038, 1.38322, 
  1.39801, 1.40645, 1.41587, 1.42669, 1.42908, 1.43156, 1.43415, 
  1.43685, 1.43969, 1.44587, 1.45293, 1.45693, 1.46137}, {0.239992, 
  0.240995, 0.242946, 0.24156, 0.241854, 0.242213, 0.242269, 0.242117,
   0.241939, 0.241855, 0.242174, 0.24208, 0.242226, 0.242298, 
  0.242438, 0.242546, 0.242783, 0.243024, 0.243261, 0.243498, 
  0.243736, 0.243974, 0.244213, 0.244811, 0.245412, 0.246621, 0.24784,
   0.249069, 0.251558, 0.254089, 0.256662, 0.259279, 0.26194, 
  0.264647, 0.2674, 0.2702, 0.273049, 0.275948, 0.278897, 0.281899, 
  0.289636, 0.297725, 0.306186, 0.315041, 0.324312, 0.334027, 0.34421,
   0.354891, 0.366101, 0.377873, 0.390245, 0.403257, 0.416954, 
  0.431386, 0.446612, 0.462696, 0.479713, 0.49775, 0.516906, 0.537303,
   0.559082, 0.582416, 0.607517, 0.634649, 0.66415, 0.696462, 
  0.732182, 0.77215, 0.789586, 0.807999, 0.827516, 0.848298, 0.870543,
   0.894505, 0.920521, 0.949042, 0.980704, 1.01643, 1.05767, 1.10684, 
  1.1357, 1.16864, 1.20737, 1.22983, 1.25522, 1.26646, 1.27846, 
  1.29139, 1.30546, 1.32098, 1.33842, 1.35857, 1.3701, 1.38301, 
  1.39787, 1.40634, 1.41579, 1.42664, 1.42903, 1.43152, 1.43411, 
  1.43682, 1.43966, 1.44585, 1.45291, 1.45692, 1.46137}, {0.21884, 
  0.216298, 0.216633, 0.217244, 0.216888, 0.217233, 0.217329, 
  0.217189, 0.217034, 0.217006, 0.217217, 0.217255, 0.217293, 
  0.217414, 0.217505, 0.217599, 0.217807, 0.218016, 0.218223, 
  0.218427, 0.218635, 0.218841, 0.21905, 0.21957, 0.220094, 0.221148, 
  0.222212, 0.223288, 0.22547, 0.227696, 0.229966, 0.232282, 0.234646,
   0.237057, 0.239517, 0.242028, 0.244591, 0.247207, 0.249877, 
  0.252604, 0.259674, 0.267128, 0.274991, 0.283291, 0.292057, 
  0.301319, 0.31111, 0.321465, 0.332419, 0.344012, 0.356285, 0.369284,
   0.383057, 0.39766, 0.413152, 0.429604, 0.447094, 0.465712, 
  0.485563, 0.506772, 0.529486, 0.553883, 0.580183, 0.608655, 
  0.639647, 0.673609, 0.711146, 0.753105, 0.771386, 0.790674, 
  0.811095, 0.832808, 0.856013, 0.880965, 0.908002, 0.937583, 
  0.970359, 1.00729, 1.04987, 1.10064, 1.13044, 1.16442, 1.20427, 
  1.22731, 1.2533, 1.26476, 1.27699, 1.29015, 1.30444, 1.32016, 
  1.3378, 1.35815, 1.36977, 1.38276, 1.3977, 1.40621, 1.4157, 1.42658,
   1.42898, 1.43147, 1.43407, 1.43678, 1.43963, 1.44582, 1.4529, 
  1.45691, 1.46135}, {0.199814, 0.201353, 0.200132, 0.200956, 
  0.200467, 0.200722, 0.200829, 0.200716, 0.200604, 0.200626, 
  0.200705, 0.200826, 0.200822, 0.200914, 0.201004, 0.201103, 
  0.201273, 0.201448, 0.201625, 0.2018, 0.201975, 0.202153, 0.202329, 
  0.202773, 0.20322, 0.204121, 0.205032, 0.205953, 0.207826, 0.209742,
   0.211702, 0.213706, 0.215757, 0.217854, 0.220001, 0.222198, 
  0.224447, 0.226749, 0.229105, 0.231518, 0.237809, 0.244492, 
  0.251597, 0.25916, 0.267216, 0.275804, 0.284965, 0.294742, 0.305178,
   0.316322, 0.328221, 0.340928, 0.354499, 0.368992, 0.384475, 
  0.401019, 0.418708, 0.437635, 0.45791, 0.479661, 0.50304, 0.528232, 
  0.555461, 0.585008, 0.617226, 0.652577, 0.691676, 0.735379, 
  0.754411, 0.774477, 0.795706, 0.818253, 0.842317, 0.868149, 
  0.896082, 0.926573, 0.960272, 0.998162, 1.04182, 1.09394, 1.1246, 
  1.15963, 1.20071, 1.22441, 1.25107, 1.2628, 1.2753, 1.28872, 
  1.30326, 1.31923, 1.3371, 1.35766, 1.36938, 1.38248, 1.39751, 
  1.40606, 1.41559, 1.4265, 1.42891, 1.43141, 1.43401, 1.43673, 
  1.43959, 1.44579, 1.45288, 1.45689, 1.46134}};

  const double imomegaqnm33[8][107] = {{0.068612, 0.0735351, 0.0774984, 0.0829193, 0.0860369, 0.0870572, 
  0.0878377, 0.0884398, 0.0889085, 0.0892767, 0.0902685, 0.0906386, 
  0.0907974, 0.0909101, 0.0909414, 0.090952, 0.0909603, 0.0909664, 
  0.0909729, 0.0909798, 0.090987, 0.0909942, 0.0910015, 0.0910197, 
  0.0910378, 0.091074, 0.09111, 0.0911457, 0.0912166, 0.0912867, 
  0.091356, 0.0914243, 0.0914917, 0.0915582, 0.0916236, 0.091688, 
  0.0917514, 0.0918136, 0.0918746, 0.0919344, 0.0920784, 0.0922137, 
  0.0923394, 0.0924547, 0.0925583, 0.0926492, 0.0927258, 0.0927867, 
  0.0928302, 0.0928541, 0.0928562, 0.0928338, 0.092784, 0.092703, 
  0.0925869, 0.0924305, 0.0922281, 0.0919726, 0.0916556, 0.0912666, 
  0.0907928, 0.0902179, 0.0895213, 0.0886763, 0.087647, 0.0863849, 
  0.0848213, 0.0828557, 0.0819248, 0.0808922, 0.0797413, 0.0784512, 
  0.0769953, 0.075339, 0.0734361, 0.0712234, 0.0686106, 0.0654629, 
  0.0615646, 0.056537, 0.0533881, 0.0496087, 0.0449046, 0.0420439, 
  0.0386881, 0.0371631, 0.0355053, 0.0336873, 0.0316714, 0.0294027, 
  0.0267968, 0.023711, 0.0219107, 0.0198652, 0.0174725, 0.0160908, 
  0.014533, 0.0127266, 0.0123252, 0.011907, 0.0114701, 0.011012, 
  0.0105299, 0.00947732, 0.00826624, 0.0075764, 
  0.00680679}, {0.277496, 0.278402, 0.278809, 0.278955, 0.278824, 
  0.278742, 0.278664, 0.278594, 0.278534, 0.278482, 0.278334, 
  0.278289, 0.278281, 0.278294, 0.278307, 0.278318, 0.278336, 
  0.278354, 0.278373, 0.278391, 0.278409, 0.278427, 0.278446, 
  0.278491, 0.278537, 0.278627, 0.278716, 0.278805, 0.27898, 0.279153,
   0.279323, 0.279489, 0.279653, 0.279813, 0.27997, 0.280123, 
  0.280273, 0.280419, 0.280561, 0.280699, 0.281026, 0.281324, 
  0.281591, 0.281825, 0.28202, 0.282175, 0.282285, 0.282345, 0.28235, 
  0.282293, 0.282168, 0.281967, 0.281681, 0.281298, 0.280807, 
  0.280191, 0.279434, 0.278515, 0.277407, 0.27608, 0.274494, 0.272602,
   0.27034, 0.26763, 0.264363, 0.260394, 0.255518, 0.249435, 0.246568,
   0.243397, 0.239871, 0.235928, 0.231489, 0.226451, 0.220675, 
  0.213971, 0.206071, 0.19657, 0.184822, 0.169692, 0.160224, 0.148867,
   0.134739, 0.12615, 0.116076, 0.111499, 0.106524, 0.101069, 
  0.0950193, 0.088212, 0.0803933, 0.0711348, 0.0657335, 0.0595967, 
  0.0524182, 0.0482728, 0.0435993, 0.03818, 0.0369756, 0.0357211, 
  0.0344103, 0.0330362, 0.0315898, 0.028432, 0.0247987, 0.0227292, 
  0.0204204}, {0.486993, 0.484741, 0.483427, 0.482389, 0.482301, 
  0.482372, 0.48246, 0.482542, 0.482611, 0.482664, 0.482744, 0.482718,
   0.482702, 0.482709, 0.482722, 0.482729, 0.482743, 0.482759, 
  0.482774, 0.482789, 0.482804, 0.482819, 0.482834, 0.482871, 
  0.482907, 0.48298, 0.483051, 0.483121, 0.483256, 0.483386, 0.483511,
   0.48363, 0.483742, 0.483849, 0.483949, 0.484043, 0.48413, 0.48421, 
  0.484283, 0.484348, 0.484477, 0.484555, 0.484575, 0.484535, 
  0.484428, 0.484247, 0.483987, 0.483638, 0.483192, 0.48264, 0.481968,
   0.481164, 0.480211, 0.479093, 0.477786, 0.476267, 0.474506, 
  0.472466, 0.470105, 0.467369, 0.464195, 0.4605, 0.456181, 0.451103, 
  0.44509, 0.437899, 0.42919, 0.418467, 0.413455, 0.407937, 0.401829, 
  0.39503, 0.387407, 0.37879, 0.36895, 0.357573, 0.344212, 0.328195, 
  0.308447, 0.283082, 0.267235, 0.248247, 0.224647, 0.210309, 
  0.193501, 0.185866, 0.177568, 0.168469, 0.158382, 0.147032, 
  0.133998, 0.118564, 0.10956, 0.0993313, 0.0873659, 0.0804562, 
  0.0726665, 0.0636339, 0.0616265, 0.0595355, 0.0573509, 0.0550606, 
  0.0526499, 0.0473868, 0.0413313, 0.0378821, 0.034034}, {0.709974, 
  0.711601, 0.713499, 0.714385, 0.713785, 0.713546, 0.713418, 
  0.713383, 0.713411, 0.713468, 0.713655, 0.713595, 0.713563, 
  0.713584, 0.713571, 0.713561, 0.71355, 0.713535, 0.713522, 0.713507,
   0.713493, 0.713479, 0.713465, 0.713428, 0.713392, 0.713317, 
  0.713241, 0.713162, 0.712999, 0.712828, 0.712649, 0.712461, 
  0.712264, 0.712058, 0.711843, 0.711618, 0.711383, 0.711138, 
  0.710882, 0.710616, 0.709899, 0.709105, 0.708229, 0.707264, 
  0.706201, 0.705033, 0.703748, 0.702336, 0.700785, 0.699079, 
  0.697202, 0.695135, 0.692855, 0.690337, 0.68755, 0.684459, 0.68102, 
  0.677183, 0.672887, 0.668056, 0.662599, 0.656401, 0.649315, 
  0.641153, 0.631664, 0.620511, 0.607221, 0.591101, 0.583641, 
  0.575473, 0.566484, 0.556529, 0.545428, 0.532943, 0.518759, 
  0.502439, 0.48336, 0.460589, 0.432625, 0.396831, 0.37452, 0.347819, 
  0.314676, 0.294558, 0.270986, 0.260282, 0.248651, 0.235902, 
  0.221769, 0.20587, 0.187614, 0.166001, 0.153393, 0.13907, 0.122317, 
  0.112642, 0.101735, 0.0890889, 0.0862783, 0.0833507, 0.080292, 
  0.0770855, 0.0737103, 0.0663418, 0.057864, 0.053035, 
  0.0476477}, {0.972477, 0.969229, 0.96766, 0.969029, 0.969442, 
  0.969155, 0.968886, 0.968746, 0.968738, 0.968815, 0.969029, 
  0.968878, 0.968919, 0.968877, 0.968857, 0.968828, 0.968768, 
  0.968706, 0.968649, 0.968589, 0.968529, 0.96847, 0.96841, 0.96826, 
  0.96811, 0.967806, 0.9675, 0.96719, 0.966561, 0.96592, 0.965265, 
  0.964597, 0.963916, 0.96322, 0.962509, 0.961784, 0.961043, 0.960287,
   0.959514, 0.958725, 0.956678, 0.954518, 0.952237, 0.949827, 
  0.947277, 0.944577, 0.941715, 0.938677, 0.935446, 0.932005, 
  0.928331, 0.924401, 0.920185, 0.915649, 0.910755, 0.905456, 
  0.899695, 0.893406, 0.88651, 0.878908, 0.870479, 0.861073, 0.850501,
   0.838516, 0.824795, 0.808901, 0.790225, 0.767883, 0.757639, 
  0.746483, 0.73427, 0.720819, 0.7059, 0.689212, 0.670355, 0.648776, 
  0.623682, 0.593881, 0.557455, 0.511016, 0.482143, 0.44764, 0.404868,
   0.378929, 0.348555, 0.334768, 0.319791, 0.303378, 0.285189, 
  0.264731, 0.241246, 0.213447, 0.197233, 0.178815, 0.157271, 
  0.144831, 0.130807, 0.114545, 0.110932, 0.107167, 0.103234, 
  0.0991113, 0.0947715, 0.0852972, 0.0743968, 0.068188, 
  0.0612614}, {1.23488, 1.2385, 1.23757, 1.23638, 1.23735, 1.23719, 
  1.23689, 1.23671, 1.23672, 1.23684, 1.23689, 1.23686, 1.23687, 
  1.23681, 1.23674, 1.23671, 1.2366, 1.23651, 1.23641, 1.23631, 
  1.23622, 1.23612, 1.23602, 1.23578, 1.23553, 1.23504, 1.23454, 
  1.23403, 1.23301, 1.23196, 1.23088, 1.22979, 1.22867, 1.22753, 
  1.22637, 1.22518, 1.22396, 1.22272, 1.22145, 1.22016, 1.21681, 
  1.21328, 1.20956, 1.20564, 1.20152, 1.19717, 1.19258, 1.18775, 
  1.18264, 1.17724, 1.17153, 1.16546, 1.15902, 1.15215, 1.14481, 
  1.13695, 1.12848, 1.11934, 1.10942, 1.09859, 1.08672, 1.07361, 
  1.05902, 1.04265, 1.02409, 1.00282, 0.978074, 0.948772, 0.935434, 
  0.920971, 0.905209, 0.887929, 0.868853, 0.847622, 0.823753, 
  0.796583, 0.765161, 0.728048, 0.682917, 0.625639, 0.590118, 
  0.547732, 0.49525, 0.463448, 0.42623, 0.409344, 0.391005, 0.370913, 
  0.348652, 0.323623, 0.294897, 0.260904, 0.241081, 0.218564, 
  0.192229, 0.177023, 0.159881, 0.140004, 0.135587, 0.130985, 
  0.126178, 0.121138, 0.115834, 0.104253, 0.09093, 0.0833412, 
  0.0748751}, {1.50812, 1.50714, 1.50886, 1.50765, 1.50833, 1.50831, 
  1.50805, 1.5079, 1.50794, 1.50808, 1.50796, 1.50805, 1.50794, 
  1.50791, 1.50787, 1.50779, 1.50768, 1.50756, 1.50745, 1.50733, 
  1.50721, 1.50709, 1.50697, 1.50667, 1.50637, 1.50577, 1.50515, 
  1.50453, 1.50326, 1.50195, 1.50061, 1.49924, 1.49784, 1.49639, 
  1.49491, 1.4934, 1.49185, 1.49025, 1.48862, 1.48695, 1.4826, 
  1.47798, 1.47308, 1.46789, 1.46239, 1.45657, 1.45042, 1.44391, 
  1.43702, 1.42973, 1.42201, 1.41383, 1.40515, 1.39591, 1.38607, 
  1.37555, 1.36428, 1.35215, 1.33904, 1.32481, 1.30928, 1.29221, 
  1.27333, 1.25226, 1.22852, 1.20145, 1.17016, 1.13335, 1.11668, 
  1.09865, 1.07906, 1.05766, 1.03413, 1.00803, 0.978819, 0.945725, 
  0.907651, 0.862934, 0.808869, 0.7406, 0.698384, 0.648072, 0.585828, 
  0.54813, 0.50403, 0.484027, 0.462309, 0.43852, 0.412171, 0.382553, 
  0.348572, 0.308374, 0.284937, 0.258319, 0.227191, 0.209218, 
  0.188957, 0.165464, 0.160243, 0.154805, 0.149123, 0.143167, 
  0.136898, 0.123211, 0.107464, 0.0984948, 0.102103}, {1.77947, 
  1.77829, 1.77828, 1.77845, 1.7786, 1.77868, 1.7785, 1.77838, 
  1.77844, 1.77854, 1.77842, 1.77841, 1.77841, 1.77834, 1.77827, 
  1.7782, 1.77808, 1.77795, 1.77783, 1.7777, 1.77757, 1.77744, 
  1.77732, 1.77699, 1.77667, 1.77601, 1.77534, 1.77466, 1.77327, 
  1.77183, 1.77035, 1.76883, 1.76726, 1.76564, 1.76398, 1.76226, 
  1.7605, 1.75868, 1.75682, 1.75489, 1.74985, 1.74445, 1.73867, 
  1.73249, 1.72589, 1.71886, 1.71137, 1.7034, 1.69493, 1.68592, 
  1.67635, 1.66618, 1.65536, 1.64384, 1.63157, 1.61846, 1.60441, 
  1.58931, 1.57303, 1.55539, 1.53617, 1.5151, 1.49185, 1.46597, 
  1.43689, 1.40384, 1.36575, 1.32111, 1.30093, 1.27916, 1.25555, 
  1.22981, 1.20156, 1.17032, 1.13545, 1.09609, 1.05101, 0.998352, 
  0.935078, 0.85569, 0.806777, 0.748562, 0.676569, 0.632966, 0.58196, 
  0.558828, 0.533715, 0.506213, 0.475757, 0.441533, 0.40228, 0.355861,
   0.328804, 0.29808, 0.262156, 0.241415, 0.218035, 0.190926, 
  0.184902, 0.178626, 0.17207, 0.165197, 0.157962, 0.142169, 0.123998,
   0.113649, 0.115718}};


  /* 4, 4 mode */
  const double reomegaqnm44[8][107] = {{0.603485, 0.613847, 0.619636, 0.623952, 0.624219, 0.623894, 
  0.623504, 0.623122, 0.62278, 0.622487, 0.621648, 0.621375, 0.621309,
   0.621365, 0.621483, 0.621613, 0.621877, 0.622141, 0.622404, 
  0.622667, 0.62293, 0.623194, 0.623458, 0.624119, 0.624781, 0.626113,
   0.627452, 0.628799, 0.631518, 0.634269, 0.637054, 0.639872, 
  0.642726, 0.645615, 0.648541, 0.651503, 0.654504, 0.657544, 
  0.660623, 0.663743, 0.671728, 0.679989, 0.688543, 0.697411, 
  0.706611, 0.716168, 0.726107, 0.736455, 0.747243, 0.758508, 
  0.770286, 0.782624, 0.795569, 0.809178, 0.823517, 0.83866, 0.854693,
   0.871718, 0.889853, 0.909242, 0.930054, 0.9525, 0.976839, 1.0034, 
  1.03259, 1.06498, 1.10131, 1.14265, 1.16092, 1.18036, 1.20114, 
  1.22345, 1.24755, 1.27374, 1.30245, 1.33422, 1.36984, 1.41042, 
  1.45773, 1.51478, 1.54862, 1.58759, 1.6339, 1.66102, 1.69194, 
  1.7057, 1.72045, 1.73641, 1.75384, 1.77314, 1.79492, 1.82022, 
  1.83474, 1.85104, 1.86985, 1.8806, 1.8926, 1.9064, 1.90945, 1.91261,
   1.91592, 1.91937, 1.92299, 1.93088, 1.93989, 1.945, 
  1.95068}, {0.59231, 0.596523, 0.598412, 0.599977, 0.600617, 
  0.600806, 0.600947, 0.601053, 0.601136, 0.601202, 0.601389, 
  0.601485, 0.601558, 0.601693, 0.601828, 0.601964, 0.602238, 
  0.602511, 0.602785, 0.60306, 0.603334, 0.603609, 0.603884, 0.604574,
   0.605265, 0.606653, 0.60805, 0.609455, 0.61229, 0.615158, 0.618061,
   0.620998, 0.623972, 0.626982, 0.630029, 0.633115, 0.63624, 
  0.639405, 0.64261, 0.645858, 0.654166, 0.662758, 0.671651, 0.680865,
   0.69042, 0.700341, 0.710653, 0.721383, 0.732565, 0.744232, 
  0.756425, 0.769188, 0.782571, 0.796632, 0.811434, 0.827055, 
  0.843581, 0.861114, 0.879773, 0.899703, 0.921074, 0.944098, 
  0.969034, 0.996209, 1.02604, 1.05909, 1.09611, 1.13816, 1.15671, 
  1.17644, 1.19752, 1.22013, 1.24453, 1.27103, 1.30005, 1.33214, 
  1.36806, 1.40896, 1.45658, 1.51396, 1.54795, 1.58707, 1.63353, 
  1.66073, 1.69172, 1.7055, 1.72028, 1.73626, 1.75372, 1.77304, 
  1.79485, 1.82017, 1.8347, 1.85101, 1.86983, 1.88058, 1.89259, 
  1.90639, 1.90944, 1.91261, 1.91591, 1.91936, 1.92299, 1.93087, 
  1.93989, 1.945, 1.95068}, {0.559386, 0.561881, 0.562985, 0.563673, 
  0.563747, 0.563732, 0.563712, 0.563695, 0.563683, 0.563678, 
  0.563718, 0.563794, 0.563872, 0.564021, 0.564168, 0.564316, 
  0.564612, 0.564908, 0.565204, 0.565501, 0.565798, 0.566096, 
  0.566394, 0.567139, 0.567887, 0.56939, 0.570901, 0.57242, 0.575484, 
  0.578584, 0.58172, 0.584892, 0.588102, 0.59135, 0.594637, 0.597965, 
  0.601333, 0.604743, 0.608195, 0.611692, 0.62063, 0.629865, 0.639415,
   0.649299, 0.65954, 0.670161, 0.681189, 0.692653, 0.704586, 
  0.717023, 0.730006, 0.743579, 0.757794, 0.77271, 0.788392, 0.804918,
   0.822375, 0.840868, 0.860518, 0.88147, 0.903897, 0.928012, 
  0.954076, 0.982419, 1.01346, 1.04777, 1.08608, 1.12947, 1.14857, 
  1.16887, 1.19051, 1.21369, 1.23867, 1.26576, 1.29538, 1.32806, 
  1.36459, 1.4061, 1.45434, 1.51232, 1.54662, 1.58604, 1.6328, 
  1.66014, 1.69127, 1.70511, 1.71994, 1.73598, 1.75348, 1.77286, 
  1.79471, 1.82008, 1.83463, 1.85096, 1.8698, 1.88055, 1.89257, 
  1.90638, 1.90943, 1.9126, 1.9159, 1.91935, 1.92298, 1.93087, 
  1.93989, 1.945, 1.95068}, {0.514416, 0.513669, 0.512866, 0.5122, 
  0.512256, 0.512336, 0.512404, 0.512455, 0.512488, 0.51251, 0.512565,
   0.51264, 0.512726, 0.512889, 0.513051, 0.513214, 0.51354, 0.513866,
   0.514193, 0.51452, 0.514848, 0.515175, 0.515503, 0.516325, 
  0.517149, 0.518804, 0.520467, 0.522139, 0.525511, 0.52892, 0.532367,
   0.535852, 0.539377, 0.542943, 0.546549, 0.550198, 0.553889, 
  0.557625, 0.561406, 0.565233, 0.575009, 0.585097, 0.595517, 
  0.606289, 0.617438, 0.628986, 0.640963, 0.653399, 0.666327, 
  0.679785, 0.693815, 0.708464, 0.723784, 0.739837, 0.75669, 0.774421,
   0.793122, 0.812899, 0.833875, 0.856197, 0.880043, 0.905628, 
  0.933216, 0.96314, 0.995828, 1.03184, 1.07193, 1.11717, 1.13703, 
  1.1581, 1.18052, 1.20451, 1.2303, 1.25821, 1.28866, 1.32219, 
  1.35958, 1.40195, 1.45106, 1.50993, 1.54467, 1.58453, 1.63172, 
  1.65927, 1.69061, 1.70453, 1.71944, 1.73555, 1.75313, 1.77258, 
  1.79451, 1.81994, 1.83452, 1.85088, 1.86975, 1.88051, 1.89254, 
  1.90636, 1.90941, 1.91258, 1.91589, 1.91934, 1.92297, 1.93086, 
  1.93988, 1.94499, 1.95068}, {0.454601, 0.454746, 0.455842, 0.456343,
   0.456046, 0.455998, 0.456018, 0.456071, 0.456127, 0.456171, 
  0.456246, 0.456327, 0.456423, 0.456596, 0.456773, 0.456949, 
  0.457302, 0.457655, 0.458009, 0.458363, 0.458717, 0.459072, 
  0.459427, 0.460317, 0.461208, 0.462999, 0.464799, 0.466609, 
  0.470257, 0.473946, 0.477675, 0.481446, 0.485259, 0.489115, 
  0.493016, 0.496962, 0.500954, 0.504993, 0.509081, 0.513219, 
  0.523786, 0.534689, 0.545948, 0.557585, 0.569625, 0.582093, 0.59502,
   0.608436, 0.622378, 0.636885, 0.652, 0.667772, 0.684256, 0.701516, 
  0.71962, 0.738651, 0.758702, 0.779882, 0.802317, 0.82616, 0.851591, 
  0.878829, 0.908143, 0.939873, 0.974449, 1.01244, 1.05461, 1.10202, 
  1.12279, 1.14477, 1.16814, 1.19308, 1.21984, 1.24875, 1.28021, 
  1.31477, 1.35322, 1.39666, 1.44687, 1.50684, 1.54214, 1.58256, 
  1.63031, 1.65813, 1.68974, 1.70376, 1.71878, 1.73499, 1.75267, 
  1.77222, 1.79424, 1.81975, 1.83438, 1.85077, 1.86967, 1.88046, 
  1.8925, 1.90633, 1.90939, 1.91256, 1.91587, 1.91932, 1.92296, 
  1.93085, 1.93987, 1.94499, 1.95067}, {0.405546, 0.404573, 0.403366, 
  0.404025, 0.404115, 0.403989, 0.40394, 0.403971, 0.404037, 0.404098,
   0.40416, 0.404265, 0.404355, 0.404535, 0.404715, 0.404897, 0.40526,
   0.405624, 0.405988, 0.406353, 0.406718, 0.407083, 0.407449, 
  0.408365, 0.409284, 0.411129, 0.412986, 0.414853, 0.418619, 0.42243,
   0.426286, 0.430188, 0.434137, 0.438134, 0.44218, 0.446277, 
  0.450425, 0.454625, 0.458878, 0.463187, 0.474206, 0.485594, 
  0.497376, 0.509573, 0.522212, 0.535321, 0.54893, 0.563075, 0.577791,
   0.593121, 0.609109, 0.625808, 0.643274, 0.661572, 0.680777, 
  0.700971, 0.722252, 0.744732, 0.768542, 0.793839, 0.820807, 
  0.849672, 0.880708, 0.914263, 0.950774, 0.990817, 1.03516, 1.0849, 
  1.10664, 1.12961, 1.15399, 1.17998, 1.20782, 1.23782, 1.27041, 
  1.30613, 1.34576, 1.39042, 1.44188, 1.50314, 1.53909, 1.58017, 
  1.62858, 1.65673, 1.68866, 1.70281, 1.71796, 1.7343, 1.7521, 
  1.77176, 1.7939, 1.81952, 1.83419, 1.85064, 1.86958, 1.88039, 
  1.89245, 1.9063, 1.90936, 1.91253, 1.91584, 1.9193, 1.92294, 
  1.93083, 1.93986, 1.94498, 1.95067}, {0.35976, 0.362048, 0.362088, 
  0.361405, 0.361838, 0.361724, 0.361635, 0.361653, 0.361731, 
  0.361798, 0.361829, 0.361951, 0.362022, 0.362204, 0.362384, 
  0.362559, 0.362914, 0.363269, 0.363624, 0.36398, 0.364336, 0.364693,
   0.36505, 0.365945, 0.366843, 0.368648, 0.370465, 0.372293, 
  0.375986, 0.379729, 0.383522, 0.387367, 0.391264, 0.395214, 0.39922,
   0.403281, 0.4074, 0.411576, 0.415813, 0.42011, 0.431129, 0.44256, 
  0.454427, 0.466756, 0.479575, 0.492914, 0.506807, 0.521289, 
  0.536401, 0.552185, 0.568691, 0.585972, 0.604089, 0.623109, 
  0.643109, 0.664176, 0.686412, 0.709933, 0.734875, 0.761399, 
  0.789696, 0.819999, 0.85259, 0.887824, 0.92615, 0.968153, 1.01462, 
  1.06665, 1.08936, 1.11334, 1.13875, 1.16581, 1.19475, 1.22589, 
  1.25966, 1.29659, 1.33748, 1.38345, 1.43625, 1.49892, 1.53559, 
  1.57741, 1.62657, 1.6551, 1.6874, 1.7017, 1.717, 1.73348, 1.75142, 
  1.77123, 1.79349, 1.81924, 1.83398, 1.85048, 1.86948, 1.8803, 
  1.89239, 1.90626, 1.90932, 1.9125, 1.91581, 1.91928, 1.92291, 
  1.93082, 1.93985, 1.94497, 1.95066},
{0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72}};

  const double imomegaqnm44[8][107] = {{0.0472887, 0.0570846, 0.0652857, 0.0769202, 0.0835313, 0.0855743, 
  0.0870596, 0.0881437, 0.0889418, 0.0895359, 0.0909619, 0.0914196, 
  0.0916056, 0.091741, 0.0917845, 0.0918026, 0.0918186, 0.0918281, 
  0.0918366, 0.0918448, 0.0918531, 0.0918614, 0.0918696, 0.0918903, 
  0.0919109, 0.0919521, 0.0919931, 0.092034, 0.0921151, 0.0921956, 
  0.0922753, 0.0923543, 0.0924324, 0.0925097, 0.0925861, 0.0926616, 
  0.0927361, 0.0928096, 0.0928821, 0.0929535, 0.0931268, 0.093292, 
  0.0934484, 0.0935949, 0.0937302, 0.0938532, 0.0939623, 0.0940559, 
  0.0941321, 0.0941887, 0.0942233, 0.094233, 0.0942145, 0.094164, 
  0.0940768, 0.0939478, 0.0937705, 0.0935374, 0.0932393, 0.092865, 
  0.0924006, 0.0918287, 0.0911274, 0.0902679, 0.0892124, 0.0879094, 
  0.086287, 0.0842401, 0.0832692, 0.0821917, 0.0809905, 0.0796442, 
  0.0781255, 0.0763992, 0.0744185, 0.0721192, 0.0694103, 0.066156, 
  0.0621396, 0.056982, 0.0537633, 0.0499111, 0.0451315, 0.0422323, 
  0.0388376, 0.037297, 0.0356236, 0.0337901, 0.0317587, 0.0294747, 
  0.0268536, 0.0237529, 0.0219453, 0.0198926, 0.0174928, 0.0161075, 
  0.0145463, 0.0127365, 0.0123344, 0.0119155, 0.0114779, 0.0110192, 
  0.0105364, 0.00948248, 0.00827008, 0.00757959, 
  0.00680933}, {0.267496, 0.27243, 0.275117, 0.277412, 0.278166, 
  0.278332, 0.278431, 0.278489, 0.278523, 0.278542, 0.278554, 
  0.278543, 0.278539, 0.278544, 0.278555, 0.278567, 0.27859, 0.278613,
   0.278636, 0.278659, 0.278682, 0.278704, 0.278727, 0.278784, 
  0.278841, 0.278955, 0.279068, 0.27918, 0.279403, 0.279623, 0.279841,
   0.280057, 0.28027, 0.28048, 0.280687, 0.280892, 0.281093, 0.281291,
   0.281485, 0.281676, 0.282137, 0.282572, 0.282979, 0.283355, 
  0.283696, 0.283998, 0.284257, 0.284467, 0.284624, 0.284721, 
  0.284749, 0.284701, 0.284567, 0.284334, 0.28399, 0.283519, 0.282901,
   0.282113, 0.281129, 0.279913, 0.278426, 0.276614, 0.274412, 
  0.271733, 0.268465, 0.264452, 0.25948, 0.253233, 0.250278, 0.247003,
   0.243357, 0.239276, 0.234678, 0.229458, 0.223476, 0.216539, 
  0.208375, 0.198576, 0.186494, 0.170992, 0.161324, 0.149756, 
  0.135408, 0.126707, 0.116519, 0.111896, 0.106875, 0.101374, 
  0.0952787, 0.0884259, 0.0805622, 0.0712595, 0.0658364, 0.0596782, 
  0.0524786, 0.0483228, 0.043639, 0.0382095, 0.0370031, 0.0357465, 
  0.0344338, 0.0330576, 0.0316092, 0.0284475, 0.0248103, 0.0227388, 
  0.020428}, {0.476992, 0.476552, 0.476059, 0.475354, 0.475023, 
  0.474945, 0.474902, 0.47488, 0.474873, 0.474873, 0.474904, 0.474923,
   0.474932, 0.474945, 0.47496, 0.474975, 0.475006, 0.475036, 
  0.475066, 0.475097, 0.475127, 0.475157, 0.475187, 0.475263, 
  0.475338, 0.475487, 0.475635, 0.475782, 0.476073, 0.476359, 
  0.476641, 0.476918, 0.47719, 0.477456, 0.477718, 0.477973, 0.478223,
   0.478467, 0.478704, 0.478936, 0.479483, 0.479984, 0.480434, 
  0.480828, 0.48116, 0.481423, 0.48161, 0.481712, 0.481721, 0.481624, 
  0.48141, 0.481064, 0.48057, 0.479908, 0.479056, 0.477986, 0.476666, 
  0.47506, 0.47312, 0.470791, 0.468003, 0.464669, 0.460679, 0.45589, 
  0.450113, 0.44309, 0.434465, 0.423712, 0.41865, 0.413055, 0.406842, 
  0.399905, 0.392108, 0.383276, 0.373175, 0.361487, 0.347759, 
  0.331313, 0.311071, 0.28514, 0.268985, 0.249669, 0.225726, 0.21121, 
  0.194221, 0.186512, 0.17814, 0.168968, 0.158807, 0.147383, 0.134275,
   0.118768, 0.109729, 0.099465, 0.0874652, 0.0805385, 0.072732, 
  0.0636828, 0.0616721, 0.0595777, 0.0573898, 0.0550962, 0.0526822, 
  0.0474125, 0.0413505, 0.037898, 0.0340467}, {0.688887, 0.687276, 
  0.68687, 0.687194, 0.687484, 0.687529, 0.687535, 0.687523, 0.687504,
   0.687487, 0.687464, 0.68748, 0.687488, 0.687497, 0.687509, 
  0.687521, 0.687544, 0.687567, 0.68759, 0.687613, 0.687636, 0.687658,
   0.687681, 0.687738, 0.687794, 0.687905, 0.688015, 0.688122, 
  0.688333, 0.688537, 0.688734, 0.688923, 0.689105, 0.689278, 
  0.689444, 0.689601, 0.689749, 0.689889, 0.690018, 0.690138, 
  0.690394, 0.690579, 0.69069, 0.690717, 0.690652, 0.690487, 0.690211,
   0.689812, 0.689276, 0.688589, 0.687732, 0.686685, 0.685425, 
  0.683924, 0.682151, 0.680068, 0.677629, 0.674782, 0.671462, 
  0.667589, 0.663066, 0.657772, 0.651551, 0.644203, 0.635463, 
  0.624972, 0.612229, 0.596503, 0.589146, 0.581044, 0.572077, 
  0.562097, 0.550915, 0.538287, 0.523888, 0.507273, 0.48781, 0.464555,
   0.436002, 0.399509, 0.376808, 0.349693, 0.316109, 0.295761, 
  0.271954, 0.261154, 0.249426, 0.236578, 0.222347, 0.206349, 
  0.187993, 0.166281, 0.153625, 0.139254, 0.122453, 0.112755, 
  0.101826, 0.0891563, 0.0863413, 0.0834091, 0.080346, 0.0771349, 
  0.0737552, 0.0663776, 0.0578907, 0.0530572, 0.0476654}, {0.918235, 
  0.920515, 0.920895, 0.920137, 0.920012, 0.920094, 0.920167, 
  0.920207, 0.920214, 0.920202, 0.920149, 0.920166, 0.920164, 0.92016,
   0.920158, 0.920156, 0.920152, 0.920148, 0.920144, 0.92014, 
  0.920136, 0.920132, 0.920127, 0.920116, 0.920104, 0.920079, 
  0.920051, 0.920022, 0.919956, 0.919881, 0.919796, 0.919702, 
  0.919597, 0.919482, 0.919356, 0.919218, 0.919069, 0.918908, 
  0.918734, 0.918547, 0.91802, 0.917401, 0.916681, 0.915851, 0.9149, 
  0.913815, 0.912584, 0.91119, 0.909617, 0.907843, 0.905847, 0.903603,
   0.901078, 0.898239, 0.895043, 0.891442, 0.887378, 0.88278, 
  0.877566, 0.871633, 0.864855, 0.857075, 0.848091, 0.837647, 0.8254, 
  0.81089, 0.79347, 0.772205, 0.762324, 0.751485, 0.739534, 0.72628, 
  0.711484, 0.694831, 0.675907, 0.654142, 0.628729, 0.598457, 0.5614, 
  0.514171, 0.484847, 0.449863, 0.406581, 0.380377, 0.349729, 
  0.335831, 0.32074, 0.304211, 0.285904, 0.265327, 0.24172, 0.213799, 
  0.197524, 0.179045, 0.157442, 0.144973, 0.13092, 0.11463, 0.111011, 
  0.107241, 0.103302, 0.0991739, 0.0948285, 0.0853429, 0.0744311, 
  0.0682165, 0.0612841}, {1.17243, 1.17025, 1.17072, 1.17153, 1.17109,
   1.17107, 1.17115, 1.17122, 1.17125, 1.17124, 1.17116, 1.17118, 
  1.17115, 1.17114, 1.17112, 1.17109, 1.17105, 1.171, 1.17096, 
  1.17092, 1.17087, 1.17083, 1.17078, 1.17067, 1.17055, 1.17033, 
  1.17009, 1.16986, 1.16938, 1.16889, 1.16838, 1.16786, 1.16733, 
  1.16679, 1.16622, 1.16565, 1.16506, 1.16445, 1.16382, 1.16318, 
  1.1615, 1.15969, 1.15776, 1.15567, 1.15344, 1.15103, 1.14843, 
  1.14562, 1.14259, 1.1393, 1.13573, 1.13184, 1.12761, 1.12298, 
  1.1179, 1.11232, 1.10615, 1.09933, 1.09174, 1.08325, 1.07371, 
  1.06294, 1.05067, 1.0366, 1.0203, 1.00121, 0.97854, 0.951144, 
  0.9385, 0.92468, 0.909499, 0.892724, 0.874063, 0.853135, 0.829435, 
  0.802272, 0.770663, 0.733138, 0.687353, 0.629184, 0.593147, 
  0.550213, 0.497162, 0.465071, 0.427558, 0.410551, 0.392089, 
  0.371872, 0.349482, 0.32432, 0.295457, 0.261323, 0.241428, 0.218839,
   0.192433, 0.177192, 0.160015, 0.140105, 0.135681, 0.131073, 
  0.126259, 0.121213, 0.115902, 0.104308, 0.0909715, 0.0833759, 
  0.0749029}, {1.43411, 1.43505, 1.43374, 1.43433, 1.43411, 1.43401, 
  1.43406, 1.43415, 1.43418, 1.43415, 1.4341, 1.43407, 1.43405, 
  1.43401, 1.43397, 1.43393, 1.43384, 1.43376, 1.43367, 1.43359, 
  1.4335, 1.43342, 1.43333, 1.43312, 1.4329, 1.43247, 1.43203, 
  1.43159, 1.43069, 1.42978, 1.42884, 1.42789, 1.42691, 1.42592, 
  1.4249, 1.42387, 1.42281, 1.42173, 1.42062, 1.41949, 1.41656, 
  1.41345, 1.41017, 1.4067, 1.40302, 1.39911, 1.39496, 1.39054, 
  1.38582, 1.38079, 1.3754, 1.36961, 1.36339, 1.35669, 1.34943, 
  1.34155, 1.33296, 1.32356, 1.31324, 1.30183, 1.28915, 1.27497, 
  1.259, 1.24085, 1.22004, 1.19589, 1.16745, 1.13339, 1.11775, 
  1.10072, 1.08207, 1.06154, 1.03877, 1.01331, 0.984584, 0.951766, 
  0.913705, 0.868675, 0.81392, 0.744589, 0.701738, 0.650764, 0.587868,
   0.549856, 0.505448, 0.485323, 0.46348, 0.439564, 0.413085, 
  0.383331, 0.349205, 0.308853, 0.285337, 0.258637, 0.227427, 
  0.209413, 0.189112, 0.16558, 0.160352, 0.154906, 0.149217, 0.143253,
   0.136976, 0.123274, 0.107512, 0.0985354, 0.0885218}, 
{0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28}};

  /* 5, 5 mode */
  const double reomegaqnm55[8][107] = {{0.752532, 0.769133, 0.779179, 0.78709, 0.786149, 0.784329, 0.782338,
   0.780485, 0.778903, 0.777615, 0.774286, 0.773282, 0.772966, 
  0.77289, 0.772999, 0.773149, 0.773476, 0.77381, 0.774145, 0.77448, 
  0.774816, 0.775152, 0.775488, 0.776331, 0.777176, 0.778874, 
  0.780582, 0.782299, 0.785766, 0.789274, 0.792825, 0.79642, 0.800059,
   0.803743, 0.807474, 0.811253, 0.81508, 0.818956, 0.822883, 
  0.826863, 0.837046, 0.847582, 0.858493, 0.869802, 0.881536, 
  0.893724, 0.906398, 0.919594, 0.93335, 0.947712, 0.962728, 0.978454,
   0.994953, 1.0123, 1.03056, 1.04985, 1.07027, 1.09194, 1.11502, 
  1.13968, 1.16614, 1.19467, 1.22558, 1.25928, 1.29631, 1.33734, 
  1.38332, 1.43555, 1.45861, 1.48313, 1.50932, 1.53742, 1.56774, 
  1.60067, 1.63671, 1.67655, 1.72115, 1.77188, 1.83092, 1.90197, 
  1.94403, 1.99239, 2.04977, 2.08333, 2.12154, 2.13852, 2.15673, 
  2.17641, 2.1979, 2.22168, 2.24849, 2.27961, 2.29746, 2.31749, 
  2.34058, 2.35377, 2.3685, 2.38542, 2.38915, 2.39303, 2.39708, 
  2.40131, 2.40575, 2.41541, 2.42645, 2.43271, 2.43967}, {0.753001, 
  0.757208, 0.757533, 0.756804, 0.756468, 0.756414, 0.756397, 
  0.756401, 0.756417, 0.756438, 0.756561, 0.756667, 0.756761, 
  0.756937, 0.757109, 0.75728, 0.757623, 0.757967, 0.758311, 0.758656,
   0.759001, 0.759346, 0.759692, 0.760558, 0.761427, 0.763172, 
  0.764927, 0.766693, 0.770255, 0.77386, 0.777508, 0.7812, 0.784938, 
  0.788722, 0.792553, 0.796433, 0.800362, 0.804341, 0.808372, 
  0.812456, 0.822904, 0.833711, 0.844898, 0.85649, 0.868513, 0.880997,
   0.893974, 0.90748, 0.921555, 0.936243, 0.951594, 0.967664, 
  0.984515, 1.00222, 1.02086, 1.04054, 1.06135, 1.08343, 1.10693, 
  1.13203, 1.15894, 1.18793, 1.21932, 1.25352, 1.29105, 1.33262, 
  1.37915, 1.43196, 1.45525, 1.48001, 1.50644, 1.53478, 1.56534, 
  1.59852, 1.63481, 1.67491, 1.71976, 1.77074, 1.83004, 1.90134, 
  1.94352, 1.992, 2.0495, 2.08311, 2.12137, 2.13838, 2.15661, 2.17631,
   2.19781, 2.22161, 2.24844, 2.27957, 2.29743, 2.31747, 2.34057, 
  2.35376, 2.36849, 2.38541, 2.38915, 2.39303, 2.39708, 2.40131, 
  2.40575, 2.41541, 2.42645, 2.43271, 2.43967}, {0.7196, 0.722567, 
  0.723924, 0.725024, 0.725386, 0.725467, 0.725518, 0.725552, 
  0.725576, 0.725595, 0.725675, 0.72576, 0.725849, 0.72603, 0.726211, 
  0.726392, 0.726754, 0.727117, 0.72748, 0.727843, 0.728207, 0.728571,
   0.728936, 0.729849, 0.730765, 0.732605, 0.734456, 0.736316, 
  0.740071, 0.743868, 0.747711, 0.751599, 0.755534, 0.759516, 
  0.763547, 0.767627, 0.771759, 0.775942, 0.780178, 0.784469, 
  0.795442, 0.806783, 0.818514, 0.830662, 0.843253, 0.856316, 
  0.869886, 0.883998, 0.898692, 0.914015, 0.930016, 0.946752, 
  0.964288, 0.982696, 1.00206, 1.02247, 1.04405, 1.06692, 1.09123, 
  1.11716, 1.14494, 1.17482, 1.20713, 1.24229, 1.28082, 1.32342, 
  1.37101, 1.42493, 1.44867, 1.47389, 1.50078, 1.5296, 1.56064, 
  1.5943, 1.63108, 1.67167, 1.71702, 1.7685, 1.82829, 1.90008, 1.9425,
   1.99122, 2.04895, 2.08267, 2.12104, 2.13809, 2.15636, 2.1761, 
  2.19764, 2.22147, 2.24834, 2.27951, 2.29738, 2.31743, 2.34054, 
  2.35374, 2.36848, 2.3854, 2.38914, 2.39302, 2.39707, 2.4013, 
  2.40575, 2.41541, 2.42645, 2.43271, 2.43966}, {0.680361, 0.681655, 
  0.681927, 0.681788, 0.681669, 0.681657, 0.681664, 0.681682, 
  0.681705, 0.681729, 0.681838, 0.681934, 0.682031, 0.682225, 0.68242,
   0.682614, 0.683004, 0.683394, 0.683784, 0.684175, 0.684566, 
  0.684957, 0.685349, 0.686331, 0.687316, 0.689293, 0.69128, 0.693279,
   0.69731, 0.701387, 0.705509, 0.709679, 0.713897, 0.718164, 
  0.722481, 0.72685, 0.731271, 0.735746, 0.740276, 0.744863, 0.756582,
   0.768683, 0.781187, 0.794122, 0.807514, 0.821395, 0.835799, 
  0.850762, 0.866327, 0.88254, 0.899452, 0.917122, 0.935614, 0.955004,
   0.975376, 0.996827, 1.01947, 1.04344, 1.06888, 1.09599, 1.12497, 
  1.1561, 1.18971, 1.22621, 1.26614, 1.31018, 1.35928, 1.41476, 
  1.43915, 1.46502, 1.49258, 1.52208, 1.5538, 1.58816, 1.62564, 
  1.66694, 1.713, 1.76521, 1.82572, 1.89823, 1.941, 1.99007, 2.04813, 
  2.08202, 2.12055, 2.13766, 2.15599, 2.17578, 2.19738, 2.22127, 
  2.24819, 2.2794, 2.2973, 2.31737, 2.34051, 2.35371, 2.36846, 
  2.38539, 2.38913, 2.39301, 2.39706, 2.40129, 2.40574, 2.4154, 
  2.42645, 2.43271, 2.43966}, {0.630475, 0.629154, 0.628652, 0.628812,
   0.628995, 0.62902, 0.629028, 0.629035, 0.629045, 0.629061, 
  0.629172, 0.629277, 0.629381, 0.629592, 0.629802, 0.630012, 
  0.630433, 0.630855, 0.631277, 0.631699, 0.632122, 0.632545, 
  0.632968, 0.634029, 0.635093, 0.637229, 0.639376, 0.641534, 
  0.645886, 0.650284, 0.65473, 0.659225, 0.66377, 0.668366, 0.673014, 
  0.677716, 0.682472, 0.687284, 0.692153, 0.69708, 0.709663, 0.722641,
   0.73604, 0.749885, 0.764206, 0.779035, 0.794407, 0.81036, 0.826938,
   0.844188, 0.862164, 0.880925, 0.900539, 0.921082, 0.94264, 
  0.965314, 0.989218, 1.01449, 1.04128, 1.06978, 1.1002, 1.13284, 
  1.168, 1.20613, 1.24774, 1.29355, 1.34449, 1.4019, 1.42708, 1.45377,
   1.48215, 1.51248, 1.54506, 1.58029, 1.61866, 1.66085, 1.70782, 
  1.76095, 1.82238, 1.89581, 1.93904, 1.98856, 2.04706, 2.08116, 
  2.1199, 2.13709, 2.1555, 2.17537, 2.19704, 2.221, 2.24799, 2.27927, 
  2.29719, 2.31729, 2.34045, 2.35367, 2.36843, 2.38537, 2.38911, 
  2.39299, 2.39704, 2.40128, 2.40573, 2.41539, 2.42644, 2.4327, 
  2.43966}, {0.571708, 0.573132, 0.573635, 0.573297, 0.573329, 
  0.573398, 0.573443, 0.573465, 0.573477, 0.573488, 0.5736, 0.573712, 
  0.573823, 0.574048, 0.574272, 0.574496, 0.574944, 0.575393, 
  0.575843, 0.576292, 0.576743, 0.577194, 0.577645, 0.578775, 
  0.579908, 0.582183, 0.58447, 0.586768, 0.591402, 0.596085, 0.600818,
   0.605603, 0.61044, 0.61533, 0.620276, 0.625278, 0.630337, 0.635455,
   0.640632, 0.645871, 0.659247, 0.673038, 0.687272, 0.701975, 
  0.717178, 0.732914, 0.74922, 0.766136, 0.783707, 0.801983, 0.821018,
   0.840873, 0.86162, 0.883336, 0.90611, 0.930046, 0.955261, 0.981893,
   1.0101, 1.04008, 1.07205, 1.10629, 1.14315, 1.18304, 1.22651, 
  1.27427, 1.32727, 1.38685, 1.41293, 1.44053, 1.46986, 1.50115, 
  1.53471, 1.57094, 1.61033, 1.65357, 1.70161, 1.75581, 1.81834, 
  1.89287, 1.93665, 1.98671, 2.04575, 2.08011, 2.1191, 2.13638, 
  2.15489, 2.17486, 2.19662, 2.22067, 2.24774, 2.2791, 2.29706, 
  2.31719, 2.34039, 2.35362, 2.36839, 2.38534, 2.38909, 2.39297, 
  2.39703, 2.40126, 2.40571, 2.41538, 2.42643, 2.4327, 
  2.43966}, {0.522196, 0.520579, 0.52057, 0.521087, 0.52093, 0.520981,
   0.521047, 0.521087, 0.521102, 0.52111, 0.521228, 0.521338, 
  0.521456, 0.521686, 0.521918, 0.522149, 0.522613, 0.523077, 
  0.523541, 0.524006, 0.524471, 0.524937, 0.525404, 0.526572, 
  0.527743, 0.530095, 0.53246, 0.534837, 0.539631, 0.544478, 0.549379,
   0.554335, 0.559348, 0.564418, 0.569547, 0.574737, 0.579987, 
  0.585301, 0.590679, 0.596123, 0.61003, 0.624382, 0.639206, 0.654532,
   0.670391, 0.686819, 0.703853, 0.721535, 0.739912, 0.759036, 
  0.778963, 0.799758, 0.821493, 0.844248, 0.868117, 0.893205, 
  0.919633, 0.947544, 0.977101, 1.0085, 1.04198, 1.07781, 1.11634, 
  1.15802, 1.20338, 1.25315, 1.30829, 1.37015, 1.39719, 1.42578, 
  1.45611, 1.48844, 1.52307, 1.56039, 1.6009, 1.64529, 1.69451, 
  1.74992, 1.81368, 1.88946, 1.93386, 1.98455, 2.0442, 2.07886, 
  2.11815, 2.13555, 2.15417, 2.17425, 2.19612, 2.22028, 2.24745, 
  2.2789, 2.29691, 2.31708, 2.34031, 2.35356, 2.36835, 2.38531, 
  2.38906, 2.39295, 2.39701, 2.40125, 2.4057, 2.41537, 2.42643, 
  2.43269, 2.43965}, 
{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}};

  const double imomegaqnm55[8][107] = {{0.0151971, 0.02987, 0.0432352, 0.064849, 0.0789971, 0.0834521, 
  0.0865117, 0.0885261, 0.0898184, 0.090638, 0.0919286, 0.0920896, 
  0.0921222, 0.0921485, 0.0921652, 0.0921765, 0.092191, 0.0922013, 
  0.0922106, 0.0922195, 0.0922283, 0.0922371, 0.0922459, 0.0922679, 
  0.0922898, 0.0923335, 0.092377, 0.0924205, 0.0925068, 0.0925926, 
  0.0926777, 0.092762, 0.0928457, 0.0929286, 0.0930106, 0.0930918, 
  0.0931722, 0.0932516, 0.09333, 0.0934074, 0.0935961, 0.0937773, 
  0.09395, 0.0941133, 0.094266, 0.0944066, 0.0945339, 0.094646, 
  0.0947411, 0.0948169, 0.0948709, 0.0949002, 0.0949014, 0.0948705, 
  0.094803, 0.0946932, 0.0945348, 0.0943198, 0.0940389, 0.0936804, 
  0.0932301, 0.09267, 0.0919774, 0.0911228, 0.0900671, 0.0887575, 
  0.0871199, 0.0850466, 0.0840612, 0.0829665, 0.0817451, 0.0803751, 
  0.0788289, 0.0770706, 0.0750526, 0.0727103, 0.0699513, 0.0666388, 
  0.0625545, 0.0573174, 0.054054, 0.0501531, 0.0453206, 0.0423931, 
  0.038969, 0.0374162, 0.0357305, 0.0338844, 0.0318402, 0.0295431, 
  0.0269089, 0.0237947, 0.0219803, 0.0199207, 0.017514, 0.0161252, 
  0.0145605, 0.0127472, 0.0123443, 0.0119247, 0.0114865, 0.011027, 
  0.0105435, 0.00948819, 0.00827437, 0.00758316, 
  0.00681219}, {0.25166, 0.264251, 0.270831, 0.275931, 0.277491, 
  0.277846, 0.278073, 0.278222, 0.278323, 0.278391, 0.278527, 
  0.278556, 0.278565, 0.278576, 0.278588, 0.2786, 0.278625, 0.27865, 
  0.278676, 0.278701, 0.278726, 0.278751, 0.278776, 0.278839, 
  0.278901, 0.279026, 0.27915, 0.279273, 0.279519, 0.279763, 0.280005,
   0.280244, 0.280481, 0.280716, 0.280948, 0.281177, 0.281404, 
  0.281627, 0.281847, 0.282065, 0.282593, 0.283097, 0.283575, 
  0.284023, 0.284438, 0.284817, 0.285154, 0.285444, 0.285682, 
  0.285862, 0.285974, 0.286011, 0.285963, 0.285817, 0.28556, 0.285175,
   0.284643, 0.28394, 0.283038, 0.281902, 0.280489, 0.278746, 
  0.276603, 0.273974, 0.27074, 0.266744, 0.261763, 0.255474, 0.252491,
   0.249179, 0.245488, 0.24135, 0.236685, 0.231383, 0.225304, 
  0.218251, 0.20995, 0.199989, 0.187715, 0.171984, 0.162185, 0.150475,
   0.135971, 0.127186, 0.116912, 0.112252, 0.107194, 0.101656, 
  0.0955223, 0.0886307, 0.0807275, 0.0713845, 0.0659411, 0.0597624, 
  0.0525421, 0.0483758, 0.0436815, 0.0382415, 0.037033, 0.0357743, 
  0.0344594, 0.0330811, 0.0316306, 0.0284646, 0.0248231, 0.0227495, 
  0.0204366}, {0.467579, 0.470028, 0.470942, 0.471348, 0.471322, 
  0.47129, 0.471261, 0.471239, 0.471223, 0.471213, 0.471201, 0.47121, 
  0.471221, 0.47124, 0.471258, 0.471277, 0.471314, 0.471351, 0.471388,
   0.471425, 0.471462, 0.471499, 0.471536, 0.471628, 0.47172, 
  0.471904, 0.472086, 0.472268, 0.472628, 0.472984, 0.473336, 
  0.473684, 0.474028, 0.474367, 0.474701, 0.475031, 0.475355, 
  0.475674, 0.475987, 0.476294, 0.477036, 0.477735, 0.478387, 
  0.478987, 0.479529, 0.480006, 0.480411, 0.480735, 0.480968, 
  0.481101, 0.481119, 0.481008, 0.480751, 0.480328, 0.479717, 
  0.478889, 0.477812, 0.476446, 0.474745, 0.472651, 0.470091, 
  0.466978, 0.463197, 0.4586, 0.452993, 0.446112, 0.437588, 0.426883, 
  0.421821, 0.416213, 0.409971, 0.402987, 0.395124, 0.386201, 
  0.375984, 0.364147, 0.350231, 0.333554, 0.313026, 0.286746, 
  0.270386, 0.250845, 0.226651, 0.212, 0.194868, 0.1871, 0.178668, 
  0.169434, 0.15921, 0.147722, 0.134549, 0.118976, 0.109903, 
  0.0996047, 0.0875705, 0.0806266, 0.0728027, 0.063736, 0.0617217, 
  0.0596238, 0.0574324, 0.0551353, 0.0527177, 0.047441, 0.0413719, 
  0.0379158, 0.034061}, {0.676801, 0.675633, 0.67493, 0.67443, 0.6744,
   0.674422, 0.674445, 0.674462, 0.674474, 0.674481, 0.67449, 
  0.674497, 0.674508, 0.674528, 0.674548, 0.674568, 0.674609, 
  0.674649, 0.674689, 0.674729, 0.674769, 0.674809, 0.674849, 
  0.674948, 0.675048, 0.675245, 0.675441, 0.675635, 0.67602, 0.676399,
   0.676771, 0.677137, 0.677496, 0.677847, 0.678192, 0.678529, 
  0.678858, 0.679179, 0.679492, 0.679795, 0.680514, 0.681168, 
  0.681753, 0.68226, 0.682681, 0.683008, 0.683228, 0.683332, 0.683303,
   0.683128, 0.682788, 0.682263, 0.681529, 0.680557, 0.679315, 
  0.677766, 0.675862, 0.673549, 0.67076, 0.667416, 0.663415, 0.658633,
   0.65291, 0.64604, 0.637751, 0.627674, 0.615292, 0.599856, 0.592588,
   0.584557, 0.575639, 0.565684, 0.554498, 0.541831, 0.527355, 
  0.510615, 0.490971, 0.467468, 0.438585, 0.401662, 0.378702, 
  0.351292, 0.317377, 0.296847, 0.272846, 0.261966, 0.250155, 
  0.237224, 0.222906, 0.20682, 0.188374, 0.16657, 0.153867, 0.139448, 
  0.1226, 0.112878, 0.101924, 0.0892306, 0.0864106, 0.0834736, 
  0.0804055, 0.0771896, 0.0738049, 0.0664175, 0.0579207, 0.0530822, 
  0.0476854}, {0.891895, 0.89147, 0.891863, 0.892311, 0.89227, 
  0.892236, 0.892217, 0.892213, 0.892217, 0.892223, 0.892242, 
  0.892246, 0.892254, 0.892269, 0.892283, 0.892298, 0.892327, 
  0.892357, 0.892386, 0.892415, 0.892444, 0.892473, 0.892502, 
  0.892574, 0.892646, 0.892788, 0.892927, 0.893065, 0.893335, 
  0.893597, 0.89385, 0.894094, 0.894329, 0.894554, 0.89477, 0.894975, 
  0.89517, 0.895353, 0.895526, 0.895686, 0.896031, 0.896291, 0.896458,
   0.896521, 0.89647, 0.896292, 0.895974, 0.895501, 0.894855, 
  0.894014, 0.892957, 0.891657, 0.890083, 0.888198, 0.88596, 0.88332, 
  0.880219, 0.876586, 0.872336, 0.867364, 0.86154, 0.854705, 0.846653,
   0.83712, 0.825754, 0.812083, 0.795442, 0.774871, 0.765236, 
  0.754619, 0.742864, 0.729775, 0.715106, 0.698537, 0.679644, 
  0.657847, 0.632324, 0.601851, 0.564477, 0.516789, 0.487173, 
  0.451846, 0.408167, 0.381741, 0.350856, 0.336857, 0.321664, 0.30503,
   0.286614, 0.265926, 0.242205, 0.214167, 0.197833, 0.179293, 
  0.15763, 0.14513, 0.131046, 0.114725, 0.1111, 0.107324, 0.103379, 
  0.099244, 0.0948923, 0.085394, 0.0744695, 0.0682486, 
  0.0613098}, {1.12614, 1.12697, 1.12642, 1.12612, 1.12629, 1.12629, 
  1.12627, 1.12625, 1.12624, 1.12624, 1.12626, 1.12625, 1.12626, 
  1.12626, 1.12626, 1.12626, 1.12626, 1.12626, 1.12627, 1.12627, 
  1.12627, 1.12627, 1.12628, 1.12628, 1.12628, 1.12629, 1.12629, 
  1.1263, 1.12629, 1.12628, 1.12625, 1.12622, 1.12617, 1.12611, 
  1.12604, 1.12595, 1.12586, 1.12574, 1.12562, 1.12548, 1.12505, 
  1.12453, 1.12388, 1.1231, 1.12219, 1.12111, 1.11986, 1.11841, 
  1.11675, 1.11485, 1.11267, 1.1102, 1.10738, 1.10418, 1.10054, 
  1.09641, 1.0917, 1.08634, 1.08021, 1.07319, 1.06513, 1.05582, 
  1.04501, 1.03238, 1.0175, 0.999787, 0.978439, 0.952276, 0.940091, 
  0.926705, 0.911928, 0.895521, 0.877184, 0.856528, 0.833035, 0.806, 
  0.774421, 0.736806, 0.690779, 0.632177, 0.595837, 0.552533, 
  0.499038, 0.466693, 0.428905, 0.411781, 0.393199, 0.372856, 
  0.350338, 0.325043, 0.296043, 0.261768, 0.241802, 0.219141, 
  0.192661, 0.177383, 0.160169, 0.140221, 0.135789, 0.131174, 
  0.126352, 0.121299, 0.11598, 0.104371, 0.0910184, 0.083415, 
  0.0749342}, {1.3749, 1.37446, 1.37529, 1.37504, 1.37506, 1.37511, 
  1.37511, 1.37508, 1.37506, 1.37506, 1.37507, 1.37505, 1.37505, 
  1.37503, 1.37501, 1.37499, 1.37496, 1.37492, 1.37488, 1.37485, 
  1.37481, 1.37478, 1.37474, 1.37465, 1.37456, 1.37437, 1.37418, 
  1.37399, 1.37359, 1.37319, 1.37277, 1.37233, 1.37188, 1.37142, 
  1.37094, 1.37045, 1.36993, 1.36941, 1.36886, 1.36829, 1.3668, 
  1.36517, 1.36339, 1.36145, 1.35934, 1.35704, 1.35452, 1.35177, 
  1.34875, 1.34544, 1.3418, 1.33779, 1.33338, 1.3285, 1.32309, 
  1.31709, 1.3104, 1.30292, 1.29452, 1.28507, 1.27436, 1.26216, 
  1.24818, 1.23202, 1.2132, 1.19101, 1.1645, 1.13228, 1.11735, 
  1.10101, 1.08301, 1.06309, 1.04089, 1.01595, 0.987663, 0.955194, 
  0.917365, 0.872418, 0.817554, 0.747869, 0.704729, 0.653377, 
  0.590005, 0.551715, 0.507, 0.486743, 0.464764, 0.440707, 0.41408, 
  0.384173, 0.34989, 0.309375, 0.285775, 0.25899, 0.227694, 0.209637, 
  0.189292, 0.165716, 0.160479, 0.155024, 0.149326, 0.143353, 
  0.137067, 0.123347, 0.107567, 0.0985815, 0.0885587}, 
{0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28}};

  /* To actually calculate the ringdown, we will use the following pointers to point to the correct mode */
  const double (*reomegaqnm)[107] = NULL;
  const double (*imomegaqnm)[107] = NULL;

  REAL8 totalMass, finalMass, finalSpin;

  /* Stuff for interpolating the data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  UINT4 i;

  /* Check we do not require more overtones than we have */
  if ( nmodes > 8 )
  {
    XLALPrintError( "Requesting more overtones than we have data to generate!\n");
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Choose the appropriate data bases on the user requested l an m */
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        reomegaqnm = reomegaqnm22;
        imomegaqnm = imomegaqnm22;
      }
      else if ( m == 1 )
      {
        reomegaqnm = reomegaqnm21;
        imomegaqnm = imomegaqnm21;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( l == 3 )
      {
        reomegaqnm = reomegaqnm33;
        imomegaqnm = imomegaqnm33;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL ); 
      }
      break;
    case 4:
      if ( l == 4 )
      {
        reomegaqnm = reomegaqnm44;
        imomegaqnm = imomegaqnm44;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( l == 5 )
      {
        reomegaqnm = reomegaqnm55;
        imomegaqnm = imomegaqnm55;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL );
      }
      break;
    default:
      XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }

  spline = gsl_spline_alloc( gsl_interp_cspline, 107 );
  acc    = gsl_interp_accel_alloc();

  totalMass = mass1 + mass2;

  /* Call XLALSimIMREOBFinalMassSpin() to get mass and spin of the final black hole */
  if ( XLALSimIMREOBFinalMassSpin(&finalMass, &finalSpin, mass1, mass2, spin1, spin2, approximant) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Now get the QNM frequencies from interpolating the above data */
  for ( i = 0; i < nmodes; i++ )
  {
    gsl_spline_init( spline, afinallist, reomegaqnm[i], 107 );
    gsl_interp_accel_reset( acc );
    
    modefreqs->data[i] = gsl_spline_eval( spline, finalSpin, acc );

    gsl_spline_init( spline, afinallist, imomegaqnm[i], 107 );
    gsl_interp_accel_reset( acc );

    modefreqs->data[i] += I*gsl_spline_eval( spline, finalSpin, acc );

    /* Scale by the appropriate mass factors */
    modefreqs->data[i] *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
  }
  /* Free memory and exit */
  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  return XLAL_SUCCESS;
}


