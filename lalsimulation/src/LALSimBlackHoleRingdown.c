#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

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
