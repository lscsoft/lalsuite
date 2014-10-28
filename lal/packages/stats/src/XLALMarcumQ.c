/*
 * Copyright (C) 2014  Kipp Cannon
 *
 * Implementation of the algorithm described in
 *
 * A. Gil, J. Segura, and N. M. Temme.  Algorithm 939:  Computation of the
 * Marcum Q-Function.  ACM Transactions on Mathematical Software (TOMS),
 * Volume 40 Issue 3, April 2014, Article No. 20. arXiv:1311.0681
 *
 * with a few modifications.  In particular, a different expression is used
 * here for the 0th term in the asymptotic expansion for large (xy) than
 * shown in Section 4.1 such that it remains numerically stable in the x ~=
 * y limit, and a special case of the recursion implemented that remains
 * valid when x = y.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#pragma GCC diagnostic ignored "-Wshadow"


/*
 * ============================================================================
 *
 *                                   Preamble
 *
 * ============================================================================
 */


/*
 * stuff from the C library
 */


#include <math.h>


/*
 * stuff from GSL
 */


#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>


/*
 * stuff from LAL
 */


#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/LALMarcumQ.h>


/*
 * ============================================================================
 *
 *                                 Support Code
 *
 * ============================================================================
 */


/*
 * equation (8), normalized incomplete Gamma function,
 *
 *	Q_{\mu}(x) = \Gamma(\mu, x) / \Gamma(\mu).
 *
 * we use the GSL implementations, but Gil et al.'s tests of their Marcum Q
 * function implemention rely on their own implementation of these
 * functions.  GSL's implementations might not work as well as theirs and
 * so the Marcum Q implemented here might not meet their accuracy claims.
 * on the other hand, GSL might be using exactly their aglorithm, or it
 * might be using an even better one, who knows.
 *
 * their algorithm is described in
 *
 * Gil , A., Segura , J., and Temme , N. M. 2012. Efficient and accurate
 * algorithms for the computation and inversion of the incomplete gamma
 * function ratios. SIAM J. Sci. Comput.  34(6), A2965–A2981.
 * arXiv:1306.1754.
 */


static double Q(double mu, double x)
{
	return gsl_sf_gamma_inc_Q(mu, x);
}


/*
 * \sqrt{x / y} * equation (14).  see
 *
 * W. Gautschi, J. Slavik, On the Computation of Modified Bessel Function
 * Ratios, Mathematics of Computation, Vol. 32, No. 143 (Jul., 1978), pp.
 * 865-875.
 *
 * for information on using continued fractions to compute the ratio
 * directly.
 */


static double c_mu(double mu, double xi)
{
	return gsl_sf_bessel_Inu_scaled(mu, xi) / gsl_sf_bessel_Inu_scaled(mu - 1., xi);
}


/*
 * log of equation (32).  computed using log \Gamma functions from GSL.
 */


static double lnA(int n, double mu)
{
	mu += 0.5;
	return gsl_sf_lngamma(mu + n) - gsl_sf_lngamma(mu - n) - n * log(2.) - gsl_sf_lngamma(n + 1);
}


/*
 * equation (84).  \zeta^{2} / 2 = phi(xi) - phi(z0), with special case for
 * y ~= x + 1 given by expansion in (85) and (86).
 */


static double half_zeta_2(double x, double y)
{
	if(fabs(y - x - 1.) < 1e-3) {
		/*
		 * we do this a little differently than implied by what's
		 * in Gil et al.  we rewrite (85) as
		 *
		 * -(2x + 1)^(3/2) * (y - x - 1) / (2x+1)^2 * \sum_{k} ...
		 * = -(2x + 1)^(3/2) * \sum_{k} c_{k} z^{k + 1}
		 *
		 * but then we're returning \zeta^{2} / 2,
		 *
		 * = (2x + 1)^3 * [\sum_{k} c[k] z^{k+1}]^2 / 2.
		 */
		double c[] = {
			+1.,	/* not used */
			-(3. * x + 1.) / 3.,
			+((72. * x + 42.) * x + 7.) / 36.,
			-(((2700. * x + 2142.) * x + 657.) * x + 73.) / 540.,
			+((((181440. * x + 177552.) * x + 76356.) * x + 15972.) * x + 1331.) / 12960.
		};
		double z = (y - x - 1.) / pow(2. * x + 1., 2.);
		double z_to_k_p_1 = z;	/* = z^{1} */
		double sum = z_to_k_p_1;	/* = c[0] * z^{1} */
		int k;
		for(k = 1; k < 5; k++) {
			z_to_k_p_1 *= z;
			sum += c[k] * z_to_k_p_1;
		}
		return pow(2. * x + 1., 3.) * sum * sum / 2.;
	} else {
		double root_1_plus_4xy = sqrt(1. + 4. * x * y);
		return x + y - root_1_plus_4xy + log((1. + root_1_plus_4xy) / (2. * y));
	}
}


/*
 * equation (56).  = \sqrt{2 (phi(xi) - phi(z0))} with the sign of x + 1 -
 * y.  phi(xi) - phi(z0) is given by the expression in equation (84).
 */


static double zeta(double x, double y)
{
	return copysign(sqrt(2. * half_zeta_2(x, y)), x + 1. - y);
}


/*
 * equation (98).  note:  rho, r, and r'sin(theta) are all even functions
 * of theta.  we assume 0 <= theta <= pi
 */


static double theta_over_sin_theta(double theta)
{
	/* Taylor series = 1 + theta^2 / 6 + 7 theta^4 / 360.  therefore,
	 * once theta is below 10^-4, the first two terms are sufficient */
	return theta < 1e-4 ? theta * theta / 6. + 1. : theta / sin(theta);
}


static double theta_over_sin_theta_primed_sin_theta(double theta)
{
	/* sin x * d/dx (x / sin x)
	 *	= (sin x - x cos x) / sin x,
	 *	= 1 - x / tan x
	 *
	 * Taylor series is
	 *
	 * 	theta^2 / 3 + theta^4 / 45 + 2 theta^6 / 945 + ...
	 *
	 * therefore once theta is below 10^-4 the first two terms are
	 * sufficient
	 */

	return theta < 1e-4 ? (theta * theta / 45. + 1. / 3.) * theta * theta : 1. - theta / tan(theta);
}


static double rho(double theta_over_sin_theta, double xi)
{
	return sqrt(theta_over_sin_theta * theta_over_sin_theta + xi * xi);
}


static double r(double theta, double y, double xi)
{
	double theta_over_sin_theta_ = theta_over_sin_theta(theta);

	xi /= theta_over_sin_theta_;
	return (1. + sqrt(1. + xi * xi)) * theta_over_sin_theta_ / (2. * y);
}


static double r_primed_sin_theta(double theta, double y, double xi)
{
	double theta_over_sin_theta_ = theta_over_sin_theta(theta);

	xi /= theta_over_sin_theta_;
	return (1. + 1. / sqrt(1. + xi * xi)) * theta_over_sin_theta_primed_sin_theta(theta) / (2. * y);
}


/*
 * equation (96).  f is an even function of theta.
 *
 * NOTE:  there is a typo in both the final printed version and in the copy
 * on arXiv, the middle term in the denominator should read "- 2 r(\theta)
 * \cos(\theta)".
 *
 * we rewrite the denominator as (r - cos(theta))^2 + 1 - cos^2(theta) =
 * ((r - cos(theta))^2 + sin^2(theta)
 */


static double f(double theta, double y, double xi)
{
	double r_ = r(theta, y, xi);
	double r_minus_cos_theta = r_ - cos(theta);
	double sin_theta = sin(theta);

	return (r_primed_sin_theta(theta, y, xi) - r_minus_cos_theta * r_) / (r_minus_cos_theta * r_minus_cos_theta + sin_theta * sin_theta);
}


/*
 * equation (97).  psi is an even function of theta.
 */


static double psi(double theta, double xi)
{
	double theta_over_sin_theta_ = theta_over_sin_theta(theta);
	double rho_ = rho(theta_over_sin_theta_, xi);
	double root_xi2_plus_1 = sqrt(1. + xi * xi);

	return cos(theta) * rho_ - root_xi2_plus_1 - log((theta_over_sin_theta_ + rho_) / (1. + root_xi2_plus_1));
}


/*
 * equation (100), unscaled variables
 */


static void f1_f2(double x, double M, double *f1, double *f2)
{
	double a = x + M;
	double b = sqrt(4. * x + 2. * M);
	*f1 = a - b;
	*f2 = a + b;
}


/*
 * ============================================================================
 *
 *                         Implementations by Domain
 *
 * ============================================================================
 */


/*
 * series expansion described in section 3.  equation (7).
 */


static double MarcumQ_small_x(double M, double x, double y)
{
	int n = 0;
	double x_to_n_over_n_factorial = 1.;
	double sum = 0.;

	do {
		double term = x_to_n_over_n_factorial * Q(M + n, y);
		sum += term;
		if(term <= 1e-17 * sum)
			break;
		n++;
		x_to_n_over_n_factorial *= x / n;
	} while(1);

	return exp(-x) * sum;
}


/*
 * asymptotic expansions in section 4.1.  equation (37) for y > x, 1 -
 * equation (39) for x < y, with modifications for x == y.
 */


static double MarcumQ_large_xy(double M, double x, double y, double xi)
{
	double root_y_minus_root_x = sqrt(y) - sqrt(x);

	/* equation (31) */
	double sigma = root_y_minus_root_x * root_y_minus_root_x / xi;
	double rho = sqrt(y / x);

	double rho_to_M_over_root_8_pi = pow(y / x, M / 2.) / sqrt(8. * LAL_PI);
	double e_to_neg_sigma_xi = exp(-root_y_minus_root_x * root_y_minus_root_x);

	/* Phi_{0}.  equation (35).  NOTE:  that in (35) they show
	 * \sqrt{\sigma \xi} being replaced by \sqrt{y} - \sqrt{x}, but it
	 * should be |\sqrt{y} - \sqrt{x}|.  the x ~= y limit is obtained
	 * from
	 *
	 *	erfc(a x) / x = 1/x - 2/\sqrt{pi} [ a / (1 0!) - a^3 x^2 / (3 1!) + a^5 x^4 / (5 2!) - a^7 x^6 / (7 3!) + ...]
	 *
	 * the ratio |2nd term| / |0th term| is < 10^-17 when ax < 3 *
	 * 10^-6.  using this,
	 *
	 *	Phi_{0} = sqrt(pi / sigma) - 2 sqrt(xi)
	 *
	 * to the limits of double precision when |sqrt(y)-sqrt(x)| < 3e-6.
	 * we cheat a bit and call that 1e-5 because our final answer isn't
	 * going to be perfect anyway;  and we compute the ratio of
	 * sqrt(pi)/sqrt(sigma), instead of sqrt(pi/sigma), to skirt
	 * numerical overflow until sigma is identically 0.
	 *
	 * the special case for \sigma = 0 is found by substituting
	 * \Phi_{0} from (35) into (36) and determining the desired value
	 * of \Phi_{1}.  when this is done we find a net positive power of
	 * \sigma in the term involving \Phi_{0} and so evidently the
	 * correct value for \Phi_{1} will be obtained by setting \Phi_{0}
	 * to any finite value and allowing the factor of \sigma appearing
	 * in (36) to eliminate it naturally when \sigma = 0 (if we let
	 * \Phi_{0} diverge we end up with inf*0 = nan and the code fails).
	 */
	double Phi = fabs(root_y_minus_root_x) < 1e-5 ? sigma == 0. ? 0. : sqrt(LAL_PI) / sqrt(sigma) - 2. * sqrt(xi) : sqrt(LAL_PI / sigma) * gsl_sf_erfc(fabs(root_y_minus_root_x));

	/* Psi_{0} from Phi_{0}.  this is a special case of equation (38)
	 * that remains valid in the small \rho-1 limit (x ~= y).  when n =
	 * 0, from equation (32) we see that A_{0}(\mu) = 1, and
	 * equation (38) can be written
	 *
	 *	\Psi_{0} = \rho^{\mu - 1} / \sqrt{8\pi} (\rho-1) * \Phi_{0}.
	 *
	 * meanwhile, from equation (31) we have that
	 *
	 *	\sqrt{\sigma} = |\rho - 1| / \sqrt{2 \rho}
	 *
	 * and so equation (35) can be written
	 *
	 *	 \Phi_{0} = \sqrt{2 \pi \rho} / |\rho - 1| erfc(\sqrt{y} - \sqrt{x}).
	 *
	 * combining these we find
	 *
	 *	\Psi_{0} = \rho^{\mu - .5} / 2 erfc(\sqrt{y} - \sqrt{x}) * sgn(\rho - 1).
	 *
	 * NOTE:  this result is only used for x < y (\rho > 1).  for x > y
	 * (\rho < 1) we compute Q=1-P where from equation (39) -P is the
	 * sum of -\tilde{\Psi}_{n} = \Psi_{n} except -\tilde{\Psi}_{0}
	 * which is given by -1 * equation (40).  comparing -1 * (40) to
	 * our result above when \rho < 1, we see that they are nearly
	 * identical except for the sign of the argument of erfc() being
	 * flipped (making it positive, again).  therefore, by taking the
	 * absolute value of the erfc() argument we obtain an expression
	 * for the 0th term in the sum valid in both cases.  the only
	 * remaining difference is that the sum needs to be initialized to
	 * +1 in the x > y case and 0 in the x < y case.
	 *
	 * if rho is exactly equal to 1, instead of relying on FPUs
	 * agreeing that 1-1=+0 and not -0 we enter the result with the
	 * correct sign by hand.  the "correct" sign is decided by the
	 * initialization of the sum.  we've chosen to initialize the x==y
	 * case to 0, so the correct sign for Psi_{0} is positive.
	 */
	double Psi = rho == 1. ? .5 : copysign(pow(rho, M - .5) * gsl_sf_erfc(fabs(root_y_minus_root_x)) / 2., rho - 1.);

	double sum = x > y ? 1. : 0.;
	int n = 0;

	do {
		/* equation (37) */
		sum += Psi;
		if(fabs(Psi) <= 1e-17 * fabs(sum))
			break;

		/* next n.  rho_to_M_over_root_8_pi is providing the factor
		 * of (-1)^n, so we flip its sign here, too.  */
		n++;
		rho_to_M_over_root_8_pi = -rho_to_M_over_root_8_pi;

		/* Phi_{n} from Phi_{n - 1}.  equation (36) */
		Phi = (e_to_neg_sigma_xi * pow(xi, .5 - n) - sigma * Phi) / (n - .5);

		/* Psi_{n} from Phi_{n}.  equation (38) */
		double lnA_n_M_minus_1 = lnA(n, M - 1.);
		Psi = rho_to_M_over_root_8_pi * exp(lnA_n_M_minus_1) * (1. - exp(lnA(n, M) - lnA_n_M_minus_1) / rho) * Phi;
	} while(1);

	/*
	 * for very small probabilities, improper cancellation can leave us
	 * with a small negative number.  nudge back into the allowed range
	 */

	if(-1e-200 < sum && sum < 0.)
		sum = 0.;

	return sum;
}


/*
 * recurrence relation in equation (14).
 */


static double MarcumQ_recurrence(double M, double x, double y, double xi)
{
	/* factor in equation (14) we've omitted from our c_mu() */
	double root_y_over_x = sqrt(y / x);

	/* where to start evaluating recurrence relation */
	double mu = sqrt(2. * xi) - 1.;
	/* make sure it differs from M by an integer */
	mu = M - ceil(M - mu);

	double Qmu_minus_1 = XLALMarcumQmodified(mu - 1., x, y);
	double Qmu = XLALMarcumQmodified(mu, x, y);

	do {
		double cmu = root_y_over_x * c_mu(mu, xi);

		/* compute Q_{\mu + 1} from Q_{\mu} and Q_{\mu - 1} */
		double Qmu_saved = Qmu;
		Qmu = (1. + cmu) * Qmu - cmu * Qmu_minus_1;
		Qmu_minus_1 = Qmu_saved;
		mu++;
	} while(mu < M);

	if(fabs(mu - M) / M > 1e-15)
		XLAL_ERROR_REAL8(XLAL_ELOSS, "recursion failed to land on correct order:  wanted M=%.16g, got M=%.16g", M, mu);

	return Qmu;
}


/*
 * asymptotic expansion in section 4.2
 */


static double MarcumQ_large_M(double M, double x, double y)
{
	/* equation (56) */
	double zeta_ = zeta(x, y);

	/* equation (67).  if size of Psi array is changed, updated safety
	 * check inside loop below */
	double e_to_neg_half_M_zeta_2 = exp(-M * half_zeta_2(x, y));
	double Psi[20];
	Psi[0] = sqrt(LAL_PI / (2. * M)) * gsl_sf_erfc(-zeta_ * sqrt(M / 2.));
	Psi[1] = e_to_neg_half_M_zeta_2 / M;

	double sum = 0.;
	int k = 1;

	do {
		/* equation (71) */
		double Bk = 0.;
		int j;
		for(j = 0; j <= k; j++)
			Bk += /* FIXME:  need f(j, k - j) * */ Psi[j] / pow(M, k - j);

		sum += Bk;
		if(Bk <= 1e-17 * sum)
			break;

		/* next k */
		k++;

		/* equation (68).  Psi_{j} from Psi_{j - 2}.  note that we
		 * have all j < k, we just need to add the j = k-th term to
		 * the sequence */
		if(k >= 20)
			/* FIXME:  allow dynamic allocation of Psi */
			XLAL_ERROR_REAL8(XLAL_EMAXITER);
		Psi[k] = (k - 1) / M * Psi[k - 2] + pow(-zeta_, k - 1) / M * e_to_neg_half_M_zeta_2;
	} while(1);

	/* equation (72) */
	return gsl_sf_erfc(-zeta_ * sqrt(M / 2.)) / 2. - sqrt(M / LAL_TWOPI) * sum;
}


/*
 * quadrature method in section 5.  the integrand is an even function of
 * theta, so we just integrate over [0, +\pi] and double the result.
 *
 * NOTE:  as theta approaches +\pi, \theta/sin(\theta) diverges in equation
 * (97) causing the cos(\theta)\rho term to go to -\infty and the argument
 * of the logarithm to go to +\infy and thus the logarithm to -\infty;
 * this makes the integrand go to 0 in (95), and so we can back the upper
 * bound away from +\pi a bit without affecting the result and in doing so
 * avoid the need to trap overflow errors inside the argument of the
 * exponential.
 *
 * NOTE:  the expression in (95) only yields Q if \zeta(x, y) < 0,
 * otherwise it yields -P.
 */


struct integrand_params {
	double M;
	double y;
	double xi;
};


static double integrand(double theta, void *integrand_params)
{
	struct integrand_params params = *(struct integrand_params *) integrand_params;

	/*
	 * equation (95)
	 */

	return exp(params.M * psi(theta, params.xi)) * f(theta, params.y, params.xi);
}


static double MarcumQ_quadrature(double M, double x, double y, double xi)
{
	struct integrand_params params = {
		.M = M,
		.y = y,
		.xi = xi
	};
	gsl_function integrand_ = {
		.function = integrand,
		.params = &params
	};
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(20);
	double integral = 0.;
	double abserr = 0.;

	/* "limit" must = size of workspace.  upper bound a bit less than
	 * +\pi (see above) */
	gsl_integration_qag(&integrand_, 0, LAL_PI - 1./512., abserr, 1e-12, 20, GSL_INTEG_GAUSS41, workspace, &integral, &abserr);
	gsl_integration_workspace_free(workspace);

	/* instead of dividing by 2*pi and doubling the integral, divide by
	 * pi */
	integral *= exp(-M * half_zeta_2(x, y)) / LAL_PI;

	if(x + 1. < y )	/* if zeta < 0 */
		return integral;
	return 1 + integral;
}


/*
 * ============================================================================
 *
 *                                 Exported API
 *
 * ============================================================================
 */


/**
 * The modified form of the Marcum Q function used by Gil <I>et al.</I> in
 *
 * A. Gil, J. Segura, and N. M. Temme.  Algorithm 939:  Computation of the
 * Marcum Q-Function.  ACM Transactions on Mathematical Software (TOMS),
 * Volume 40 Issue 3, April 2014, Article No. 20. arXiv:1311.0681
 *
 * The relationship between this function and the standard Marcum Q
 * function is
 *
 *	XLALMarcumQmodified(M, x, y) = XLALMarcumQ(M, sqrt(2. * x), sqrt(2. * y)).
 *
 * The function is defined for \f$1 \leq M\f$, \f$0 \leq x\f$, \f$0 \leq
 * y\f$.  Additionally, the implementation here becomes inaccurate when
 * \f$M\f$, \f$x\f$, or \f$y\f$ is \f$\geq 10000\f$.
 */


double XLALMarcumQmodified(double M, double x, double y)
{
	double xi;
	double f1, f2;
	double Q;

	/*
	 * check input
	 */

	if(M < 1.)
		XLAL_ERROR_REAL8(XLAL_EDOM, "require 1 <= M: M=%.16g", M);
	if(x < 0.)
		XLAL_ERROR_REAL8(XLAL_EDOM, "0 <= x: x=%.16g", x);
	if(y < 0.)
		XLAL_ERROR_REAL8(XLAL_EDOM, "0 <= y: y=%.16g", y);

	if(M > 10000.)
		XLAL_ERROR_REAL8(XLAL_ELOSS, "require M <= 10000: M=%.16g", M);
	if(x > 10000.)
		XLAL_ERROR_REAL8(XLAL_ELOSS, "require x <= 10000: x=%.16g", x);
	if(y > 10000.)
		XLAL_ERROR_REAL8(XLAL_ELOSS, "require y <= 10000: y=%.16g", y);

	/*
	 * compute some parameters
	 */

	xi = 2. * sqrt(x * y);
	f1_f2(x, M, &f1, &f2);

	/*
	 * selection scheme from section 6
	 */

	if(x < 30.) {
		/*
		 * use the series expansion in section 3
		 */

		Q = MarcumQ_small_x(M, x, y);
	} else if(xi > 30. && M * M < 2. * xi) {
		/*
		 * use the asymptotic expansion in section 4.1
		 */

		Q = MarcumQ_large_xy(M, x, y, xi);
	} else if(f1 < y && y < f2 && M < 135.) {
		/*
		 * use the recurrence relations in (14)
		 */

		Q = MarcumQ_recurrence(M, x, y, xi);
	} else if(f1 < y && y < f2 && M >= 135.) {
		/*
		 * use the asymptotic expansion in section 4.2
		 */

		/* FIXME:  missing an implementation of f_{jk} needed by
		 * this code path */
		XLAL_ERROR_REAL8(XLAL_EDOM, "not implemented: %s(%.16g, %.16g, %.16g)", __func__, M, x, y);

		Q = MarcumQ_large_M(M - 1., x / (M - 1.), y / (M - 1.));
	} else {
		/*
		 * use the integral representation in section 5
		 */

		Q = MarcumQ_quadrature(M, x / M, y / M, xi / M);
	}

	/*
	 * we're generally only correct to 12 digits, and the error can
	 * result in an impossible probability.  nudge back into the
	 * allowed range.
	 */

	if(1. < Q && Q < 1. + 1e-12)
		Q = 1.;

	if(Q < 0. || Q > 1.)
		XLAL_ERROR_REAL8(XLAL_ERANGE, "%s(%.17g, %.17g, %.17g) = %.17g", __func__, M, x, y, Q);

	return Q;
}


/**
 * The function defined by J.\ Marcum,
 * \f[
 *	Q_{M}(a, b) = \int_{b}^{\infty} x \left( \frac{x}{a} \right)^{M - 1} \exp \left( -\frac{x^{2} + a^{2}}{2} \right) I_{M - 1}(a x) \,\mathrm{d}x,
 * \f]
 * where \f$I_{M - 1}\f$ is the modified Bessel function of order \f$M -
 * 1\f$.
 *
 * The CCDF for the random variable \f$x\f$ distributed according to the
 * noncentral \f$\chi^{2}\f$ distribution with \f$k\f$ degrees-of-freedom
 * and noncentrality parameter \f$\lambda\f$ is \f$Q_{k/2}(\sqrt{\lambda},
 * \sqrt{x})\f$.
 *
 * The CCDF for the random variable \f$x\f$ distributed according to the
 * Rice distribution with noncentrality parameter \f$\nu\f$ and width
 * \f$\sigma\f$ is \f$Q_{1}(\nu/\sigma, x/\sigma)\f$.
 *
 * The probability that a signal that would be seen in a two-phase matched
 * filter with SNR \f$\rho_{0}\f$ is seen to have matched filter SNR
 * \f$\geq \rho\f$ in stationary Gaussian noise is \f$Q_{1}(\rho_{0},
 * \rho)\f$.
 *
 * This function is implemented by computing the modified form used by Gil
 * <I>et al.</I>,
 *
 *	XLALMarcumQ(M, a, b) = XLALMarcumQmodified(M, a * a / 2., b * b / 2.).
 */


double XLALMarcumQ(double M, double a, double b)
{
	/*
	 * inverse of equation (5).  note:  all error checking handled
	 * inside XLALMarcumQmodified().
	 */

	return XLALMarcumQmodified(M, a * a / 2., b * b / 2.);
}
