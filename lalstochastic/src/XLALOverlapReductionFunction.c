/*
*  Copyright (C) 2007 Jolien Creighton, Robert Adam Mercer, John Whelan
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


#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <math.h>
#include <string.h>
#include <lal/StochasticCrossCorrelation.h>

static void evaluateBessels(REAL4 rho[3], REAL4 alpha);
static REAL4 cartesianInnerProduct(REAL4 a[3], REAL4 b[3]);

/**
 * \ingroup OverlapReductionFunction_c
 * \author UTB Relativity Group; contact whelan@phys.utb.edu
 *
 * \brief Calculates the values of the overlap reduction function for a pair
 * of gravitational wave detectors.
 *
 * XLALOverlapReductionFunction() calculates the values of the overlap reduction:
 *
 * \anchor stochastic_e_gamma \f{equation}{
 * \gamma(f):=\frac{5}{8\pi}\sum_A\int_{S^2}d\hat\Omega\
 * e^{i2\pi f\hat\Omega\cdot\Delta \vec x/c}\
 * F_1^A(\hat\Omega)F_2^A(\hat\Omega)\ ,
 * \tag{stochastic_e_gamma}
 * \f}
 *
 * where \f$\hat \Omega\f$ is a unit vector specifying a direction on the
 * two-sphere, \f$\Delta\vec x:=\vec x_1-\vec x_2\f$ is the separation vector
 * between the two detectors, and
 *
 * \anchor stochastic_e_F_i \f{equation}{
 * F_i^A(\hat\Omega):=e_{ab}^A(\hat\Omega)\ d_i^{ab}
 * \tag{stochastic_e_F_i}
 * \f}
 *
 * is the response of the \f$i\f$th detector \f$(i=1,2)\f$ to the \f$A=+,\times\f$
 * polarization.
 * Here \f$d_i^{ab}\f$ is the response tensor for the \f$i\f$th detector, which
 * relates the "strain" \f$h\f$ measured by the detector to the metric
 * perturbation \f$h_{ab}\f$ due to gravitational waves by
 *
 * \f{equation}{
 * h = d_i^{ab} h_{ab}
 * \f}
 *
 * The Cartesian components of \f$d_{ab}\f$ are constant in an earth-fixed
 * rotating coördinate system.  \f$\{e_{ab}^A(\hat\Omega)|A=+,\times\}\f$
 * are the spin-2 polarization tensors for the "plus" and "cross"
 * polarization states, normalized so that \f$e_{ab}^A e^{Bab}=2\delta^{AB}\f$.
 * With this definition,
 *
 * \f{equation}{
 * \gamma(f)=d_{1ab}d_2^{cd}\frac{5}{4\pi}\int_{S^2}d\hat\Omega\
 * e^{i2\pi f\hat\Omega\cdot\Delta \vec x/c}\
 * P^{ab}_{cd}(\hat\Omega)
 * \f}
 *
 * where \f$P^{ab}_{cd}(\hat\Omega)\f$ is a projection operator onto the space
 * of symmetric, traceless second-rank tensors orthogonal to \f$\hat\Omega\f$.
 *
 * The overlap reduction function for a pair of identical detectors is a
 * maximum when they are coïncident and coäligned; it
 * decreases when the detectors are shifted apart (so there is a phase
 * shift between the signals in the two detectors), or rotated out of
 * coälignment (so the detectors are sensitive to different
 * polarizations).  The overlap reduction function arises naturally when
 * calculating the cross-correlated signal due to an isotropic and
 * unpolarized stochastic gravitational-wave background.
 *
 * Given a choice of two detector sites, a frequency spacing \f$\delta f\f$,
 * a start frequency \f$f_0\f$, and the number of desired values \f$N\f$, <tt>
 * XLALOverlapReductionFunction()\/</tt> calculates the values of \f$\gamma(f)\f$ at the discrete
 * frequencies \f$f_i=f_0 + i\Delta f\f$, \f$i=0,1,\cdots, N-1\f$.
 *
 * ### Algorithm ###
 *
 * As shown in Appendix B of [\ref Flanagan1993] and Sec.~III.B of
 * [\ref Allen1999], the overlap reduction function can be written in
 * closed form in terms of the traceless response tensor
 * \f$D_{ab}=d_{ab}-\delta_{ab} d^c_c/3\f$
 * as sum of three spherical Bessel functions:
 *
 * \anchor stochastic_e_closed1 \f{equation}{
 * \gamma(f)=\rho_1(\alpha)\ D_1^{ab}D_{2ab}
 * +\rho_2(\alpha)\ D_1^{ab}D_{2a}{}^c s_b s_c
 * +\rho_3(\alpha)\ D_1^{ab}D_2^{cd}s_a s_b s_c s_d\ ,
 * \tag{stochastic_e_closed1}
 * \f}
 *
 * where
 *
 * \anchor stochastic_e_closed2 \f{equation}{
 * \left[
 * \begin{array}{c}
 * \rho_1\\
 * \rho_2\\
 * \rho_3
 * \end{array}
 * \right]
 * =
 * \frac{1}{2\alpha^2}
 * \left[
 * \begin{array}{rrr}
 * 10\alpha^2 & -20\alpha   & 10\\
 * -20\alpha^2 &  80\alpha   & -100\\
 * 5\alpha^2 & -50\alpha   & 175
 * \end{array}
 * \right]\
 * \left[
 * \begin{array}{c}
 * j_0\\
 * j_1\\
 * j_2
 * \end{array}
 * \right]\ ,
 * \tag{stochastic_e_closed2}
 * \f}
 *
 * \f$j_0\f$, \f$j_1\f$, and \f$j_2\f$ are the standard spherical Bessel functions:
 *
 * \f{eqnarray*}{
 * j_0(\alpha)&=&\frac{\sin\alpha}{\alpha} ,\\
 * j_1(\alpha)&=&\frac{\sin\alpha}{\alpha^2}-\frac{\cos\alpha}{\alpha}\ ,\\
 * j_2(\alpha)&=&3\ \frac{\sin\alpha}{\alpha^3}-3\ \frac{\cos\alpha}{\alpha^2} -\frac{\sin\alpha}{\alpha}\ ,
 * \f}
 *
 * \f$\vec s\f$ is a unit vector pointing in the direction of
 * \f$\Delta \vec x:=\vec x_1-\vec x_2\f$, and \f$\alpha:=2\pi f|\Delta\vec x|/c\f$.
 *
 * <tt>XLALOverlapReductionFunction()\/</tt> calculates the values of \f$\gamma(f)\f$
 * as follows:
 *
 * <ol>
 *
 * <li> Gets the locations and response tensors for the two detectors
 * from the \c LALDetector structures in the input.</li>
 *
 * <li> Constructs the traceless parts \f$D_{iab}\f$ of the two detector
 * response tensors and finds the distance \f$|\Delta\vec x|\f$ and
 * direction \f$s^a\f$ between the sites.</li>
 *
 * <li> Calculates the frequency-independent coëfficients
 * \f$D_1^{ab}D_{2ab}\f$, \f$D_1^{ab}D_{2a}{}^c s_b s_c\f$, and
 * \f$D_1^{ab}D_2^{cd}s_a s_b s_c s_d\f$ that appear in
 * Eq.\eqref{stochastic_e_closed1}.</li>
 *
 * <li> Calculates \f$\gamma(f)\f$ at each discrete frequency
 * \f$f_i:=f_0+i\Delta f\f$, \f$i=0,1,\cdots N-1\f$, using the power series
 * expansion
 * \f{eqnarray}{
 * j_0(\alpha) &=& 1 - \frac{\alpha^2}{6} + \frac{\alpha^4}{120} + \mathcal{O}(\alpha^6) \\
 * \frac{j_1(\alpha)}{\alpha} &=& \frac{1}{3} - \frac{\alpha^2}{30} + \frac{\alpha^4}{840} + \mathcal{O}(\alpha^6) \\
 * \frac{j_2(\alpha)}{\alpha^2} &=& \frac{1}{15} - \frac{\alpha^2}{210} + \frac{\alpha^4}{7560} + \mathcal{O}(\alpha^6)
 * \f}
 * for the spherical Bessel functions \f$j_0(\alpha_i)\f$,
 * \f$j_a(\alpha_i)\f$, \f$j_2(\alpha_i)\f$ when \f$\alpha_i=2\pi f_i
 * |\Delta\vec x|/c<0.01\f$.
 * </li>
 * </ol>
 *
 * ### Uses ###
 *
 * \code
 * LALUnitRaise()
 * sin()
 * cos()
 * sqrt()
 * strncpy()
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 *
 * <li> The \f$\gamma(f)\f$ here is related to the unnormalized \f$\Gamma(f)\f$
 * defined by Maggiore
 * [\ref Maggiore2000a,\ref Maggiore2000b] by
 * \f$\gamma(f) = \frac{5}{2}\Gamma(f)\f$.  This normalization, which
 * agrees with the literature
 * [\ref Flanagan1993,\ref Allen1997,\ref Allen1999]
 * on interferometers, is chosen so that \f$\gamma(f)\equiv 1\f$ for a pair
 * of coïncident, coäligned interferometers with perpendicular
 * arms.  It means that, for combinations other than a pair of
 * interferometers, our \f$\gamma(f)\f$ is \e not equal to the
 * generalization of \f$\gamma(f)\f$ defined by Maggiore, whose
 * relationship to \f$\Gamma(f)\f$ depends on the type of detector.
 * Defining \f$\gamma(f)\f$ as we do allows us to use the formulae from,
 * e.g., [\ref Allen1999], irrespective of the detector
 * type in question.</li>
 *
 * <li> While \f$\gamma(f)\f$ is usually considered to be dimensionless,
 * this routine attaches to it units of strain\f$^2\f$.  This is because it
 * contains two powers of the response tensor \f$d^{ab}\f$, which converts
 * the dimensionless metric perturbation \f$h_{ab}\f$ to \f$h=h_{ab}d^{ab}\f$,
 * which has units of strain.
 * </li>
 * </ul>
 *
 */

void
XLALOverlapReductionFunction(
    LALStatus                                *status,
    REAL4FrequencySeries                     *output,
    const LALDetectorPair                    *detectors,
    const XLALOverlapReductionFunctionParameters *parameters)

{
  UINT4 length;
  REAL8 deltaF;
  REAL8 f0;
  UINT4 i, j;
  REAL4 trace1, trace2;
  REAL4 d1DotS[3], d2DotS[3];
  REAL4 s[3];
  REAL4 distance;
  REAL4 c1, c2, c3;
  REAL4 alpha, alpha0, deltaAlpha;
  REAL4 rho[3];
  REAL4 d1[3][3], d2[3][3];
  RAT4 power;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* check that pointer to parameters is not null */
  ASSERT(parameters!=NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that specified length of output vector is > 0 */
  length = parameters->length;
  ASSERT(length > 0, status, \
			STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that frequency spacing is > 0 */
  deltaF = parameters->deltaF;
  ASSERT(deltaF > 0, status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that minimum frequency is >= 0 */
  f0 = parameters->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
				STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that pointer to output frequency series is not null */
  ASSERT(output != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of output frequency series is
	 * not null */
  ASSERT(output->data != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of the data member of output frequency series
	 * agrees with length specified in input parameters */
  if(output->data->length != length)
	{
		ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
				STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of output vector is not null */
  ASSERT(output->data->data != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to input structure is not null */
  ASSERT(detectors != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* everything okay here -------------------------------------------- */

  /* the overlap reduction function has units of strain^2 */
  power.numerator = 2;
  power.denominatorMinusOne = 0;
  TRY(LALUnitRaise(status->statusPtr, &(output->sampleUnits), \
				&lalStrainUnit, &power), status);

  /* set parameters for output */
  strncpy(output->name, "Overlap reduction function", LALNameLength);
  output->epoch.gpsSeconds = 0.0;
  output->epoch.gpsNanoSeconds = 0.0;
  output->deltaF = parameters->deltaF;
  output->f0 = parameters->f0;

  /* calculate separation vector between sites */
  for (i = 0; i < 3; i++)
	{
		s[i] = (REAL4)(detectors->detectorOne.location[i] - \
				detectors->detectorTwo.location[i]);
		d1[i][i] = detectors->detectorOne.response[i][i];
    d2[i][i] = detectors->detectorTwo.response[i][i];

		for (j = i; j<3; j++) {
      d1[i][j] = d1[j][i] = detectors->detectorOne.response[i][j];
      d2[i][j] = d2[j][i] = detectors->detectorTwo.response[i][j];

			/* check for non symmetric response tensor */
      ASSERT(d1[j][i] == detectors->detectorOne.response[j][i], status, \
					STOCHASTICCROSSCORRELATIONH_ENONSYMDIJ, \
					STOCHASTICCROSSCORRELATIONH_MSGENONSYMDIJ);
      ASSERT(d2[j][i] == detectors->detectorTwo.response[j][i], status, \
					STOCHASTICCROSSCORRELATIONH_ENONSYMDIJ, \
					STOCHASTICCROSSCORRELATIONH_MSGENONSYMDIJ);
		}
	}

	/* calculate distance between sites, in meters */
  distance = sqrt(cartesianInnerProduct(s,s));

  /* calculate unit separation vector */
  if (distance != 0)
	{
		for (i = 0; i < 3; i++)
		{
      s[i] /= distance;
    }
  }

  trace1 = d1[0][0] + d1[1][1] + d1[2][2];
  if (trace1)
	{
    trace1 /= 3.0;
    for (i = 0; i < 3; i++)
		{
      d1[i][i] -= trace1;
    }
  }
  trace2 = d2[0][0] + d2[1][1] + d2[2][2];
  if (trace2)
	{
    trace2 /= 3.0;
    for (i = 0; i < 3; i++)
		{
      d2[i][i] -= trace2;
    }
  }

  /* calculate coefficients c1, c2, c3 for overlap reduction funtion */

  /* c1 = d1 : d2 */
  c1 = 0;
  for (i = 0; i < 3; i++)
	{
    for (j = 0; j < 3; j++)
		{
      c1 += d1[i][j] * d2[i][j];
    }
  }

  for (i = 0; i < 3; i++)
	{
    d1DotS[i] = cartesianInnerProduct(d1[i], s);
    d2DotS[i] = cartesianInnerProduct(d2[i], s);
  }

  /* c2 = s . d1 . d2 . s */
  c2 = cartesianInnerProduct(d1DotS, d2DotS);

  /* c3 = (s . d1 . s)(s . d2 . s) */
  c3 = cartesianInnerProduct(s, d1DotS) * cartesianInnerProduct(s, d2DotS);

  distance *= (2 * LAL_PI / LAL_C_SI);
  deltaAlpha = deltaF * distance;
  alpha0 = f0 * distance;

  if (f0 == 0)
  {
    for (i = 0; i < length; ++i)
    {
      alpha = deltaAlpha * (REAL4)i;
      evaluateBessels(rho, alpha);
      output->data->data[i] = c1 * rho[0] + c2 * rho[1] + c3 * rho[2];
    }
  }
  else
  {
    for (i = 0; i < length; ++i)
    {
      alpha = alpha0 + deltaAlpha * (REAL4)i;
      evaluateBessels(rho, alpha);
      output->data->data[i] = c1 * rho[0] + c2 * rho[1] + c3 * rho[2];
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

static void evaluateBessels(REAL4 rho[3], REAL4 alpha)
{
  REAL8 alpha2, alpha4;
  REAL8 s, c;
  REAL8 b0, b1, b2;

  alpha2 = alpha * alpha;

  /* if the argument is close to zero, use power series */
  if (alpha < 0.01)
	{
		alpha4 = alpha2 * alpha2;
    b0 = 1.0 - alpha2/6.0 + alpha4/120.0;
    b1 = 1.0/3.0 - alpha2/30.0 + alpha4/840.0;
    b2 = 1.0/15.0 - alpha2/210.0 + alpha4/7560.0;
  }
  else
	{
		s = sin(alpha);
    c = cos(alpha);

    /* define spherical bessel functions j0, j1, j2 */

    b0 = s/alpha; /* = j0 */
    b1 = (b0 - c)/alpha;
    b2 = ((3.0 * b1) - s) / alpha;
    b1 /= alpha; /* = j1/alpha */
    b2 /= alpha2; /* = j2/alpha2 */
  }

  rho[0] = 5.0*b0 - 10.0*b1 + 5.0*b2;
  rho[1] = -10.0*b0 + 40.0*b1 - 50.0*b2;
  rho[2] = 2.5*b0 - 25.0*b1 + 87.5*b2;

  return;
}

static REAL4 cartesianInnerProduct(REAL4 a[3], REAL4 b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
