/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, Peter Shawhan, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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

/*
	Created: 7.9.96.
	Author: B.S.Sathyaprakash, Caltech, Cardiff University.
	Purpose: To compute the metric and the template bank
              parameters, corresponding to 2PN chirps.
	Revision History: Updates 15.7.97; 31.8.97.; 18.3.99; 25.05.02
        First C version October 2000.
        Major revision in May 2002: Direct computation in (tau0, tau3)
        space avoid (m, eta) space. This means the code won't create
        template parameters in (tau0, tau2) space. An interface to the
        old code will have to be provided, if that is necessary.

	Dependencies: moments.f, inverse.f transform.f (not needed in the version of 02/05)
	Outputs:
	    g00:
	    g11:
	  theta: Angle which the t0-axis makes with semi-major (dx0) axis.
         srate: The minimal sampling rate required (computed but not outputted.
	Notes: Owen and Sathyaprakash (Caltech collaboration notes).
              Also Sathyaprakash, note from October 2000, May 24/25, 2002.
*/

#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

#define METRIC_DIMENSION 2
#define METRIC_ORDER 5

static void
InspiralComputeMetricGetPsiCoefficients (
    REAL8              Psi[METRIC_DIMENSION][METRIC_ORDER],
    InspiralMomentsEtc *moments,
    REAL8 fLower,
    REAL8 t0,
    REAL8 t3
    );

/**
 * \ingroup LALInspiralBank_h
 * \deprecated Use XLALInspiralComputeMetric() instead.
 */
void
LALInspiralComputeMetric (
    LALStatus          *status,		/**< LAL status pointer */
    InspiralMetric     *metric,		/**< [out] the metric at the lattice point defined by \c params */
    InspiralTemplate   *params,		/**< [in] the parameters where metric must be computed */
    InspiralMomentsEtc *moments		/**< [in] moments \f$J(1), \ldots, J(17),\f$ of the PSD and other constants needed in the computation of the metric */
    )
{
  INITSTATUS(status);

  XLALPrintDeprecationWarning("LALInspiralComputeMetric", "XLALInspiralComputeMetric");

  if (XLALInspiralComputeMetric(metric, moments, params->fLower, params->order, params->t0, params->t3) == XLAL_FAILURE){
    ABORTXLAL( status );
  };

  RETURN( status );
}


/**
 * \ingroup LALInspiralBank_h
 * \brief Function to compute the components of the metric which is used to describe distances on the signal manifold.
 * \author Churches, D. K., Cokelaer, T., Sathyaprakash, B. S.
 *
 * We calculate the components of the metric using the procedure outlined
 * in Owen [\ref Owen_96].
 * This uses the moments of the noise curve,
 * \f{equation}{
 * I(q) \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{S_{h}(x)}
 * \, dx
 * \f}
 * and
 * \f{equation}{
 * J(q) \equiv \frac{I(q)}{I(7)} \,.
 * \f}
 * (Please note that the function \c LALInspiralMoments doesn't compute \f$I(q)\f$ defined
 * here; the index \f$q\f$ is definted differently there and the normalisation is supplied by the user.
 * For ease of writing, here we shall follow the standard notation.)
 * Then the moment functional \f$\mathcal{J}\f$ is defined such that, for a function \f$a\f$,
 * \f{equation}{
 * \mathcal{J} [a] \equiv \frac{1}{I(7)} \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-7/3}}{S_{h}(x)}
 * a(x) \, dx
 * \f}
 * which gives us
 * \f{equation}{
 * \mathcal{J} = \left[ \sum_{n} a_{n} x^{n} \right] = \sum_{n} a_{n} J(7-3n).
 * \f}
 * The above equation is used to calculate the components of the metric using the following formula:
 * \f{equation}{
 * \gamma_{\alpha \beta} = \frac{1}{2} \left( \mathcal{J} [ \psi_{\alpha} \psi_{\beta} ] -
 * \mathcal{J} [ \psi_{\alpha} ] \mathcal{J} [ \psi_{\beta} ] \right)
 * \f}
 * where \f$\psi_\alpha\f$ is the derivative of the Fourier phase of the inspiral waveform
 * with respect to the parameter \f$\lambda^\alpha,\f$ that is
 * \f$\psi_\alpha \equiv \Psi_{,\alpha}.\f$  Writing the derivative index as
 * \f$\alpha=0,j,\f$ with \f$j=1,\ldots,n\f$ we have
 * \f{equation}{
 * \psi_{0} \equiv 2 \pi f, \ \ \
 * \psi_{j} \equiv \frac{\partial \Delta \Psi}{\partial \Delta \lambda^{j}}.
 * \f}
 * The phase \f$\Psi\f$ is that which appears in the usual stationary-phase formula for
 * the Fourier transform:
 * \f{equation}{
 * \tilde{h}(f)  \propto f^{-7/6} e^{i[-\pi/4 - \Phi_{0} + 2 \pi f t_{0} + \Psi(f;\vec{\lambda}]}.
 * \f}
 * If we take the usual chirp times and multiply each by \f$(\pi f_{0})\f$ then we get
 * <em>dimensionless
 * chirp times</em> \f$\tau_{k}\f$, and then the phase \f$\Psi\f$ may be written in the form
 * \f{equation}{
 * \Psi = 2 \pi f t_{c} + \sum_{k} \Psi_{k}(f) \tau_{k}
 * \f}
 * where, defining \f$v_0 = (\pi m f_0)^{1/3}\f$ (\f$m\f$ being total mass and \f$f_0\f$ a
 * fiducial starting frequency), the chirptimes \f$\tau_{k},\f$ up to 2nd PN order,
 * are given by
 * \f[
 * \tau_{0} = \frac{5}{256 \eta v_{0}^{5}},\ \
 * \tau_{2} = \frac{5}{192 \eta v_{0}^{3}} \left( \frac{743}{336} + \frac{11}{4} \eta \right),
 * \f]
 * \f{equation}{
 * \tau_{3} = \frac{\pi}{8 \eta v_{0}^{2}},\ \
 * \tau_{4} = \frac{5}{128 \eta v_{0}} \left( \frac{3\,058\,673}{1\,016\,064} + \frac{5429}{1008}
 * \eta +
 * \frac{617}{144} \eta^{2} \right).
 * \f}
 * Up to second post-Newtonian approximation the \f$\Psi_{k}\f$ are given by
 * \f{equation}{
 * \Psi_{0} = \frac{6}{5 \nu^{5/3}},\ \
 * \Psi_{2} = \frac{2}{\nu},\ \
 * \Psi_{3} = - \frac{3}{\nu^{2/3}},\ \
 * \Psi_{4} = \frac{6}{\nu^{1/3}}.
 * \f}
 * where \f$\nu = f/f_{0}\f$.
 *
 * If we now make the substitution \f$f = v^{3}/\pi m\f$ we then the find that the phase may be
 * expressed in a simpler form
 * \f{equation}{
 * \Psi(v) = 2 \pi f t_{c} + \sum_{k} \theta_{k} v^{k-5}
 * \f}
 * where the <em>chirp parameters</em> \f$\theta_{k}\f$ are given by
 * \f[
 * \theta_{0} = \frac{3}{128 \eta}, \ \
 * \theta_{2} = \frac{5}{96 \eta} \left( \frac{743}{336} + \frac{11}{4} \eta \right),
 * \f]
 * \f{equation}{
 * \theta_{3} = - \frac{3 \pi}{8 \eta},\ \
 * \theta_{4} = \frac{15}{64 \eta} \left( \frac{3\,058\,673}{1\,016\,064} + \frac{5429}{1008} \eta
 * + \frac{617}{144}
 * \eta^{2} \right).
 * \f}
 *
 * If we want to express \f$\Psi\f$ in terms of \f$f\f$ rather than \f$v\f$ we simply substitute \f$v = (\pi m
 * f)^{1/3}\f$
 * to obtain
 * \anchor phaselabel \f{equation}{
 * \Psi(f) = 2 \pi f t_{c} + \sum_{k} \theta^{\prime}_{k} f^{(k-5)/3}
 * \tag{phaselabel}
 * \f}
 * where
 * \f{equation}{
 * \theta^{\prime}_{k} = (\pi m)^{(k-5)/3} \theta_{k}.
 * \f}
 *
 * We are now in a position to start calculating components of \f$\gamma_{\alpha \beta}\f$. We had
 * \f{equation}{
 * \psi_{j} \equiv \frac{\partial \Delta \Psi}{\partial \Delta \lambda^{j}}
 * \f}
 * where \f$\Psi\f$ is given by Eq.\eqref{phaselabel}. Therefore we may write
 * \f{equation}{
 * \Delta \Psi = \Delta \theta^{\prime}_{0} f^{-5/3} + \Delta \theta^{\prime}_{2} f^{-1} + \Delta
 * \theta^{\prime}_{3} f^{-2/3} + \Delta \theta^{\prime}_{4} f^{-1/3}
 * \f}
 * All we need to do now is specify the coordinates \f$\lambda^{j}\f$ with respect to which the
 * derivatives
 * will be taken. In general, the template placement algorithm works in \f$(\tau_{0},\tau_{3})\f$
 * coordinates. It is simplest for us to calculate the components of \f$\gamma_{\alpha \beta}\f$ in
 * the \f$(m,\eta)\f$ coordinate system and then perform a coordinate transformation to get the
 * components in
 * the \f$(\tau_{0},\tau_{3})\f$ system.
 * So, we first of all calculate the components of \f$\gamma_{\alpha \beta}\f$ in the \f$(m,\eta)\f$
 * system.
 *
 * This involves calculating the following:
 * \f{equation}{
 * \frac{\partial \Delta \Psi}{\partial \Delta m} = \frac{\Delta \theta^{\prime}_{0}}{\Delta m}
 * f^{-5/3} +
 * \frac{\Delta \theta^{\prime}_{2}}{\Delta m} f^{-1} - \frac{\Delta \theta^{\prime}_{3}}{\Delta m}
 * f^{-2/3} + \frac{\delta \theta^{\prime}_{4}}{\Delta m} f^{-1/3}
 * \f}
 * and
 * \f{equation}{
 * \frac{\partial \Delta \Psi}{\partial \Delta \eta} = \frac{\Delta \theta^{\prime}_{0}}{\Delta
 * \eta} f^{-5/3} +
 * \frac{\Delta \theta^{\prime}_{2}}{\Delta \eta} f^{-1} - \frac{\Delta \theta^{\prime}_{3}}{\Delta
 * \eta}
 * f^{-2/3} + \frac{\delta \theta^{\prime}_{4}}{\Delta \eta} f^{-1/3}
 * \f}
 * where all of the derivatives are easily calculable. This gives us the terms \f$\psi_{j}\f$ as a
 * power
 * series in \f$f\f$. These are then used in the formula
 * \f{equation}{
 * \gamma_{\alpha \beta} = \frac{1}{2} \left( \mathcal{J} [ \psi_{\alpha} \psi_{\beta} ] -
 * \mathcal{J} [
 * \psi_{\alpha}] \mathcal{J} [\psi_{\beta}] \right)
 * \f}
 * to calculate the components of \f$\gamma_{\alpha \beta}\f$. The fact that each of the \f$\psi_{j}\f$ is
 * in the
 * form of a power series in \f$f\f$ allows us to calculate \f$\gamma_{\alpha \beta}\f$ using
 * \f{equation}{
 * \mathcal{J} \left[ \sum_{n} a_{n} x^{n} \right] = \sum_{n} J(7-3n).
 * \f}
 * i.e.\ we can express all the \f$\mathcal{J}[]\f$ in terms of the integral \f$J(q)\f$ which we calculate
 * numerically at the outset for the required values of \f$q\f$.
 *
 * Once we have obtained \f$\gamma_{\alpha \beta}\f$ in this way, we take the inverse of this matrix to
 * give us \f$\gamma^{\alpha \beta}\f$ in the \f$(m,\eta)\f$ system. Then
 * we perform the following coordinate transformation to give us the components of
 * \f$\gamma^{\alpha^{\prime}
 * \beta^{\prime}}\f$ in our chosen system,
 * \f{equation}{
 * \gamma^{\alpha^{\prime} \beta^{\prime}} = \Lambda^{\alpha^{\prime}}_{\,\,\sigma}
 * \Lambda^{\beta^{\prime}}_{\,\,\delta} \gamma^{\sigma \delta}
 * \f}
 * where the transformation matrix \f$\Lambda^{\alpha^{\prime}}_{\,\,\beta}\f$ is defined by
 * \f{equation}{
 * \Lambda^{\alpha^{\prime}}_{\,\,\beta} = \frac{\partial x^{\alpha^{\prime}}}{\partial x^{\beta}}
 * \f}
 * Finally, we take the inverse of this matrix to obtain \f$\gamma_{\alpha^{\prime} \beta^{\prime}}\f$
 * in the
 * chosen system.
 * Since the unprimed system corresponds to \f$(t_{c},m,\eta)\f$ coordinates and the primed system to
 * \f$(t_{c},\tau_{0},\tau_{3})\f$ coordinates, the matrix
 * \f$\Lambda^{\alpha^{\prime}}_{\,\,\beta^{\prime}}\f$ has
 * element
 * \f{equation}{
 * \Lambda^{\alpha^{\prime}}_{\,\,\beta} = \left(
 * \begin{array}{ccc}
 * 1  &  0  &  0  \\
 * 0  & \frac{\partial \tau_{0}}{\partial m}  &  \frac{\partial \tau_{0}}{\partial \eta}  \\
 * 0  & \frac{\partial \tau_{3}}{\partial m}  &  \frac{\partial \tau_{3}}{\partial \eta}
 * \end{array}
 * \right) = \left(
 * \begin{array}{ccc}
 * 1  &  0  &  0  \\
 * 0  &  -\frac{5 \tau_{0}}{3m}  &  -\frac{\tau_{0}}{\eta}  \\
 * 0  &  -\frac{2 \tau_{3}}{3m}  &  -\frac{\tau_{3}}{\eta}
 * \end{array}
 * \right)
 * \f}
 *
 * Finally, what is needed in laying a lattice in the space of dynamical parameters
 * (also referred to as intrinsic parameters) is the metric with the kinematical
 * parameter (also called extrinsic parameter) projected out: In other words one defines
 * the 2-dimensional metric \f$g_{mn}\f$ by
 * \f{equation}{
 * g_{mn} = \gamma_{mn} - \frac{\gamma_{0m} \gamma_{0n}}{\gamma_{00}}.
 * \f}
 *
 * ### Metric computation in the \f$\tau_0-\tau_3\f$ space ###
 *
 * The metric cannot be directly computed in the \f$(\tau_0,\tau_2)\f$ space.
 * Therefore, in the previous Section we first computed the metric
 * in the \f$(m,\eta)\f$ space and then transformed to \f$(\tau_0,\tau_2)\f$ space.
 * The same method can also be used to find the metric in the \f$(\tau_0,\tau_3)\f$ space.
 * However, in \f$(\tau_0,\tau_3)\f$ space one can directly compute the
 * metric without recourse to \f$(m,\eta)\f$ coordinates. It is of interest to see
 * whether this yields the same results as the previous method.
 *
 * The starting point of our derivation is Eq.\ (3.7) of Owen and Sathyaprakash
 * (Phys. Rev. D 60, 022002, 1999) for the Fourier domain phase which we shall
 * write as:
 * \f{eqnarray}{
 * \Psi(f; \theta_1, \theta_2) & = & a_{01}\theta_1 v^{-5}
 * + \left [a_{21} \frac {\theta_1}{\theta_2} + a_{22} \left ( \theta_1 \theta_2^2 \right )^{1/3} \right ] v^{-3}
 * + a_{31} \theta_2 v^{-2} \nonumber \\
 * & + & \left [a_{41} \frac {\theta_1}{\theta_2^2} + a_{42} \left ( \frac {\theta_1}{\theta_2} \right )^{1/3}
 * + a_{43} \left ( \frac{\theta_2^4}{\theta_1} \right )^{1/3} \right ] v^{-1},
 * \f}
 * to 2nd post-Newtonain order.  Here \f$v=(f/f_0)^{1/3},\f$ \f$\theta_1\f$ and \f$\theta_2\f$ are
 * identical to the  \f$\theta^1\f$ and \f$\theta^2\f$ parameters
 * of Owen and Sathyaprakash defined in Eq.\ (3.3) there and the \f$a\f$ coefficients are given by:
 * \f{eqnarray}{
 * a_{01} = \frac{3}{5}, \ \ a_{21} = \frac{11\pi}{12}, \ \
 * a_{22} = \frac{743}{2016} \left ( \frac {25}{2\pi^2} \right )^{1/3}, \ \ a_{31} = -\frac{3}{2}, \nonumber \\
 * a_{41} = \frac {617}{384} \pi^2, \ \ a_{42} = \frac{5429}{5376} \left ( \frac{25 \pi}{2} \right )^{1/3},\ \
 * a_{43} = \frac {15293365}{10838016} \left ( \frac{5}{4\pi^4} \right )^{1/3}.
 * \f}
 * Differentials of the phase with respect to the coordinates \f$\theta_1\f$ and \f$\theta_2\f$ appear in the
 * metric which we write as:
 * \f{equation}{
 * \psi_m \equiv \frac{\partial \Psi}{\partial \theta_m} = \sum_0^N \Psi_{mk} v^{k-5}.
 * \f}
 * where \f$N\f$ is the post-Newtonian order up to which the phase is known, or the post-Newtonian
 * order at which the metric is desired.
 * Expansion coefficients \f$\Psi_{mn}\f$ can be considered be \f$(2\times N)\f$ matrix which to
 * second post-Newtonian order is given by:
 * \f{equation}{
 * \Psi =
 * \left [ \begin{matrix}
 * a_{01}
 * & 0
 * & {a_{21}}/{\theta_2} + ({a_{22}}/{3}) \left ( {\theta_2}/{\theta_1} \right )^{2/3}
 * & 0
 * & {a_{41}}/{\theta_2^2} + {a_{42}}/\left ({3 \left ( \theta_1^2\theta_2 \right )^{1/3} } \right )
 * - ({a_{43}}/{3}) \left ( {\theta_2}/{\theta_1} \right )^{4/3} \cr
 * 0
 * & 0
 * & - {a_{21}\theta_1}/{\theta_2^2} + (2 {a_{22}}/{3}) \left ( {\theta_1}/{\theta_2} \right )^{1/3}
 * & a_{31}
 * & - {2a_{41} \theta_1}/{\theta_2^3} - ({a_{42}}/{3}) \left ( {\theta_1}/{\theta_2^4} \right )^{1/3}
 * + ({4a_{43}}/{3}) \left ( {\theta_2}/{\theta_1} \right )^{1/3}
 * \end{matrix}
 * \right ].
 * \f}
 *
 * Using the definition of the
 * metric introduced earlier and projecting out the \f$t_c\f$ coordinate, one finds that
 * \f{eqnarray}{
 * g_{mn}  & = & \frac{1}{2}\sum_{k,l=0}^N \Psi_{mk} \Psi_{nl}
 * \biggl  [ J(17-k-l) - J(12-k) J(12-l) \biggr . \nonumber \\
 * & - & 	\biggl . \frac { \left ( J(9-k) - J(4)J(12-k) \right )
 * \left ( J(9-l) - J(4)J(12-l) \right )} {\left (J(1) - J(4)^2 \right)}
 * \biggr ]
 * \f}
 * where \f$J\f$'s are the moments introduced earlier.
 */
int
XLALInspiralComputeMetric (
    InspiralMetric     *metric,
    InspiralMomentsEtc *moments,
    REAL8 fLower,
    LALPNOrder order,
    REAL8 t0,
    REAL8 t3
    )

{
  static REAL8 Psi[METRIC_DIMENSION][METRIC_ORDER];
  static REAL8 g[METRIC_DIMENSION][METRIC_DIMENSION];

  REAL8 a, b, c, q;
  UINT4 PNorder, m, n;

  const REAL8 two_pi_flower = LAL_TWOPI * fLower;
  const REAL8 two_pi_flower_sq = two_pi_flower * two_pi_flower;

  /* check inputs */
  if (!metric || !moments){
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  };

  if (t0 <= 0.L || t3 <= 0.L){
    XLALPrintError("t0 and t3 must be positive");
    XLAL_ERROR(XLAL_EDOM);
  };

  if (fLower <= 0){
    XLALPrintError("fLower must be positive");
    XLAL_ERROR(XLAL_EDOM);
  }

  /* use the order of the waveform to compute the metric */
  /* summation below will be carried out up to PNorder   */
  if ( order != LAL_PNORDER_ONE &&
      order != LAL_PNORDER_ONE_POINT_FIVE &&
      order != LAL_PNORDER_TWO )
  {
    /* Let us force the order to be twoPN because that is the only order
     * available for the template bank anyway. */
    XLALPrintWarning("forcing LAL_PNORDER_TWO");
    PNorder = LAL_PNORDER_TWO;
  }
  else
  {
    PNorder = (UINT4) order;
  }


  /* Setting up \Psi_{mn} coefficients  */
  InspiralComputeMetricGetPsiCoefficients( Psi, moments, fLower, t0, t3 );

  for ( m = 0; m < METRIC_DIMENSION; m++ )
  {
    for ( n = m; n < METRIC_DIMENSION; n++ )
    {
      UINT4 k, l;
      g[m][n] = 0.L;

      for ( k = 0 ; k < PNorder; k++ )
      {
        for ( l = 0; l < PNorder; l++ )
        {
          g[m][n] += Psi[m][k] * Psi[n][l] * (
              moments->j[17-k-l] - moments->j[12-k] * moments->j[12-l]
              - ( moments->j[9-k] - moments->j[4] * moments->j[12-k] )
              * ( moments->j[9-l] - moments->j[4] * moments->j[12-l] )
              / ( moments->j[1]   - moments->j[4] * moments->j[4]    )
              );
        }
      }
      g[m][n] /= 2.;
      g[n][m] = g[m][n];
    }
  }

#if 0
  The minimum sampling rate for given MM is
    srate =
    2 * LAL_PI * f0 sqrt( (moments.j[1] - moments.j[4]*moments.j[4]) /
        (2.0*(1.-MM)));
#endif

  /* The calculation above gives the metric in coordinates   */
  /* (t0=2\pi f_0 \tau0, t3=2\pi f_0 \tau3). Re-scale metric */
  /* coefficients to get metric in (tau0, tau3) coordinates  */
  a = g[0][0] * two_pi_flower_sq;
  b = g[0][1] * two_pi_flower_sq;
  c = g[1][1] * two_pi_flower_sq;


  /* The metric in tau0-tau2,3 space. */
  metric->G00 = a;
  metric->G01 = b;
  metric->G11 = c;

  /* Diagonalize the metric. */
  q = sqrt( (a-c)*(a-c) + 4. * b*b );

  metric->g00 = 0.5 * (a + c - q);
  metric->g11 = 0.5 * (a + c + q);

  if ( a == c )
  {
    metric->theta = LAL_PI/2.;
  }
  else
  {
    /* metric->theta = 0.5 * atan(2.*b/(a-c));                  */
    /* We want to always measure the angle from the             */
    /* semi-major axis to the tau0 axis which is given by       */
    /* the following line as opposed to the line above          */
    metric->theta = atan( b / (metric->g00 - c) );
  }

  /* memset the metric->Gamma array to zero before populating them with correct
   * values. This prevents junk getting stored in unused fields */
  memset (metric->Gamma, 0, 10*sizeof(REAL4));

  /* Now we compute the 3d metric in tc,\tau_0,\tau_3 co-ordinates */
  /* We only need metric->Gamma[0,...,5].                          */
  metric->Gamma[0] = 0.5*two_pi_flower_sq*
          ( moments->j[1] - (moments->j[4]*moments->j[4]) );

  metric->Gamma[1] = 0.5*two_pi_flower_sq*
          ( Psi[0][0]*(moments->j[9] - (moments->j[4]*moments->j[12]) )
          + Psi[0][2]*(moments->j[7] - (moments->j[4]*moments->j[10]) )
          + Psi[0][4]*(moments->j[5] - (moments->j[4]*moments->j[8]) ));

  metric->Gamma[2] = 0.5*two_pi_flower_sq*
          ( Psi[1][2]*(moments->j[7] - (moments->j[4]*moments->j[10]) )
          + Psi[1][3]*(moments->j[6] - (moments->j[4]*moments->j[9])  )
          + Psi[1][4]*(moments->j[5] - (moments->j[4]*moments->j[8])  ));


  metric->Gamma[3] = metric->G00 + metric->Gamma[1]*metric->Gamma[1]/metric->Gamma[0];
  metric->Gamma[4] = metric->G01 + metric->Gamma[1]*metric->Gamma[2]/metric->Gamma[0];
  metric->Gamma[5] = metric->G11 + metric->Gamma[2]*metric->Gamma[2]/metric->Gamma[0];


  return XLAL_SUCCESS;
}

static void
InspiralComputeMetricGetPsiCoefficients (
    REAL8              Psi[METRIC_DIMENSION][METRIC_ORDER],
    InspiralMomentsEtc *moments,
    REAL8 fLower,
    REAL8 t0,
    REAL8 t3
    )
{
  REAL8 t1 = LAL_TWOPI * fLower * t0;
  REAL8 t2 = LAL_TWOPI * fLower * t3;

  Psi[0][0] = moments->a01;
  Psi[0][1] = 0.L;
  Psi[0][2] = moments->a21/t2 + moments->a22/3.L * cbrt(t2 * t2 / (t1 * t1));
  Psi[0][3] = 0.L;
  Psi[0][4] = moments->a41/(t2*t2) + moments->a42/(3.L* cbrt(t1*t1*t2))
    - moments->a43/3.L * t2 / t1 * cbrt(t2 / t1);

  Psi[1][0] = 0.L;
  Psi[1][1] = 0.L;
  Psi[1][2] = -moments->a21*t1/(t2*t2) + 2.L *
    moments->a22/3.L * cbrt(t1/t2);
  Psi[1][3] =  moments->a31;
  Psi[1][4] = - 2.L * moments->a41*t1 / (t2*t2*t2) -
    moments->a42/3.L * cbrt(t1/(t2*t2*t2*t2)) +
    4.L * moments->a43/3.L * cbrt(t2/t1);
}

/**
 * UNDOCUMENTED
 * \see See LALInspiralComputeMetric() for documentation
 */
void
LALInspiralComputeMetricBCV (
    LALStatus             *status,
    InspiralMetric        *metric,
    REAL8FrequencySeries  *psd,
    InspiralTemplate      *params
    )
{
  REAL8 g[METRIC_DIMENSION][METRIC_DIMENSION];
  InspiralMomentsEtcBCV moments;
  REAL8 num;
  REAL8 a, b, c, q;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  moments.alpha = params->alpha;
  moments.n0 = 5.L/3.L;
  moments.n15 = 2.L/3.L;

  LALGetInspiralMomentsBCV( status->statusPtr, &moments, psd, params );
  CHECKSTATUSPTR( status );

  num =  moments.M3[0][0] *moments.M3[1][1]
    - moments.M3[0][1] * moments.M3[1][0];

  g[0][0] =moments.M2[0][0]*(moments.M3[1][1]*moments.M2[0][0]
      -moments.M3[0][1]*moments.M2[0][1])
    +moments.M2[0][1]*(-moments.M3[0][1]*moments.M2[0][0]
        +moments.M3[0][0]*moments.M2[0][1]);
  g[0][0] /= num;


  g[1][1] =moments.M2[0][1]*(moments.M3[1][1]*moments.M2[0][1]
      -moments.M3[0][1]*moments.M2[1][1])
    +moments.M2[1][1]*(-moments.M3[0][1]*moments.M2[0][1]
        +moments.M3[0][0]*moments.M2[1][1]);
  g[1][1] /= num;


  g[0][1] = moments.M2[0][0]*(moments.M3[1][1]*moments.M2[0][1]
      -moments.M3[0][1]*moments.M2[1][1])
    +moments.M2[0][1]*(-moments.M3[0][1]*moments.M2[0][1]
        +moments.M3[0][0]*moments.M2[1][1]);
  g[0][1] /= num ;

  metric->G00 = .5 *(moments.M1[0][0] - g[0][0] );
  metric->G01 = .5 *(moments.M1[0][1] - g[0][1] );
  metric->G11 = .5 *(moments.M1[1][1] - g[1][1] );

  a = metric->G00;
  b = metric->G01;
  c = metric->G11;
  q = sqrt( (a-c)*(a-c) + 4. * b*b );
  metric->g00 = 0.5 * (a + c - q);
  metric->g11 = 0.5 * (a + c + q);
  if ( a == c )
  {
    metric->theta = LAL_PI/2.;
  }
  else
  {
    /* metric->theta = 0.5 * atan(2.*b/(a-c));                  */
    /* We want to always measure the angle from the             */
    /* semi-major axis to the tau0 axis which is given by       */
    /* the following line as opposed to the line above          */
    metric->theta = atan(b/(metric->g00 - c));
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

#undef METRIC_ORDER
#undef METRIC_DIMENSION
