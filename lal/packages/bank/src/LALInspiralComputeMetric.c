/*  <lalVerbatim file="LALInspiralComputeMetricCV">
Author: Churches, D. K., Cokelaer, T., Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralComputeMetric.c}}

Module to compute the components of the metric which is used 
to describe distances on the signal manifold.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralComputeMetricCP}
\idx{LALInspiralComputeMetric()}
\begin{itemize}
   \item \texttt{metric,} Output, the metric at the lattice point defined by \texttt{params}
   \item \texttt{params,} Input, the parameters where metric must be computed
   \item \texttt{moments,} Input, moments $J(1), \ldots, J(17),$ of the PSD and other constants needed
   in the computation of the metric.
\end{itemize}

\subsubsection*{Description}
We calculate the components of the metric using the procedure outlined 
in Owen \cite{Owen:96}. 
This uses the moments of the noise curve,
\begin{equation}
I(q) \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{S_{h}(x)}
\, dx
\end{equation}
and
\begin{equation}
J(q) \equiv \frac{I(q)}{I(7)} \,.
\end{equation}
(Please note that the function {\tt LALInspiralMoments} doesn't compute $I(q)$ defined 
here; the index $q$ is definted differently there and the normalisation is supplied by the user.
For ease of writing, here we shall follow the standard notation.)
Then the moment functional $\mathcal{J}$ is defined such that, for a function $a$,
\begin{equation}
\mathcal{J} [a] \equiv \frac{1}{I(7)} \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-7/3}}{S_{h}(x)}
a(x) \, dx
\end{equation}
which gives us
\begin{equation}
\mathcal{J} = \left[ \sum_{n} a_{n} x^{n} \right] = \sum_{n} a_{n} J(7-3n).
\end{equation}
The above equation is used to calculate the components of the metric using the following formula:
\begin{equation}
\gamma_{\alpha \beta} = \frac{1}{2} \left( \mathcal{J} [ \psi_{\alpha} \psi_{\beta} ] -
\mathcal{J} [ \psi_{\alpha} ] \mathcal{J} [ \psi_{\beta} ] \right)
\end{equation}
where $\psi_\alpha$ is the derivative of the Fourier phase of the inspiral waveform
with respect to the parameter $\lambda^\alpha,$ that is
$\psi_\alpha \equiv \Psi_{,\alpha}.$  Writing the derivative index as
$\alpha=0,j,$ with $j=1,\ldots,n$ we have 
\begin{equation}
\psi_{0} \equiv 2 \pi f, \ \ \ 
\psi_{j} \equiv \frac{\partial \Delta \Psi}{\partial \Delta \lambda^{j}}.
\end{equation}
The phase $\Psi$ is that which appears in the usual stationary-phase formula for
the Fourier transform:
\begin{equation}
\tilde{h}(f)  \propto f^{-7/6} e^{i[-\pi/4 - \Phi_{0} + 2 \pi f t_{0} + \Psi(f;\vec{\lambda}]}.
\end{equation}
If we take the usual chirp times and multiply each by $(\pi f_{0})$ then we get
\emph{dimensionless
chirp times} $\tau_{k}$, and then the phase $\Psi$ may be written in the form
\begin{equation}
\Psi = 2 \pi f t_{c} + \sum_{k} \Psi_{k}(f) \tau_{k}
\end{equation}
where, defining $v_0 = (\pi m f_0)^{1/3}$ ($m$ being total mass and $f_0$ a
fiducial starting frequency), the chirptimes $\tau_{k},$ up to 2nd PN order,
are given by
$$
\tau_{0} = \frac{5}{256 \eta v_{0}^{5}},\ \ 
\tau_{2} = \frac{5}{192 \eta v_{0}^{3}} \left( \frac{743}{336} + \frac{11}{4} \eta \right),
$$
\begin{equation}
\tau_{3} = \frac{\pi}{8 \eta v_{0}^{2}},\ \ 
\tau_{4} = \frac{5}{128 \eta v_{0}} \left( \frac{3\,058\,673}{1\,016\,064} + \frac{5429}{1008}
\eta +
\frac{617}{144} \eta^{2} \right).
\end{equation}
Up to second post-Newtonian approximation the $\psi_{k}$ are given by
\begin{equation}
\psi_{0} = \frac{6}{5 \nu^{5/3}},\ \ 
\psi_{2} = \frac{2}{\nu},\ \ 
\psi_{3} = - \frac{3}{\nu^{2/3}},\ \ 
\psi_{4} = \frac{6}{\nu^{1/3}}.
\end{equation}
where $\nu = f/f_{0}$.

If we now make the substitution $f = v^{3}/\pi m$ we then the find that the phase may be
expressed in a simpler form
\begin{equation}
\Psi(v) = 2 \pi f t_{c} + \sum_{k} \theta_{k} v^{k-5}
\end{equation}
where the \emph{chirp parameters} $\theta_{k}$ are given by
$$
\theta_{0} = \frac{3}{128 \eta}, \ \ 
\theta_{2} = \frac{5}{96 \eta} \left( \frac{743}{336} + \frac{11}{4} \eta \right),
$$
\begin{equation}
\theta_{3} = - \frac{3 \pi}{8 \eta},\ \ 
\theta_{4} = \frac{15}{64 \eta} \left( \frac{3\,058\,673}{1\,016\,064} + \frac{5429}{1008} \eta
+ \frac{617}{144}
\eta^{2} \right).
\end{equation}
 
If we want to express $\Psi$ in terms of $f$ rather than $v$ we simply substitute $v = (\pi m
f)^{1/3}$
to obtain
\begin{equation}
\Psi(f) = 2 \pi f t_{c} + \sum_{k} \theta^{\prime}_{k} f^{(k-5)/3}
\label{phaselabel}
\end{equation}
where
\begin{equation}
\theta^{\prime}_{k} = (\pi m)^{(k-5)/3} \theta_{k}.
\end{equation}
 
We are now in a position to start calculating components of $\gamma_{\alpha \beta}$. We had
\begin{equation}
\psi_{j} \equiv \frac{\partial \Delta \Psi}{\partial \Delta \lambda^{j}}
\end{equation}
where $\Psi$ is given by Eq.(\ref{phaselabel}). Therefore we may write
\begin{equation}
\Delta \Psi = \Delta \theta^{\prime}_{0} f^{-5/3} + \Delta \theta^{\prime}_{2} f^{-1} + \Delta
\theta^{\prime}_{3} f^{-2/3} + \Delta \theta^{\prime}_{4} f^{-1/3}
\end{equation}
All we need to do now is specify the coordinates $\lambda^{j}$ with respect to which the
derivatives
will be taken. In general, the template placement algorithm works in $(\tau_{0},\tau_{3})$
coordinates. It is simplest for us to calculate the components of $\gamma_{\alpha \beta}$ in
the $(m,\eta)$ coordinate system and then perform a coordinate transformation to get the
components in
the $(\tau_{0},\tau_{3})$ system.
So, we first of all calculate the components of $\gamma_{\alpha \beta}$ in the $(m,\eta)$
system.
 
This involves calculating the following:
\begin{equation}
\frac{\partial \Delta \Psi}{\partial \Delta m} = \frac{\Delta \theta^{\prime}_{0}}{\Delta m}
f^{-5/3} +
\frac{\Delta \theta^{\prime}_{2}}{\Delta m} f^{-1} - \frac{\Delta \theta^{\prime}_{3}}{\Delta m}
f^{-2/3} + \frac{\delta \theta^{\prime}_{4}}{\Delta m} f^{-1/3}
\end{equation}
and
\begin{equation}
\frac{\partial \Delta \Psi}{\partial \Delta \eta} = \frac{\Delta \theta^{\prime}_{0}}{\Delta
\eta} f^{-5/3} +
\frac{\Delta \theta^{\prime}_{2}}{\Delta \eta} f^{-1} - \frac{\Delta \theta^{\prime}_{3}}{\Delta
\eta}
f^{-2/3} + \frac{\delta \theta^{\prime}_{4}}{\Delta \eta} f^{-1/3}
\end{equation}
where all of the derivatives are easily calculable. This gives us the terms $\psi_{j}$ as a
power
series in $f$. These are then used in the formula
\begin{equation}
\gamma_{\alpha \beta} = \frac{1}{2} \left( \mathcal{J} [ \psi_{\alpha} \psi_{\beta} ] -
\mathcal{J} [
\psi_{\alpha}] \mathcal{J} [\psi_{\beta}] \right)
\end{equation}
to calculate the components of $\gamma_{\alpha \beta}$. The fact that each of the $\psi_{j}$ is
in the
form of a power series in $f$ allows us to calculate $\gamma_{\alpha \beta}$ using
\begin{equation}
\mathcal{J} \left[ \sum_{n} a_{n} x^{n} \right] = \sum_{n} J(7-3n).
\end{equation}
i.e.\ we can express all the $\mathcal{J}[]$ in terms of the integral $J(q)$ which we calculate
numerically at the outset for the required values of $q$.
 
 
Once we have obtained $\gamma_{\alpha \beta}$ in this way, we take the inverse of this matrix to
give us $\gamma^{\alpha \beta}$ in the $(m,\eta)$ system. Then
we perform the following coordinate transformation to give us the components of
$\gamma^{\alpha^{\prime}
\beta^{\prime}}$ in our chosen system,
\begin{equation}
\gamma^{\alpha^{\prime} \beta^{\prime}} = \Lambda^{\alpha^{\prime}}_{\,\,\sigma}
\Lambda^{\beta^{\prime}}_{\,\,\delta} \gamma^{\sigma \delta}
\end{equation}
where the transformation matrix $\Lambda^{\alpha^{\prime}}_{\,\,\beta}$ is defined by
\begin{equation}
\Lambda^{\alpha^{\prime}}_{\,\,\beta} = \frac{\partial x^{\alpha^{\prime}}}{\partial x^{\beta}}
\end{equation}
Finally, we take the inverse of this matrix to obtain $\gamma_{\alpha^{\prime} \beta^{\prime}}$
in the
chosen system.
Since the unprimed system corresponds to $(t_{c},m,\eta)$ coordinates and the primed system to
$(t_{c},\tau_{0},\tau_{3})$ coordinates, the matrix
$\Lambda^{\alpha^{\prime}}_{\,\,\beta^{\prime}}$ has
element
\begin{equation}
\Lambda^{\alpha^{\prime}}_{\,\,\beta} = \left(
\begin{array}{ccc}
1  &  0  &  0  \\
0  & \frac{\partial \tau_{0}}{\partial m}  &  \frac{\partial \tau_{0}}{\partial \eta}  \\
0  & \frac{\partial \tau_{3}}{\partial m}  &  \frac{\partial \tau_{3}}{\partial \eta}
\end{array}
\right) = \left(
\begin{array}{ccc}
1  &  0  &  0  \\
0  &  -\frac{5 \tau_{0}}{3m}  &  -\frac{\tau_{0}}{\eta}  \\
0  &  -\frac{2 \tau_{3}}{3m}  &  -\frac{\tau_{3}}{\eta}
\end{array}
\right)
\end{equation}

Finally, what is needed in laying a lattice in the space of dynamical parameters
(also referred to as intrinsic parameters) is the metric with the kinematical
parameter (also called extrinsic parameter) projected out: In other words one defines
the 2-dimensional metric $g_{mn}$ by
\begin{equation}
g_{mn} = \gamma_{mn} - \frac{\gamma_{0m} \gamma_{0n}}{\gamma_{00}}.
\end{equation}

\subsubsection*{Metric computation in the $\tau_0-\tau_3$ space}

The metric cannot be directly computed in the $(\tau_0,\tau_2)$ space.
Therefore, in the previous Section we first computed the metric
in the $(m,\eta)$ space and then transformed to $(\tau_0,\tau_2)$ space.
The same method was also be used to find the metric in the $(\tau_0,\tau_3)$ space.
However, in $(\tau_0,\tau_3)$ space one can directly compute the
metric without recourse to $(m,\eta)$ coordinates. It is of interest to see
whether this yields the same results as the previous method.

The starting point of our derivation is Eq.~(3.7) of Owen and Sathyaprakash
(Phys. Rev. D 60, 022002, 1999) for the Fourier domain phase which we shall
write as:

\begin{eqnarray}
\psi(f; \theta_1, \theta_2) & = & a_{01}\theta_1 v^{-5}  
+ \left [a_{21} \frac {\theta_1}{\theta_2} + a_{22} \left ( \theta_1 \theta_2^2 \right )^{1/3} \right ] v^{-3}
+ a_{31} \theta_2 v^{-2} \nonumber \\
& + & \left [a_{41} \frac {\theta_1}{\theta_2^2} + a_{42} \left ( \frac {\theta_1}{\theta_2} \right )^{1/3} 
+ a_{43} \left ( \frac{\theta_2^4}{\theta_1} \right )^{1/3} \right ] v^{-1},
\end{eqnarray}
to 2nd post-Newtonain order.  Here $v=(f/f_0)^{1/3},$ $\theta_1$ and $\theta_2$ are 
identical to the  $\theta^1$ and $\theta^2$ parameters
of Owen and Sathyaprakash defined in Eq.~(3.3) there and the $a$ coefficients are given by:
\begin{eqnarray}
a_{01} = \frac{3}{5}, \ \ a_{21} = \frac{11\pi}{12}, \ \ 
a_{22} = \frac{743}{2016} \left ( \frac {25}{2\pi^2} \right )^{1/3}, \ \ a_{31} = -\frac{3}{2}, \nonumber \\
a_{41} = \frac {617}{384} \pi^2, \ \ a_{42} = \frac{5429}{5376} \left ( \frac{25 \pi}{2} \right )^{1/3},\ \ 
a_{43} = \frac {15293365}{10838016} \left ( \frac{5}{4\pi^4} \right )^{1/3}.
\end{eqnarray}
Differentials of the phase with respect to the coordinates $\theta_1$ and $\theta_2$ appear in the
metric which we write as:
\begin{equation}
\psi_m = \frac{\partial \psi}{\partial \theta_m} = \sum_0^N \Psi_{mk} v^{k-5}.
\end{equation}
where $N$ is the post-Newtonian order up to which the phase is known, or the post-Newtonian
order at which the metric is desired.
Expansion coefficients $\Psi_{mn}$ can be considered be $(2\times N)$ matrix which to
second post-Newtonian order is given by: 
\begin{equation}
\Psi = 
\left [ \matrix { 
	  a_{01} 
	& 0 
	& {a_{21}}/{\theta_2} + ({a_{22}}/{3}) \left ( {\theta_2}/{\theta_1} \right )^{2/3} 
	& 0 
	& {a_{41}}/{\theta_2^2} + {a_{42}}/\left ({3 \left ( \theta_1^2\theta_2 \right )^{1/3} } \right ) 
	- ({a_{43}}/{3}) \left ( {\theta_2}/{\theta_1} \right )^{4/3} \cr 
	  0
	& 0 
	& - {a_{21}\theta_1}/{\theta_2^2} + (2 {a_{22}}/{3}) \left ( {\theta_1}/{\theta_2} \right )^{1/3} 
	& a_{31} 
	& - {2a_{41} \theta_1}/{\theta_2^3} - ({a_{42}}/{3}) \left ( {\theta_1}/{\theta_2^4} \right )^{1/3}  
	+ ({4a_{43}}/{3}) \left ( {\theta_2}/{\theta_1} \right )^{1/3} }
\right ].
\end{equation}

Using the definition of the
metric introduced earlier and projecting out the $t_c$ coordinate, one finds that
\begin{eqnarray}
g_{mn}  & = & \frac{1}{2}\sum_{k,l=0}^N \Psi_{mk} \Psi_{nl} 
\biggl  [ J(17-k-l) - J(12-k) J(12-l) \biggr . \nonumber \\
	& - & 	\biggl . \frac { \left ( J(9-k) - J(4)J(12-k) \right )
		\left ( J(9-l) - J(4)J(12-l) \right )} {\left (J(1) - J(4)^2 \right)}
\biggr ]
\end{eqnarray}
where $J$'s are the moments introduced earlier. 


\subsubsection*{Algorithm}
 
 
\subsubsection*{Uses}
\begin{verbatim}
LALMAlloc
LALInspiralMoments
LALInverse3
LALMatrixTransform
LALFree
\end{verbatim}
 
\subsubsection*{Notes}
 
\vfill{\footnotesize\input{LALInspiralComputeMetricCV}}
 
</lalLaTeX>  */

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
	    det: Determinant of the metric. (not in versions after 02/05)
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
    InspiralTemplate   *params, 
    InspiralMomentsEtc *moments
    );

NRCSID(LALINSPIRALCOMPUTEMETRICC, "$Id$");

/* <lalVerbatim file="LALInspiralComputeMetricCP">  */
void 
LALInspiralComputeMetric (
    LALStatus          *status,
    InspiralMetric     *metric,
    InspiralTemplate   *params,
    InspiralMomentsEtc *moments
    )
/* </lalVerbatim> */
{
  static REAL8 Psi[METRIC_DIMENSION][METRIC_ORDER];
  static REAL8 g[METRIC_DIMENSION][METRIC_DIMENSION];

  REAL8 a, b, c, q, det;
  UINT4 PNorder, m, n;

  INITSTATUS(status, "LALInspiralComputeMetric", LALINSPIRALCOMPUTEMETRICC );

  ASSERT( metric, status, 
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( params, status, 
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( moments, status, 
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( params->t0 > 0.L, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT( params->t3 > 0.L, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  /* use the order of the waveform to compute the metric */
  /* summation below will be carried out up to PNorder   */
  if ( params->order != onePN && 
      params->order != onePointFivePN &&
      params->order != twoPN )
  {
    ABORT( status, LALINSPIRALBANKH_EORDER, LALINSPIRALBANKH_MSGEORDER );
  }
  PNorder = (UINT4) params->order;

  /* Setting up \Psi_{mn} coefficients  */
  InspiralComputeMetricGetPsiCoefficients( Psi, params, moments );

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
  a = g[0][0] * pow(2.*LAL_PI*params->fLower,2.); 
  b = g[0][1] * pow(2.*LAL_PI*params->fLower,2.); 
  c = g[1][1] * pow(2.*LAL_PI*params->fLower,2.); 


  /* The metric in tau0-tau2,3 space. */
  metric->G00 = a;
  metric->G01 = b;
  metric->G11 = c;

  /* Diagonalize the metric. */
  det = a * c - b * b;
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

  RETURN( status );
}

static void 
InspiralComputeMetricGetPsiCoefficients (
    REAL8              Psi[METRIC_DIMENSION][METRIC_ORDER], 
    InspiralTemplate   *params, 
    InspiralMomentsEtc *moments
    )
{
  REAL8 t1 = 2.L * LAL_PI * params->fLower * params->t0; 
  REAL8 t2 = 2.L * LAL_PI * params->fLower * params->t3; 

  Psi[0][0] = moments->a01;
  Psi[0][1] = 0.L;
  Psi[0][2] = moments->a21/t2 + moments->a22/3.L * pow(t2/t1,2.L/3.L);
  Psi[0][3] = 0.L;
  Psi[0][4] = moments->a41/(t2*t2) + moments->a42/(3.L* pow(t1*t1*t2,1.L/3.L)) 
    - moments->a43/3.L * pow(t2/t1,4.L/3.L);

  Psi[1][0] = 0.L;
  Psi[1][1] = 0.L;
  Psi[1][2] = -moments->a21*t1/pow(t2,2.L) + 2.L * 
    moments->a22/3.L * pow(t1/t2,1.L/3.L);
  Psi[1][3] =  moments->a31;
  Psi[1][4] = - 2.L * moments->a41*t1 / pow(t2,3.L) - 
    moments->a42/3.L * pow(t1/pow(t2,4.L),1.L/3.L) + 
    4.L * moments->a43/3.L * pow(t2/t1,1.L/3.L);
}

/* <lalVerbatim file="LALInspiralComputeMetricCP">  */
void 
LALInspiralComputeMetricBCV (
    LALStatus             *status,
    InspiralMetric        *metric,
    REAL8FrequencySeries  *psd,
    InspiralTemplate      *params
    )
{ /* </lalVerbatim> */
  REAL8 g[METRIC_DIMENSION][METRIC_DIMENSION];
  InspiralMomentsEtcBCV moments;
  REAL8 num;
  REAL8 a, b, c, q, det;

  INITSTATUS( status, 
      "LALInspiralComputeMetricBCV ", LALINSPIRALCOMPUTEMETRICC );
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
  det = a * c - b * b;
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
