/*  <lalVerbatim file="LALInspiralComputeMetricCV">
Author: Churches, D. K., Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralComputeMetric.c}}

Module to compute the components of the metric which is used 
to describe distances on the signal manifold.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralComputeMetricCP}
\index{\verb&LALInspiralComputeMetric()&}

\subsubsection*{Description}

We calculate the components of the metric using the procedure outlined 
in Owen \cite{Owen:96}. 
This uses the moments of the noise curve,
\begin{equation}
I(q) \equiv s_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{S_{h}(xf_{0})}
\, dx
\end{equation}
and
\begin{equation}
J(q) \equiv \frac{I(q)}{I(7)} \,.
\end{equation}
Then the moment functional $\mathcal{J}$ is defined such that, for a function $a$,
\begin{equation}
\mathcal{J} [a] = \frac{1}{I(7)} \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-7/3}}{S_{h}(xf_{0}}
a(x) \, dx
\end{equation}
which gives us
\begin{equation}
\mathcal{J} \left[ \sum_{n} a_{n} x^{n} \right] = \sum_{n} a_{n} J(7-3n)
\end{equation}
This is used to calculate the components of the metric using the following formula:
\begin{equation}
\gamma_{\alpha \beta} = \frac{1}{2} \left( \mathcal{J} [ \psi_{\alpha} \psi_{\beta} ] -
\mathcal{J} [ \psi_{\alpha} ] \mathcal{J} [ \psi_{\beta} ] \right)
\end{equation}
where we have
\begin{equation}
\psi_{0} \equiv 2 \pi f
\end{equation}
and
\begin{equation}
\psi_{j} \equiv \frac{\partial \Delta \Psi}{\partial \Delta \lambda^{j}}
\end{equation}
where
$\Psi$ is the phase of the waveform which appears in the usual stationary phase formula:
\begin{equation}
\tilde{u}  = f^{-7/6} e^{i[-pi/4 - \Phi_{0} + 2 \pi f t_{0} + \Psi(f;\vec{\lambda}]}
\end{equation}
If we take the usual chirp times and multiply each by $(\pi f_{0})$ then we get
\emph{dimensionless
chirp times} $\tau_{k}$, and then the phase $\Psi$ may be written in the form
 
\begin{equation}
\Psi = 2 \pi f t_{c} + \sum_{k} \Psi_{k}(f) \tau_{k}
\end{equation}
where the $\tau_{k}$ are given by
\begin{equation}
\tau_{0} = \frac{5}{256} \eta v_{0}^{5}
\end{equation}
\begin{equation}
\tau_{2} = \frac{5}{192} \eta v_{0}^{3} \left( \frac{743}{336} + \frac{11}{4} \eta \right)
\end{equation}
\begin{equation}
\tau_{3} = \frac{\pi}{8 \eta v_{0}^{2}}
\end{equation}
\begin{equation}
\tau_{4} = \frac{5}{128 \eta v_{0}} \left( \frac{3\,058\,673}{1\,016\,064} + \frac{5429}{1008}
\eta +
\frac{617}{144} \eta^{2} \right)
\end{equation}
and the $\Psi_{k}$ are given by
\begin{equation}
\Psi_{0} = \frac{6}{5 \nu^{5/3}}
\end{equation}
\begin{equation}
\Psi_{2} = \frac{2}{\nu}
\end{equation}
\begin{equation}
\Psi_{3} = - \frac{3}{\nu^{2/3}}
\end{equation}
\begin{equation}
\Psi_{4} = \frac{6}{\nu^{1/3}}
\end{equation}
where $\nu = f/f_{0}$.

If we now make the substitution $f = v^{3}/\pi m$ we then the find that the phase may be
expressed in the
simpler form
\begin{equation}
\Psi(v) = 2 \pi f t_{c} + \sum_{k} \theta_{k} v^{k-5}
\end{equation}
where the \emph{chirp parameters} $\theta_{k}$ are given by
\begin{equation}
\theta_{0} = \frac{3}{128 \eta}
\end{equation}
\begin{equation}
\theta_{2} = \frac{5}{96 \eta} \left( \frac{743}{336} + \frac{11}{4} \eta \right)
\end{equation}
\begin{equation}
\theta_{3} = - \frac{3 \pi}{8 \eta}
\end{equation}
\begin{equation}
\theta_{4} = \frac{15}{64 \eta} \left( \frac{3\,058\,673}{1\,016\,064} + \frac{5429}{1008} \eta
+ \frac{617}{144}
\eta^{2} \right)
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




#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

/*
*	Created: 7.9.96.
*	Author: B.S.Sathyaprakash, Caltech, Cardiff University.
*	Revision History: Updates 15.7.97; 31.8.97.; 18.3.99
*       First C version October 2000.
*	Purpose: To compute the metric and the template bank
*              parameters, corresponding to 2PN chirps.
*	Dependencies: moments.f, inverse.f transform.f
*	Outputs:
*	    det: Determinant of the metric.
*	    g00: 
*	    g11: 
*	  theta: Angle which the t0-axis makes with semi-major (dx0) axis.
*         srate: The minimal sampling rate required (computed but not outputted.
*	Notes: Owen and Sathyaprakash (Caltech collaboration notes).
*              Also Sathyaprakash, note from October 2000.
*/

NRCSID(LALINSPIRALCOMPUTEMETRICC, "$Id$");

/* <lalVerbatim file="LALInspiralComputeMetricCP">  */

void LALInspiralComputeMetric(LALStatus        *status,
                              InspiralMetric   *metric,
                              InspiralTemplate params,
                              INT4             pass)
{ /* </lalVerbatim> */

   INT4 i, Dim=3;
   REAL8 **trans=NULL, **mm3=NULL, **tm3=NULL, *dummy;
   REAL8 t_0, t_2, t_3, t_4, s0, s2, s3, s4, tm11, tm12, tm22, m2;
   REAL8 t_02,t22,t32,t42,s02,s22,s32,s42,eta2,k0sq,k1sq,k2sq;
   REAL8 k0, k1, k2, k00, k01, k02, k11, k12, k22, flso;
   static REAL8 i7,j1,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j17;
   static REAL8 a1, a2, a3, f1, f2, c4, c0, c2, c3, fr;
   REAL8 bsq, rtbsqm4ac, totmass, eta, det;
   InspiralMomentsIn in;

   INITSTATUS (status, "LALInspiralComputeMetric", LALINSPIRALCOMPUTEMETRICC);
   ATTATCHSTATUSPTR(status);

   ASSERT (metric,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (params.totalMass > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (params.eta > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (params.eta <= 0.25, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (params.fCutoff > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (params.t0 > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (params.t2 > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   totmass = params.totalMass * LAL_MTSUN_SI;
   eta = params.eta;
   flso = 1/(LAL_PI * totmass * pow(6.,1.5));
/* Allocating space for three, (Dim x Dim) matrices and then point the
   mm3, tm3 and trans arrays to dummy */
/* Arrays for the metric (mm3 -> m and eta, tm3->chirp times) and 
transformation   matrix */

   if (!(dummy = LALMalloc(sizeof(REAL8) * Dim * Dim * 3))) {
      ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
   }
   if (!(trans = (REAL8 **) LALMalloc(sizeof(REAL8*) * Dim * 3))) {
      LALFree(dummy);
      ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
   }

   mm3 = trans+Dim;
   tm3 = trans+2*Dim;
  
   for (i=0; i<Dim; i++) {
      trans[i] = &dummy[Dim*i]; 
      mm3[i] = &dummy[Dim*Dim + Dim*i]; 
      tm3[i] = &dummy[2*Dim*Dim + Dim*i]; 
   } 
   
/* If first time set all static variables and compute moments else go
   straight to the metric calculation */
   if (pass==1) {
/* The factor that converts chirp times from fLower to f0 and vice versa */
   fr = params.fLower;

   c0 = pow(fr,-8./3.);
   c2 = pow(fr,-6./3.);
   c3 = pow(fr,-5./3.);
   f1 = 2.*LAL_PI;
   f2 = f1*f1;
   c4 = 3.0586730/1.0160640;
   a1 = 924.0/743.0;
   a2 = 5.4290/(1.0080 * c4);
   a3 = 6.170/(1.440 * c4);

   in.NoisePsd = metric->NoisePsd;
   in.xmin = params.fLower;
   in.xmax = params.fCutoff;
   in.norm = 1.;

/* Normalised moments of the noise PSD from 1/3 to 17/3. */
   in.ndx = 7.00/3.0; LALInspiralMoments(status->statusPtr,&i7,in); CHECKSTATUSPTR(status);
   in.norm = i7;
   j7 = 1.0;
   in.ndx = 1.00/3.0; LALInspiralMoments(status->statusPtr,&j1,in);  CHECKSTATUSPTR(status);
   in.ndx = 4.00/3.0; LALInspiralMoments(status->statusPtr,&j4,in);  CHECKSTATUSPTR(status);
   in.ndx = 5.00/3.0; LALInspiralMoments(status->statusPtr,&j5,in);  CHECKSTATUSPTR(status);
   in.ndx = 6.00/3.0; LALInspiralMoments(status->statusPtr,&j6,in);  CHECKSTATUSPTR(status);
   in.ndx = 8.00/3.0; LALInspiralMoments(status->statusPtr,&j8,in);  CHECKSTATUSPTR(status);
   in.ndx = 9.00/3.0; LALInspiralMoments(status->statusPtr,&j9,in);  CHECKSTATUSPTR(status);
   in.ndx = 10.0/3.0; LALInspiralMoments(status->statusPtr,&j10,in); CHECKSTATUSPTR(status);
   in.ndx = 11.0/3.0; LALInspiralMoments(status->statusPtr,&j11,in); CHECKSTATUSPTR(status);
   in.ndx = 12.0/3.0; LALInspiralMoments(status->statusPtr,&j12,in); CHECKSTATUSPTR(status);
   in.ndx = 13.0/3.0; LALInspiralMoments(status->statusPtr,&j13,in); CHECKSTATUSPTR(status);
   in.ndx = 14.0/3.0; LALInspiralMoments(status->statusPtr,&j14,in); CHECKSTATUSPTR(status);
   in.ndx = 15.0/3.0; LALInspiralMoments(status->statusPtr,&j15,in); CHECKSTATUSPTR(status);
   in.ndx = 17.0/3.0; LALInspiralMoments(status->statusPtr,&j17,in); CHECKSTATUSPTR(status);
   }
   
/* Rescale the chirp times to begin from f0=1Hz */
   t_0 = params.t0 * pow(fr,8./3.);
   t_2 = params.t2 * pow(fr,6./3.);
   t_3 = params.t3 * pow(fr,5./3.);
   t_4 = params.t4 * pow(fr,4./3.);
   
/* Compute the various terms and factors involved in the definition of the metric */
   s0 = 0.6*t_0;
   s2 = t_2/(1.0+ a1*eta);
   s3 = 1.50*t_3;
   s4 = 3.0*t_4*(1.0- a3*eta*eta)/(1.0 + a2*eta + a3*eta*eta);

   k0 = f1*j4;
   k1 = -(f1/totmass)*(t_0*j12 + t_2*j10 - t_3*j9 + t_4*j8);
   k2 = -(f1/eta)*(s0*j12 + s2*j10 - s3*j9 + s4*j8);

   m2 = totmass*totmass;
   t_02 = t_0*t_0;
   t22 = t_2*t_2;
   t32 = t_3*t_3;
   t42 = t_4*t_4;
   s02 = s0*s0;
   s22 = s2*s2;
   s32 = s3*s3;
   s42 = s4*s4;
   eta2 = eta*eta;

   k00 = f2*j1;
   k01 = - (f2/totmass)*(t_0*j9 + t_2*j7 - t_3*j6 + t_4*j5);
   k02 = - (f2/eta)*(s0*j9 + s2*j7 - s3*j6 + s4*j5);
   k11 = (f2/m2)*(t_02*j17 + t22*j13 + t32*j11 + t42*j9 +
       2.*t_0*t_2*j15 - 2.*t_0*t_3*j14 + 2.*t_0*t_4*j13 -
       2.*t_2*t_3*j12 + 2.*t_2*t_4*j11 - 2.*t_3*t_4*j10);
   k22 = (f2/eta2)*(s02*j17 + s22*j13 + s32*j11 + s42*j9 +
      2.*s0*s2*j15 - 2.*s0*s3*j14 + 2.*s0*s4*j13 -
      2.*s2*s3*j12 + 2.*s2*s4*j11 - 2.*s3*s4*j10);
   k12 = (f2/eta/totmass)*(t_0*s0*j17 + t_2*s2*j13 + t_3*s3*j11 + t_4*s4*j9 +
       (t_0*s2 + s0*t_2)*j15 - (t_0*s3 + s0*t_3)*j14 + 
       (t_0*s4 + s0*t_4)*j13 - (t_2*s3 + s2*t_3)*j12 + 
       (t_2*s4 + s2*t_4)*j11 - (t_3*s4 + s3*t_4)*j10);
   
/*   if (lalDebugLevel==1) {
      fprintf(stderr, "%e %e %e %e\n", t_0, t_2, t_3, t_4);
      fprintf(stderr, "%e %e %e\n", k0, k1, k2);
      fprintf(stderr, "%e %e %e\n", k00, k01, k02);
      fprintf(stderr, "%e %e %e\n", k11, k22, k12);
   }
*/

   k0sq = k0*k0;
   k1sq = k1*k1;
   k2sq = k2*k2;
   mm3[0][0] = 0.50*(k00 - k0sq) ;
   mm3[0][1] = 0.50*(k01 - k0*k1);
   mm3[0][2] = 0.50*(k02 - k0*k2);
   mm3[1][0] = mm3 [0][1];
   mm3[1][1] = 0.50*(k11 - k1sq);
   mm3[1][2] = 0.50*(k12 - k1*k2);
   mm3[2][0] = mm3 [0][2];
   mm3[2][1] = mm3 [1][2];
   mm3[2][2] = 0.50*(k22 - k2sq);

/*   if (lalDebugLevel==1) {
      fprintf(stderr, "%e %e %e\n", mm3[0][0], mm3[0][1], mm3[0][2]);
      fprintf(stderr, "%e %e %e\n", mm3[1][0], mm3[1][1], mm3[1][2]);
      fprintf(stderr, "%e %e %e\n", mm3[2][0], mm3[2][1], mm3[2][2]);
   }
*/

   LALInverse3 (status->statusPtr, tm3, mm3);
   CHECKSTATUSPTR(status);

/*   if (lalDebugLevel==1) {
      fprintf(stderr, "%e %e %e\n", tm3[0][0], tm3[0][1], tm3[0][2]);
      fprintf(stderr, "%e %e %e\n", tm3[1][0], tm3[1][1], tm3[1][2]);
      fprintf(stderr, "%e %e %e\n", tm3[2][0], tm3[2][1], tm3[2][2]);
   }
*/

   trans[0][0] = 1.0;
   trans[0][1] = 0.0;
   trans[0][2] = 0.0;
   trans[1][0] = 0.0;
   trans[2][0] = 0.0;
   trans[1][1] = -5.0*t_0*c0/(3.0*totmass);
   trans[1][2] = -t_0*c0/eta;
   ASSERT (metric->space>=0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (metric->space<=1, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   switch (metric->space) {
      case Tau0Tau2:
         trans[2][1] = -t_2*c2/(totmass);
         trans[2][2] = -t_2*c2/(eta*(1.0+a1*eta));
         break;
      case Tau0Tau3:
         trans[2][1] = -2.0*t_3*c3/(3.0*totmass);
         trans[2][2] = -t_3*c3/eta;
         break;
   }


/*   if (lalDebugLevel==1) {
      fprintf(stderr, "%e %e %e\n", trans[0][0], trans[0][1], trans[0][2]);
      fprintf(stderr, "%e %e %e\n", trans[1][0], trans[1][1], trans[1][2]);
      fprintf(stderr, "%e %e %e\n", trans[2][0], trans[2][1], trans[2][2]);
   }
*/

   LALMatrixTransform (status->statusPtr, Dim, trans, tm3, mm3);
   CHECKSTATUSPTR(status);
   LALInverse3 (status->statusPtr, tm3, mm3);
   CHECKSTATUSPTR(status);

/*   if (lalDebugLevel==1) {
      fprintf(stderr, "%e %e %e\n", tm3[0][0], tm3[0][1], tm3[0][2]);
      fprintf(stderr, "%e %e %e\n", tm3[1][0], tm3[1][1], tm3[1][2]);
      fprintf(stderr, "%e %e %e\n", tm3[2][0], tm3[2][1], tm3[2][2]);
   }
*/

/*
   The minimum sampling rate for an MM=0.970 is
   srate = sqrt(tm3[0][0]/(2.0*(1.0-0.970)));
*/
     
   tm11 = tm3[1][1] - tm3[0][1] * tm3[0][1] / tm3[0][0];
   tm12 = tm3[1][2] - tm3[0][1] * tm3[0][2] / tm3[0][0];
   tm22 = tm3[2][2] - tm3[0][2] * tm3[0][2] / tm3[0][0];

/*   if (lalDebugLevel==1) {
      fprintf(stderr, "%e %e\n", tm11, tm12 );
      fprintf(stderr, "%e %e\n", tm12, tm22);
   }
*/
   
   det = tm11 * tm22 - tm12 * tm12;
   
   bsq = (tm11-tm22) * (tm11-tm22);
   rtbsqm4ac = sqrt(bsq + 4.0 * tm12 * tm12);
   metric->g00 = 0.50*(tm11+tm22 - rtbsqm4ac);
   metric->g11 = 0.50*(tm11+tm22 + rtbsqm4ac);
   metric->theta = 0.50 * atan(2.0*tm12/(tm11-tm22));


/*   if (lalDebugLevel==1) {
      fprintf(stderr, "det=%e g00=%e g11=%e theta=%e\n", det, metric->g00, metric->g11, metric->theta);
   }
*/
   LALFree(dummy); 
   LALFree(trans);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
