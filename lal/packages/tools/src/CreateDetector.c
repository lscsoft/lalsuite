/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Peter Shawhan, John Whelan
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

/************************************ <lalVerbatim file="CreateDetectorCV">
Author: J. T. Whelan <john.whelan@ligo.org>
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{CreateDetector.c}}
\label{tools:ss:CreateDetector.c}

Creates a \texttt{LALDetector} structure from a \texttt{LALFrDetector}
structure and the type of detector.

\subsubsection*{Prototypes}
\input{CreateDetectorCP}
\idx{XLALCreateDetector()}

\subsubsection*{Description}

This routine takes the site geometry described in the
\texttt{LALFrDetector} structure, along with a
\texttt{LALDetectorType} parameter, and constructs the Cartesian
detector location and response tensor needed to fill the
\texttt{LALDetector} output.

The detector type is needed because different types of detectors have
different response tensors.  In each case the response tensor is
determined by the unit vectors $\hat{u}_X$ and $\hat{u}_Y$
which are constant in an Earth-fixed rotating reference frame and
point in the ``X arm'' and ``Y arm'' directions, respectively; the
headings of these directions in a local frame at the detector are
stored in the \texttt{LALFrDetector} structure.

The detector types recognized are (all names are prefaced by
\texttt{LALDETECTORTYPE\_}):
\begin{itemize}
\item[\texttt{IFODIFF}] An interferometer in differential mode.  The
response tensor is given by $d^{ab}=\frac{1}{2} (u_X^au_X^b-u_Y^au_Y^b)$.
Note that this is the preferred form even in the two arms of the
detector are not perpendicular (e.g., at the GEO600 site).
\item[\texttt{IFOXARM}] An interferometer in one-armed mode with the
X arm active.  The response tensor is given by
$d^{ab}=\frac{1}{2}u_X^au_X^b$.
\item[\texttt{IFOYARM}] An interferometer in one-armed mode with the
Y arm active.  The response tensor is given by
$d^{ab}=\frac{1}{2}u_Y^au_Y^b$.
\item[\texttt{IFOCOMM}] An interferometer in common mode.  The
response tensor is given by $d^{ab}=\frac{1}{2} (u_X^au_X^b+u_Y^au_Y^b)$.
\item[\texttt{CYLBAR}] A cylindrical bar detector.  In this case the
``X arm'' is actually the symmetry axis of the bar, and the ``Y arm''
is ignored.  The response tensor is
$d^{ab}=u_X^au_X^b$.
\end{itemize}

In each of these cases, the basic transformation needed is to express a
unit vector $\hat{u}$ in terms of its
components in the Earth-fixed basis
$\{\hat{e}_1,\hat{e}_2,\hat{e}_3\}$.  The altitude angle ${\mathcal{A}}$ and
azimuth angle $\zeta$ allow us to express the unit vector  $\hat{u}$
corresponding to a direction in terms of an orthonormal basis consisting
of a vector $\hat{e}_{\scriptstyle{\rm E}}$ pointing due East within the
local horizontal, a vector $\hat{e}_{\scriptstyle{\rm N}}$ pointing due
North within the local horizontal, and an upward-pointing vector
$\hat{e}_{\scriptstyle{\rm U}}$ normal to the local horizontal
plane.\footnote{These form a right-handed basis, providing an answer to
the age-old question ``What's Up?'': ``East cross North.''}  The relationship
is
\begin{equation}
\hat{u} =   ( \hat{e}_{\scriptstyle{\rm E}}\sin\zeta
                           + \hat{e}_{\scriptstyle{\rm N}}\cos\zeta )
                            \cos{\mathcal{A}}
          + \hat{e}_{\scriptstyle{\rm U}} \sin{\mathcal{A}}
\ .
\end{equation}
Since the local horizontal is defined as the tangent plane to the
reference ellipsoid at the point with the detector's latitude $\beta$
and longitude $\lambda$, the local basis is related to the orthonormal
basis
$\{\hat{e}_\rho,\hat{e}_\lambda,\hat{e}_z\}$ of a cylindrical
co\"{o}rdinate system (related to the Earth-fixed Cartesian
co\"{o}rdinates by $x^1=\rho\cos\lambda$, $x^2=\rho\sin\lambda$, $x^3=z$,
so that $\hat{e}_\rho$ points away from the Earth's axis,
$\hat{e}_\lambda$ points in the direction of increasing longitude, and
$\hat{e}_z$ points in the direction of increasing $x^3$)
by
\begin{eqnarray}
\hat{e}_{\scriptstyle{\rm E}} &=& \hat{e}_\lambda \\
\hat{e}_{\scriptstyle{\rm N}} &=& - \hat{e}_\rho \sin\beta
                                  + \hat{e}_z \cos\beta \\
\hat{e}_{\scriptstyle{\rm U}} &=&   \hat{e}_\rho \cos\beta
                                  + \hat{e}_z \sin\beta
\end{eqnarray}
It is then straightforward to relate the cylindrical basis vectors to
those in the Earth-fixed Cartesian system by
\begin{eqnarray}
\hat{e}_\rho    &=&  \hat{e}_1\cos\lambda  +  \hat{e}_2\sin\lambda  \\
\hat{e}_\lambda &=& -\hat{e}_1\sin\lambda  +  \hat{e}_2\cos\lambda  \\
\hat{e}_z       &=& \hat{e}_3
\end{eqnarray}

To express $\hat{u}$ in the Cartesian basis, we need
$\hat{u}\cdot\hat{e}_1$, $\hat{u}\cdot\hat{e}_2$, and
$\hat{u}\cdot\hat{e}_3$.  We first observe that
\begin{eqnarray}
\label{tools:e:eE}
\hat{u}\cdot\hat{e}_{\scriptstyle{\rm E}} &=& \cos{\mathcal{A}}\,\sin\zeta \\
\hat{u}\cdot\hat{e}_{\scriptstyle{\rm N}} &=& \cos{\mathcal{A}}\,\cos\zeta \\
\hat{u}\cdot\hat{e}_{\scriptstyle{\rm U}} &=& \sin{\mathcal{A}}
\end{eqnarray}
then that
\begin{eqnarray}
\hat{u}\cdot\hat{e}_\rho &=& (\hat{u}\cdot\hat{e}_{\scriptstyle{\rm N}})
                             (\hat{e}_{\scriptstyle{\rm N}}\cdot\hat{e}_\rho)
                           + (\hat{u}\cdot\hat{e}_{\scriptstyle{\rm U}})
                             (\hat{e}_{\scriptstyle{\rm U}}\cdot\hat{e}_\rho)
= -(\hat{u}\cdot\hat{e}_{\scriptstyle{\rm N}}) \sin\beta
 +(\hat{u}\cdot\hat{e}_{\scriptstyle{\rm U}}) \cos\beta \\
\hat{u}\cdot\hat{e}_\lambda &=& \hat{u}\cdot\hat{e}_{\scriptstyle{\rm E}}\\
\hat{u}\cdot\hat{e}_z &=& (\hat{u}\cdot\hat{e}_{\scriptstyle{\rm N}})
                             (\hat{e}_{\scriptstyle{\rm N}}\cdot\hat{e}_z)
                           + (\hat{u}\cdot\hat{e}_{\scriptstyle{\rm U}})
                             (\hat{e}_{\scriptstyle{\rm U}}\cdot\hat{e}_z)
= (\hat{u}\cdot\hat{e}_{\scriptstyle{\rm N}}) \cos\beta
 +(\hat{u}\cdot\hat{e}_{\scriptstyle{\rm U}}) \sin\beta
\end{eqnarray}
and finally that
\begin{eqnarray}
\hat{u}\cdot\hat{e}_1 &=& (\hat{u}\cdot\hat{e}_\rho)
                             (\hat{e}_\rho\cdot\hat{e}_1)
                           + (\hat{u}\cdot\hat{e}_\lambda)
                             (\hat{e}_\lambda\cdot\hat{e}_1)
= (\hat{u}\cdot\hat{e}_\rho) \cos\lambda
 -(\hat{u}\cdot\hat{e}_\lambda) \sin\lambda \\
\hat{u}\cdot\hat{e}_2 &=& (\hat{u}\cdot\hat{e}_\rho)
                             (\hat{e}_\rho\cdot\hat{e}_2)
                           + (\hat{u}\cdot\hat{e}_\lambda)
                             (\hat{e}_\lambda\cdot\hat{e}_2)
= (\hat{u}\cdot\hat{e}_\rho) \sin\lambda
 +(\hat{u}\cdot\hat{e}_\lambda) \cos\lambda \\
\hat{u}\cdot\hat{e}_3 &=& \hat{u}\cdot\hat{e}_z
\label{tools:e:e3ez}
\end{eqnarray}

\subsubsection*{Cached Detectors}
\idx[Variable]{lalCachedDetectors[]}

To avoid repeatedly calculating the Cartesian co\"{o}rdinates and
response tensor of known detectors, the constant array
\texttt{lalCachedDetectors[]} contains the site geometry and
response tensors of the most commonly used detectors.  These are
defined in this file and listed in Table~\ref{tools:tab:cached}.
\input{CreateDetectorCT}

\subsubsection*{Algorithm}
\texttt{XLALCreateDetector()} first checks the
\texttt{lalCachedDetectors[]} array to see if the specified type and
the name in the input \texttt{LALFrDetector} match any of the
predefined constant detectors.  If so, it returns a copy of the
constant detector (not just a pointer to the constant).

If not, it calculates the Cartesian co\"{o}rdinates $\{x^1,x^2,x^3\}$
of the detector location defined by (\ref{tools:e:cart1}--\ref{tools:e:cart3});
 in
particular, it calculates the denominator
$\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}$ and the distance from the axis
\begin{equation}
\rho = \left(\frac{a^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}} + h \right)
\cos\beta
\end{equation}
as intermediate steps.

It then calculates the Cartesian components of the unit vectors
$\hat{u}_X$ and $\hat{u}_Y$ in the arm directions from the altitude
and azimuth angles by use of a \texttt{static} function which
implements (\ref{tools:e:eE}--\ref{tools:e:e3ez}).  (Depending on the detector
type specified, only the unit vector(s) which are actually needed are
calculated.)  Using this components it constructs $d^{ab}$ according
to the formula appropriate to the detector type.

The calculation of $x^a$ is done to double precision, that of $d^{ab}$
to single precision.

\subsubsection*{Uses}

\begin{verbatim}
XLALCreateDetector()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item The conventions in the \texttt{LALFrDetector} structure are based
  on version 6 of the frame specification \cite{tools:LIGOVIRGO:2000}.
\item If the location and response tensor information for a
\texttt{LALDetector} are filled in by hand (e.g., for testing
purposes), the \texttt{type} field should be set to
\texttt{LALDETECTORTYPE\_ABSENT}.
\item The range of \texttt{LALDetectorType}s could be expanded to
include the  monopole and five quadrupole modes for a spherical
resonant detector
\cite{tools:Maggiore:2000b,tools:Zhou:1995,tools:Bianchi:1998,tools:Maggiore:2000a}.
\end{itemize}

\vfill{\footnotesize\input{CreateDetectorCV}}

******************************************************* </lalLaTeX> */

/**************************************** <lalLaTeX file="CreateDetectorCB">
\bibitem{tools:Anderson:2000}
  W.~G.~Anderson, J.~T.~Whelan, P.~R.~Brady, J.~D.~E.~Creighton,
  D.~Chin, and K.~Riles, ``Beam Pattern Response Functions and Times
  of Arrival for Earthbound Interferometers'', Maple worksheet,
  pdf version at
  \href{http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/T010110.pdf}
  {\texttt{http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/T010110.pdf}}
\bibitem{tools:Maggiore:2000b}
  M.~Maggiore, ``Gravitational Wave Experiments and Early Universe Cosmology'',
  Phys.\ Rept.\ \textbf{331}, 283-367 (2000);
  \href{http://www.arXiv.org/abs/gr-qc/9909001}{gr-qc/9909001}
\bibitem{tools:Zhou:1995} C.~Z.~Zhou and P.~F.~Michelson, ``Spherical
  resonant-mass gravitational wave detectors'', Phys.\ Rev.\ D.\ {\bf
  51}, 2517-2545 (1995).
\bibitem{tools:Bianchi:1998} M.~Bianchi, M.~Brunetti, E.~Coccia, F.~Fucito,
 and J.~A.~Lobo,
``Cross section of a resonant-mass detector for scalar gravitational waves''
  Phys.\ Rev.\ D.\ {\bf 57}, 4525--4534 (1998);
  \href{http://www.arXiv.org/abs/gr-qc/9709045}{gr-qc/9709045}.
\bibitem{tools:Maggiore:2000a}
  M.~Maggiore and A.~Nicholis, ``Detection strategies for scalar
  gravitational waves with interferometers and resonant spheres'',
  Phys.\ Rev.\ D \textbf{62}, 024004 (2000);
  \href{http://www.arXiv.org/abs/gr-qc/9907055}{gr-qc/9907055}
******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <math.h>
#include <string.h>
#include <lal/DetectorSite.h>

NRCSID( CREATEDETECTORC, "$Id$" );


/**************************************** <lalLaTeX file="CreateDetectorCT">

\begin{table}[htbp]
  \begin{center}
    \begin{tabular}{|c|c|c|}
\hline
      index & \texttt{LAL\_LHO\_4K\_DETECTOR}
            & \texttt{LAL\_LLO\_4K\_DETECTOR}
\\ \hline
      prefix & \texttt{H1}
            & \texttt{L1}
\\ \hline
      $x^a$
      &
      $
      \left(
        \begin{array}{c}
        -2.1614149\times 10^6 \\
        -3.8346952\times 10^6 \\
        4.6003502\times 10^6
        \end{array}
      \right)
      $
      &
      $
      \left(
        \begin{array}{c}
        -74276.044 \\
        -5.496283721\times 10^6 \\
        3.224257018\times 10^6
        \end{array}
      \right)
      $
\\ \hline
      $d^{ab}$
      &
      $
      \left(
        \begin{array}{ccc}
          -0.3926141& -0.0776130& -0.2473886 \\
          -0.0776130&  0.3195244&  0.2279981 \\
          -0.2473886&  0.2279981&  0.0730903
        \end{array}
      \right)
      $
      &
      $
      \left(
        \begin{array}{ccc}
          0.4112809&  0.1402097&  0.2472943 \\
          0.1402097& -0.1090056& -0.1816157 \\
          0.2472943& -0.1816157& -0.3022755
        \end{array}
      \right)
      $
\\ \hline
     type & \texttt{LALDETECTORTYPE\_IFODIFF} & \texttt{LALDETECTORTYPE\_IFODIFF}
\\ \hline
     name & LHO\_4k & LLO\_4k
\\ \hline
     $(\lambda,\beta,h)$
     & $(-(119^\circ24'27''\!\!.5657),46^\circ27'18''\!\!.528,
         142.544\,\textrm{m})$
     & $(-(90^\circ46'27''\!\!.2654),30^\circ33'46\!\!.4196,
           -6.574\,\textrm{m})$
\\ \hline
      $({\mathcal{A}}_X,\zeta_X)$
      & $(         -6.195\times 10^{-4},      324^\circ\!\!.0006)$
      & $(          -3.121\times 10^{-4},     252^\circ\!\!.2835)$
\\ \hline
      $({\mathcal{A}}_Y,\zeta_Y)$
      & $(           1.25\times 10^{-5},      234^\circ\!\!.0006)$
      & $(          -6.107\times 10^{-4},     162^\circ\!\!.2835)$
\\ \hline
     $(L_X/2,L_Y/2)$
     & $(2000\,\textrm{m},2000\,\textrm{m})$
     & $(2000\,\textrm{m},2000\,\textrm{m})$
\\ \hline
\hline
      index & \texttt{LAL\_VIRGO\_DETECTOR}
            & \texttt{LAL\_GEO\_600\_DETECTOR}
\\ \hline
      $x^a$
      &
      $
      \left(
        \begin{array}{c}
         4.54637409863\times 10^6 \\
         842989.697467\\
         4.37857696275\times 10^6
        \end{array}
      \right)
      $
      &
      $
      \left(
        \begin{array}{c}
        3.85630994953\times 10^6 \\
        666598.956352 \\
        5.01964141692\times 10^6
        \end{array}
      \right)
      $
\\ \hline
      $d^{ab}$
      &
      $
      \left(
        \begin{array}{ccc}
        0.2438740 &  -0.0990838 & -0.2325762 \\
        -0.0990838 &  -0.4478258 &   0.1878331 \\
        -0.2325762 &   0.1878331 &   0.2039518
        \end{array}
      \right)
      $
      &
      $
      \left(
        \begin{array}{ccc}
        -0.0968250 &  -0.3657823 &    0.1221373 \\
        -0.3657823 &   0.2229681 &    0.2497174 \\
         0.1221373 &   0.2497174 &   -0.1261431
        \end{array}
      \right)
      $
\\ \hline
     type & \texttt{LALDETECTORTYPE\_IFODIFF} & \texttt{LALDETECTORTYPE\_IFODIFF}
\\ \hline
     name & VIRGO & GEO\_600
\\ \hline
     $(\lambda,\beta,h)$
     & $(10^\circ30'16''\!\!.1878,43^\circ37'\!\!53''.0921,
         51.884\,\textrm{m})$
     & $(9^\circ48'25''\!\!.894,52^\circ14'42''\!\!.528,
         114.425\,\textrm{m})$
\\ \hline
      $({\mathcal{A}}_X,\zeta_X)$
        & $( 0,          19^\circ\!\!.4326)$
        & $( 0,          68^\circ\!\!.3883)$
\\ \hline
      $({\mathcal{A}}_Y,\zeta_Y)$
      & $( 0,           289^\circ\!\!.4326)$
      & $( 0,           334^\circ\!\!.0569)$
\\ \hline
     $(L_X/2,L_Y/2)$
     & $(1500\,\textrm{m},1500\,\textrm{m})$
     & $(300\,\textrm{m},300\,\textrm{m})$
\\ \hline
\hline
      index & \texttt{LAL\_TAMA\_300\_DETECTOR}
            & \texttt{LAL\_CIT\_40\_DETECTOR}
\\ \hline
      $x^a$
      &
      $
      \left(
        \begin{array}{c}
        -3.94640898771\times 10^6 \\
         3.36625903242\times 10^6\\
         3.69915069189\times 10^6
        \end{array}
      \right)
      $
      &
      $
      \left(
        \begin{array}{c}
        -2.49064958399\times 10^6 \\
        -4.65869968229\times 10^6 \\
         3.56206411337\times 10^6
        \end{array}
      \right)
      $
\\ \hline
      $d^{ab}$
      &
      $
      \left(
        \begin{array}{ccc}
         0.1121397 & 0.3308421 & -0.1802193 \\
         0.3308421 & 0.2177940 &  0.1537258 \\
        -0.1802193 & 0.1537258 & -0.3299337
        \end{array}
      \right)
      $
      &
      $
      \left(
        \begin{array}{ccc}
        -0.3537959 &  0.2734713 &  0.1095458 \\
         0.2734713 &  0.0115214 &  0.2049027 \\
         0.1095458 &  0.2049027 &  0.3422745
        \end{array}
      \right)
      $
\\ \hline
     type & \texttt{LALDETECTORTYPE\_IFODIFF} & \texttt{LALDETECTORTYPE\_IFODIFF}
\\ \hline
     name & TAMA\_300 & CIT\_40
\\ \hline
     $(\lambda,\beta,h)$
     & $(139^\circ32'9''\!\!.8,35^\circ40'35''\!\!.6,
         90\,\textrm{m})$
     & $(-118^\circ\!\!.13,34^\circ\!\!.17,
         0\,\textrm{m})$
\\ \hline
      $({\mathcal{A}}_X,\zeta_X)$
        & $( 0,          270^\circ)$
        & $( 0,          180^\circ)$
\\ \hline
      $({\mathcal{A}}_Y,\zeta_Y)$
      & $( 0,         180^\circ)$
      & $(0,          90^\circ)$
\\ \hline
     $(L_X/2,L_Y/2)$
     & $(150\,\textrm{m},150\,\textrm{m})$
     & $(20\,\textrm{m},20\,\textrm{m})$
\\ \hline

    \end{tabular}
    \caption{Selected redefined gravitational wave detectors, contained in
      the \texttt{lalCachedDetectors[]} array.
      Not shown in the table are the LHO 2\,km detector (H2) and the bar
      detectors ALLEGRO, AURIGA, EXPLORER, NIOBE and NAUTILUS.
      The LIGO site data
      come directly from \cite{tools:Althouse:1999}, including the
      Cartesian position vectors $x^a$ and the response tensor
      $d^{ab}$, which was dermined from the quoted components of the
      detector frame basis vectors $\hat{x}_G\equiv\hat{u}_X$ and
      $\hat{y}_G\equiv\hat{u}_Y$.  The data on the other detectors
      comes from \cite{tools:Anderson:2000}.}
    \label{tools:tab:cached}
  \end{center}
\end{table}

******************************************************* </lalLaTeX> */

/*  { name,
      vertexLatitiudeRadians,
      vertexLongitudeRadians,
      vertexElevation,
      xArmAltitudeRadians, xArmAzimuthRadians,
      yArmAltitudeRadians, yArmAzimuthRadians }   */

/* New method for creating cached detectors:
 *
 * Construct the detector structures from the macros describing the
 * detectors.
 *
 * Use a bit of macro magic to do this.
 */

#define LAL_CAT(x,y) x ## y
#define LAL_XCAT(x,y) LAL_CAT(x,y)

/* expands to constant c of detector d */
#define LAL_DETECTOR_CONSTANT(d,c) LAL_XCAT(LAL_XCAT(LAL_,d),LAL_XCAT(_,c))

/* initializer for detector location vector */
#define LAL_DETECTOR_LOCATION(d) \
{ \
  LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_X_SI),\
  LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_Y_SI),\
  LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_Z_SI) \
}

/* expands to component c (X,Y,Z) of arm X of detector d */
#define LAL_ARM_X(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(ARM_X_DIRECTION_,c))

/* expands to component c (X,Y,Z) of arm Y of detector d */
#define LAL_ARM_Y(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(ARM_Y_DIRECTION_,c))

/* expands to component c (X,Y,Z) of axis of detector d */
#define LAL_AXIS(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(AXIS_DIRECTION_,c))

/* expands to a 3x3 matix initializer for the response for IFODIFF detector d */
#define LAL_DETECTOR_RESPONSE_IFODIFF(d) \
{ \
  { \
    0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,X) - LAL_ARM_Y(d,X) * LAL_ARM_Y(d,X) ), \
    0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Y) - LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Y) ), \
    0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Z) - LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Z) )  \
  }, \
  { \
    0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,X) - LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,X) ), \
    0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Y) - LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Y) ), \
    0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Z) - LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Z) )  \
  }, \
  { \
    0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,X) - LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,X) ), \
    0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Y) - LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Y) ), \
    0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Z) - LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Z) )  \
  } \
}

/* expands to a 3x3 matix initializer for the response for IFOCOMM detector d */
#define LAL_DETECTOR_RESPONSE_IFOCOMM(d) \
{ \
  { \
    0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,X) + LAL_ARM_Y(d,X) * LAL_ARM_Y(d,X) ), \
    0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Y) + LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Y) ), \
    0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Z) + LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Z) )  \
  }, \
  { \
    0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,X) + LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,X) ), \
    0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Y) + LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Y) ), \
    0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Z) + LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Z) )  \
  }, \
  { \
    0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,X) + LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,X) ), \
    0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Y) + LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Y) ), \
    0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Z) + LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Z) )  \
  } \
}

/* expands to a 3x3 matix initializer for the response for IFOXARM detector d */
#define LAL_DETECTOR_RESPONSE_IFOXARM(d) \
{ \
  { \
    0.5 * LAL_ARM_X(d,X) * LAL_ARM_X(d,X), \
    0.5 * LAL_ARM_X(d,X) * LAL_ARM_X(d,Y), \
    0.5 * LAL_ARM_X(d,X) * LAL_ARM_X(d,Z)  \
  }, \
  { \
    0.5 * LAL_ARM_X(d,Y) * LAL_ARM_X(d,X), \
    0.5 * LAL_ARM_X(d,Y) * LAL_ARM_X(d,Y), \
    0.5 * LAL_ARM_X(d,Y) * LAL_ARM_X(d,Z)  \
  }, \
  { \
    0.5 * LAL_ARM_X(d,Z) * LAL_ARM_X(d,X), \
    0.5 * LAL_ARM_X(d,Z) * LAL_ARM_X(d,Y), \
    0.5 * LAL_ARM_X(d,Z) * LAL_ARM_X(d,Z)  \
  } \
}

/* expands to a 3x3 matix initializer for the response for IFOYARM detector d */
#define LAL_DETECTOR_RESPONSE_IFOYARM(d) \
{ \
  { \
    0.5 * LAL_ARM_Y(d,X) * LAL_ARM_Y(d,X), \
    0.5 * LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Y), \
    0.5 * LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Z)  \
  }, \
  { \
    0.5 * LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,X), \
    0.5 * LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Y), \
    0.5 * LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Z)  \
  }, \
  { \
    0.5 * LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,X), \
    0.5 * LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Y), \
    0.5 * LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Z)  \
  } \
}

/* expands to a 3x3 matix initializer for the response for CYLBAR detector d */
#define LAL_DETECTOR_RESPONSE_CYLBAR(d) \
{ \
  { \
    LAL_AXIS(d,X) * LAL_AXIS(d,X), \
    LAL_AXIS(d,X) * LAL_AXIS(d,Y), \
    LAL_AXIS(d,X) * LAL_AXIS(d,Z)  \
  }, \
  { \
    LAL_AXIS(d,Y) * LAL_AXIS(d,X), \
    LAL_AXIS(d,Y) * LAL_AXIS(d,Y), \
    LAL_AXIS(d,Y) * LAL_AXIS(d,Z)  \
  }, \
  { \
    LAL_AXIS(d,Z) * LAL_AXIS(d,X), \
    LAL_AXIS(d,Z) * LAL_AXIS(d,Y), \
    LAL_AXIS(d,Z) * LAL_AXIS(d,Z)  \
  } \
}

#define LAL_FR_DETECTOR_STRUCT(d) \
{ \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_NAME), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_PREFIX), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_LONGITUDE_RAD), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_LATITUDE_RAD), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ELEVATION_SI), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_X_ALTITUDE_RAD), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_X_AZIMUTH_RAD), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_Y_ALTITUDE_RAD), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_Y_AZIMUTH_RAD), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_X_MIDPOINT_SI), \
  LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_Y_MIDPOINT_SI) \
}

#define LAL_DETECTOR_RESPONSE(d,t) \
  LAL_XCAT( LAL_DETECTOR_RESPONSE_, t )(d)

#define LAL_DETECTOR_STRUCT(d,t) \
{ \
  LAL_DETECTOR_LOCATION(d),      \
  LAL_DETECTOR_RESPONSE(d,t),    \
  LAL_XCAT(LALDETECTORTYPE_,t),  \
  LAL_FR_DETECTOR_STRUCT(d)      \
}

/** Pre-existing detectors. */
const LALDetector lalCachedDetectors[LAL_NUM_DETECTORS] = {
  LAL_DETECTOR_STRUCT( TAMA_300, IFODIFF ),
  LAL_DETECTOR_STRUCT( VIRGO, IFODIFF ),
  LAL_DETECTOR_STRUCT( GEO_600, IFODIFF ),
  LAL_DETECTOR_STRUCT( LHO_2K, IFODIFF ),
  LAL_DETECTOR_STRUCT( LHO_4K, IFODIFF ),
  LAL_DETECTOR_STRUCT( LLO_4K, IFODIFF ),
  LAL_DETECTOR_STRUCT( CIT_40, IFODIFF ),
  LAL_DETECTOR_STRUCT( ALLEGRO_320, CYLBAR ),
  LAL_DETECTOR_STRUCT( AURIGA, CYLBAR ),
  LAL_DETECTOR_STRUCT( EXPLORER, CYLBAR ),
  LAL_DETECTOR_STRUCT( NIOBE, CYLBAR ),
  LAL_DETECTOR_STRUCT( NAUTILUS, CYLBAR )
};


static
void getCartesianComponents( REAL4 u[3],
                             REAL8 cosAlt, REAL8 sinAlt,
                             REAL8 cosAz,  REAL8 sinAz,
                             REAL8 cosLat, REAL8 sinLat,
                             REAL8 cosLon, REAL8 sinLon )
{
  REAL8 uNorth = cosAlt * cosAz;
  REAL8 uEast = cosAlt * sinAz;
  /* uUp == sinAlt */
  REAL8 uRho = - sinLat * uNorth + cosLat * sinAlt;
  /* uLambda == uEast */

#if LALDETECTORSH_PRINTF
  printf("uNorth = %g\n",uNorth);
  printf("uEast = %g\n",uEast);
  printf("uUp = %g\n",sinAlt);
  printf("uRho = %g\n",uRho);
#endif

  u[0] = cosLon * uRho - sinLon * uEast;
  u[1] = sinLon * uRho + cosLon * uEast;
  u[2] = cosLat * uNorth + sinLat * sinAlt;

  return;
}


/* <lalVerbatim file="CreateDetectorCP"> */
LALDetector * XLALCreateDetector( LALDetector *detector,
    const LALFrDetector *frDetector, LALDetectorType type )
/* </lalVerbatim> */
{
  static const char *func = "XLALCreateDetector";
  INT2                i, j;
  REAL8               latRad, lonRad;
  REAL8               cosLat, sinLat, cosLon, sinLon;
  REAL8               locationRho, ellipsoidalDenominator;
  REAL4               xArm[3], yArm[3];
  const LALDetector  *detectorPtr, *detectorStopPtr;

  /* if detector is NULL, we are to allocate memory for it */
  if ( ! detector )
    detector = LALCalloc( 1, sizeof( *detector ) );

  if ( ! detector )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  /* if frDetector is NULL, just return a blank detector structure,
   * but set the type */
  if ( ! frDetector )
  {
    detector->type = type;
    return detector;
  }

  /* Check to see if this is a cached detector */
  detectorStopPtr = lalCachedDetectors + LAL_NUM_DETECTORS;

  for ( detectorPtr = lalCachedDetectors;
        detectorPtr < detectorStopPtr;
        ++detectorPtr )
  {
    if (  type == detectorPtr->type
          && !strncmp(detectorPtr->frDetector.name, frDetector->name,
                  LALNameLength)
          )
      {
        *detector = *detectorPtr;
        return detector;
      }
  }

  /* If it's not, construct Cartesian position vector and response tensor */

  latRad = frDetector->vertexLatitudeRadians;
  lonRad = frDetector->vertexLongitudeRadians;

#if LALDETECTORSH_PRINTF
  printf("LAT = %g radians, LON = %g radians\n", latRad, lonRad);
#endif

  cosLat = cos(latRad); sinLat = sin(latRad);
#if LALDETECTORSH_PRINTF
  printf("cos(LAT) = %g, sin(LAT) = %g\n", cosLat, sinLat);
#endif
  cosLon = cos(lonRad); sinLon = sin(lonRad);
#if LALDETECTORSH_PRINTF
  printf("cos(LON) = %g, sin(LON) = %g\n", cosLon, sinLon);
#endif

  ellipsoidalDenominator = sqrt( (LAL_AWGS84_SI * LAL_AWGS84_SI)
                            * (cosLat * cosLat)
                            + (LAL_BWGS84_SI * LAL_BWGS84_SI)
                            * (sinLat * sinLat) );

  locationRho
    = cosLat * ( (LAL_AWGS84_SI * LAL_AWGS84_SI) / ellipsoidalDenominator
                 + (REAL8) frDetector->vertexElevation );
  detector->location[0] = locationRho * cosLon;
  detector->location[1] = locationRho * sinLon;
  detector->location[2]
    = sinLat * ( (LAL_BWGS84_SI * LAL_BWGS84_SI) / ellipsoidalDenominator
                 + (REAL8) frDetector->vertexElevation );

#if LALDETECTORSH_PRINTF
  printf("%d %d\n", type, LALDETECTORTYPE_IFODIFF);
#endif

  if (type != LALDETECTORTYPE_IFOYARM)
  {
    getCartesianComponents ( xArm,
                             cos(frDetector->xArmAltitudeRadians),
                             sin(frDetector->xArmAltitudeRadians),
                             cos(frDetector->xArmAzimuthRadians),
                             sin(frDetector->xArmAzimuthRadians),
                             cosLat, sinLat, cosLon, sinLon );

#if LALDETECTORSH_PRINTF
    printf("xArm = (%g, %g, %g)\n", xArm[0], xArm[1], xArm[2]);
#endif
  }

  if (type != LALDETECTORTYPE_IFOXARM && type != LALDETECTORTYPE_CYLBAR)
  {
    getCartesianComponents ( yArm,
                             cos(frDetector->yArmAltitudeRadians),
                             sin(frDetector->yArmAltitudeRadians),
                             cos(frDetector->yArmAzimuthRadians),
                             sin(frDetector->yArmAzimuthRadians),
                             cosLat, sinLat, cosLon, sinLon );

#if LALDETECTORSH_PRINTF
    printf("yArm = (%g, %g, %g)\n", yArm[0], yArm[1], yArm[2]);
#endif
  }


  switch (type)
  {
    case LALDETECTORTYPE_IFODIFF:
      for ( i=0; i<3; ++i )
      {
        detector->response[i][i]
          = ( xArm[i] * xArm[i] - yArm[i] * yArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          detector->response[i][j] = detector->response[j][i]
            = ( xArm[i] * xArm[j] - yArm[i] * yArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_IFOXARM:
      for ( i=0; i<3; ++i )
      {
        detector->response[i][i]
          = ( xArm[i] * xArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          detector->response[i][j] = detector->response[j][i]
            = ( xArm[i] * xArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_IFOYARM:
      for ( i=0; i<3; ++i )
      {
        detector->response[i][i]
          = ( yArm[i] * yArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          detector->response[i][j] = detector->response[j][i]
            = ( yArm[i] * yArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_IFOCOMM:
      for ( i=0; i<3; ++i )
      {
        detector->response[i][i]
          = ( xArm[i] * xArm[i] + yArm[i] * yArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          detector->response[i][j] = detector->response[j][i]
            = ( xArm[i] * xArm[j] + yArm[i] * yArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_CYLBAR:
      for ( i=0; i<3; ++i )
      {
        detector->response[i][i]
          = xArm[i] * xArm[i];
        for ( j=i+1; j<3; ++j )
        {
          detector->response[i][j] = detector->response[j][i]
            = xArm[i] * xArm[j];
        }
      }
      break;
    default:
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
  } /* switch (type) */

  detector->frDetector = *frDetector;
  detector->type = type;
  return detector;
}


void LALCreateDetector( LALStatus             *status,
                        LALDetector           *output,
                        const LALFrDetector   *input,
                        const LALDetectorType  type )
{
  INITSTATUS( status, "LALCreateDetector", CREATEDETECTORC );

  ASSERT( input != NULL, status, LALDETECTORSH_ENULLP,
          LALDETECTORSH_MSGENULLP );

  ASSERT( output != NULL, status, LALDETECTORSH_ENULLP,
          LALDETECTORSH_MSGENULLP );

  output = XLALCreateDetector( output, input, type );
  if ( ! output )
    switch ( XLALClearErrno() )
    {
      case XLAL_EINVAL:
        ABORT( status, LALDETECTORSH_ETYPE, LALDETECTORSH_MSGETYPE );
        break;
      default:
        ABORTXLAL( status );
    }

  RETURN(status);
}
