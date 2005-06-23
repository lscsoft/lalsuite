/************************************ <lalVerbatim file="CreateDetectorCV">
Author: J. T. Whelan <jtwhelan@loyno.edu>
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{CreateDetector.c}}
\label{tools:ss:CreateDetector.c}

Creates a \texttt{LALDetector} structure from a \texttt{LALFrDetector}
structure and the type of detector.

\subsubsection*{Prototypes}
\input{CreateDetectorCP}
\idx{LALCreateDetector()}

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
\texttt{LALCreateDetector()} first checks the
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
LALCreateDetector()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item The conventions in the \texttt{LALFrDetector} structure are an
  extension of those used in version 4 of the frame specification
  \cite{tools:LIGOVIRGO:2000} and do not agree with the conventions
  used by LDAS, which are an anticipation of the yet-to-be-released
  version 5 frame spec.  The \texttt{LALFrDetector} structure will
  change once a final version 5 frame spec is available.
\item If the location and response tensor information for a
\texttt{LALDetector} are filled in by hand (e.g., for testing
purposes), the \texttt{type} field should be set to
\texttt{LALDETECTORTYPE\_ABSENT}.
\item The range of \texttt{LALDetectorType}s could be expanded to
include the  monopole and five quadrupole modes for a spherical
resonant detector
\cite{tools:Maggiore:2000b,tools:Zhou:1995,tools:Bianchi:1998,tools:Maggiore:2000a}.
\item At the moment, this code still writes some diagnostics to
standard output.  These are supposed to be removed once it's been tested.
\end{itemize}

\vfill{\footnotesize\input{CreateDetectorCV}}

******************************************************* </lalLaTeX> */

/**************************************** <lalLaTeX file="CreateDetectorCB">
\bibitem{tools:Anderson:2000}
  W.~G.~Anderson, J.~T.~Whelan, P.~R.~Brady, J.~D.~E.~Creighton,
  D.~Chin, and K.~Riles, ``Beam Pattern Response Functions and Times
  of Arrival for Earthbound Interferometers'', Maple worksheet,
  \href{http://dirac.phys.utb.edu/~warren/unprot/beam\_patterns\_260401.tar.gz}
  {\texttt{http://dirac.phys.utb.edu/\~ warren/unprot/beam\_patterns\_260401.tar.gz}}
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
      index & \texttt{LALDetectorIndexLHODIFF}
            & \texttt{LALDetectorIndexLLODIFF}
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
     name & LIGO Hanford Observatory & LIGO Livingston Observatory
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
\hline
      index & \texttt{LALDetectorIndexVIRGODIFF}
            & \texttt{LALDetectorIndexGEO600DIFF}
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
     name & VIRGO Interferometer & GEO-600 Interferometer
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
\hline
      index & \texttt{LALDetectorIndexTAMA300DIFF}
            & \texttt{LALDetectorIndexCIT40DIFF}
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
     name & TAMA-300 Interferometer & Caltech-40 Interferometer
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

    \end{tabular}
    \caption{Predefined gravitational wave detectors, contained in
      the \texttt{lalCachedDetectors[]} array.  The LIGO site data
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

#if 0
const LALDetector lalCachedDetectors[LALNumCachedDetectors]
= { { { -2.1614149e+06L, -3.8346952e+06L,   4.6003502e+06L },
      { { -0.3926141, -0.0776130, -0.2473886 },
        { -0.0776130,  0.3195244,  0.2279981 },
        { -0.2473886,  0.2279981,  0.0730903 }
      },
      LALDETECTORTYPE_IFODIFF,
      { "LIGO Hanford Observatory",
	- (119.0L + 24.0L/60.0L + 27.5657L/3600.0L ) * LAL_PI_180,
	(46.0L + 27.0L/60.0L + 18.528L/3600.0L) * LAL_PI_180,
	142.554,
	-6.195e-4,     (360.0L + 90.0L - 125.9994L) * LAL_PI_180,
	1.25e-5,       (360.0L + 90.0L - 215.9994L) * LAL_PI_180
      }
    },
    { { -74276.044L,     -5.496283721e+06L, 3.224257018e+06L },
      { {  0.4112809,  0.1402097,  0.2472943 },
        {  0.1402097, -0.1090056, -0.1816157 },
        {  0.2472943, -0.1816157, -0.3022755 },
      },
      LALDETECTORTYPE_IFODIFF,
      { "LIGO Livingston Observatory",
	- (90.0L + 46.0L/60.0L + 27.2654L/3600.0L) * LAL_PI_180,
	(30.0L + 33.0L/60.0L + 46.4196L/3600.0L) * LAL_PI_180,
	-6.574,
	-3.121e-4,   (360.0L + 90.0L - 197.7165L) * LAL_PI_180, 
	-6.107e-4,   (360.0L + 90.0L - 287.7165L) * LAL_PI_180, 
      }
    },
    { { 4.54637409863e+06L, 842989.697467L, 4.37857696275e+06L},
      { {  0.2438740,  -0.0990838, -0.2325762},
        { -0.0990838, - 0.4478258,  0.1878331},
        { -0.2325762,   0.1878331,  0.2039518},
      },
      LALDETECTORTYPE_IFODIFF,
      { "VIRGO Interferometer",
	(10.0L + 30.0L/60.0L + 16.1878L/3600.0L) * LAL_PI_180,
	(43.0L + 37.0L/60.0L + 53.0921L/3600.0L) * LAL_PI_180,
	51.884,
        0.0,           (90.0L - 70.5674L) * LAL_PI_180,
        0.0,           (360.0L + 90.0L - 160.5674L) * LAL_PI_180
      }
    },
    { { 3.85630994953e+06L,   666598.956352L,   5.01964141692e+06L},
      { { -0.0968250, -0.3657823,   0.1221373},
        { -0.3657823,  0.222968,    0.2497174},
        {  0.1221373,  0.2497174,  -0.1261431},
      },
      LALDETECTORTYPE_IFODIFF,
      { "GEO-600 Interferometer",
	(9.0L + 48.0L/60.0L + 25.894L/3600.0L) * LAL_PI_180, 
	(52.0L + 14.0L/60.0L + 42.528L/3600.0L) * LAL_PI_180,
	114.425,
	0.0,         (90.0L - 21.6117L) * LAL_PI_180,
	0.0,         (360.0L + 90.0L - 115.9431L) * LAL_PI_180
      }
    },
    { { -3.94640898771e+06L,  3.36625903242e+06L, 3.69915069189e+06L},
      { {  0.1121397,  0.3308421, -0.1802193},
        {  0.3308421,  0.2177940,  0.1537258},
        { -0.1802193,  0.1537258, -0.3299337},
      },
      LALDETECTORTYPE_IFODIFF,
      { "TAMA-300 Interferometer",
	(139.0L + 32.0L/60.0L + 9.8L/3600.0L) * LAL_PI_180,
	(35.0L + 40.0L/60.0L + 35.6L/3600.0L) * LAL_PI_180,
	90,
	0.0, (360.0L + 90.0L - 180.0L) * LAL_PI_180,
	0.0, (360.0L + 90.0L - 270.0L) * LAL_PI_180
      }
    },
    { { -2.49064958399e+06L,  -4.65869968229e+06L,  3.56206411337e+06L},
      { { -0.3537959,  0.2734713, 0.1095458},
        {  0.2734713,  0.0115214, 0.2049027},
        {  0.1095458,  0.2049027, 0.3422745},
      },
      LALDETECTORTYPE_IFODIFF,
      { "Caltech-40 Interferometer",
	-118.13 * LAL_PI_180,
	34.17 * LAL_PI_180, 
	0,  
	0.0,            (360.0L + 90.0L - 270.0L) * LAL_PI_180,
	0.0,            (90.0L - 0.0L) * LAL_PI_180
      }
    }
};
#endif

const LALDetector lalCachedDetectors[LALNumCachedDetectors] =
{
  /* Cached Detector 0: LIGO Hanford Observatory */
  {
    /* detector position vector (m) */
    { -2.16141490000e+06, -3.83469520000e+06, +4.60035020000e+06 },
    /* detector response matrix */
    {
      { -0.39261409640, -0.07761300355, -0.24738860130 },
      { -0.07761300355, +0.31952440739, +0.22799809277 },
      { -0.24738860130, +0.22799809277, +0.07309029996 },
    },
    LALDETECTORTYPE_IFODIFF, 	/* detector type */
    /* LALFrDetector structure */
    {
      "LIGO Hanford Observatory", 	/* detector name */
      -2.08405676917, 	/* vertex longitude (rad) */
      +0.81079526383, 	/* vertex latitude (rad) */
      +142.554      , 	/* vertex elevation (m) */
      -0.00061950000, 	/* x arm altitude (rad) */
      +5.65487718582, 	/* x arm azimuth (rad) */
      +0.00001250000, 	/* y arm altitude (rad) */
      +4.08408069611  	/* y arm azimuth (rad) */
    },
  },
  /* Cached Detector 1: LIGO Livingston Observatory */
  {
    /* detector position vector (m) */
    { -7.42760440000e+04, -5.49628372100e+06, +3.22425701800e+06 },
    /* detector response matrix */
    {
      { +0.41128090024, +0.14020970464, +0.24729430676 },
      { +0.14020970464, -0.10900560021, -0.18161569536 },
      { +0.24729430676, -0.18161569536, -0.30227550864 },
    },
    LALDETECTORTYPE_IFODIFF, 	/* detector type */
    /* LALFrDetector structure */
    {
      "LIGO Livingston Observatory", 	/* detector name */
      -1.58430937078, 	/* vertex longitude (rad) */
      +0.53342313506, 	/* vertex latitude (rad) */
      -6.574        , 	/* vertex elevation (m) */
      -0.00031209999, 	/* x arm altitude (rad) */
      +4.40317773819, 	/* x arm azimuth (rad) */
      -0.00061069999, 	/* y arm altitude (rad) */
      +2.83238148689  	/* y arm azimuth (rad) */
    },
  },
  /* Cached Detector 2: VIRGO Interferometer */
  {
    /* detector position vector (m) */
    { +4.54637409863e+06, +8.42989697467e+05, +4.37857696275e+06 },
    /* detector response matrix */
    {
      { +0.24387399852, -0.09908380359, -0.23257620633 },
      { -0.09908380359, -0.44782578945, +0.18783310056 },
      { -0.23257620633, +0.18783310056, +0.20395180583 },
    },
    LALDETECTORTYPE_IFODIFF, 	/* detector type */
    /* LALFrDetector structure */
    {
      "VIRGO Interferometer", 	/* detector name */
      +0.18333805213, 	/* vertex longitude (rad) */
      +0.76151183984, 	/* vertex latitude (rad) */
      +51.884       , 	/* vertex elevation (m) */
      +0.00000000000, 	/* x arm altitude (rad) */
      +0.33916285634, 	/* x arm azimuth (rad) */
      +0.00000000000, 	/* y arm altitude (rad) */
      +5.05155181885  	/* y arm azimuth (rad) */
    },
  },
  /* Cached Detector 3: GEO-600 Interferometer */
  {
    /* detector position vector (m) */
    { +3.85630994953e+06, +6.66598956352e+05, +5.01964141692e+06 },
    /* detector response matrix */
    {
      { -0.09682500362, -0.36578229070, +0.12213730067 },
      { -0.36578229070, +0.22296799719, +0.24971739948 },
      { +0.12213730067, +0.24971739948, -0.12614309788 },
    },
    LALDETECTORTYPE_IFODIFF, 	/* detector type */
    /* LALFrDetector structure */
    {
      "GEO-600 Interferometer", 	/* detector name */
      +0.17116780435, 	/* vertex longitude (rad) */
      +0.91184982752, 	/* vertex latitude (rad) */
      +114.425      , 	/* vertex elevation (m) */
      +0.00000000000, 	/* x arm altitude (rad) */
      +1.19360101223, 	/* x arm azimuth (rad) */
      +0.00000000000, 	/* y arm altitude (rad) */
      +5.83039283752  	/* y arm azimuth (rad) */
    },
  },
  /* Cached Detector 4: TAMA-300 Interferometer */
  {
    /* detector position vector (m) */
    { -3.94640898771e+06, +3.36625903242e+06, +3.69915069189e+06 },
    /* detector response matrix */
    {
      { +0.11213970184, +0.33084210753, -0.18021929264 },
      { +0.33084210753, +0.21779400110, +0.15372580290 },
      { -0.18021929264, +0.15372580290, -0.32993370295 },
    },
    LALDETECTORTYPE_IFODIFF, 	/* detector type */
    /* LALFrDetector structure */
    {
      "TAMA-300 Interferometer", 	/* detector name */
      +2.43536359469, 	/* vertex longitude (rad) */
      +0.62267336022, 	/* vertex latitude (rad) */
      +90.000       , 	/* vertex elevation (m) */
      +0.00000000000, 	/* x arm altitude (rad) */
      +4.71238899231, 	/* x arm azimuth (rad) */
      +0.00000000000, 	/* y arm altitude (rad) */
      +3.14159274101  	/* y arm azimuth (rad) */
    },
  },
  /* Cached Detector 5: Caltech-40 Interferometer */
  {
    /* detector position vector (m) */
    { -2.49064958399e+06, -4.65869968229e+06, +3.56206411337e+06 },
    /* detector response matrix */
    {
      { -0.35379588604, +0.27347129583, +0.10954579711 },
      { +0.27347129583, +0.01152139995, +0.20490269363 },
      { +0.10954579711, +0.20490269363, +0.34227448702 },
    },
    LALDETECTORTYPE_IFODIFF, 	/* detector type */
    /* LALFrDetector structure */
    {
      "Caltech-40 Interferometer", 	/* detector name */
      -2.06175744538, 	/* vertex longitude (rad) */
      +0.59637900541, 	/* vertex latitude (rad) */
      +0.000        , 	/* vertex elevation (m) */
      +0.00000000000, 	/* x arm altitude (rad) */
      +3.14159274101, 	/* x arm azimuth (rad) */
      +0.00000000000, 	/* y arm altitude (rad) */
      +1.57079637051  	/* y arm azimuth (rad) */
    },
  }
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

#if DETECTORSITEH_PRINTF
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
void LALCreateDetector( LALStatus             *status,
                        LALDetector           *output,
                        const LALFrDetector   *input,
                        const LALDetectorType  type )
/* </lalVerbatim> */
{
  INT2                i, j;
  REAL8               latRad, lonRad;
  REAL8               cosLat, sinLat, cosLon, sinLon;
  REAL8               locationRho, ellipsoidalDenominator;
  REAL4               xArm[3], yArm[3];
  const LALDetector  *detectorPtr, *detectorStopPtr;

  INITSTATUS( status, "LALCreateDetector", CREATEDETECTORC );

  ASSERT( input != NULL, status, DETECTORSITEH_ENULLP,
          DETECTORSITEH_MSGENULLP );

  ASSERT( output != NULL, status, DETECTORSITEH_ENULLP,
          DETECTORSITEH_MSGENULLP );

  /* Check to see if this is a cached detector */

  detectorStopPtr = lalCachedDetectors + LALNumCachedDetectors;

  for ( detectorPtr = lalCachedDetectors;
        detectorPtr < detectorStopPtr;
        ++detectorPtr )
  {
    if (  type == detectorPtr->type
          && !strncmp(detectorPtr->frDetector.name, input->name,
                  LALNameLength)
          )
      {
        *output = *detectorPtr;
        RETURN(status);
      }
  }

  /* If it's not, construct Cartesian position vector and response tensor */

  latRad = input->vertexLatitudeRadians;
  lonRad = input->vertexLongitudeRadians;

#if DETECTORSITEH_PRINTF
  printf("LAT = %g radians, LON = %g radians\n", latRad, lonRad);
#endif

  cosLat = cos(latRad); sinLat = sin(latRad);
#if DETECTORSITEH_PRINTF
  printf("cos(LAT) = %g, sin(LAT) = %g\n", cosLat, sinLat);
#endif
  cosLon = cos(lonRad); sinLon = sin(lonRad);
#if DETECTORSITEH_PRINTF
  printf("cos(LON) = %g, sin(LON) = %g\n", cosLon, sinLon);
#endif

  ellipsoidalDenominator = sqrt( (LAL_AWGS84_SI * LAL_AWGS84_SI)
                            * (cosLat * cosLat)
                            + (LAL_BWGS84_SI * LAL_BWGS84_SI)
                            * (sinLat * sinLat) );

  locationRho
    = cosLat * ( (LAL_AWGS84_SI * LAL_AWGS84_SI) / ellipsoidalDenominator
                 + (REAL8) input->vertexElevation );
  output->location[0] = locationRho * cosLon;
  output->location[1] = locationRho * sinLon;
  output->location[2]
    = sinLat * ( (LAL_BWGS84_SI * LAL_BWGS84_SI) / ellipsoidalDenominator
                 + (REAL8) input->vertexElevation );

#if DETECTORSITEH_PRINTF
  printf("%d %d\n", type, LALDETECTORTYPE_IFODIFF);
#endif

  if (type != LALDETECTORTYPE_IFOYARM)
  {
    getCartesianComponents ( xArm,
                             cos(input->xArmAltitudeRadians),
                             sin(input->xArmAltitudeRadians),
                             cos(input->xArmAzimuthRadians),
                             sin(input->xArmAzimuthRadians),
                             cosLat, sinLat, cosLon, sinLon );

#if DETECTORSITEH_PRINTF
    printf("xArm = (%g, %g, %g)\n", xArm[0], xArm[1], xArm[2]);
#endif
  }

  if (type != LALDETECTORTYPE_IFOXARM && type != LALDETECTORTYPE_CYLBAR)
  {
    getCartesianComponents ( yArm,
                             cos(input->yArmAltitudeRadians),
                             sin(input->yArmAltitudeRadians),
                             cos(input->yArmAzimuthRadians),
                             sin(input->yArmAzimuthRadians),
                             cosLat, sinLat, cosLon, sinLon );

#if DETECTORSITEH_PRINTF
    printf("yArm = (%g, %g, %g)\n", yArm[0], yArm[1], yArm[2]);
#endif
  }


  switch (type)
  {
    case LALDETECTORTYPE_IFODIFF:
      for ( i=0; i<3; ++i )
      {
        output->response[i][i]
          = ( xArm[i] * xArm[i] - yArm[i] * yArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          output->response[i][j] = output->response[j][i]
            = ( xArm[i] * xArm[j] - yArm[i] * yArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_IFOXARM:
      for ( i=0; i<3; ++i )
      {
        output->response[i][i]
          = ( xArm[i] * xArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          output->response[i][j] = output->response[j][i]
            = ( xArm[i] * xArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_IFOYARM:
      for ( i=0; i<3; ++i )
      {
        output->response[i][i]
          = ( yArm[i] * yArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          output->response[i][j] = output->response[j][i]
            = ( yArm[i] * yArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_IFOCOMM:
      for ( i=0; i<3; ++i )
      {
        output->response[i][i]
          = ( xArm[i] * xArm[i] + yArm[i] * yArm[i] ) / 2;
        for ( j=i+1; j<3; ++j )
        {
          output->response[i][j] = output->response[j][i]
            = ( xArm[i] * xArm[j] + yArm[i] * yArm[j] ) / 2;
        }
      }
      break;
    case LALDETECTORTYPE_CYLBAR:
      for ( i=0; i<3; ++i )
      {
        output->response[i][i]
          = xArm[i] * xArm[i];
        for ( j=i+1; j<3; ++j )
        {
          output->response[i][j] = output->response[j][i]
            = xArm[i] * xArm[j];
        }
      }
      break;
    default:
      ABORT( status, DETECTORSITEH_ETYPE, DETECTORSITEH_MSGETYPE );
  } /* switch (type) */

  output->frDetector = *input;
  output->type = type;

  RETURN(status);
}





