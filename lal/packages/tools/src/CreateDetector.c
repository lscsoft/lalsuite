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

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <math.h>
#include <string.h>
#include <lal/DetectorSite.h>

/**
 * \defgroup CreateDetector_c Module CreateDetector.c
 * \ingroup LALDetectors_h
 *
 * \author J. T. Whelan <john.whelan@ligo.org>
 *
 * \brief Creates a \c LALDetector structure from a \c LALFrDetector structure and the type of detector.
 *
 * This routine takes the site geometry described in the
 * \c LALFrDetector structure, along with a
 * \c LALDetectorType parameter, and constructs the Cartesian
 * detector location and response tensor needed to fill the
 * \c LALDetector output.
 *
 * The detector type is needed because different types of detectors have
 * different response tensors.  In each case the response tensor is
 * determined by the unit vectors \f$\hat{u}_X\f$ and \f$\hat{u}_Y\f$
 * which are constant in an Earth-fixed rotating reference frame and
 * point in the "X arm" and "Y arm" directions, respectively; the
 * headings of these directions in a local frame at the detector are
 * stored in the \c LALFrDetector structure.
 *
 * The detector types recognized are (all names are prefaced by \c LALDETECTORTYPE_):
 * <ul>
 * <li>[\c IFODIFF] An interferometer in differential mode.  The
 * response tensor is given by \f$d^{ab}=\frac{1}{2} (u_X^au_X^b-u_Y^au_Y^b)\f$.
 * Note that this is the preferred form even in the two arms of the
 * detector are not perpendicular (e.g., at the GEO600 site).</li>
 * <li>[\c IFOXARM] An interferometer in one-armed mode with the
 * X arm active.  The response tensor is given by
 * \f$d^{ab}=\frac{1}{2}u_X^au_X^b\f$.</li>
 * <li>[\c IFOYARM] An interferometer in one-armed mode with the
 * Y arm active.  The response tensor is given by
 * \f$d^{ab}=\frac{1}{2}u_Y^au_Y^b\f$.</li>
 * <li>[\c IFOCOMM] An interferometer in common mode.  The
 * response tensor is given by \f$d^{ab}=\frac{1}{2} (u_X^au_X^b+u_Y^au_Y^b)\f$.</li>
 * <li>[\c CYLBAR] A cylindrical bar detector.  In this case the
 * "X arm" is actually the symmetry axis of the bar, and the "Y arm"
 * is ignored.  The response tensor is
 * \f$d^{ab}=u_X^au_X^b\f$.</li>
 * </ul>
 *
 * In each of these cases, the basic transformation needed is to express a
 * unit vector \f$\hat{u}\f$ in terms of its
 * components in the Earth-fixed basis
 * \f$\{\hat{e}_1,\hat{e}_2,\hat{e}_3\}\f$.  The altitude angle \f${\mathcal{A}}\f$ and
 * azimuth angle \f$\zeta\f$ allow us to express the unit vector  \f$\hat{u}\f$
 * corresponding to a direction in terms of an orthonormal basis consisting
 * of a vector \f$\hat{e}_{\scriptstyle\textrm{E}}\f$ pointing due East within the
 * local horizontal, a vector \f$\hat{e}_{\scriptstyle\textrm{N}}\f$ pointing due
 * North within the local horizontal, and an upward-pointing vector
 * \f$\hat{e}_{\scriptstyle\textrm{U}}\f$ normal to the local horizontal
 * plane. [These form a right-handed basis, providing an answer to
 * the age-old question "What's Up?": "East cross North."]
 * The relationship is
 * \f{equation}{
 * \hat{u} =   ( \hat{e}_{\scriptstyle\textrm{E}}\sin\zeta
 * + \hat{e}_{\scriptstyle\textrm{N}}\cos\zeta )
 * \cos{\mathcal{A}}
 * + \hat{e}_{\scriptstyle\textrm{U}} \sin{\mathcal{A}}
 * \ .
 * \f}
 * Since the local horizontal is defined as the tangent plane to the
 * reference ellipsoid at the point with the detector's latitude \f$\beta\f$
 * and longitude \f$\lambda\f$, the local basis is related to the orthonormal
 * basis
 * \f$\{\hat{e}_\rho,\hat{e}_\lambda,\hat{e}_z\}\f$ of a cylindrical
 * coördinate system (related to the Earth-fixed Cartesian
 * coördinates by \f$x^1=\rho\cos\lambda\f$, \f$x^2=\rho\sin\lambda\f$, \f$x^3=z\f$,
 * so that \f$\hat{e}_\rho\f$ points away from the Earth's axis,
 * \f$\hat{e}_\lambda\f$ points in the direction of increasing longitude, and
 * n\f$\hat{e}_z\f$ points in the direction of increasing \f$x^3\f$)
 * by
 * \f{align}{
 * \hat{e}_{\scriptstyle\textrm{E}} &= \hat{e}_\lambda \\
 * \hat{e}_{\scriptstyle\textrm{N}} &= - \hat{e}_\rho \sin\beta
 * + \hat{e}_z \cos\beta \\
 * \hat{e}_{\scriptstyle\textrm{U}} &=   \hat{e}_\rho \cos\beta
 * + \hat{e}_z \sin\beta
 * \f}
 * It is then straightforward to relate the cylindrical basis vectors to
 * those in the Earth-fixed Cartesian system by
 * \f{align}{
 * \hat{e}_\rho    &=  \hat{e}_1\cos\lambda  +  \hat{e}_2\sin\lambda  \\
 * \hat{e}_\lambda &= -\hat{e}_1\sin\lambda  +  \hat{e}_2\cos\lambda  \\
 * \hat{e}_z       &= \hat{e}_3
 * \f}
 *
 * To express \f$\hat{u}\f$ in the Cartesian basis, we need
 * \f$\hat{u}\cdot\hat{e}_1\f$, \f$\hat{u}\cdot\hat{e}_2\f$, and
 * \f$\hat{u}\cdot\hat{e}_3\f$.  We first observe that
 * \anchor tools_e_eE \f{align}{
 * \hat{u}\cdot\hat{e}_{\scriptstyle\textrm{E}} &= \cos{\mathcal{A}}\,\sin\zeta \tag{tools_e_eE}\\
 * \hat{u}\cdot\hat{e}_{\scriptstyle\textrm{N}} &= \cos{\mathcal{A}}\,\cos\zeta \\
 * \hat{u}\cdot\hat{e}_{\scriptstyle\textrm{U}} &= \sin{\mathcal{A}}
 * \f}
 * then that
 * \f{align}{
 * \hat{u}\cdot\hat{e}_\rho &= (\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{N}})
 * (\hat{e}_{\scriptstyle\textrm{N}}\cdot\hat{e}_\rho)
 * + (\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{U}})
 * (\hat{e}_{\scriptstyle\textrm{U}}\cdot\hat{e}_\rho)
 * = -(\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{N}}) \sin\beta
 * +(\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{U}}) \cos\beta \\
 * \hat{u}\cdot\hat{e}_\lambda &=& \hat{u}\cdot\hat{e}_{\scriptstyle\textrm{E}}\\
 * \hat{u}\cdot\hat{e}_z &= (\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{N}})
 * (\hat{e}_{\scriptstyle\textrm{N}}\cdot\hat{e}_z)
 * + (\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{U}})
 * (\hat{e}_{\scriptstyle\textrm{U}}\cdot\hat{e}_z)
 * = (\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{N}}) \cos\beta
 * +(\hat{u}\cdot\hat{e}_{\scriptstyle\textrm{U}}) \sin\beta
 * \f}
 * and finally that
 * \anchor tools_e_e3ez \f{align}{
 * \hat{u}\cdot\hat{e}_1 &= (\hat{u}\cdot\hat{e}_\rho)
 * (\hat{e}_\rho\cdot\hat{e}_1)
 * + (\hat{u}\cdot\hat{e}_\lambda)
 * (\hat{e}_\lambda\cdot\hat{e}_1)
 * = (\hat{u}\cdot\hat{e}_\rho) \cos\lambda
 * -(\hat{u}\cdot\hat{e}_\lambda) \sin\lambda \\
 * \hat{u}\cdot\hat{e}_2 &= (\hat{u}\cdot\hat{e}_\rho)
 * (\hat{e}_\rho\cdot\hat{e}_2)
 * + (\hat{u}\cdot\hat{e}_\lambda)
 * (\hat{e}_\lambda\cdot\hat{e}_2)
 * = (\hat{u}\cdot\hat{e}_\rho) \sin\lambda
 * +(\hat{u}\cdot\hat{e}_\lambda) \cos\lambda \\
 * \hat{u}\cdot\hat{e}_3 &= \hat{u}\cdot\hat{e}_z
 * \tag{tools_e_e3ez}
 * \f}
 *
 * ### Cached Detectors ###
 *
 * To avoid repeatedly calculating the Cartesian coördinates and
 * response tensor of known detectors, the constant array
 * <tt>lalCachedDetectors[]</tt> contains the site geometry and
 * response tensors of the most commonly used detectors.  These are
 * defined in this file and listed in Table.\tableref{tools_tab_cached}.
 *
 * ### Algorithm ###
 *
 * <tt>XLALCreateDetector()</tt> first checks the
 * <tt>lalCachedDetectors[]</tt> array to see if the specified type and
 * the name in the input \c LALFrDetector match any of the
 * predefined constant detectors.  If so, it returns a copy of the
 * constant detector (not just a pointer to the constant).
 *
 * If not, it calculates the Cartesian coördinates \f$\{x^1,x^2,x^3\}\f$
 * of the detector location defined by\eqref{tools_e_cart1}-\eqref{tools_e_cart3};
 * in
 * particular, it calculates the denominator
 * \f$\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}\f$ and the distance from the axis
 * \f{equation}{
 * \rho = \left(\frac{a^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}} + h \right)
 * \cos\beta
 * \f}
 * as intermediate steps.
 *
 * It then calculates the Cartesian components of the unit vectors
 * \f$\hat{u}_X\f$ and \f$\hat{u}_Y\f$ in the arm directions from the altitude
 * and azimuth angles by use of a \c static function which
 * implements\eqref{tools_e_eE}-\eqref{tools_e_e3ez}.  (Depending on the detector
 * type specified, only the unit vector(s) which are actually needed are
 * calculated.)  Using this components it constructs \f$d^{ab}\f$ according
 * to the formula appropriate to the detector type.
 *
 * The calculation of \f$x^a\f$ is done to double precision, that of \f$d^{ab}\f$
 * to single precision.
 *
 * ### Notes ###
 *
 * <ul>
 * <li> The conventions in the \c LALFrDetector structure are based
 * on version 6 of the frame specification [\ref LIGOVIRGO_2000].</li>
 * <li> If the location and response tensor information for a
 * \c ::LALDetector are filled in by hand (e.g., for testing
 * purposes), the \c type field should be set to
 * \c #LALDETECTORTYPE_ABSENT.</li>
 * <li> The range of \c ::LALDetectorType could be expanded to
 * include the  monopole and five quadrupole modes for a spherical
 * resonant detector
 * [\ref Maggiore_2000b, \ref Zhou_1995, \ref Bianchi_1998, \ref Maggiore_2000a].</li>
 * </ul>
 *
 * <b>Table tools_tab_cached:</b> Selected redefined gravitational wave detectors, contained in the ::lalCachedDetectors array.
 * Not shown in the table are the LHO 2km detector (H2) and the bar detectors ALLEGRO, AURIGA, EXPLORER, NIOBE and NAUTILUS.
 * The LIGO site data come directly from [\ref Althouse_1999], including the Cartesian position vectors \f$x^a\f$ and the response tensor
 * \f$d^{ab}\f$, which was dermined from the quoted components of the detector frame basis vectors \f$\hat{x}_G\equiv\hat{u}_X\f$ and
 * \f$\hat{y}_G\equiv\hat{u}_Y\f$.  The data on the other detectors comes from [\ref Anderson_2000].
 * \anchor tools_tab_cached
 * <table>
 * <tr><th>index</th><th>\c LAL_LHO_4K_DETECTOR</th><th>\c LAL_LLO_4K_DETECTOR</th></tr>
 * <tr><td>      prefix</td><td>\c H1                 </td><td>\c L1 </td></tr>
 * <tr><td>   \f$x^a\f$</td><td>
 * \f$
 * \left(
 * \begin{array}{c}
 * -2.1614149 \times 10^6 \\
 * -3.8346952 \times 10^6 \\
 * 4.6003502 \times 10^6 \\
 * \end{array}
 * \right)
 * \f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{c}
 * -74276.044\\
 * -5.496283721\times 10^6\\
 * 3.224257018\times 10^6\\
 * \end{array}
 * \right)
 * \f$
 * </td></tr>
 * <tr><td>  \f$d^{ab}\f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{ccc}
 * v          -0.3926141  & -0.0776130  & -0.2473886 \\
 * -0.0776130  &  0.3195244  &  0.2279981 \\
 * -0.2473886  &  0.2279981  &  0.0730903 \\
 * \end{array}
 * \right)
 * \f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{ccc}
 * 0.4112809  &  0.1402097  &  0.2472943 \\
 * 0.1402097  & -0.1090056  & -0.1816157 \\
 * 0.2472943  & -0.1816157  & -0.3022755 \\
 * \end{array}
 * \right)
 * \f$
 * </td></tr>
 * <tr><td>     type</td><td>\c LALDETECTORTYPE_IFODIFF</td><td>\c LALDETECTORTYPE_IFODIFF</td></tr>
 * <tr><td>     name</td><td>LHO_4k</td><td>LLO_4k</td></tr>
 * <tr><td>     \f$(\lambda,\beta,h)\f$
 * </td><td>\f$(-(119^\circ24'27"\!\!.5657),46^\circ27'18"\!\!.528, 142.544\,\textrm{m})\f$
 * </td><td>\f$(-(90^\circ46'27"\!\!.2654),30^\circ33'46\!\!.4196,   -6.574\,\textrm{m})\f$
 * </td></tr>
 * <tr><td>      \f$({\mathcal{A}}_X,\zeta_X)\f$
 * </td><td>\f$(         -6.195\times 10^{-4},      324^\circ\!\!.0006)\f$
 * </td><td>\f$(          -3.121\times 10^{-4},     252^\circ\!\!.2835)\f$
 * </td></tr>
 * <tr><td>      \f$({\mathcal{A}}_Y,\zeta_Y)\f$
 * </td><td>\f$(           1.25\times 10^{-5},      234^\circ\!\!.0006)\f$
 * </td><td>\f$(          -6.107\times 10^{-4},     162^\circ\!\!.2835)\f$
 * </td></tr>
 * <tr><td>     \f$(L_X/2,L_Y/2)\f$ </td><td>\f$(2000\,\textrm{m},2000\,\textrm{m})\f$ </td><td>\f$(2000\,\textrm{m},2000\,\textrm{m})\f$ </td></tr>
 * </table>
 *
 * <br>
 * <table class="doxtable">
 * <tr><th> index</th><th>\c LAL_VIRGO_DETECTOR </th><th>\c LAL_GEO_600_DETECTOR</th></tr>
 * <tr><td>      \f$x^a\f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{c}
 * 4.54637409863 \times 10^6 \\
 * 842989.697467  \\
 * 4.37857696275\times 10^6 \\
 * \end{array}
 * \right)
 * \f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{c}
 * 3.85630994953 \times 10^6 \\
 * 666598.956352 \\
 * 5.01964141692 \times 10^6 \\
 * \end{array}
 * \right)
 * \f$
 * </td></tr>
 * <tr><td>      \f$d^{ab}\f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{ccc}
 * 0.2438740  & -0.0990838  & -0.2325762 \\
 * -0.0990838 & -0.4478258  & 0.1878331  \\
 * -0.2325762 & 0.1878331   & 0.2039518  \\
 * \end{array}
 * \right)
 * \f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{ccc}
 * -0.0968250  & -0.3657823  &  0.1221373  \\
 * -0.3657823  &  0.2229681  &  0.2497174  \\
 * 0.1221373  &  0.2497174  & -0.1261431	\\
 * \end{array}
 * \right)
 * \f$
 * </td></tr>
 * <tr><td>     type</td><td>\c LALDETECTORTYPE_IFODIFF</td><td>\c LALDETECTORTYPE_IFODIFF </td></tr>
 * <tr><td>     name</td><td>VIRGO</td><td>GEO_600 </td></tr>
 * <tr><td>     \f$(\lambda,\beta,h)\f$
 * </td><td>\f$(10^\circ30'16"\!\!.1878,43^\circ37'\!\!53".0921, 51.884\,\textrm{m})\f$
 * </td><td>\f$(9^\circ48'25"\!\!.894,52^\circ14'42"\!\!.528, 114.425\,\textrm{m})\f$
 * </td></tr>
 * <tr><td>      \f$({\mathcal{A}}_X,\zeta_X)\f$
 * </td><td>\f$( 0,          19^\circ\!\!.4326)\f$
 * </td><td>\f$( 0,          68^\circ\!\!.3883)\f$
 * </td></tr>
 * <tr><td>      \f$({\mathcal{A}}_Y,\zeta_Y)\f$
 * </td><td>\f$( 0,           289^\circ\!\!.4326)\f$
 * </td><td>\f$( 0,           334^\circ\!\!.0569)\f$
 * </td></tr>
 * <tr><td>     \f$(L_X/2,L_Y/2)\f$ </td><td>\f$(1500\,\textrm{m},1500\,\textrm{m})\f$ </td><td>\f$(300\,\textrm{m},300\,\textrm{m})\f$ </td></tr>
 * </table>
 *
 * <br>
 * <table class="doxtable">
 * <tr><th>index</th><th>\c LAL_TAMA_300_DETECTOR </th><th>\c LAL_CIT_40_DETECTOR </th></tr>
 * <tr><td>      \f$x^a\f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{c}
 * -3.94640898771 \times 10^6 \\
 * 3.36625903242 \times 10^6 \\
 * 3.69915069189 \times 10^6 \\
 * \end{array}
 * \right)
 * \f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{c}
 * -2.49064958399 \times 10^6 \\
 * -4.65869968229 \times 10^6 \\
 * 3.56206411337 \times 10^6 \\
 * \end{array}
 * \right)
 * \f$
 * </td></tr>
 * <tr><td>      \f$d^{ab}\f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{ccc}
 * 0.1121397  & 0.3308421  & -0.1802193 \\
 * 0.3308421  & 0.2177940  &  0.1537258 \\
 * -0.1802193  & 0.1537258  & -0.3299337 \\
 * \end{array}
 * \right)
 * \f$
 * </td><td>
 * \f$
 * \left(
 * \begin{array}{ccc}
 * -0.3537959  & 0.2734713  & 0.1095458  \\
 * 0.2734713   & 0.0115214  & 0.2049027  \\
 * 0.1095458   & 0.2049027  & 0.3422745  \\
 * \end{array}
 * \right)
 * \f$
 * </td></tr>
 * <tr><td>     type</td><td>\c LALDETECTORTYPE_IFODIFF</td><td>\c LALDETECTORTYPE_IFODIFF </td></tr>
 * <tr><td>     name</td><td>TAMA_300</td><td>CIT_40 </td></tr>
 * <tr><td>     \f$(\lambda,\beta,h)\f$
 * </td><td>\f$(139^\circ32'9"\!\!.8,35^\circ40'35"\!\!.6, 90\,\textrm{m})\f$
 * </td><td>\f$(-118^\circ\!\!.13,34^\circ\!\!.17, 0\,\textrm{m})\f$
 * </td></tr>
 * <tr><td>      \f$({\mathcal{A}}_X,\zeta_X)\f$
 * </td><td>\f$( 0,          270^\circ)\f$
 * </td><td>\f$( 0,          180^\circ)\f$
 * </td></tr>
 * <tr><td>      \f$({\mathcal{A}}_Y,\zeta_Y)\f$
 * </td><td>\f$( 0,         180^\circ)\f$
 * </td><td>\f$(0,          90^\circ)\f$
 * </td></tr>
 * <tr><td>     \f$(L_X/2,L_Y/2)\f$
 * </td><td>\f$(150\,\textrm{m},150\,\textrm{m})\f$
 * </td><td>\f$(20\,\textrm{m},20\,\textrm{m})\f$
 * </td></tr>
 * </table>
 */
/*@{*/

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

/** expands to constant c of detector d */
#define LAL_DETECTOR_CONSTANT(d,c) LAL_XCAT(LAL_XCAT(LAL_,d),LAL_XCAT(_,c))

/** initializer for detector location vector */
#define LAL_DETECTOR_LOCATION(d) \
{ \
  LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_X_SI),\
  LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_Y_SI),\
  LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_Z_SI) \
}

/** expands to component c (X,Y,Z) of arm X of detector d */
#define LAL_ARM_X(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(ARM_X_DIRECTION_,c))

/** expands to component c (X,Y,Z) of arm Y of detector d */
#define LAL_ARM_Y(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(ARM_Y_DIRECTION_,c))

/** expands to component c (X,Y,Z) of axis of detector d */
#define LAL_AXIS(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(AXIS_DIRECTION_,c))

/** expands to a 3x3 matix initializer for the response for IFODIFF detector d */
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

/** expands to a 3x3 matix initializer for the response for IFOCOMM detector d */
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

/** expands to a 3x3 matix initializer for the response for IFOXARM detector d */
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

/** expands to a 3x3 matix initializer for the response for IFOYARM detector d */
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

/** expands to a 3x3 matix initializer for the response for CYLBAR detector d */
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

#define LAL_FR_STREAM_DETECTOR_STRUCT(d) \
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
  LAL_FR_STREAM_DETECTOR_STRUCT(d)      \
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


/** UNDOCUMENTED */
LALDetector * XLALCreateDetector( LALDetector *detector,
    const LALFrDetector *frDetector, LALDetectorType type )

{
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
    XLAL_ERROR_NULL( XLAL_ENOMEM );

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
      XLAL_ERROR_NULL( XLAL_EINVAL );
  } /* switch (type) */

  detector->frDetector = *frDetector;
  detector->type = type;
  return detector;
}

/**
 * DEPRECATED.
 * \deprecated Use XLALCreateDetector() instead.
 */
void LALCreateDetector( LALStatus             *status,
                        LALDetector           *output,
                        const LALFrDetector   *input,
                        const LALDetectorType  type )
{
  INITSTATUS(status);

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
/*@}*/
