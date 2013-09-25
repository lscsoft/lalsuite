/*
 * Copyright (C) 2007 Bruce Allen, Duncan Brown, Jolien Creighton, Kipp
 * Cannon, Teviet Creighton
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#ifndef _WINDOW_H
#define _WINDOW_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup Window_h Header Window.h
 * \ingroup pkg_window
 * \brief This header file provides routines and structures to create and store window functions (also called a taper,
 * lag window, or apodization function).
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Window.h>
 * \endcode
 *
 * ### Description ###
 *
 * These functions create or destroy a time-domain window function in a vector
 * of specified length.  If you wish to construct a custom window, call
 * <tt>XLALCreateRectangularREAL8Window()</tt> (or the \c REAL4
 * version), then replace the samples inside it with your own, and update the
 * \c sumofsquares and \c sum elements.  If the window function
 * proves useful, consider adding it here so that others can benefit.
 *
 * It is convenient to describe the windows as functions on the normalized
 * domain \f$y \in [-1, 1]\f$.  The window is zero outside this domain.  The
 * window functions defined in this package are as follows.
 *
 * ### Rectangle ###
 *
 * \f{equation}{
 * w(y)
 * = 1.
 * \f}</dd>
 *
 * ### Hann ###
 *
 * \f{equation}{
 * w(y)
 * = \cos^2 \frac{\pi}{2} y.
 * \f}</dd>
 *
 * ### Welch ###
 *
 * \f{equation}{
 * w(y)
 * = 1 - y^2.
 * \f}</dd>
 *
 * ### Bartlett ###
 *
 * \f{equation}{
 * w(y)
 * = 1 - |y|.
 * \f}</dd>
 *
 * ### Parzen ###
 *
 * \f{equation}{
 * w(y)
 * = \left\{ \begin{array}{ll}
 * 1 - 6 y^2 (1 - |y|) & |y| \leq 1 / 2, \\
 * 2 (1 - |y|)^3 & |y| > 1 / 2.
 * \end{array}\right.
 * \f}</dd>
 *
 * ### Papoulis ###
 *
 * \f{equation}{
 * w(y)
 * = \frac{1}{\pi} \sin \pi |y| + (1 - |y|) \cos \pi |y|.
 * \f}</dd>
 *
 * ### Hamming ###
 *
 * \f{equation}{
 * w(y)
 * = 0.08 + 0.92 \cos^{2} \frac{\pi}{2} y.
 * \f}
 * This is the same as the Hann window, but with an additional DC bias, or
 * &quot;foot,&quot; of 0.08.</dd>
 *
 * ### Kaiser ###
 *
 * \f{equation}{
 * w(y)
 * = I_0 \left( \beta \sqrt{1-y^2} \right) / I_0(\beta),
 * \f}
 * where \f$I_0(x)\f$ is the \f$0\f$th order, modified Bessel function of the first
 * kind.  The shape parameter \f$\beta \in [0, \infty]\f$ sets the sharpness of
 * the central peak.  \f$\beta = 0\f$ yields the rectangle window, \f$\beta
 * \rightarrow \infty\f$ yields a \f$\delta\f$ function with a single non-zero
 * sample in the middle.  This window is difficult to compute for large
 * \f$\beta\f$, and an asymptotic approximation is used for \f$\beta \ge 705\f$.  A
 * linearly-interpolated transition occurs between \f$\beta = 695\f$ and \f$\beta =
 * 705\f$.  Finite-difference derivatives of the window with respect to \f$\beta\f$
 * are unlikely to work well in this regime.</dd>
 *
 * ### Creighton ###
 *
 * \f{equation}{
 * w(y)
 * = \exp \left( -\beta \frac{y^2}{1 - y^2} \right).
 * \f}
 * This window function is based on a fairly standard \f$C_{\infty}\f$ test
 * function used in distribution theory, e.g.\ <em>Green's Functions and
 * Boundary Value Problems</em> \cite stakgold79, by Stakgold.  The shape parameter
 * \f$\beta \in [0, \infty]\f$ sets the sharpness of the central peak.  \f$\beta =
 * 0\f$ yields the rectangle window, \f$\beta \rightarrow \infty\f$ yields a
 * \f$\delta\f$ function with a single non-zero sample in the middle.</dd>
 *
 * ### Tukey ###
 *
 * \f{equation}{
 * w(y)
 * = \left\{ \begin{array}{ll}
 * \sin^2 \left[ \frac{\pi}{2} (|y| - 1) / \beta \right] & |y| \geq 1 -
 * \beta, \\
 * 1 & |y| < 1 - \beta.
 * \end{array} \right.
 * \f}
 * The shape parameter \f$\beta \in [0, 1]\f$ sets what fraction of the window is
 * spanned by the tapers.  \f$\beta = 0\f$ yields the rectangle window, \f$\beta =
 * 1\f$ yields the Hann window.</dd>
 *
 * ### Gauss ###
 *
 * \f{equation}{
 * w(y)
 * = \exp \left( -\frac{1}{2} \beta^{2} y^{2} \right).
 * \f}
 * The shape parameter \f$\beta \in [0, \infty]\f$ sets the sharpness of the
 * central peak.  \f$\beta = 0\f$ yields the rectangle window, \f$\beta \rightarrow \infty\f$ yields a \f$\delta\f$ function with a single non-zero sample in the
 * middle.
 *
 * These window functions are shown in \figref{window_t}, showing various windows as functions of the normalized
 * independend variable \f$y\f$, choosing \f$\beta = 6\f$ for the Kaiser window, \f$\beta = 2\f$ for the Creighton window,
 * \f$\beta = 0.5\f$ for the Tukey window, and \f$\beta = 3\f$ for the Gauss window.
 *
 * \figure{window_t,pdf,0.6,Various windows as functions of the normalized independend variable y}
 *
 * For a vector of length \f$L\f$ (an integer), the mapping from integer array
 * index \f$i\f$ to normalized co-ordinate \f$y\f$ is
 * \f{equation}{
 * y(i)
 * = \left\{ \begin{array}{ll}
 * 0 & L \le 1, \\
 * 2 i / (L - 1) - 1 & L > 1,
 * \end{array} \right.
 * \f}
 * where \f$0 \le i < L\f$, and floating-point division is used.  This agrees with
 * J. G. Proakis and D. G. Manolakis, <em>Digital Signal Processing</em>
 * \cite pm95, and \c MatLab.  The first sample is \f$y = -1\f$, the last
 * sample is \f$y = +1\f$.  For odd-lengthed vectors, the middle sample is \f$y =
 * 0\f$, while for even-lengthed vectors \f$y = 0\f$ occurs half-way between the two
 * middle samples.  Substituting \f$y(i)\f$ into the definitions of the window
 * functions above yields \f$w(i)\f$, the value of the window function at the
 * integer sample \f$i\f$.
 *
 * The Fourier transforms of the windows are shown as functions of \f$1 / y\f$ in
 * \figref{window_f}, showing frequency behaviour of various windows as functions
 * of the inverse of the normalized independend variable \f$y\f$, choosing \f$\beta = 6\f$
 * for the Kaiser window, \f$\beta = 2\f$ for the Creighton window, \f$\beta = 0.5\f$ for
 * the Tukey window, and \f$\beta = 3\f$ for the Gauss window.
 *
 * \figure{window_f,pdf,0.6,Frequency behaviour of various windows as functions of the inverse of the normalized independend variable y}
 *
 * Since the Fourier transform of windowed data is the Fourier transform of
 * the data convolved with the Fourier transform of the window,
 * \figref{window_f} is the major guideline for selecting a window.  One
 * can see that windows with a narrow central lobe tend to have higher
 * sidelobes, and windows which suppress their low-order sidelobes tend to
 * have more power in the high-order sidelobes.  The choice of window thus
 * depends on whether one is trying to resolve nearby spectral features of
 * comparable magnitude (suggesting a rectangular or a Welch window), to
 * reduce spectral bias and low-order sidelobes (a Hamming or Kaiser window),
 * or to measure a broad spectrum with a large dynamical range (a Creighton or
 * a Papoulis window).
 *
 */
/*@{*/


/**
 * Structure for storing REAL4 window function data, providing storage for a sequence of samples
 * as well as metadata about the window such as the sum-of-squarse of the samples
 */
typedef struct tagREAL4Window {
        REAL4Sequence *data;		/**< The window function samples */
	REAL8          sumofsquares;	/**< The sum of the squares of the window function samples */
        REAL8          sum;		/**< The sum of the window function samples */
} REAL4Window;


/**
 * Structure for storing REAL8 window function data, providing storage for a sequence of samples
 * as well as metadata about the window such as the sum-of-squarse of the samples
 */
typedef struct tagREAL8Window {
	REAL8Sequence *data;		/**< The window function samples */
	REAL8          sumofsquares;	/**< The sum of the squares of the window function samples */
	REAL8          sum;		/**< The sum of the window function samples */
} REAL8Window;


#ifdef SWIG   // SWIG interface directives
SWIGLAL(ACQUIRES_OWNERSHIP(REAL4Sequence*, sequence));
SWIGLAL(ACQUIRES_OWNERSHIP(REAL8Sequence*, sequence));
#endif
REAL4Window *XLALCreateREAL4WindowFromSequence(REAL4Sequence *sequence);
REAL8Window *XLALCreateREAL8WindowFromSequence(REAL8Sequence *sequence);
#ifdef SWIG   // SWIG interface directives
SWIGLAL_CLEAR(ACQUIRES_OWNERSHIP(REAL4Sequence*, sequence));
SWIGLAL_CLEAR(ACQUIRES_OWNERSHIP(REAL8Sequence*, sequence));
#endif

REAL4Window *XLALCreateRectangularREAL4Window(UINT4 length);
REAL4Window *XLALCreateHannREAL4Window(UINT4 length);
REAL4Window *XLALCreateWelchREAL4Window(UINT4 length);
REAL4Window *XLALCreateBartlettREAL4Window(UINT4 length);
REAL4Window *XLALCreateParzenREAL4Window(UINT4 length);
REAL4Window *XLALCreatePapoulisREAL4Window(UINT4 length);
REAL4Window *XLALCreateHammingREAL4Window(UINT4 length);
REAL4Window *XLALCreateKaiserREAL4Window(UINT4 length, REAL4 beta);
REAL4Window *XLALCreateCreightonREAL4Window(UINT4 length, REAL4 beta);
REAL4Window *XLALCreateTukeyREAL4Window(UINT4 length, REAL4 beta);
REAL4Window *XLALCreateGaussREAL4Window(UINT4 length, REAL4 beta);

REAL8Window *XLALCreateRectangularREAL8Window(UINT4 length);
REAL8Window *XLALCreateHannREAL8Window(UINT4 length);
REAL8Window *XLALCreateWelchREAL8Window(UINT4 length);
REAL8Window *XLALCreateBartlettREAL8Window(UINT4 length);
REAL8Window *XLALCreateParzenREAL8Window(UINT4 length);
REAL8Window *XLALCreatePapoulisREAL8Window(UINT4 length);
REAL8Window *XLALCreateHammingREAL8Window(UINT4 length);
REAL8Window *XLALCreateKaiserREAL8Window(UINT4 length, REAL8 beta);
REAL8Window *XLALCreateCreightonREAL8Window(UINT4 length, REAL8 beta);
REAL8Window *XLALCreateTukeyREAL8Window(UINT4 length, REAL8 beta);
REAL8Window *XLALCreateGaussREAL8Window(UINT4 length, REAL8 beta);

void XLALDestroyREAL4Window(REAL4Window *window);
void XLALDestroyREAL8Window(REAL8Window *window);

REAL4Sequence *XLALUnitaryWindowREAL4Sequence(REAL4Sequence *sequence, const REAL4Window *window);
COMPLEX8Sequence *XLALUnitaryWindowCOMPLEX8Sequence(COMPLEX8Sequence *sequence, const REAL4Window *window);
REAL8Sequence *XLALUnitaryWindowREAL8Sequence(REAL8Sequence *sequence, const REAL8Window *window);
COMPLEX16Sequence *XLALUnitaryWindowCOMPLEX16Sequence(COMPLEX16Sequence *sequence, const REAL8Window *window);

REAL8Window *XLALCreateNamedREAL8Window ( const char *windowName, REAL8 beta, UINT4 length );
REAL4Window *XLALCreateNamedREAL4Window ( const char *windowName, REAL8 beta, UINT4 length );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _WINDOW_H */
