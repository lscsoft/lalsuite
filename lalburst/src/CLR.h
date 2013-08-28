/*
*  Copyright (C) 2007 Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: CLR.h
 *
 * Author: Sintes, A. M.
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * CLR.h
 *
 * SYNOPSIS
 * #include "CLR.h"
 *
 * DESCRIPTION
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

/*
 * 2. include-loop protection (see below). Note the naming convention!
 */

#ifndef _CLR_H
#define _CLR_H

/*
 * 3. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \author Sintes, A. M.
 * \defgroup CLR_h		Header CLR.h
 * \ingroup pkg_clremoval
 *
 * Provides routines for finding line harmonics,
 * generating a reference interference signal, and removing
 * all the interference harmonics.
 *
 * ### Synopsis ###
 *
 * \code
 * #include "CLR.h"
 * \endcode
 *
 * The principal of  \c clr is the following:
 *
 * We assume that the interference has the form
 * \anchor e3 \f{equation}{
 * y(t)=\sum_n a_n m(t)^n + \left( a_n m(t)^n\right)^* \ ,
 * \tag{e3}
 * \f}
 * where \f$a_n\f$ are
 * complex amplitudes and  \f$m(t)\f$ is a nearly monochromatic function
 * near
 * a frequency \f$f_0\f$.
 * The idea is to
 * use the information in the different harmonics of the interference
 * to construct a function \f$M(t)\f$ that is  as close a replica as possible of
 * \f$m(t)\f$ and then construct a function close to \f$y(t)\f$ which
 * is subtracted from the output  of the system cancelling the interference.
 * The key is that real gravitational wave signals will not be present
 * with multiple harmonics and that \f$M(t)\f$ is constructed from many
 * frequency bands with independent noise. Hence, \c clr will little
 * affect the statistics of the noise in any one band and
 * any gravitational wave signal masked by the interference can
 * be recovered without any disturbance.
 *
 * We assume that the data produced by the system
 * is just the sum of the interference plus  noise
 * \anchor e6 \f{equation}{
 * x(t)=y(t)+n(t) \ ,
 * \tag{e6}
 * \f}
 * where \f$y(t)\f$ is given by Eq.\eqref{e3} and the noise \f$n(t)\f$ in the
 * detector  is a zero-mean stationary
 * stochastic process.
 * The procedure consists in   defining a set of functions \f$\tilde z_k(\nu)\f$
 * in the frequency domain as
 * \anchor e8 \f{equation}{
 * \tilde z_k(\nu)\equiv \left\{
 * \begin{array}{cc}
 * \tilde x(\nu) & \nu_{ik}<\nu <\nu_{fk}\\
 * 0 & \mbox{elsewhere}\ ,
 * \end{array}
 * \right.
 * \tag{e8}
 * \f}
 * where  \f$(\nu_{ik}, \nu_{fk})\f$ correspond to the upper and lower frequency
 * limits of the harmonics of the interference
 * and \f$k\f$ denotes the harmonic considered.
 * These functions are equivalent to
 * \f{equation}{
 * \tilde z_k(\nu)= a_k \widetilde{m^k}(\nu) +\tilde n_k(\nu) \ ,
 * \f}
 * where \f$ \tilde n_k(\nu)\f$ is the noise in the frequency band of the
 * harmonic considered.  Their inverse  Fourier transforms yield
 * \f{equation}{
 * z_k(t)=a_k m(t)^k +n_k(t) \ .
 * \f}
 * Since  \f$m(t)\f$ is supposed to be a narrow-band function near a frequency \f$f_0\f$,
 * each \f$z_k(t)\f$ is a  narrow-band function near \f$kf_0\f$.
 * Then, we  define
 * \anchor e10a \f{equation}{
 * B_k(t)\equiv \left[ z_k(t)\right]^{1/k}\ ,\tag{e10a}
 * \f}
 * that can be rewritten as
 * \anchor e10 \f{equation}{
 * B_k(t)= (a_k)^{1/k}m(t) \beta_k(t) \ , \qquad
 * \beta_k(t)=\left[ 1+ \frac{n_k(t)}{a_k m(t)^k}\right]^{1/k} \ .
 * \tag{e10}
 * \f}
 * All these  functions, \f$\{B_k(t)\}\f$, are almost monochromatic around the
 * fundamental frequency, \f$f_0\f$, but they differ basically by a certain
 * complex amplitude. These factors, \f$\Gamma_k\f$, can easily be  calculated,
 * and  we can construct a set of functions   \f$\{b_k(t)\}\f$
 * \f{equation}{
 * b_k(t)=\Gamma_k B_k(t)\ ,
 * \f}
 * such that, they all have the same mean value. Then,   \f$M(t)\f$ can be
 * constructed  as a function of all \f$\{b_k(t)\}\f$
 * in such a way
 * that it has the same mean and minimum variance.
 * If
 * \f$M(t)\f$ is linear with \f$\{b_k(t)\}\f$, the statistically the best is
 * \f{equation}{
 * M(t)=\left(\sum_k \frac{b_k(t)}{\textrm{Var}[\beta_k(t)]} \right) {\Big
 * { /}}
 * \left( \sum_k \frac{1}{\mathrm{Var}[\beta_k(t)]}\right) \ ,
 * \f}
 * where
 * \f{equation}{
 * \textrm{Var}[\beta_k(t)]= \frac{\langle n_k(t) n_k(t)^*\rangle}{k^2
 * \vert a_k m(t)^k\vert^2}+ \mbox{corrections} \ .
 * \f}
 * In practice,
 * we  approximate
 * \f{equation}{
 * \vert a_k m(t)^k\vert^2 \approx \vert z_k(t)\vert^2 \ ,
 * \f}
 * and we assume  stationary noise. Therefore,
 * \f{equation}{
 * \langle n_k(t) n_k(t)^*\rangle= \int_{\nu_{ik}}^{\nu_{fk}} S(\nu) d\nu \ ,
 * \f}
 * where \f$S(\nu)\f$ is the power spectral density of the noise.
 *
 * Finally, it only remains to determine the amplitude of the different
 * harmonics, which can be obtained applying a least square method.
 */

/*@{*/


/*
 * 5. Macros. But, note that macros are deprecated.
 */

/*
 * 8. Structure, enum, union, etc., typdefs.
 */

/**\name Error Codes */
/*@{*/
#define CLRH_ENULL 1	/**< Null pointer */
#define CLRH_ESIZE 2	/**< Invalid input size */
#define CLRH_ESZMM 4	/**< Size mismatch */
#define CLRH_EINT  6	/**< Invalid interval */
#define CLRH_ESAME 8	/**< Input/Output data vectors are the same */
#define CLRH_EFREQ 10	/**< Invalid frequency */
/*@}*/

/** \cond DONT_DOXYGEN */
#define CLRH_MSGENULL "Null pointer"
#define CLRH_MSGESIZE "Invalid input size"
#define CLRH_MSGESZMM "Size mismatch"
#define CLRH_MSGEINT  "Invalid interval"
#define CLRH_MSGESAME "Input/Output data vectors are the same"
#define CLRH_MSGEFREQ "Invalid frequency"
/** \endcond */


/* Available CLR Time Vector  */

/**
 * This structure stores the time domain data and other
 * information needed in order to remove all the line interference
 * harmonics.  The fields are:
 *
 * <dl>
 * <dt><tt>UINT4 length</tt></dt><dd>  The number of elements in data.</dd>
 * <dt><tt>REAL4  *data</tt></dt><dd> The (real) time domain data.</dd>
 * <dt><tt>REAL4  deltaT</tt></dt><dd> The sample spacing in time (in seconds).</dd>
 * <dt><tt>REAL4  fLine</tt></dt><dd> The interference fundamental frequency \f$f_0\f$
 * (in Hz), e.g., 60 Hz.</dd>
 * </dl>
 */
typedef struct tagREAL4TVectorCLR{
  UINT4   length;  /**< number of elements in data */
  REAL4   *data;   /**< time domain data */
  REAL4   deltaT;  /**< sample spacing in time (in seconds) */
  REAL4   fLine;   /**< interference fundamental frequency (e.g., 50 Hz) */
} REAL4TVectorCLR;


  /* Available CLR Frequency Vector  */

/**
 * This structure stores the spectrum, \f$\vert\tilde x(\nu)\vert^2\f$,
 * and other information
 * needed to find the location of the line harmonics.
 * The fields are:
 *
 * <dl>
 * <dt><tt>UINT4 length</tt></dt><dd>  The number of elements in data.</dd>
 * <dt><tt>REAL4   *data</tt></dt><dd> The (real)  frequency domain data, e.g.,
 * \f$\vert\tilde x(\nu)\vert^2\f$.</dd>
 * <dt><tt>REAL4   deltaF</tt></dt><dd> The \f$\Delta\f$F offset between samples (in Hz).</dd>
 * <dt><tt>REAL4   fLine</tt></dt><dd> The interference fundamental frequency \f$f_0\f$
 * (in Hz), e.g., 60 Hz.</dd>
 * </dl>
 *
 */
typedef struct tagREAL4FVectorCLR{
  UINT4   length;  /**< number of elements in data */
  REAL4   *data;   /**< frequency domain data */
  REAL4   deltaF;  /**< delta F offset between samples (in Hz)*/
  REAL4   fLine;   /**< interference fundamental frequency (e.g., 50 Hz) */
} REAL4FVectorCLR;


  /*
   * 9. Functions Declarations (i.e., prototypes).
   */

void LALRefInterference ( LALStatus          *status,
			  COMPLEX8Vector     *out, /*  M(t), size n */
			  COMPLEX8Vector     *in,  /* x(nu), size n/2+1 */
			  INT4Vector         *par
			  );

void LALCleanAll ( LALStatus        *status,
		   REAL4Vector      *out,  /* clean data */
		   COMPLEX8Vector   *in2,  /* M(t), ref. interference */
		   REAL4TVectorCLR  *in1   /* x(t), data + information */
		 );


void LALHarmonicFinder (
         LALStatus          *status,
         INT4Vector         *out,   /* harmonic index, bin location (3*l) */
         REAL4FVectorCLR    *in2,   /* |x(f)|^2, data + inform.  */
         INT4Vector         *in1    /* the harmonic index (l) */
		 );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _CLR_H */
