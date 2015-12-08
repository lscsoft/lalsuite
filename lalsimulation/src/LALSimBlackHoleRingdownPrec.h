/*
 * Copyright (C) 2011 J. Creighton
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

#ifndef _LALSIMBLACKHOLERINGDOWNPREC_H
#define _LALSIMBLACKHOLERINGDOWNPREC_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * Computes the final mass and spin of the black hole resulting from merger.
 * They are given by fittings of NR simulations results. Specifically,
 * for EOBNR, Table I of Buonanno et al. PRD76, 104049;
 * for EOBNRv2 and EOBNRv2HM, Eqs. 29a and 29b of Pan et al. PRD84, 124052;
 * for SEOBNRv1, Eq. 8 of Tichy and Marronetti PRD78, 081501 and
 * Eqs. 1 and 3 of Barausse and Rezzolla ApJ704, L40.
 */
INT4 XLALSimIMREOBFinalMassSpinPrec(
	REAL8    *finalMass,  /**<< OUTPUT, the final mass (scaled by original total mass) */
	REAL8    *finalSpin,  /**<< OUTPUT, the final spin (scaled by final mass) */
  const REAL8     mass1,      /**<< The mass of the 1st component of the system */
  const REAL8     mass2,      /**<< The mass of the 2nd component of the system */
  const REAL8     spin1[3],   /**<< The spin of the 1st object; only needed for spin waveforms */
  const REAL8     spin2[3],   /**<< The spin of the 2nd object; only needed for spin waveforms */
  Approximant     approximant /**<< The waveform approximant being used */
);

/**
 * This function generates the quasinormal mode frequencies for a black
 * hole ringdown. At present, this function works for the 22, 21, 33, 44
 * and 55 modes, and includes 8 overtones. The final frequencies are
 * computed by interpolating the data found on the webpage of
 * Vitor Cardoso, http://centra.ist.utl.pt/~vitor/?page=ringdown
 * In this page, frequecy data are given for positive final spins only.
 * For a negative final spin chi<0 case, the (l,m) mode frequency is given by
 * the (l,-m) mode frequency of the positive final spin -chi case.
 */
INT4 XLALSimIMREOBGenerateQNMFreqV2Prec(
  COMPLEX16Vector *modefreqs, /**<< OUTPUT, complex freqs of overtones (scaled by total mass) */
  const REAL8      mass1,     /**<< The mass of the 1st component (in Solar masses) */
  const REAL8      mass2,     /**<< The mass of the 2nd component (in Solar masses) */
  const REAL8      spin1[3],  /**<< The spin of the 1st object; only needed for spin waveforms */
  const REAL8      spin2[3],  /**<< The spin of the 2nd object; only needed for spin waveforms */
  UINT4            l,         /**<< The l value of the mode in question */
  INT4             m,         /**<< The m value of the mode in question */
  UINT4            nmodes,    /**<< The number of overtones that should be included (max 8) */
  Approximant      approximant/**<< The waveform approximant being used */
  );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMBLACKHOLERINGDOWNPREC_H */
