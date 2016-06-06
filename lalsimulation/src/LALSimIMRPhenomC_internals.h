#ifndef _LALSIM_IMR_PHENOMC_INTERNALS_H
#define _LALSIM_IMR_PHENOMC_INTERNALS_H

/*
 * Copyright (C) 2012 Prayush Kumar, Frank Ohme
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


/* The paper refered to here as the Main paper, is Phys. Rev. D 82, 064016 (2010)
 * */

#include <math.h>
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/StringInput.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>


/*********************************************************************/
/* This structure stores the PN coefficients used to calculate flux  */
/* and waveform amplitude, and Fourier phase. It also stores some    */
/* frequently used expressions which are constant during waveform    */
/* generation.                                                       */
/*********************************************************************/

// MP: could we move this into the header file?
typedef struct tagBBHPhenomCParams{
  REAL8 piM;
  REAL8 m_sec;

  REAL8 fmin;
  REAL8 fCut;
  REAL8 df;

  REAL8 f0;
  REAL8 f1;
  REAL8 f2;
  REAL8 d0;
  REAL8 d1;
  REAL8 d2;

  REAL8 afin;
  REAL8 fRingDown;
  REAL8 MfRingDown;
  REAL8 Qual;

  REAL8 pfaN;
  REAL8 pfa1;
  REAL8 pfa2;
  REAL8 pfa3;
  REAL8 pfa4;
  REAL8 pfa5;
  REAL8 pfa6;
  REAL8 pfa6log;
  REAL8 pfa7;

  REAL8 xdotaN;
  REAL8 xdota2;
  REAL8 xdota3;
  REAL8 xdota4;
  REAL8 xdota5;
  REAL8 xdota6;
  REAL8 xdota6log;
  REAL8 xdota7;

  REAL8 AN;
  REAL8 A2;
  REAL8 A3;
  REAL8 A4;
  REAL8 A5;
  REAL8 A5imag;
  REAL8 A6;
  REAL8 A6log;
  REAL8 A6imag;

  REAL8 a1;
  REAL8 a2;
  REAL8 a3;
  REAL8 a4;
  REAL8 a5;
  REAL8 a6;
  REAL8 g1;
  REAL8 del1;
  REAL8 del2;
  REAL8 b1;
  REAL8 b2;
}
BBHPhenomCParams;

/**
 *
 * private function prototypes; all internal functions use solar masses.
 *
 */

static BBHPhenomCParams *ComputeIMRPhenomCParamsSPA( const REAL8 m1, const REAL8 m2, const REAL8 chi, const LALSimInspiralTestGRParam *extraParams  );
static BBHPhenomCParams *ComputeIMRPhenomCParams( const REAL8 m1, const REAL8 m2, const REAL8 chi , const LALSimInspiralTestGRParam *extraParams );
static REAL8 wPlus( const REAL8 f, const REAL8 f0, const REAL8 d, const BBHPhenomCParams *params );
static REAL8 wMinus( const REAL8 f, const REAL8 f0, const REAL8 d, const BBHPhenomCParams *params );

UNUSED static size_t NextPow2_PC(const size_t n);
static REAL8 IMRPhenomCGeneratePhasePM( REAL8 f, REAL8 eta, const BBHPhenomCParams *params );
UNUSED static int IMRPhenomCGenerateAmpPhase( REAL8 *amplitude, REAL8 *phasing, REAL8 f, REAL8 eta, const BBHPhenomCParams *params);

#endif	// of #ifndef _LALSIM_IMR_PHENOMC_INTERNALS_H
