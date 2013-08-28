/*
*  Copyright (C) 2010 Craig Robinson
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


#include <lal/LALInspiral.h>


#ifndef _LALEOBNR_H
#define _LALEOBNR_H

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define LALEOB_MAX_MULTIPOLE 8

/**
 * \defgroup LALEOBNRv2Waveform_h Header LALEOBNRv2Waveform.h
 * \ingroup pkg_CBC_NEW
 * \author Craig Robinson
 *
 * \brief %Header containing most of the structures and prototypes which are
 * used in the generation of the EOBNRv2 waveform.
 */
/*@{*/

/**
 * Structure containing the coefficients for EOBNRv2 A potential function.
 * The elements in the structure are labelled as follows:
 * aN, where a is denotes whether the parameter is in the numerator (n)
 * or denominator (d); and N is the power of r which will multiply this
 * coefficient. For example, the coefficient of r^5 in the numerator
 * will be called n5.
 */
typedef struct
tagEOBACoefficients
{
  REAL8 n4;
  REAL8 n5;
  REAL8 d0;
  REAL8 d1;
  REAL8 d2;
  REAL8 d3;
  REAL8 d4;
  REAL8 d5;
}
EOBACoefficients;

/**
 * Structure containing the coefficients for calculating the factorized
 * waveform. The coefficients are precomputed in the function
 * XLALCalcFacWaveformCoefficients()
 */
typedef struct
tagFacWaveformCoeffs
{
  REAL8 delta22vh3;
  REAL8 delta22vh6;
  REAL8 delta22vh8;
  REAL8 delta22vh9;
  REAL8 delta22v5;

  REAL8 rho22v2;
  REAL8 rho22v3;
  REAL8 rho22v4;
  REAL8 rho22v5;
  REAL8 rho22v6;
  REAL8 rho22v6l;
  REAL8 rho22v7;
  REAL8 rho22v8;
  REAL8 rho22v8l;
  REAL8 rho22v10;
  REAL8 rho22v10l;

  REAL8 delta21vh3;
  REAL8 delta21vh6;
  REAL8 delta21vh7;
  REAL8 delta21vh9;
  REAL8 delta21v5;
  REAL8 delta21v7;

  REAL8 rho21v1;
  REAL8 rho21v2;
  REAL8 rho21v3;
  REAL8 rho21v4;
  REAL8 rho21v5;
  REAL8 rho21v6;
  REAL8 rho21v6l;
  REAL8 rho21v7;
  REAL8 rho21v7l;
  REAL8 rho21v8;
  REAL8 rho21v8l;
  REAL8 rho21v10;
  REAL8 rho21v10l;

  REAL8 delta33vh3;
  REAL8 delta33vh6;
  REAL8 delta33vh9;
  REAL8 delta33v5;
  REAL8 delta33v7;

  REAL8 rho33v2;
  REAL8 rho33v3;
  REAL8 rho33v4;
  REAL8 rho33v5;
  REAL8 rho33v6;
  REAL8 rho33v6l;
  REAL8 rho33v7;
  REAL8 rho33v8;
  REAL8 rho33v8l;

  REAL8 delta32vh3;
  REAL8 delta32vh4;
  REAL8 delta32vh6;
  REAL8 delta32vh9;

  REAL8 rho32v;
  REAL8 rho32v2;
  REAL8 rho32v3;
  REAL8 rho32v4;
  REAL8 rho32v5;
  REAL8 rho32v6;
  REAL8 rho32v6l;
  REAL8 rho32v8;
  REAL8 rho32v8l;

  REAL8 delta31vh3;
  REAL8 delta31vh6;
  REAL8 delta31vh7;
  REAL8 delta31vh9;
  REAL8 delta31v5;

  REAL8 rho31v2;
  REAL8 rho31v3;
  REAL8 rho31v4;
  REAL8 rho31v5;
  REAL8 rho31v6;
  REAL8 rho31v6l;
  REAL8 rho31v7;
  REAL8 rho31v8;
  REAL8 rho31v8l;

  REAL8 delta44vh3;
  REAL8 delta44vh6;
  REAL8 delta44v5;

  REAL8 rho44v2;
  REAL8 rho44v3;
  REAL8 rho44v4;
  REAL8 rho44v5;
  REAL8 rho44v6;
  REAL8 rho44v6l;

  REAL8 delta43vh3;
  REAL8 delta43vh4;
  REAL8 delta43vh6;

  REAL8 rho43v;
  REAL8 rho43v2;
  REAL8 rho43v4;
  REAL8 rho43v5;
  REAL8 rho43v6;
  REAL8 rho43v6l;

  REAL8 delta42vh3;
  REAL8 delta42vh6;

  REAL8 rho42v2;
  REAL8 rho42v3;
  REAL8 rho42v4;
  REAL8 rho42v5;
  REAL8 rho42v6;
  REAL8 rho42v6l;

  REAL8 delta41vh3;
  REAL8 delta41vh4;
  REAL8 delta41vh6;

  REAL8 rho41v;
  REAL8 rho41v2;
  REAL8 rho41v4;
  REAL8 rho41v5;
  REAL8 rho41v6;
  REAL8 rho41v6l;

  REAL8 delta55vh3;
  REAL8 delta55v5;
  REAL8 rho55v2;
  REAL8 rho55v3;
  REAL8 rho55v4;
  REAL8 rho55v5;
  REAL8 rho55v6;

  REAL8 delta54vh3;
  REAL8 delta54vh4;
  REAL8 rho54v2;
  REAL8 rho54v3;
  REAL8 rho54v4;

  REAL8 delta53vh3;
  REAL8 rho53v2;
  REAL8 rho53v3;
  REAL8 rho53v4;
  REAL8 rho53v5;

  REAL8 delta52vh3;
  REAL8 delta52vh4;
  REAL8 rho52v2;
  REAL8 rho52v3;
  REAL8 rho52v4;

  REAL8 delta51vh3;
  REAL8 rho51v2;
  REAL8 rho51v3;
  REAL8 rho51v4;
  REAL8 rho51v5;

  REAL8 delta66vh3;
  REAL8 rho66v2;
  REAL8 rho66v3;
  REAL8 rho66v4;

  REAL8 delta65vh3;
  REAL8 rho65v2;
  REAL8 rho65v3;

  REAL8 delta64vh3;
  REAL8 rho64v2;
  REAL8 rho64v3;
  REAL8 rho64v4;

  REAL8 delta63vh3;
  REAL8 rho63v2;
  REAL8 rho63v3;

  REAL8 delta62vh3;
  REAL8 rho62v2;
  REAL8 rho62v3;
  REAL8 rho62v4;

  REAL8 delta61vh3;
  REAL8 rho61v2;
  REAL8 rho61v3;

  REAL8 delta77vh3;
  REAL8 rho77v2;
  REAL8 rho77v3;

  REAL8 rho76v2;

  REAL8 delta75vh3;
  REAL8 rho75v2;
  REAL8 rho75v3;

  REAL8 rho74v2;

  REAL8 delta73vh3;
  REAL8 rho73v2;
  REAL8 rho73v3;

  REAL8 rho72v2;

  REAL8 delta71vh3;
  REAL8 rho71v2;
  REAL8 rho71v3;

  REAL8 rho88v2;
  REAL8 rho87v2;
  REAL8 rho86v2;
  REAL8 rho85v2;
  REAL8 rho84v2;
  REAL8 rho83v2;
  REAL8 rho82v2;
  REAL8 rho81v2;
}
FacWaveformCoeffs;

/**
 * Structure containing all the terms of the Newtonian multipole which
 * are constant over the course of the evolution, and can therefore be
 * pre-computed. They are stored in a two-dimensional array, which is
 * indexed as values[l][m]. Since m has to be \<= l, this structure
 * is larger than it needs to be; but it makes the coding a bit neater...
 */
typedef
struct tagNewtonMultipolePrefixes
{
  COMPLEX16 values[LALEOB_MAX_MULTIPOLE+1][LALEOB_MAX_MULTIPOLE+1];
}
NewtonMultipolePrefixes;

/**
 * Structure containing all the parameters needed for the EOB waveform.
 * It contains eta, the pre-computed parameters for the A potential function,
 * and the pre-computed parameters for the factorized waveform
 */
typedef
struct tagEOBParams
{
  REAL8 eta;
  REAL8 omega;
  REAL8 m1;
  REAL8 m2;
  EOBACoefficients        *aCoeffs;
  FacWaveformCoeffs       *hCoeffs;
  EOBNonQCCoeffs          *nqcCoeffs;
  NewtonMultipolePrefixes *prefixes;
}
EOBParams;


/*---------------------------------------------------------------- */
/* Functions used in calculating the Newtonian multipolar waveform */
/* Defined in Damour et al Phys.Rev.D79:064004,2009 */
int
XLALCalculateNewtonianMultipole(
                            COMPLEX16 *multipole,
                            REAL8 x,
                            REAL8 r,
                            REAL8 phi,
                            UINT4  l,
                            INT4  m,
                            EOBParams    *params
                            );

REAL8
XLALInspiralFactorizedFlux(
                      REAL8Vector           *values,
                      const REAL8           omega,
                      EOBParams             *ak,
                      const INT4             lMax
                     );

INT4
XLALGetFactorizedWaveform(
                      COMPLEX16             *hlm,
                      REAL8Vector           *values,
                      const REAL8           v,
                      const INT4            l,
                      const INT4            m,
                      EOBParams             *params
                     );

int XLALGetCalibratedNQCCoeffs( EOBNonQCCoeffs *coeffs,
                                INT4            l,
                                INT4            m,
                                REAL8           eta
                                );

int  XLALEOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,
                      REAL8Vector    * restrict values,
                      const REAL8               omega,
                      EOBNonQCCoeffs * restrict coeffs
                     );

int XLALCalculateNQCCoefficients(
                 REAL8Vector    * restrict amplitude,
                 REAL8Vector    * restrict phase,
                 REAL8Vector    * restrict q1,
                 REAL8Vector    * restrict q2,
                 REAL8Vector    * restrict q3,
                 REAL8Vector    * restrict p1,
                 REAL8Vector    * restrict p2,
                 INT4                      l,
                 INT4                      m,
                 REAL8                     timePeak,
                 REAL8                     deltaT,
                 REAL8                     eta,
                 EOBNonQCCoeffs * restrict coeffs );

REAL8 XLALGetNRPeakDeltaT( INT4 l, INT4 m, REAL8 eta );

int XLALCalcFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs,
          const REAL8               eta
          );

int XLALModifyFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs,
          const REAL8               eta
          );

int XLALComputeNewtonMultipolePrefixes(
                NewtonMultipolePrefixes *prefix,
                const REAL8             m1,
                const REAL8             m2 );

/*@}*/ /* end:LALEOBNRv2Waveform_h */

#ifdef  __cplusplus
}
#endif

#endif /* _LALEOBNR_H */
