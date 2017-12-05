#ifndef _LALSimIMRCalculateSpinPrecEOBHCoeffs_C
#define _LALSimIMRCalculateSpinPrecEOBHCoeffs_C

/**
 * \author Craig Robinson, Yi Pan, Stas Babak, Prayush Kumar, Andrea Taracchini
 *
 * This function was originally part of LALSimIMRSpinEOBHamiltonianPrec.c,
 * and moved here during the development of v3_opt.  Function relocation
 * implemented by R. Devine, Z. Etienne, D. Buch, and T. Knowles.  In comments,
 * R.H. refers to Roland Hass.
 */

#include <stdio.h>
#include <math.h>
#include <LALSimIMRSpinEOB.h>

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static int XLALSimIMRCalculateSpinPrecEOBHCoeffs(
        SpinEOBHCoeffs *coeffs,
        const REAL8    eta,
        const REAL8    a,
        const UINT4    SpinAlignedEOBversion
        );


/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 *
 * This function is used to calculate some coefficients which will be used in the
 * spinning EOB Hamiltonian. It takes the following inputs:
 *
 * coeffs - a (non-null) pointer to a SpinEOBParams structure. This will be populated
 * with the output.
 * eta - the symmetric mass ratio.
 * sigmaKerr - the spin of the effective Kerr background (a combination of the individual spins).
 *
 * If all goes well, the function will return XLAL_SUCCESS. Otherwise, XLAL_FAILURE is returned.
 */
static int XLALSimIMRCalculateSpinPrecEOBHCoeffs(
        SpinEOBHCoeffs *coeffs, /**<< OUTPUT, EOB parameters including pre-computed coefficients */
        const REAL8    eta,     /**<< symmetric mass ratio */
        const REAL8    a,       /**<< Normalized deformed Kerr spin */
        const UINT4    SpinAlignedEOBversion  /**<< 1 for SEOBNRv1; 2 for SEOBNRv2 */
        )
{

  REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
  REAL8 m1PlusEtaKK;

  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }

  coeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;

  const int debugPK = 0;
  if( debugPK )
  {
    XLAL_PRINT_INFO("In XLALSimIMRCalculateSpinPrecEOBHCoeffs: SpinAlignedEOBversion = %d,%d\n",
        (int) SpinAlignedEOBversion, (int) coeffs->SpinAlignedEOBversion );
    fflush( NULL );
  }
  /* Constants are fits taken from PRD 86, 024011 (2012) Eq. 37 */
  static const REAL8 c0  = 1.4467; /* needed to get the correct self-force results */
  static const REAL8 c1  = -1.7152360250654402;
  static const REAL8 c2  = -3.246255899738242;

  static const REAL8 c20  = 1.712;
  static const REAL8 c21  = -1.803949138004582;
  static const REAL8 c22  = -39.77229225266885;
  static const REAL8 c23  = 103.16588921239249;

  static const REAL8 third = 1./3.;
  static const REAL8 fifth = 1./5.;
  static const REAL8 ln2 = 0.6931471805599453094172321214581765680755; // log(2.)

  // RH: this assumes that SpinAlignedEOBversion is either 1 or 2
  // RH: the ifthenelse macros return their ifvalue if cond>=0 (specifically
  // the positive sign is set) and elsevalue otherwise. So:
  // RH: 1.5-SpinAlignedEOBversion is positive for SpinAlignedEOBversion==1 and
  // RH: negative for SpinAlignedEOBversion==2
  // RH: SpinAlignedEOBversion-1.5 is reversed

  // RH: TODO check if b3 can ever be non-zero. If not remove all terms using b3.
  coeffs->b3  = 0.;
  coeffs->bb3 = 0.;
#define ifthenelse(cond, ifvalue, elsevalue) ((elsevalue) + (0.5 + copysign(0.5, cond)) * ((ifvalue)-(elsevalue)))
  coeffs->KK = KK = ifthenelse(1.5-SpinAlignedEOBversion,
                               c0 + c1*eta + c2*eta*eta,
                               c20 + c21*eta + c22*(eta*eta) + c23*(eta*eta)*eta);
  m1PlusEtaKK = -1. + eta*KK;
  const REAL8 invm1PlusEtaKK = 1./m1PlusEtaKK;
  /* Eqs. 5.77 - 5.81 of PRD 81, 084024 (2010) */
  coeffs->k0 = k0 = KK*(m1PlusEtaKK - 1.);
  coeffs->k1 = k1 = - 2.*(k0 + KK)*m1PlusEtaKK;
  k1p2= k1*k1;
  k1p3= k1*k1p2;
  coeffs->k2 = k2 = (k1 * (k1 - 4.*m1PlusEtaKK)) * 0.5 - a*a*k0*m1PlusEtaKK*m1PlusEtaKK;
  coeffs->k3 = k3 = -(k1*k1)*k1 * third + k1*k2 + (k1*k1)*m1PlusEtaKK - 2.*(k2 - m1PlusEtaKK)*m1PlusEtaKK - a*a*k1*(m1PlusEtaKK*m1PlusEtaKK);
  coeffs->k4 = k4 = ((24./96.)*(k1*k1)*(k1*k1) - (96./96.)*(k1*k1)*k2 + (48./96.)*k2*k2 - (64./96.)*(k1*k1)*k1*m1PlusEtaKK
      + (48./96.)*(a*a)*(k1*k1 - 2.*k2)*(m1PlusEtaKK*m1PlusEtaKK) +
      (96./96.)*k1*(k3 + 2.*k2*m1PlusEtaKK) - m1PlusEtaKK*((192./96.)*k3 + m1PlusEtaKK*(-(3008./96.) + (123./96.)*LAL_PI*LAL_PI)));
#define ifthenelsezero(cond, ifvalue) ((0.5 + copysign(0.5, cond)) * (ifvalue))
  coeffs->k5 = k5 = ifthenelsezero(SpinAlignedEOBversion-1.5,
                     m1PlusEtaKK*m1PlusEtaKK
                     * (-4237./60.+128./5.*LAL_GAMMA+2275.*LAL_PI*LAL_PI/512.
                     - third*(a*a)*(k1p3-3.*(k1*k2)+3.*k3)
                     - ((k1p3*k1p2)-5.*(k1p3*k2)+5.*k1*k2*k2+5.*k1p2*k3-5.*k2*k3-5.*k1*k4)*fifth*invm1PlusEtaKK*invm1PlusEtaKK
                     + ((k1p2*k1p2)-4.*(k1p2*k2)+2.*k2*k2+4.*k1*k3-4.*k4)*0.5*invm1PlusEtaKK+(256./5.)*ln2)
                    );
  coeffs->k5l= k5l= ifthenelsezero(SpinAlignedEOBversion-1.5, (m1PlusEtaKK*m1PlusEtaKK) * (64./5.));

  /* Now calibrated parameters for spin models */
  coeffs->d1 = ifthenelsezero(1.5-SpinAlignedEOBversion, -69.5);
  coeffs->d1v2 = ifthenelsezero(SpinAlignedEOBversion-1.5, -74.71 - 156.*eta + 627.5*eta*eta);
  coeffs->dheffSS = ifthenelsezero(1.5-SpinAlignedEOBversion, 2.75);
  coeffs->dheffSSv2 = ifthenelsezero(SpinAlignedEOBversion-1.5, 8.127 - 154.2*eta + 830.8*eta*eta);

  return XLAL_SUCCESS;
}

#endif // _LALSimIMRCalculateSpinPrecEOBHCoeffs_C
