/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RingUtils.h>

static REAL4 ring_spin_factor( REAL4 a )
{
  /* Cardoso's equation from Berti et al. (2008) */
  /* Define 3 parameters for l=m=2, n=0 */
  const REAL4 fparam1 = 1.5251;
  const REAL4 fparam2 = -1.1568;
  const REAL4 fparam3 = 0.1292;

  return fparam1 + fparam2 * pow( 1.0 - a, fparam3);

  /* Echeverria's equation */
  /* return 1.0 - 0.63 * pow( 1.0 - a, 3.0/10.0 ); */
}

static REAL4 ring_quality_fn( REAL4 Q )
{
    return  1.0 + 7.0/(24.0*Q*Q);
}


/* FIXME Is there already a function*/
/* in LAL that will return the sign */
static int SIGN( REAL4 x )
{
    if( x > 0 )

        return 1;

    else if( x < 0 )

        return -1;

    else

        return 0;
}

/*
 *
 * XLAL Routines.
 *
 */



REAL4 XLALBlackHoleRingSpin( REAL4 Q )

{

  /* Cardoso's equation from Berti et al. (2008) */
  /* Define 3 parameters for l=m=2, n=0 */
  const REAL4 qparam1 = 0.7000;
  const REAL4 qparam2 = 1.4187;
  const REAL4 qparam3 = -0.4990;

  return 1.0 - pow( qparam2/(Q - qparam1), -1.0/qparam3 );

  /* Echeverria's equation */
  /* return 1.0 - pow( 2.0/Q, 20.0/9.0 ); */
}


REAL4 XLALBlackHoleRingMass( REAL4 f, REAL4 Q )

{
  const REAL4 c = LAL_C_SI;
  const REAL4 a = XLALBlackHoleRingSpin( Q );
  const REAL4 g = ring_spin_factor( a );
  return (c * c * c * g) / ( LAL_TWOPI * LAL_G_SI * LAL_MSUN_SI * f );
}


REAL4 XLALBlackHoleRingQuality( REAL4 a )

{

  /* Cardoso's equation from Berti et al. (2008) */
  /* Define 3 parameters for l=m=2, n=0 */
  const REAL4 qparam1 = 0.7000;
  const REAL4 qparam2 = 1.4187;
  const REAL4 qparam3 = -0.4990;

  return qparam1 + (qparam2 * pow( (1.0 - a), qparam3 ));

  /* Echeverria's equation */
  /*  return 2.0 * pow( ( 1.0 - a ), -0.45 ); */
}


REAL4 XLALBlackHoleRingFrequency( REAL4 M, REAL4 a )

{
  const REAL4 c = LAL_C_SI;
  const REAL4 g = ring_spin_factor( a );
  return (c * c * c * g) / ( LAL_TWOPI * LAL_G_SI * LAL_MSUN_SI * M );
}

/* Formulas for final mass and spin of a non-spinning binary */
/* Spin equation is not used. See XLALSpinBinaryFinalBHSpin() below. */
/* Buonanno et al arxiv:0706.3732v3 */

REAL4 XLALNonSpinBinaryFinalBHSpin( REAL4 eta )

{
  return sqrt(12.0) * eta - 2.9 * eta *eta;
}


REAL4 XLALNonSpinBinaryFinalBHMass( REAL4 eta, REAL4 mass1, REAL4 mass2 )

{
  return ( 1 + ( sqrt(8.0/9.0) - 1) * eta - 0.498 * eta * eta) * (mass1 + mass2);
}

/* Formulas for final spin of an initially spinning or initially */
/* non-spinning binary. */
/* Barausse & Rezzolla arxiv:0904.2577 */
/* Equations 6, 7, and 10 */

REAL4 XLALSpinBinaryFinalBHSpin( REAL4 eta, REAL4 mass1, REAL4 mass2, REAL4 spin1x, REAL4 spin2x,
   REAL4 spin1y, REAL4 spin2y, REAL4 spin1z, REAL4 spin2z )

{
  REAL4 norml;
  REAL4 cosalpha;
  REAL4 cosbeta;
  REAL4 cosgamma;
  const REAL4 t0 = -2.8904;
  const REAL4 t2 = -3.5171;
  const REAL4 t3 = 2.5763;
  const REAL4 s4 = -0.1229;
  const REAL4 s5 = 0.4537;
  const REAL4 q = mass2 / mass1;
  const REAL4 q2 = q * q;
  const REAL4 q4 = q2 * q2;
  const REAL4 q2sum = 1.0 + q2;
  const REAL4 norma1 = sqrt( spin1x * spin1x + spin1y * spin1y + spin1z * spin1z );
  const REAL4 norma2 = sqrt( spin2x * spin2x + spin2y * spin2y + spin2z * spin2z );
  const REAL4 norma12 = norma1 * norma1;
  const REAL4 norma22 = norma2 * norma2;

  /* FIXME The following calculation for Eq. 7 of arxiv:0904.2577 */
  /* assumes that we have spin (anti-)aligned PhenomB waveforms */
  /* aligned in the \hat{z} direction. Furthermore, we assume that */
  /* the orbital angular momentum L is along the positive \hat{z} */
  /* direction. For the completely general case, we will need to store */
  /* the components of \hat{L}. */
  const int a1sign = SIGN( spin1z );
  const int a2sign = SIGN( spin2z );
  cosalpha = (REAL4) a1sign * a2sign;
  cosbeta = (REAL4) a1sign;
  cosgamma = (REAL4) a2sign;

  norml = 2.0 * sqrt(3.0) + t2 * eta + t3 * eta * eta + s4  / (q2sum * q2sum) *
     (norma12 + norma22 * q4 + 2.0 * norma1 * norma2 * q2 * cosalpha) +
     (s5 * eta + t0 + 2.0) / q2sum * (norma1 * cosbeta + norma2 * q2 * cosgamma);

  return 1.0 / ((1.0 + q) * (1.0 + q)) * sqrt(norma12 + norma22 * q4 + 2.0 * norma1 *
     norma2 * q2 * cosalpha + 2.0 * (norma1 * cosbeta + norma2 * q2 * cosgamma) * norml *
     q + norml * norml * q2);
}


REAL4 XLALBlackHoleRingAmplitude( REAL4 f, REAL4 Q, REAL4 r, REAL4 epsilon )

{
  const REAL4 c = LAL_C_SI;
  const REAL4 M = XLALBlackHoleRingMass( f, Q );
  const REAL4 a = XLALBlackHoleRingSpin( Q );
  const REAL4 g = ring_spin_factor( a );
  const REAL4 F = ring_quality_fn( Q );

  return sqrt(5.0/2.0 * epsilon) *
    ( (LAL_G_SI * M * LAL_MSUN_SI) / ( c * c * r * 1.0e6 * LAL_PC_SI) ) *
    (1.0 / sqrt( Q * F * g) );
}


REAL4 XLALBlackHoleRingEpsilon( REAL4 f, REAL4 Q, REAL4 r, REAL4 amplitude )

{
  const REAL4 c = LAL_C_SI;
  const REAL4 M = XLALBlackHoleRingMass( f, Q );
  const REAL4 a = XLALBlackHoleRingSpin( Q );
  const REAL4 g = ring_spin_factor( a );
  const REAL4 F = ring_quality_fn( Q );

  return (2.0/5.0 * amplitude * amplitude) *
    pow( ( c * c * r * 1.0e6 * LAL_PC_SI) / (LAL_G_SI * M * LAL_MSUN_SI), 2 ) *
    ( Q * F * g );
}


REAL4 XLALBlackHoleRingHRSS( REAL4 f, REAL4 Q, REAL4 amplitude, REAL4 plus, REAL4 cross )

{
  return amplitude * sqrt( Q * ( (0.5 + Q*Q ) * plus*plus + Q * plus*cross + Q*Q * cross*cross ) / ( 1.0 + 4.0 * Q*Q) / LAL_PI / f );
}


REAL8 XLAL2DRingMetricDistance( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb )

{
  REAL8 Q2 = Qa*Qa;
  REAL8 gQQ;
  REAL8 gff;
  REAL8 gQf;

  gQQ = ( 3.0 + 16.0 * Q2 * Q2) / ( Q2 * ( 1.0 + 4.0 * Q2 ) * ( 1.0 + 4.0 * Q2 ) );
  gff = ( 3.0 + 8.0 * Q2) / ( fa * fa);
  gQf = - 2.0 * ( 3.0 + 4.0 * Q2 ) / ( Qa * fa * ( 1.0 + 4.0 * Q2 ));

  return ( 1.0/8.0 * ( gQQ * pow(Qb-Qa,2) + gQf * (Qb-Qa) * (fb-fa) + gff * pow(fb-fa,2) ) );
}


REAL8 XLAL3DRingMetricDistance( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb, REAL8 dt )

{
  REAL8 gQQ, gff, gtt;
  REAL8 gQf, gtf, gtQ;
  REAL8 df, dQ, ds2;
  REAL8 f = (fa+fb)/2.;
  REAL8 Q = (Qa+Qb)/2.;
  REAL8 Q2 = Q*Q;

  gQQ = ( 1. + 28.*Q2*Q2 + 128.*Q2*Q2*Q2 + 64.*Q2*Q2*Q2*Q2) / ( 4. * Q2 * ( 1. + 6.*Q2 + 8.*Q2*Q2 )*( 1. + 6.*Q2 + 8.*Q2*Q2 ) );
  gff = ( 1. + 6.*Q2 + 16.*Q2*Q2) / ( 4. * f*f * ( 1. + 2.*Q2 ) );
  gtt = ( LAL_PI*LAL_PI * f*f ) * ( 1. + 4.*Q2 ) / ( Q2 );
  gQf = - ( 1. + 2.*Q2 + 8.*Q2*Q2 ) / ( 4.*Q*f * ( 1. + 6.*Q2 + 8.*Q2*Q2 ) );
  gtf = - ( LAL_PI * Q ) * ( 1. + 4.*Q2) / ( 1. + 2.*Q2 );
  gtQ = ( LAL_PI * f ) * ( 1. - 2.*Q2 ) / ( ( 1. + 2.*Q2 )*( 1. + 2.*Q2 ) );

  df = fa - fb;
  dQ = Qa - Qb;

  ds2 = ( gQQ * dQ*dQ + gff * df*df + gtt * dt*dt + gQf * 2.*dQ*df + gtf * 2.*dt*df + gtQ * 2.*dt*dQ );

  return ( ds2 );
}



REAL8 XLAL3DRingTimeMinimum( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb)

{
  REAL8 gtt;
  REAL8 gtf, gtQ;
  REAL8 df, dQ, dt;
  REAL8 f = (fa+fb)/2.;
  REAL8 Q = (Qa+Qb)/2.;
  REAL8 Q2 = Q*Q;

  gtt = ( LAL_PI*LAL_PI * f*f ) * ( 1. + 4.*Q2 ) / ( Q2 );
  gtf = - ( LAL_PI * Q ) * ( 1. + 4.*Q2) / ( 1. + 2.*Q2 );
  gtQ = ( LAL_PI * f ) * ( 1. - 2.*Q2 ) / ( ( 1. + 2.*Q2 )*( 1. + 2.*Q2 ) );

  df = fa - fb;
  dQ = Qa - Qb;

  dt = -(gtf * df + gtQ * dQ)/gtt;

  return ( dt );
}


REAL8 XLALRingdownTimeError( const SnglRingdownTable *table,  REAL8 lal_ring_ds_sq )
{
  REAL8 gtt;
  REAL8 f = table->frequency;
  REAL8 Q = table->quality;
  REAL8 Q2 = Q*Q;

  gtt = ( LAL_PI*LAL_PI * f*f ) * ( 1. + 4.*Q2 ) / ( Q2 );

  return ( sqrt( lal_ring_ds_sq / gtt ) );
}

/**
 * This routine computes the ringdown waveform
 * \f{equation}{
 * r(t) = \left\{
 * \begin{array}{ll}
 * e^{-\pi ft/Q}\cos(2\pi ft) & \mbox{for} t\ge0 \\
 * 0 & \mbox{for} t<0
 * \end{array}
 * \right.
 * \f}
 * where the parameters \f$f\f$ and \f$Q\f$ are specified in the input structure.
 * The output must have an appropriate amount of memory allocated, and must
 * have the desired temporal spacing set.  Note: Ref.\ \cite JDECreighton99
 * used a different convention for the ringdown normlization: there the
 * ringdown waveform was taken to be \f$q(t)=(2\pi)^{1/2}r(t)\f$.
 */
int XLALComputeRingTemplate( REAL4TimeSeries *output, SnglRingdownTable *input )

{
  const REAL8 efolds = 10;
  REAL8 amp = 1.0;
  REAL8 fac;
  REAL8 a;
  REAL8 y;
  REAL8 yy;
  UINT4 i;
  UINT4 n;

  if ( ! output || ! output->data || ! input )
    XLAL_ERROR( XLAL_EFAULT );

  if ( ! output->data->length )
    XLAL_ERROR( XLAL_EBADLEN );

  if ( output->deltaT <= 0 || input->quality <= 0 || input->frequency <= 0 )
    XLAL_ERROR( XLAL_EINVAL );

  /* exponential decay variables */
  fac = exp( - LAL_PI * input->frequency * output->deltaT / input->quality );
  n = ceil( - efolds / log( fac ) );

  /* oscillator variables */
  a = 2 * cos( 2 * LAL_PI * input->frequency * output->deltaT );
  y = sin( -2 * LAL_PI * input->frequency * output->deltaT +
      0.5 * LAL_PI + input->phase );
  yy = sin( -4 * LAL_PI * input->frequency * output->deltaT +
      0.5 * LAL_PI + input->phase );

  if ( n < output->data->length )
    memset( output->data->data + n, 0,
        ( output->data->length - n ) * sizeof( *output->data->data ) );
  else
    n = output->data->length;

  for ( i = 0; i < n; ++i )
  {
    REAL4 tmp = a * y - yy;
    yy = y;
    output->data->data[i] = amp * ( y = tmp );
    amp *= fac;
  }

  return 0;
}

/**
 * This routine computes a waveform for a
 * black hole with the specified physical parameters (in the input structure).
 * The parameters are the black hole mass \f$M\f$ (in solar masses \f$M_\odot\f$), the
 * spin \f$S={\hat{a}}GM^2/c\f$ expressed in terms of the dimensionless mass
 * parameter \f${\hat{a}}\f$, the fractional mass lost \f$\epsilon\f$ in ringdown
 * radiation expressed as a percent, and the distance to the typical source
 * (angle-averaged waveform) \f$r\f$ given in megaparsecs (Mpc).  The central
 * frequency and quality of the ringdown are approximated
 * as\ \cite EWLeaver85,\cite EBertiEtAl06:
 * \f{equation}{
 * f \simeq 32\,\textrm{kHz}\times[f_1+f_2(1-{\hat{a}})^{f_3}](M_\odot/M)
 * \f}
 * and
 * \f{equation}{
 * Q \simeq q_1+q_2(1-{\hat{a}})^{q_3},
 * \f}
 * where the values of the constants (f_1,f_2,f_3) and (q_1,q_2,q_3) are
 * given for each of (l,m,n) in\ \cite EBertiEtAl06.
 * The strain waveform produced is \f$h(t)=A_q q(t)\f$ where the amplitude factor
 * is\ \cite JDECreighton99
 * \f{equation}{
 * A_q = 2.415\times10^{-21}Q^{-1/2}[1-0.63(1-{\hat{a}})^{3/10}]^{-1/2}
 * \left(\frac{\textrm{Mpc}}{r}\right)
 * \left(\frac{M}{M_\odot}\right)
 * \left(\frac{\epsilon}{0.01}\right)^{1/2}.
 * \f}
 * Note that this is written \f$A_q\f$ to emphasize that it is the amplitude
 * factor for \f$q(t)\f$ rather than \f$r(t)\f$.
 */
int XLALComputeBlackHoleRing(
    REAL4TimeSeries     *output,
    SnglRingdownTable   *input,
    REAL4                dynRange
    )

{
  const REAL4 amp = dynRange *
    XLALBlackHoleRingAmplitude(
        input->frequency, input->quality, input->eff_dist, input->epsilon );
  UINT4 i;

  if ( XLALComputeRingTemplate( output, input ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );

  for ( i = 0; i < output->data->length; ++i )
    output->data->data[i] *= amp;

  return 0;
}


static int MakeBank( SnglRingdownTable *tmplt, RingTemplateBankInput *input )
{
  UINT4 count = 0;
  REAL4 dseff = 4 * sqrt( input->maxMismatch );
  REAL4 minlogf = log( input->minFrequency );
  REAL4 maxlogf = log( input->maxFrequency );
  REAL4 q = input->minQuality;

  while ( q < input->maxQuality )
  {
    REAL4 q2 = q * q;
    REAL4 logfreq = minlogf;

    while ( logfreq < maxlogf )
    {
      if ( tmplt )
      {
        tmplt[count].quality   = q;
        tmplt[count].frequency = exp( logfreq );
        tmplt[count].phase     = input->templatePhase;
        tmplt[count].epsilon   = input->templateEpsilon;
        tmplt[count].eff_dist  = input->templateDistance;
      }
      ++count;
      logfreq += dseff / sqrt( 3 + 8 * q2 );
    }

    q += dseff * q * ( 1 + 4 * q2 ) / sqrt( 3 + 16 * q2 * q2 );
  }

  return count;
}

/**
 * This routine creates a bank of ringdown
 * templates that cover a set range in the parameters \f$f\f$ and \f$Q\f$.
 */
RingTemplateBank *XLALCreateRingTemplateBank( RingTemplateBankInput *input )

{
  UINT4 i;
  RingTemplateBank *bank;

  if ( ! input )
    XLAL_ERROR_NULL( XLAL_EFAULT );

  bank = LALCalloc( 1, sizeof( *bank ) );
  if ( ! bank )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

  bank->numTmplt = MakeBank( NULL, input );
  bank->tmplt = LALCalloc( bank->numTmplt, sizeof( *bank->tmplt ) );
  if ( ! bank->tmplt )
  {
    LALFree( bank );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }
  for ( i = 0; i < bank->numTmplt - 1; ++i )
    bank->tmplt[i].next = bank->tmplt + i + 1;

  MakeBank( bank->tmplt, input );
  return bank;
}


/** Destroys a RingTemplateBank */
void XLALDestroyRingTemplateBank( RingTemplateBank *bank )

{
  if ( bank )
  {
    if ( bank->tmplt )
      LALFree( bank->tmplt );
    LALFree( bank );
  }
  return;
}
