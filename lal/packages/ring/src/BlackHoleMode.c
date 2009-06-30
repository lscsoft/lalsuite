/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>
#include <lal/BlackHoleMode.h>

NRCSID (BLACKHOLEMODEC,"$Id$");


static int mysignbit( REAL8 x )
{
  if ( x == 0 )
  {
    REAL8 y = 0;
    if ( memcmp( &x, &y, sizeof(x) ) )
      return 1;
    return 0;
  }
  return x < 0.0;
}


/*
 *
 * NOTE: STRANGE CONVENTIONS ...
 *
 * Use sign of a to indicate branch: a > 0 means +ve omega branch.
 *                                   a < 0 means -ve omega branch.
 *
 */


/* preprocessor macros */
#define TINY LAL_REAL8_MIN
#define EPS LAL_REAL4_EPS
#define SMALL LAL_REAL4_EPS
#define MAXITER 10000
#define ctiny crect(TINY,0)
#define RCmul(x,a) cmulr((a),(x))

enum { LAL_BLACK_HOLE_MODE_DATA_TABLE_SIZE = 10 };

/*
 *
 * XLAL Routines
 *
 */



static int XLALBlackHoleModeEigenSolveContinuedFraction(
    COMPLEX16 *result,
    int (*coef)(COMPLEX16 *, COMPLEX16 *, COMPLEX16 *, int n, void *),
    void *params
    )
{
  static const char *func = "XLALBlackHoleModeEigenSolveContinuedFraction";
  int n = 0;
  COMPLEX16 a, b, alp, gam;
  COMPLEX16 C, D, Delta;

  a = czero;
  coef( &alp, &b, &gam, n, params );
  a = cneg(cmul(a,gam));
  if ( cabs(b) < TINY )
    *result = ctiny;
  else
    *result = b;
  C = *result;
  D = czero;

  while ( n++ < MAXITER )
  {
    COMPLEX16 prev = *result;
    a = alp;
    coef( &alp, &b, &gam, n, params );
    a = cneg(cmul(a,gam));
    C = cadd(b,cdiv(a,C));
    D = cadd(b,cmul(a,D));
    if (cabs(D) < TINY)
      D = ctiny;
    if (cabs(C) < TINY)
      C = ctiny;
    D = cinv(D);
    Delta = cmul(C,D);
    *result = cmul(*result,Delta);
    if (cabs(csub(*result,prev))<SMALL)
      return 0;
  }
  XLAL_ERROR( func, XLAL_EMAXITER );
  return -1;
}



static int XLALBlackHoleModeSchwarzschildCoefficients( COMPLEX16 *alp, COMPLEX16 *bet, COMPLEX16 *gam, int n, void *params )
{
  struct tagBlackHoleMode *p = params;
  int l, s;
  COMPLEX16 omega, rho, A;

  s = p->s;
  l = p->l;
  omega = p->omega;
  A = crect(l*(l+1)-s*(s+1),0.0);
  rho = cneg(cmul(I,omega));
  *alp = cadd(RCmul(n,cadd(csetr(n),RCmul(2,cadd(cunit,rho)))),
	      cadd(RCmul(2,rho),cunit));
  *bet = cneg(cadd(RCmul(2*n,cadd(csetr(n),cadd(RCmul(4,rho),cunit))),
    cadd(cmul(RCmul(4,rho),cadd(RCmul(2,rho),cunit)),cadd(A,csetr(1+s)))));
  *gam = cadd(RCmul(n,cadd(csetr(n),RCmul(4,rho))),
              cadd(RCmul(4,cmul(rho,rho)),csetr(-s*s)));

  return 0;
}

static int XLALBlackHoleModeKerrRadialCoefficients( COMPLEX16 *alp, COMPLEX16 *bet, COMPLEX16 *gam, int n, void *params )
{
  struct tagBlackHoleMode *p = params;
  int m=p->m, s=p->s;
  REAL8 a=p->a;
  REAL8 b=sqrt(1-4*a*a);
  COMPLEX16 omega=p->omega, A=p->A;
  COMPLEX16 c0,c1,c2,c3,c4,cfact,rho,omega2;

  rho = cneg(cmul(I,omega));
  omega2 = cmul(omega,omega);

  cfact = RCmul(1/b,csub(omega,csetr(2*a*m)));
  c0 = cadd(csetr(1-s),cadd(rho,cneg(cmul(I,cfact))));
  c1 = cadd(csetr(-4),cadd(RCmul(-2*(2+b),rho),cmul(cseti(2),cfact)));
  c2 = cadd(csetr(s+3),cadd(RCmul(3,rho),cneg(cmul(I,cfact))));
  c3 = cadd(RCmul(4+2*b-a*a,omega2),
            cadd(RCmul(-2*m*a,omega),
                 cadd(csetr(-s-1),
                      cadd(RCmul(-2-b,rho),
                           cadd(cneg(A),
                                cmul(cadd(RCmul(2,omega),I),cfact))))));
  c4 = cadd(csetr(s+1),
            cadd(RCmul(-2,omega2),
                 cadd(RCmul(2*s+3,rho),
                      cneg(cmul(cadd(RCmul(2,omega),I),cfact)))));
  *alp = cadd(RCmul(n,cadd(csetr(n),cadd(c0,cunit))),c0);
  *bet = cadd(RCmul(n,cadd(csetr(-2*n),cadd(c1,csetr(2)))),c3);
  *gam = cadd(RCmul(n,cadd(csetr(n),cadd(c2,csetr(-3)))),
              cadd(c4,cadd(cneg(c2),csetr(2))));
  return 0;
}

static int XLALBlackHoleModeKerrAngularCoefficients( COMPLEX16 *alp, COMPLEX16 *bet, COMPLEX16 *gam, int n, void *params )
{
  struct tagBlackHoleMode *p = params;
  int m=p->m, s=p->s;
  REAL8 a=p->a;
  COMPLEX16 omega=p->omega, A=p->A;
  COMPLEX16 w;
  REAL8 k1=0.5*abs(m-s),k2=0.5*abs(m+s);

  w = RCmul(a,omega);
  *alp = csetr(-2*(n+1)*(n+2*k1+1));
  *bet = cadd(csetr(n*(n+2*k1+2*k2+1)+(k1+k2)*(k1+k2+1)-s*(s+1)),
	      cadd(RCmul(-2*(2*n+2*k1+s+1),w),
                   cadd(cneg(cmul(w,w)),cneg(A))));
  *gam = RCmul(2*(n+k1+k2+s),w);
  return 0;
}

static int XLALBlackHoleModeEigenSolveSchwarzschildResid( const gsl_vector *x, void *params, gsl_vector *f )
{
  COMPLEX16 cf;
  struct tagBlackHoleMode *p = params;
  creal(p->omega) = gsl_vector_get(x,0);
  cimag(p->omega) = gsl_vector_get(x,1);

  XLALBlackHoleModeEigenSolveContinuedFraction( &cf, XLALBlackHoleModeSchwarzschildCoefficients, p );

  gsl_vector_set(f, 0, creal(cf));
  gsl_vector_set(f, 1, cimag(cf));
  return 0;
}

static int XLALBlackHoleModeEigenSolveKerrResid( const gsl_vector *x, void *params, gsl_vector *f )
{
  COMPLEX16 cf1, cf2;
  struct tagBlackHoleMode *p = params;
  creal(p->A) = gsl_vector_get(x,0);
  cimag(p->A) = gsl_vector_get(x,1);
  creal(p->omega) = gsl_vector_get(x,2);
  cimag(p->omega) = gsl_vector_get(x,3);

  XLALBlackHoleModeEigenSolveContinuedFraction( &cf1, XLALBlackHoleModeKerrRadialCoefficients, p );
  XLALBlackHoleModeEigenSolveContinuedFraction( &cf2, XLALBlackHoleModeKerrAngularCoefficients, p );

  gsl_vector_set(f, 0, creal(cf1));
  gsl_vector_set(f, 1, cimag(cf1));
  gsl_vector_set(f, 2, creal(cf2));
  gsl_vector_set(f, 3, cimag(cf2));
  return 0;
}


static int XLALBlackHoleModeEigenSolveSchwarzschild( COMPLEX16 *omega, int l, int m, int s )
{
  static const char *func = "XLALBlackHoleModeEigenSolveSchwarzschild";
  enum { ndim = 2 };
  struct tagBlackHoleMode p;
  const gsl_multiroot_fsolver_type *solverType;
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function f;
  gsl_vector *x;
  size_t iter = 0;
  int status;

  f.f = &XLALBlackHoleModeEigenSolveSchwarzschildResid;
  f.n = ndim;
  f.params = &p;
  x = gsl_vector_alloc( ndim );
  gsl_vector_set( x, 0, omega->re );
  gsl_vector_set( x, 1, omega->im );
  p.a = 0.0;
  p.l = l;
  p.m = m;
  p.s = s;

  solverType = gsl_multiroot_fsolver_hybrids;
  solver = gsl_multiroot_fsolver_alloc( solverType, ndim );
  gsl_multiroot_fsolver_set(solver, &f, x);

  do
  {
    ++iter;
    status = gsl_multiroot_fsolver_iterate( solver );
    if ( status )
      break;
    status = gsl_multiroot_test_residual(solver->f, EPS);
  }
  while ( status == GSL_CONTINUE && iter < MAXITER );
  if ( iter >= MAXITER )
    XLAL_ERROR( func, XLAL_EMAXITER );

  omega->re = gsl_vector_get(solver->x, 0);
  omega->im = gsl_vector_get(solver->x, 1);

  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);
  return 0;
}

static int XLALBlackHoleModeEigenSolveKerr( COMPLEX16 *A, COMPLEX16 *omega, REAL8 a, int l, int m, int s )
{
  static const char *func = "XLALBlackHoleModeEigenSolveKerr";
  enum { ndim = 4 };
  struct tagBlackHoleMode p;
  const gsl_multiroot_fsolver_type *solverType;
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function f;
  gsl_vector *x;
  size_t iter = 0;
  int status;

  f.f = &XLALBlackHoleModeEigenSolveKerrResid;
  f.n = ndim;
  f.params = &p;
  x = gsl_vector_alloc( ndim );
  gsl_vector_set( x, 0, A->re );
  gsl_vector_set( x, 1, A->im );
  gsl_vector_set( x, 2, omega->re );
  gsl_vector_set( x, 3, omega->im );
  p.a = a;
  p.l = l;
  p.m = m;
  p.s = s;

  solverType = gsl_multiroot_fsolver_hybrids;
  solver = gsl_multiroot_fsolver_alloc( solverType, ndim );
  gsl_multiroot_fsolver_set(solver, &f, x);

  do
  {
    ++iter;
    status = gsl_multiroot_fsolver_iterate( solver );
    if ( status )
      break;
    status = gsl_multiroot_test_residual(solver->f, EPS);
  }
  while ( status == GSL_CONTINUE && iter < MAXITER );
  if ( iter >= MAXITER )
    XLAL_ERROR( func, XLAL_EMAXITER );

  A->re = gsl_vector_get(solver->x, 0);
  A->im = gsl_vector_get(solver->x, 1);
  omega->re = gsl_vector_get(solver->x, 2);
  omega->im = gsl_vector_get(solver->x, 3);

  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);
  return 0;
}

/* Note: A and omega are input and output.
 * On input, they are an initial guess for the solver */
static int XLALBlackHoleModeEigenSolve( COMPLEX16 *A, COMPLEX16 *omega, REAL8 a, int l, int m, int s )
{
  static const char *func = "XLALBlackHoleModeEigenSolve";
  int aIsNegative;
  int status;

  /* if a is negative, use an identity to obtain eigenvalues from -m mode */
  if ( ( aIsNegative = mysignbit( a ) ) )
  {
    a = fabs(a);
    m = -m;
    /* now use identity to modify the guessed eigenvalue */
    A->im = -A->im;
    omega->re = -omega->re;
  }

  if ( a >= 0.5 || l < 0 || abs(m) > l || s > 0 || s < -2 )
    /* EINVAL */
    return -1;

  if ( a == 0 )
    status = XLALBlackHoleModeEigenSolveSchwarzschild( omega, l, m, s );
  else
    status = XLALBlackHoleModeEigenSolveKerr( A, omega, a, l, m, s );

  if ( status < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  /* if a was negative, apply the identity to restore this eigenvalue */
  if ( aIsNegative )
  {
    A->im = -A->im;
    omega->re = -omega->re;
  }

  return 0;
}


void XLALDestroyBlackHoleModeTable( BlackHoleModeTable *mode )
{
  if ( mode )
  {
    LALFree( mode->positiveEigenTable );
    LALFree( mode->negativeEigenTable );
    LALFree( mode );
  }
  return;
}

static int XLALBlackHoleModeMakeEigenTable( BlackHoleModeTable *mode )
{
  static const char *func = "XLALBlackHoleModeMakeEigenTable";
  const REAL8 fact = 0.19245008972987525483; /* 1.0/sqrt(27.0) */
  size_t nobj = LAL_BLACK_HOLE_MODE_DATA_TABLE_SIZE;
  int l,m,s;
  REAL8 a;
  COMPLEX16 A;
  COMPLEX16 omega;
  size_t i;

  if ( mode->positiveEigenTable || mode->negativeEigenTable )
    XLAL_ERROR( func, XLAL_EINVAL );

  mode->eigenTableSize = nobj;
  mode->positiveEigenTable = LALCalloc(nobj, sizeof(*mode->positiveEigenTable));
  mode->negativeEigenTable = LALCalloc(nobj, sizeof(*mode->negativeEigenTable));
  if ( ! mode->positiveEigenTable || ! mode->negativeEigenTable )
    XLAL_ERROR( func, XLAL_ENOMEM );

  l = mode->l;
  m = mode->m;
  s = mode->s;

  /* do the positive eigenvalues first */
  A = crect(l*(l+1) - s*(s+1),0);
  omega = RCmul(fact,crect(2*l+1,-1));
  a = 0;
  for ( i = 0; i < nobj; ++i )
  {
    struct tagBlackHoleModeEigenvalues *eigen = mode->positiveEigenTable + i;
    a = 0.5 * tanh(0.5*i);
    if ( XLALBlackHoleModeEigenSolve( &A, &omega, a, l, m, s ) < 0 )
    {
      XLALDestroyBlackHoleModeTable( mode );
      XLAL_ERROR( func, XLAL_EFUNC );
    }
    eigen->a = a;
    eigen->A = A;
    eigen->omega = omega;
    a += 0.1 * (0.5 - a);
  }

  /* now do the negative eigenvalues */
  A = crect(l*(l+1) - s*(s+1),0);
  omega = RCmul(fact,crect(-2*l+1,-1));
  a = 0;
  for ( i = 0; i < nobj; ++i )
  {
    struct tagBlackHoleModeEigenvalues *eigen = mode->negativeEigenTable + i;
    a = 0.5 * tanh(0.5*i);
    if ( XLALBlackHoleModeEigenSolve( &A, &omega, -a, l, m, s ) < 0 )
    {
      XLALDestroyBlackHoleModeTable( mode );
      XLAL_ERROR( func, XLAL_EFUNC );
    }
    eigen->a = -a;
    eigen->A = A;
    eigen->omega = omega;
    a += 0.1 * (0.5 - a);
  }

  return 0;
}


BlackHoleModeTable * XLALCreateBlackHoleModeTable( int l, int m, int s )
{
  static const char *func = "XLALCreateBlackHoleModeTable";
  BlackHoleModeTable *mode;

  if ( l < 0 || abs(m) > l || s > 0 || s < -2 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );

  mode = LALCalloc( 1, sizeof( *mode ) );
  if ( ! mode )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  mode->l = l;
  mode->m = m;
  mode->s = s;

  if ( XLALBlackHoleModeMakeEigenTable( mode ) < 0 )
  {
    XLALDestroyBlackHoleModeTable( mode );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }

  return mode;
}


int XLALBlackHoleModeEigenvaluesLeaverT( COMPLEX16 *A, COMPLEX16 *omega, REAL8 a, BlackHoleModeTable *mode )
{
  static const char *func = "XLALBlackHoleModeEigenvaluesLeaverT";
  struct tagBlackHoleModeEigenvalues *eigen;
  int l, m, s;
  size_t i;

  if ( ! (fabs(a) < 0.5) )
    XLAL_ERROR( func, XLAL_EINVAL );

  l = mode->l;
  m = mode->m;
  s = mode->s;

  /* locate the value of a in the table that is closer to zero than
   * the requested value */
  eigen = mysignbit(a) ? mode->negativeEigenTable : mode->positiveEigenTable;
  for ( i = 1; i < mode->eigenTableSize; ++i )
    if ( fabs( a ) < fabs( eigen[i].a ) )
      break;

  /* here is the guess */
  *A = eigen[i-1].A;
  *omega = eigen[i-1].omega;

  if ( XLALBlackHoleModeEigenSolve( A, omega, a, l, m, s ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  return 0;
}


int XLALBlackHoleModeEigenvaluesLeaver( COMPLEX16 *A, COMPLEX16 *omega, REAL8 a, int l, int m, int s )
{
  static const char *func = "XLALBlackHoleModeEigenvaluesLeaver";
  const REAL8 fact = 0.19245008972987525483; /* 1.0/sqrt(27.0) */
  size_t nobj = LAL_BLACK_HOLE_MODE_DATA_TABLE_SIZE;
  size_t i;
  int sign;

  if ( ! (fabs(a) < 0.5) )
    XLAL_ERROR( func, XLAL_EINVAL );
  sign = mysignbit(a) ? -1 : 1;

  /* step towards desired value of a */
  *A = crect(l*(l+1) - s*(s+1),0);
  *omega = RCmul(fact,crect(sign*2*l+1,-1));
  for ( i = 0; i < nobj; ++i )
  {
    REAL8 atry = sign * 0.5 * tanh(0.5*i);
    if ( fabs(atry) > fabs(a) )
      break;
    if ( XLALBlackHoleModeEigenSolve( A, omega, atry, l, m, s ) < 0 )
      XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* we now have our guess */

  if ( XLALBlackHoleModeEigenSolve( A, omega, a, l, m, s ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  return 0;
}


int XLALSpheroidalWaveFunction1( COMPLEX16 *result, REAL8 mu, struct tagBlackHoleMode *mode )
{
#if 0
  int s = mode->s, m = mode->m;
  REAL8 k1 = 0.5*abs(m - s), k2 = 0.5*abs(m + s);
  REAL8 fact = 1 + mu;
  REAL8 prod = fact;
  COMPLEX16 alp, bet, gam;
  COMPLEX16 sum, a, a_;
  size_t n = 0;

  sum = a_ = cunit;
  XLALBlackHoleModeKerrAngularCoefficients(&alp, &bet, &gam, n, mode);
  a = cdiv(cmul(bet,a_),cneg(alp));
  while ( n++ < MAXITER )
  {
    COMPLEX16 aa, delta;
    delta = cmulr(a,prod);
    sum = cadd(sum,delta);
    if ( cabs2(delta) < SMALL )
      break;
    XLALBlackHoleModeKerrAngularCoefficients(&alp, &bet, &gam, n, mode);
    aa = a;
    a = cdiv(cadd(cmul(bet,a),cmul(gam,a_)),cneg(alp));
    a_ = aa;
    prod *= fact;
  }
  if ( n >= MAXITER)
    /* EMAXIT */
    abort();

  *result = cmulr(sum,pow(1+mu,k1)*pow(1-mu,k2));
  *result = cmul(*result,cexp(cmulr(mode->omega,mu*mode->a)));
  return 0;

#else
  static const char *func = "XLALSpheroidalWaveFunction1";
  struct tagBlackHoleMode params;
  int s=mode->s, m=mode->m;
  REAL8 k1=0.5*abs(m-s), k2=0.5*abs(m+s);
  REAL8 prod, fact=1+mu;
  COMPLEX16 sum, term;
  COMPLEX16 aprev, a, anext;
  COMPLEX16 alp, bet, gam;
  size_t n = 0;

  params = *mode;
  /* params.a = fabs(params.a); */
  if ( XLALBlackHoleModeKerrAngularCoefficients(&alp, &bet, &gam, n, &params) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  prod  = 1;
  sum   = cunit;
  a     = cunit;
  anext = cneg(cdiv(cmul(bet,a),alp)); /* Eq. 19, first line */

  while ( n++ < MAXITER )
  {
    aprev = a;
    a     = anext;
    prod *= fact;
    term  = cmulr(a,prod);
    sum   = cadd(sum,term);
    /*
    if ( cabs2(term) < SMALL )
      break;
      */
    if ( cabs2(term) < 1e-4 )
      break;
    if ( XLALBlackHoleModeKerrAngularCoefficients(&alp, &bet, &gam, n, &params) < 0 )
      XLAL_ERROR( func, XLAL_EFUNC );
    /* Eq. 19, second line: */
    anext = cneg(cdiv(cadd(cmul(bet,a),cmul(gam,aprev)),alp));
  }
  if ( n >= MAXITER)
    XLAL_ERROR( func, XLAL_EMAXITER );

  *result = cmulr(sum,pow(1+mu,k1)*pow(1-mu,k2));
  *result = cmul(*result,cexp(cmulr(mode->omega,mu*mode->a)));
  return 0;
#endif
}


static REAL8 XLALSpheroidalWaveFunction1Abs2( REAL8 mu, void *params )
{
  static const char *func = "XLALSpheroidalWaveFunction1Abs2";
  COMPLEX16 swf;
  REAL8 result;
  if ( XLALSpheroidalWaveFunction1( &swf, mu, params ) )
    XLAL_ERROR_REAL8( func, XLAL_EFUNC );
  result = cabs2(swf);
  return result;
}

int XLALSpheroidalWaveFunctionNorm( COMPLEX16 *norm, struct tagBlackHoleMode *params )
{
  static const char *func = "XLALSpheroidalWaveFunctionNorm";
  enum { WORKSPACESIZE = 1000 };
  gsl_integration_workspace *work;
  gsl_function function;
  REAL8 integral, error;
  COMPLEX16 swf;
  int signneg;

  work = gsl_integration_workspace_alloc( WORKSPACESIZE );
  function.function = XLALSpheroidalWaveFunction1Abs2;
  function.params   = params;
  /* gsl_integration_qags( &function, -1, 1, 0, LAL_REAL8_EPS, WORKSPACESIZE, work, &integral, &error ); */
  gsl_integration_qags( &function, -1, 1, 0, 1e-6, WORKSPACESIZE, work, &integral, &error );
  gsl_integration_workspace_free( work );

  /* get complex part so that spheroidal wave function is real at mu=0 */
  if ( XLALSpheroidalWaveFunction1( &swf, 0.0, params ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );
  *norm = conj(cdivr(swf,cabs(swf)));

  /* sign convention: to agree with sw spherical harmonics */
  /* TODO: CHECKME */

  if ( XLALSpheroidalWaveFunction1( &swf, -1.0, params ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );
  swf = cmul(swf,*norm);
  signneg = ( params->l - ( params->m > params->s ? params->m : params->s ) ) % 2 ? 0 : 1;
  if ( signneg && creal(swf) > 0 || ! signneg && creal(swf) < 0 )
    *norm = cneg(*norm);

  *norm = cdivr(*norm,sqrt(integral));

  return 0;
}

int XLALSpheroidalWaveFunction( COMPLEX16 *result, REAL8 mu, struct tagBlackHoleMode *params )
{
  static const char *func = "XLALSpheroidalWaveFunction";
  COMPLEX16 norm;
  if ( XLALSpheroidalWaveFunctionNorm( &norm, params ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );
  if ( XLALSpheroidalWaveFunction1( result, mu, params ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );
  *result = cmul( *result, norm );
  return 0;
}


/* Here a is in standard conventions, not Leaver's */
int XLALSetBlackHoleModeParams( struct tagBlackHoleMode *params, REAL8 a, int l, int m, int s )
{
  static const char *func = "XLALSetBlackHoleModeParams";
  if ( l < 0 || abs(m) > l || s > 0 || s < -2 )
    XLAL_ERROR( func, XLAL_EINVAL );
  params->a = 0.5 * a; /* to convert from standard conventions to Leaver conventions */
  params->l = l;
  params->m = m;
  params->s = s;
  XLALBlackHoleModeEigenvaluesLeaver( &params->A, &params->omega, params->a, params->l, params->m, params->s );
  return 0;
}

/* Here a is in standarad conventions, not Leaver's */
int XLALSpinWeightedSpheroidalHarmonic( COMPLEX16 *result, REAL8 mu, REAL8 a, int l, int m, int s )
{
  struct tagBlackHoleMode params;
  XLALSetBlackHoleModeParams( &params, a, l, m, s );
  XLALSpheroidalWaveFunction( result, mu, &params );
  return 0;
}

int XLALBlackHoleRingdownAmplitude(
    COMPLEX16 *amplitudePlus,
    COMPLEX16 *amplitudeCross,
    REAL8 massSolar,
    REAL8 fracMassLoss,
    REAL8 distanceMpc,
    REAL8 inclinationRad,
    REAL8 azimuthRad,
    BlackHoleMode *mode
    )
{
  COMPLEX16 swf_1, swf_2;
  COMPLEX16 omega, amp;
  REAL8 mu; /* cosine of inclinaion */

  mu = cos( inclinationRad );
  XLALSpheroidalWaveFunction( &swf_1, mu, mode );
  XLALSpheroidalWaveFunction( &swf_2, -mu, mode );

  /* change from Leaver's conventions to usual conventions */
  omega = cmulr( mode->omega, 0.5 );

  amp = cexp(cmulr(I, mode->m*azimuthRad));
  amp = cmulr( amp, -4.0*massSolar*sqrt(-1.0*cimag(omega)*0.5*fracMassLoss/cabs2(omega)));
  amp = cmulr( amp, LAL_MRSUN_SI/(distanceMpc*1e6*LAL_PC_SI));

  *amplitudePlus = cmul(cadd(swf_1,swf_2),amp);
  *amplitudeCross = cmul(cmul(I,csub(swf_1,swf_2)),amp);

  return 0;
}

int XLALBlackHoleRingdownWaveform(
    REAL4TimeSeries *plus,
    REAL4TimeSeries *cross,
    REAL8 massSolar,
    REAL8 spinParam,
    REAL8 fracMassLoss,
    REAL8 distanceMpc,
    REAL8 inclinationRad,
    REAL8 azimuthRad,
    int l,
    int m
    )
{
  static const char *func = "XLALBlackHoleRingdownWaveform";
  const int s = -2;
  struct tagBlackHoleMode params;
  COMPLEX16 swf_1, swf_2, omega;
  REAL8 mu; /* cosine of inclinaion */
  REAL8 dt; /* time step in seconds */
  COMPLEX16 amplitude_1, I_omega_dt_1;
  COMPLEX16 amplitude_2, I_omega_dt_2;
  INT4 ndat;
  INT4 i;

  if ( ! plus || ! cross )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( plus->deltaT != cross->deltaT )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! plus->data || ! cross->data
      || plus->data->length != cross->data->length )
    XLAL_ERROR( func, XLAL_EINVAL );

  mu = cos( inclinationRad );
  dt = plus->deltaT;

  XLALSetBlackHoleModeParams( &params, spinParam, l, m, s );
  XLALSpheroidalWaveFunction( &swf_1, mu, &params );
  XLALSpheroidalWaveFunction( &swf_2, -mu, &params );

  /* change from Leaver's conventions to usual conventions */
  omega = cmulr( params.omega, 0.5 );

  /* how many points in the waveform to compute */
  ndat = ceil( log( LAL_REAL4_EPS ) * massSolar * LAL_MTSUN_SI / ( cimag(omega) * dt ) );
  if ( ndat < 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( ndat > (INT4)plus->data->length )
  {
    /* warning: waveform will be truncated */
    ndat = plus->data->length;
  }

  /* compute complex amplitude and phase factors */
  amplitude_1 = cmul(swf_1, cexp(cmulr(I, m*azimuthRad)));
  amplitude_1 = cmulr(amplitude_1, -4.0*massSolar*sqrt(-1.0*cimag(omega)*0.5*fracMassLoss/cabs2(omega)));
  amplitude_1 = cmulr(amplitude_1, LAL_MRSUN_SI/(distanceMpc*1e6*LAL_PC_SI));
  I_omega_dt_1 = cmulr(cmul(I,omega),-1.0*dt/(massSolar*LAL_MTSUN_SI));

  amplitude_2 = conj(cmul(swf_2, cexp(cmulr(I, m*azimuthRad))));
  amplitude_2 = cmulr(amplitude_2, -4.0*massSolar*sqrt(-1.0*cimag(omega)*0.5*fracMassLoss/cabs2(omega)));
  amplitude_2 = cmulr(amplitude_2, LAL_MRSUN_SI/(distanceMpc*1e6*LAL_PC_SI));
  I_omega_dt_2 = cmulr(cmul(I,conj(omega)),dt/(massSolar*LAL_MTSUN_SI));

  /* zero the data */
  memset( plus->data->data, 0, plus->data->length * sizeof( *plus->data->data ) );
  memset( cross->data->data, 0, cross->data->length * sizeof( *cross->data->data ) );

  for ( i = 0; i < ndat; ++i )
  {
    COMPLEX16 strain_1, strain_2, strain;
    strain_1 = cmul( amplitude_1, cexp( cmulr( I_omega_dt_1, i ) ) );
    strain_2 = cmul( amplitude_2, cexp( cmulr( I_omega_dt_2, i ) ) );
    strain = cadd( strain_1, strain_2 );
    plus->data->data[i] = creal( strain );
    cross->data->data[i] = -1.0 * cimag( strain );
  }

  return ndat;
}
