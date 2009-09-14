/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Julien Sylvestre
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

#include <lal/LPC.h>
#include "lal/LALRCSID.h"
#include <lal/LALError.h>
#include <lal/LALStatusMacros.h>
#include <lal/RealFFT.h>
#include <string.h>

NRCSID (LPCC, "$Id$");

static INT4 balanc(REAL4 **a, INT4 n);
static INT4 toeplz(REAL4 r[], REAL4 x[], REAL4 y[], INT4 n);
static INT4 zrhqr(REAL4 a[], INT4 m, REAL4 rtr[], REAL4 rti[]);
static INT4 hqr(REAL4 **a, INT4 n, REAL4 wr[], REAL4 wi[]);

/******** <lalLaTeX file="LPCC"> ********
\noindent
Compute the coefficients of a linear predictor filter.
\subsubsection*{Prototype}
********* </lalLaTeX> ********/
/* <lalVerbatim> */
void LALLPC(LALStatus *status,
	    REAL4Vector *aout,    /* returned filter coefficients */
	    REAL4Vector *x,    /* training data */
	    UINT4 p            /* filter order */
	    ) {
/* </lalVerbatim> */
/******** <lalLaTeX file="LPCC"> ********
\subsubsection*{Description}
Train a FIR filter aout of order p on the data x.

\subsubsection*{Uses}
\begin{verbatim}
...
\end{verbatim}
********* </lalLaTeX> ********/
  UINT4 i,/*j,*/npad;
  REAL4 *r, *R, *a, *y;
  RealFFTPlan *pfwd = NULL, *pinv = NULL;
  COMPLEX8Vector *Hvec = NULL;
  REAL4Vector rv;

  INITSTATUS (status, "LALLPC", LPCC);
  ATTATCHSTATUSPTR (status);

  ASSERT(x->length>=p+1, status, LPCH_EIN, LPCH_MSGEIN);

  /* compute correlations */
  npad = (UINT4)pow(2.0,ceil(log((REAL4)(x->length))/log(2.0)));

  pfwd = XLALCreateForwardREAL4FFTPlan( npad, 0 );
  if (pfwd == NULL)
    ABORTXLAL(status);

  pinv = XLALCreateReverseREAL4FFTPlan( npad, 0 );
  if (pinv == NULL)
    ABORTXLAL(status);

  LALCCreateVector( status->statusPtr, &Hvec, npad/2 + 1);
  CHECKSTATUSPTR (status);

  r = (REAL4 *)LALMalloc(npad * sizeof(REAL4));
  if(!r) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

  rv.length = npad;
  rv.data = r;

  memcpy(r,x->data,x->length * sizeof(REAL4));
  bzero(r + x->length, (npad - x->length) * sizeof(REAL4));

  if (XLALREAL4ForwardFFT(Hvec, &rv, pfwd) != 0)
    ABORTXLAL(status);

  for(i=0;i<Hvec->length;i++) {
    Hvec->data[i].re = (Hvec->data[i].re*Hvec->data[i].re + Hvec->data[i].im*Hvec->data[i].im)/((REAL4)(npad*(x->length-1)));
    Hvec->data[i].im = 0.0;
  }

  if (XLALREAL4ReverseFFT(&rv, Hvec, pinv) != 0)
    ABORTXLAL(status);

  XLALDestroyREAL4FFTPlan(pinv);

  XLALDestroyREAL4FFTPlan(pfwd);

  LALCDestroyVector( status->statusPtr, &Hvec);
  CHECKSTATUSPTR (status);


  /* pack correlations */
  R = (REAL4 *)LALMalloc((2*p-1) * sizeof(REAL4));
  if(!R) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

  R[p-1] = r[0];
  for(i=0;i<p-1;i++) {
    R[p-i-2] = R[i+p] = r[i+1];
  }

  y = (REAL4 *)LALMalloc(p * sizeof(REAL4));
  if(!y) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

  for(i=0;i<p;i++) {
    y[i] = -r[i+1];
  }

  a = (REAL4 *)LALMalloc((1+p) * sizeof(REAL4));
  if(!a) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

  a[0] = 1.0;

  /* solve toeplitz system */
  if(toeplz(R,a+1,y,p)) {
    ABORT(status, LPCH_ENUM, LPCH_MSGENUM);
  }

  /* set output */
  if(aout->data) {
    memcpy(aout->data, a, (1+p)*sizeof(REAL4));
    aout->length = p+1;
    LALFree(a);
  } else {
    aout->data = a;
    aout->length = p+1;
  }


  /* stabilize polynomial */
  /*
  {
    REAL4 a0;

    LALPolystab(status->statusPtr, aout);
    CHECKSTATUSPTR (status);

    a0 = aout->data[0];

    for(i=0;i<aout->length;i++) {
      aout->data[0] /= a0;
    }
  }
  */

  /* clean up */
  LALFree(y);
  LALFree(R);
  LALFree(r);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/******** <lalLaTeX file="LPCC"> ********
\noindent
Stabilizes a polynomial.
\subsubsection*{Prototype}
********* </lalLaTeX> ********/
/* <lalVerbatim> */
void LALPolystab(LALStatus *status,
		 REAL4Vector *a) {
/* </lalVerbatim> */
/******** <lalLaTeX file="LPCC"> ********
\subsubsection*{Description}
Reflects poles and zeroes of a inside the complex unit circle.

\subsubsection*{Uses}
\begin{verbatim}
...
\end{verbatim}
********* </lalLaTeX> ********/

  INITSTATUS (status, "LALPolystab", LPCC);
  ATTATCHSTATUSPTR (status);

  if(a->length > 1) {

    REAL4 *rtr, *rti, *ari;
    REAL4 aind1;
    INT4 i, j, m = a->length - 1;

    ari = (REAL4 *)LALCalloc(a->length, sizeof(REAL4));
    if(!ari) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

    rtr = (REAL4 *)LALCalloc(m, sizeof(REAL4));
    if(!rtr) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

    rti = (REAL4 *)LALCalloc(m, sizeof(REAL4));
    if(!rti) {ABORT(status, LPCH_EMEM, LPCH_MSGEMEM);}

    for(i=0;i<m+1;i++) {
      if(a->data[i] != 0) {break;
      }
    }
    aind1 = a->data[i];

    for(i=0; i<(int)a->length/2; i++) {
      REAL4 tmp;
      tmp = a->data[i];
      a->data[i] = a->data[a->length-i-1];
      a->data[a->length-i-1] = tmp;
    }

    rtr--;
    rti--;

    if(zrhqr(a->data,m,rtr,rti)) {
      ABORT(status, LPCH_ENUM, LPCH_MSGENUM);
    }

    for(i=1;i<=m;i++) {
      if(rtr[i] != 0 || rti[i] != 0) {
	REAL4 vs;
	REAL4 absv2 = rtr[i]*rtr[i] + rti[i]*rti[i];
	REAL4 absv = sqrt(absv2);
	if(absv >= 1.0) {vs=1.0;} else {vs=0.0;}
	rtr[i] = (1.0-vs)*rtr[i] + vs*rtr[i]/absv2;
	rti[i] = (1.0-vs)*rti[i] + vs*rti[i]/absv2;
      }
    }

    bzero(a->data,(m+1)*sizeof(REAL4));
    a->data[0]=1.0;

    for(j=1;j<=m;j++) {
      for(i=j;i>=1;i--) {
	a->data[i] -= rtr[j]*a->data[i-1] - rti[j]*ari[i-1];
	ari[i] -= rtr[j]*ari[i-1] + rti[j]*a->data[i-1];
      }
    }

    for(i=0; i<(int)a->length; i++) {
      a->data[i] *= aind1;
    }

    LALFree(ari);
    LALFree(rtr+1);
    LALFree(rti+1);

  }

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




#define FREERETURN {LALFree(h+1);LALFree(g+1); return 0;}

static INT4 toeplz(REAL4 r[], REAL4 x[], REAL4 y[], INT4 n) {

  INT4 j,k,m,m1,m2;
  REAL4 pp,pt1,pt2,qq,qt1,qt2,sd,sgd,sgn,shn,sxn;
  REAL4 *g,*h;

  r--;
  y--;
  x--;

  if (r[n] == 0.0) {
    fprintf(stderr,"toeplz-1 singular principal minor");
    return 1;
  }

  g = (REAL4 *)LALCalloc(n, sizeof(REAL4));
  h = (REAL4 *)LALCalloc(n, sizeof(REAL4));

  g--;
  h--;

  x[1]=y[1]/r[n];
  if (n == 1) FREERETURN
		g[1]=r[n-1]/r[n];
  h[1]=r[n+1]/r[n];
  for (m=1;m<=n;m++) {
    m1=m+1;
    sxn = -y[m1];
    sd = -r[n];
    for (j=1;j<=m;j++) {
      sxn += r[n+m1-j]*x[j];
      sd += r[n+m1-j]*g[m-j+1];
    }
    if (sd == 0.0) {
      fprintf(stderr,"toeplz-2 singular principal minor");
      return 1;
    }
    x[m1]=sxn/sd;
    for (j=1;j<=m;j++) x[j] -= x[m1]*g[m-j+1];
    if (m1 == n) FREERETURN
		   sgn = -r[n-m1];
    shn = -r[n+m1];
    sgd = -r[n];
    for (j=1;j<=m;j++) {
      sgn += r[n+j-m1]*g[j];
      shn += r[n+m1-j]*h[j];
      sgd += r[n+j-m1]*h[m-j+1];
    }
    if (sd == 0.0 || sgd == 0.0) {
      fprintf(stderr,"toeplz-3 singular principal minor");
      return 1;
    }
    g[m1]=sgn/sgd;
    h[m1]=shn/sd;
    k=m;
    m2=(m+1) >> 1;
    pp=g[m1];
    qq=h[m1];
    for (j=1;j<=m2;j++) {
      pt1=g[j];
      pt2=g[k];
      qt1=h[j];
      qt2=h[k];
      g[j]=pt1-pp*qt2;
      g[k]=pt2-pp*qt1;
      h[j]=qt1-qq*pt2;
      h[k--]=qt2-qq*pt1;
    }
  }

  fprintf(stderr,"toeplz - should not arrive here!");
  return 1;
}
#undef FREERETURN



/*
#define MAXM 50
*/

static INT4 zrhqr(REAL4 a[], INT4 m, REAL4 rtr[], REAL4 rti[])
{
	INT4 j,k;
	REAL4 **hess,xr,xi;
	INT4 MAXM = m+1;

	hess=(REAL4 **)LALCalloc(1+MAXM, sizeof(REAL4 *));
	for(j=0;j<1+MAXM;j++) {
	  hess[j] = (REAL4 *)LALCalloc(1+MAXM, sizeof(REAL4));
	}

	if (m > MAXM || a[m] == 0.0) {fprintf(stderr,"bad args in zrhqr"); return 1;}
	for (k=1;k<=m;k++) {
		hess[1][k] = -a[m-k]/a[m];
		for (j=2;j<=m;j++) hess[j][k]=0.0;
		if (k != m) hess[k+1][k]=1.0;
	}
	if(balanc(hess,m)) {return 1;}
	if(hqr(hess,m,rtr,rti)) {return 1;}
	for (j=2;j<=m;j++) {
		xr=rtr[j];
		xi=rti[j];
		for (k=j-1;k>=1;k--) {
			if (rtr[k] <= xr) break;
			rtr[k+1]=rtr[k];
			rti[k+1]=rti[k];
		}
		rtr[k+1]=xr;
		rti[k+1]=xi;
	}

	for(j=0;j<1+MAXM;j++) {
	  LALFree(hess[j]);
	}
	LALFree(hess);

	return 0;
}
#undef MAXM

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static INT4 hqr(REAL4 **a, INT4 n, REAL4 wr[], REAL4 wi[])
{
	INT4 nn,m,l,k,j,its,i,mmin;
	REAL4 z,y,x,w,v,u,t,s,r,q,p,anorm;

        r=q=p=0;

	anorm=0.0;
	for (i=1;i<=n;i++)
		for (j=IMAX(i-1,1);j<=n;j++)
			anorm += fabs(a[i][j]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if ((REAL4)(fabs(a[l][l-1]) + s) == s) break;
			}
			x=a[nn][nn];
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					}
					nn -= 2;
				} else {
				  if (its == 30) {fprintf(stderr,"Too many iterations in hqr"); return 1;}
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i][i] -= x;
						s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if ((REAL4)(u+v) == v) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i][i-2]=0.0;
						if (i != (m+2)) a[i][i-3]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}

	return 0;

}
#undef IMAX
#undef SIGN

#define RADIX 2.0

static INT4 balanc(REAL4 **a, INT4 n)
{
	INT4 last,j,i;
	REAL4 s,r,g,f,c,sqrdx;

	sqrdx=RADIX*RADIX;
	last=0;
	while (last == 0) {
		last=1;
		for (i=1;i<=n;i++) {
			r=c=0.0;
			for (j=1;j<=n;j++)
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			if (c && r) {
				g=r/RADIX;
				f=1.0;
				s=c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i][j] *= g;
					for (j=1;j<=n;j++) a[j][i] *= f;
				}
			}
		}
	}

	return 0;
}
#undef RADIX

