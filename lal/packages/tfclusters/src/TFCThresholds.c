/*-----------------------------------------------------------------------*
 *
 * File Name: TFCThresholds.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------*/
 

/******** <lalVerbatim file="TFCThresholdsCV"> ********
Author: Sylvestre, J
$Id$
********* </lalVerbatim> ********/

#include "lal/LALRCSID.h"

NRCSID (TFCLUSTERSC, "$Id$");


#include <lal/TFCThresholds.h>
#include <math.h>

/* local util functions */
static REAL8 RiceThreshold(REAL8 bpp, REAL8 Q, REAL8 P0, REAL4 eGoal);

static REAL4 NewtonRoot(void (*funcd)(REAL4, REAL4 *, REAL4 *, REAL8, REAL8, REAL8), REAL4, REAL4, REAL4, REAL8, REAL8, REAL8);

static void RiceWrapper(REAL4, REAL4 *, REAL4 *, REAL8, REAL8, REAL8);

static REAL4 Rice(REAL4, REAL8, REAL8);

static REAL4 RiceInt(REAL4, REAL4, REAL8, REAL8);

static REAL8 i0(REAL8 x);

static REAL4 RombergInt(REAL4 (*)(REAL4, REAL8, REAL8), REAL4, REAL4, REAL8, REAL8);

static void PolynomialInterpolation(REAL4 xa[], REAL4 ya[], INT4, REAL4, REAL4 *, REAL4 *);

static REAL4 Trapezoid(REAL4 (*)(REAL4, REAL8, REAL8), REAL4, REAL4, INT4, REAL8, REAL8);

static REAL4 *vector(long nl, long nh);

static void free_vector(REAL4 *v, long nl);

/******** <lalLaTeX file="TFCThresholdsC"> ********
\noindent
Computes thresholds on the power from a best fit to the Rice distribution.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void
LALTFCRiceThreshold ( LALStatus *status,
		      REAL4* rho,
		      RiceThresholdParams* thr
		      );
}
\idx{LALTFCRiceThreshold()}

\subsubsection*{Description}
Computes the thresholds on the power at every frequency, using a best fit to the Rice distribution. 

\subsubsection*{Notes}
\begin{itemize}
\item \texttt{rho} must be pointing to allocated memory. 
\end{itemize}
\vfill{\footnotesize\input{TFCThresholdsCV}}
********* </lalLaTeX> ********/

#define EPS 1.0e-6 /* error on Romberg integration */
static BOOLEAN numError = 0;

void
LALTFCRiceThreshold ( LALStatus *status,
		      REAL4* rho,
		      RiceThresholdParams* thr
		      )
{
  UINT4 i;
  

  INITSTATUS (status, "LALComputeSpectrogram", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);


  ASSERT(rho, status, TFCTHRESHOLDSH_ENULLP, TFCTHRESHOLDSH_MSGENULLP);
  ASSERT(thr, status, TFCTHRESHOLDSH_ENULLP, TFCTHRESHOLDSH_MSGENULLP);

  ASSERT(thr->meanRe, status, TFCTHRESHOLDSH_ENULLP, TFCTHRESHOLDSH_MSGENULLP);
  ASSERT(thr->meanIm, status, TFCTHRESHOLDSH_ENULLP, TFCTHRESHOLDSH_MSGENULLP);
  ASSERT(thr->varRe, status, TFCTHRESHOLDSH_ENULLP, TFCTHRESHOLDSH_MSGENULLP);
  ASSERT(thr->varIm, status, TFCTHRESHOLDSH_ENULLP, TFCTHRESHOLDSH_MSGENULLP);

  ASSERT(thr->bpp >= 0.0 && thr->bpp <= 1.0, status, TFCTHRESHOLDSH_E01, TFCTHRESHOLDSH_MSGE01);

  ASSERT(thr->eGoal > EPS, status, TFCTHRESHOLDSH_EEGOAL, TFCTHRESHOLDSH_MSGEEGOAL);

  for(i=0; i<thr->nFreq; i++) {

    rho[i] = RiceThreshold(thr->bpp, 
			   thr->meanRe[i]*thr->meanRe[i] + 
			   thr->meanIm[i]*thr->meanIm[i], 
			   thr->varRe[i] + thr->varIm[i],
			   thr->eGoal);

    if(numError) {ABORT(status,TFCTHRESHOLDSH_ENERR,TFCTHRESHOLDSH_MSGENERR);}
  }

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




/**********************UTILITY FUNCTIONS******************************/




static REAL8 RiceThreshold(REAL8 bpp, REAL8 Q, REAL8 P0, REAL4 eGoal) {
  /* computes threshold to apply on power to have black pixel probability
     bpp, for a Rice distribution with scale P0 and non-centrality Q */

  REAL8 P;

  P = NewtonRoot(RiceWrapper, P0+Q, P0+Q-2.0*log(bpp)*sqrt(P0*(P0+2.0*Q)), eGoal, P0, Q, bpp);

  return P;

}



static void RiceWrapper(REAL4 P, REAL4 *f, REAL4 *df, REAL8 P0, REAL8 Q, REAL8 bpp) {

  *f = RiceInt(bpp, P, Q, P0);
  *df = Rice(P, Q, P0);
}

static REAL4 Rice(REAL4 Pt, REAL8 Q, REAL8 P0) {
  /* computes the value of the Rice distribution */
  REAL8 p, P = (REAL8)Pt;

  p = i0(2.0*sqrt(P*Q)/P0) * exp(-(P+Q)/P0) / P0;

  return p;

}

static REAL4 RiceInt(REAL4 bpp, REAL4 P, REAL8 Q, REAL8 P0) {
  /* computes the error from goal bpp */

  REAL4 er;

  er = 1.0 - bpp - RombergInt(Rice, 0, P, Q, P0);

  return er;

}




#define JMAX 30
#define JMAXP (JMAX+1)
#define K 5

static REAL4 RombergInt(REAL4 (*func)(REAL4, REAL8, REAL8), REAL4 a, REAL4 b,
            REAL8 P0, REAL8 Q)
{
        REAL4 ss,dss;
        REAL4 s[JMAXP],h[JMAXP+1];
        INT4 j;

        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
                s[j]=Trapezoid(func,a,b,j,P0,Q);
                if (j >= K) {
                        PolynomialInterpolation(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        if (fabs(dss) <= EPS*fabs(ss)) return ss;
                }
                h[j+1]=0.25*h[j];
        }
	numError = 1;
        return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K


#define FUNC(x) ((*func)(x, P0, Q))

static REAL4 Trapezoid(REAL4 (*func)(REAL4, REAL8, REAL8), REAL4 a, REAL4 b, INT4 n, REAL8 P0, REAL8 Q)
{
        REAL4 x,tnm,sum,del;
        static REAL4 s;
        INT4 it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}
#undef FUNC



static void PolynomialInterpolation(REAL4 xa[], REAL4 ya[], INT4 n, REAL4 x, REAL4 *y, REAL4 *dy)
{

        INT4 i,m,ns=1;
        REAL4 den,dif,dift,ho,hp,w;
        REAL4 *c,*d;

        dif=fabs(x-xa[1]);
        c=vector(1,n);
        d=vector(1,n);
        for (i=1;i<=n;i++) {
                if ( (dift=fabs(x-xa[i])) < dif) {
                        ns=i;
                        dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
        }
        *y=ya[ns--];
        for (m=1;m<n;m++) {
                for (i=1;i<=n-m;i++) {
                        ho=xa[i]-x;
                        hp=xa[i+m]-x;
                        w=c[i+1]-d[i];
                        if ( (den=ho-hp) == 0.0) {
			  numError = 1;
			  return;
			}
                        den=w/den;
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
        }
        free_vector(d,1);
        free_vector(c,1);
}


static REAL4 *vector(long nl, long nh)
/* allocate a REAL4 vector with subscript range v[nl..nh] */
{
        REAL4 *v;

        v=(REAL4 *)malloc((size_t) ((nh-nl+2)*sizeof(REAL4)));
        return v-nl+1;
}

static void free_vector(REAL4 *v, long nl)
/* free a REAL4 vector allocated with vector() */
{
        free((char *) (v+nl-1));
}





#define MAXIT 1000

static REAL4 NewtonRoot(void (*funcd)(REAL4, REAL4 *, REAL4 *, REAL8, REAL8, REAL8), 
             REAL4 x1, REAL4 x2,
             REAL4 xacc,
             REAL8 P0, REAL8 Q, REAL8 bpp)
{
        INT4 j;
        REAL4 df,dx,dxold,f,fh,fl;
        REAL4 temp,xh,xl,rts;

        (*funcd)(x1,&fl,&df,P0,Q,bpp);
        (*funcd)(x2,&fh,&df,P0,Q,bpp);
        if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
	  numError = 1;
	  return 0.0;
	}

        if (fl == 0.0) return x1;
        if (fh == 0.0) return x2;
        if (fl < 0.0) {
                xl=x1;
                xh=x2;
        } else {
                xh=x1;
                xl=x2;
        }
        rts=0.5*(x1+x2);
        dxold=fabs(x2-x1);
        dx=dxold;
        (*funcd)(rts,&f,&df,P0,Q,bpp);
        for (j=1;j<=MAXIT;j++) {
                if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
                        || (fabs(2.0*f) > fabs(dxold*df))) {
                        dxold=dx;
                        dx=0.5*(xh-xl);
                        rts=xl+dx;
                        if (xl == rts) return rts;
                } else {
                        dxold=dx;
                        dx=f/df;
                        temp=rts;
                        rts -= dx;
                        if (temp == rts) return rts;
                }
                if (fabs(dx) < xacc) return rts;
                (*funcd)(rts,&f,&df,P0,Q,bpp);
                if (f < 0.0)
                        xl=rts;
                else
                        xh=rts;
        }
	numError = 1;
        return 0.0;
}
#undef MAXIT



/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */

static REAL8 A[] =
{
-4.41534164647933937950E-18,
 3.33079451882223809783E-17,
-2.43127984654795469359E-16,
 1.71539128555513303061E-15,
-1.16853328779934516808E-14,
 7.67618549860493561688E-14,
-4.85644678311192946090E-13,
 2.95505266312963983461E-12,
-1.72682629144155570723E-11,
 9.67580903537323691224E-11,
-5.18979560163526290666E-10,
 2.65982372468238665035E-9,
-1.30002500998624804212E-8,
 6.04699502254191894932E-8,
-2.67079385394061173391E-7,
 1.11738753912010371815E-6,
-4.41673835845875056359E-6,
 1.64484480707288970893E-5,
-5.75419501008210370398E-5,
 1.88502885095841655729E-4,
-5.76375574538582365885E-4,
 1.63947561694133579842E-3,
-4.32430999505057594430E-3,
 1.05464603945949983183E-2,
-2.37374148058994688156E-2,
 4.93052842396707084878E-2,
-9.49010970480476444210E-2,
 1.71620901522208775349E-1,
-3.04682672343198398683E-1,
 6.76795274409476084995E-1
};


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */


static REAL8 B[] =
{
-7.23318048787475395456E-18,
-4.83050448594418207126E-18,
 4.46562142029675999901E-17,
 3.46122286769746109310E-17,
-2.82762398051658348494E-16,
-3.42548561967721913462E-16,
 1.77256013305652638360E-15,
 3.81168066935262242075E-15,
-9.55484669882830764870E-15,
-4.15056934728722208663E-14,
 1.54008621752140982691E-14,
 3.85277838274214270114E-13,
 7.18012445138366623367E-13,
-1.79417853150680611778E-12,
-1.32158118404477131188E-11,
-3.14991652796324136454E-11,
 1.18891471078464383424E-11,
 4.94060238822496958910E-10,
 3.39623202570838634515E-9,
 2.26666899049817806459E-8,
 2.04891858946906374183E-7,
 2.89137052083475648297E-6,
 6.88975834691682398426E-5,
 3.36911647825569408990E-3,
 8.04490411014108831608E-1
};

static REAL8 chbevl(REAL8 x, REAL8 array[], INT4 n);

static REAL8 i0(REAL8 x)
{
REAL8 y;

if( x < 0 )
        x = -x;
if( x <= 8.0 )
        {
        y = (x/2.0) - 2.0;
        return( exp(x) * chbevl( y, A, 30 ) );
        }

return(  exp(x) * chbevl( 32.0/x - 2.0, B, 25 ) / sqrt(x) );

}




static REAL8 i0e(REAL8 x)
{
REAL8 y;

if( x < 0 )
        x = -x;
if( x <= 8.0 )
        {
        y = (x/2.0) - 2.0;
        return( chbevl( y, A, 30 ) );
        }

return(  chbevl( 32.0/x - 2.0, B, 25 ) / sqrt(x) );

}


static REAL8 chbevl(REAL8 x, REAL8 array[], INT4 n)
{
REAL8 b0, b1, b2, *p;
INT4 i;

p = array;
b0 = *p++;
b1 = 0.0;
i = n - 1;

do
        {
        b2 = b1;
        b1 = b0;
        b0 = x * b1  -  b2  + *p++;
        }
while( --i );

return( 0.5*(b0-b2) );
}

#undef EPS
