/******************************* <lalVerbatim file="TBaryPtolemaicCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TBaryPtolemaic.c}}
\label{ss:TBaryPtolemaic.c}

Computes the barycentric arrival time of an incoming wavefront using a
circular model of the Earth's orbit.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TBaryPtolemaicCP}
\idx{LALTBaryPtolemaic()}
\idx{LALDTBaryPtolemaic()}

\subsubsection*{Description}

These routines compute the barycentric time transformation and its
derivatives.  That is, if a signal originating from a right ascension
$\alpha$ and declination $\delta$ on the sky and arrives at the
detector at a time $t$, then it will pass the centre of the solar
system at a time $t_b(t,\alpha,\delta)$.

The routines obey the calling convention presented in the header
\verb@PulsarTimes.h@, with $n=2$ variable parameters
$\lambda^1=\alpha$, $\lambda^2=\delta$ (both measured in radians).
The constant parameter fields used by these routines are
\verb@constants->tAutumn@, \verb@constants->tMidnight@,
\verb@constants->latitude@, and \verb@constants->longitude@.

Note that \verb@*variables@ must have a length of at least 3, and can
be longer, but values beyond the third are ignored.  
\verb@*dtBary@ must be at least of length 1, and the number of derivatives 
computed is determined by \verb@*dtBary->length-1@. 
All elements beyond the fourth will be set to zero. 

\subsubsection*{Algorithm}

Let $\mathbf{\hat{n}}(\alpha,\delta)$ be the unit vector to the source
on the sky, and $\mathbf{x}(t)$ be the position of the detector
relative to the solar system barycentre.  Then, ignoring relativistic
effects and considering only the Roemer time delay, the barycentred
time is:
$$
t_b(t,\alpha,\delta) = t + \frac{\mathbf{x}(t)}{c}\cdot
	\mathbf{\hat{n}}(\alpha,\delta) \; .
$$
Choosing a right-handed coordinate system with $\mathbf{\hat{e}}_x$
toward the vernal equinox and $\mathbf{\hat{e}}_z$ toward celestial
North, the components of $\mathbf{\hat{n}}(\alpha,\delta)$ can be
written as:
\begin{eqnarray}
n_x & = & \cos\alpha\cos\delta \; , \nonumber\\
n_y & = & \sin\alpha\cos\delta \; , \nonumber\\
n_z & = & \sin\delta           \; .
\label{eq:n-alphadelta}
\end{eqnarray}

In a first-order Ptolemaic model of the solar system, the rotation of
the Earth \emph{and} its orbit can be treated as simple circular
motions, inclined to one another at an angle $i=23.5^\circ$.  This
approximation is insufficient for actual signal demodulation, since it
can result in timing errors of several seconds due to the Earth's
orbital eccentricity.  However, the \emph{derivatives} of $t_b$ will
not be off by more than a percent or so, which is good enough for,
say, computing a parameter-space metric.  We define angles of rotation
and revolution
$\theta_\mathrm{rot}(t)=2\pi(t-t_\mathrm{midnight})/P_\mathrm{rot}$,
$\theta_\mathrm{rev}(t)=2\pi(t-t_\mathrm{autumn})/P_\mathrm{rev}$,
where $P_\mathrm{rot}$ is a sidereal day, $P_\mathrm{rev}$ a sidereal
year, $t_\mathrm{midnight}$ the time of a sidereal midnight at the
prime meridian, and $t_\mathrm{autumn}$ is the time of an autumnal
equinox.  In other words, $\theta_\mathrm{rot}=0$ corresponds to the
Earth lying on the positive $\mathbf{\hat{e}}_x$ axis and
$\theta_\mathrm{rev}=0$ corresponds to a point with zero longitude
lying in the $\mathbf{\hat{e}}_x$ direction from the Earth's axis.  If
the detector has a north latitude $\lambda$ and an east longitude $l$,
then the detector position is:
\begin{eqnarray}
x & = & R\cos\theta_\mathrm{rev} +
	r\cos(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
y & = & R\sin\theta_\mathrm{rev}\cos i +
	r\sin(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
z & = & R\sin\theta_\mathrm{rev}\sin i +
	r\sin\lambda \; , \nonumber
\end{eqnarray}
where $R$ is the Earth-Sun distance and $r$ is the Earth's radius.
The time dependence of $\mathbf{x}$ is implicit in the time dependence
of $\theta_\mathrm{rot}$ and $\theta_\mathrm{rev}$.

Differentiating with respect to $t$, we obtain:
$$
\frac{\partial t_b(t,\alpha,\delta)}{\partial t} = 1 +
	\frac{\mathbf{v}(t)}{c}\cdot
	\mathbf{\hat{n}}(\alpha,\delta) \; ,
$$
where:
\begin{eqnarray}
v_x & = & -V\sin\theta_\mathrm{rev} -
	v\sin(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
v_y & = & V\cos\theta_\mathrm{rev}\cos i +
	v\cos(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
v_z & = & V\cos\theta_\mathrm{rev}\sin i \; , \nonumber
\end{eqnarray}
$V=2\pi R/P_\mathrm{rev}$ is the Earth's orbital velocity, and $v=2\pi
r/P_\mathrm{rot}$ is the Earth's equatorial rotational velocity.

Differentiating with respect to $\alpha$ gives:
$$
\frac{\partial t_b(t,\alpha,\delta)}{\partial\alpha} =
	\frac{\mathbf{x}(t)}{c}\cdot
	\frac{\partial\mathbf{\hat{n}}(\alpha,\delta)}
		{\partial\alpha} \; ,
$$
where $\partial\mathbf{n}/\partial\alpha$ is easily obtained from
Eqs.~\ref{eq:n-alphadelta}.  Similarly for $\delta$.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TBaryPtolemaicCV}}

******************************************************* </lalLaTeX> */

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>
#include<lal/PulsarTimes.h>


/* Uncomment to remove orbital motion.
#define LAL_AU_SI 0.0 */

NRCSID(TBARYPTOLEMAICC,"$Id$");

/* <lalVerbatim file="TBaryPtolemaicCP"> */
void
LALTBaryPtolemaic( LALStatus             *stat,
		   REAL8                 *tBary,
		   REAL8Vector           *variables,
		   PulsarTimesParamStruc *constants )
{ /* </lalVerbatim> */
  INT4 i;      /* An index. */
  REAL8 t;     /* A time variable. */
  REAL8 ra;    /* Source right ascension, in radians. */
  REAL8 dec;   /* Source declination, in radians. */
  REAL8 rot;   /* Earth rotation angle, in radians. */
  REAL8 rev;   /* Earth revolution angle, in radians. */
  REAL8 *data; /* Pointer to vector data. */
  REAL8 x[3];  /* Detector position vector, in s. */
  REAL8 n[3];  /* Source position unit vector, or derivative. */

  /* Some local constants. */
  REAL8 tRot=LAL_REARTH_SI/LAL_C_SI;
  REAL8 tRev=LAL_AU_SI/LAL_C_SI;
  REAL8 cosLat;
  REAL8 cosi=cos(LAL_IEARTH);
  REAL8 sini=sin(LAL_IEARTH);

  INITSTATUS(stat,"TBaryPtolemaic",TBARYPTOLEMAICC);

  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if(lalDebugLevel){

    /* Make sure parameter structures and their fields exist. */
    ASSERT(tBary,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);

    /* Make sure array sizes are consistent. */
    ASSERT(variables->length>2,stat,
	   PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
  }
#endif

  /* Set some temporary variables. */
  data=variables->data;
  t=*(data++);
  ra=*(data++);
  dec=*(data);
  rot=LAL_TWOPI*(t-constants->tMidnight)/LAL_DAYSID_SI
    +constants->longitude;
  rev=LAL_TWOPI*(t-constants->tAutumn)/LAL_YRSID_SI;
  cosLat=cos(constants->latitude);

  /* Get detector position. */
  x[0]=tRot*cosLat*cos(rot)+tRev*cos(rev);
  x[1]=tRot*cosLat*sin(rot)+tRev*sin(rev)*cosi;
  x[2]=tRot*sin(constants->latitude)+tRev*sin(rev)*sini;

  /* Get source position unit vector. */
  n[0]=cos(ra)*cos(dec);
  n[1]=sin(ra)*cos(dec);
  n[2]=sin(dec);

  /* Compute barycentred time. */
  i=2;
  while(i--)
    t+=x[i]*n[i];
  *tBary=t;

  RETURN(stat);
}


/* <lalVerbatim file="TBaryPtolemaicCP"> */
void
LALDTBaryPtolemaic( LALStatus             *stat,
		    REAL8Vector           *dtBary,
		    REAL8Vector           *variables,
		    PulsarTimesParamStruc *constants )
{ /* </lalVerbatim> */
  INT4 i;      /* An index. */
  REAL8 t;     /* A time variable. */
  REAL8 ra;    /* Source right ascension, in radians. */
  REAL8 dec;   /* Source declination, in radians. */
  REAL8 rot;   /* Earth rotation angle, in radians. */
  REAL8 rev;   /* Earth revolution angle, in radians. */
  REAL8 *data; /* Pointer to vector data. */
  REAL8 x[3];  /* Detector position vector, in s. */
  REAL8 v[3];  /* Detector velocity vector, in c. */
  REAL8 n[3];  /* Source position unit vector, or derivative. */

  /* Some local constants. */
  REAL8 tRot=LAL_REARTH_SI/LAL_C_SI;
  REAL8 tRev=LAL_AU_SI/LAL_C_SI;
  REAL8 vRot=LAL_TWOPI*tRot/LAL_DAYSID_SI;
  REAL8 vRev=LAL_TWOPI*tRev/LAL_YRSID_SI;
  REAL8 cosLat;
  REAL8 cosi=cos(LAL_IEARTH);
  REAL8 sini=sin(LAL_IEARTH);

  UINT4 numDeriv; /* number of variables wrt which to compute derivatives */

  INITSTATUS(stat,"DTBaryPtolemaic",TBARYPTOLEMAICC);

  numDeriv = dtBary->length - 1;

  /* Make sure parameter structures and their fields exist. */
  ASSERT(dtBary,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(dtBary->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(variables,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(variables->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(constants,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  
  /* Make sure array sizes are consistent. */
  ASSERT( dtBary->length >= 1,stat, PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
  /* need at least (and exactly) [t, alpha, delta] */
  ASSERT( variables->length >= 3, stat, PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);

  /* Set some temporary variables. */
  data=variables->data;
  t=*(data++);
  ra=*(data++);
  dec=*(data);
  rot=LAL_TWOPI*(t-constants->tMidnight)/LAL_DAYSID_SI
    +constants->longitude;
  rev=LAL_TWOPI*(t-constants->tAutumn)/LAL_YRSID_SI;
  cosLat=cos(constants->latitude);
  vRot*=cosLat;

  /* Get detector position. */
  x[0]=tRot*cosLat*cos(rot)+tRev*cos(rev);
  x[1]=tRot*cosLat*sin(rot)+tRev*sin(rev)*cosi;
  x[2]=tRot*sin(constants->latitude)+tRev*sin(rev)*sini;

  /* Get detector velocity. */
  v[0]=-vRot*sin(rot)-vRev*sin(rev);
  v[1]=vRot*cos(rot)+vRev*cos(rev)*cosi;
  v[2]=vRev*cos(rev)*sini;

  /* Get source position unit vector. */
  n[0]=cos(ra)*cos(dec);
  n[1]=sin(ra)*cos(dec);
  n[2]=sin(dec);

  /* Initialize dtBary. */
  data=dtBary->data;
  memset(data,0,dtBary->length*sizeof(REAL8));

  /* Compute barycentred time. */
  i=3;
  while(i--)
    t+=x[i]*n[i];
  *(data++)=t;

  /* only compute the following derivatives if requested */
  if ( numDeriv >= 1 )
    {
      /* Compute time derivative. */
      t=1.0;
      i=3;
      while(i--)
	t+=v[i]*n[i];
      *(data++)=t;
    }
  if ( numDeriv >= 2 )
    {
      /* Get the right-ascension-derivative of the source position and
	 barycentred time. */
      n[0]=-sin(ra)*cos(dec);
      n[1]=cos(ra)*cos(dec);
      t=0;
      i=2;
      while(i--)
	t+=x[i]*n[i];
      *(data++)=t;
    }
  if ( numDeriv >= 3 )
    {
      /* Get the declination-derivative of the source position and
	 barycentred time. */
      n[0]=-cos(ra)*sin(dec);
      n[1]=-sin(ra)*sin(dec);
      n[2]=cos(dec);
      t=0;
      i=3;
      while(i--)
	t+=x[i]*n[i];
      *data=t;
    }

  RETURN(stat);

} /* LALDTBaryPtolemaic() */
