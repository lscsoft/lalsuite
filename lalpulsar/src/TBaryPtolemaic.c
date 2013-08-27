/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton
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

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>
#include<lal/PulsarTimes.h>


/* Uncomment to remove orbital motion.
#define LAL_AU_SI 0.0 */

/** \file
    \author Creighton, T. D.
    \ingroup PulsarTimes_h
    \brief Compute the barycentric arrival time of an incoming wavefront using a circular model of the Earth's orbit.


These routines compute the barycentric time transformation and its derivatives.
That is, if a signal originating from a right ascension
\f$\alpha\f$ and declination \f$\delta\f$ on the sky and arrives at the
detector at a time \f$t\f$, then it will pass the centre of the solar
system at a time \f$t_b(t,\alpha,\delta)\f$.

The routines obey the calling convention presented in the module
\ref PulsarTimes_h, with \f$n=2\f$ variable parameters
\f$\lambda^1=\alpha\f$, \f$\lambda^2=\delta\f$ (both measured in radians).
The constant parameter fields used by these routines are
<tt>constants->tAutumn</tt>, <tt>constants->tMidnight</tt>,
<tt>constants->latitude</tt>, and <tt>constants->longitude</tt>.

Note that <tt>*variables</tt> must have a length of at least 3, and can
be longer, but values beyond the third are ignored.
<tt>*dtBary</tt> must be at least of length 1, and the number of derivatives
computed is determined by <tt>*dtBary->length-1</tt>.
All elements beyond the fourth will be set to zero.

<b>Algorithm</b>

Let \f$\mathbf{\hat{n}}(\alpha,\delta)\f$ be the unit vector to the source
on the sky, and \f$\mathbf{x}(t)\f$ be the position of the detector
relative to the solar system barycentre.  Then, ignoring relativistic
effects and considering only the Roemer time delay, the barycentred
time is:
\f[
t_b(t,\alpha,\delta) = t + \frac{\mathbf{x}(t)}{c}\cdot
	\mathbf{\hat{n}}(\alpha,\delta) \; .
\f]
Choosing a right-handed coordinate system with \f$\mathbf{\hat{e}}_x\f$
toward the vernal equinox and \f$\mathbf{\hat{e}}_z\f$ toward celestial
North, the components of \f$\mathbf{\hat{n}}(\alpha,\delta)\f$ can be
written as:
\anchor eq_n-alphadelta
\f{eqnarray}{
n_x & = & \cos\alpha\cos\delta \; , \nonumber\\
n_y & = & \sin\alpha\cos\delta \; , \nonumber\\
n_z & = & \sin\delta           \; .
\tag{eq_n-alphadelta}
\f}

In a first-order Ptolemaic model of the solar system, the rotation of
the Earth \e and its orbit can be treated as simple circular
motions, inclined to one another at an angle \f$i=23.5^\circ\f$.  This
approximation is insufficient for actual signal demodulation, since it
can result in timing errors of several seconds due to the Earth's
orbital eccentricity.  However, the \e derivatives of \f$t_b\f$ will
not be off by more than a percent or so, which is good enough for,
say, computing a parameter-space metric.  We define angles of rotation
and revolution
\f$\theta_\mathrm{rot}(t)=2\pi(t-t_\mathrm{midnight})/P_\mathrm{rot}\f$,
\f$\theta_\mathrm{rev}(t)=2\pi(t-t_\mathrm{autumn})/P_\mathrm{rev}\f$,
where \f$P_\mathrm{rot}\f$ is a sidereal day, \f$P_\mathrm{rev}\f$ a sidereal
year, \f$t_\mathrm{midnight}\f$ the time of a sidereal midnight at the
prime meridian, and \f$t_\mathrm{autumn}\f$ is the time of an autumnal
equinox.  In other words, \f$\theta_\mathrm{rot}=0\f$ corresponds to the
Earth lying on the positive \f$\mathbf{\hat{e}}_x\f$ axis and
\f$\theta_\mathrm{rev}=0\f$ corresponds to a point with zero longitude
lying in the \f$\mathbf{\hat{e}}_x\f$ direction from the Earth's axis.  If
the detector has a north latitude \f$\lambda\f$ and an east longitude \f$l\f$,
then the detector position is:
\f{eqnarray}{
x & = & R\cos\theta_\mathrm{rev} +
	r\cos(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
y & = & R\sin\theta_\mathrm{rev}\cos i +
	r\sin(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
z & = & R\sin\theta_\mathrm{rev}\sin i +
	r\sin\lambda \; , \nonumber
\f}
where \f$R\f$ is the Earth-Sun distance and \f$r\f$ is the Earth's radius.
The time dependence of \f$\mathbf{x}\f$ is implicit in the time dependence
of \f$\theta_\mathrm{rot}\f$ and \f$\theta_\mathrm{rev}\f$.

Differentiating with respect to \f$t\f$, we obtain:
\f[
\frac{\partial t_b(t,\alpha,\delta)}{\partial t} = 1 +
	\frac{\mathbf{v}(t)}{c}\cdot
	\mathbf{\hat{n}}(\alpha,\delta) \; ,
\f]
where:
\f{eqnarray}{
v_x & = & -V\sin\theta_\mathrm{rev} -
	v\sin(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
v_y & = & V\cos\theta_\mathrm{rev}\cos i +
	v\cos(\theta_\mathrm{rot}+l)\cos\lambda \; , \nonumber\\
v_z & = & V\cos\theta_\mathrm{rev}\sin i \; , \nonumber
\f}
\f$V=2\pi R/P_\mathrm{rev}\f$ is the Earth's orbital velocity, and \f$v=2\pi
r/P_\mathrm{rot}\f$ is the Earth's equatorial rotational velocity.

Differentiating with respect to \f$\alpha\f$ gives:
\f[
\frac{\partial t_b(t,\alpha,\delta)}{\partial\alpha} =
	\frac{\mathbf{x}(t)}{c}\cdot
	\frac{\partial\mathbf{\hat{n}}(\alpha,\delta)}
		{\partial\alpha} \; ,
\f]
where \f$\partial\mathbf{n}/\partial\alpha\f$ is easily obtained from
Eqs.\eqref{eq_n-alphadelta}.  Similarly for \f$\delta\f$.

\heading{Uses}
\code
lalDebugLevel
\endcode
*/
/*@{*/

void
LALTBaryPtolemaic( LALStatus             *stat,
		   REAL8                 *tBary,
		   REAL8Vector           *variables,
		   PulsarTimesParamStruc *constants )
{
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

  INITSTATUS(stat);

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


void
LALDTBaryPtolemaic( LALStatus             *stat,
		    REAL8Vector           *dtBary,
		    REAL8Vector           *variables,
		    PulsarTimesParamStruc *constants )
{
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

  INITSTATUS(stat);

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
/*@}*/
