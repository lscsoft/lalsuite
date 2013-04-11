/*
*  Copyright (C) 2007 Matt Pitkin
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

/* Write function definitions for reading in the Crab ephemeris and */
/* computing the coefficients of the spin down model for f3 and f4  */

/* Matt Pitkin 10/03/04 */
/* LALTimingNoiseHeterodyne function has now been changed so that it
   doesn't have to heterodyne at the SSB, cos the data has already been
	 transformed to the SSB - soooo much easier! */

/* Matt Pitkin 26/03/04 */
/* The function's inputs params have be swithced about to conform to
	 the LAL Spec */

/**
         \file
	 \author M. D. Pitkin
         \brief
	 Functions for reading the Crab Ephemeris, interpolating between these values and
	 heterodyning the timing noise phase difference.

	 \heading{Description}

	 The first function will read the data from a file containing the Crab
	 pulsar ephemeris (<tt>crab_ephemeris.txt</tt>) which has four columns with
	 the Modified Julian Date (MJD) in Barycentre Dynamical Time (TDB),
	 the arrival time of the first pulse after
	 midnight on the given MJD in seconds, the frequency \f$f\f$ (Hz) at that time
	 and \f$\dot{f}\f$ (\f$10^{-10}\,\textrm{Hz}^2\f$) CrabEphem. The ephemeris values
	 are given on a roughly monthly basis. The function converts the MJD
	 times into GPS seconds including adding leap seconds. The ephemeris values
	 are available from 15 February 1982 so all leap seconds from that point are
	 added at the correct times - this is overkill as there is no data to be
	 analysed prior 2001, but I've included it just for completeness.

	 The second function will use the ephemeris values and create a fourth order
	 spline interpolation between consecutive monthly points. The simple phase
	 evolution of a pulsar signal can be described by the truncated Talyor series

	 \f{equation}{
	 \phi(t) = \phi_0 + f_0(t-t_0) + \frac{\dot{f_0}(t-t_0)^2}{2} + \dots,
	 \f}

	 To completely model the phase evolution of the Crab pulsar the higher order
	 frequency derivatives of the spin down model (\f$\ddot{f}\f$, \f$\dddot{f}\f$,
	 and eventually \f$\ddddot{f}\f$) need to be calculated for each point in the
	 ephemeris. This is done by fitting the spin-down model between consecutive
	 points using the known values of \f$f\f$,\f$\dot{f}\f$ and the phase \f$\phi\f$ at
	 either point as boundary conditions. The evolution of the phase, frequency
	 and higher order derivatives are given by the following matrix:

	 \anchor numatrix \f{equation}{\label{numatrix}
	 \left( \begin{array}{c} \phi \\ f \\ \dot{f} \\ \ddot{f} \\
	 \dddot{f} \\ \ddddot{f}
	 \end{array} \right) =
	 \left( \begin{array}{cccccc} 1 & t & \frac{1}{2}t^2 & \frac{1}{6}t^3 &
	 \frac{1}{24}t^4 & \frac{1}{120}t^5 \\
	 0 & 1 & t & \frac{1}{2}t^2 & \frac{1}{6}t^3 & \frac{1}{24}t^4 \\
	 0 & 0 & 1 & t & \frac{1}{2}t^2 & \frac{1}{6}t^3 \\
	 0 & 0 & 0 & 1 & t & \frac{1}{2}t^2 \\
	 0 & 0 & 0 & 0 & 1 & t \\
	 0 & 0 & 0 & 0 & 0 & 1 \\
	 \end{array} \right)
	 \left( \begin{array}{c} \phi_0 \\ f_0 \\ \dot{f_0} \\ \ddot{f_0} \\
	 \dddot{f_0} \\ \ddddot{f_0}
	 \end{array} \right)
	 \f}

	 Using the relation \f$\vec{a} = \mathbf{X}\vec{b}\f$, where \f$\vec{a}\f$ and
	 \f$\vec{b}\f$ are the vectors of \f$\phi\f$ and \f$\phi_0\f$ derivatives respectively
	 and \f$\mathbf{X}\f$ is the matrix of coefficients, it follows that \f$\vec{b} =
	 \mathbf{X}^{-1}\vec{a}\f$. The inverse matrix \f$\mathbf{X}^{-1}\f$ is

	 \anchor nu2matrix \f{equation}{\label{nu2matrix}
	 \mathbf{X}^{-1} =
	 \left( \begin{array}{cccccc} 1 & -t & \frac{1}{2}t^2 & -\frac{1}{6}t^3 &
	 \frac{1}{24}t^4 & -\frac{1}{120}t^5 \\
	 0 & 1 & -t & \frac{1}{2}t^2 & -\frac{1}{6}t^3 & \frac{1}{24}t^4 \\
	 0 & 0 & 1 & -t & \frac{1}{2}t^2 & -\frac{1}{6}t^3 \\
	 0 & 0 & 0 & 1 & -t & \frac{1}{2}t^2 \\
	 0 & 0 & 0 & 0 & 1 & -t \\
	 0 & 0 & 0 & 0 & 0 & 1 \\
	 \end{array} \right).
	 \f}

	 The values of \f$\phi\f$ and \f$\phi_0\f$ are unknown at each data point, but it is
	 known that they must be integer values. The phase at the first point
	 can be set to zero with the consecutive points having phases with integer
	 values. The exact number of then needs to be calculated. To get a constraint
	 on the phase at each data point, the boundary conditions for the known
	 parameters can be used to calculate the values of the unknowns \f$\ddot{f}\f$,
	 \f$\ddot{f_0}\f$, \f$\dddot{f}\f$, and \f$\dddot{f_0}\f$. Taking the \f$4\times 4\f$
	 matrices in the centres of\eqref{numatrix} and\eqref{nu2matrix} and vectors
	 of \f$f\f$ and its first three derivatives, one can then construct 4 equations
	 containing the 4 unknowns by equating the equations, when worked out at a
	 time \f$t\f$ halfway between consecutive points and solving them
	 simultaneously using the boundary conditions (by inverting the matrix in
	 \eqref{simult1} to give\eqref{inverse1} and solving as above). The values
	 could then be used to calculate a value of \f$\phi\f$ from\eqref{phi}.

	 \anchor simult1 \f{equation}{\label{simult1}
	 \left( \begin{array}{c} f_0 + \dot{f_0}t - f + \dot{f}t \\
	 \dot{f_0} - \dot{f} \\ 0 \\ 0 \end{array} \right) = \left(
	 \begin{array}{cccc}
	 -\frac{1}{2}t^2 & -\frac{1}{6}t^3 & \frac{1}{2}t^2 & -\frac{1}{6}t^3 \\
	 -t & -\frac{1}{2}t^2 & -t & \frac{1}{2}t^2 \\
	 -1 & -t & 1 & -t \\
	 0 & -1 & 0 & 1 \end{array} \right) \left( \begin{array}{c}
	 \ddot{f_0} \\ \dddot{f_0} \\ \ddot{f} \\ \dddot{f}
	 \end{array} \right).
	 \f}

	 \anchor inverse1 \f{equation}{\label{inverse1}
	 \mathbf{X}^{-1}
	 \left( \begin{array}{cccc}
	 -\frac{1}{2}t^2 & -\frac{1}{6}t^3 & \frac{1}{2}t^2 & -\frac{1}{6}t^3 \\
	 -t & -\frac{1}{2}t^2 & -t & \frac{1}{2}t^2 \\
	 -1 & -t & 1 & -t \\
	 0 & -1 & 0 & 1 \end{array} \right)
	 \f}

	 Equating the whole of equations\eqref{numatrix} and\eqref{nu2matrix} at a
	 point halfway inbetween consecutive points, gave the six equations necessary
	 to work out the six unknowns, which are shown in matrix form in equation
	 \eqref{simult2}.

	 \anchor simult2 \f{equation}{\label{simult2}
	 \left( \begin{array}{c} \phi_0 + f_0 t + \frac{\dot{f_0}}{2}t^2 - \phi + f t - \frac{\dot{f}}{2} \\ f_0 + \dot{f_0}t - f + \dot{f}t \\ \dot{f_0} - \dot{f} \\ 0 \\ 0 \\ 0 \end{array} \right) = \left( \begin{array}{cccccc}
	 -\frac{1}{6}t^3 & -\frac{1}{24}t^4 & -\frac{1}{120}t^5 & -\frac{1}{6}t^3 & \frac{1}{24}t^4 & \frac{1}{120}t^5 \\
	 -\frac{1}{2}t^2 & -\frac{1}{6}t^3 & -\frac{1}{24}t^4 & \frac{1}{2}t^2 & -\frac{1}{6}t^3 & \frac{1}{24}t^4 \\
	 -t & -\frac{1}{2}t^2 & -\frac{1}{6}t^3 & -t & \frac{1}{2}t^2 & -\frac{1}{6}t^3 \\
	 -1 & -t & -\frac{1}{2}t^2 & 1 & -t & \frac{1}{2}t^2 \\
	 0 & -1 & -t & 0 & 1 & -t \\
	 0 & 0 & -1 & 0 & 0 & 1 \end{array} \right) \left( \begin{array}{c}
	 \ddot{f_0} \\ \dddot{f_0} \\ \ddddot{f_0} \\ \ddot{f} \\ \dddot{f} \\ \ddddot{f}
	 \end{array} \right).
	 \f}

	 The matrix has an inverse,

	 \anchor inverse2 \f{equation}{\label{inverse2}
	 \mathbf{X}^{-1} =
	 \left( \begin{array}{cccccc}
	 -\frac{15}{2}\frac{1}{t^3} & -\frac{3}{2}\frac{1}{t^2} &
	 \frac{3}{4}\frac{1}{t} & \frac{1}{4}& -\frac{1}{16}t & -\frac{1}{16}t^2 \\
	 \frac{45}{2}\frac{1}{t^4} & \frac{3}{2}\frac{1}{t^3} &
	 -\frac{15}{4}\frac{1}{t^2} & -\frac{3}{4}\frac{1}{t} & \frac{7}{16} &
	 \frac{5}{16}t \\
	 -\frac{42}{2}\frac{1}{t^5} & 0 & \frac{15}{4}\frac{1}{t^3} & 0 &
	 -\frac{15}{16}\frac{1}{t} & -\frac{1}{2} \\
	 -\frac{15}{2}\frac{3}{2} & \frac{3}{2}\frac{1}{t} & \frac{3}{4}\frac{1}{t} &
	 -\frac{1}{4} & -\frac{1}{16}t & \frac{1}{16}t^2 \\
	 -\frac{45}{2}\frac{1}{t^4} & \frac{3}{2}\frac{1}{t^3} &
	 \frac{15}{4}\frac{1}{t^2} & -\frac{3}{4}\frac{1}{t} & -\frac{7}{16} &
	 \frac{5}{16}t \\
	 -\frac{45}{2}\frac{1}{t^5} & 0 & \frac{15}{4}\frac{1}{t^3} & 0 &
	 -\frac{15}{16}\frac{1}{t} & \frac{1}{2} \end{array} \right)
	 \f}

	 This leaves three equations for the values of the unknowns \f$\ddot{f_0}\f$,
	 \f$\dddot{f_0}\f$ and \f$\ddddot{f_0}\f$,

	 \f{equation}{
	 \ddot{f_0} =  -\frac{15}{2t^3}a - \frac{3}{2t^2}b + \frac{3}{4t}c,
	 \f}

	 \f{equation}{
	 \dddot{f_0} = \frac{45}{2t^4}a + \frac{3}{2t^3}b - \frac{15}{4t^2}c,
	 \f}

	 \f{equation}{
	 \ddddot{f_0} = -\frac{45}{2t^5}a + \frac{15}{4t^3}c,
	 \f}
	 where \f$a = \phi_0 + f_0 t + \frac{\dot{f_0}}{2}t^2 - \phi + f t -
	 \frac{\dot{f}}{2}\f$, \f$b = f_0 + \dot{f_0}t - f + \dot{f}t\f$, \f$c =
	 \dot{f_0} - \dot{f}\f$ and \f$t\f$ is the time halfway between ephemeris points.

	 The third function sets the values of the Crab pulsar spin-down to be used
	 to heterodyne a particular data point. The values are set to those that are
	 the first values prior to the time of the data point.

	 The fourth function calculates the phase difference between that used to
	 initially heterodyne the data and the phase including timing noise as
	 calculated using the interpolated ephemeris values. It then
	 removes this phase difference via another complex heterodyne,

	 \f{equation}{
	 B_\textrm{k final} = B_\textrm{k initial}e^{-i\Delta\phi}.
	 \f}

*/

#include <math.h>
#include <stdio.h>

#include <lal/LALConstants.h>
#include <lal/BinaryPulsarTiming.h>
#include "HeterodyneCrabPulsar.h"

void
LALGetCrabEphemeris	( LALStatus			*status,
			  CrabSpindownParamsInput	*output,
				GetCrabEphemerisInput *input )
{
  UINT4 i=0, j=0;

  REAL8 MJD;
  REAL8 MJDVec[1000];
	REAL8 GPStemp;

  /* array of the time of the pulse from the ephemeris  */
  /* minus the time of the first ephemeris data point */
  REAL8 GPStArr;
  REAL8 tArrVec[1000];

  REAL8 f0;
  REAL8 f1;

  FILE *fp;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Check Input Arguments */
  ASSERT(output != (CrabSpindownParamsInput *)NULL, status, HETERODYNECRABPULSARH_ENULLOUTPUT,
  HETERODYNECRABPULSARH_MSGENULLOUTPUT);

  ASSERT(input->filename != NULL, status, HETERODYNECRABPULSARH_EEPHEMERISFILENAME,
  HETERODYNECRABPULSARH_MSGEEPHEMERISFILENAME);

  if((fp = fopen(input->filename,"r")) == NULL){
    fprintf(stderr,"Unable to open file!");
    exit(1);
  }

  /* read in ephermeris data */
  while(!feof(fp)){
    fscanf(fp,"%lf%lf%lf%lf",&MJD,&GPStArr,&f0,&f1);
    output->f0->data[j] = f0;
    output->f1->data[j] = f1;
    MJDVec[j] = MJD;
    tArrVec[j] = GPStArr;
    j++;
  }

  fclose(fp);

  ASSERT(j != 0, status, HETERODYNECRABPULSARH_ENUMEPHEMERISDATA,
  HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA);

  output->numOfData = j;

/* check that the number of times in the ephemeris file (N2) is equal to the */
  /* number of lines read in (i-1) if not give error and exit                */
  /*if ( N2 != i-1 ){
    fprintf(stderr,"Error in ephemeris file, check or update\n");
    exit(1);
  }*/

  /* convert time in MJD to secs in TDB */
  for(i=0;i<j;i++){
    GPStemp = XLALTTMJDtoGPS(MJDVec[i]);

    MJDVec[i] = GPStemp + tArrVec[i];

    /*printf("%lf\n",fmod(MJDVec[i],50000));*/
    output->tArr->data[i] = MJDVec[i];
  }

  /* LALCheckMemoryLeaks();*/

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALComputeFreqDerivatives	( LALStatus			*status,
				  CrabSpindownParamsOutput	*output,
					CrabSpindownParamsInput	*input )
{
  REAL8 t, t1; /* temp time values to calculate freq derivs */
  REAL8 a, b, c; /* values used to calculate frequency derivs */
  REAL8 tempf2, tempf3; /* temp values of f2 and f3 */
  UINT4 i;
  REAL8 f0;
  REAL8 f1;
  REAL8 f2;  /* second deriv of nu */
  REAL8 f3;  /* third deriv of nu */
  REAL8 f4;  /* fourth deriv of nu */
  REAL8 tArr;	/* arrival time of pulse (GPS) */
  REAL8 phase0, phase; /* phase of signal */
  REAL8 phaseNew;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Check Input Arguments */
  ASSERT(input != (CrabSpindownParamsInput *)NULL, status, HETERODYNECRABPULSARH_ENULLINPUT,
  HETERODYNECRABPULSARH_MSGENULLINPUT);

  ASSERT(output != (CrabSpindownParamsOutput *)NULL, status, HETERODYNECRABPULSARH_ENULLOUTPUT,
  HETERODYNECRABPULSARH_MSGENULLOUTPUT);

  ASSERT(input->numOfData != 0, status, HETERODYNECRABPULSARH_ENUMEPHEMERISDATA,
  HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA);

  input->f1->data[0] = input->f1->data[0]*1.e-15;  /* convert f1 to correct units */
  phase0 = 0.0;
  phaseNew = 0.0;

  for (i=1;i<input->numOfData;i++){
    tArr = input->tArr->data[i-1];
    t = input->tArr->data[i]-input->tArr->data[i-1];
    t1 = t/2.0;
    input->f1->data[i] = input->f1->data[i]*1.0e-15;

    /* calculate values for spline */
    b = input->f0->data[i-1] + input->f1->data[i-1]*t1 - input->f0->data[i] +
    input->f1->data[i]*t1;
    c = input->f1->data[i-1] - input->f1->data[i];

    /* calculate second and third derivatives of f with constrained f0 and f1 */

    tempf2 = -(3.0/(2.0*t1*t1))*b - (1.0/(2.0*t1))*c;
    tempf3 = (3.0/(2.0*t1*t1*t1))*b;

    /* calculate the phase */
    phase = phase0 +input->f0->data[i-1]*t + (input->f1->data[i-1]/2.0)*t*t +
    (tempf2/6.0)*t*t*t + (tempf3/24.0)*t*t*t*t;

    /* round phase to nearest integer */
    /* (phase is tracked well enough using the values of f0, f1, f2 and f3
      calculated above to be able to constrain the phase after time t to the nearest
      integer) */
    phaseNew = floor(phase + 0.5);
    /* fprintf(stdout,"%lf\t%lf\n",phase,phaseNew); */

    /* recalculate spindown params with constrained phase */
    a = phase0 + input->f0->data[i-1]*t1 + (input->f1->data[i-1]/2.0)*t1*t1 -
    phaseNew + input->f0->data[i]*t1 - (input->f1->data[i]/2.0)*t1*t1;

    /* calculate the second, third and fourth derivatives of frequency */
    f0 = input->f0->data[i-1];
    f1 = input->f1->data[i-1];
    f2 = -(15.0/(2.0*t1*t1*t1))*a - (3.0/(2.0*t1*t1))*b + (3.0/(4.0*t1))*c;
    f3 = (45.0/(2.0*t1*t1*t1*t1))*a + (3.0/(2.0*t1*t1*t1))*b - (15.0/(4.0*t1*t1))*c;
    f4 = -(45.0/(2.0*t1*t1*t1*t1*t1))*a + (15.0/(4.0*t1*t1*t1))*c;

    output->f0->data[i-1] = f0;
    output->f1->data[i-1] = f1;
    output->f2->data[i-1] = f2;
    output->f3->data[i-1] = f3;
    output->f4->data[i-1] = f4;
    output->tArr->data[i-1] = tArr;
    output->numOfData = input->numOfData;
  }

  /* LALCheckMemoryLeaks(); */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALSetSpindownParams	( LALStatus  *status,
        ParamsForHeterodyne   *output,
        CrabSpindownParamsOutput  *input,
        LIGOTimeGPS   epoch)
{
  UINT4 i = 0;
  REAL8 dataEpoch;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Check Input Arguments */
  ASSERT(input != (CrabSpindownParamsOutput *)NULL, status, HETERODYNECRABPULSARH_ENULLINPUT,
  HETERODYNECRABPULSARH_MSGENULLINPUT);

  ASSERT(input->numOfData != 0, status, HETERODYNECRABPULSARH_ENUMEPHEMERISDATA,
  HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA);

  dataEpoch = (REAL8)epoch.gpsSeconds + ((REAL8)epoch.gpsNanoSeconds/1.0e9);

  for(i=0;i<input->numOfData;i++){
    if((dataEpoch<input->tArr->data[i+1])&&(dataEpoch>=input->tArr->data[i])){
      output->f0 = input->f0->data[i];
      output->f1 = input->f1->data[i];
      output->f2 = input->f2->data[i];
      output->f3 = input->f3->data[i];
      output->f4 = input->f4->data[i];

      output->epoch = input->tArr->data[i];

      break;
    }
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* function that performs the timing noise heterodyne but no longer at SSB */
/* also it now does one point at a time */
void
LALTimingNoiseHeterodyne	( LALStatus  *status,
          TNHeterodyneOutput  *output,
          TNHeterodyneInput *input,
          ParamsForHeterodyne *params,
          BarycenterInput baryinput,
          EarthState earth )
{
  REAL8 DPhi; /* phase difference to be removed */
  REAL8 t1;
  REAL8 t;
  REAL8 phi, phi0; /* phases with and without timing noise */

  EmissionTime emit;

  /* find the SSB barycenter time delay */
  XLALBarycenter(&emit, &baryinput, &earth);

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Check Input Arguments */
  ASSERT(input != (TNHeterodyneInput *)NULL, status, HETERODYNECRABPULSARH_ENULLINPUT,
  HETERODYNECRABPULSARH_MSGENULLINPUT);

  ASSERT(output != (TNHeterodyneOutput *)NULL, status, HETERODYNECRABPULSARH_ENULLOUTPUT,
  HETERODYNECRABPULSARH_MSGENULLOUTPUT);

  ASSERT(params->f0 > 0, status, HETERODYNECRABPULSARH_EINVALIDF0,
  HETERODYNECRABPULSARH_MSGEINVALIDF0);

  t = (REAL8)input->epoch.gpsSeconds+(REAL8)input->epoch.gpsNanoSeconds/1.e9;

  /* Dt between time of data point and the epoch of the FIRST data point */
  /* epoch of the FIRST data point is always the GPS time of that point  */
  t = t - (REAL8)input->t0 + emit.deltaT;

  /* Dt between time of data point and the epoch of that data point */
  t1 = t + (REAL8)input->t0 - (REAL8)params->epoch;

  /* calculate phase difference between signal with timing noise and one without */
  /* phi - phi0 */
  /* input f0 is already at GW freq (i.e. spin freq * 2) */
  phi0 = input->f0*t + 0.5*input->f1*t*t + (1.0/6.0)*input->f2*t*t*t +
    (input->f3/24.)*t*t*t*t;

  phi = 2.0*(params->f0*t1 + 0.5*params->f1*t1*t1 + (params->f2/6.0)*t1*t1*t1+
    (params->f3/24.0)*t1*t1*t1*t1 + (params->f4/120.0)*t1*t1*t1*t1*t1);

  DPhi = phi-phi0;
  DPhi = 2.0*(REAL8)LAL_PI*fmod(DPhi,1.0);

  output->Dphase = DPhi;
  output->phi0 = 2.*LAL_PI*fmod(phi0,1.0);
  output->phi1 = 2.*LAL_PI*fmod(phi,1.0);

  /* Heterodyne to remove timing noise */
  output->Vh = (REAL8)creal(input->Vh)*cos(-DPhi) - (REAL8)cimag(input->Vh)*sin(-DPhi)
    + I *((REAL8)creal(input->Vh)*sin(-DPhi) + (REAL8)cimag(input->Vh)*cos(-DPhi));

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
