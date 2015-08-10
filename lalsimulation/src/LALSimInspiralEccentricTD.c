/*
 * Copyright (C) 2015 M. Haney, A. Gopakumar
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

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 * This structure contains the intrinsic parameters and angular velocity
 * co-efficient for the evolution equations for angular velocity, 
 * orbital eccentricity and the mean anomaly of the quasi-Keplerian orbit.
 * These are computed by XLALSimInspiralEccTaylorT4Setup routine.
 */

typedef struct
tagexpnCoeffsEccTaylorT4 {
   
   /* angular velocity coefficient*/
   REAL8 av;

   /* static (mass) parameters in orbital evolution equations */
   REAL8 mt,eta,ST;

   /* symmetric mass ratio, total mass, component masses*/
   REAL8 nu,m,m1,m2,mu;
}expnCoeffsEccTaylorT4;

/*
 * 
 */

typedef REAL8 (SimInspiralEvolutionEquations4)(
   REAL8 v,			/* post-Newtonian parameter */
   REAL8 et, 			/* orbital eccentricity */
   expnCoeffsEccTaylorT4 *ak	
);

/*
 * This strucuture contains pointers to the functions for calculating
 * the post-Newtonian accurate evolution equations at the desired order. They can be set by
 * XLALSimInspiralEccTaylorT4Setup by passing an appropriate PN order.
 */

typedef struct
tagexpnFuncEccTaylorT4
{
   SimInspiralEvolutionEquations4 *orbvel4;
   SimInspiralEvolutionEquations4 *orbecc4;
   SimInspiralEvolutionEquations4 *meananom4;
} expnFuncEccTaylorT4;

/*
 * Computes the rate of increase of the orbital velocity for a post-Newtonian
 * inspiral in an eccentric orbit.
 *
 * Implements Equations (3.12a), (3.14a), (3.15a), (B9a) and (B9b) of: 
 * Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * https://dcc.ligo.org/P1500148-v1
 * 
 * Note that above equation uses x = v^2 as the PN expansion parameter.
 */

static REAL8 
XLALSimInspiralOrbitalVelocityEvolution4_0PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdvdt;
	REAL8 POLdvdt0;
	
	OTS = sqrt(1. - et * et);

	CFdvdt = ak->eta * pow(v, 9.0) / (2. * ak->mt * ak->ST);
	POLdvdt0 = (192. + 584. * pow(et, 2.0) + 74 * pow(et, 4.0)) / (15. * pow(OTS, 7.0));

	ans = CFdvdt * POLdvdt0;

	return ans;
}

static REAL8 
XLALSimInspiralOrbitalVelocityEvolution4_2PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdvdt;
	REAL8 POLdvdt0;
	REAL8 POLdvdt2;

	OTS = sqrt(1. - et * et);

	CFdvdt = ak->eta * pow(v, 9.0) / (2. * ak->mt * ak->ST);
	POLdvdt0 = (192. + 584. * pow(et, 2.0) + 74 * pow(et, 4.0)) / (15. * pow(OTS, 7.0));
	POLdvdt2 = ((-11888. + pow(et, 2.0) * (87720. - 159600. * ak->eta) + pow(et, 4.0) * (171038. - 141708. * ak->eta)
		    + pow(et, 6.0) * (11717. - 8288. * ak->eta) - 14784. * ak->eta) * pow(v, 2.0)) / (420. * pow(OTS, 9.0));

	ans = CFdvdt * (POLdvdt0 + POLdvdt2);

	return ans;
}

static REAL8 
XLALSimInspiralOrbitalVelocityEvolution4_3PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdvdt;
	REAL8 POLdvdt0;
	REAL8 POLdvdt2;
	REAL8 FACdvdt3;
	REAL8 RFTdvdt3;

	OTS = sqrt(1. - et * et);

	CFdvdt = ak->eta * pow(v, 9.0) / (2. * ak->mt * ak->ST);
	POLdvdt0 = (192. + 584. * pow(et, 2.0) + 74 * pow(et, 4.0)) / (15. * pow(OTS, 7.0));
	POLdvdt2 = ((-11888. + pow(et, 2.0) * (87720. - 159600. * ak->eta) + pow(et, 4.0) * (171038. - 141708. * ak->eta)
		    + pow(et, 6.0) * (11717. - 8288. * ak->eta) - 14784. * ak->eta) * pow(v, 2.0)) / (420. * pow(OTS, 9.0));
	FACdvdt3 = (256./5.) * LAL_PI * pow(v, 3.0);
	RFTdvdt3 = (1. + 7.260831042 * pow(et, 2.0) + 5.844370473 * pow(et, 4.0) + 0.8452020270 * pow(et, 6.0) 
		    + 0.07580633432 * pow(et, 8.0) + 0.002034045037 * pow(et, 10.0)) / (1. - 4.900627291 * pow(et, 2.0) 
		    + 9.512155497 * pow(et, 4.0) - 9.051368575 * pow(et, 6.0) + 4.096465525 * pow(et, 8.0) 
		    - 0.5933309609 * pow(et, 10.0) - 0.05427399445 * pow(et, 12.0) - 0.009020225634 * pow(et, 14.0));
	
	ans = CFdvdt * (POLdvdt0 + POLdvdt2 + (FACdvdt3 * RFTdvdt3));
	
	return ans;
}

static REAL8 
XLALSimInspiralOrbitalVelocityEvolution4_4PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdvdt;
	REAL8 POLdvdt0;
	REAL8 POLdvdt2;
	REAL8 FACdvdt3;
	REAL8 RFTdvdt3;
	REAL8 POLdvdt4;
	
	OTS = sqrt(1. - et * et);

	CFdvdt = ak->eta * pow(v, 9.0) / (2. * ak->mt * ak->ST);
	POLdvdt0 = (192. + 584. * pow(et, 2.0) + 74 * pow(et, 4.0)) / (15. * pow(OTS, 7.0));
	POLdvdt2 = ((-11888. + pow(et, 2.0) * (87720. - 159600. * ak->eta) + pow(et, 4.0) * (171038. - 141708. * ak->eta)
		    + pow(et, 6.0) * (11717. - 8288. * ak->eta) - 14784. * ak->eta) * pow(v, 2.0)) / (420. * pow(OTS, 9.0));
	FACdvdt3 = (256./5.) * LAL_PI * pow(v, 3.0);
	RFTdvdt3 = (1. + 7.260831042 * pow(et, 2.0) + 5.844370473 * pow(et, 4.0) + 0.8452020270 * pow(et, 6.0) 
		    + 0.07580633432 * pow(et, 8.0) + 0.002034045037 * pow(et, 10.0)) / (1. - 4.900627291 * pow(et, 2.0) 
		    + 9.512155497 * pow(et, 4.0) - 9.051368575 * pow(et, 6.0) + 4.096465525 * pow(et, 8.0) 
		    - 0.5933309609 * pow(et, 10.0) - 0.05427399445 * pow(et, 12.0) - 0.009020225634 * pow(et, 14.0));
	POLdvdt4 = ((-360224. + 4514976. * ak->eta + 1903104. * pow(ak->eta, 2.0) 
		    + pow(et, 8.0) * (3523113. - 3259980. * ak->eta + 1964256. * pow(ak->eta, 2.0))
		    + pow(et, 2.0) * (-92846560. + 15464736. * ak->eta + 61282032. * pow(ak->eta, 2.0)) 
		    + pow(et, 6.0) * (83424402. - 123108426. * ak->eta + 64828848. * pow(ak->eta, 2.0))
		    + pow(et, 4.0) * (783768. - 207204264. * ak->eta + 166506060. * pow(ak->eta, 2.0)) 
		    - 3024. * (96. + 4268. * pow(et, 2.0) + 4386. * pow(et, 4.0)
		    + 175. * pow(et, 6.0))*(-5. + 2. * ak->eta) * OTS) * pow(v, 4.0)) / (45360. * pow(OTS, 11.0));

	ans = CFdvdt * (POLdvdt0 + POLdvdt2 + (FACdvdt3 * RFTdvdt3) + POLdvdt4);

	return ans;
}

/*
 * Computes the rate of increase of the orbital eccentricity for a post-Newtonian
 * inspiral in an eccentric orbit.
 *
 * Implements Equations (3.12b), (3.14a), (3.14b), (3.15b), (3.16), (B9c) and (B9d) of: 
 * Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * https://dcc.ligo.org/P1500148-v1
 * 
 * Note that above equation uses x = v^2 as the PN expansion parameter.
 */

static REAL8 
XLALSimInspiralOrbitalEccentricityEvolution4_0PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdedt;
	REAL8 POLdedt0;
	
	OTS = sqrt(1. - et * et);

	CFdedt = - et * ak->eta * pow(v, 8.0) / (ak->mt * ak->ST);
	POLdedt0 = (304. + 121. * pow(et, 2.0)) / (15. * pow(OTS, 5.0));

	ans = CFdedt * POLdedt0;

	return ans;
}

static REAL8 
XLALSimInspiralOrbitalEccentricityEvolution4_2PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdedt;
	REAL8 POLdedt0;
	REAL8 POLdedt2;
	
	OTS = sqrt(1. - et * et);    

	CFdedt = - et * ak->eta * pow(v, 8.0) / (ak->mt * ak->ST);
	POLdedt0 = (304. + 121. * pow(et, 2.0)) / (15. * pow(OTS, 5.0));
	POLdedt2 = - ((67608. + 228704. * ak->eta + pow(et, 4.0) * (-125361. + 93184. * ak->eta) 
		    + pow(et, 2.0) * (-718008. + 651252. * ak->eta)) * pow(v, 2.0)) / (2520. * pow(OTS, 7.0));

	ans = CFdedt * (POLdedt0 + POLdedt2);

	return ans;
}

static REAL8 
XLALSimInspiralOrbitalEccentricityEvolution4_3PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdedt;
	REAL8 POLdedt0;
	REAL8 POLdedt2;
	REAL8 FACdedt3;
	REAL8 RFTdedt3;
	
	OTS = sqrt(1. - et * et);

	CFdedt = - et * ak->eta * pow(v, 8.0) / (ak->mt * ak->ST);
	POLdedt0 = (304. + 121. * pow(et, 2.0)) / (15. * pow(OTS, 5.0));
	POLdedt2 = - ((67608. + 228704. * ak->eta + pow(et, 4.0) * (-125361. + 93184. * ak->eta) 
		    + pow(et, 2.0) * (-718008. + 651252. * ak->eta)) * pow(v, 2.0)) / (2520. * pow(OTS, 7.0));
	FACdedt3 = (394./3.) * LAL_PI * pow(v, 3.0);
	RFTdedt3 = 0.1949238579 * OTS * (OTS * (1. + 7.260831042 * pow(et, 2.0) + 5.844370473 * pow(et, 4.0) 
		    + 0.8452020270 * pow(et, 6.0) + 0.07580633432 * pow(et, 8.0) 
		    + 0.002034045037 * pow(et, 10.0)) / (1. - 4.900627291 * pow(et, 2.0) + 9.512155497 * pow(et, 4.0) 
		    - 9.051368575 * pow(et, 6.0) + 4.096465525 * pow(et, 8.0) - 0.5933309609 * pow(et, 10.0) 
		    - 0.05427399445 * pow(et, 12.0) - 0.009020225634 * pow(et, 14.0)) - (1. + 1.893242666 * pow(et, 2.0) 
		    - 2.708117333 * pow(et, 4.0) + 0.6192474531 * pow(et, 6.0) + 0.0500847462 * pow(et, 8.0) 
		    - 0.01059040781 * pow(et, 10.0)) / (1. - 4.638007334 * pow(et, 2.0) + 8.716680569 * pow(et, 4.0) 
		    - 8.451197591 * pow(et, 6.0) + 4.435922348 * pow(et, 8.0) - 1.199023304 * pow(et, 10.0) 
		    + 0.1398678608 * pow(et, 12.0) - 0.004254544193 * pow(et, 14.0))) / pow(et, 2.0);

	ans = CFdedt * (POLdedt0 + POLdedt2 + (FACdedt3 * RFTdedt3));

	return ans;
}

static REAL8 
XLALSimInspiralOrbitalEccentricityEvolution4_4PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdedt;
	REAL8 POLdedt0;
	REAL8 POLdedt2;
	REAL8 FACdedt3;
	REAL8 RFTdedt3;
	REAL8 POLdedt4;
	
	OTS = sqrt(1. - et * et);

	CFdedt = - et * ak->eta * pow(v, 8.0) / (ak->mt * ak->ST);
	POLdedt0 = (304. + 121. * pow(et, 2.0)) / (15. * pow(OTS, 5.0));
	POLdedt2 = - ((67608. + 228704. * ak->eta + pow(et, 4.0) * (-125361. + 93184. * ak->eta) 
		    + pow(et, 2.0) * (-718008. + 651252. * ak->eta)) * pow(v, 2.0)) / (2520. * pow(OTS, 7.0));
	FACdedt3 = (394./3.) * LAL_PI * pow(v, 3.0);
	RFTdedt3 = 0.1949238579 * OTS * (OTS * (1. + 7.260831042 * pow(et, 2.0) + 5.844370473 * pow(et, 4.0) 
		    + 0.8452020270 * pow(et, 6.0) + 0.07580633432 * pow(et, 8.0) 
		    + 0.002034045037 * pow(et, 10.0)) / (1. - 4.900627291 * pow(et, 2.0) + 9.512155497 * pow(et, 4.0) 
		    - 9.051368575 * pow(et, 6.0) + 4.096465525 * pow(et, 8.0) - 0.5933309609 * pow(et, 10.0) 
		    - 0.05427399445 * pow(et, 12.0) - 0.009020225634 * pow(et, 14.0)) - (1. + 1.893242666 * pow(et, 2.0) 
		    - 2.708117333 * pow(et, 4.0) + 0.6192474531 * pow(et, 6.0) + 0.0500847462 * pow(et, 8.0) 
		    - 0.01059040781 * pow(et, 10.0)) / (1. - 4.638007334 * pow(et, 2.0) + 8.716680569 * pow(et, 4.0) 
		    - 8.451197591 * pow(et, 6.0) + 4.435922348 * pow(et, 8.0) - 1.199023304 * pow(et, 10.0) 
		    + 0.1398678608 * pow(et, 12.0) - 0.004254544193 * pow(et, 14.0))) / pow(et, 2.0);
	POLdedt4 = ((-15238352. + 12823920. * ak->eta + 4548096. * pow(ak->eta, 2.0) 
		    + pow(et, 6.0) * (3786543. - 4344852. * ak->eta + 2758560. * pow(ak->eta, 2.0)) 
		    + pow(et, 4.0) * (46566110. - 78343602. * ak->eta + 42810096 * pow(ak->eta, 2.0)) 
		    + pow(et, 2.0) * (-37367868 - 41949252 * ak->eta + 48711348 * pow(ak->eta, 2.0)) - 1008. * (2672. 
		    + 6963. * pow(et, 2.0) + 565. * pow(et, 4.0)) * (-5. + 2. * ak->eta) * OTS) * pow(v, 4.0)) / (30240. * pow(OTS, 9.0));

	ans = CFdedt * (POLdedt0 + POLdedt2 + (FACdedt3 * RFTdedt3) + POLdedt4);

	return ans;
}

/*
 * Computes the rate of increase of the mean anomaly for a post-Newtonian
 * inspiral in an eccentric orbit.
 *
 * Implements Equations (3.17a), (B9e) and (B9f) of: 
 * Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * https://dcc.ligo.org/P1500148-v1
 * 
 * Note that above equation uses x = v^2 as the PN expansion parameter.
 */

static REAL8 
XLALSimInspiralMeanAnomalyEvolution4_0PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 UNUSED et,	/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;

	REAL8 CFdldt;
	REAL8 POLdldt0;
	
	CFdldt = pow(v, 3.0) / (ak->mt * ak->ST);
	POLdldt0 = 1.;

	ans = CFdldt * POLdldt0;

	return ans;
}

static REAL8 
XLALSimInspiralMeanAnomalyEvolution4_2PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdldt;
	REAL8 POLdldt0;
	REAL8 POLdldt2;

	OTS = sqrt(1. - et * et);

	CFdldt = pow(v, 3.0) / (ak->mt * ak->ST);
	POLdldt0 = 1.;
	POLdldt2 = - (3. * pow(v, 2.0)) / pow(OTS, 2.0);

	ans = CFdldt * (POLdldt0 + POLdldt2);

	return ans;
}

static REAL8 
XLALSimInspiralMeanAnomalyEvolution4_4PN(
	REAL8 v,		/* post-Newtonian parameter */
	REAL8 et,		/* orbital eccentricity */
	expnCoeffsEccTaylorT4 *ak	/* PN co-efficients and intrinsic parameters */
	)
{
	REAL8 ans;
	REAL8 OTS;

	REAL8 CFdldt;
	REAL8 POLdldt0;
	REAL8 POLdldt2;
	REAL8 POLdldt4;
	
	OTS = sqrt(1. - et * et);

	CFdldt = pow(v, 3.0) / (ak->mt * ak->ST);
	POLdldt0 = 1.;
	POLdldt2 = - (3. * pow(v, 2.0)) / pow(OTS, 2.0);
	POLdldt4 = ((-18. + 28. * ak->eta + pow(et, 2.0)*(-51. + 26. * ak->eta)) * pow(v, 4.0)) / (4. * pow(OTS, 4.0));

	ans = CFdldt * (POLdldt0 + POLdldt2 + POLdldt4);

	return ans;
}

/**
 * 
 */

typedef struct
{
    REAL8 (*funcv)(REAL8 v, REAL8 et, expnCoeffsEccTaylorT4 *ak);
    REAL8 (*funcet)(REAL8 v, REAL8 et, expnCoeffsEccTaylorT4 *ak);
    REAL8 (*funcl)(REAL8 v, REAL8 et, expnCoeffsEccTaylorT4 *ak);
    expnCoeffsEccTaylorT4 ak;
}XLALSimInspiralEccTaylorT4PNEvolveOrbitParams;

/*
 * This function is used in the call to the GSL integrator.
 * 
 * Implements Equation (3.17b) of: Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 */

static int 
XLALSimInspiralEccTaylorT4PNEvolveOrbitIntegrand(double UNUSED t, const double y[], double ydot[], void *params)
{
	XLALSimInspiralEccTaylorT4PNEvolveOrbitParams* p = (XLALSimInspiralEccTaylorT4PNEvolveOrbitParams*)params;
	ydot[0] = p->funcv(y[0],y[1],&p->ak);
	ydot[1] = p->funcet(y[0],y[1],&p->ak);
	ydot[2] = p->funcl(y[0],y[1],&p->ak);
	ydot[3] = y[0]*y[0]*y[0]*p->ak.av;		
	t = 0.0;
	return GSL_SUCCESS;
}

/*
 * Set up the expnCoeffsEccTaylorT4 and expnFuncEccTaylorT4 structures for
 * generating the eccentric TaylorT4 waveform and select the post-Newtonian
 * functions corresponding to the desired order.
 *
 * Inputs given in SI units.
 */

static int 
XLALSimInspiralEccTaylorT4Setup(
    expnCoeffsEccTaylorT4 *ak,         /* coefficients for TaylorT4 evolution [modified] */
    expnFuncEccTaylorT4 *f,            /* functions for TaylorT4 evolution [modified] */
    REAL8 m1,                       /* mass of companion 1 */
    REAL8 m2,                       /* mass of companion 2 */
    int O                           /* twice post-Newtonian order */
)
{
    ak->m1 = m1;
    ak->m2 = m2;
    ak->m = ak->m1 + ak->m2;
    ak->mu = m1 * m2 / ak->m;
    ak->nu = ak->mu / ak->m;
    
    ak->eta = m1 * m2 / pow(ak->m, 2.0);
    ak->mt = ak->m / LAL_MSUN_SI;
    ak->ST = LAL_MTSUN_SI;

    /* angular velocity co-efficient */
    ak->av = 1./(ak->ST * ak->mt);

    switch (O)
    {
	case 0:
            f->orbvel4 = &XLALSimInspiralOrbitalVelocityEvolution4_0PN;
	    f->orbecc4 = &XLALSimInspiralOrbitalEccentricityEvolution4_0PN;
	    f->meananom4 = &XLALSimInspiralMeanAnomalyEvolution4_0PN;
            break;
	case 1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
	case 2:
            f->orbvel4 = &XLALSimInspiralOrbitalVelocityEvolution4_2PN;
	    f->orbecc4 = &XLALSimInspiralOrbitalEccentricityEvolution4_2PN;
	    f->meananom4 = &XLALSimInspiralMeanAnomalyEvolution4_2PN;
            break;
	case 3:
            f->orbvel4 = &XLALSimInspiralOrbitalVelocityEvolution4_3PN;
	    f->orbecc4 = &XLALSimInspiralOrbitalEccentricityEvolution4_3PN;
	    f->meananom4 = &XLALSimInspiralMeanAnomalyEvolution4_2PN;
            break;
	case -1: // Highest implemented PN order - move if higher terms added!
        case 4:
            f->orbvel4 = &XLALSimInspiralOrbitalVelocityEvolution4_4PN;
	    f->orbecc4 = &XLALSimInspiralOrbitalEccentricityEvolution4_4PN;
	    f->meananom4 = &XLALSimInspiralMeanAnomalyEvolution4_4PN;
            break;
	case 5:
	    XLALPrintError("XLAL Error - %s: PN approximant not yet implemented for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
	case 6:
	    XLALPrintError("XLAL Error - %s: PN approximant not yet implemented for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
	case 7:
	    XLALPrintError("XLAL Error - %s: PN approximant not yet implemented for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
	case 8:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
	default:
            XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
    }
  
  return 0;
}

/*
 * Given time series of the mean anomaly the orbital eccentricity of a quasi-Keplerian orbit,
 * solves Kepler equation using a modified version of Mikkola's method.
 * 
 * See:
 * Manuel Tessmer, and Achamveedu Gopakumar,
 * \"Accurate and efficient gravitational waveforms for certain galactic
 * compact binaries\", Mon. Not. R. Astron. Soc. 374, 721–728 (2007);
 * arXiv:gr-qc/0610139
 * 
 * and references therein.
 */

static int XLALSimInspiralKeplerEquationSolver(
		REAL8TimeSeries **V,		/* post-Newtonian parameter */
		REAL8TimeSeries **ECC,		/* orbital eccentricity */
		REAL8TimeSeries *TIME,		/* mean anomaly of quasi-Keplerian orbit */
		REAL8TimeSeries *U		/* eccentric anomaly of quasi-Keplerian orbit [returned] */
	   ){
  
  UINT4 len, k = 0;
  len = (*V)->data->length;
  
  REAL8 l, et;
  
  REAL8 ND, SIGN, SIGN2, alpha, beta1, zplus, zminus, SIGN3, z, s, w, E0, u;
  REAL8 f, f1, f2, f3, f4, u1, u2, u3, u4, xi, sol;
  
  for(k = 0; k < len; k++){
    
	l = TIME->data->data[k];	et = (*ECC)->data->data[k];
	
	/* mapping of l into interval -Pi<=l<=Pi to make auxiliary variable s as small as possible */
	if(l > 0.0 ||l == 0.0){ ND = floor(l/(2.*LAL_PI));
				      SIGN = 1.; }
	else{ ND = ceil(l/(2.*LAL_PI)); 
	      SIGN = -1.; }
	
	l = SIGN*l;
	
	if(l > (2.*LAL_PI)){ l = (l - floor(l/(2.*LAL_PI))*2.*LAL_PI); }
	
	if(l < LAL_PI){ SIGN2 = 1.; }
	else{ SIGN2 = -1.; }

	if(SIGN2 == 1.){ ; }
	else{ l = (2.*LAL_PI - l); }
	
	/* root-finding method to solve KE as cubic equation in terms of s=sin(u/3) */
	alpha  = (1. - et)/(4.*et + 1./2.);
	beta1  = (l/2.)/(4.*et + 1./2.);
	zplus  = pow(beta1 + sqrt(pow(alpha, 3.0) + pow(beta1, 2.0)), (1./3.));
	zminus = pow(beta1 - sqrt(pow(alpha, 3.0) + pow(beta1, 2.0)), (1./3.));
	
	if(beta1 > 0.0 ||beta1 == 0.0){ SIGN3 = 1.; }
	else{ SIGN3 = -1.; }
	
	if(SIGN3 > 0.0){ z = zplus; }
	else{ z = zminus; }
	
	z = zplus;
	s = (z - alpha/z);
	
	/* negate largest error occuring at l=Pi by correction term */
	w = (s - (0.078*pow(s, 5.0))/(1. + et));
	
	/* initial guess for u */
	E0 = (l + et*(3.*w - 4.*pow(w, 3.0)));
	u = E0 ;

	/* define f(u) and its first four derivatives with respect to u */
	f = u - et*sin(u) - l;
	f1 = 1. - et*cos(u);
	f2 = et*sin(u);
	f3 = et*cos(u);
	f4 = -et*sin(u);
	
	/* apply a fourth-order extension of Newton’s method to improve on initial guess */
	u1 = -f/f1;
	u2 = -f/(f1 + (1./2.)*f2*u1);
	u3 = -f/(f1 + (1./2.)*f2*u2 + (1./6.)*f3*(pow(u2, 2.0)));
	u4 = -f/(f1 + (1./2.)*f2*u3 + (1./6.)*f3*(pow(u3, 2.0)) + (1./24.)*f4*(pow(u3, 3.0)));
	xi = (E0 + u4);
	
	if(SIGN2 > 0.0){ sol = xi; }
	else{ sol = (2.*LAL_PI - xi); }
	
	if(SIGN < 0.0){ sol = -sol; }
	
	u = (sol + ND*2.*LAL_PI);
	
	U->data->data[k] = u;
  }
  
  return XLAL_SUCCESS; 

}

/*
 * Given time series for orbital velocity, orbital eccentricity as well as mean
 * and eccentric anomaly of the quasi-Keplerian orbit, computes the PN-contributions
 * to the 2PN-accurate Kepler equation.
 * 
 * Implements Equation (3.11) of: Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * 
 * Note that above equation uses x = v^2 as the PN expansion parameter.
 */

static int XLALSimInspiralKeplerEquationLHS_PN(
		expnCoeffsEccTaylorT4 *ak,
		REAL8TimeSeries **V,		/* post-Newtonian parameter */
		REAL8TimeSeries *U,		/* eccentric anomaly of quasi-Keplerian orbit */
		REAL8TimeSeries **ECC,		/* orbital eccentricity */
		REAL8TimeSeries **TIME,		/* mean anomaly l of quasi-Keplerian orbit */
		REAL8TimeSeries *L2		/* LHS of 2PN-accurate Kepler equation, i.e., l - l_2PN [returned] */
	    ){
  
  UINT4 len, k = 0;
  len = (*V)->data->length;
  
  REAL8 u, et, v, xp;
  REAL8 dt;
  REAL8 OTS, betaN, vmuN;
  REAL8 L2PN;
  
  for(k = 0; k < len; k++){
    
    /* Abbreviated names in lower case for time series at this sample */
    u = U->data->data[k];	et = (*ECC)->data->data[k];
    v = (*V)->data->data[k];	xp = (*TIME)->data->data[k];
    
    dt    = 1. - et*cos(u);
    
    OTS   = sqrt(1. - et*et);
    betaN = (1. - OTS)/et;
    vmuN = 2*atan(betaN*sin(u)/(1. - betaN*cos(u)));
    
    L2PN = pow(v, 4.0)*(dt*(60. - 24.*ak->eta)*vmuN 
	      + et*(15.*ak->eta - pow(ak->eta, 2.0))*OTS*sin(u))/(8.*dt*OTS);
	      
    L2->data->data[k] = xp - L2PN;
    
  }
	      
 return XLAL_SUCCESS;
  
}

/*
 * Given time series for orbital velocity, orbital eccentricity and mean 
 * anomaly of the quasi-Keplerian orbit, computes the periodic contribution W 
 * to the phase split, i.e., phi = lambda + W.
 * 
 * Implements Equations (3.8) - (3.10) and (B4c), (B4d), (B7a), (B7b) of: 
 * Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * 
 * Note that above equation uses x = v^2 as the PN expansion parameter.
 */

static int XLALSimInspiralPhaseContributionW(
		expnCoeffsEccTaylorT4 *ak,
		REAL8TimeSeries **V,		/* post-Newtonian parameter */
		REAL8TimeSeries **ECC,		/* orbital eccentricity */
		REAL8TimeSeries **U,		/* eccentric anomaly of quasi-Keplerian orbit */
		REAL8TimeSeries *W		/* periodic contribution to phase split [returned] */
	   ){
    
  UINT4 len, k = 0;
  len = (*V)->data->length;
  
  REAL8 v, et, u, w;
  REAL8 dt, OTS;
  REAL8 beta1PN, beta2PN;
  REAL8 vmu1PN, vmu2PN;
  
  for(k = 0; k < len; k++){
    
    v = (*V)->data->data[k];	et = (*ECC)->data->data[k];
    u = (*U)->data->data[k];
    
    dt    = 1. - et * cos(u);
    OTS   = sqrt(1. - et * et);
    
    beta1PN = (1. - OTS) / et
	      //=========================================================================================================
	      + ((-4. + pow(et, 2.0) * (8. - 2. * ak->eta) + ak->eta + (4. - ak->eta) * OTS) * pow(v, 2.0)) / (et * OTS);
	      //=========================================================================================================
    beta2PN = beta1PN
	      //=========================================================================================================
	      + (-528. - 220. * ak->eta + 4. * pow(ak->eta, 2.0) + pow(et, 4.0) * (-3840. + 2086. * ak->eta 
	      - 178. * pow(ak->eta, 2.0)) + pow(et, 2.0) * (5232. - 1659. * ak->eta + 177. * pow(ak->eta, 2.0)) 
	      + (528. + 220. * ak->eta - 4. * pow(ak->eta, 2.0) + pow(et, 2.0) * (288. + 83. * ak->eta 
	      - 41. * pow(ak->eta, 2.0))) * OTS) * pow(v, 4.0) / (96. * et * pow(OTS, 3.0));
	      //=========================================================================================================
	    
    vmu1PN = 2. * atan(beta1PN * sin(u) / (1. - beta1PN * cos(u)));
    vmu2PN = 2. * atan(beta2PN * sin(u) / (1. - beta2PN * cos(u)));
    
    w	= et * sin(u) + vmu2PN
	  //=========================================================================================================
	  + 3. * (et * sin(u) + vmu1PN) * pow(v, 2.0) / pow(OTS, 2.0)
	  //=========================================================================================================
	  + et * sin(u) * (4. * pow(dt, 2.0) * (dt * (108. + pow(et, 2.0) * (102. - 52. * ak->eta) - 56. * ak->eta)
	  - 15. * ak->eta + pow(ak->eta, 2.0) + pow(et, 2.0) * (30. * ak->eta - 2. * pow(ak->eta, 2.0))
	  + pow(et, 4.0) * (-15. * ak->eta + pow(ak->eta, 2.0))) + (4. * ak->eta - 12. * pow(ak->eta, 2.0)
	  + dt * (8. + pow(et, 2.0) * (-8. - 144. * ak->eta) + 144 * ak->eta) + pow(et, 4.0) * (4. * ak->eta 
	  - 12. * pow(ak->eta, 2.0)) + pow(et, 2.0) * (-8. * ak->eta + 24. * pow(ak->eta, 2.0)) 
	  + pow(dt, 2.0) * (-8. - 148. * ak->eta + 12. * pow(ak->eta, 2.0) + pow(et, 2.0) * (- ak->eta 
	  + 3. * pow(ak->eta, 2.0)))) * OTS) * pow(v, 4.0) / (32. * pow(dt, 3.0) * pow(OTS, 4.0));
	  //=========================================================================================================
	      
    W->data->data[k] = w;
    
  }
	      
 return XLAL_SUCCESS;
  
}

/**
 * Evolves a post-Newtonian orbit using the eccentric Taylor T4 method.
 *
 * See:
 * Sashwat Tanay, Maria Haney, and Achamveedu Gopakumar,
 * \"Frequency and time domain inspiral waveforms from comparable 
 * mass compact binaries in eccentric orbits\", (2015);
 * arXiv:TBD
 * https://dcc.ligo.org/P1500148-v1
 */

int XLALSimInspiralEccentricTDPNEvolveOrbit(
		REAL8TimeSeries **v,            /**< post-Newtonian parameter [returned] */
		REAL8TimeSeries **et,		/**< orbital eccentricity [returned] */
		REAL8TimeSeries **l,		/**< mean anomaly of quasi-Keplerian orbit [returned] */
		REAL8TimeSeries **lambda,       /**< secular contribution to orbital phase [returned] */
		REAL8TimeSeries **u,		/**< eccentric anomaly of quasi-Keplerian orbit [returned] */
		REAL8TimeSeries **phi,		/**< orbital phase [returned] */
		REAL8 phiRef,                   /**< reference orbital phase (rad) */
		REAL8 deltaT,                   /**< sampling interval (s) */
		REAL8 m1,                       /**< mass of companion 1 (kg) */
		REAL8 m2,                       /**< mass of companion 2 (kg) */
		REAL8 f_min,                    /**< start frequency (Hz) */
		REAL8 fRef,                     /**< reference frequency (Hz) */
		REAL8 e_min,			/**< initial orbital eccentricity at f_min */
		int O                           /**< twice post-Newtonian order */
		)
{
	const UINT4 blocklen = 1024;
	const REAL8 visco = 1./sqrt(6.);
	REAL8 VRef = 0.;
	XLALSimInspiralEccTaylorT4PNEvolveOrbitParams params;
	expnFuncEccTaylorT4 expnfunc;
	expnCoeffsEccTaylorT4 ak;

	if(XLALSimInspiralEccTaylorT4Setup(&ak,&expnfunc,m1,m2,O))
		XLAL_ERROR(XLAL_EFUNC);

	params.funcv=expnfunc.orbvel4;
	params.funcet=expnfunc.orbecc4;
	params.funcl=expnfunc.meananom4;
	params.ak=ak;

	UINT4 j, len, idxRef = 0;
	LIGOTimeGPS tc = LIGOTIMEGPSZERO;
	double y[4];
	double yerr[4];
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
	gsl_odeiv_step *s;
	gsl_odeiv_system sys;

	/* setup ode system */
	sys.function = XLALSimInspiralEccTaylorT4PNEvolveOrbitIntegrand;
	sys.jacobian = NULL;
	sys.dimension = 4;
	sys.params = &params;

	/* allocate memory */
	*v = XLALCreateREAL8TimeSeries( "ORBITAL_VELOCITY_PARAMETER", &tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	*et = XLALCreateREAL8TimeSeries( "ORBITAL_ECCENTRICITY", &tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	*l = XLALCreateREAL8TimeSeries( "ORBITAL_MEAN_ANOMALY", &tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	*lambda = XLALCreateREAL8TimeSeries( "SECULAR_ORBITAL_PHASE", &tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	if ( !v || !et || !l || !lambda )
		XLAL_ERROR(XLAL_EFUNC);

	y[0] = (*v)->data->data[0] = cbrt(LAL_PI*LAL_G_SI*ak.m*f_min)/LAL_C_SI;
	y[1] = (*et)->data->data[0] = e_min;
	y[2] = (*l)->data->data[0] = 0.;
	y[3] = (*lambda)->data->data[0] = 0.;

	j = 0;

	s = gsl_odeiv_step_alloc(T, 4);
	while (1) {
		++j;
		gsl_odeiv_step_apply(s, j*deltaT, deltaT, y, yerr, NULL, NULL, &sys);
		/* ISCO termination condition for quadrupole, 1PN, 2.5PN */
		if ( y[0] > visco ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", __func__);
			break;
		}
		if ( j >= (*v)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*v, 0, (*v)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*et, 0, (*et)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*l, 0, (*l)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*lambda, 0, (*lambda)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
		}
		(*v)->data->data[j] = y[0];
		(*et)->data->data[j] = y[1];
		(*l)->data->data[j] = y[2];
		(*lambda)->data->data[j] = y[3];
	}
	gsl_odeiv_step_free(s);	

	/* make the correct length */
	if ( ! XLALResizeREAL8TimeSeries(*v, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*et, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*l, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*lambda, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);

	/* adjust to correct time */
	XLALGPSAdd(&(*v)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*et)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*l)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*lambda)->epoch, -1.0*j*deltaT);

	UINT4 LEN, k = 0;
	LEN = (*v)->data->length;
	
	/* 2-step solution of 1PN- and 2PN-accurate Kepler equation to obtain u(l) */
	/* 1PN KE: l1 = u1 - et*sin(u1) */
	/* 2PN KE: l2 = l1 - l_2PN = u2 - et*sin(u2) */
	
	REAL8TimeSeries *l1;
	REAL8TimeSeries *u1;
	REAL8TimeSeries *l2;
	REAL8TimeSeries *u2;
	
	/* allocate memory for auxiliary time series */
	l1 = XLALCreateREAL8TimeSeries( "MEAN_ANOMALY_1PN", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, (*v)->data->length );	
	memset(l1->data->data, 0, l1->data->length
			* sizeof(*l1->data->data));
	u1 = XLALCreateREAL8TimeSeries( "ECCENTRIC_ANOMALY_1PN", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, (*v)->data->length );	
	memset(u1->data->data, 0, u1->data->length
			* sizeof(*u1->data->data));
	l2 = XLALCreateREAL8TimeSeries( "2PN_accurate_l", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, (*v)->data->length );
	memset(l2->data->data, 0, l2->data->length
			* sizeof(*l2->data->data));
	u2 = XLALCreateREAL8TimeSeries( "ECCENTRIC_ANOMALY_2PN", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, (*v)->data->length );	
	memset(u2->data->data, 0, u2->data->length
			* sizeof(*u2->data->data));
	
	
	for(k = 0; k < LEN; k++){
	  l1->data->data[k] = (*l)->data->data[k];
	}
	if(XLALSimInspiralKeplerEquationSolver(v,et,l1,u1))
		XLAL_ERROR(XLAL_EFUNC);
	if(XLALSimInspiralKeplerEquationLHS_PN(&ak,v,u1,et,l,l2))
		XLAL_ERROR(XLAL_EFUNC);
	if(XLALSimInspiralKeplerEquationSolver(v,et,l2,u2))
		XLAL_ERROR(XLAL_EFUNC);
	
	/* allocate memory for u [returned] */
	*u = XLALCreateREAL8TimeSeries( "ECCENTRIC_ANOMALY", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, LEN );
	
	for(k = 0; k < LEN; k++){
	  (*u)->data->data[k] = u2->data->data[k];
	}
	
	XLALDestroyREAL8TimeSeries(l1);
	XLALDestroyREAL8TimeSeries(u1);
	XLALDestroyREAL8TimeSeries(l2);
	XLALDestroyREAL8TimeSeries(u2);
	
	/* compute full orbital phase variable phi entering h+/x */
	
	/* allocate memory for periodic contribution to phase split phi = lambda + W */
	REAL8TimeSeries *w;
	w = XLALCreateREAL8TimeSeries( "PERIODIC_PHASE_CONTRIBUTION", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, (*v)->data->length );	
	memset(w->data->data, 0, w->data->length
			* sizeof(*w->data->data));
	
	if(XLALSimInspiralPhaseContributionW(&ak,v,et,u,w))
	  XLAL_ERROR(XLAL_EFUNC);
	
	*phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", &(*v)->epoch, 0., (*v)->deltaT, &lalDimensionlessUnit, LEN );
	
	for(k = 0; k < LEN; k++){
	  (*phi)->data->data[k] = (*lambda)->data->data[k] + w->data->data[k];
	}
	
	XLALDestroyREAL8TimeSeries(w);
	
	/* Do a constant phase shift to get desired value of phiRef */
	len = (*phi)->data->length;
	/* For fRef==0, phiRef is phase of last sample */
	if( fRef == 0. )
		phiRef -= (*phi)->data->data[len-1];
	/* For fRef==fmin, phiRef is phase of first sample */
	else if( fRef == f_min )
		phiRef -= (*phi)->data->data[0];
	/* phiRef is phase when f==fRef */
	else
	{
		VRef = pow(LAL_PI * LAL_G_SI*(m1+m2) * fRef, 1./3.) / LAL_C_SI;
		j = 0;
		do {
			idxRef = j;
			j++;
		} while ((*v)->data->data[j] <= VRef);
		phiRef -= (*phi)->data->data[idxRef];
	}
	for (j = 0; j < len; ++j)
		(*phi)->data->data[j] += phiRef;

	return (int)(*v)->data->length;
		
}

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different PN orders
 * for phasing calculation vs. amplitude calculations.
 */

int XLALSimInspiralEccentricTDPNGenerator(
		REAL8TimeSeries **hplus,        /**< +-polarization waveform */
		REAL8TimeSeries **hcross,       /**< x-polarization waveform */
		REAL8 phiRef,                   /**< reference orbital phase (rad) */
		REAL8 deltaT,                   /**< sampling interval (s) */
		REAL8 m1,                       /**< mass of companion 1 (kg) */
		REAL8 m2,                       /**< mass of companion 2 (kg) */
		REAL8 f_min,                    /**< start frequency (Hz) */
		REAL8 fRef,                     /**< reference frequency (Hz) */
		REAL8 r,                        /**< distance of source (m) */
		REAL8 i,                        /**< inclination of source (rad) */
		REAL8 e_min,			/**< initial orbital eccentricity at f_min */
		int amplitudeO,                 /**< twice post-Newtonian amplitude order */
		int phaseO                      /**< twice post-Newtonian phase order */
		)
{

	/* The Schwarzschild ISCO frequency - for sanity checking fRef */
	REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

	/* Sanity check fRef value */
	if( fRef < 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n",
				__func__, fRef);
		XLAL_ERROR(XLAL_EINVAL);
	}
	if( fRef != 0. && fRef < f_min )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
				__func__, fRef, f_min);
		XLAL_ERROR(XLAL_EINVAL);
	}
	if( fRef >= fISCO )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
				__func__, fRef, fISCO);
		XLAL_ERROR(XLAL_EINVAL);
	}


	REAL8TimeSeries *v;
	REAL8TimeSeries *et;
	REAL8TimeSeries *l;
	REAL8TimeSeries *lambda;
	
	REAL8TimeSeries *u;
	REAL8TimeSeries *phi;
	
	int n;	
	n = XLALSimInspiralEccentricTDPNEvolveOrbit(&v, &et, &l, &lambda, &u, &phi, phiRef, deltaT,
			m1, m2, f_min, fRef, e_min, phaseO);
	if ( n < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	
	int status;
	status = XLALSimInspiralPNPolarizationWaveformsEccentric(hplus, hcross, v, et, u, phi,
			m1, m2, r, i, amplitudeO, phaseO);
	
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(u);

	XLALDestroyREAL8TimeSeries(lambda);
	XLALDestroyREAL8TimeSeries(l);
	XLALDestroyREAL8TimeSeries(et);
	XLALDestroyREAL8TimeSeries(v);
	
	if ( status < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	
	return n;
}

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 * 
 * Note that at present phasing calculations for eccentric waveforms are implemented 
 * up to 2PN order, while amplitude calculations are accurate up to Newtonian 
 * (quadrupolar) order.
 * Higher-order contributions to phasing and amplitude calculations will be implemented
 * in the near future.
 */

int XLALSimInspiralEccentricTDPN(
		REAL8TimeSeries **hplus,        /**< +-polarization waveform */
		REAL8TimeSeries **hcross,       /**< x-polarization waveform */
		REAL8 phiRef,                   /**< reference orbital phase (rad) */
		REAL8 deltaT,                   /**< sampling interval (Hz) */
		REAL8 m1,                       /**< mass of companion 1 (kg) */
		REAL8 m2,                       /**< mass of companion 2 (kg) */
		REAL8 f_min,                    /**< start frequency (Hz) */
		REAL8 fRef,                     /**< reference frequency (Hz) */
		REAL8 r,                        /**< distance of source (m) */
		REAL8 i,                        /**< inclination of source (rad) */
		REAL8 e_min,			/**< initial orbital eccentricity at f_min */
		int O                           /**< twice post-Newtonian order */
		)
{
	return XLALSimInspiralEccentricTDPNGenerator(hplus, hcross, phiRef,
			deltaT, m1, m2, f_min, fRef, r, i, e_min, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 * 
 * Note that amplitudes of eccentric waveforms are at present Newtonian order by default.
 * Amplitude corrections will be implemented in the near future.
 */
int XLALSimInspiralEccentricTDPNRestricted(
		REAL8TimeSeries **hplus,        /**< +-polarization waveform */
		REAL8TimeSeries **hcross,       /**< x-polarization waveform */
		REAL8 phiRef,                   /**< reference orbital phase (rad) */
		REAL8 deltaT,                   /**< sampling interval (s) */
		REAL8 m1,                       /**< mass of companion 1 (kg) */
		REAL8 m2,                       /**< mass of companion 2 (kg) */
		REAL8 f_min,                    /**< start frequency (Hz) */
		REAL8 fRef,                     /**< reference frequency (Hz) */
		REAL8 r,                        /**< distance of source (m) */
		REAL8 i,                        /**< inclination of source (rad) */
		REAL8 e_min,			/**< initial orbital eccentricity at f_min */
		int O                           /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	return XLALSimInspiralEccentricTDPNGenerator(hplus, hcross, phiRef, 
			deltaT, m1, m2, f_min, fRef, r, i, e_min, 0, O);
}


#if 0
#include <lal/PrintFTSeries.h>
#include <lal/PrintFTSeries.h>
int main(void)
{
	LIGOTimeGPS tc = { 888888888, 222222222 };
	REAL8 phic = 1.0;
	REAL8 deltaT = 1.0/16384.0;
	REAL8 m1 = 1.4*LAL_MSUN_SI;
	REAL8 m2 = 1.4*LAL_MSUN_SI;
	REAL8 r = 1e6*LAL_PC_SI;
	REAL8 i = 0.5*LAL_PI;
	REAL8 f_min = 100.0;
	REAL8 fRef = 0.;
	int O = -1;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	XLALSimInspiralEccentricTDPNRestricted(&hplus, &hcross, &tc, phic, deltaT, m1, m2, f_min, fRef, r, i, e_min, O);
	LALDPrintTimeSeries(hplus, "hp.dat");
	LALDPrintTimeSeries(hcross, "hc.dat");
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	LALCheckMemoryLeaks();
	return 0;
}
#endif
