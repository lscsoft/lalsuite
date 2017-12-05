/*
 *  Copyright (C) 2016 Riccardo Sturani
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

#include <stdlib.h>
#include <math.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

#define EPSILON 1.e-11

static int compare(
    REAL8 val1,
    REAL8 val2)
{
    if (fabs(val1 - val2) > EPSILON)
        return 1;
    else
        return 0;
}

int main (int argc, char **argv)
{
    /* Ignore unused parameters. */
    (void)argc;
    (void)argv;

    const UINT4 Ntest=8;
    REAL8 idxr;
    int ret = 0;
    int errCode=0;
    REAL8 incl, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z;
    REAL8 thetaJN, phiJL, theta1, theta2, phi12, chi1, chi2, phiRef;
    REAL8 inclIn,s1xIn, s1yIn, s1zIn, s2xIn, s2yIn, s2zIn;
    REAL8 m1, m2, fRef;
    const REAL8 m1min=1.4;
    const REAL8 m1max=100.;
    const REAL8 m2min=1.4;
    const REAL8 m2max=50.;

    LALSimInspiralWaveformFlags *waveFlag=XLALSimInspiralCreateWaveformFlags();

    for (UINT4 idx=0;idx<=Ntest;idx++) {
      idxr=((REAL8)idx)/((REAL8) Ntest);
      thetaJN = idxr*LAL_PI;
      phiJL   = idxr*2.*LAL_PI;
      theta1  = idxr*LAL_PI;
      theta2  = (1.-idxr)*LAL_PI;
      phi12   = idxr*2.*LAL_PI;
      chi1    = idxr;
      chi2    = (1.-idxr);
      m1=m1min+idxr*(m1max-m1min);
      m2=m2min+idxr*(m2max-m2min);
      fRef=0.5/(sqrt(6.)*LAL_PI*(m1+m2)*LAL_MTSUN_SI);
      phiRef=idxr*LAL_PI;

      s1xIn=chi1*sin(theta1);
      s1yIn=0.;
      s1zIn=chi1*cos(theta1);
      s2xIn=chi2*sin(theta2)*cos(phi12);
      s2yIn=chi2*sin(theta2)*sin(phi12);
      s2zIn=chi2*cos(theta2);
      inclIn=thetaJN;

      errCode=XLALSimInspiralTransformPrecessingNewInitialConditions(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, thetaJN, phiJL, theta1, theta2, phi12, chi1, chi2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fRef, phiRef);
      ret += compare(chi1*chi1, spin1x*spin1x+spin1y*spin1y+spin1z*spin1z);
      ret += compare(chi2*chi2, spin2x*spin2x+spin2y*spin2y+spin2z*spin2z);
      ret += compare(s1zIn, spin1z);
      ret += compare(s2zIn, spin2z);
      ret += compare(s1xIn*s2xIn+s1yIn*s2yIn+s1zIn*s2zIn, spin1x*spin2x+spin1y*spin2y+spin1z*spin2z);
      ret += compare(s1xIn*s2xIn+s1yIn*s2yIn, spin1x*spin2x+spin1y*spin2y );

      XLALSimInspiralSetFrameAxis(waveFlag, LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT);
      errCode+=XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclIn, s1xIn, s1yIn, s1zIn, s2xIn, s2yIn, s2zIn, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fRef, phiRef,  XLALSimInspiralGetFrameAxis(waveFlag));
      ret += compare(s1xIn*s2xIn+s1yIn*s2yIn+s1zIn*s2zIn,spin1x*spin2x+spin1y*spin2y+spin1z*spin2z);
      ret += compare(chi1*chi1, spin1x*spin1x+spin1y*spin1y+spin1z*spin1z);
      ret += compare(chi2*chi2, spin2x*spin2x+spin2y*spin2y+spin2z*spin2z);
      ret += compare(s1zIn,spin1x*sin(incl)+spin1z*cos(incl));
      ret += compare(s2zIn,spin2x*sin(incl)+spin2z*cos(incl));
      ret += compare(inclIn,incl);

      XLALSimInspiralSetFrameAxis(waveFlag, LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW);
      errCode+=XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclIn, s1xIn, s1yIn, s1zIn, s2xIn, s2yIn, s2zIn, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fRef, phiRef, XLALSimInspiralGetFrameAxis(waveFlag));
      ret += compare(chi1*chi1, spin1x*spin1x+spin1y*spin1y+spin1z*spin1z);
      ret += compare(chi2*chi2, spin2x*spin2x+spin2y*spin2y+spin2z*spin2z);
      ret += compare(s1xIn*s2xIn+s1yIn*s2yIn+s1zIn*s2zIn,spin1x*spin2x+spin1y*spin2y+spin1z*spin2z);
      ret += compare(s1xIn*sin(incl)+s1zIn*cos(incl),spin1x*sin(inclIn)+spin1z*cos(inclIn));
      ret += compare(s2xIn*sin(incl)+s2zIn*cos(incl),spin2x*sin(inclIn)+spin2z*cos(inclIn));
      ret += compare(inclIn,incl);

      XLALSimInspiralSetFrameAxis(waveFlag, LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J);
errCode+=XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclIn, s1xIn, s1yIn, s1zIn, s2xIn, s2yIn, s2zIn, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fRef, phiRef, XLALSimInspiralGetFrameAxis(waveFlag));
      ret += compare(chi1*chi1, spin1x*spin1x+spin1y*spin1y+spin1z*spin1z);
      ret += compare(chi2*chi2, spin2x*spin2x+spin2y*spin2y+spin2z*spin2z);
      ret += compare(s1xIn*s2xIn+s1yIn*s2yIn+s1zIn*s2zIn,spin1x*spin2x+spin1y*spin2y+spin1z*spin2z);
    }

    if ( (ret == 0) && (errCode == 0) )
    {
        fprintf(stdout, "\nInitial condition rotation test passed.\n");
    }
    else
    {
        fprintf(stderr, "\nFAILURE: %u Initial condition rotation test failed.\n", ret);
    }

    return ret;
}
