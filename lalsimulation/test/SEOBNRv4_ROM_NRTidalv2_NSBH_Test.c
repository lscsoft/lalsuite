 /*
 *  Copyright (C) 2019 Andrew Matas
 *  Test code for LALSimIMRSEOBNRv4ROM_NSBHAmplitudeCorrection
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/LALConstants.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>


#include <lal/LALSimNoise.h>
#include <lal/ComplexFFT.h>

#include <lal/ComplexFFT.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include "../lib/LALSimIMRSEOBNRv4ROM_NSBHAmplitudeCorrection.c"

#define MYUNUSED(expr) do { (void)(expr); } while (0)


bool approximatelyEqual(REAL8 a, REAL8 b, REAL8 epsilon);
bool approximatelyEqual(REAL8 a, REAL8 b, REAL8 epsilon) {
  if (a == 0)
    return fabs(b) < epsilon;
  else if (b == 0)
    return fabs(a) < epsilon;
  else
    return !gsl_fcmp(a, b, epsilon);
}
bool approximatelyEqualREAL8Sequence(REAL8Sequence * a, REAL8Sequence * b, REAL8 epsilon);
bool approximatelyEqualREAL8Sequence(REAL8Sequence * a, REAL8Sequence * b, REAL8 epsilon) {
  if (a->length != b->length)
     return false;
  for (UINT4 i=0; i<a->length; i++){
     if (!approximatelyEqual(a->data[i],b->data[i],epsilon))
         return false;
  }
  return true;
}


void print_difference(const char *name, REAL8 u, REAL8 u_expected);
void print_difference(const char *name, REAL8 u, REAL8 u_expected) {
  printf("%s: %-20.17g\t%-20.17g\t%-20.17g\n", name, u, u_expected, u - u_expected);
}

void print_difference_FrequencySequence(const char *name, REAL8Sequence *u, REAL8Sequence *u_expected,REAL8Sequence *freqs);
void print_difference_FrequencySequence(const char *name, REAL8Sequence *u, REAL8Sequence *u_expected,REAL8Sequence *freqs){
    for (UINT4 i=0; i<freqs->length; i++){
  printf("%s [%2.5f Hz]: %-20.17g\t%-20.17g\t%-20.17g\n", name, freqs->data[i], u->data[i], u_expected->data[i], u->data[i] - u_expected->data[i]);
    }
}
static void Test_amplitude_correction_disruptive(void);
static void Test_amplitude_correction_disruptive(void){
    REAL8 m1_SI=1.4*LAL_MSUN_SI;
    REAL8 m2_SI=1.4*LAL_MSUN_SI;
    REAL8 chi1=0;
    REAL8 lambda2=1000;

    REAL8 deltaF=100;
    REAL8 Fmax=1000;
    int nfreqs = (int)(Fmax/deltaF);

    REAL8Sequence *freqs=NULL;
    freqs = XLALCreateREAL8Sequence(nfreqs);
    for (int i=0; i<nfreqs; i++){
        freqs->data[i]=i*deltaF;
    }

    REAL8Sequence *amp_tidal = NULL;
    amp_tidal = XLALCreateREAL8Sequence(freqs->length);
    int ret=XLALSEOBNRv4ROMNSBHAmplitudeCorrectionFrequencySeries(
            amp_tidal, freqs,
            m1_SI, m2_SI, chi1, lambda2);
    MYUNUSED(ret);
    
    REAL8Sequence *expected_amp_tidal=NULL;
    expected_amp_tidal=XLALCreateREAL8Sequence(freqs->length);
    expected_amp_tidal->data[0]=0.98956307697285817;
    expected_amp_tidal->data[1]=0.98623701074185499;
    expected_amp_tidal->data[2]=0.98187040424514516;
    expected_amp_tidal->data[3]=0.97615189747808495;
    expected_amp_tidal->data[4]=0.96868716084731876;
    expected_amp_tidal->data[5]=0.95898405100871376;
    expected_amp_tidal->data[6]=0.94644042313762244;
    expected_amp_tidal->data[7]=0.93033932007923714;
    expected_amp_tidal->data[8]=0.90985889942914522;
    expected_amp_tidal->data[9]=0.8841072573277694;

    print_difference_FrequencySequence("Amplitude Correction",amp_tidal,expected_amp_tidal,freqs);
    assert(approximatelyEqualREAL8Sequence(amp_tidal,expected_amp_tidal,1e-9));
}

int main(int argc, char *argv[]) {
  MYUNUSED(argc);
  MYUNUSED(argv);

  Test_amplitude_correction_disruptive();

  return 0;
}

