/* [PURPOSE]
 * 	test file for the LALSTPNWaveform function. 
 * [USAGE]
 * 	change paramters in the file, compile and launch the executable
 * [INPUT]
 * 	none	
 * [OUTPUT]
 * 	write the waveform in the file wave1.dat
 * [COMPILATION]
 *      gcc  LALSTPNWaveformTest LALSTPNWaveformTest.c -L<lib lal directory> 
 *      -llal -I<include lal directory>
 * [AUTHORS]
 * 	Michele Vallisneri, Thomas Cokelaer, April 2004
 *
 * [COMMENTS]
 * 	we don't put that file into the Makefile.am. It is indeed a temporary 
 * 	test file dedicated to one function namely LALSTPNWaveformForInjection
 * */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

int main(int argc, char **argv) {
    static LALStatus mystatus;
    CoherentGW thewaveform;
    InspiralTemplate parameters;
    FILE *outputfile;
    int i,length;
    double dt;
    double a1,a2,phi,shift;

    /* --- first we fill the InspiralTemplate structure --- */
    parameters.mass1 	= 14.0;
    parameters.mass2 	= 10.0;
    parameters.fCutoff 	= 1000;
    parameters.fLower 	= 100.0;
    parameters.approximant = SpinTaylorPN;
    parameters.massChoice  = m1Andm2;
    /* possible orders: newtonian, oneHalfPN, onePN, onePointFivePN,
       twoPN, twoPointFivePN, threePN, threePointFivePN */

    /* twoPointFivePN seems to have termination problems as well as 1 and .5*/

    parameters.order = twoPN;
    
    parameters.distance 	= 1e10;
    parameters.tSampling 	= 4000;
    parameters.orbitTheta0 	= 1.5;
    parameters.orbitPhi0 	= 0.0;

    parameters.spin1[0] 	= 0.0/sqrt(3.0);
    parameters.spin1[1] 	= 0.0/sqrt(3.0);
    parameters.spin1[2] 	= 0.0/sqrt(3.0);

    parameters.spin2[0] 	= 0.0;
    parameters.spin2[1] 	= 0.0;
    parameters.spin2[2] 	= 0.0;

    /* --- then the coherent structure --- */
    thewaveform.a = 0;
    thewaveform.f = 0;
    thewaveform.phi = 0;
    thewaveform.shift = 0;

    /* --- now we can call the injection function --- */
    LALSTPNWaveformForInjection(&mystatus,&thewaveform,&parameters);
    REPORTSTATUS(&mystatus);


    /* --- and finally save in a file --- */
    outputfile = fopen("wave1.dat","a");

    length 	= thewaveform.phi->data->length;
    dt 		= thewaveform.phi->deltaT;

    for(i = 0; i < length; i++) {

	a1 	= thewaveform.a->data->data[2*i];
	a2 	= thewaveform.a->data->data[2*i+1];
	phi 	= thewaveform.phi->data->data[i];
	shift 	= thewaveform.shift->data->data[i];

	fprintf(outputfile,"%le\t%le\t%le\n",
		i*dt,
		a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
		a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
    }

    fclose(outputfile);
return 0;
}
