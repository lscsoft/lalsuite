/* <lalVerbatim file="LALSTPNWaveformTestCV">
Author: Cokelaer, T. and Vallisneri M.
$Id: LALSTPNWaveformTest.c,v 1.1 2004/05/05 20:06:23 thomas Exp
</lalVerbatim>*/



/* <lalLaTeX>
  

\subsection{Test program \texttt{LALSTPNWaveformTest.c}}
\label{ss:LALSTPNWaveformTest.c}

Create a waveform based on SpinTaylor model (LALSTPNWaveform).

\subsection*{Usage}
Just change the InspiralTemplate strcuture parameters in the c file. 
</lalLaTeX> */






#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID( LALSTPNWaveformTestC, "$Id: LALSTPNWaveformTest.c,v 1.1 2004/05/05 20:06:23 thomas Exp");


int main() {
    static LALStatus 	mystatus;
    CoherentGW 		thewaveform;
    InspiralTemplate 	parameters;
    PPNParamStruc       ppnParams;
    FILE 		*outputfile;
    INT4 		i,length;
    REAL8		dt;
    REAL8 		a1, a2, phi, shift;

    /* --- first we fill the InspiralTemplate structure --- */
    parameters.mass1 		= 14.0;
    parameters.mass2 		= 10.0;
    parameters.fCutoff 		= 1000.0;
    parameters.fLower 		= 40.0;
    parameters.approximant 	= SpinTaylor;
    parameters.massChoice  	= m1Andm2;
    /* possible orders: newtonian, oneHalfPN, onePN, onePointFivePN,
       twoPN, twoPointFivePN, threePN, threePointFivePN */

    /* twoPointFivePN seems to have termination problems as well as 1 and .5*/

    parameters.order 		= twoPN;
    
    parameters.distance 	= 1e10;
    parameters.tSampling 	= 4000;
    parameters.orbitTheta0 	= 1.5;
    parameters.orbitPhi0 	= 0.0;

    parameters.spin1[0] 	= 0.3/sqrt(3.0);
    parameters.spin1[1] 	= 0.86/sqrt(3.0);
    parameters.spin1[2] 	= 1.53/sqrt(3.0);

    parameters.spin2[0] 	= 0.7;
    parameters.spin2[1] 	= 0.86;
    parameters.spin2[2] 	= 1.53;

    /* --- then the coherent structure --- */
    thewaveform.a 	= 0;
    thewaveform.f 	= 0;
    thewaveform.phi 	= 0;
    thewaveform.shift 	= 0;

    /* --- now we can call the injection function --- */
    LALSTPNWaveformForInjection(&mystatus,&thewaveform,&parameters, &ppnParams);



    /* --- and finally save in a file --- */
    outputfile = fopen("wave1.dat","a");

    length 	= thewaveform.phi->data->length;
    dt 		= thewaveform.phi->deltaT;

    for(i = 0; i < length; i++) {

	a1 	= thewaveform.a->data->data[2*i];
	a2 	= thewaveform.a->data->data[2*i+1];
	phi 	= thewaveform.phi->data->data[i];
	shift 	= thewaveform.shift->data->data[i];

	fprintf(outputfile,"%e\t%e\t%e\n",
		i*dt,
		a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
		a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
    }

    fclose(outputfile);

    return 0;
}
