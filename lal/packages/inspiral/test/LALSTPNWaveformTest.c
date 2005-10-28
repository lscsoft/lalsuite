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
    SimInspiralTable    injParams;
    PPNParamStruc       ppnParams;

    FILE 		*outputfile;
    INT4 		i,length;
    REAL8		dt;
    REAL8 		a1, a2, phi, shift;

    memset( &mystatus, 0, sizeof(LALStatus) );
    memset( &thewaveform, 0, sizeof(CoherentGW) );
    memset( &injParams, 0, sizeof(SimInspiralTable) );
    memset( &ppnParams, 0, sizeof(PPNParamStruc) );

    /* --- first we fill the SimInspiral structure --- */
    injParams.mass1 		= 100.0;
    injParams.mass2 		= 1.4;
    injParams.f_final           = 500.0;
    injParams.f_lower 		= 10.0;
    
    LALSnprintf( injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinTaylortwoPN" );
    
    injParams.distance 	        = 1e10;
    injParams.theta0            = 1.5;
    injParams.phi0              = 0.0;

    injParams.spin1x 	= 0.3/sqrt(3.0);
    injParams.spin1y 	= 0.86/sqrt(3.0);
    injParams.spin1z 	= 1.53/sqrt(3.0);

    injParams.spin2x 	= 0.7;
    injParams.spin2y 	= 0.86;
    injParams.spin2z 	= 1.53;
    
    ppnParams.deltaT = 1.0 / 4096.0;

    /* --- now we can call the injection function --- */
    LALGenerateInspiral( &mystatus, &thewaveform, &injParams, &ppnParams );
    if ( mystatus.statusCode )
    {
      fprintf( stderr, "error generating waveform\n" );
      exit( 1 );
    }

    /* --- and finally save in a file --- */
    outputfile = fopen("wave1.dat","w");

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
