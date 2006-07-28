/* <lalVerbatim file="LALSTPNWaveformTestCV">
Author: Cokelaer, T. and Vallisneri, M.
$Id: LALSTPNWaveformTest.c,v 1.1 2004/05/05 20:06:23 thomas Exp
</lalVerbatim>*/

/* <lalLaTeX>

\subsection{Test program \texttt{LALSTPNWaveformTest.c}}
\label{ss:LALSTPNWaveformTest.c}

Create a waveform based on SpinTaylor model (LALSTPNWaveform).
Outputs a file with three columns corresponding to time (in seconds),
$h_+$, and $h_\times$.

\subsection*{Usage}

\texttt{LALSTPNWaveformTest m1 m2 S1x S1y S1z S2x S2y S2z theta0 phi0 finit distance PNorder [outputfile]}

The masses are given in solar masses.
The spins are given as \texttt{chi1} and \texttt{chi2} times the unit vector;
the direction of the initial orbital angular momentum is given with the standard polar angles;
the final frequency is given in Hz, the distance in Mpc.
Legal values for \texttt{PNorder} include the strings
\texttt{newtonian}, \texttt{oneHalfPN}, \texttt{onePN}, \texttt{onePointFivePN},
\texttt{twoPN}, \texttt{twoPointFivePN}, \texttt{threePN}, \texttt{threePointFivePN}.
If outputfile is not given it is \texttt{wave1.dat} in the current directory.

</lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID(LALSTPNWaveformTestC, "$Id: LALSTPNWaveformTest.c,v 1.1 2004/05/05 20:06:23 thomas Exp");

int main(int argc,char **argv) {
    static LALStatus    mystatus;
    CoherentGW      thewaveform;
    SimInspiralTable    injParams;
    PPNParamStruc       ppnParams;

    const char        *filename = "wave1.dat";
    FILE        *outputfile;
    INT4        i,length;
    REAL8       dt;
    REAL8       a1, a2, phi, shift, phi0;

    memset( &mystatus, 0, sizeof(LALStatus) );
    memset( &thewaveform, 0, sizeof(CoherentGW) );
    memset( &injParams, 0, sizeof(SimInspiralTable) );
    memset( &ppnParams, 0, sizeof(PPNParamStruc) );

    if(argc < 14) {
        printf("Usage: %s m1 m2 S1x S1y S1z S2x S2y S2z theta0 phi0 finit distance PNorder [outputfile]\n\n",argv[0]);
        printf("The spins are given as chi1 and chi2 times the unit vector;\n");
        printf("the initial orbital angular momentum is given as standard polar angles,\n");
        printf("the initial frequency in Hz, the distance in Mpc.\n");
        printf("Legal values for PNorder include newtonian, oneHalfPN, onePN,\n");
        printf("onePointFivePN, twoPN, twoPointFivePN, threePN, threePointFivePN.\n");
        printf("If outputfile is not given it is wave1.dat in the current directory.\n");
        exit(1);
    }
    
    if(argc > 14) {
        filename = argv[14];
    }

    /* --- first we fill the SimInspiral structure --- */
    
    injParams.mass1 = atof(argv[1]);
    injParams.mass2 = atof(argv[2]);

    /* MV-20060224: I believe this is not used in the SpinTaylor code! */

    injParams.f_final = 500.0;

    injParams.f_lower = atof(argv[11]);
    
    LALSnprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),argv[13]);
    
    /* MV-20060224: this is given in Mpc */
    
    injParams.distance = atof(argv[12]);

    injParams.theta0 = atof(argv[9]);
    injParams.phi0   = atof(argv[10]);

    injParams.spin1x = atof(argv[3]);
    injParams.spin1y = atof(argv[4]);
    injParams.spin1z = atof(argv[5]);
    
    injParams.spin2x = atof(argv[6]);
    injParams.spin2y = atof(argv[7]);
    injParams.spin2z = atof(argv[8]);

    ppnParams.deltaT = 1.0 / 4096.0;

    /* --- now we can call the injection function --- */
    LALGenerateInspiral( &mystatus, &thewaveform, &injParams, &ppnParams );
    if ( mystatus.statusCode )
    {
      fprintf( stderr, "LALSTPNWaveformTest: error generating waveform\n" );
      exit( 1 );
    }

    /* --- and finally save in a file --- */

    outputfile = fopen(filename,"w");

    length  = thewaveform.phi->data->length;
    dt      = thewaveform.phi->deltaT;

    phi0    = thewaveform.phi->data->data[0];

    for(i = 0; i < length; i++) {
        a1  = thewaveform.a->data->data[2*i];
        a2  = thewaveform.a->data->data[2*i+1];
        phi     = thewaveform.phi->data->data[i] - phi0;
        shift   = thewaveform.shift->data->data[i];

        fprintf(outputfile,"%e\t%e\t%e\n",
            i*dt,
            a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
            a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
    }

    fclose(outputfile);

    return 0;
}
