/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Thomas Cokelaer, Michele Vallisneri, Riccardo Sturani
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
*  
* Compile with 
* gcc `pkg-config --cflags lal lalsupport lalinspiral` LALSpinInspiralRDTest.c `pkg-config --libs lal lalsupport lalinspiral` -o LALSpinInspiralRDTest.out
*
*/

/**** <lalVerbatim file="LALSpinInspiralRDTestCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Test program \texttt{LALSpinInspiralRDTest.c}}
 * \label{ss:LALSpinInspiralRDTest.c}
 *
 * Create a waveform based on HSpinTaylorRD model (LALSpinInspiralRD).
 * Outputs a file with three columns corresponding to time (in seconds),
 * $h_+$, and $h_\times$.
 * 
 * \subsection*{Usage}
 * 
 * \texttt{LALSpinInspiralRDTest m1 m2 S1x S1y S1z S2x S2y S2z theta0 phi0 finit distance PNorder [outputfile]}
 *
 * The masses are given in solar masses. 
 * The spins are given as \texttt{chi1} and \texttt{chi2} times the unit vector;
 * the direction of the initial orbital angular momentum is assumed along the z axisg;
 * the final frequency is given in Hz, the distance in Mpc.
 * Legal values for \texttt{PNorder} include the strings
 * \texttt{newtonian}, \texttt{oneHalfPN}, \texttt{onePN}, \texttt{onePointFivePN},
 * \texttt{twoPN}, \texttt{twoPointFivePN}, \texttt{threePN}, \texttt{threePointFivePN}.
 * If outputfile is not given it is \texttt{wave1.dat} in the current directory.

 **** </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/NRWaveInject.h>

NRCSID(LALSpinInspiralRDTestC, "$Id: LALSpinInspiralRDTest.c,v 1.1 2009/07/14 09:47:23 ");

int main() {
    static LALStatus    mystatus;

    CoherentGW          thewaveform;
    SimInspiralTable    injParams; // Class defined in LIGOMetadataTables.h
    PPNParamStruc       ppnParams; // Defined in GeneratePPNInspiral.h

    const char  *filename = "wave1.dat";
    FILE        *outputfile;
    INT4        i,length;
    REAL8       dt;
    REAL8       a1, a2, om;
    
    REAL8 max=0.;

    const REAL8 omf=0.065;
    const REAL8 fi =40.;
    REAL8 omi,ff;

    CHAR message[256];

    memset( &mystatus, 0, sizeof(LALStatus) );
    memset( &thewaveform, 0, sizeof(CoherentGW) );
    memset( &injParams, 0, sizeof(SimInspiralTable) );
    memset( &ppnParams, 0, sizeof(PPNParamStruc) );

    /* --- first we fill the SimInspiral structure --- */

    injParams.mass1 = 10.;
    injParams.mass2 = 10.;

    snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"HSpinTaylorRDthreePointFivePN");

    /* this is given in Mpc */
    injParams.distance = 1.;

    /* this is given in Hz*/
    injParams.f_lower  = fi;
    ff=omf/(injParams.mass1+injParams.mass2)/LAL_MTSUN_SI/LAL_PI;

    //injParams.f_final  = ff;
    omi=injParams.f_lower*(injParams.mass1+injParams.mass2)*LAL_MTSUN_SI*LAL_PI;

    fprintf(stdout,"**** Test:  %11.4e < om < %11.4e\n",omi,omf);
    fprintf(stdout,"**** Test:  %11.4e < f  < %11.4e\n",fi,ff);

    injParams.theta0 = 0.;
    injParams.phi0   = 0.1;

    injParams.inclination = 0.0;
    injParams.polarization   = 0.0;

    injParams.spin1x = 0.68;
    injParams.spin1y = 0.;
    injParams.spin1z = -0.68;

    injParams.spin2x = 0.68;
    injParams.spin2y = 0.;
    injParams.spin2z = -0.68;

    /*Spin units are such that multiplying spini by m_i^2 one obtains the physical spin */

    ppnParams.deltaT = 1.0 / 4096.0;
    ppnParams.fStop  = ff;

    /* --- now we can call the injection function --- */
    LALGenerateInspiral( &mystatus, &thewaveform, &injParams, &ppnParams );
    if ( mystatus.statusCode )
    {
      fprintf( stderr, "LALSpinInspiralRDTest: error generating waveform %d\n",mystatus.statusCode );
      fprintf(stderr,"*** Test: h defined %s %11.3e\n",thewaveform.h->name,thewaveform.h->data->data[0]);
      exit( 1 );
    }

    /* --- and finally save in a file --- */

    outputfile = fopen(filename,"w");

    length  = thewaveform.h->data->length;
    
    fprintf(stderr,"*** Test: Length=%d\n",length);

    if (thewaveform.h) {
      fprintf(stderr,"*** Test: h defined %s %11.3e\n",thewaveform.h->name,thewaveform.h->data->data[0]);
      fprintf(stderr,"          phi: %11.3e  shift %11.3e\n",thewaveform.phi->data->data[0],thewaveform.shift->data->data[0]);
    }
    dt      = thewaveform.phi->deltaT;

    //    phi0    = thewaveform.phi->data->data[0];


    for(i = 0; i < length; i++) {
        a1  = thewaveform.h->data->data[2*i];
        a2  = thewaveform.h->data->data[2*i+1];
	om  = thewaveform.f->data->data[i];

	if ( fabs(a1) > max) max = a1;
	if ( fabs(a2) > max) max = a2;

	//phi     = thewaveform.phi->data->data[i] - phi0;
        //shift   = thewaveform.shift->data->data[i];

        fprintf(outputfile,"%e\t%e\t%e\t%e\n",i*dt,a1,a2,om);
    }

    fclose(outputfile);
    fprintf(stderr,"*** Test: Max = %11.3e\n",max);
    fprintf(stderr,"*** Test: waveform saved in %s\n          Final time %11.3e\n", filename, i*dt);

    return 0;
}
