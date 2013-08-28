/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Thomas Cokelaer, Michele Vallisneri
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

/**
 * \author Cokelaer, T. and Vallisneri, M.
 *
 * \brief Create a waveform based on SpinTaylor model (LALSTPNWaveform),
 * switching to LALSTPNWaveform2 if anything is given as argument.
 *
 * Outputs a file with three columns corresponding to time (in seconds),
 * \f$h_+\f$, and \f$h_\times\f$. The outputfile is <tt>wave1.dat</tt> in the
 * current directory. [This is an obsolete test program that should
 * be refactored with access to the new frameless LALSTPNWaveform.]
 *
 * \heading{Usage}
 *
 * <tt>LALSTPNWaveformTest [switch]</tt>
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

#include <lal/LALSTPNWaveform2.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

extern int newswitch;

int main(int argc,char UNUSED **argv) {
    static LALStatus    mystatus;

    CoherentGW      thewaveform;
    SimInspiralTable    injParams;
    PPNParamStruc       ppnParams;

    const char  *filename = "wave1.dat";
    FILE        *outputfile;
    INT4        i,length;
    REAL8       dt;
    REAL8       a1, a2, phi, shift, phi0;

    memset( &mystatus, 0, sizeof(LALStatus) );
    memset( &thewaveform, 0, sizeof(CoherentGW) );
    memset( &injParams, 0, sizeof(SimInspiralTable) );
    memset( &ppnParams, 0, sizeof(PPNParamStruc) );

    /* --- first we fill the SimInspiral structure --- */

    injParams.mass1 = 10.;
    injParams.mass2 = 10.;

    /* MV-20060224: I believe this is not used in the SpinTaylor code! */
    injParams.f_final = 500.0;
    injParams.f_lower = 40.;

    snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"SpinTaylortwoPN");

    /* this is given in Mpc */
    injParams.distance = 100.;
    /* this should not be zero*/
    injParams.theta0 = 0.01;
    injParams.phi0   = 0.5;

    injParams.inclination = 0.8;
    injParams.polarization   = 0.9;

    injParams.spin1x = 0.1;
    injParams.spin1y = 0.2;
    injParams.spin1z = 0.3;

    injParams.spin2x = 0.4;
    injParams.spin2y = 0.5;
    injParams.spin2z = 0.6;

    ppnParams.deltaT = 1.0 / 4096.0;

    fprintf(stderr, "Lower cut-off frequency used will be %fHz\n", injParams.f_lower);

    if(argc > 1) newswitch=1;
    /* --- now we can call the injection function --- */
    LALGenerateInspiral( &mystatus, &thewaveform, &injParams, &ppnParams );
    if ( mystatus.statusCode )
    {
      fprintf( stderr, "LALSTPNWaveformTest: error generating waveform %d\n", mystatus.statusCode );
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
    fprintf(stderr,"waveform saved in wave1.dat\n" );
    return 0;
}
