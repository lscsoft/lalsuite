/*
*  Copyright (C) 2007 Patrick Brady, Reinhard Prix, Thomas Cokelaer
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
#include <string.h>
#include <math.h>
#include <lal/LALConstants.h>

#define NLINES 16000
#define SNR 8.0

#define USAGE "Usage: %s [options]\n"\
        "--help                           Print this help message\n" \
        "--masses m1 m1       Masses of binary elements\n" \
        "--spectrum filename  File containing amplitude spectrum m/rtHz \n" \
        "--length arm_length  Length of the arm in meters \n" \

int main ( int argc, char *argv[])
{
    FILE *fpin = NULL;
    int i,nlines=NLINES,arg=1,record=0;
    float f[NLINES],amp[NLINES],dum1;
    float m1, m2, eta, mtot, mu, arm_length=1.0;
    float ins_amp,d_optimal;

    m1 = m2 = 1.4;

    /*******************************************************************
     * PARSE ARGUMENTS (arg stores the current position)               *
     *******************************************************************/

    if (argc <= 1){
        fprintf(stderr,  USAGE, *argv );
        return 0;
    }

    while ( arg < argc ) {
        /*********************************************************
         * File containing veto information metadata
         *********************************************************/
        if ( !strcmp( argv[arg], "--help" ) ) {
            fprintf(stderr,  USAGE, *argv );
            return 0;
        }
        /*********************************************************
         * Masses
         *********************************************************/
        else if ( !strcmp( argv[arg], "--masses" ) ) {
                arg++;
                m1 = atof(argv[arg++]);
                m2 = atof(argv[arg++]);
        }
        /*********************************************************
         * Armlength
         *********************************************************/
        else if ( !strcmp( argv[arg], "--length" ) ) {
                arg++;
                arm_length = atof(argv[arg++]);
        }
        /*********************************************************
         * Spectrum [m / sqert(Hz)]
         *********************************************************/
        else if ( !strcmp( argv[arg], "--spectrum" ) ) {
                arg++;
                fpin=fopen(argv[arg++],"r");
        }
        /*********************************************************
         * Switch to generate output for plotting later
         *********************************************************/
        else if ( !strcmp( argv[arg], "--record" ) ) {
                arg++;
                record=1;
        }
        /* Check for unrecognized options. */
        else if ( argv[arg][0] == '-' ) {
            fprintf(stderr,  USAGE, *argv );
            return 1;
        }
    } /* End of argument parsing loop. */


    /* Compute the inspiral parameters */
    mtot = m1+m2;
    eta = m1*m2/(mtot*mtot);
    mu = eta * mtot;
    ins_amp = (LAL_MTSUN_SI * LAL_C_SI / (1.0e6 *  LAL_PC_SI))
        * sqrt( 5.0*mu / 96.0 ) * ( pow( mtot/(LAL_PI*LAL_PI) , 0.33333 ) /
            pow(LAL_MTSUN_SI, 1.0 / 6.0) ) ;

    /* Do the integral over frequency */
    dum1=0;
    for(i=1;i<nlines;i++){
        if ( fscanf(fpin,"%f %f",&f[i],&amp[i]) == EOF )
            break;
        dum1+=pow(f[i],-7.0/3.0)/(amp[i]*amp[i]);
        /* Assume ISCO is at 200 (20 Msun / mtot) Hz */
        if (f[i] * mtot > 4000)
            break;
    }

    /* print the information to stderr */
    d_optimal= arm_length*ins_amp*sqrt(2.0 * dum1)/SNR;
    fprintf(stderr,"Sensitive to optimally oriented binary at %0.4f Mpc\n",
            d_optimal);
    fprintf(stderr,"Effective sensitivity over sky %0.4f Mpc\n",
            d_optimal/sqrt(5.0));

    /* and to stdout if desired */
    if ( record ){
        fprintf(stdout,"%f %f\n",mtot,d_optimal);
    }

    return 1;
}
