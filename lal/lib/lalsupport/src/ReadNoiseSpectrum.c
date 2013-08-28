/*
*  Copyright (C) 2007 Jolien Creighton, Patrick Brady
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
#include <errno.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/Interpolate.h>
#include <lal/ReadNoiseSpectrum.h>

/*********************************************************************
 * STATIC FUNCTION to locate the point nearest to the desired
 * frequency for the interpolation step below
 *********************************************************************/
static UINT4 mylocate(REAL8 *farray, REAL8 target, UINT4 n)
{
  UINT4 i=0;
  UINT4 jmin=0;

  for (i=0 ; i<n ; i++){
    if ( (farray[i]-target) > 0 ){
      break;
    }
  }

  jmin=i;

  if (i>0){
     if ( (target - farray[i-1]) < (farray[i] - target) ){
         jmin =  i-1;
     }
  }

  return jmin;
}


/*********************************************************************
 * MAIN FUNCTION to get populate the amplitude spectrum from the data
 * file specified.
 *********************************************************************/


/**
 * \ingroup ReadNoiseSpectrum_h
 * \brief Function to read in noise spectrum from a formatted ascii file and return the
 * amplitude noise spectrum in \f$\textrm{strain}/\sqrt{\textrm{Hz}}\f$.
 *
 * ### Description ###
 *
 * <tt>LALReadNoiseSpectrum()</tt> fills the contents of the REAL4FrequencySeries
 * \c spectrum from data read from formatted ascii file with name \c fname.
 * The ascii file should have a header (greater than or equal to one line) which is
 * indicated by a \f$\#\f$ at the start of the line.   The first line of the file must
 * contain the number of points at which the spectrum is sampled.  If the spectrum
 * is sampled at 500 different points,  then the first line would be
 * \code
 * # npoints=500
 * \endcode
 * Replace 500 by the number of sample points in your particular data.
 *
 * The REAL4FrequencySeries \c spectrum should be allocated before calling the
 * routine which uses the \c length and metadata information to determine the
 * evenly sampled output that is reqruired.   The function does nearest neighbor
 * interpolation to get the points in the outpu frequency series.
 *
 */
    void
LALReadNoiseSpectrum(LALStatus *stat, REAL4FrequencySeries *spectrum, CHAR *fname)
{

    FILE *fp=NULL; /* where the spectrum data is stored */
    CHAR line[LALREADNOISESPECTRUM_MAXLINELENGTH];
    INT4 npoints;   /* number of points in the specrtum */
    UINT4 j;         /* dummy index*/
    REAL8 *f=NULL;  /* dummy variable for frequency values */
    REAL8 *s=NULL;  /* dummy variable for spectrum values */
    REAL4 freq, myfmin, df;
    UINT4 location;
    DInterpolateOut  intOutput;
    DInterpolatePar  intParams;

    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* this is the file containing the spectrum data */
    if ( !(fp = LALOpenDataFile( fname )) )
    {
        ABORT(stat, LALREADNOISESPECTRUMH_EOPEN, LALREADNOISESPECTRUMH_MSGEOPEN);
    }

    /* read in the first line of the file */
    if (fgets(line,sizeof(line),fp) == NULL)
    {
        ABORT(stat, LALREADNOISESPECTRUMH_EPARS, LALREADNOISESPECTRUMH_MSGEPARS);
    }

    /* check that it has the right format */
    if ( (strstr(line,"# npoints=")) == NULL)
    {
        ABORT(stat, LALREADNOISESPECTRUMH_EPARS, LALREADNOISESPECTRUMH_MSGEPARS);
    }
    sscanf(line,"# npoints=%i", &npoints);

    /* memory for the input data */
    f = (REAL8 *) LALMalloc( npoints * sizeof(REAL8) );
    s = (REAL8 *) LALMalloc( npoints * sizeof(REAL8) );

    /* read data into arrays */
    j=0;
    while (1) {
        if (fgets(line,sizeof(line),fp)==NULL) {
            LALFclose(fp);
            break;
        }
        if (line[0] != '#'){
            sscanf(line,"%lf %lf\n",&f[j],&s[j]);
            j++;
        }
    }

    /* populate the frequency series */
    intParams.n=4;
    location = 0;
    myfmin = spectrum->f0;
    df = spectrum->deltaF;
    for(j=0 ; j < spectrum->data->length  ; j++) {
        freq = myfmin + ((REAL4) j)*df;

        /* if the frequency is above lowest in noise file ... */
        if ( freq >= f[0] )
        {
            /* ... interpolate to desired frequency ... */
            location=mylocate(f,(REAL8)freq,npoints);
            if (location > (npoints-intParams.n)){
                location = (npoints-intParams.n);
            }
            else if ( location < (intParams.n/2) ){
                location = 0;
            }
            else {
                location-=(intParams.n/2);
            }
            intParams.x = &f[location];
            intParams.y = &s[location];
            LALDPolynomialInterpolation(stat->statusPtr, &intOutput, freq, &intParams);
            CHECKSTATUSPTR (stat);
            spectrum->data->data[j] = (REAL4)(intOutput.y);
        }

        /* ...... otherwise,  just fill be lowest available value */
        else {
            spectrum->data->data[j] = s[0];
        }
    }

    /* set the units to strain / sqrt(Hz) */
    memset( &(spectrum->sampleUnits), 0, sizeof(LALUnit) );
    spectrum->sampleUnits.unitNumerator[LALUnitIndexStrain] = 1;
    spectrum->sampleUnits.unitNumerator[LALUnitIndexSecond] = 1;
    spectrum->sampleUnits.unitDenominatorMinusOne[LALUnitIndexSecond] = 1;

    /* .... and clean up */
    LALFree( f );
    LALFree( s );

    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}
