/*
*  Copyright (C) 2008 Yi Pan
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

#if 0
<lalVerbatim file="TestRingdownCV">
Author: Yi Pan
$Id$
</lalVerbatim>

<lalLaTeX>

\subsection{Program \texttt{TestRingdown.c}}
\label{ss:TestRingdown.c}

Generate a full waveform from inspiral to ring-down.

\subsubsection*{Usage}

\subsubsection*{Description}

\subsubsection*{Exit codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TestRingdownCV}}

</lalLaTeX>
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

NRCSID(METRICTESTC,"$Id$");

int lalDebugLevel = 1;

int main( int argc, char *argv[] )
{
  UINT4 i, j;
  LALStatus status;
  INT4 errcode;
  InspiralTemplate tmplt;
  REAL8 deltaT = 1.0/32768.0;
  UINT4 N = 65536;
  REAL8 deltaF = 1.0/((REAL8)N * deltaT);

  memset( &status, 0, sizeof(LALStatus) );
  memset( &tmplt, 0, sizeof(InspiralTemplate) );

  tmplt.approximant = EOB;
  tmplt.order = 8;

  tmplt.OmegaS = 0.;
  tmplt.Theta = 0.;
  tmplt.ieta = 1;
  tmplt.signalAmplitude = 1.0;

  tmplt.massChoice = m1Andm2;
  tmplt.mass1 = 10.0;
  tmplt.mass2 = 10.0;
  tmplt.chi = 0.000001;
  tmplt.kappa = 0.000001;

  tmplt.fLower= 40.0;
  tmplt.fCutoff = 1000;
  tmplt.tSampling = 1.0 / deltaT;

  tmplt.sourceTheta = LAL_PI / 3.;
  tmplt.sourcePhi = LAL_PI / 6.;
  tmplt.polarisationAngle = LAL_PI / 4.;
  tmplt.startPhase = 0.;
  tmplt.startTime = 0.;

  LALInspiralWaveLength(&status, &N, tmplt);
  LALInspiralParameterCalc(&status, &tmplt);
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  /* Create memory for the inspiral waveform */
  REAL4Vector *omega;
  REAL4Vector *inspwave1;
  REAL4Vector *inspwave2;
  omega		= XLALCreateREAL4Vector( N );
  inspwave1 = XLALCreateREAL4Vector( N );
  inspwave2 = XLALCreateREAL4Vector( N );

  /* Generate T4pn waveforms
  errcode = XLALInspiralComputePTFWaveform( inspwave1, &tmplt);
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALInspiralComputePTFWaveform failed\n" );
    exit( 1 );
  }
  tmplt.startPhase = LAL_PI / 2.;
  errcode = XLALInspiralComputePTFWaveform( inspwave2, &tmplt);
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALInspiralComputePTFWaveform failed\n" );
    exit( 1 );
  } */

  /* Generate inspiral waveforms */
  LALEOBWaveformTemplates(&status, omega, inspwave1, inspwave2, &tmplt);

  /* Attaching position set by omega_match */
  /* Omega_match is given by Eq.(37) of PRD 76, 104049 (2007) */
  /* -0.01 because the current EOBNR 4PN setting can't reach omega_match */
  REAL8 omegamatch = -0.01 + 0.133 + 0.183 * tmplt.eta + 0.161 * tmplt.eta * tmplt.eta;
  UINT4 attpos = 0;
  for (j = 0; j < N; ++j)
  {
    if(omega->data[j] > omegamatch)
	{
	  attpos = j - 1;
	  break;
	}
  }

  FILE *omegafile;
  omegafile = fopen( "myomega.dat", "w" );
  if( omegafile != NULL )
  {
	for (j = 0; j < N; ++j)
	{
	  fprintf( omegafile, "%f\n", omega->data[j]);
	}
  }
  fclose(omegafile);

  /* Create memory for the QNM frequencies */
  COMPLEX8Vector *modefreqs;
  modefreqs = XLALCreateCOMPLEX8Vector( 3 );
  errcode = XLALGenerateQNMFreq( modefreqs, &tmplt, 2, 2, 3 );
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALGenerateQNMFreq failed\n" );
    exit( 1 );
  }

  /* Ringdown signal length: 10 times the decay time of the n=0 mode */
  UINT4 Nrdwave = 10 / modefreqs->data[0].im / deltaT;
  /* Patch length, centered around the matching point "attpos" */
  UINT4 Npatch = 11;

  /* Create memory for the ring-down and full waveforms, and derivatives of inspirals */
  REAL4Vector			*rdwave1;
  REAL4Vector			*rdwave2;
  REAL4Vector			*fullwave1;
  REAL4Vector			*fullwave2;
  REAL4Vector			*inspwave;
  REAL4Vector			*dinspwave;
  REAL4Vector			*ddinspwave;
  REAL4VectorSequence	*inspwaves1;
  REAL4VectorSequence	*inspwaves2;

  rdwave1 = XLALCreateREAL4Vector( Nrdwave );
  rdwave2 = XLALCreateREAL4Vector( Nrdwave );
  fullwave1 = XLALCreateREAL4Vector( attpos + Nrdwave - 1 );
  fullwave2 = XLALCreateREAL4Vector( attpos + Nrdwave - 1 );
  inspwave = XLALCreateREAL4Vector( Npatch );
  dinspwave = XLALCreateREAL4Vector( Npatch );
  ddinspwave = XLALCreateREAL4Vector( Npatch );
  inspwaves1 = XLALCreateREAL4VectorSequence( 3, (Npatch + 1) / 2 );
  inspwaves2 = XLALCreateREAL4VectorSequence( 3, (Npatch + 1) / 2 );

  /* Generate derivatives of the last part of inspiral waves */
  /* Take the last part of inspwave1 */
  for (j = 0; j < Npatch; j++)
  {
	inspwave->data[j] = inspwave1->data[attpos - (Npatch + 1) / 2 + j];
  }
  /* Get derivatives of inspwave1 */
  errcode = XLALGenerateWaveDerivatives( dinspwave, ddinspwave, inspwave, &tmplt );
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALGenerateQNMFreq failed\n" );
    exit( 1 );
  }
  for (j = 0; j < (Npatch + 1) / 2; j++)
  {
	inspwaves1->data[j] = inspwave->data[j];
	inspwaves1->data[j + (Npatch + 1) / 2] = dinspwave->data[j];
	inspwaves1->data[j + 2 * (Npatch + 1) / 2] = ddinspwave->data[j];
  }
  /* Take the last part of inspwave2 */
  for (j = 0; j < Npatch; j++)
  {
	inspwave->data[j] = inspwave2->data[attpos - (Npatch + 1) / 2 + j];
  }
  /* Get derivatives of inspwave2 */
  errcode = XLALGenerateWaveDerivatives( dinspwave, ddinspwave, inspwave, &tmplt );
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALGenerateQNMFreq failed\n" );
    exit( 1 );
  }
  for (j = 0; j < (Npatch + 1) / 2; j++)
  {
	inspwaves2->data[j] = inspwave->data[j];
	inspwaves2->data[j + (Npatch + 1) / 2] = dinspwave->data[j];
	inspwaves2->data[j + 2 * (Npatch + 1) / 2] = ddinspwave->data[j];
  }


  /* Generate ring-down waveforms */
  errcode = XLALInspiralRingdownWave( rdwave1, rdwave2, &tmplt, inspwaves1, inspwaves2,
								  modefreqs, 3);
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALInspiralComputePTFWDeriv failed\n" );
    exit( 1 );
  }
  /* Generate full waveforms, by stitching inspiral and ring-down waveforms */
  for (j = 0; j < attpos; ++j)
  {
	fullwave1->data[j] = inspwave1->data[j];
	fullwave2->data[j] = inspwave2->data[j];
  }
  for (j = 1; j < Nrdwave; ++j)
  {
	fullwave1->data[j + attpos - 1] = rdwave1->data[j];
	fullwave2->data[j + attpos - 1] = rdwave2->data[j];
  }


  /* Write waveforms to file */
  FILE *inspfile;
  FILE *rdfile;
  FILE *fullfile;
  inspfile = fopen( "myinspwave.dat", "w" );
  if( inspfile != NULL )
  {
	for (j = 0; j < attpos; ++j)
	{
	  fprintf( inspfile, "%e      %e\n", inspwave1->data[j], inspwave2->data[j]);
	}
  }
  fclose(inspfile);
  rdfile = fopen( "myrdwave.dat", "w" );
  if( rdfile != NULL )
  {
	for (j = 0; j < Nrdwave; ++j)
	{
	  fprintf( rdfile, "%e      %e\n", rdwave1->data[j], rdwave2->data[j]);
	}
  }
  fclose(rdfile);
  fullfile = fopen( "myfullwave.dat", "w" );
  if( fullfile != NULL )
  {
	for (j = 0; j < attpos + Nrdwave -1; ++j)
	{
	  fprintf( fullfile, "%e      %e\n", fullwave1->data[j], fullwave2->data[j]);
	}
  }
  fclose(fullfile);

  /* Free memory */
  XLALDestroyCOMPLEX8Vector( modefreqs );
  XLALDestroyREAL4Vector( fullwave1 );
  XLALDestroyREAL4Vector( fullwave2 );
  XLALDestroyREAL4Vector( rdwave1 );
  XLALDestroyREAL4Vector( rdwave2 );
  XLALDestroyREAL4VectorSequence( inspwaves1 );
  XLALDestroyREAL4VectorSequence( inspwaves2 );
  XLALDestroyREAL4Vector( omega );
  XLALDestroyREAL4Vector( inspwave1 );
  XLALDestroyREAL4Vector( inspwave2 );

  return 0;
}
