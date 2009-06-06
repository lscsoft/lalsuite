/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Thomas Cokelaer
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

/* <lalVerbatim file="GetOrientationEllipseCV">
Author: Thomas Cokelaer
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{GetOrientationEllipse.c}}
\label{ss:GetOrientationEllipse.c}

Test code for the \texttt{bank} modules.


\subsubsection*{Usage}
\begin{verbatim}
./getOrientationEllipse
\end{verbatim}

\subsubsection*{Description}

This code illustrates the use of several functions such as
\texttt{LALInspiralParameterCalc}, \texttt{LALGetInspiralMoments},
and \texttt{LALInspiralComputeMetric}. It shows how to defined
a suitable \texttt{InspiralCoarseBankIn} structure so as to extract
the metric components for a set of binary parameters. In this
example, we first declare all the relevant parameter needed
(minimum  and maximum mass, fLower, design sensitivity curve
and so on), which can be changed by the user before compilation.

Then, a loop spans a square parameter space defined by tau0 in the
range [.1,40] seconds and tau3 in [1, 2] seconds. For each set of
parameter, the metric is computed and the code prints on stdout
the value of the coordinate used (tau0, tau3) and the orientation
of the metric in degrees. We do not check whether a template is valid
or not in this code but one could have use a function such as
\texttt{LALInspiralValidtemplate} to do so.

\subsubsection*{Notes}
\vfill{\footnotesize\input{GetOrientationEllipseCV}}
</lalLaTeX> */


#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

#include <lal/LALRCSID.h>
NRCSID (GETORIENTATIONELLIPSEC,"$Id$");



INT4 lalDebugLevel=33;

int
main(int argc, char **argv)
{
  INT4 arg;
  /* top-level status structure */
  static LALStatus status;
  /* Structure specifying the nature of the bank needed */
  static InspiralCoarseBankIn coarseIn;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  UINT4   j, numPSDpts=262144/4/4; /*Size of the vectors*/
  InspiralTemplate tempPars;
  InspiralMetric metric;
  InspiralMomentsEtc moments;


  /* In order to use any functions related to the metric, we need to fill the
   * coarseIn structure which is defined here. The most important are the
   * mMin, mMax and fLower parameters which influence the size of the
   * parameter space and metric. Of course, mmCoarse, the minimal match is
   * also relevant. */
  coarseIn.LowGM        = -2;
  coarseIn.HighGM       = 6;
  coarseIn.fLower       = 40.L;
  coarseIn.fUpper       = 2047.L;
  coarseIn.tSampling    = 4096.L;
  coarseIn.order        = LAL_PNORDER_TWO;
  coarseIn.space        = Tau0Tau3;
  coarseIn.mmCoarse     = 0.97;
  coarseIn.mmFine       = 0.97;
  coarseIn.iflso        = 0.0L;
  coarseIn.mMin         = 3;
  coarseIn.mMax         = 30.0;
  coarseIn.MMax         = coarseIn.mMax * 2.;
  coarseIn.massRange    = MinMaxComponentMass;
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin       = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
/*  coarseIn.psi0Min      = 1.e0;
  coarseIn.psi0Max      = 2.5e4;
  coarseIn.psi3Min      = -1e4;
  coarseIn.psi3Max      = -10;
  coarseIn.alpha        = 0.L;
  coarseIn.numFcutTemplates = 4;
*/

  /* init the vector for the PSD */
  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
  coarseIn.shf.deltaF = coarseIn.tSampling / (2.*(REAL8) coarseIn.shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status,
			   coarseIn.shf.data,
			   noisemodel, coarseIn.shf.deltaF );



  tempPars.totalMass = coarseIn.MMax;
  tempPars.eta = 0.25;
  tempPars.ieta = 1.L;
  tempPars.fLower = coarseIn.fLower;
  tempPars.massChoice = totalMassAndEta;

  metric.space = coarseIn.space;
  LALInspiralSetParams( &status, &tempPars, coarseIn );


  /* */

  {
    REAL4 tau0 = 0.;
    REAL4 tau3 = 0.;
    INT4 valid = 0;
    InspiralBankParams bankPars, bankParsOld;

    /* Get the limit of the parameter space*/
    LALInspiralSetSearchLimits( &status, &bankPars, coarseIn );
    /* Loop over a rectangle with pre-defined range in tau_0 and tau_3. We will
     * then check if this is a valid template. */
    for (tau0=.8; tau0<10; tau0+=0.01){
      for (tau3=.3; tau3<2; tau3+=0.01){
        valid = 1;
        bankPars.x0 = tau0;
        bankPars.x1 = tau3;

        LALInspiralValidParams(&status, &valid, bankPars, coarseIn );
	if  (valid==1)
        {
	  tempPars.t0 = tau0;
  	  tempPars.t3 = tau3;
	  LALInspiralParameterCalc( &status, &tempPars );
          /* Even for non physical template, we can get the metric componenents*/
  	  LALGetInspiralMoments( &status, &moments, &coarseIn.shf, &tempPars );
	  LALInspiralComputeMetric( &status, &metric, &tempPars, &moments );
          /*Now, we simply ouput the values we are interested in : t0, t3 and
           * angle. */
          printf("%f  %f %f\n", tempPars.t0, tempPars.t3, metric.theta*180/3.14159);
  	  fflush(stdout);
        }

      }
    }
  }


  /* Free the list, and exit. */

  LALDDestroyVector( &status, &(coarseIn.shf.data) );
  LALCheckMemoryLeaks();
  return(0);
}



