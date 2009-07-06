/*
*  Copyright (C) 2007 Chad Hanna, Benjamin Owen
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

/*_______________________________________________________________________________________
 *
 * File Name: InspiralSpinBankwNDTemplateBankTest.c
 *
 * Author: Hanna C. R.
 *
 * Revision: $Id$
 *
 *_______________________________________________________________________________________
 */


/* ------------------------------------ AUTO-DOC ------------------------------------- */
/* ----------------------------------------------------------------------------------- */

/*<lalVerbatim file="InspiralSpinBankwNDTemplateBankTestCV">
  Author: Hanna, C.R.
  $Id$
  </lalVerbatim>*/

/*SUBSECTION - PROGRAM - "InspiralSpinBankwNDTemplateBankTest.c" ---------------------------- <lalLaTeX>
  \subsection{Program \texttt{InspiralSpinBankwNDTemplateBankTest.c}}
  \label{s:InspiralSpinBankwNDTemplateBankTest.c}
  \providecommand{\MATHEMATICA}{$M\scriptstyle{ATHEMATICA}^{\textrm{{\small\textregistered} }}$}
* Tests InpsiralSpinBankwNDTemplateBank().
  </lalLaTeX>*/

  /*SUBSUBSECTION - USAGE - "InspiralSpinBankwNDTemplateBankTest.c" ------------------------- <lalLaTeX>
    \begin{verbatim}
    InspiralSpinBankwNDTemplateBankTest
    \end{verbatim}
    </lalLaTeX>
    END SUBSUBSECTION - USAGE - "InspiralSpinBankwNDTemplateBankTest.c" ----------------------------- */

  /*SUBSUBSECTION - DESCRIPTION - "InspiralSpinBankwNDTemplateBankTest.c" ------------------- <lalLaTeX>
    \subsubsection{Description}
  * This program uses InspiralSpinBankwNDTemplateBank() to generate a template bank from command line
  * parameter input.  It also has the option to make a \MATHEMATICA notebook using
  * LALMath3DPlot() which will plot the 3D template bank.
    </lalLaTeX>
    END SUBSUBSECTION - DESCRIPTION - "InspiralSpinBankTest.c" ----------------------- */

  /*SUBSUBSECTION - OPTIONS - "InspiralSpinBankTest.c" ----------------------- <lalLaTeX>
    \subsubsection{Command line options}
    \begin{description}
    \item[-n]
  * Specifies the minimum smaller mass between 0 and 5.0 $M\odot$.
    \item[-x]
  * Specifies the maximum smaller mass between 0 and 5.0 $M\odot$.
    \item[-m]
  * Specifies the minimum mismatch threshold (typically 0.03) but for the sake of testing
  * it is best to pick a value $O[1]$ to save compiling time.
    \item[-p]
  * Specifies that the program should generate a \MATHEMATICA notebook ``Math3DNotebook.nb".
    \end{description}
    </lalLaTeX>
    END SUBSUBSECTION - OPTIONS - "InspiralSpinBankTest.c" --------------------------- */

  /*SUBSUBSECTION - EXIT CODES ------------------------------------------------------- */
    /* <lalLaTeX>
    \subsubsection*{Exit codes}
    \input{InspiralSpinBankwNDTemplateBankTestCE}
    </lalLaTeX>
    END - SUBSUBSECTION - EXIT CODES ------------------------------------------------- */

  /*SUBSUBSECTION - NOTES - "InspiralSpinBankTest.c" ------------------------- <lalLaTeX>
    \subsubsection{Notes}
    \begin{itemize}
  * \item The metric used in InspiralSpinBank() is only valid for binary systems with a
  * total mass $<15M\odot$ where the minimum larger mass is at least twice the maximum
  * smaller mass.  Choosing mass values that violate these conditions will cause an
  * error message.
  * \item It is unlikely that you will be able to run a \MATHEMATICA notebook that
  * contains more than 10,000 tiles.  Adjust your parameters accordingly if you plan to
  * view a plot.
    \end{itemize}
    </lalLaTeX>
    END - SUBSUBSECTION - NOTES - "InspiralSpinBankTest.c" --------------------------- */
  /*<lalLaTeX>
  \vfill{\footnotesize\input{InspiralSpinBankTestCV}}
  </lalLaTeX>
  END SUBSECTION - PROGRAM - "InspiralSpinBankTest.c" -------------------------------- */

/* ----------------------------------------------------------------------------------- */
/* ----------------------------------- END AUTO-DOC ---------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALConfig.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LALMathematica.h>
#include <lal/LALNoiseModels.h>

extern char *optarg;

/*<lalErrTable file="InspiralSpinBankwNDTemplateBankTestCE">*/
#define INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_ENORM     0
#define INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_EMEM      1
#define INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_ESUB      2

#define INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGENORM  "Normal exit"
#define INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGEMEM   "Memory allocation error"
#define INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGESUB   "Subroutine error"
/*</lalErrTable>*/

NRCSID(LALINSPIRALSPINBANKWNDTEMPLATEBANKTESTC, "$Id$");

int main(int argc, char *argv[]){
  static LALStatus stat;
  INT4 loop = 0; 			/* loop counter */
  Math3DPointList *list = NULL; 	/* Pointer to structure for mathematica plot */
  Math3DPointList *first = NULL;
  InspiralTemplateList *Tiles = NULL;
  InspiralCoarseBankIn CoarseIn;	/* this input stucture is superfluous */
  INT4 ntiles = 0;			/* number of tiles */
  INT2 Math3DPlot = 0;			/* option flag for Mathematica plot */
  INT4 opt = 0;				/* returning value of getopt() */
  INT4 optflag = -1;			/* Command Line option */
  REAL4 PtSize = 0.02;
  REAL8Vector *psd = NULL;
  REAL8 df = 1.0;
  InspiralTemplate inspiralTemplate;
  INT2 printMoments = 0;
  InspiralMomentsEtc moments;
  REAL4 F0 = 164.0; /* Roughly the minimum of the noise curve */
                    /* This should be calculated */
  REAL8 minfreq = 1;
  FILE *plot;
  plot = fopen("plot.dat", "w");

  if ((list = (Math3DPointList *) LALCalloc(1, sizeof(Math3DPointList))) == NULL){
    LALError(&stat, INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGEMEM);
    printf(INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGEMEM);
    return INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_EMEM;
    }

  first = list;

  /* Stuff for calculating the PSD and Noise Moments */
  CoarseIn.fLower = 30;
  CoarseIn.fUpper = 3300;
  CoarseIn.iflso = 0;
  inspiralTemplate.fLower = 30;
  inspiralTemplate.fCutoff = 2000;

 /* Parse options. */
  do{
    optflag++;
    switch (opt) {
      case 'n':
        CoarseIn.mMin = atof( optarg );
        break;
      case 'm':
        CoarseIn.mmCoarse = atof( optarg );
        break;
      case 'x':
        CoarseIn.MMax = atof( optarg );
        break;
      case 'p':
         Math3DPlot = 1;
         break;
      case 's':
         printMoments = 1;
         break;
      default:
         CoarseIn.mmCoarse = 1.0;
         CoarseIn.mMin = 1.0;
         CoarseIn.MMax = 3.0;
         Math3DPlot = 0;
         printMoments = 0;
         break;
      }
    } while ((opt = getopt( argc, argv, "n:m:x:ps" )) != -1);

  CoarseIn.shf.data = NULL;
  memset( &(CoarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  CoarseIn.shf.f0 = 0;
  LALDCreateVector(&stat, &psd, 3300);
  df = 1.0;
  LALNoiseSpectralDensity(&stat, psd, &LALLIGOIPsd, df);
  CoarseIn.shf.data = psd;
  CoarseIn.shf.deltaF = df;

  for(loop = 0; loop < psd->length; loop++){
    if ((psd->data[loop]) && (psd->data[loop] <= minfreq)){
      F0 = (REAL4) CoarseIn.shf.deltaF * loop;
      minfreq = psd->data[loop];
      }
    }
  printf("\nF0=%f\n",F0);


  printf("\n\n\n____________________________________________________________________");
  if (optflag) printf("\nThe parameters you specified (along with defaults) are:\n");
  else {
    printf("\nNOTE: you can provide input parameters in the command line");
    printf("\n\nThe options are:\n");
    printf("\t-x argument = maximum mass of the smaller body in Solar Mass\n");
    printf("\t-n argument = minimum mass of the same\n");
    printf("\t-m argument = Mismatch threshold (e.g. 0.03 is typical");
    printf("\n\t\t      but takes a while to run)\n");
    printf("\t-p (no arg) - calling this option makes a mathematica plot\n");
    printf("\n\t\t      of the tiles\n");
    printf("\t-s (no arg) - this prints the noise moments");
    printf("\nThe default parameters are:\n");
    }
  printf("\nMismatch = %f\nMathematica Plot = %i\nMass two min = %f\nMass two max = %f\nPrint Moments = %i\n",
    CoarseIn.mmCoarse, Math3DPlot, CoarseIn.mMin, CoarseIn.MMax, printMoments);
  printf("______________________________________________________________________\n");

  if (printMoments){
    LALGetInspiralMoments(&stat, &moments, &CoarseIn.shf, &inspiralTemplate);
    printf("\nThe noise moments are:\n");
    for(loop = 1; loop <=17; loop++){
      moments.j[loop] *= pow((inspiralTemplate.fLower/F0), ((7.0-(REAL4) loop)/3.0));
      printf("J[%i] = %f\n", loop, moments.j[loop]);
      }
    }

  printf("\n\nCalling LALInspiralSpinBankwNDTemplateBank()......\n");
  LALInspiralSpinBankwNDTemplateBank(&stat, &Tiles, &ntiles, CoarseIn);
  REPORTSTATUS(&stat);
  if (stat.statusCode){
    LALError(&stat, INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGESUB);
    printf(INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGESUB);
    return INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_ESUB;
  }

  printf("\nThese parameters yeild %i tiles.\n\n\n", ntiles);

  for(loop=0; loop < ntiles; loop++){
    fprintf(plot, "%f\t%f\n", Tiles[loop].params.mass1/LAL_MTSUN_SI, Tiles[loop].params.mass2/LAL_MTSUN_SI);
    }

  /* Mathematica Plot Stuff */
  if (Math3DPlot){
    for(loop=0; loop < (ntiles); loop++){
      list->x = Tiles[loop].params.psi0;
      list->y = Tiles[loop].params.psi3;
      list->z = Tiles[loop].params.beta;
      list->GrayLevel = 1.0*loop/ntiles;
      if ((list = list->next = (Math3DPointList *) LALCalloc(1, sizeof(Math3DPointList))) == NULL){
        LALError(&stat, INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGEMEM);
        printf(INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGEMEM);
        return INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_EMEM;
      }
    }
    list->next = NULL;
    printf("\nCalling LALMath3DPlot()......\n");
    LALMath3DPlot(&stat, first, &ntiles, NULL);
    REPORTSTATUS(&stat);
    if (stat.statusCode){
      LALError(&stat, INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGESUB);
      printf(INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_MSGESUB);
      return INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_ESUB;
    }


    /* Clean Up the memory from the MathPlot3D structure */
    list = first;
    while(list->next){
      first = list->next;
      LALFree(list);
      list = first;
    }
  }

  /* free the last (first?) memory allocated for Math3DPlot. */
  LALFree(list);

  if (stat.statusCode)
    return INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_ESUB;
  else
    return INSPIRALSPINBANKWNDTEMPLATEBANKTESTC_ENORM;

}

