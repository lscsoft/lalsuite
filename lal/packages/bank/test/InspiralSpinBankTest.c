/* <lalVerbatim file="InspiralSpinBankTestCV">
 * Author: Hanna, C.R. and Owen, B. J.
 * $Id$
 * </lalVerbatim> */ 

/* <lalLaTeX>
 * \subsection{Program \texttt{InspiralSpinBankTest}}
 * \label{s:InspiralSpinBankTest.c}
 * \providecommand{\MATHEMATICA}{$M\scriptstyle{ATHEMATICA}^{\textrm{{\small\textregistered} }}$}
 * Tests InpsiralSpinBank(). 
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * InspiralSpinBankTest
 * \end{verbatim} 
 *
 * This program uses InspiralSpinBank() to generate a template bank from
 * command line parameter input.  It also has the option to make a
 * \MATHEMATICA notebook using LALMath3DPlot() which will plot the 3D
 * template bank. 
 *
 * \subsubsection{Command line options}
 * \begin{description}
 * \item[-n] 
 * Specifies the minimum smaller mass between 0 and 5.0 $M\odot$.
 * \item[-x]
 * Specifies the maximum smaller mass between 0 and 5.0 $M\odot$.
 * \item[-m] 
 * Specifies the minimum mismatch threshold (typically 0.03) but for the
 * sake of testing it is best to pick a value $O[1]$ to save compiling time.
 * \item[-p]
 * Specifies that the program should generate a \MATHEMATICA notebook
 * ``Math3DNotebook.nb''.
 * \end{description}
 *
 * \subsubsection*{Exit codes}
 * \input{InspiralSpinBankTestCE}
 *
 * \subsubsection*{Notes}
 * \begin{itemize}
 * \item The metric used in InspiralSpinBank() is only valid for binary
 * systems with a total mass $<15M\odot$ where the minimum larger mass is at
 * least twice the maximum smaller mass.  Choosing mass values that violate
 * these conditions will cause an error message.
 * \item It is unlikely that you will be able to run a \MATHEMATICA notebook
 * that contains more than 10,000 tiles.  Adjust your parameters accordingly
 * if you plan to view a plot.
 * \end{itemize}
 *
 * \vfill{\footnotesize\input{InspiralSpinBankTestCV}}
 * </lalLaTeX> */


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
#include <lal/LALMathematica.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>


extern char *optarg;


/* <lalErrTable file="InspiralSpinBankTestCE"> */
#define INSPIRALSPINBANKTESTC_ENORM     0
#define INSPIRALSPINBANKTESTC_EMEM      1
#define INSPIRALSPINBANKTESTC_ESUB      2

#define INSPIRALSPINBANKTESTC_MSGENORM  "Normal exit"
#define INSPIRALSPINBANKTESTC_MSGEMEM   "Memory allocation error"
#define INSPIRALSPINBANKTESTC_MSGESUB   "Subroutine error"
/* </lalErrTable> */

NRCSID(LALINSPIRALSPINBANKTESTC, "$Id$");

int main( int argc, char *argv[] )
{
  static LALStatus stat;
  INT4 loop = 0; /* loop counter */
  Math3DPointList *list = NULL;    /* structure for mathematica plot */
  Math3DPointList *first = NULL;
  SnglInspiralTable *tiles = NULL;
  SnglInspiralTable *tmplt = NULL;
  InspiralCoarseBankIn coarseIn;
  INT4 ntiles = 0;                 /* number of tiles */
  INT2 Math3DPlot = 0;             /* option flag for Mathematica plot */
  INT4 opt = 0;                    /* returning value of getopt() */
  INT4 optflag = -1;               /* Command Line option */
  REAL4 PtSize = 0.02;
  REAL8Vector *psd = NULL;
  REAL8 df = 1.0;
  InspiralTemplate inspiralTemplate;
  INT2 printMoments = 0;
  InspiralMomentsEtc moments;
  REAL4 F0 = 0;
  REAL4 noiseMin = 1;
 
  if( (list = (Math3DPointList *) LALCalloc( 1, sizeof( Math3DPointList )))
      == NULL )
  {
    LALError( &stat, INSPIRALSPINBANKTESTC_MSGEMEM );
    printf( INSPIRALSPINBANKTESTC_MSGEMEM );
    return INSPIRALSPINBANKTESTC_EMEM;
  }

  first = list;

  /* Stuff for calculating the PSD and Noise Moments */
  coarseIn.fLower = inspiralTemplate.fLower = 30;
  coarseIn.fUpper = inspiralTemplate.fCutoff = 2000;
  coarseIn.iflso = 0;

 /* Parse options. */
  do
  {
    optflag++;  
    switch (opt)
    {
      case 'm':
        coarseIn.mmCoarse = atof( optarg );       
        break;
      case 'n':
        coarseIn.mMin = atof( optarg );
        break;
      case 'p':
        Math3DPlot = 1;
        break;
      case 's':
        printMoments = 1;
        break;
      case 'x':
        coarseIn.MMax = atof( optarg );
        break;
      default:
        coarseIn.mmCoarse = 0.1;
        coarseIn.mMin = 1.0;
        coarseIn.MMax = 2*3.0;
        Math3DPlot = 0;
        printMoments = 0;
        break;
    }
  }
  while( (opt = getopt( argc, argv, "n:m:x:ps" )) != -1 );
  
  coarseIn.shf.data = NULL;
  memset( &(coarseIn.shf), 0, sizeof( REAL8FrequencySeries ) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &stat, &psd, coarseIn.fUpper ); 
  df = 1.0;
  LALNoiseSpectralDensity( &stat, psd, &LALLIGOIPsd, df );
  coarseIn.shf.data = psd;
  coarseIn.shf.deltaF = df;
  
  for( loop = 0; loop < psd->length; loop++ )
  {
    if( psd->data[loop] > 0 && psd->data[loop] < noiseMin )
    { 
      F0 = (REAL4) coarseIn.shf.deltaF * loop;
      noiseMin = psd->data[loop];
    }
  }

  LALInspiralSpinBank( &stat, &tiles, &ntiles, coarseIn );
  REPORTSTATUS(&stat);
  if( stat.statusCode )
  {
    LALError( &stat, INSPIRALSPINBANKTESTC_MSGESUB );
    printf( INSPIRALSPINBANKTESTC_MSGESUB );
    return INSPIRALSPINBANKTESTC_ESUB;
  }

  /* Mathematica Plot Stuff */
  if( Math3DPlot )
  {
    for( tmplt = tiles; tmplt != NULL; tmplt = tmplt->next )
    {
      list->x = tmplt->psi0;
      list->y = tmplt->psi3;
      list->z = tmplt->beta;
      list->grayLevel = 1.0*loop/ntiles;
      if( (list = list->next = (Math3DPointList *) LALCalloc( 1, sizeof(
          Math3DPointList ))) == NULL ) {
        LALError( &stat, INSPIRALSPINBANKTESTC_MSGEMEM );
        printf( INSPIRALSPINBANKTESTC_MSGEMEM );
        return INSPIRALSPINBANKTESTC_EMEM;
      }
    }
    list->next = NULL;
    LALMath3DPlot( &stat, first, &ntiles, NULL );
    if( stat.statusCode )
    {
      LALError( &stat, INSPIRALSPINBANKTESTC_MSGESUB );
      printf( INSPIRALSPINBANKTESTC_MSGESUB );
      return INSPIRALSPINBANKTESTC_ESUB;
    }

    /* Clean Up the memory from the MathPlot3D structure */
    list = first;
    while( list->next )
    {
      first = list->next;
      LALFree( list );
      list = first;
    }
  }
  
  /* free the last (first?) memory allocated for Math3DPlot. */
  LALFree( list ); 
 
  if (stat.statusCode)
    return INSPIRALSPINBANKTESTC_ESUB;
  else
    return INSPIRALSPINBANKTESTC_ENORM;
}
