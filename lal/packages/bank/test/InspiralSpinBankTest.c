/*_______________________________________________________________________________________
 * 
 * File Name: InspiralSpinBankTest.c
 *
 * Author: Hanna C. R.
 * 
 * Revision: $Id$
 * 
 *_______________________________________________________________________________________
 */


/* ------------------------------------ AUTO-DOC ------------------------------------- */
/* ----------------------------------------------------------------------------------- */

/*<lalVerbatim file="InspiralSpinBankTestCV">
  Author: Hanna, C.R.
  $Id$
  </lalVerbatim>*/ 

/*SUBSECTION - PROGRAM - "InspiralSpinBankTest.c" ---------------------------- <lalLaTeX>
  \subsection{Program \texttt{InspiralSpinBankTest.c}}
  \label{s:InspiralSpinBankTest.c}
  \providecommand{\MATHEMATICA}{$M\scriptstyle{ATHEMATICA}^{\textrm{{\small\textregistered} }}$}
* Tests InpsiralSpinBank(). 
  </lalLaTeX>*/
  
  /*SUBSUBSECTION - USAGE - "InspiralSpinBankTest.c" ------------------------- <lalLaTeX>
    \begin{verbatim}
    InspiralSpinBankTest
    \end{verbatim} 
    </lalLaTeX>
    END SUBSUBSECTION - USAGE - "InspiralSpinBankTest.c" ----------------------------- */

  /*SUBSUBSECTION - DESCRIPTION - "InspiralSpinBankTest.c" ------------------- <lalLaTeX>
    \subsubsection{Description}
  * This program uses InspiralSpinBank() to generate a template bank from command line
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
    \input{InspiralSpinBankTestCE}
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
#include <lal/LALConfig.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LALMathematica.h>

/*<lalErrTable file="InspiralSpinBankTestCE">*/
#define INSPIRALSPINBANKTESTC_ENORM     0
#define INSPIRALSPINBANKTESTC_EMEM      1
#define INSPIRALSPINBANKTESTC_ESUB      2
                                                                                                                                                              
#define INSPIRALSPINBANKTESTC_MSGENORM  "Normal exit"
#define INSPIRALSPINBANKTESTC_MSGEMEM   "Memory allocation error"
#define INSPIRALSPINBANKTESTC_MSGESUB   "Subroutine error"
/*</lalErrTable>*/

NRCSID(LALINSPIRALSPINBANKTESTC, "$Id$");

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

  if ((list = (Math3DPointList *) LALCalloc(1, sizeof(Math3DPointList))) == NULL){
    LALError(&stat, INSPIRALSPINBANKTESTC_MSGEMEM);
    printf(INSPIRALSPINBANKTESTC_MSGEMEM);
    return INSPIRALSPINBANKTESTC_EMEM;
    }

  first = list;
 
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
      default:
        CoarseIn.mmCoarse = 1.0;
        CoarseIn.mMin = 1.0;
        CoarseIn.MMax = 3.0;
        Math3DPlot = 0;
        break;
      }
    } while ((opt = getopt( argc, argv, "n:m:x:p" )) != -1);

  printf("\n\n\n____________________________________________________________________");
  if (optflag) printf("\nThe parameters you specified (along with defaults) are:\n");
  else {
    printf("\nNOTE: you can provide input parameters in the command line");
    printf("\n\nThe options are:\n");
    printf("\t-x argument = maximum mass of the smaller body in Solar Mass\n");
    printf("\t-n argument = minimum mass of the same\n");
    printf("\t-m argument = Mismatch threshold (e.g. 0.03 is typical");
    printf("\n\t\t      but takes a while to run)\n");
    printf("\t-p (no arg) - calling this option makes a mathematica plot");
    printf("\n\t\t      of the tiles");
    printf("\nThe default parameters are:\n");
    } 
  printf("\nMismatch = %f\nMathematica Plot = %i\nMass two min = %f\nMass two max = %f\n", 
    CoarseIn.mmCoarse, Math3DPlot, CoarseIn.mMin, CoarseIn.MMax);
  printf("______________________________________________________________________\n");
    
  printf("\n\nCalling LALInspiralSpinBank()......\n"); 
  LALInspiralSpinBank(&stat, &Tiles, &ntiles, CoarseIn);  
  REPORTSTATUS(&stat);
  if (stat.statusCode){
    LALError(&stat, INSPIRALSPINBANKTESTC_MSGESUB);
    printf(INSPIRALSPINBANKTESTC_MSGESUB);
    return INSPIRALSPINBANKTESTC_ESUB;
  }
    
  printf("\nThese parameters yeild %i tiles.\n\n\n", ntiles);

  /* Mathematica Plot Stuff */
  if (Math3DPlot){
    for(loop=0; loop < (ntiles); loop++){
      list->x = Tiles[loop].params.psi0;
      list->y = Tiles[loop].params.psi3;
      list->z = Tiles[loop].params.beta;
      list->GrayLevel = 1.0*loop/ntiles;
      if ((list = list->next = (Math3DPointList *) LALCalloc(1, sizeof(Math3DPointList))) == NULL){
        LALError(&stat, INSPIRALSPINBANKTESTC_MSGEMEM);
        printf(INSPIRALSPINBANKTESTC_MSGEMEM);
        return INSPIRALSPINBANKTESTC_EMEM;
      }
    }
    list->next = NULL;
    printf("\nCalling LALMath3DPlot()......\n");
    LALMath3DPlot(&stat, first, &ntiles, NULL);
    REPORTSTATUS(&stat);
    if (stat.statusCode){
      LALError(&stat, INSPIRALSPINBANKTESTC_MSGESUB);
      printf(INSPIRALSPINBANKTESTC_MSGESUB);
      return INSPIRALSPINBANKTESTC_ESUB;
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
    return INSPIRALSPINBANKTESTC_ESUB;
  else
    return INSPIRALSPINBANKTESTC_ENORM;

}

