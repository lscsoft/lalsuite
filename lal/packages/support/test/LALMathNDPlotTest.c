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

/* ______________________________________________________________________________________
 * File name: LALMathNDPlotTest.c
 *
 * Author: Hanna C. R.
 *
 * Revision: $Id$
 * ______________________________________________________________________________________
 */

/* ------------------------------------ AUTO-DOC ------------------------------------- */
/* ----------------------------------------------------------------------------------- */

/*<lalVerbatim file="LALMathNDPlotTestCV">
  Author: Hanna, C.R.
  $Id$
  </lalVerbatim>*/

/*SUBSECTION - PROGRAM - "LALMathNDPlotTest.c" ------------------------------- <lalLaTeX>
  \subsection{Program \texttt{LALMathNDPlotTest.c}}
  \label{s:LALMathNDPlotTest.c}
  \providecommand{\MATHEMATICA}{$M\scriptstyle{ATHEMATICA}^{\textrm{{\small\textregistered} }}$}
* Tests LALMathNDPlot().
  </lalLaTeX>*/

  /*SUBSUBSECTION - USAGE - "LALMathNDPlotTest.c" ---------------------------- <lalLaTeX>
    \begin{verbatim}
    LALMathNDPlotTest
    \end{verbatim}
    </lalLaTeX>
    END SUBSUBSECTION - USAGE - "LALMathNDPlotTest.c" -------------------------------- */

  /*SUBSUBSECTION - DESCRIPTION - "LALMathNDPlotTest.c" ---------------------- <lalLaTeX>
    \subsubsection{Description}
  * This program generates a set of points simulating a 4-D template bank and calls
  * LALMathNDPlot() to generate a \MATHEMATICA notebook to display the permutations of
  * 3D projections of the template bank.  Instructions on how to evaluate the notebook
  * appear when it is opened.
    </lalLaTeX>
    END SUBSUBSECTION - DESCRIPTION - "LALMathNDPlotTest.c" -------------------------- */

  /*SUBSUBSECTION - EXIT CODES ------------------------------------------------------- */
    /* <lalLaTeX>
    \subsubsection*{Exit codes}
    \input{LALMathNDPlotTestCE}
    </lalLaTeX>
    END - SUBSUBSECTION - EXIT CODES ------------------------------------------------- */

  /*SUBSUBSECTION - NOTES - "LALMathNDPlotTest.c" ------------------------- <lalLaTeX>
    \subsubsection{Notes}
    \begin{itemize}
  * \item No notes yet.
    \end{itemize}
    </lalLaTeX>
    END - SUBSUBSECTION - NOTES - "LALMathNDPlotTest.c" --------------------------- */

  /*<lalLaTeX>
  \vfill{\footnotesize\input{LALMathNDPlotTestCV}}
  </lalLaTeX>
  END SUBSECTION - PROGRAM - "LALMathNDPlotTest.c" -------------------------------- */

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
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LALMathematica.h>

/*<lalErrTable file="LALMathNDPlotTestCE">*/
#define LALMATHNDPLOTTESTC_ENORM        0
#define LALMATHNDPLOTTESTC_EMEM         1
#define LALMATHNDPLOTTESTC_ESUB         2


#define LALMATHNDPLOTTESTC_MSGENORM     "Normal exit"
#define LALMATHNDPLOTTESTC_MSGEMEM      "Memory allocation error"
#define LALMATHNDPLOTTESTC_MSGESUB      "Subroutine error"
/*</lalErrTable>*/

NRCSID(LALMATHNDPLOTTESTC, "$Id:");

int main(void){
  static LALStatus status;
  INT4 loopx = 0;                       /* loop counters */
  INT4 loopy = 0;
  INT4 loopz = 0;
  INT4 loopw = 0;
  INT4 ntiles = 0;
  MathNDPointList *list = NULL;         /* Pointer to structure for mathematica plot */
  MathNDPointList *first = NULL;
  UINT4 dim = 4;

  if ((list = (MathNDPointList *) LALCalloc(1, sizeof(MathNDPointList))) == NULL){
    LALError(&status, LALMATHNDPLOTTESTC_MSGEMEM);
    printf(LALMATHNDPLOTTESTC_MSGEMEM);
    return LALMATHNDPLOTTESTC_EMEM;
  }
  first=list;

  for(loopx=1; loopx <= 20; loopx++){
    for(loopy=1; loopy <= 20; loopy++){
      for(loopz=0; loopz <= 1; loopz++){
        for(loopw=0; loopw <= 1; loopw++){
          LALSCreateVector(&status, &list->coordinates, dim);
          list->coordinates->data[0] = loopx;
          list->coordinates->data[1] = loopy;
          list->coordinates->data[2] = loopz;
          list->coordinates->data[3] = loopw;
          list->grayLevel = 0.0;
          ntiles++;
          if ((list = list->next = (MathNDPointList *) LALCalloc(1, sizeof(MathNDPointList))) == NULL){
            LALError(&status, LALMATHNDPLOTTESTC_MSGEMEM);
            printf(LALMATHNDPLOTTESTC_MSGEMEM);
            return LALMATHNDPLOTTESTC_EMEM;
            }
          }
        }
      }
    }

  for(loopx=1; loopx <= 20; loopx++){
    for(loopy=1; loopy <= 20; loopy++){
      for(loopw=0; loopw <= 1; loopw++){
        LALSCreateVector(&status, &list->coordinates, dim);
        list->coordinates->data[0] = loopx;
        list->coordinates->data[1] = loopy;
        list->coordinates->data[2] = 2;
        list->coordinates->data[3] = loopw;
        list->grayLevel = 1.0;
        ntiles++;
        if ((list = list->next = (MathNDPointList *) LALCalloc(1, sizeof(MathNDPointList))) == NULL){
          LALError(&status, LALMATHNDPLOTTESTC_MSGEMEM);
          printf(LALMATHNDPLOTTESTC_MSGEMEM);
          return LALMATHNDPLOTTESTC_EMEM;
          }
        }
      }
    }

  /* LAL!!! */
  for(loopx=1; loopx <= 20; loopx++){
    for(loopy=1; loopy <= 20; loopy++){
      for(loopz=3; loopz <= 4; loopz++){
        for(loopw=0; loopw <=1; loopw++){
          if( ((loopx==6)||(loopx==19)) && (loopy<16) && (loopy>5)) continue;
          if((loopy==15)&&(((loopx<20)&&(loopx>14))||((loopx<7)&&(loopx>1)))) continue;
          if((loopx>9)&&(loopx<12)&&(((loopy>6)&&(loopy<10))||(loopy==12))) continue;
          if(((loopx==9)||(loopx==12)) && ((loopy>9)&&(loopy<13))) continue;
          if(((loopx==8)||(loopx==13)) && ((loopy>12)&&(loopy<16))) continue;
          LALSCreateVector(&status, &list->coordinates, dim);
          if (status.statusCode){
            LALError(&status, LALMATHNDPLOTTESTC_MSGESUB);
            printf(LALMATHNDPLOTTESTC_MSGESUB);
            return LALMATHNDPLOTTESTC_ESUB;
            }
          list->coordinates->data[0] = loopx;
          list->coordinates->data[1] = loopy;
          list->coordinates->data[2] = loopz;
          list->coordinates->data[3] = loopw;
          list->grayLevel = 0.0;
          ntiles++;
          if ((list = list->next = (MathNDPointList *) LALCalloc(1, sizeof(MathNDPointList))) == NULL){
            LALError(&status, LALMATHNDPLOTTESTC_MSGEMEM);
            printf(LALMATHNDPLOTTESTC_MSGEMEM);
            return LALMATHNDPLOTTESTC_EMEM;
            }
          }
        }
      }
    }


  list->next = NULL;
  printf("\nCalling LALMathNDPlot()......\n");
  LALMathNDPlot(&status, first, &ntiles, NULL);
  REPORTSTATUS(&status);
  if (status.statusCode){
    LALError(&status, LALMATHNDPLOTTESTC_MSGESUB);
    printf(LALMATHNDPLOTTESTC_MSGESUB);
    return LALMATHNDPLOTTESTC_ESUB;
  }


  /* Clean Up the memory from the MathPlot3D structure */
  list = first;
  while(list->next){
    first = list->next;
    if (list->coordinates)
      LALSDestroyVector(&status, &list->coordinates);
    if (status.statusCode){
      LALError(&status, LALMATHNDPLOTTESTC_MSGESUB);
      printf(LALMATHNDPLOTTESTC_MSGESUB);
      return LALMATHNDPLOTTESTC_ESUB;
      }
    LALFree(list);
    if (status.statusCode){
      LALError(&status, LALMATHNDPLOTTESTC_MSGEMEM);
      printf(LALMATHNDPLOTTESTC_MSGEMEM);
      return LALMATHNDPLOTTESTC_EMEM;
      }
    list = first;
  }


  /* free the last (first?) memory allocated for Math3DPlot. */
  if(list) LALFree(list);

  if (status.statusCode)
    return LALMATHNDPLOTTESTC_ESUB;
  else
    return LALMATHNDPLOTTESTC_ENORM;

}

