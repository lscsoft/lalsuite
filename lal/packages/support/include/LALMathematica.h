/*_______________________________________________________________________________________ 
 * 
 * File Name: LALMathematica.h
 *
 * Author: Hanna C. R.
 * 
 * Revision: $Id$
 * 
 *_______________________________________________________________________________________
 */

/* -------------------------------------- AUTO-DOC ----------------------------------- */
/* ----------------------------------------------------------------------------------- */


/*<lalVerbatim file="LALMathematicaHV">
  Author: Hanna, C. R.
  $Id$
  </lalVerbatim> */ 

/*SECTION - HEADER - "LalMathematica.h" -------------------------------------- <lalLaTeX>
  \section{Header \texttt{LALMathematica.h}}
  \label{s:LALMathematica.h}
* Provides structures, functions and macro definitions for modules
* that generate Mathematica notebooks.  Currently, the only module
* that takes advantage of this header file is \texttt{Math3DPlot.c}.
* \texttt{Math3DPlot.c} generates 3D animated plots of template banks
* with three parameters.  See the \texttt{Math3DPlot.c} documenation
* for more detail.  
  </lalLaTeX> */
 
  /*SUBSECTION - SYNOPSIS ---------------------------------------------------- <lalLaTeX>
    \subsection*{Synopsis}
    \begin{verbatim}
    #include <lal/LALMathematica.h>
    \end{verbatim}
  * This header file defines macros containing Mathematica syntax that 
  * is otherwise messy to implement into C source files.  The result is 
  * a very readable module that is similar in style to most markup 
  * languages (e.g. XML, HTML, etc).  One simply has to use the macros
  * as tags for inserting their plain text or Mathematica commands.
  * These insertions will appear in the generated notebook according to
  * the user's specification.  It makes subsequent editing of source
  * files very easy.  One doesn't need to worry excessively about 
  * trying to combine C and Mathematica syntax.  
    </lalLaTeX> 
    END SUBSECTION - SYNOPSIS -------------------------------------------------------- */

  /*<lalLaTeX> 
    \vfill{\footnotesize\input{LALMathematicaHV}} 
  </lalLaTeX> */

  /*NEWPAGE - SUBSECTION - ERROR CODES --------------------------------------- <lalLaTeX>
    \newpage
    \subsection*{Error codes} 
    \input{LALMathematicaHE}
    </lalLaTeX> 
    END SUBSECTION - ERROR CODES ----------------------------------------------------- */

  /*SUBSECTION - MACROS ------------------------------------------------------ <lalLaTeX>
    \subsection*{Macros} 
    \input{LALMathematicaHM}
    \noindent
  * See the source file \texttt{Math3DPlot.c} for an example of how to use these macros
  * to generate a Mathematica notebook in your own programs.  
  * When grouping several cell objects together the last object in the list should have 
  * an???????tag instead of the END tags without underscores.  The underscore
  * implies that it is the end of a list and syntactically is simply missing a comma 
  * delimiter.  When in doubt use the one without the underscore!  The notebook will 
  * compile most of the time if you always use the tags without an ending underscore.  
  * However, the notebook will certainly not compile if you misuse the ending underscore 
  * tag.  When you view the notebook in mathematica you may see the word "NULL" printed on 
  * a line.  That is an indication that you should use the underscore version of the tag 
  * which preceeded the "NULL" statement.  That should fix the problem. 
    </lalLaTeX> 
    END SUBSECTION - MACROS ---------------------------------------------------------- */

  /*SUBSECTION - TYPES ----------------------------------------------------------------*/
    /* <lalLaTeX> 
    \subsection*{Types} 
    \input{LALMathematicaHT}
    \idx[Type]{Math3DPointList}
    \noindent
  * The Math3DPointList type is used by \texttt{Math3DPlot.c} as an input structure to 
  * plot 3-dimensional template banks.  It is a linked list and has parameters for each 
  * coordinate x,y,z and the next pointer.  It also has a parameter called GrayLevel 
  * which must be of the set ?????.  It specifies the shading of the point in the final 
  * plot.  By creatively assigning its value the intensity of the points may convey 
  * additional information.  If for no other reason it usually makes the plot more 
  * readable. 
    </lalLaTeX> 
    END SUBSECTION TYPES --------------------------------------------------------------*/

  /*NEWPAGE - MODULES - "LALMath3DPlot.c" ---------------------------------------------*/
    /* <lalLaTeX> 
    \newpage\input{LALMath3DPlotC} 
    </lalLaTeX> 
    END - MODULES - "LALMath3DPlot.c" ------------------------------------------------ */
/*END SECTION - HEADER - "LalMathematica.h" -------------------------------------------*/


/* ----------------------------------------------------------------------------------- */
/* ---------------------------------- END AUTO-DOC ----------------------------------- */




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>

/* <lalErrTable file="LALMathematicaHE"> */
#define LALMATHEMATICAH_ENULL 1 
#define LALMATHEMATICAH_MSGENULL "NULL pointer to a LALMathematica.h input structure"
#define LALMATHEMATICAH_EFILE 2 
#define LALMATHEMATICAH_MSGEFILE "Could not open file to write a Mathematica Notebook"
/* </lalErrTable> */


/* <lalVerbatim file="LALMathematicaHM"> */
#define BEG_NOTEBOOK 		fprintf(nb, "Notebook[{\n")
#define END_NOTEBOOK		fprintf(nb, "}]\n")
#define BEG_TITLECELL		fprintf(nb, "Cell[\"")
#define END_TITLECELL		fprintf(nb, "\", \"Title\"],\n")
#define END_TITLECELL_ 		fprintf(nb, "\", \"Title\"]\n")
#define BEG_GROUPCELL  		fprintf(nb, "Cell[CellGroupData[{\n")
#define END_GROUPCELLC  	fprintf(nb, "}, Closed  ]],\n")
#define END_GROUPCELLC_ 	fprintf(nb, "}, Closed  ]]\n")
#define BEG_SECTIONCELL		fprintf(nb, "Cell[\"")
#define END_SECTIONCELL		fprintf(nb, "\", \"Section\"],\n")
#define END_SECTIONCELL_	fprintf(nb, "\", \"Section\"]\n")
#define BEG_INPUTCELL		fprintf(nb, "Cell[BoxData[\\(")
#define END_INPUTCELL   	fprintf(nb, "\\)], \"Input\"],\n")
#define END_INPUTCELL_  	fprintf(nb, "\\)], \"Input\"]\n")
#define BEG_TEXTCELL		fprintf(nb, "Cell[\"\\<")
#define END_TEXTCELL		fprintf(nb, "\\>\", \"Text\"],\n")
#define END_TEXTCELL_		fprintf(nb, "\\>\", \"Text\"]\n")
/* </lalVerbatim> */

NRCSID (LALMATHEMATICAH, "$Id$");

/* <lalVerbatim file="LALMathematicaHT"> */
typedef struct Math3DPointList{
  struct Math3DPointList *next;
  REAL4 x;			
  REAL4 y;			
  REAL4 z;  		
  REAL4 GrayLevel;		
  }Math3DPointList;  
/* </lalVerbatim> */

void  
LALMath3DPlot( LALStatus *stat, 
               Math3DPointList *first, 
               INT4 *ntiles );
