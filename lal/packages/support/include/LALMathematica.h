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
  \providecommand{\MATHEMATICA}{$M\scriptstyle{ATHEMATICA}^{\textrm{{\small\textregistered} }}$}
* Provides structures, functions and macro definitions for modules
* that generate \MATHEMATICA notebooks. Currently, the only module using this header file 
* is \texttt{Math3DPlot.c}, which generates 3D animated plots of template banks having
* three parameters.  See the \texttt{Math3DPlot.c} documenation for more detail.  
  </lalLaTeX> */
 
  /*SUBSECTION - SYNOPSIS ---------------------------------------------------- <lalLaTeX>
    \subsection*{Synopsis}
    \begin{verbatim}
    #include <lal/LALMathematica.h>
    \end{verbatim}
  * This header file defines macros containing \MATHEMATICA syntax that is otherwise messy 
  * to implement into C source files.  With these macros one can easily write a function 
  * that will output a \MATHEMATICA notebook. Simply\ldots
    \begin{enumerate}
  * \item Open a file with a pointer named ``nb" and a file extenstion ``.nb".  
  * \item Use BEG\_NOTEBOOK to start the notebook file.  
  * \item Use the appropriate BEG and END macros with fprint(nb, ``Your Text") in 
  * between to write your text to the cells of the notebook.  If you are writing 
  * \MATHEMATICA commands use the INPUT macros; for plain text, use TEXT Macros.  
  * \item Denote titles and sections with the appropriate macros.  
  * \item Use END\_NOTEBOOK to end the notebook and 
  * use fclose(nb) to close the file ``nb". 
    \end{enumerate}
  * The result is very readable/changeable source similar in style to most markup 
  * languages. An example program might look like:
    \begin{verbatim}
  FILE *nb;
  nb = fopen("YourFileName.nb", "rw");
  BEG_NOTEBOOK;
  BEG_TITLECELL;
   fprintf(nb, "Sample Program Title");
  END_TITLECELL_;
  BEG_SECTIONCELL;
    fprintf(nb, "Sample Program Section Name");
  END_SECTIONCELL;
  END_NOTEBOOK;
  fclose(nb);
  \end{verbatim}
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
  * to generate a \MATHEMATICA notebook in your own program.  
    \begin{itemize}
  * \item{NOTEBOOK} Denotes the beginning and ending of the notebook file.  A 
  * BEG\_NOTEBOOK tag must start the file and an END\_NOTEBOOK tag must end it.
  * \item{TITLE} Placing an fprint(nb, "Your Title") inbetween these BEG and END tags 
  * will place a \emph{title font} cell in the notebook.
  * \item{GROUP} Cells placed in between these tags will be grouped together
  * \item{SECTION} Same as title except the text printed will be in \emph{section font}.
  * Subsequent input and text cells following the END\_SECTIONCELL tag will be grouped
  * with that section until a new BEG\_SECTIONCELL tag is encountered.
  * \item{INPUT} provides cells to input \MATHEMATICA  commands.
  * \item{TEXT} provides cells to input plain text.
    \end{itemize}
  * Notice that the file 
  * pointer must be named ``nb" in order to use the macros defined in this header.  
  * When grouping several cell objects together the last object in the list should have 
  * an underscored tag instead of an END tag without an underscore.  
  * Although the notebook will compile (usually) if you use the tags without an ending 
  * underscore, the dangling comma is taken as a null member of the list of grouped cells  
  * and when you view the notebook in \MATHEMATICA  you may see the word 
  * ``NULL" printed on a line.  That is an indication that you should use the underscore 
  * version of the tag which preceeded the ``NULL" statement.   
    </lalLaTeX> 
    END SUBSECTION - MACROS ---------------------------------------------------------- */

  /*<lalLaTeX> 
    \vfill{\footnotesize\input{LALMathematicaHV}} 
  </lalLaTeX> */

  /*NEWPAGE - SUBSECTION - TYPES ----------------------------------------------<lalLaTeX>
    \newpage
    \subsection*{Types} 
    \input{LALMathematicaHT}
    \idx[Type]{Math3DPointList}
    \noindent
  * The Math3DPointList type is used by \texttt{Math3DPlot.c} as an input structure to 
  * plot 3-dimensional template banks.  It is a linked list and has parameters for each 
  * coordinate x,y,z and the next pointer.  It also has a parameter called GrayLevel 
  * which must be $\epsilon [0,1]$.  It specifies the shading of the point in the final 
  * plot with 0 representing black and 1 representing white.  By creatively assigning 
  * its value the grayscale shade of the points may convey additional information.  
    </lalLaTeX> 
    END SUBSECTION TYPES --------------------------------------------------------------*/

  /*SUBSECTION - NOTES --------------------------------------------------------<lalLaTeX>
    \subsection*{Notes}
    \noindent\begin{itemize}
  * \item Obviously the definitions and functions associated with 
  * this header are NOT lal complient and thus do not belong in any lal routines except 
  * test programs.
  * \item There are many more commands to manipulate \MATHEMATICA notebooks that are not
  * included in this header.  The macros are only what is necessary for a 
  * \emph{bare minimum} interface.   
    \end{itemize} 
    </lalLaTeX>
    END SUBSECTION NOTES --------------------------------------------------------------*/

  /*<lalLaTeX> 
    \vfill{\footnotesize\input{LALMathematicaHV}} 
  </lalLaTeX> */

  /*NEWPAGE - MODULES INPUT - "LALMath3DPlot.c" ---------------------------------------*/
    /* <lalLaTeX> 
    \newpage\input{LALMath3DPlotC}
    </lalLaTeX> 
    END - MODULES INPUT - "LALMath3DPlot.c" -------------------------------------------*/ 
  
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
