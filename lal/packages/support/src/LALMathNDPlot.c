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
 * File Name: LALMathNDPlot.c
 *
 * Author: Hanna C. R.
 *
 * Revision: $Id$
 *
 *_______________________________________________________________________________________
 */

/* ------------------------------------- AUTO-DOC ------------------------------------ */
/* ----------------------------------------------------------------------------------- */


/*<lalVerbatim file="LALMathNDPlotCV">
  Author: Hanna, C. R.
  $Id$
  </lalVerbatim>*/

/*SUBSECTION - MODULE - "LALMathNDPlot.c" ------------------------------------ <lalLaTeX>
  \subsection{Module \texttt{LALMathNDPlot.c}}
  \label{ss:LALMathNDPlot}
  </lalLaTeX> */

  /* SUBSUBSECTION - PROTOTYPES - "LALMathNDPlot()" -------------------------- <lalLaTeX>
     \subsubsection{Prototypes}
     \input{LALMathNDPlotCP}
     \idx{LALMathNDPlot()}
     \noindent\texttt{*stat} LALStatus structure pointer\\
     \\\texttt{*first} MathNDPointList stucture pointer\\
     \\\texttt{*ntiles} INT4 pointer to the number of templates you \emph{plan} to plot.
     * This may be called as NULL.  If it is called with a value this function will check
     * to see if the MathNDPointList has the correct number of templates.  If it does not
     * a warning will be printed. \\
     \\\texttt{pointSize} REAL4 $\epsilon[0,1]$ which specifies the relative size of each
     * point to the final display area.  (e.g. 1 would fill the enire plot.)  This may be
     * called as NULL and a calculated value will be assigned.  (Its only a rough guess)
     </lalLaTeX>
     END SUBSUBSECTION - PROTOTYPES "LALMathNDPlot()" -------------------------------- */

  /* SUBSUBSECTION - DESCRIPTION --------------------------------------------- <lalLaTeX>
     \subsubsection{Description}
   * This module contains a function for plotting N-Dimensional template banks by creating
   * a \MATHEMATICA notebook.  The notebook renders the templates as points in all of the
   * 3-Dimensional projection permutations. Each projection may be animated so the user
   * can see the template bank from different perspectives.
     </lalLaTeX>
     END SUBSUBSECTION - DESCRIPTION ------------------------------------------------- */


  /* SUBSUBSECTION - NOTES --------------------------------------------------- <lalLaTeX>
     \subsubsection{Notes}
     \begin{itemize}
   * \item The output of this function is ``MathNDNotebook.nb" and will appear in the
   * directory of the program that called this function.
   * \item Exported \MATHEMATICA graphics  will appear in
   * your home directory for unix users and in the $\backslash$Mathematica directory for
   * Windows users unless you have another path configured in your \MATHEMATICA
   * installation. It is necessary to change the file name within the notebook to avoid
   * overwriting previous files.
   * \item The number of projections is N!/(3!(N-3)!).  Thus plotting 6 dimensions would
   * produce 20 projections, 7 dimensions would yeilds 35 amd 8 gives 56.
     \end{itemize}
     </lalLaTeX>
     END SUBSUBSECTION - NOTES ------------------------------------------------------- */

/*END - SUBSECTION - MODULE - LALMathNDPlot.c" --------------------------------------- */

/*<lalLaTeX>
\vfill{\footnotesize\input{LALMathNDPlotCV}}
</lalLaTeX>*/


/* -------------------------------------END AUTO DOC --------------------------------- */
/* ----------------------------------------------------------------------------------- */


#include <lal/LALConfig.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LALMathematica.h>

#define INSTRUCTIONS    fprintf(nb, "Running this entire notebook using ctrl+A and shift+enter may crash your computer.  Evaluate each section as needed.  The Initialization and User Variables sections must be evaluated first.  The 3-dimensional projections are represented in the sections below User Varibles as PointList (x1, x2, x3) etc.  Evaluating the entire Image Generation sections creates animated plots (if AnimationPlot := True).  If (AnimationPlot := False) you get only still plots, saving time and memory.")

NRCSID(LALMATHNDPLOTC, "$Id$");

/* <lalVerbatim file="LALMathNDPlotCP"> */
void
LALMathNDPlot( LALStatus *stat,
               MathNDPointList *first,
               INT4 *ntiles,
               REAL4 *pointSize)
/* </lalVerbatim>*/
{
  FILE *nb;                             /* pointer to the notebook file */
  INT4 jflag = 0;                       /* flag to justify the output data */
  MathNDPointList *list;                /* loop counter */
  REAL4 PtSize = 0.02;
  INT4 counter = 0;
  INT4 dim = first->coordinates->length;
  INT4 x = 0;
  INT4 y = 0;
  INT4 z = 0;


  INITSTATUS( stat, "LALMathNDPlot", LALMATHNDPLOTC );

  /* Check that the PointList isn't NULL */
  if (!first) {
    ABORT(stat, LALMATHEMATICAH_ENULL, LALMATHEMATICAH_MSGENULL);
  }

  /* Open a file for writing a notebook */
  if ((nb = fopen("MathNDNotebook.nb", "w")) == NULL) {
    ABORT(stat, LALMATHEMATICAH_EFILE, LALMATHEMATICAH_MSGEFILE);
  }

  /* Appropriately handle the inputs for ntiles and pointsize to assure
     that a propter pointSize is chosen.  Also print a warning if the
     length of the MathNDPointList is not equal in length to the parameter
     ntiles passed to this function. */
  if (!pointSize){
    if (!ntiles){
      list=first;
      while(list->next){
        counter++;
        list=list->next;
        }
      ntiles = &counter;
      }
    else{
      list=first;
      while(list->next){
        counter++;
        list=list->next;
        }
      if (*ntiles != counter)
        printf("\nWARNING!!! The value of argument ntiles (%i) != the MathNDPointList length (%i)\n",
               *ntiles, counter);
      }
    if (*ntiles <=0) {
      ABORT(stat, LALMATHEMATICAH_EVAL, LALMATHEMATICAH_MSGEVAL);
    }
    PtSize = 0.50*(1.0/(pow((*ntiles),0.333333)));
    if (*ntiles > 10000)
      printf("\nWARNING!!! More than 10,000 tiles may crash Mathematica:)\n");
  }

  else{
    if ((*pointSize <= 0.0) || (*pointSize >= 1.0)) {
      printf("\nIllegal value of pointSize; it must be between 0 and 1.\n");
      printf("The default value of 0.02 will be used");
      PtSize = 0.02;
    }
  }

  /* The code that generates the notebook */
  BEG_NOTEBOOK;
    BEG_TITLECELL;
      fprintf(nb, "LALMath3D Output");
    END_TITLECELL;
    BEG_SECTIONCELL;
      fprintf(nb, "Instructions");
    END_SECTIONCELL;
    BEG_TEXTCELL;
      INSTRUCTIONS;
    END_TEXTCELL;
    BEG_GROUPCELL;
      BEG_SECTIONCELL;
        fprintf(nb, "Initialization");
      END_SECTIONCELL;
      BEG_INPUTCELL;
        fprintf(nb, "Off[General::spell];");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "Off[General::spell1];");
      END_INPUTCELL_;
    END_GROUPCELLC;
    BEG_SECTIONCELL;
      fprintf(nb, "User Variables");
    END_SECTIONCELL;
    BEG_GROUPCELL;
      BEG_INPUTCELL;
        fprintf(nb, "AnimationPlot\t:= False;");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "AnimationSize\t= {400,400};");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "StillSize\t= {600,600};");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "PtSize\t= %f;", PtSize);
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "AnimationName\t:= \"AnimationTilePlot\""); /* dont forget to giv a unique name */
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "StillName\t:= \"StillTilePlot\""); /* dont forget to give a unique name for each section */
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "StillType\t:=\"EPS\"");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "frames\t= 30;");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "FrameTime\t= 0.2;");
      END_INPUTCELL;
      for(counter = 1; counter <= dim; counter++){
        BEG_INPUTCELL; /* Create an axis label for each dimension */
          fprintf(nb, "X%iAxisLabel := \"x%i\"", counter, counter);
        END_INPUTCELL;
        }
      BEG_TEXTCELL;
        fprintf(nb, "AnimationPlot:\tFlag to set for generating animations (may take a while to run)\n");
        fprintf(nb, "AnimationSize:\tThe size of the final animation in PIXELS x PIXELS\n");
        fprintf(nb, "StillSize:\t\tThe size of the final still image in PIXELS x PIXELS\n");
        fprintf(nb, "PtSize:\t\tThe relative size of the template points to the final display width.\n");
        fprintf(nb, "\t\t\tIt is given as a decimal part of one. (e.g. PtSize=0.02 is 1/20 of the display width)\n");
        fprintf(nb, "AnimationName:\tWhat to name the final animation.\n");
        fprintf(nb, "StillName:\t\tWhat to name the final still image - extension determined by StillType\n");
        fprintf(nb, "StillType:\t\tThe file type and extension for the still image\n");
        fprintf(nb, "\t\t\tChoose any standard format (e.g. JPG, GIF, PDF, EPS, etc.)\n");
        fprintf(nb, "frames:\t\tThe number of frames for each rotation of the image.\n");
        fprintf(nb, "\t\t\tThe final image will have 2 times the number frames\n");
        fprintf(nb, "FrameTime:\t\tSets the delay time in seconds between each frame in the animated gif\n");
        fprintf(nb, "\t\t\tApplications seem to interpret this differently.  You may have to adjust this setting\n");
        fprintf(nb, "\t\t\tbased on the intended application of the animations.\n");
        fprintf(nb, "XiAxisLabel:\t\tSets the Xi-axis label\n");
      END_TEXTCELL_;
    END_GROUPCELLC;
    /* begin iterations of projection permutations */
    for(x = 0; x < dim; x++){
      for(y = (x+1); y < dim; y++){
        for(z = (y+1); z < dim; z++){
    BEG_GROUPCELL;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Point List (x%i,x%i,x%i)", (x+1), (y+1), (z+1));
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "TILES%i%i%i  = \n", (x+1), (y+1), (z+1));
          fprintf(nb, "Graphics3D[{PointSize[PtSize]");
          list = first;
          while(list->next)
          {
            fprintf(nb, ",{GrayLevel[%f], Point[{%f,%f,%f}]}",
                    list->grayLevel, list->coordinates->data[x], list->coordinates->data[y], list->coordinates->data[z]);
            if (jflag%2) fprintf(nb,"\n");
            ++jflag;
            list = list->next;
          }
          fprintf(nb, "}]");
        END_INPUTCELL_;
      END_GROUPCELLC;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Image generation (x%i,x%i,x%i)", (x+1), (y+1), (z+1));
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "still%i%i%i = Show[TILES%i%i%i, Background-> RGBColor[.93, .91, .89], ViewPoint -> {1, 1.3, 2.4}, ",
                  (x+1), (y+1), (z+1), (x+1), (y+1), (z+1));
          fprintf(nb, "ImageSize->StillSize, Axes->True, AxesLabel->{X%iAxisLabel, X%iAxisLabel, X%iAxisLabel}];\n",
                  (x+1), (y+1), (z+1));
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "If[AnimationPlot,{Do[tile%i%i%i[T]=Show[TILES%i%i%i, Background -> RGBColor[.93, .91, .89], ",
                  (x+1), (y+1), (z+1), (x+1), (y+1), (z+1));
          fprintf(nb, "ViewPoint -> {1-(.99 T/frames)^2, T/(4 frames), 2 (T/frames)^2},ImageSize->AnimationSize], {T, 0, frames, 1}],\n");
          fprintf(nb, "Do[tile%i%i%i[frames+T]=Show[TILES%i%i%i, Background -> RGBColor[.93, .91, .89], ViewPoint -> {.005+(T/frames)^2, ",
                  (x+1), (y+1), (z+1), (x+1), (y+1), (z+1));
          fprintf(nb, "0.25-T/(4 frames), 2-2 (.99 T/frames)^2},ImageSize->AnimationSize], {T, 0, frames, 1}]}];\n");
        END_INPUTCELL_;
      END_GROUPCELLC;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Image Export (x%i,x%i,x%i)", (x+1), (y+1), (z+1));
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "If[AnimationPlot, images%i%i%i = Evaluate[Table[tile%i%i%i[j], {j, 0, 2 frames, 1}]]];\n",
                  (x+1), (y+1), (z+1), (x+1), (y+1), (z+1));
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "Export[StillName<>\"%i%i%i\"<>\".\"<>ToLowerCase[StillType], still%i%i%i, StillType, ImageSize->StillSize, ",
                  (x+1), (y+1), (z+1), (x+1), (y+1), (z+1));
          fprintf(nb, "ConversionOptions->{\"ColorReductionDither\" -> False}]");
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "If[AnimationPlot, Export[AnimationName<>\"%i%i%i\"<>\".gif\", images%i%i%i, \"GIF\", ImageSize -> AnimationSize, ",
                  (x+1), (y+1), (z+1), (x+1), (y+1), (z+1));
          fprintf(nb, "ConversionOptions -> {\"Loop\" -> True,\"AnimationDisplayTime\" -> FrameTime, ");
          fprintf(nb, "\"ColorReductionDither\" -> False}]]");
        END_INPUTCELL_;
      END_GROUPCELLC_;
    END_GROUPCELLC;
	  }
        }
      }
  END_NOTEBOOK;
  fclose(nb);
  RETURN(stat);
}/* END - LALMathNDPlot() */




