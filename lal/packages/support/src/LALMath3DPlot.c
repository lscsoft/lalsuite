/*_______________________________________________________________________________________
 * 
 * File Name: LALMath3DPlot.c
 *
 * Author: Hanna C. R.
 * 
 * Revision: $Id$
 * 
 *_______________________________________________________________________________________
 */


/* ------------------------------------- AUTO-DOC ------------------------------------ */
/* ----------------------------------------------------------------------------------- */


/*<lalVerbatim file="LALMath3DPlotCV">
Author: Hanna, C. R.
$Id$
  </lalVerbatim>*/ 

/*SUBSECTION - MODULE - "LALMath3DPlot.c" ------------------------------------ <lalLaTeX>
  \subsection{Module \texttt{LALMath3DPlot.c}}
  \label{ss:LALMath3DPlot}
  </lalLaTeX> */
 
  /* SUBSUBSECTION - PROTOTYPES - "LALMath3DPlot()" -------------------------- <lalLaTeX>
     \subsubsection{Prototypes}
     \input{LALMath3DPlotCP}
     \idx{LALMath3DPlot()}
     </lalLaTeX>
     END SUBSUBSECTION - PROTOTYPES "LALMath3DPlot()" -------------------------------- */

  /* SUBSUBSECTION - DESCRIPTION --------------------------------------------- <lalLaTeX>
     \subsubsection{Description}
   * This module contains a function for plotting 3D template banks by creating a 
   * \MATHEMATICA notebook.  The notebook renders the templates as points in a three 
   * dimensional lattice.  The plot is animated so the user can see the template bank from 
   * different perspectives.  See figure 1.1.
     </lalLaTeX> 
     END SUBSUBSECTION - DESCRIPTION ------------------------------------------------- */

  /* FIGURE - "LALMathematicaHplot1" ----------------------------------------- <lalLaTeX>
     \begin{figure}
     \begin{center}
     \includegraphics[width=0.5\textwidth]{LALMathematicaHplot1}
     \end{center}
     \caption{Here is an example template bank produced by running InspiralSpinBankTest.c
              to generate roughly 5000 templates.  Currently the plot doesn't show the 
              contour of the templates; it renders them as spheres.  In the case of 
              metrics with disimilar scales along the principle directions you will 
 	      notice considerable space between points accordingly.}
     \end{figure}
     </lalLaTeX>
     END - FIGURE -------------------------------------------------------------------- */

  /* SUBSUBSECTION - NOTES --------------------------------------------------- <lalLaTeX>
     \subsubsection{Notes}
     \begin{itemize}
   * \item The output of this function is ``Math3DNotebook.nb" and will appear in the 
   * directory of the program that called this function.
   * \item Exported \MATHEMATICA graphics have no preferred directory and will appear in 
   * your home directory for unix users and in the $\backslash$Mathematica directory for 
   * Windows users unless you have another path configured in your \MATHEMATICA 
   * installation. It is necessary to change the file name within the notebook to avoid 
   * overwriting previous files.
     \end{itemize}
     </lalLaTeX>
     END SUBSUBSECTION - NOTES ------------------------------------------------------- */
     
/*END - SUBSECTION - MODULE - LALMath3DPlot.c" --------------------------------------- */

/*<lalLaTeX>  
\vfill{\footnotesize\input{LALMath3DPlotCV}}
</lalLaTeX>*/


/* -------------------------------------END AUTO DOC --------------------------------- */
/* ----------------------------------------------------------------------------------- */


#include <lal/LALConfig.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LALMathematica.h>

#define INSTRUCTIONS 	fprintf(nb, "This notebook will produce an animated 3D plot of your template bank.  See the next section to change any user variables before evaluating.  The cells of this notebook must be evaluated sequentially.  If you wish to evaluate the entire notebook at once press Ctrl+A then press Shift+Enter in most operating systems.")


NRCSID(LALMATH3DPLOTC, "$Id$");

/* <lalVerbatim file="LALMath3DPlotCP"> */
void
LALMath3DPlot( LALStatus *stat, 
               Math3DPointList *first, 
               INT4 *ntiles) 
/* </lalVerbatim>*/
{
  FILE *nb; 				/* pointer to the notebook file */
  INT4 jflag = 0;			/* flag to justify the output data */
  Math3DPointList *list;		/* loop counter */
  INT4 loop = 0;
  
  INITSTATUS( stat, "LALMath3DPlot", LALMATH3DPLOTC ); 

  if (!first)
    ABORT(stat, LALMATHEMATICAH_ENULL, LALMATHEMATICAH_MSGENULL);

  if ((nb = fopen("Math3DNotebook.nb", "w")) == NULL)
    ABORT(stat, LALMATHEMATICAH_EFILE, LALMATHEMATICAH_MSGEFILE);
  

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
        fprintf(nb, "AnimationSize\t= {400,400};");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "StillSize\t= {600,600};");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "AnimationName\t:= \"AnimationTilePlot.gif\"");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "StillName\t:= \"StillTilePlot\"");
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
      BEG_INPUTCELL;
        fprintf(nb, "XAxisLabel := \"Psi01\"");
      END_INPUTCELL;
      BEG_INPUTCELL; 
        fprintf(nb, "YAxisLabel := \"Psi03\"");
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "ZAxisLabel := \"Beta\"");
      END_INPUTCELL;
      BEG_TEXTCELL;
        fprintf(nb, "AnimationSize:\tThe size of the final animation in PIXELS x PIXELS\n");
        fprintf(nb, "StillSize:\t\tThe size of the final still image in PIXELS x PIXELS\n");
        fprintf(nb, "AnimationName:\tWhat to name the final animation - must use .gif extension\n");
        fprintf(nb, "StillName:\t\tWhat to name the final still image - extension determined by StillType\n");
        fprintf(nb, "StillType:\t\tThe file type and extension for the still image\n");
        fprintf(nb, "\t\t\tChoose any standard format (e.g. JPG, GIF, PDF, EPS, etc.)\n");
        fprintf(nb, "frames:\t\tThe number of frames for each rotation of the image.\n"); 
        fprintf(nb, "\t\t\tThe final image will have 2 times the number frames\n");
        fprintf(nb, "FrameTime:\t\tSets the delay time in seconds between each frame in the animated gif\n");
        fprintf(nb, "\t\t\tApplications seem to interpret this differently.  You may have to adjust this setting\n");
        fprintf(nb, "\t\t\tbased on the intended application of the animations.\n");
        fprintf(nb, "XAxisLabel:\t\tSets the X-axis label\n");
        fprintf(nb, "XAxisLabel:\t\tSets the Y-axis label\n");
        fprintf(nb, "XAxisLabel:\t\tSets the Z-axis label");
      END_TEXTCELL_;
    END_GROUPCELLC;
    BEG_GROUPCELL;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Point List");
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "TILES  = \n"); 
          fprintf(nb, "Graphics3D[{PointSize[0.02],"); 
          list = first;
          while(list->next) 
          {
            fprintf(nb, "{GrayLevel[%f], Point[{%f,%f,%f}]},",
                    list->GrayLevel, list->x, list->y, list->z);
            if (jflag%2) fprintf(nb,"\n");
            ++jflag;
            list = list->next;
          }
          fprintf(nb, "{GrayLevel[%f], Point[{%f, %f, %f}]}}]", list->GrayLevel, list->x, list->y, list->z);
        END_INPUTCELL_;
      END_GROUPCELLC;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Image generation");
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "still = Show[TILES, Background-> RGBColor[.93, .91, .89], ViewPoint -> {1, 1.3, 2.4}, ");
          fprintf(nb, "ImageSize->StillSize, Axes->True, AxesLabel->{XAxisLabel, YAxisLabel, ZAxisLabel}];\n");
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "Do[tile[T]=Show[TILES, Background -> RGBColor[.93, .91, .89], ");
          fprintf(nb, "ViewPoint -> {1-(.99 T/frames)^2, T/(4 frames), (T/frames)^2},ImageSize->AnimationSize], {T, 0, frames, 1}];\n");
          fprintf(nb, "Do[tile[frames+T]=Show[TILES, Background -> RGBColor[.93, .91, .89], ViewPoint -> {.005+(T/frames)^2, ");
          fprintf(nb, "0.25-T/(4 frames), 1-(.99 T/frames)^2},ImageSize->AnimationSize], {T, 0, frames, 1}];\n");
        END_INPUTCELL_;  
      END_GROUPCELLC;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Animation Generation");
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "images = Evaluate[Table[tile[j], {j, 0, 2 frames, 1}]];\n");
        END_INPUTCELL;
        BEG_INPUTCELL; 
          fprintf(nb, "Export[StillName<>\".\"<>ToLowerCase[StillType], still, StillType, ImageSize->StillSize, ");
          fprintf(nb, "ConversionOptions->{\"ColorReductionDither\" -> False}]");
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "Export[AnimationName, images, \"GIF\", ImageSize -> AnimationSize, ");
          fprintf(nb, "ConversionOptions -> {\"Loop\" -> True,\"AnimationDisplayTime\" -> FrameTime, ");
          fprintf(nb, "\"ColorReductionDither\" -> False}]");
        END_INPUTCELL_;   
      END_GROUPCELLC_;
    END_GROUPCELLC_;    
  END_NOTEBOOK;
  fclose(nb);
  RETURN(stat);
}

