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

#include <lal/LALConfig.h>
#include <lal/LALMalloc.h>
#include <lal/LALMathematica.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>

#define INSTRUCTIONS    fprintf(nb, "This notebook will produce an animated 3D plot of your template bank.  See the next section to change any user variables before evaluating.  The cells of this notebook must be evaluated sequentially.  If you wish to evaluate the entire notebook at once press Ctrl+A then press Shift+Enter in most operating systems.")

/**
 * \brief This function is for plotting 3D template banks by creating a MATHEMATICA notebook.
 * \ingroup LALMathematica_h
 * \author Hanna, C. R.
 *
 * The notebook renders the templates as
 * points in a three dimensional lattice.  The plot is animated so the user
 * can see the template bank from different perspectives.  See \figref{LALMathematicaHplot1}.
 *
 * \figure{LALMathematicaHplot1,eps,0.6,An example template bank produced by running InspiralSpinBankTest.c to generate roughly 5000 templates}
 *
 * Currently the plot doesn't show the contour of the templates; it renders them as
 * spheres.  In the case of metrics with disimilar scales along the
 * principle directions you will notice considerable space between points
 * accordingly.
 *
 * ### Notes ###
 *
 * <ul>
 * <li> The output of this function is &quot;Math3DNotebook.nb&quot; and will appear
 * in the directory of the program that called this function.</li>
 * <li> Exported MATHEMATICA graphics  will appear in your home directory
 * for unix users and in the Mathematica directory for Windows
 * users unless you have another path configured in your MATHEMATICA
 * installation. It is necessary to change the file name within the notebook
 * to avoid overwriting previous files.
 * </li>
 * </ul>
 *
 */
void
LALMath3DPlot ( LALStatus *stat,        /**< LALStatus structure pointer */
               Math3DPointList *first,  /**< Math3DPointList stucture pointer */
               INT4 *ntiles,            /**< pointer to the number of templates you \e plan to plot.
                                         * This may be called as NULL.  If it is called with a
                                         * value this function will check to see if the Math3DPointList has the
                                         * correct number of templates.  If it does not a warning will be printed. */
                REAL4 *pointSize        /**< \f$\epsilon[0,1]\f$ which specifies the relative
                                         * size of each point to the final display area.  (e.g. 1 would fill the
                                         * entire plot.)  This may be called as NULL and a calculated value will be
                                         * assigned.  (It's only a rough guess)
                                         */
                )

{
  FILE *nb; 				/* pointer to the notebook file */
  INT4 jflag = 0;			/* flag to justify the output data */
  Math3DPointList *list;		/* loop counter */
  REAL4 PtSize = 0.02;
  INT4 counter = 0;
  REAL4 xmax, ymax, zmax; /* maximum values plotted */
  INT2 xlog, ylog, zlog;  /* log10 of axis scaling factors */

  INITSTATUS(stat);

  if (!first) {
    ABORT(stat, LALMATHEMATICAH_ENULL, LALMATHEMATICAH_MSGENULL);
  }

  if ((nb = LALFopen("Math3DNotebook.nb", "w")) == NULL) {
    ABORT(stat, LALMATHEMATICAH_EFILE, LALMATHEMATICAH_MSGEFILE);
  }

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
        printf("\nWARNING!!! The value of argument ntiles (%i) != the Math3DPointList length (%i)\n",
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

  /* Find scaling factors for axes from maximum absolute values. */
  xmax = ymax = zmax = 0;
  for( list = first; list != NULL; list = list->next )
  {
    if( fabs( list->x ) > xmax )
      xmax = fabs( list->x );
    if( fabs( list->y ) > ymax )
      ymax = fabs( list->y );
    if( fabs( list->z ) > zmax )
      zmax = fabs( list->z );
  }
  xlog = (INT2)(log(xmax)/log(10));
  ylog = (INT2)(log(ymax)/log(10));
  zlog = (INT2)(log(zmax)/log(10));

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
        fprintf(nb, "AnimationPlot\t:=True");
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
        fprintf(nb, "XAxisLabel = \"Psi0 / 1e%d\"", xlog );
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "YAxisLabel = \"Psi3 / 1e%d\"", ylog );
      END_INPUTCELL;
      BEG_INPUTCELL;
        fprintf(nb, "ZAxisLabel = \"Beta / 1e%d\"", zlog );
      END_INPUTCELL;
      BEG_TEXTCELL;
        fprintf(nb, "AnimationPlot:\tFlag that determines whether to generate animations");
        fprintf(nb, "AnimationSize:\tThe size of the final animation in PIXELS x PIXELS\n");
        fprintf(nb, "StillSize:\t\tThe size of the final still image in PIXELS x PIXELS\n");
        fprintf(nb, "PtSize:\t\tThe relative size of the template points to the final display width.\n");
        fprintf(nb, "\t\t\tIt is given as a decimal part of one. (e.g. PtSize=0.02 is 1/20 of the display width)\n");
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
          fprintf(nb, "Graphics3D[{PointSize[PtSize]");
          list = first;
          while(list->next)
          {
            fprintf( nb, ",{GrayLevel[%f], Point[{%f,%f,%f}]}",
                     list->grayLevel, list->x/pow(10,xlog),
                     list->y/pow(10,ylog), list->z/pow(10,zlog) );
            if (jflag%2) fprintf(nb,"\n");
            ++jflag;
            list = list->next;
          }
          fprintf(nb, "}]");
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
          fprintf(nb, "If[AnimationPlot,{");
          fprintf(nb, "Do[tile[T]=Show[TILES, Background -> RGBColor[.93, .91, .89], ");
          fprintf(nb, "ViewPoint -> {1-(.99 T/frames)^2, T/(4 frames), 2 (T/frames)^2},ImageSize->AnimationSize], {T, 0, frames, 1}],\n");
          fprintf(nb, "Do[tile[frames+T]=Show[TILES, Background -> RGBColor[.93, .91, .89], ViewPoint -> {.005+(T/frames)^2, ");
          fprintf(nb, "0.25-T/(4 frames), 2-2 (.99 T/frames)^2},ImageSize->AnimationSize], {T, 0, frames, 1}]}];\n");
        END_INPUTCELL_;
      END_GROUPCELLC;
      BEG_GROUPCELL;
        BEG_SECTIONCELL;
          fprintf(nb, "Animation Generation");
        END_SECTIONCELL;
        BEG_INPUTCELL;
          fprintf(nb, "If[AnimationPlot,{images = Evaluate[Table[tile[j], {j, 0, 2 frames, 1}]]}];\n");
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "Export[StillName<>\".\"<>ToLowerCase[StillType], still, StillType, ImageSize->StillSize, ");
          fprintf(nb, "ConversionOptions->{\"ColorReductionDither\" -> False}]");
        END_INPUTCELL;
        BEG_INPUTCELL;
          fprintf(nb, "If[AnimationPlot,");
          fprintf(nb, "Export[AnimationName, images, \"GIF\", ImageSize -> AnimationSize, ");
          fprintf(nb, "ConversionOptions -> {\"Loop\" -> True,\"AnimationDisplayTime\" -> FrameTime, ");
          fprintf(nb, "\"ColorReductionDither\" -> False}]]");
        END_INPUTCELL_;
      END_GROUPCELLC_;
    END_GROUPCELLC_;
  END_NOTEBOOK;
  fclose(nb);
  RETURN(stat);
}
