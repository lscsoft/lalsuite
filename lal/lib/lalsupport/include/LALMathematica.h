/*
*  Copyright (C) 2007 Chad Hanna, Jolien Creighton, Benjamin Owen
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

#ifndef _LALMATHEMATICA_H
#define _LALMATHEMATICA_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * \addtogroup LALMathematica_h
 * \author Hanna, C. R.
 *
 * \brief Provides structures, functions and macro definitions for modules that
 * generate MATHEMATICA notebooks.
 *
 * Currently, the only modules using this
 * header file are <tt>LALMath3DPlot()</tt>, which generates 3D animated
 * plots of template banks having three parameters and
 * <tt>LALMathNDPlot()</tt> which plots the 3-dimensional projections of a
 * bank that is N-dimensional.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALMathematica.h>
 * \endcode
 *
 * This header file defines macros containing MATHEMATICA syntax that is
 * otherwise messy to implement into C source files.  Here is how to use
 * these macros to make your own program generate a MATHEMATICA notebook.
 *
 * <ol>
 * <li> Open a file with a pointer named &quot;nb&quot; and a file extenstion &quot;.nb&quot;.
 * </li><li> Use BEG_NOTEBOOK to start the notebook file.
 * </li><li> Use the appropriate BEG and END macros with fprint(nb, &quot;Your
 * Text&quot) in between to write your text to the cells of the notebook.  If
 * you are writing MATHEMATICA commands use the INPUT macros; for plain
 * text, use TEXT Macros.
 * </li><li> Denote titles and sections with the appropriate macros.
 * </li><li> Use END_NOTEBOOK to end the notebook and use fclose(nb) to close
 * the file &quot;nb&quot;.
 * </li></ol>
 *
 * The result is very readable/changeable source similar in style to most
 * markup languages. An example program might look like:
 * \code
 * FILE *nb;
 * nb = fopen("YourFileName.nb", "rw");
 * BEG_NOTEBOOK;
 * BEG_TITLECELL;
 * fprintf(nb, "Sample Program Title");
 * END_TITLECELL_;
 * BEG_SECTIONCELL;
 * fprintf(nb, "Sample Program Section Name");
 * END_SECTIONCELL;
 * END_NOTEBOOK;
 * fclose(nb);
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 * <li> Obviously the definitions and functions associated with this header
 * are NOT LAL compliant and thus do not belong in any lal routines except
 * test programs.</li>
 * <li> There are many more commands to manipulate the structure of
 * MATHEMATICA notebooks that are not included in this header.  The macros
 * are only what is necessary for a <em>bare minimum</em> interface.</li>
 * </ul>
 *
 */
/*@{*/

/**\name Error Codes */ /*@{*/
#define LALMATHEMATICAH_ENULL 1         /**< NULL pointer to a LALMathematica.h input structure */
#define LALMATHEMATICAH_EFILE 2         /**< Could not open file to write a Mathematica Notebook */
#define LALMATHEMATICAH_EVAL  3         /**< Invalid parameter value */
/*@}*/

/** \cond DONT_DOXYGEN */
#define LALMATHEMATICAH_MSGENULL        "NULL pointer to a LALMathematica.h input structure"
#define LALMATHEMATICAH_MSGEFILE        "Could not open file to write a Mathematica Notebook"
#define LALMATHEMATICAH_MSGEVAL         "Invalid parameter value"
/** \endcond */

/**
 * \name Macros
 * See the source file \ref LALMath3DPlot.c for an example of how to use
 * these macros to generate a MATHEMATICA notebook in your own program.
 *
 * <dl>
 * <dt>NOTEBOOK</dt><dd>Denotes the beginning and ending of the notebook file.  A
 * BEG_NOTEBOOK tag must start the file and an END_NOTEBOOK tag must end
 * it.</dd>
 * <dt>TITLE</dt><dd>Placing an fprint(nb, &quot;Your Title&quot;) between BEG and END
 * tags will place a <em>title font</em> cell in the notebook.</dd>
 * <dt>GROUP</dt><dd>Cells placed in between these tags will be grouped together</dd>
 * <dt>SECTION</dt><dd>Same as title except the text printed will be in
 * <em>section font</em>.  Subsequent input and text cells following the
 * END_SECTIONCELL tag will be grouped with that section until a new
 * BEG_SECTIONCELL tag is encountered.</dd>
 * <dt>INPUT</dt><dd>provides cells to input MATHEMATICA  commands.</dd>
 * <dt>TEXT</dt><dd>provides cells to input plain text.</dd>
 * </dl>
 *
 * Notice that the file pointer must be named &quot;nb&quot; in order to use the
 * macros defined in this header.  When grouping several cell objects
 * together the last object in the list should have an underscored END tag
 * instead of an END tag without an underscore.  Although the notebook will
 * compile (usually) if you use the tags without an ending underscore, the
 * dangling comma is taken as a null member of the list of grouped cells.
 * Therefore, when you view the notebook in MATHEMATICA  you may see the
 * word &quot;NULL&quot; printed on a line.  That is an indication that you should
 * use the underscore version of the tag which preceeded the &quot;NULL&quot;
 * statement.
 */
/*@{*/
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
/*@}*/


/**
 * This type is used by \ref LALMath3DPlot.c as an input structure to plot 3-dimensional template banks.
 * It is a linked list with
 * parameters for each coordinate x,y,z and a next pointer.  It also has a
 * parameter called grayLevel which must be \f$\epsilon [0,1]\f$.  It specifies
 * the shading of the point in the final plot with 0 representing black and
 * 1 representing white.  By creatively assigning its value the grayscale
 * shade of the points may convey additional information.
 */
typedef struct tagMath3DPointList{
  struct tagMath3DPointList *next;
  REAL4 x;
  REAL4 y;
  REAL4 z;
  REAL4 grayLevel;
  }Math3DPointList;

/**
 * This type is similar to Math3DPointList except the coordinates are
 * stored as data in the REAL4Vector \c coordinates.
 */
typedef struct tagMathNDPointList{
  struct tagMathNDPointList *next;
  REAL4Vector *coordinates;
  INT4 dimension;
  REAL4 grayLevel;
  } MathNDPointList;

/*@}*/

void
LALMath3DPlot( LALStatus *status,
               Math3DPointList *first,
               INT4 *ntiles,
               REAL4 *pointSize);

void
LALMathNDPlot( LALStatus *status,
               MathNDPointList *first,
               INT4 *ntiles,
               REAL4 *pointSize );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALMATHEMATICA_H */
