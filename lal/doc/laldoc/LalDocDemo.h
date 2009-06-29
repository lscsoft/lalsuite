/*
*  Copyright (C) 2007 Alan Wiseman
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


(1)  LaTeX:

/*
 * <lalLaTeX file="LalDocDemo_LaTeX_This_File">
 *
 * \documentclass[oneside]{book}
 * \setlength{\textheight}{9.0in}
 * \setlength{\textwidth}{6.0in}
 * \setlength{\topmargin}{-0.00in}
 * \setlength{\oddsidemargin}{-0.25in}
 * \sloppy
 *
 * \begin{document}
 *
 * \chapter{The \LaTeX environment}
 *
 * The \LaTeX source that made this file was extracted from the
 * LalDocDemo.h. It was surrounded by the keyword, extracted and stuffed
 * into LalDocDemo\_LaTeX\_This\_File.tex.
 *
 * When programming in c, many programmers like to put *'s down the
 * side of their comments. If you look at the input file
 * ({\tt LalDocDemo.h}), you will see that {\tt laldoc} is smart
 * enough to strip them off when in the \LaTeX enviroment.
 *
 * \chapter{The Verbatim environment}
 *
 * The  following is the result of verbatim extraction.  This was extracted
 * and stored in the default file LalDocDemoH.tex. The name having
 * been constructed from the input file name.  It was then stuffed in
 * with an \verb@\input{}@ command.
 *
 * \input{LalDocDemoH}
 *
 * \vspace*{1.5in}
 * The following was also extracted from the input file, but it was
 * stuffed into the file {\tt LalDocDemoHMyVerbatim.tex}
 *
 * \vspace*{0.5in}
 * \input{LalDocDemoHMyVerbatim}
 *
 * \chapter{A Good Error Table}
 * \input{LalDocDemoHErrCodeGood}
 *
 * \vspace*{0.5in}
 * The caption below the table is also supplied the {\tt laldoc}.
 *
 *
 * \chapter{A  Not So Good Error Table}
 * \input{LalDocDemoHErrCodeBad}
 *
 * \vspace*{0.5in}
 * The caption below the table was also  supplied the {\tt laldoc}.
 * If you look at the source code from which this table was generated,
 * you will notice that laldoc was smart enough to notice that the error
 * codes in this table do not follow the naming convention and flagged
 * the problem in the caption.  This problem was also flagged in
 * the error file Error.out.
 *
 * \end{document}
 * </lalLaTeX>
 *
*/

/* ---------------------------------------------------------------------- */

(2) Verbatim with default output file name:

<lalVerbatim>

This junk was extracted and stuffed into the default file:
{\tt LalDocDemoH.tex}.  Also note the margin par on the right. This
tells you exactly where this came from.  You can put all the garbage
you want in here.

</lalVerbatim>

/* ---------------------------------------------------------------------- */

(3) Verbatim with specified output file name:

<lalVerbatim file=LalDocDemoHMyVerbatim>
This junk was extracted and stuffed into the file LalDocDemoHMyVerbatim.tex.
Again, the margin par tells where it came from.

</lalVerbatim>

/* ---------------------------------------------------------------------- */

(4) An Error Table that follows the convention:

/* <lalErrTable file="LalDocDemoHErrCodeGood"> */

#define LALDOCDEMOH_EANERROR 1
#define LALDOCDEMOH_MSGEADIFFERENTERROR "This is a different error message"
/* Comments in the middle of the Error Table enviroment will be ignored. */
#define LALDOCDEMOH_EALONGERROR  3

#define LALDOCDEMOH_MSGEANERROR "This is an error message"
/* putting the message before the code is ok. */
#define LALDOCDEMOH_EADIFFERENTERROR  2
#define LALDOCDEMOH_MSGEALONGERROR "This is very very long error message that
goes on for several lines.  There is a length limit, but the limit is
set by {\tt MAXSTR}. Currently it is 1024 characters."

/* </lalErrTable> */

/* ---------------------------------------------------------------------- */

(4) An Error Table that follows the convention:

/* <lalErrTable file="LalDocDemoHErrCodeBad"> */

#define INPUTH_EANERROR 1
/* Comments in the middle will be ignored. */
#define INPUTH_EADIFFERENTERROR  2
#define INPUTH_EALONGERROR  3

#define BADNAME_MSGEANERROR "This is an error message"
#define BADNAME_MSGEADIFFERENTERROR "This is a different error message"
#define BADNAME_MSGEALONGERROR "This is very very long error message that
goes on for several lines.  There is a length limit, but the limit is
set by {\tt MAXSTR}. Currently it is 1024 characters."


/* </lalErrTable> */

