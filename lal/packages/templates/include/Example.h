/*
*  Copyright (C) 2007 Jolien Creighton
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

/**** <lalVerbatim file="ExampleHV">
 * Author: Al A. Lal
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{Example.h}}
 *
 * %% One sentence briefly defining scope of the header
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Example.h>
 * \end{verbatim}
 *
 * %% Generic documentation on the header; this is the main place to
 * %% document any stuff not specific to the module
 *
 * \subsection*{Error conditions}
 * \input{ExampleHE}
 *
 * \subsection*{Structures}
 *
 * %% Document here any structures defined in the header.
 * %% Also include any of them in the index; e.g.:
 * %% \idx[Type]{ExampleOutput}
 * %% \idx[Type]{ExampleInput}
 * %% \idx[Type]{ExampleParams}
 *
 * \vfill{\footnotesize\input{ExampleHV}}
 * \newpage\input{ExampleC}
 * \newpage\input{ExampleTestC}
 *
 **** </lalLaTeX> */

#ifndef _EXAMPLE_H
#define _EXAMPLE_H

/** INCLUDE LAL HEADERS NEEDED FOR HEADER (NOT MODULE) **/
#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#pragma } /** to match the previous brace **/
#endif

/** DEFINE RCS ID STRING **/
NRCSID( EXAMPLEH, "$Id$" );

/** DEFINE ERROR CODES AND MESSAGES **/
/**** <lalErrTable file="ExampleHE"> */
#define EXAMPLEH_ENULL 1
#define EXAMPLEH_EALOC 2
#define EXAMPLEH_EIEIO 4
#define EXAMPLEH_EMOTE 8
#define EXAMPLEH_MSGENULL "Null pointer"
#define EXAMPLEH_MSGEALOC "Memory allocation error"
#define EXAMPLEH_MSGEIEIO "Have bought the barn"
#define EXAMPLEH_MSGEMOTE "I can't believe this is happening!"
/**** </lalErrTable> */

/** DEFINE NEW STRUCTURES AND TYPES **/
typedef struct tagExampleOutput { /** contents **/ } ExampleOutput;
typedef struct tagExampleInput  { /** contents **/ } ExampleInput;
typedef struct tagExampleParams { /** contents **/ } ExampleParams;

/** FUNCTION PROTOTYPES **/
void
LALExample(
    LALStatus     *status,
    ExampleOutput *output,
    ExampleInput  *input,
    ExampleParams *params
    );

#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _EXAMPLE_H */
