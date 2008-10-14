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

/********************************** <lalVerbatim file="StringTokenCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StringToken.c}}
\label{ss:StringToken.c}

Converts a string into a series of tokens, for use by other routines.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StringTokenCP}
\idx{LALCreateTokenList()}
\idx{LALDestroyTokenList()}

\subsubsection*{Description}

The routine \verb@LALCreateTokenList()@ parses \verb@*string@ as a
sequence of tokens (substrings of non-null characters that do not
appear in \verb@delimiters@), separated by delimiters (substrings
consisting only of characters that appear in \verb@delimiters@), and
terminated by the null character \verb@'\0'@.  The structure
\verb@**list@ is created, storing the sequence of tokens as a list
null-terminated character strings.

The output \verb@list@ should be a non-\verb@NULL@ handle that points
to the value \verb@NULL@ (i.e.\ \verb@list@$\neq$\verb@NULL@ but
\verb@*list@=\verb@NULL@).  Even if no tokens were found, \verb@*list@
will be created, but will have \verb@(*list)->nTokens@=0,
\verb@(*list)->tokens[0]@=\verb@NULL@, and
\verb@(*list)->list@=\verb@NULL@.  Note that this is \emph{not} an
error, so the calling routine need not guarantee in advance that
\verb@string@ contain any non-delimiter characters.

The routine \verb@LALDestroyTokenList()@ destroys a list of tokens as
created by \verb@LALCreateTokenList()@, setting \verb@*list@ to
\verb@NULL@.

\subsubsection*{Algorithm}

The \verb@LALCreateTokenList()@ function is not particularly
memory-efficient, requiring internal storage up to twice the length of
\verb@*string@.  It first creates a working copy of
\verb@string->data@, and replaces all occurences of characters
appearing in \verb@*delimiters@ with \verb@'\0'@, while at the same
time keeping track of the number and total length of all tokens.  It
then allocates a contiguous block of memory to store all the tokens
(separated by and terminated with single \verb@'\0'@ characters), and
a set of \verb@CHAR *@ pointers to point to the individual tokens in
this block.  Then the routine proceeds through the working copy one
last time, copying tokens into the token list and setting the token
pointers accordingly, before destroying the working copy.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                     LALCHARCreateVector()
LALFree()                       LALCHARDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StringTokenCV}}

******************************************************* </lalLaTeX> */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>

NRCSID( STRINGTOKENC, "$Id$" );

/* <lalVerbatim file="StringTokenCP"> */
void
LALCreateTokenList( LALStatus  *stat,
		    TokenList  **list,
		    const CHAR *string,
		    const CHAR *delimiters )
{ /* </lalVerbatim> */
  BOOLEAN delimiter = 1; /* whether current character is a delimiter */
  UINT4 i = 0, j = 0;    /* indecies */
  UINT4 nTokens = 0;     /* number of tokens */
  UINT4 sLength;         /* length of string */
  UINT4 tLength = 0;     /* length of token list */
  CHAR *copy;            /* working copy of token list */

  INITSTATUS( stat, "LALCreateTokenList", STRINGTOKENC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( list, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( delimiters, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( !*list, stat, STRINGINPUTH_EOUT, STRINGINPUTH_MSGEOUT );

  /* Create working copy of token list. */
  sLength = strlen( string ) + 1;
  if ( !( copy = (CHAR *)LALMalloc( sLength*sizeof(CHAR) ) ) ) {
    ABORT( stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  }
  for ( i = 0; i < sLength; i++ ) {
    CHAR c = string[i];
    if ( strchr( delimiters, c ) ) {
      copy[i] = '\0';
      delimiter = 1;
    } else {
      copy[i] = c;
      tLength++;
      if ( delimiter ) {
	delimiter = 0;
	nTokens++;
      }
    }
  }

  /* Create the token list. */
  if ( !( *list = (TokenList *)LALMalloc( sizeof(TokenList) ) ) ) {
    LALFree( copy );
    ABORT( stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  }
  if ( !( (*list)->tokens =
	  (CHAR **)LALMalloc( ( nTokens + 1 )*sizeof(CHAR *) ) ) ) {
    LALFree( *list );
    *list = NULL;
    LALFree( copy );
    ABORT( stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  }
  (*list)->nTokens = nTokens;
  (*list)->list = NULL;


  /* If tokens were found, copy them over and set up pointers. */
  if ( nTokens ) {
    CHAR *listData; /* pointer to token list data */
    LALCHARCreateVector( stat->statusPtr, &((*list)->list),
			 nTokens + tLength );
    BEGINFAIL( stat ) {
      LALFree( (*list)->tokens );
      LALFree( *list );
      *list = NULL;
      LALFree( copy );
    } ENDFAIL( stat );
    listData = (*list)->list->data;
    i = 0;
    while ( i < sLength ) {
      if ( copy[i] ) {
	tLength = strlen( copy + i ) + 1;
	memcpy( listData, copy + i, tLength*sizeof(CHAR) );
	(*list)->tokens[j++] = listData;
	i += tLength;
	listData += tLength;
      } else
	i++;
    }
  }
  (*list)->tokens[j] = NULL;

  /* Clean up and exit. */
  LALFree( copy );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}



/* <lalVerbatim file="StringTokenCP"> */
int
XLALCreateTokenList(TokenList  **list,
		    const CHAR *string,
		    const CHAR *delimiters )
{ /* </lalVerbatim> */
  BOOLEAN delimiter = 1; /* whether current character is a delimiter */
  UINT4 i = 0, j = 0;    /* indecies */
  UINT4 nTokens = 0;     /* number of tokens */
  UINT4 sLength;         /* length of string */
  UINT4 tLength = 0;     /* length of token list */
  CHAR *copy;            /* working copy of token list */

  /* Check for valid input arguments. */
  if (!list || !string || !delimiters || *list != NULL) {
      fprintf(stderr, STRINGINPUTH_MSGENUL);
      return STRINGINPUTH_ENUL;
  }

  /* Create working copy of token list. */
  sLength = strlen( string ) + 1;
  if ( !( copy = (CHAR *)LALMalloc( sLength*sizeof(CHAR) ) ) ) {
      fprintf(stderr, STRINGINPUTH_MSGENUL);
      return STRINGINPUTH_ENUL;
  }
  for ( i = 0; i < sLength; i++ ) {
    CHAR c = string[i];
    if ( strchr( delimiters, c ) ) {
      copy[i] = '\0';
      delimiter = 1;
    } else {
      copy[i] = c;
      tLength++;
      if ( delimiter ) {
	delimiter = 0;
	nTokens++;
      }
    }
  }

  /* Create the token list. */
  if ( !( *list = (TokenList *)LALMalloc( sizeof(TokenList) ) ) ) {
    LALFree( copy );
    fprintf(stderr, STRINGINPUTH_MSGENUL);
    return STRINGINPUTH_ENUL;
  }
  if ( !( (*list)->tokens =
	  (CHAR **)LALMalloc( ( nTokens + 1 )*sizeof(CHAR *) ) ) ) {
    LALFree( *list );
    *list = NULL;
    LALFree( copy );
    fprintf(stderr, STRINGINPUTH_MSGENUL);
    return STRINGINPUTH_ENUL;
  }
  (*list)->nTokens = nTokens;
  (*list)->list = NULL;


  /* If tokens were found, copy them over and set up pointers. */
  if ( nTokens ) {
    CHAR *listData; /* pointer to token list data */
    (*list)->list = XLALCreateCHARVector(nTokens + tLength);

    if ( ! (*list)->list ) {
      LALFree( (*list)->tokens );
      LALFree( *list );
      *list = NULL;
      LALFree( copy );
    }
    listData = (*list)->list->data;
    i = 0;
    while ( i < sLength ) {
      if ( copy[i] ) {
	tLength = strlen( copy + i ) + 1;
	memcpy( listData, copy + i, tLength*sizeof(CHAR) );
	(*list)->tokens[j++] = listData;
	i += tLength;
	listData += tLength;
      } else
	i++;
    }
  }
  (*list)->tokens[j] = NULL;

  /* Clean up and exit. */
  LALFree( copy );

  return 0;
}



/* <lalVerbatim file="StringTokenCP"> */
void
LALDestroyTokenList( LALStatus *stat,
		     TokenList **list )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "LALCreateTokenList", STRINGTOKENC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( list, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( *list, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Free everything and exit. */
  if ( (*list)->list ) {
    TRY( LALCHARDestroyVector( stat->statusPtr, &((*list)->list) ),
	 stat );
  }
  LALFree( (*list)->tokens );
  LALFree( *list );
  *list = NULL;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StringTokenCP"> */
int
XLALDestroyTokenList( TokenList **list )
{ /* </lalVerbatim> */
  /* Check for valid input arguments. */
  if ( !list || !*list ) {
    fprintf(stderr, STRINGINPUTH_MSGENUL);
    return STRINGINPUTH_ENUL;
  }
  /* Free everything and exit. */
  if ( (*list)->list )
    XLALDestroyCHARVector( (*list)->list );

  LALFree( (*list)->tokens );
  LALFree( *list );
  *list = NULL;

  return 0;
}
