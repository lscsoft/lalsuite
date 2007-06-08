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

/********************************** <lalVerbatim file="LALTemplateHV">
Author: Creighton, T. D.
$Id$
*********************************** </lalVerbatim> */

/* Description goes here */

#ifndef _IIRFILTER_H    /* Protect against double-inclusion */
#define _IIRFILTER_H

#include <lal/LALStdlib.h>  /* Include any other headers */

#ifdef __cplusplus
extern "C" {            /* Protect against C++ name mangling */
#endif

/* Define the RCS ID string */
NRCSID(LALTEMPLATEH,"$Id$");

/* Define error codes and messages.  These must be auto-extracted for
   inclusion into the documentation */
/***************************************************** <lalErrTable> */
#define LALTEMPLATEH_EONE 1
#define LALTEMPLATEH_ETWO 2

#define LALTEMPLATEH_MSGEONE "An error condition"
#define LALTEMPLATEH_MSGETWO "Another error condition"
/*************************************************** </lalErrTable> */

/* Define other global constants or macros */

/* Define new structures and types */

/* Include external global variables */

/* Declare global function prototypes */
void
LALTemplate( LALStatus *stat );

#ifdef __cplusplus
}                   /* Close C++ protection */
#endif

#endif              /* Close double-include protection */
