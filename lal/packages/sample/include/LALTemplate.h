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
