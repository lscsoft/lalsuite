/*----------------------------------------------------------------------- 
 * 
 * File Name: LALHello.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALHELLO_H
#define _LALHELLO_H

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

NRCSID (LALHELLOH, "$Id$");

#define LALHELLO_EOPEN  1
#define LALHELLO_EWRITE 2
#define LALHELLO_ECLOSE 4
#define LALHELLO_EFLUSH 8

#define LALHELLO_MSGEOPEN  "Could not open file"
#define LALHELLO_MSGEWRITE "Error in writing to file"
#define LALHELLO_MSGECLOSE "Error in closing file"
#define LALHELLO_MSGEFLUSH "Error in flushing stdout"

void
LALHello (Status *status, const CHAR *fileName);

#endif
