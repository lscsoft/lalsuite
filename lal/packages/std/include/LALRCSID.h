/*----------------------------------------------------------------------- 
 * 
 * File Name: LALRCSID.h
 * 
 * Author: Unknown?  Provided by Stuart Anderson
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * LALRCSID.h
 * 
 * SYNOPSIS 
 * #include "LALRCSID.h" 
 *
 * DESCRIPTION
 * Macros for defining RCS Id variables.  In each C or H file use:
 *
 *   NRCSID (FILENAMEC, "$Id$");
 *
 * or
 *
 *   NRCSID (FILENAMEH, "$Id$");
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALRCSID_H
#define _LALRCSID_H

#if !defined(lint)
#  ifndef __GNUC__
#    define RCSID(id)       static volatile const char *rcsid = (id)
#    define NRCSID(name,id) static volatile const char *name  = (id)
#  else
#    define RCSID(id) \
       static volatile const char __attribute__ ((unused)) *rcsid = (id)
#    define NRCSID(name,id) \
       static volatile const char __attribute__ ((unused)) *name  = (id)
#  endif /* !__GNUC__ */
#else
#  define RCSID(id)	   typedef void useless
#  define NRCSID(name,id)  typedef void useless
#endif /* !lint */

NRCSID (LALRCSIDH, "$Id$");

#endif
