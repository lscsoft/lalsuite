/*----------------------------------------------------------------------- 
 * 
 * File Name: CommTestGlobal.h
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _COMMTESTGLOBAL_H
#define _COMMTESTGLOBAL_H

#ifndef _LALDATATYPES_H
#include <lal/LALDatatypes.h>
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

NRCSID (COMMTESTGLOBALH, "$Id$");

#ifdef COMMTESTGLOBAL_INIT
#define DECLARE(type,name,value) type name = (value)
#else
#define DECLARE(type,name,value) extern type name
#endif
DECLARE (INT4, numPoints, 16);
#undef DECLARE

#endif
