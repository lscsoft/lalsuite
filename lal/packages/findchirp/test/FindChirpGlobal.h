#ifndef _FINDCHIRPGLOBAL_H
#define _FINDCHIRPGLOBAL_H

#ifndef _LALDATATYPES_H
#include <lal/LALDatatypes.h>
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

NRCSID (COMMTESTGLOBALH, "$Id$");

#ifdef FINDCHIRPGLOBAL_INIT
#define DECLARE(type,name,value) type name = (value)
#else
#define DECLARE(type,name,value) extern type name
#endif
DECLARE (INT4, numPoints, 1024);
DECLARE (INT4, numSegments, 2);
DECLARE (INT4, maxNumTemplates, 128);
#undef DECLARE

#endif
