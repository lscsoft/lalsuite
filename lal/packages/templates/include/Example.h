/*-----------------------------------------------------------------------
 *
 * File Name: Example.h
 *
 * Author: Hacker, J. R.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Example.h
 *
 * SYNOPSIS
 * #include "Example.h"
 *
 * DESCRIPTION
 * Example header file prolog.
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

/* 
 * Header contents go here, in order specified:
 * 
 * 1. Prolog (Comment field with file name, author, revision, etc., as 
 *    specified above)
 * 2. include-loop protection (see below). Note the naming convention!
 */

#ifndef _EXAMPLE_H
#define _EXAMPLE_H

/*
 * 3. Includes. This header may include others; if so, they go immediately 
 *    after include-loop protection. Includes should appear in the following 
 *    order: 
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 *    Includes should be double-guarded!
 */

#ifndef _STDLIB_H
#include <stdlib.h>
#ifndef _STDLIB_H
#define _STDLIB_H
#endif
#endif

#ifndef _LDAS_CONSTANTS_H
#include "LDAS_CONSTANTS.h"
#ifndef _LDAS_CONSTANTS_H
#define _LDAS_CONSTANTS_H
#endif
#endif

/*
 * 4. Header file version string (from CVS; see below). Note the string name. 
 */

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
#endif

NRCSID (EXAMPLEH, "$Id$");

/*
 * 5. Macros. But, note that macros are deprecated. 
 */

/* 
 * 6. Extern Constant Declarations. These should not be present unless a
 *    specific waiver has been granted. 
 */

/* 
 * 7. Extern Global Variables. These should also not be present unless a 
 *    specific waiver has been granted. 
 */ 

/*
 * 8. Structure, enum, union, etc., typdefs.
 */

/*
 * 9. Functions Declarations (i.e., prototypes).
 */

#endif
