/*----------------------------------------------------------------------- 
 * 
 * File Name: Example.c 
 * 
 * Author: J. Random Hacker
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * Example
 * 
 * SYNOPSIS 
 * (void) Example()
 * 
 * DESCRIPTION 
 * Example source file prolog. 
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 * (list of LLAL, LDAS, other non-system functions/procedures called. 
 * 
 * NOTES
 * (Other notes)
 * 
 *-----------------------------------------------------------------------
 */

/* 1. Prolog: an extended comment field containing summary information
 *    about the source code module (see above);
 * 2. Source file version string (from CVS). Note the string name.
 */

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
#endif

NRCSID (EXAMPLEC, "$Id$");

/*
 * 3. Include directives. These should be guarded and appear in the
 *    following order:
 *    a. Standard library includes;
 *    b. LDAS includes;
 *    c. LAL includes.
 *    Each source file must include at least its own header. Includes should 
 *    be guarded, as in 
 */

#ifndef _STDLIB_H
#include <stdlib.h>
#ifndef _STDLIB_H
#define _STDLIB_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LDASDCAPI_H
#include "LDASDCAPI.h"
#ifndef _LDASDCAPI_H
#define _LDASDCAPI_H
#endif
#endif

#ifndef _EXAMPLE_H
#include "Example.h"
#ifndef _EXAMPLE_H
#define _EXAMPLE_H
#endif
#endif

/* 
 * 4. Constants, enumerated types, structures, etc., used only internally;
 */

/*
 * 5. Type declarations ({\em i.e.,} {\tt typedefs\/}) used only
 *    internally; 
 */

/*
 * 6. Function macros for which a waiver has been granted;
 */

/*
 * 7. Extern global variable declarations for which a waiver has been
 *    granted;
 */

/*
 * 8. Static global variables for which a waiver has been granted; 
 */

/*
 * 9. Static function declarations for which a waiver has been granted; 
 */

/*
 * 10. Function prototypes
 */
