/*
 * Copyright (C) 2016 Karl Wette
 * Copyright (C) 2015 Reinhard Prix
 * Copyright (C) 2010 Reinhard Prix (xlalified)
 * Copyright (C) 2004, 2005 Reinhard Prix
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

#ifndef _USERINPUT_H  /* Double-include protection. */
#define _USERINPUT_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

#include <lal/ConfigFile.h>
#include <lal/UserInputParse.h>

/**
 * \defgroup UserInput_h Header UserInput.h
 * \ingroup lal_support
 * \author Reinhard Prix
 * \brief Module for simple unified handling of user-input from config-file and/or command-line.
 *
 *
 * ### Description ###
 *
 * This module provides functions and macros to 'register' a set off C-variables as 'User Variables',
 * which can then be read in from the commandline and/or an input config file, as parsed by \ref ConfigFile_h.
 *
 * The module also handles generating and outputting a help-string on the available input options when requested, and
 * can deal with enforcing input of required options and using default values.
 *
 * ### Usage ###
 *
 * The general approach consists of the following steps:
 * <ol>
 * <li> set default-values for optional user-variables as appropriate</li>
 * <li> \c register all user-variables using calls to \c XLALRegister<TYPE>UserVar(), or more conveniently, using the shortcut-macros
 * XLALreg<TYPE>UserStruct() that assume a pointer named 'uvar' to a struct containing all user-variables in the form 'uvar->UserVariable'.</li>
 * <li> parse user-input using XLALUserVarReadAllInput()</li>
 * </ol>
 *
 * One can use XLALUserVarWasSet() to determine whether a given user-input option has been set by the user.
 *
 * The function XLALUserVarGetLog() can be used to obtain a log-string containing the full user-input, either in \c commandline- or \c ConfigFile format.
 *
 * Here is a worked simple example of its recommended use:
 * \code
 * #include <stdio.h>
 * #include <lal/XLALError.h>
 * #include <lal/LALDatatypes.h>
 *
 * #include <lal/UserInput.h>
 *
 * // these are the C-variables we want to read in from user-input
 * typedef struct {
 *   BOOLEAN help;                // did user request help-output?
 *   INT4 anInteger;
 *   REAL8 aDoubleVar;
 *   CHAR *andAString;
 *   REAL8 specialGeekSwitch;
 *   LIGOTimeGPS someEpoch;
 *   REAL8 RA;
 *   REAL8 DEC;
 * } UserInput_t;
 *
 *
 * int main(int argc,char *argv[])
 * {
 *   UserInput_t XLAL_INIT_DECL(UserVariables); // initializes this struct to {0}
 *   UserInput_t *uvar = &UserVariables;          // struct-pointer allows us to use the XLALreg<TYPE>UserStruct() macros...
 *
 *   // 1. step: set default-values for optional user-input variables
 *   uvar->anInteger = 0;
 *   uvar->andAString = NULL;     // Note: need to assign allocated strings here as default!!
 *
 *   // 2. step: Register all user-variables using the shortcut macros:
 *   XLALRegisterUvarMember( help,               BOOLEAN, 'h',  HELP,      "Output this help-message");
 *   XLALRegisterUvarMember( anInteger,          INT4,    'i',  OPTIONAL,  "An example user-variable of an optional integer");
 *   XLALRegisterUvarMember( aDoubleVar,         REAL8,   'r',  REQUIRED,  "This REAL8 user-variable is required");
 *   XLALRegisterUvarMember( andAString,         STRING,   0,   OPTIONAL,  "Optional string-input, has no short-option");
 *   XLALRegisterUvarMember( someEpoch,          EPOCH,    0,   OPTIONAL,  "Reference epoch (format 'xx.yy[GPS|MJD]')");
 *   XLALRegisterUvarMember( RAJ,                RAJ,      0,   OPTIONAL,  "Sky location: equatorial right ascension in [0,2pi] (in radians or hours:minutes:seconds)");
 *   XLALRegisterUvarMember( DEC,                DECJ,     0,   OPTIONAL,  "Sky location: equatorial declination [-pi/2,pi/2] (in radians or degrees:minutes:seconds)");
 *   XLALRegisterUvarMember( specialGeekSwitch,  REAL8,    'g', DEVELOPER, "This REAL8 user-variable may not be relevant for standard usage");
 *
 *   // 3. step: parse all user-input, from either config-file if given, or commandline (overloads config-file values)
 *   XLAL_CHECK ( XLALUserVarReadAllInput ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC);
 *
 *   if (uvar->help){      // if user had requested help, then we're already done here
 *     return 0;
 *   }
 *
 *   printf ("User-input was: anInteger = %d, aDoubleVar = %f, andAString = %s\n", uvar->anInteger, uvar->aDoubleVar, uvar->andAString );
 *   printf ("someEpoch = {%d s, %d ns}, RA = %f rad, DEC = %f rad\n", uvar->someEpoch.gpsSeconds, uvar->someEpoch.gpsNanoSeconds, uvar->RA, uvar->DEC );
 *
 *   // 4. step: free user-input module memory
 *   XLALDestroyUserVars();
 *
 *   LALCheckMemoryLeaks();
 *   return 0;
 * } // main()
 * \endcode
 *
 * \note This code can be compiled <b>as is</b> within lalapps, and yields
 *
 * \verbatim
$ LAL_DEBUG_LEVEL=1 ./tmp --help

Usage: ./tmp [@ConfigFile] [options], where options are:

-h, --help                BOOLEAN    Output this help-message []
-i, --anInteger           INT4       An example user-variable of an optional integer [0]
-r, --aDoubleVar          REAL8      This REAL8 user-variable is required [REQUIRED]
    --andAString          STRING     Optional string-input, has no short-option [NULL]
    --someEpoch           EPOCH      Reference epoch (format 'xx.yy[GPS|MJD]') [0.000000000GPS]
    --RA                  RAJ        Sky location: right ascension in [0,2pi] (in radians or hours:minutes:seconds) [0.0]
    --DEC                 DECJ       Sky location: declination [-pi/2,pi/2] (in radians or degrees:minutes:seconds) [0.0]

---------- The following are 'Developer'-options not useful for most users:----------

-g, --specialGeekSwitch   REAL8      This REAL8 user-variable may not be relevant for standard usage [0.0]

\endverbatim
 *
 * And if called
 * \verbatim
$ ./tmp -r 3.1415 --andAString="stupid example" --someEpoch=5555MJD --RA=10:25:10.123 --DEC=-30:0:0
User-input was: anInteger = 0, aDoubleVar = 3.141500, andAString = stupid example
someEpoch = {2147483596 s, 816000000 ns}, RA = 2.727813 rad, DEC = -0.523599 rad
\endverbatim
 *
 * \note For a real-world example of usage, see various codes under lalapps/src/pulsar/{Injections,Fstatistic}
 *
 */
/*@{*/

/**
 * Shortcut Macro for registering new user variables, which are assumed
 * to be accessible via the \e struct-pointer '*uvar'
 */
#define XLALRegisterUvarMember(name,type,option,category,...)             \
  XLALRegister ##type## UserVar( &(uvar-> name), #name, option, UVAR_CATEGORY_ ## category, __VA_ARGS__)

/// (mutually exclusive) UserVariable categories: optional, required, help, developer, ...
typedef enum {
  UVAR_CATEGORY_START = 0,	///< internal start marker for range checking

  UVAR_CATEGORY_OPTIONAL,	///< optional
  UVAR_CATEGORY_REQUIRED,	///< required
  UVAR_CATEGORY_DEVELOPER,	///< optional and hidden in help-output until lalDebugLevel>=warning
  UVAR_CATEGORY_DEPRECATED,	///< optional and hidden until lalDebugLevel>=info; still supported but output warning if used
  UVAR_CATEGORY_DEFUNCT,	///< hidden completely from help output; not supported, will output error + help-string if used
  UVAR_CATEGORY_SPECIAL,	///< optional and *turns off* all checking of required variables, useful for output of version info
  UVAR_CATEGORY_NODEFAULT,	///< optional and supresses printing the default value in the help, where it doesn't make sense
  UVAR_CATEGORY_END		///< internal end marker for range checking
} UserVarCategory;

// ***** the following are provided only for backwards compatibility until full transition to new XLAL interface *****
#define UVAR_OPTIONAL  	UVAR_CATEGORY_OPTIONAL
#define UVAR_REQUIRED  	UVAR_CATEGORY_REQUIRED
#define UVAR_DEVELOPER	UVAR_CATEGORY_DEVELOPER
#define UVAR_SPECIAL	UVAR_CATEGORY_SPECIAL
// **********

/**
 * Format for logging User-input: configFile- or cmdLine-style.
 * This determines the format of the string returned from XLALUserVarGetLog().
 */
typedef enum {
  UVAR_LOGFMT_CFGFILE,		/**< return UserVars as a config-file */
  UVAR_LOGFMT_CMDLINE,		/**< return UserVars as a command-line */
  UVAR_LOGFMT_PROCPARAMS, 	/**< return UserVars suitable for filling in process-params struct */
  UVAR_LOGFMT_LAST
} UserVarLogFormat;

/* Global variables */
#ifndef SWIG /* exclude from SWIG interface */
extern const char *lalUserVarHelpBrief;
#endif /* SWIG */

/* Function prototypes */
void XLALDestroyUserVars( void );
int XLALUserVarReadCmdline( BOOLEAN *should_exit, int argc, char *argv[] );
int XLALUserVarReadCfgfile( BOOLEAN *should_exit, const CHAR *cfgfile );
int XLALUserVarReadAllInput( BOOLEAN *should_exit, int argc, char *argv[] );
int XLALUserVarWasSet( const void *cvar );
void XLALUserVarCheck( BOOLEAN *should_exit, const int assertion, const CHAR *fmt, ... ) _LAL_GCC_PRINTF_FORMAT_(3,4);
CHAR * XLALUserVarGetLog ( UserVarLogFormat format );

/**
 * \name Convenience macros for checking how many of a set of user input variables were set
 */
/*@{*/
#define UVAR_SET(n)                             (XLALUserVarWasSet(&(uvar-> n)) ? 1 : 0)
#define UVAR_SET2(n1,n2)                        (UVAR_SET(n1) + UVAR_SET(n2))
#define UVAR_SET3(n1,n2,n3)                     (UVAR_SET2(n1,n2) + UVAR_SET(n3))
#define UVAR_SET4(n1,n2,n3,n4)                  (UVAR_SET3(n1,n2,n3) + UVAR_SET(n4))
#define UVAR_SET5(n1,n2,n3,n4,n5)               (UVAR_SET4(n1,n2,n3,n4) + UVAR_SET(n5))
#define UVAR_SET6(n1,n2,n3,n4,n5,n6)            (UVAR_SET5(n1,n2,n3,n4,n5) + UVAR_SET(n6))
/*@}*/

/**
 * \name Convenience macros for checking whether all of a set of user input variables were set
 */
/*@{*/
#define UVAR_ALLSET2(n1,n2)                     (UVAR_SET2(n1,n2) == 2)
#define UVAR_ALLSET3(n1,n2,n3)                  (UVAR_SET3(n1,n2,n3) == 3)
#define UVAR_ALLSET4(n1,n2,n3,n4)               (UVAR_SET4(n1,n2,n3,n4) == 4)
#define UVAR_ALLSET5(n1,n2,n3,n4,n5)            (UVAR_SET5(n1,n2,n3,n4,n5) == 5)
#define UVAR_ALLSET6(n1,n2,n3,n4,n5,n6)         (UVAR_SET6(n1,n2,n3,n4,n5,n6) == 6)
/*@}*/

/**
 * \name Convenience macros for checking whether any of a set of user input variables were set
 */
/*@{*/
#define UVAR_ANYSET2(n1,n2)                     (UVAR_SET2(n1,n2) > 0)
#define UVAR_ANYSET3(n1,n2,n3)                  (UVAR_SET3(n1,n2,n3) > 0)
#define UVAR_ANYSET4(n1,n2,n3,n4)               (UVAR_SET4(n1,n2,n3,n4) > 0)
#define UVAR_ANYSET5(n1,n2,n3,n4,n5)            (UVAR_SET5(n1,n2,n3,n4,n5) > 0)
#define UVAR_ANYSET6(n1,n2,n3,n4,n5,n6)         (UVAR_SET6(n1,n2,n3,n4,n5,n6) > 0)
/*@}*/

/**
 * \name Convenience macros for printing user input variables in error messages
 */
/*@{*/
#define UVAR_FMT                                "`--%s'"
#define UVAR_STR(n)                             "`--"#n"'"
#define UVAR_STR2AND(n1,n2)                     "`--"#n1"' and `--"#n2"'"
#define UVAR_STR2OR(n1,n2)                      "`--"#n1"' or `--"#n2"'"
#define UVAR_STR3AND(n1,n2,n3)                  "`--"#n1"', `--"#n2"', and `--"#n3"'"
#define UVAR_STR3OR(n1,n2,n3)                   "`--"#n1"', `--"#n2"', or `--"#n3"'"
#define UVAR_STR4AND(n1,n2,n3,n4)               "`--"#n1"', `--"#n2"', `--"#n3"', and `--"#n4"'"
#define UVAR_STR4OR(n1,n2,n3,n4)                "`--"#n1"', `--"#n2"', `--"#n3"', or `--"#n4"'"
#define UVAR_STR5AND(n1,n2,n3,n4,n5)            "`--"#n1"', `--"#n2"', `--"#n3"', `--"#n4"', and `--"#n5"'"
#define UVAR_STR5OR(n1,n2,n3,n4,n5)             "`--"#n1"', `--"#n2"', `--"#n3"', `--"#n4"', or `--"#n5"'"
#define UVAR_STR6AND(n1,n2,n3,n4,n5,n6)         "`--"#n1"', `--"#n2"', `--"#n3"', `--"#n4"', `--"#n5"', and `--"#n6"'"
#define UVAR_STR6OR(n1,n2,n3,n4,n5,n6)          "`--"#n1"', `--"#n2"', `--"#n3"', `--"#n4"', `--"#n5"', or `--"#n6"'"
/*@}*/

// declare type-specific wrappers to XLALRegisterUserVar() to allow for strict C type-checking!
#define DECL_REGISTER_UVAR(UTYPE,CTYPE)                                 \
  int XLALRegister ##UTYPE## UserVar ( CTYPE *cvar, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *fmt, ... ) _LAL_GCC_PRINTF_FORMAT_(5,6)

// ------ declare registration functions
DECL_REGISTER_UVAR(REAL8,REAL8);
DECL_REGISTER_UVAR(INT4,INT4);
DECL_REGISTER_UVAR(INT8,INT8);
DECL_REGISTER_UVAR(BOOLEAN,BOOLEAN);
DECL_REGISTER_UVAR(EPOCH,LIGOTimeGPS);
DECL_REGISTER_UVAR(RAJ,REAL8);
DECL_REGISTER_UVAR(DECJ,REAL8);
DECL_REGISTER_UVAR(STRING,CHAR*);

DECL_REGISTER_UVAR(REAL8Range,REAL8Range);
DECL_REGISTER_UVAR(EPOCHRange,LIGOTimeGPSRange);
DECL_REGISTER_UVAR(RAJRange,REAL8Range);
DECL_REGISTER_UVAR(DECJRange,REAL8Range);

DECL_REGISTER_UVAR(STRINGVector,LALStringVector*);
DECL_REGISTER_UVAR(REAL8Vector,REAL8Vector*);
DECL_REGISTER_UVAR(INT4Vector,INT4Vector*);

/* ========== Deprecated LAL interface wrappers ========== */

/**
 * \name Deprecated LAL-interface
 * These functions and macros are deprecated, and you should user their XLAL-equivalents instead.
 */
/*@{*/

void LALRegisterREALUserVar(LALStatus *, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *helpstr, REAL8 *cvar);
void LALRegisterINTUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *helpstr, INT4 *cvar);
void LALRegisterBOOLUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *helpstr, BOOLEAN *cvar);
void LALRegisterSTRINGUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *helpstr, CHAR **cvar);
void LALRegisterLISTUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *helpstr, LALStringVector **cvar);
void LALDestroyUserVars (LALStatus *);
void LALUserVarReadAllInput(LALStatus *, BOOLEAN *should_exit, int argc, char *argv[]);
void LALUserVarReadCmdline (LALStatus *, BOOLEAN *should_exit, int argc, char *argv[]);
void LALUserVarReadCfgfile (LALStatus *, BOOLEAN *should_exit, const CHAR *cfgfile);
INT4 LALUserVarWasSet (const void *cvar);
void LALUserVarGetLog (LALStatus *, CHAR **logstr,  UserVarLogFormat format);

/*@}*/

/*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
