/*
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
#if 0
#include <lal/LIGOMetadataTables.h>
#endif

/**
 * \addtogroup UserInput_h
 * \author Reinhard Prix
 * \brief Module for simple unified handling of user-input from config-file and/or command-line.
 *
 * \heading{Description}
 *
 * This module provides simple function and macros to 'register' a set of C-variables as 'User Variables',
 * which can then be read in from the commandline and/or an input config file, as parsed by \ref ConfigFile_h.
 *
 * The module handles generating and outputting a help-string on the available inputs when requested, and
 * can deal with required inputs and providing defaults.
 *
 * \heading{Usage}
 *
 * The general approach consists of these steps
 * <ol>
 * <li> set default-value for optional user-variables</li>
 * <li> \c register all user-variables using calls to \c XLALRegister<TYPE>UserVar(), or more conveniently, using the shortcut-macros
 * XLALreg<TYPE>UserStruct() that assume a struct-pointer named 'uvar' containing all user-variables as 'uvar->UserVariable'.</li>
 * <li> parse all user-input using XLALUserVarReadAllInput()</li>
 * <li> At the end, free user-input structure</li>
 * </ol>
 *
 * One can use XLALUserVarWasSet() to determine wheter the user specified input for a given (optional) variable, or if it still has just its default value.
 *
 * The function XLALUserVarGetLog() can be used to obtain a log-string containing the full user-input, either in \c commandline- or \c ConfigFile format.
 *
 * Here is a worked simple example of its recommended use:
 * \code
 * #include <stdio.h>
 * #include <lal/UserInput.h>
 *
 * // these are the C-variables we want to read in from user-input
 * typedef struct {
 *   BOOLEAN help;                // did user request help-output?
 *   INT4 anInteger;
 *   REAL8 aDoubleVar;
 *   CHAR *andAString;
 *   REAL8 specialGeekSwitch;
 * } UserInput_t;
 *
 * UserInput_t empty_UserInput;    // this is zero-intialized!
 *
 * int main(int argc,char *argv[])
 * {
 *   UserInput_t UserVariables = empty_UserInput; // initializes this struct to {0}
 *   UserInput_t *uvar = &UserVariables;          // struct-pointer allows us to use the XLALreg<TYPE>UserStruct() macros...
 *
 *   // 1. step: set default-values for optional user-input variables
 *   uvar->anInteger = 0;
 *   uvar->andAString = NULL;     // Note: need to assign allocated strings here as default!!
 *
 *   // 2. step: Register all user-variables using the shortcut macros:
 *   XLALregBOOLUserStruct  ( help,               'h',  UVAR_HELP,     "Output this help-message");
 *   XLALregINTUserStruct   ( anInteger,          'i',  UVAR_OPTIONAL, "An example user-variable of an optional integer");
 *   XLALregREALUserStruct  ( aDoubleVar,         'r',  UVAR_REQUIRED, "This REAL8 user-variable is required");
 *   XLALregSTRINGUserStruct( andAString,           0,  UVAR_OPTIONAL, "Optional string-input, has no short-option");
 *   XLALregREALUserStruct  ( specialGeekSwitch,   'g',  UVAR_DEVELOPER, "This REAL8 user-variable is required");
 *
 *   // 3. step: parse all user-input, from either config-file if given, or commandline (overloads config-file values)
 *   if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS )
 *     XLAL_ERROR ( XLAL_EFUNC );
 *
 *   if (uvar->help)      // if user had requested help, then we're already done here
 *     return 0;
 *
 *   printf ("User-input was: anInteger = %d, aDoubleVar = %f, andAString = %s\n", uvar->anInteger, uvar->aDoubleVar, uvar->andAString );
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
 * \code
 * $ ./testUserInput -v1 --help
 *
 * Usage: testUserInput [@ConfigFile] [options], where options are:
 *
 * -v                        INT      set lalDebugLevel [0]
 * -h, --help                BOOL     Output this help-message []
 * -i, --anInteger           INT      An example user-variable of an optional integer [0]
 * -r, --aDoubleVar          REAL     This REAL8 user-variable is required [REQUIRED]
 * --andAString          STRING   Optional string-input, has no short-option [NULL]
 *
 * ---------- The following are 'Developer'-options not useful for most users:----------
 *
 * -g, --specialGeekSwitch   REAL     This REAL8 user-variable is required [0.0]
 * \endcode
 *
 * And if called correctly:
 * \code
 * $ ./testUserInput -r 3.1415 --andAString="stupid example"
 * User-input was: anInteger = 0, aDoubleVar = 3.141500, andAString = stupid example
 * \endcode
 *
 * \note For a real-world example of usage, see various codes in lalapps/src/pulsar, notably synthesizeLVStats.c
 *
 */
/*@{*/

/**
 * \name Shortcut Macros
 * With this family of short-cut macros one can conveniently register User-variables
 * that are accessible via a \e struct-pointer &quot;uvar->&quot;
 */
/*@{*/
#define XLALregREALUserStruct(name,option,flag,help) \
  XLALRegisterREALUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregINTUserStruct(name,option,flag,help) \
  XLALRegisterINTUserVar(#name, option,flag, help, &(uvar-> name))

#define XLALregBOOLUserStruct(name,option,flag,help) \
  XLALRegisterBOOLUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregSTRINGUserStruct(name,option,flag,help) \
  XLALRegisterSTRINGUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregLISTUserStruct(name,option,flag,help)                    \
  XLALRegisterLISTUserVar(#name, option, flag, help, &(uvar-> name))
/*@}*/

/** State-flags: variable is optional, required, help, developer or was_set */
typedef enum {
  UVAR_OPTIONAL		= 0,	/**< not required, and hasn't been set */
  UVAR_REQUIRED 	= 1<<0,	/**< we require the user to set this variable */
  UVAR_HELP		= 1<<1,	/**< special variable: trigger output of help-string */
  UVAR_DEVELOPER	= 1<<2,	/**< OPTIONAL and hidden in help-output at lalDebugLevel==0 */
  UVAR_SPECIAL		= 1<<3,	/**< OPTIONAL and *turns off* checking of required variables (LALUserVarCheckRequired) */
  UVAR_WAS_SET 		= 1<<7	/**< flag that this user-var has been set by user */
} UserVarState;

/**
 * Format for logging User-input: configFile- or cmdLine-style.
 * This determines the format of the string returned from LALLogUserInput().
 */
typedef enum {
  UVAR_LOGFMT_CFGFILE,	/**< return UserVars as a config-file */
  UVAR_LOGFMT_CMDLINE,	/**< return UserVars as a command-line */
  UVAR_LOGFMT_PROCPARAMS, /**< return UserVars suitable for filling in process-params struct */
  UVAR_LOGFMT_LAST
} UserVarLogFormat;

/* Function prototypes */
void XLALDestroyUserVars( void );
int XLALUserVarReadCmdline (int argc, char *argv[]);
int XLALUserVarReadCfgfile ( const CHAR *cfgfile );
CHAR *XLALUserVarHelpString ( const CHAR *progname );
int XLALUserVarReadAllInput ( int argc, char *argv[] );
int XLALUserVarCheckRequired( void );
int XLALUserVarWasSet (const void *cvar);
CHAR * XLALUserVarGetLog ( UserVarLogFormat format );

/* type-specific wrappers to XLALRegisterUserVar() to allow type-checking! */
int XLALRegisterREALUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, REAL8 *cvar );
int XLALRegisterINTUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, INT4 *cvar );
int XLALRegisterBOOLUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, BOOLEAN *cvar );
int XLALRegisterSTRINGUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, CHAR **cvar );
int XLALRegisterLISTUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, LALStringVector **cvar);


/* ========== Deprecated LAL interface wrappers ========== */

/** \name Error codes */
/*@{*/
#define USERINPUTH_ENULL        1       /**< Arguments contained an unexpected null pointer. */
#define USERINPUTH_ENONULL      2       /**< Output pointer is not NULL */
#define USERINPUTH_EMEM         3       /**< Out of memory */
#define USERINPUTH_EOPT         4       /**< Unknown command-line option encountered */
#define USERINPUTH_ENOUVARS     5       /**< No user-variables have been registered! */
#define USERINPUTH_ECMDLARG     6       /**< Illegal command-line argument */
#define USERINPUTH_EUNKNOWN     7       /**< Unknown user-variable */
#define USERINPUTH_ENOTSET      8       /**< Required user-variable was not set */
#define USERINPUTH_EDEBUG       9       /**< lalDebugLevel can only be read before ANY mallocs(), even hidden */
#define USERINPUTH_EONECONFIG   10      /**< Currently one ONE config-file can be specified using '\@' */
#define USERINPUTH_ERECFORMAT   11      /**< Unknown format for recording user-input */
#define USERINPUTH_EXLAL        12      /**< Failure in XLAL function */
#define USERINPUTH_ENAMECOLL    13      /**< Commandline option assigned more than once */
/*@}*/

/** \cond DONT_DOXYGEN */
#define USERINPUTH_MSGENULL 	"Arguments contained an unexpected null pointer."
#define USERINPUTH_MSGENONULL	"Output pointer is not NULL"
#define USERINPUTH_MSGEMEM	"Out of memory"
#define USERINPUTH_MSGEOPT	"Unknown command-line option encountered"
#define USERINPUTH_MSGENOUVARS	"No user-variables have been registered!"
#define USERINPUTH_MSGECMDLARG   "Illegal command-line argument"
#define USERINPUTH_MSGEUNKNOWN	"Unknown user-variable"
#define USERINPUTH_MSGENOTSET	"Required user-variable was not set"
#define USERINPUTH_MSGEDEBUG	"lalDebugLevel can only be read before ANY mallocs(), even hidden.."
#define USERINPUTH_MSGEONECONFIG "Currently one ONE config-file can be specified using '@'"
#define USERINPUTH_MSGERECFORMAT   "Unknown format for recording user-input"
#define USERINPUTH_MSGEXLAL	"Failure in XLAL function"
#define USERINPUTH_MSGENAMECOLL "Commandline option assigned more than once"
/** \endcond */

/**
 * \name Deprecated LAL-interface
 * These functions and macros are deprecated, and you should user their XLAL-equivalents instead.
 */
/*@{*/
/*----- short-cut macros to register global "uvar_" User-Variables ----- */
#define LALregREALUserVar(status,name,option,flag,help) \
TRY(LALRegisterREALUserVar((status)->statusPtr, #name, option, flag, help,&(uvar_ ## name)), status)

#define LALregINTUserVar(status,name,option,flag,help) \
TRY(LALRegisterINTUserVar((status)->statusPtr, #name, option,flag, help,&(uvar_ ## name)), status)

#define LALregBOOLUserVar(status,name,option,flag,help) \
TRY(LALRegisterBOOLUserVar((status)->statusPtr, #name, option, flag, help, &(uvar_ ## name)),status)

#define LALregSTRINGUserVar(status,name,option,flag,help) \
TRY(LALRegisterSTRINGUserVar((status)->statusPtr, #name, option, flag, help, &(uvar_ ## name)),status)

#define LALregLISTUserVar(status,name,option,flag,help) \
TRY(LALRegisterLISTUserVar((status)->statusPtr, #name, option, flag, help, &(uvar_ ## name)),status)

/*----- another family of short-cut macros: register _struct-pointer_ "uvar->" User-Variables ----- */
#define LALregREALUserStruct(status,name,option,flag,help) \
TRY(LALRegisterREALUserVar((status)->statusPtr, #name, option, flag, help, &(uvar-> name)), status)

#define LALregINTUserStruct(status,name,option,flag,help) \
TRY(LALRegisterINTUserVar((status)->statusPtr, #name, option,flag, help, &(uvar-> name)), status)

#define LALregBOOLUserStruct(status,name,option,flag,help) \
TRY(LALRegisterBOOLUserVar((status)->statusPtr, #name, option, flag, help, &(uvar-> name)),status)

#define LALregSTRINGUserStruct(status,name,option,flag,help) \
TRY(LALRegisterSTRINGUserVar((status)->statusPtr, #name, option, flag, help, &(uvar-> name)),status)

#define LALregLISTUserStruct(status,name,option,flag,help) \
TRY(LALRegisterLISTUserVar((status)->statusPtr, #name, option, flag, help, &(uvar-> name)),status)


void LALRegisterREALUserVar(LALStatus *, const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, REAL8 *cvar);
void LALRegisterINTUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, INT4 *cvar);
void LALRegisterBOOLUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, BOOLEAN *cvar);
void LALRegisterSTRINGUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, CHAR **cvar);
void LALRegisterLISTUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, LALStringVector **cvar);

void LALDestroyUserVars (LALStatus *);

void LALUserVarReadAllInput(LALStatus *, int argc, char *argv[]);
void LALUserVarReadCmdline (LALStatus *, int argc, char *argv[]);
void LALUserVarReadCfgfile (LALStatus *, const CHAR *cfgfile);

void LALUserVarHelpString (LALStatus *, CHAR **helpstring, const CHAR *progname);
void LALUserVarCheckRequired (LALStatus *);
INT4 LALUserVarWasSet (const void *cvar);
void LALUserVarGetLog (LALStatus *, CHAR **logstr,  UserVarLogFormat format);
#if 0
void LALUserVarGetProcParamsTable (LALStatus *status, ProcessParamsTable **out, CHAR *progname);
#endif

/*@}*/

/*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
