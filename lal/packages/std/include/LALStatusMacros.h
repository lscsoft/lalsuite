/*----------------------------------------------------------------------- 
 * 
 * File Name: LALStatusMacros.h
 * 
 * Author: Creighton, J. D. E. and Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * LALStatusMacros.h
 * 
 * SYNOPSIS 
 * #include "LALStatusMacros.h"
 *
 * DESCRIPTION
 *
 * This header file defines programming macros for handling LAL status
 * pointers.  The intent is simultaneously to standardize the error
 * reporting, and to make the reporting as transparent as possible to
 * people coding individual routines.
 *
 * The following summarized everything the common programmer needs to
 * know in order to follow LAL standard error reporting.
 *
 * 0. The Status structure
 *
 *    All error reporting is handled by a structure type named Status.
 *    This structure has the following fields:
 *
 *    INT4 statusCode: A code indicating the exit status of a
 *      function.  0 represents a normal exit.  Negative values are
 *      reserved for certain standard error types.  The authors of
 *      individual functions should assign positive values to the
 *      various ways in which their code can fail.
 *
 *    const CHAR *statusDescription: An explanatory string
 *      corresponding to the numerical status code.
 *
 *    const CHAR *Id: A character string identifying the source file
 *      and version number of the function being reported on.
 *
 *    const CHAR *file: The file name of the actual .c file containing
 *      the function code.
 *
 *    INT4 line: The line number in the .c file of the instruction
 *      where any error was reported.
 *
 *    Status *statusPtr: A recursive pointer to another status
 *      pointer.  This structure is used to report an error in a
 *      subroutine of the current function.  Thus if an error occurs
 *      in a deeply-nested routine, the status structure returned to
 *      the main program will be the head of a linked list of status
 *      structures, one for each nested level, with the tail structure
 *      reporting the actual error that caused the overlying routines
 *      to fail.
 *
 *    INT4 level: The nested-function level where any error was reported.
 *
 *    The standard status codes are as follows:
 *
 *     0: Nominal execution: the function returned successfully.
 *
 *    -1: Recursive error: the function aborted due to failure of a
 *        subroutine.
 *
 *    -2: The status structure passed to the function had a non-null
 *        statusPtr field, which blocks the function from calling
 *        subroutines (it is symptomatic of something screwy going on
 *        in the calling routine).
 *
 *    -4: Memory allocation error: the function was unable to allocate
 *        the statusPtr field to pass down to a subroutine.
 *
 *    -8: The statusPtr could not be deallocated at the end of all
 *        subroutine calls; one of the subroutines must have lost it
 *        or set it to null.
 *
 *
 * 1. Every source file should have a unique character string
 *    identifying that version of that file.  The standard convention,
 *    for a file Myfunk.c, is to declare a string at the top of the
 *    module using the macro NRCSID (defined in the include file LALRCSID.h):
 *
 *      NRCSID (MYFUNKC, "\044Id\044");
 *
 *    where \044Id\044 is expanded by RCS to give the full name and
 *    version number of the source file.
 *
 * 2. All functions should have return type void.  The first argument
 *    of any function should be a pointer to a structure of type
 *    Status.  Thus:
 *
 *      void Myfunk(Status *stat, ... )
 *
 *    Since the function has no return code, it must report all errors
 *    or failure through the status structure.  A function that is
 *    passed a NULL pointer in place of the status pointer should
 *    abort the program, as this is its only way to report the error.
 *    However, this is the only circumstance under which a function
 *    sould normally abort; other errors should be trapped, reported
 *    in the status structure, and control returned to the calling
 *    routine.
 *
 * 3. The first instruction in any function, after variable
 *    declarations, should be the macro INITSTATUS(), which takes two
 *    arguments: the function's status pointer, and the module's RCS
 *    Id string.
 *
 *        INITSTATUS(stat, MYFUNKC);
 *
 *    This macro checks that a valid status pointer has been passed to
 *    the function, and if so, initializes the other fields to
 *    indicate (by default) nominal execution.  If stat is null, the
 *    macro causes the program to abort, this being its only way to
 *    report that something went wrong.
 *
 * 4. Upon completion, the function should issue the macro RETURN(),
 *    which takes one argument: the function's status pointer.
 *
 *        RETURN(stat);
 *
 *    This takes the place of any return statements, and may also log
 *    status reports to some suitable log file (often stderr),
 *    depending on implementation and the value of a global debuglevel
 *    parameter.  Typically RETURN() is used only for successful
 *    completion, with other macros ABORT(), ASSERT(), and TRY() being
 *    used to report failure.  However, it is possible for the
 *    programmer to assign the fields of stat by hand, and then issue
 *    RETURN(stat).
 *
 * 5. The standard method to terminate a function unsuccessfully is
 *    with the ABORT() macro, which takes three arguments: the status
 *    pointer, the status code, and the status description string.
 *    Normally the various error codes and descriptions will be
 *    constants defined in the function's header file:
 *
 *        ABORT(stat, MYFUNK_EMYERR, MYFUNK_MSGEMYERR);
 *
 *    where the error code MYFUNK_EMYERR and the error message
 *    MYFUNK_MSGEMYERR are defined in the function's header file.
 *    Like RETURN(), ABORT() correctly handles any status logging
 *    required by the implementation and the debuglevel.  Notably,
 *    ABORT() does not raise a SIGABRT flag, but instead returns
 *    control to the calling routine.
 *
 * 6. Another way to indicate an unsuccessful termination is with the
 *    macro ASSERT(), which takes as arguments a test statement, a
 *    status pointer, a status code, and a status description.  The
 *    statement ASSERT(assertion,...); is in all ways equivalent to
 *    the statement if(!assertion) ABORT(...);.  For instance, in the
 *    above example, one might have:
 *
 *        ASSERT(assertion, stat, MYFUNK_EMYERR, MYFUNK_MSGEMYERR)
 *
 *    Coding purists may argue that ASSERT() should be used only to
 *    trap coding errors rather than runtime errors, which would be
 *    trapped using ABORT().  In other words, the assertion should
 *    always test true in the final debugged program.  At present,
 *    however, this coding practice is not enforced by the LAL
 *    standard.
 *
 * 7. If the function is to call other LAL functions as subroutines,
 *    four more macros must be used to report possible errors arising
 *    in these routines.  The macros are ATTATCHSTATUSPTR(),
 *    DETATCHSTATUSPTR(), CHECKSTATUSPTR(), and TRY().  The usage of
 *    these macros is as follows.
 *
 *    a. First, before any subroutines are called, the function must
 *       call the macro ATTATCHSTATUSPTR() which takes as its argument
 *       the status pointer of the current function:
 *
 *           ATTATCHSTATUSPTR(stat);
 *
 *       This allocates stat->statusPtr, which is the status pointer
 *       that will be handed down into any and all subroutines.  If
 *       the pointer has already been allocated, ATTATCHSTATUSPTR()
 *       will abort, as this is symptomatic of a coding error.  Note
 *       that ATTATCHSTATUSPTR() need only be called once in a given
 *       function, no matter how many subroutine calls that function
 *       makes.  Normally it should be called immediately after
 *       INITSTATUS().
 *
 *    b. When a subroutine is called, it should be handed the
 *       statusPtr field of the calling functions status structure, to
 *       report its own errors.  The calling function should test the
 *       returned status code, and either attempt to deal
 *       with any abnormal returns, or abort with status code -1.  The
 *       macro CHECKSTATUSPTR() helps simplify this procedure.  It takes one
 *       arguments: the status pointer of the current function (not
 *       the subroutine).
 *
 *           mysubroutine(stat->statusPtr,...);
 *           CHECKSTATUSPTR(stat);
 *
 *       The TRY() macro is a somewhat more streamlined approach but
 *       with equivalent results.  It takes two arguments.  The first
 *       is the subroutine call, and the second is the status pointer.
 *       Thus:
 *
 *           TRY(mysubroutine(stat->statusPtr,...), stat);
 *
 *       The only practical difference between these two approaches is
 *       that the TRY() macro takes up only one line, and can also
 *       report the name of the failed subroutine call when logging
 *       errors.
 *
 *    c. After all subroutines have been called, but before any
 *       RETURN() statement, the function must call the
 *       DETATCHSTATUSPTR() macro, with the status pointer of the
 *       current function (not the subroutines) as its argument:
 *
 *           DETATCHSTATUSPTR(stat);
 *
 *       This simply deallocates stat->statusPtr and sets it to NULL.
 *       It is an error to exit the function with non-NULL statusPtr,
 *       unless the exit was due to a subroutine failure.  ABORT() and
 *       ASSERT() check for this automatically; the only place you
 *       need to call DETATCHSTATUSPTR() is immediately before
 *       RETURN().
 *
 * 8. The REPORTSTATUS() macro is used to issue a current status
 *    report from the current function; the report is printed to
 *    stderr or some other implementation-specific stream used for
 *    error reporting.  REPORTSTATUS() takes the current status
 *    pointer as its argument, and iteratively reports down any chain
 *    of non-NULL statusPtr structures passed back by subroutines.
 *
 * 9. The top-level main function should start with an empty (all fields set
 *    to zero) status structure, as in the following example:
 *
 *      int
 *      main ()
 *      {
 *        static Status status;
 *        MYFUNK (&status);
 *        REPORTSTATUS (&status);
 *        return 0;
 *      }
 *
 *
 * Non-Conformant Functions:
 *
 * These standards apply only to functions that will be publicly
 * available in the LAL libraries.  Within a module, a programmer may
 * define and use subroutines that do not conform to the LAL function
 * standards, provided these routines are only visible within that
 * module.  Such functions should be declared as static to ensure
 * this.  A publicly-visible non-conformant function requires special
 * dispensation.
 *
 *
 * Example:
 *
 * As an example to illustrate these standards, here is a pair of
 * (rather silly) routines to perform inversion and division of real
 * numbers.  In the file division.h one might have the following error
 * codes and error messages defined:
 *
 *   #define DIVISION_ENULL 1
 *   #define DIVISION_EDIV0 2
 *   #define DIVISION_MSGEDIV0 "Null pointer"
 *   #define DIVISION_MSGEDIV0 "Division by zero"
 *
 * followed by function prototypes.  The file division.c might contain
 * the following code:
 *
 *   static const char *DIVISIONC="\044Id\044";
 *
 *   ReturnCode InvertREAL4(Status *stat,
 *                          REAL4  *output,
 *                          REAL4  input)
 *   {
 *     INITSTATUS(stat,DIVISIONC);
 *     ASSERT(output!=NULL,stat,DIVISION_ENULL,DIVISION_MSGENULL);
 *     if(input==0.0)
 *       ABORT(stat,DIVISION_EDIV0,DIVISION_MSGEDIV0);
 *     *output = 1.0/input;
 *     RETURN(stat);
 *   }
 *
 *   ReturnCode DivideREAL4(Status *stat,
 *                          REAL4  *output,
 *                          REAL4  numerator,
 *                          REAL4  denominator)
 *   {
 *     REAL4 invDenom;
 *
 *     INITSTATUS(stat,DIVISIONC);
 *     ATTATCHSTATUSPTR(stat);
 *     TRY(InvertREAL4(stat->statusPtr,&invDenom,denominator),stat);
 *     *output = numerator*invDenom;
 *     DETATCHSTATUSPTR(stat);
 *     RETURN(stat);
 *   }
 *
 *
 * NOTES
 *
 * Why are the status handling routines written as macros rather than
 * functions?  There are three good reasons.
 *
 * First, many of the handling routines must be able to force an exit
 * from the function calling them.  This cannot be done if the routine
 * is in its own function, except by raising signal flags (which is a
 * Bad Thing according to LAL standards).
 *
 * Second, it is useful for these routines to assign a status
 * structure's file and line fields using the __FILE__ and __LINE__
 * macros.  If the routine is its own function, then these will just
 * give the file and line number where the error handling routine is
 * defined.  If the routine is a macro, then these will give the file
 * and line number where the macro was called, which is much more
 * interesting.
 *
 * Third, by expanding macros at compile time, the runtime performance
 * of the resulting code is marginally better.  Most of these macros
 * will, under nominal conditions, reduce to a single conditional test
 * of an integer value, with no additional overhead from function
 * calling and parameter passing.  Thus programmers can be encouraged
 * to include extensive error trapping in all their routines, without
 * having to worry about compromising performance.
 *
 *
 *----------------------------------------------------------------------- 
 */


#ifndef _LALSTATUSMACROS_H
#define _LALSTATUSMACROS_H

#ifndef _LALMALLOC_H
#include "LALMalloc.h"
#ifndef _LALMALLOC_H
#define _LALMALLOC_H
#endif
#endif

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

#ifndef _LALERROR_H
#include "LALError.h"
#ifndef _LALERROR_H
#define _LALERROR_H
#endif
#endif

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
#endif

NRCSID (LALSTATUSMACROSH, "$Id$");

extern int debuglevel;

#define INITSTATUS(statusptr, id)                                     \
do                                                                    \
{                                                                     \
  INT4 level;                                                         \
  if(!(statusptr))                                                    \
    {                                                                 \
      CHAR msg[1024];                                                 \
      sprintf(msg,"Abort: line %d, file %s, %s\n"                     \
	      "       Null status pointer passed to function\n",      \
	      __LINE__,__FILE__,(id));                                \
      LALAbort(msg);                                                  \
    }                                                                 \
  level = (statusptr)->level;                                         \
  memset((statusptr),0,sizeof(Status));                               \
  (statusptr)->level = level > 0 ? level : 1 ;                        \
  (statusptr)->Id    = (id);                                          \
} while (0)

#define RETURN(statusptr)                                             \
do                                                                    \
{                                                                     \
  (statusptr)->file=__FILE__;                                         \
  (statusptr)->line=__LINE__;                                         \
  if(debuglevel==0 ||((debuglevel==1)&&((statusptr)->statusCode==0))) \
    {                                                                 \
      return;                                                         \
    }                                                                 \
  else if((statusptr)->statusCode==0)                                 \
    {                                                                 \
      LALPrintError("Nominal[%d]: line %d, file %s, %s\n",            \
	      (statusptr)->level,(statusptr)->line,(statusptr)->file, \
              (statusptr)->Id);                                       \
      return;                                                         \
    }                                                                 \
  else                                                                \
    {                                                                 \
      LALPrintError("Error[%d] %d: line %d, file %s, %s\n"            \
              "          %s\n",(statusptr)->level,                    \
	      (statusptr)->statusCode,(statusptr)->line,              \
              (statusptr)->file,(statusptr)->Id,                      \
              (statusptr)->statusDescription);                        \
      return;                                                         \
    }                                                                 \
} while (0)

#define ABORT(statusptr,code,mesg)                                    \
do                                                                    \
{                                                                     \
  (statusptr)->file=__FILE__;                                         \
  (statusptr)->line=__LINE__;                                         \
  (statusptr)->statusCode=(code);                                     \
  (statusptr)->statusDescription=(mesg);                              \
  if((statusptr)->statusPtr)                                          \
    {                                                                 \
      LALFree((statusptr)->statusPtr);                                \
      (statusptr)->statusPtr=NULL;                                    \
    }                                                                 \
  if(debuglevel==0 || ((debuglevel==1)&&((code)==0)))                 \
    {                                                                 \
      return;                                                         \
    }                                                                 \
  else if((code)==0)                                                  \
    {                                                                 \
      LALPrintError("Nominal[%d]: line %d, file %s, %s\n",            \
	      (statusptr)->level,(statusptr)->line,(statusptr)->file, \
              (statusptr)->Id);                                       \
      return;                                                         \
    }                                                                 \
  else                                                                \
    {                                                                 \
      LALPrintError("Error[%d] %d: line %d, file %s, %s\n"            \
              "         %s\n",(statusptr)->level,(code),              \
              (statusptr)->line,(statusptr)->file,(statusptr)->Id,    \
              (mesg));                                                \
      return;                                                         \
    }                                                                 \
} while (0)

#define ASSERT(assertion,statusptr,code,mesg)                         \
do                                                                    \
{                                                                     \
  if(!(assertion))                                                    \
    {                                                                 \
      (statusptr)->file=__FILE__;                                     \
      (statusptr)->line=__LINE__;                                     \
      (statusptr)->statusCode=(code);                                 \
      (statusptr)->statusDescription=(mesg);                          \
      if((statusptr)->statusPtr)                                      \
	{                                                             \
	  LALFree((statusptr)->statusPtr);                            \
	  (statusptr)->statusPtr=NULL;                                \
	}                                                             \
      if(debuglevel==0 || ((debuglevel==1)&&((code)==0)))             \
	{                                                             \
	  return;                                                     \
	}                                                             \
      else if((code)==0)                                              \
	{                                                             \
	  LALPrintError("Nominal[%d]: line %d, file %s, %s\n",        \
		  (statusptr)->level,(statusptr)->line,               \
                  (statusptr)->file,(statusptr)->Id);                 \
	  return;                                                     \
	}                                                             \
      else                                                            \
	{                                                             \
	  LALPrintError("Error[%d] %d: line %d, file %s, %s\n"        \
		  "         Assertion %s failed: %s\n",               \
                  (statusptr)->level,(code),(statusptr)->line,        \
                  (statusptr)->file,(statusptr)->Id,#assertion,       \
                  (mesg));                                            \
	  return;                                                     \
	}                                                             \
    }                                                                 \
} while (0)

#define ATTATCHSTATUSPTR(statusptr)                                   \
do                                                                    \
{                                                                     \
  ASSERT(!(statusptr)->statusPtr,statusptr,-2,                        \
	 "ATTATCHSTATUSPTR: non-null status pointer");                \
  (statusptr)->statusPtr=(Status *)LALCalloc(1,sizeof(Status));       \
  ASSERT((statusptr)->statusPtr,statusptr,-4,                         \
	 "ATTATCHSTATUSPTR: memory allocation error");                \
  (statusptr)->statusPtr->level=(statusptr)->level + 1;               \
} while (0)

#define DETATCHSTATUSPTR(statusptr)                                   \
do                                                                    \
{                                                                     \
  ASSERT((statusptr)->statusPtr,statusptr,-8,                         \
	 "DETATCHSTATUSPTR: null status pointer");                    \
  LALFree((statusptr)->statusPtr);                                    \
  (statusptr)->statusPtr=NULL;                                        \
} while (0)

#define TRY(function,statusptr)                                       \
do                                                                    \
{                                                                     \
  (function);                                                         \
  if((statusptr)->statusPtr->statusCode)                              \
    {                                                                 \
      (statusptr)->file=__FILE__;                                     \
      (statusptr)->line=__LINE__;                                     \
      (statusptr)->statusCode=-1;                                     \
      (statusptr)->statusDescription="Recursive error";               \
      if(debuglevel>0)                                                \
	{                                                             \
	  LALPrintError("Error[%d] %d: line %d, file %s, %s\n"        \
		  "          Function call %s failed\n",              \
                  (statusptr)->level,-1,(statusptr)->line,            \
                  (statusptr)->file,(statusptr)->Id,#function);       \
	}                                                             \
      return;                                                         \
    }                                                                 \
} while (0)

#define CHECKSTATUSPTR(statusptr)                                     \
do                                                                    \
{                                                                     \
  if((statusptr)->statusPtr->statusCode)                              \
    {                                                                 \
      (statusptr)->file=__FILE__;                                     \
      (statusptr)->line=__LINE__;                                     \
      (statusptr)->statusCode=-1;                                     \
      (statusptr)->statusDescription="Recursive error";               \
      if(debuglevel>0)                                                \
	{                                                             \
	  LALPrintError("Error[%d] %d: line %d, file %s, %s\n"        \
		  "          Function call failed\n",                 \
                  (statusptr)->level,-1,(statusptr)->line,            \
                  (statusptr)->file,(statusptr)->Id);                 \
	}                                                             \
      return;                                                         \
    }                                                                 \
} while (0)

#define REPORTSTATUS(statusptr)                                       \
do                                                                    \
{                                                                     \
  Status *ptr;                                                        \
  for(ptr=(statusptr);ptr;ptr=(ptr->statusPtr))                       \
    {                                                                 \
      LALPrintError("\nLevel %i: %s\n",ptr->level,ptr->Id);           \
      if (ptr->statusCode)                                            \
      {                                                               \
        LALPrintError("\tStatus code %i: %s\n",ptr->statusCode,       \
	              ptr->statusDescription);                        \
      }                                                               \
      else                                                            \
      {                                                               \
        LALPrintError("\tStatus code 0: Nominal\n");                  \
      }                                                               \
      LALPrintError("\tfile %s, line %i\n",ptr->file,ptr->line);      \
    }                                                                 \
} while (0)

#endif
