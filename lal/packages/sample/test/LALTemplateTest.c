/****************************** <lalVerbatim file="LALTemplateTestCV">
Author: Creighton, T. D.
$Id$
******************************* </lalVerbatim> */

/* Usage and description */

#include <lal/LALTemplate.h>  /* Include any required headers */

/* Define RCS ID string */
NRCSID(LALTEMPLATETESTC,"$Id$");

/* Define local constants and macros */

/* Declare and set global lalDebugLevel */
int lalDebugLevel = 0;

/* Declare local (static) functions (definitions can go here or at the
   end of the file) */

/* Define main function */
int main( int argc, char **argv )
{
  /* Variable declarations */
  static LALStatus stat;

  /* Code */
  LALTemplate( &stat );
  return 0;
}
