#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALVersion.h>

int lalDebugLevel = 0;

int main( void )
{
  static LALStatus status;
  char msg[1024];
  int verbose = 1;

  if ( strcmp( LAL_VERSION, lalVersion ) ||
       strcmp( LAL_CONFIGURE_ARGS, lalConfigureArgs ) ||
       strcmp( LAL_CONFIGURE_DATE, lalConfigureDate ) )
  {
    fputs( "LAL Version Mismatch!\n\n", stderr );
    fputs( "Header Version ",           stderr );
    fputs( LAL_VERSION,                 stderr );
    fputs( "\nCompiled on ",            stderr );
    fputs( LAL_CONFIGURE_DATE,          stderr );
    fputs( "\nWith arguments ",         stderr );
    fputs( LAL_CONFIGURE_ARGS,          stderr );
    fputs( "\n\n",                      stderr );
    fputs( "Library Version ",          stderr );
    fputs( lalVersion,                  stderr );
    fputs( "\nCompiled on ",            stderr );
    fputs( lalConfigureDate,            stderr );
    fputs( "\nWith arguments ",         stderr );
    fputs( lalConfigureArgs,            stderr );
    fputs( "\n",                        stderr );
    return 1;
  }

  LALVersion( &status, msg, sizeof( msg ), verbose );

  if ( status.statusCode )
  {
    LALStatus *next = &status;
    do
    {
      fputs( next->statusDescription, stderr );
      fputs( "\n", stderr );
      next = next->statusPtr;
    }
    while ( next );
    return 2;
  }

  puts( msg );

  return 0;
}
