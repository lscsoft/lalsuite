/******************************** <lalVerbatim file="DetectorSiteTestCV">
Author: Whelan, J. T.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{DetectorSiteTest.c}}
\label{tools:ss:DetectorSiteTest.c}

Tests the detector response and site parameter structures and the
routine to create one from the other.

\subsubsection*{Usage}
\begin{verbatim}

\subsubsection*{Usage}
\begin{verbatim}
DetectorSiteTest [options]
Options:
  -h         print help
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\end{verbatim}

\subsubsection*{Description}

Right now the test routine does very little.  It contains a static
function \texttt{PrintLALDetector()} which will print the fields of a
\texttt{LALDetector} to standard output in the same format that would
be used for a C initialization.  This function is not currently
called.  It also contains a static function \texttt{CheckDetector()}
which extracts the \texttt{LALFrDetector} and type from a
\texttt{LALDetector}, changes the name of the \texttt{LALFrDetector}
(in case it's one of the predefined constant detectors), constructs a
new \texttt{LALDetector} and compares the values of the fields of the
old and new structures.  The program currently performs this check
for the two LIGO sites.


\subsubsection*{Exit codes}
\input{DetectorSiteTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALCreateDetector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DetectorSiteTestCV}}
******************************************************* </lalLaTeX> */

#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

int lalDebugLevel = LALMSGLVL1;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

NRCSID( DETECTORSITETESTC, "$Id$" );

#define DETECTORSITETESTC_LOCTOL 3e-2
#define DETECTORSITETESTC_RESTOL 1e-06

/***************************** <lalErrTable file="DetectorSiteTestCE"> */
#define DETECTORSITETESTC_ENOM 0
#define DETECTORSITETESTC_ECHK 1
#define DETECTORSITETESTC_EFLS 2
#define DETECTORSITETESTC_MSGENOM "Nominal exit"
#define DETECTORSITETESTC_MSGECHK "Error checking failed to catch bad data"
#define DETECTORSITETESTC_MSGEFLS "Incorrect answer for valid data"
/***************************** </lalErrTable> */

static void
PrintLALDetector(LALDetector *detector)
{
  printf( "  = { { %.15e, %.15e, %.15e },\n", detector->location[0],
          detector->location[1], detector->location[2] );
  printf( "      { { %.15g, %.15g, %.15g },\n", detector->response[0][0],
          detector->response[0][1], detector->response[0][2]);
  printf( "        { %.15g, %.15g, %.15g },\n", detector->response[1][0],
          detector->response[1][1], detector->response[1][2]);
  printf( "        { %.15g, %.15g, %.15g } },\n", detector->response[2][0],
          detector->response[2][1], detector->response[2][2]);
  printf("      { \"%s\",\n", detector->frDetector.name);
  printf("        %.15g, %.15g, %.15g,\n",
         detector->frDetector.vertexLongitudeDegrees,
         detector->frDetector.vertexLatitudeDegrees,
         detector->frDetector.vertexElevation);
  printf("        %.15e, %.15g,\n",
         detector->frDetector.xArmAltitudeRadians,
         detector->frDetector.xArmAzimuthRadians);
  printf("        %.15e, %.15g } }\n",
         detector->frDetector.yArmAltitudeRadians,
         detector->frDetector.yArmAzimuthRadians);
  return;
}


static int
CheckDetector(LALStatus *status, LALDetector *cachedDetector)
{
  LALDetector          constructedDetector;
  LALDetectorType      type = cachedDetector->type;
  LALFrDetector        frDetector = cachedDetector->frDetector;
  INT2                 i, j;

  printf("  Name is \"%s\",\n", frDetector.name);

  frDetector.name[1] -= 4;

  printf("  Changing to \"%s\",\n", frDetector.name);

  LALCreateDetector( status, &constructedDetector, &frDetector, type );
  TestStatus( status, CODES(0), DETECTORSITETESTC_EFLS );

  printf("  Enum'ed type is %d (LALDETECTORTYPE_IFODIFF=%d)\n",
         constructedDetector.type, LALDETECTORTYPE_IFODIFF);

  for (i=0; i<3; ++i)
  {
    printf("  x[%d]:\n    cached: %.15g calc: %.15g diff: %g\n", i+1,
           cachedDetector->location[i],
           constructedDetector.location[i],
           cachedDetector->location[i] - constructedDetector.location[i]);
    if ( cachedDetector->location[i] - constructedDetector.location[i]
         >= DETECTORSITETESTC_LOCTOL )
    {
      return DETECTORSITETESTC_EFLS;
    }
  }

  for (i=0; i<3; ++i)
  {
    for (j=0; j<3; ++j)
    {
      printf("  d[%d,%d]:\n    cached: %g calc: %g diff: %g\n", i+1, j+1,
             cachedDetector->response[i][j],
             constructedDetector.response[i][j],
             cachedDetector->response[i][j]
             - constructedDetector.response[i][j]);
      if ( cachedDetector->response[i][j]
           - constructedDetector.response[i][j] >= DETECTORSITETESTC_RESTOL )
      {
        return DETECTORSITETESTC_EFLS;
      }
    }
  }

  printf("  Latitude:\n    cached: %g calc: %g diff: %g\n",
         cachedDetector->frDetector.vertexLatitudeDegrees,
         constructedDetector.frDetector.vertexLatitudeDegrees,
         cachedDetector->frDetector.vertexLatitudeDegrees -
         constructedDetector.frDetector.vertexLatitudeDegrees);

  if ( cachedDetector->frDetector.vertexLatitudeDegrees !=
         constructedDetector.frDetector.vertexLatitudeDegrees)
  {
    return DETECTORSITETESTC_EFLS;
  }

  printf("  Longitude:\n    cached: %g calc: %g diff: %g\n",
         cachedDetector->frDetector.vertexLongitudeDegrees,
         constructedDetector.frDetector.vertexLongitudeDegrees,
         cachedDetector->frDetector.vertexLongitudeDegrees -
         constructedDetector.frDetector.vertexLongitudeDegrees);

  if ( cachedDetector->frDetector.vertexLongitudeDegrees !=
         constructedDetector.frDetector.vertexLongitudeDegrees)
  {
    return DETECTORSITETESTC_EFLS;
  }

  printf("  X Arm altitide:\n    cached: %g calc: %g diff: %g\n",
         cachedDetector->frDetector.xArmAltitudeRadians,
         constructedDetector.frDetector.xArmAltitudeRadians,
         cachedDetector->frDetector.xArmAltitudeRadians -
         constructedDetector.frDetector.xArmAltitudeRadians);

  if ( cachedDetector->frDetector.xArmAltitudeRadians !=
         constructedDetector.frDetector.xArmAltitudeRadians)
  {
    return DETECTORSITETESTC_EFLS;
  }

  printf("  X Arm azimuth:\n    cached: %g calc: %g diff: %g\n",
         cachedDetector->frDetector.xArmAzimuthRadians,
         constructedDetector.frDetector.xArmAzimuthRadians,
         cachedDetector->frDetector.xArmAzimuthRadians -
         constructedDetector.frDetector.xArmAzimuthRadians);

  if ( cachedDetector->frDetector.xArmAzimuthRadians !=
         constructedDetector.frDetector.xArmAzimuthRadians)
  {
    return DETECTORSITETESTC_EFLS;
  }

  printf("  Y Arm altitide:\n    cached: %g calc: %g diff: %g\n",
         cachedDetector->frDetector.yArmAltitudeRadians,
         constructedDetector.frDetector.yArmAltitudeRadians,
         cachedDetector->frDetector.yArmAltitudeRadians -
         constructedDetector.frDetector.yArmAltitudeRadians);
  if ( cachedDetector->frDetector.yArmAltitudeRadians !=
         constructedDetector.frDetector.yArmAltitudeRadians)
  {
    return DETECTORSITETESTC_EFLS;
  }

  printf("  Y Arm azimuth:\n    cached: %g calc: %g diff: %g\n",
         cachedDetector->frDetector.yArmAzimuthRadians,
         constructedDetector.frDetector.yArmAzimuthRadians,
         cachedDetector->frDetector.yArmAzimuthRadians -
         constructedDetector.frDetector.yArmAzimuthRadians);

  if ( cachedDetector->frDetector.yArmAzimuthRadians !=
         constructedDetector.frDetector.yArmAzimuthRadians)
  {
    return DETECTORSITETESTC_EFLS;
  }

  return DETECTORSITETESTC_ENOM;
}

/* The main function */
int main( int argc, char *argv[] )
{
  static LALStatus     status;
  int                  code;

  /*
  LALFrDetector testFrDetector = {"Test Site",
                                  0.0, 0.0, 0.0,
                                  0.0, 0.0,
                                  0.0, LAL_PI_2 };
  LALFrDetector lhoFrDetector
    = {"LIGO Hanford Observatory",
       -119.40765714,   46.4551467,           142.544,
         -6.195e-4,      2.199104,
          1.25e-5,       3.769901 };
  LALFrDetector lloFrDetector
    = {"LIGO Livingston Observatory",
        -90.77424039,   30.56289433,           -6.574,
         -3.121e-4,      3.4508039,
         -6.107e-4,      5.021600 };
  */
  /*
  LALDetector          testDetector
    = { {  4510094.634714941322397L,   4510094.634714941322397L,      0.0L },
      { {  0.0,        0.0,        0.0       },
        {  0.0,        0.0,        0.0       },
        {  0.0,        0.0,        0.0       }
      },
      LALDETECTORTYPE_IFODIFF,
      { "Test Interferometer",
          45.0L,          0.0L,           100.0,
           0.0,           0.0,
           0.0,           0.0
      }
    };
    = { {  4517590.8789357600195508L,   0.0L,    4487348.40860601410662157L },
      { {  0.0,        0.0,        0.0       },
        {  0.0,        0.0,        0.0       },
        {  0.0,        0.0,        0.0       }
      },
      LALDETECTORTYPE_IFODIFF,
      { "Test Interferometer",
           0.0L,         45.0L,           0.0,
           0.0,           0.0,
           0.0,           0.0
      }
    };
  */

  LALDetector          cachedDetector;

  ParseOptions( argc, argv );

  printf("Checking LHO...\n");
  cachedDetector = lalCachedDetectors[LALDetectorIndexLHODIFF];

  if ( code = CheckDetector( &status, &cachedDetector ) )
  {
    return code;
  }

  printf("Checking LLO...\n");
  cachedDetector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if ( code = CheckDetector( &status, &cachedDetector ) )
  {
    return code;
  }

  printf("Checking VIRGO...\n");
  cachedDetector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  if ( code = CheckDetector( &status, &cachedDetector ) )
  {
    return code;
  }

  printf("Checking GEO600...\n");
  cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if ( code = CheckDetector( &status, &cachedDetector ) )
  {
    return code;
  }


  printf("Checking TAMA300...\n");
  cachedDetector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  if ( code = CheckDetector( &status, &cachedDetector ) )
  {
    return code;
  }

  printf("Checking CIT40...\n");
  cachedDetector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  if ( code = CheckDetector( &status, &cachedDetector ) )
  {
    return code;
  }

  /*
  printf("Checking trivial detector...\n");

  if ( code = CheckDetector( &status, &testDetector ) )
  {
    return code;
  }

  */
  /*

  PrintLALDetector(&lalCachedDetector[LALDetectorIndexLHODIFF]);
  PrintLALDetector(&lalCachedDetector[LALDetectorIndexLLODIFF]);

  LALCreateDetector( &status, &lalDetector, &testFrDetector, &type );
  PrintLALDetector( &lalDetector );

  LALCreateDetector( &status, &lalDetector, &lloFrDetector, &type );
  PrintLALDetector( &lalDetector );

  LALCreateDetector( &status, &lalDetector, &lhoFrDetector, &type );
  PrintLALDetector( &lalDetector );

  printf("%15g %15g %15g\n\n", lalDetector.location[0],
         lalDetector.location[1], lalDetector.location[2]);

  for (i=0; i<3; ++i)
  {
    printf("%15g %15g %15g\n", lalDetector.response[i][0],
           lalDetector.response[i][1], lalDetector.response[i][2]);
  }
  */

  LALCheckMemoryLeaks();

  printf("PASS: All tests\n");

  return DETECTORSITETESTC_ENOM;
}


/*----------------------------------------------------------------------*/


/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (LALStatus *status, const char *ignored, int exitcode)
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (LALStatus *status)
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}


/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -q         quiet: run silently\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  fprintf (stderr, "  -d level   set lalDebugLevel to level\n");
  exit (exitcode);
}


/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}

