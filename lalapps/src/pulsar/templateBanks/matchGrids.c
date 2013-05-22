/*
 * Copyright (C) 2005 Reinhard Prix
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

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup pulsarApps
 * \brief Match two "equivalent" (i.e. corresponding to same parameters) sky-grids 
 *  at a given tolerance. This is to only meaningful for sky-grids which are identical
 * up to small numerical errors plus allowing for 'outlier'-points, which don't match up.
 * Return two 'matched'-grids (match1.txt, match2.txt), i.e. the subsets of the two grids 
 * that match within the given tolerance , and two lists of the  respective 'skipped points'
 * (skip1.txt, skip2.txt).
 *
 * This was originally implemented as an octave-script, which turned out to be way too slow
 * for large grids of order 1e5 points.
 */


/*---------- INCLUDES ----------*/
#include <math.h>
#include <stdlib.h>

#include <lal/UserInput.h>

#include <lalapps.h>

#include <lal/DopplerScan.h>

/*---------- DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

#define myMAX(x,y) ( (x) > (y) ? (x) : (y) )

#define SQ(x) ((x)*(x))
#define modERR(x,y) ( SQ((x)->Alpha - (y)->Alpha) + SQ((x)->Delta - (y)->Delta) )

/*---------- User variables ---------- */
BOOLEAN uvar_help;
REAL8 uvar_tolerance;
CHAR *uvar_grid1;
CHAR *uvar_grid2;

/*---------- global vars ---------- */
extern int vrbflg;


/* ---------- prototypes ---------- */
int main(int argc, char *argv[]);

void initUserVars (LALStatus *status);
void loadSkyGridFile (LALStatus *, DopplerSkyGrid **grid, const CHAR *fname);

/*--------------------------------------------------*/
/* Test function(s) */
/*--------------------------------------------------*/
int main(int argc, char *argv[])
{
  UINT4 len1, len2;
  UINT4 i;
  UINT4 l1, l2;
  FILE *fpMatch1, *fpMatch2;
  FILE *fpSkip1, *fpSkip2;
  CHAR fnMatch1[512], fnMatch2[512];
  CHAR fnSkip1[512], fnSkip2[512];

  LALStatus status = blank_status;  

  DopplerSkyGrid *grid1 = NULL;
  DopplerSkyGrid *grid2 = NULL;

  DopplerSkyGrid *node1, *node2;

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  if ( uvar_help )
    {
      printf ("\nDescription: Outputs the matched subset-grids into 'match1.txt' and 'match2.txt',\n"
	      "the skipped points are output into 'skip1.txt' and 'skip2.txt' including\n"
	      "the original point-index in the first column.\n\n");

      return 0;
    }

  /* load grid1 */
  LAL_CALL ( loadSkyGridFile (&status, &grid1, uvar_grid1), &status);

  /* load grid2 */
  LAL_CALL ( loadSkyGridFile (&status, &grid2, uvar_grid2), &status);

  /* count lengths */
  len1 = 0;
  node1 = grid1;
  while (node1)
    {
      len1 ++;
      node1 = node1->next;
    }

  len2 = 0;
  node2 = grid2;
  while (node2)
    {
      len2 ++;
      node2 = node2->next;
    }

  /* open output-files */
  sprintf (fnMatch1, "match1.txt");
  sprintf (fnMatch2, "match2.txt");

  sprintf (fnSkip1, "skip1.txt");
  sprintf (fnSkip2, "skip2.txt");

  fpMatch1 = fopen(fnMatch1, "wb");
  fpMatch2 = fopen(fnMatch2, "wb");
  fpSkip1 = fopen(fnSkip1, "wb");
  fpSkip2 = fopen(fnSkip2, "wb");

  if ( !fpMatch1 || !fpMatch2 || !fpSkip1 || !fpSkip2 )
    {
      XLALPrintError("\nFailed to open of of ('%s','%s','%s','%s' for writing\n\n",
		    fnMatch1, fnMatch2, fnSkip1, fnSkip2);
      return 1;
    }


  /* step trough both grids: */
  l1 = 0;
  node1 = grid1;
  l2 = 0;
  node2 = grid2;
  for (i=0; i < myMAX(len1,len2); i++ )
    {
      BOOLEAN haveMatch = FALSE;

      /* case one: direct match */
      if ( modERR(node1,node2) < uvar_tolerance )
	haveMatch = TRUE;
      /* case two: 1-neighbor match */ 
      else if ( node1->next && (modERR(node1->next, node2) < uvar_tolerance ) )
	{
	  XLALPrintError ( "Skipped index=%d in '%s': [ %f, %f ]\n", 
			  l1, uvar_grid1, node1->Alpha, node1->Delta);
	  fprintf (fpSkip1, "%d %f %f\n", l1, node1->Alpha, node1->Delta );
	  node1 = node1->next;
	  l1 ++;
	  haveMatch = TRUE;
	}
      else if ( node2->next && (modERR(node1, node2->next) < uvar_tolerance) )
	{
	  XLALPrintError ( "Skipped index=%d in '%s': [ %f, %f ]\n", 
			  l2, uvar_grid2, node2->Alpha, node2->Delta);
	  fprintf (fpSkip2, "%d %f %f\n", l2, node2->Alpha, node2->Delta );
	  node2 = node2->next;
	  l2 ++;
	  haveMatch = TRUE;
	}
      /* case three: no neighbor-match found */
      else
	{
	  XLALPrintError ( "Skipped index=%d in '%s': [ %f, %f ]\n", 
			  l1, uvar_grid1, node1->Alpha, node1->Delta);
	  XLALPrintError ( "Skipped index=%d in '%s': [ %f, %f ]\n", 
			  l2, uvar_grid2, node2->Alpha, node2->Delta);


	  fprintf (fpSkip1, "%d %f %f\n", l1, node1->Alpha, node1->Delta );
	  fprintf (fpSkip2, "%d %f %f\n", l2, node2->Alpha, node2->Delta );

	  haveMatch = FALSE;
	}

      if ( haveMatch )
	{
	  fprintf (fpMatch1, "%f %f\n", node1->Alpha, node1->Delta);
	  fprintf (fpMatch2, "%f %f\n", node2->Alpha, node2->Delta);
	}

      node1 = node1->next;
      node2 = node2->next;
      l1 ++;
      l2 ++;

      if ( !node1 || !node2 )
	break;
      
    } /* for i < max(len1,len2) */

  if ( node1 )
    {
      XLALPrintError("\nHave %d leftover-nodes in '%s':\n", len1 - l1, uvar_grid1);
      while ( node1 )
	{
	  XLALPrintError ( "%d: [%f %f]\n", l1, node1->Alpha, node1->Delta );
	  fprintf (fpSkip1, "%d %f %f\n", l1, node1->Alpha, node1->Delta );
	  node1 = node1->next;
	  l1 ++;
	}
    }

  if ( node2 )
    {
      XLALPrintError("\nHave %d leftover-nodes in '%s':\n", len2 - l2, uvar_grid2);
      while ( node2 )
	{
	  XLALPrintError ( "%d: [%f %f]\n", l2, node2->Alpha, node2->Delta );
	  fprintf (fpSkip2, "%d %f %f\n", l2, node2->Alpha, node2->Delta );
	  node2 = node2->next;
	  l2 ++;
	}
    }

  /* done */
  fclose(fpMatch1);
  fclose(fpMatch2);
  fclose(fpSkip1);
  fclose(fpSkip2);

  return 0;
}

void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR ( status );

  uvar_tolerance = 1e-5;

  LALregBOOLUserVar(status,	help,		'h', UVAR_HELP, "Print this help/usage message");
  LALregREALUserVar(status,	tolerance,	'e', UVAR_OPTIONAL, "Tolerance in |error|^2 for grid-point matching");

  LALregSTRINGUserVar(status,	grid1,		'1', UVAR_REQUIRED, "Filename of grid 1");
  LALregSTRINGUserVar(status,	grid2,		'2', UVAR_REQUIRED, "Filename of grid 2");

  DETATCHSTATUSPTR ( status );
  RETURN(status);

} /* initUserVars() */
