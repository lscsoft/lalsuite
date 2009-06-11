/*
*  Copyright (C) 2007 Cristina Valeria Torres, Jolien Creighton
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

/*
 *  Example program for the use of the the "TrackSearch() routine;
 *  New Maintainer: Torres, C (Univ TX at Browsville)
 *  Original Author:  R. Balasubramanian
 */

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/TrackSearch.h>

NRCSID (TESTSTSC,"$Id$");

static void OutputInvertedMap(TimeFreqRep, TrackSearchParams);

INT4 lalDebugLevel=0;
int main( void )
{
  static LALStatus status;
  static TrackSearchOut out;
  static TimeFreqRep in;
  static TrackSearchParams params;
  INT4 i,j,cnt,maxVal;
  FILE *fp;
  unsigned char dummy;
  char stringd[256];


  /* set the parameters */
  params.sigma=2;  /* 2 */
  params.high = 1; /* 3.0 */
  params.low= 3; /* 1 */
  params.low = params.high/3; /* ? */
  /* open an input pgm file */
  fp = LALOpenDataFile("a.pgm");
  fgets(stringd,255,fp);
  fgets(stringd,255,fp);
  /* read the height and width of the image */
  fscanf(fp,"%d ",&params.height);
  fscanf(fp,"%d ",&params.width);
  fscanf(fp,"%d ",&maxVal);
  /* Allocate space for the input array */
  in.map=LALMalloc(params.height*sizeof(REAL4 *));
  for(i=0;i<params.height;i++)
    in.map[i] = LALMalloc(params.width*sizeof(REAL4));
  /* Read the image */
  for(j=0;j<params.width;j++){
    for(i=0;i<params.height;i++){
      fscanf(fp,"%c",&dummy);
      /* bright parts of the image should have higer values
	 hence the inversion */
      *(in.map[i] + j) = maxVal - dummy;
    }
  }
  /*
   * Use information about image to setup TF Rep tCol and fRow
   * values which must be
   * fRow/2+1 == width
   * tCol == height
   */
  in.fRow=(params.width-1)*2;
  in.tCol=params.height;
  /* Added static function to write out in.map structure */
  /* in is of type timefreqrep */
  printf("MaxVal = %d \n",maxVal);
  fclose(fp);
  /* set the allocFlag so that space can be allocated */
  params.allocFlag=1;
  /* Search for curves */
  LALSignalTrackSearch(&status, &out, &in, &params);
  REPORTSTATUS(&status);
  /* Output the details of the curves found.*/
    for(i=0;i<params.height;i++)
      for(j=0;j<params.width;j++)
	*(in.map[i] + j)=0;
  printf("number of curves %d\n",out.numberOfCurves);
  cnt=0;
  for(i=0;i<out.numberOfCurves;i++){
    if(out.curves[i].n>5){
      cnt++;
      printf(" curve number = %d, length =  %d, power = %4.2e\n",cnt,out.curves[i].n,out.curves[i].totalPower);
      for(j=0;j<out.curves[i].n;j++){
	printf("%d %d %4.4e\n",out.curves[i].row[j],out.curves[i].col[j],out.curves[i].depth[j]);
	*(in.map[out.curves[i].row[j]] + out.curves[i].col[j]) = cnt;
      }
    }
  }

  /* Free the space allocated for output structures */
  for(i=0;i<out.numberOfCurves;i++){
    LALFree(out.curves[i].row);
    LALFree(out.curves[i].col);
    LALFree(out.curves[i].depth); /*added*/
  }
  if(out.curves!=NULL)
    LALFree(out.curves);
  out.curves=NULL;
  /* set params.allocFlag=2 and call the routine again to free space*/
  params.allocFlag=2;
  LALSignalTrackSearch(&status, &out, &in, &params);
  REPORTSTATUS(&status);
  return(0);
}

