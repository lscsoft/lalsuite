/*
*  Copyright (C) 2007 Badri Krishnan
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
 * \file
 * \ingroup pulsarApps
 */

#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
int main(int argc, char *argv[]){
  
  int i, arg;
  float patchlist[23][4] = {
    {0.0000,  1.5708,  0.9020,  0.9020},
    {0.0000,  0.7122,  0.8976,  0.8151},
    {0.8976,  0.7122,  0.8976,  0.8151}, 
    {1.7952,  0.7122,  0.8976,  0.8151},
    {2.6928,  0.7122,  0.8976,  0.8151},
    {3.5904,  0.7122,  0.8976,  0.8151},
    {4.4880,  0.7122,  0.8976,  0.8151},
    {5.3856,  0.7122,  0.8976,  0.8151},
    {0.0000,  0.0000,  0.8976,  0.6094},
    {0.8976,  0.0000,  0.8976,  0.6094},
    {1.7952,  0.0000,  0.8976,  0.6094},
    {2.6928,  0.0000,  0.8976,  0.6094},
    {3.5904,  0.0000,  0.8976,  0.6094},
    {4.4880,  0.0000,  0.8976,  0.6094},
    {5.3856,  0.0000,  0.8976,  0.6094},
    {0.0000, -0.7122,  0.8976,  0.8151},
    {0.8976, -0.7122,  0.8976,  0.8151},
    {1.7952, -0.7122,  0.8976,  0.8151},
    {2.6928, -0.7122,  0.8976,  0.8151},
    {3.5904, -0.7122,  0.8976,  0.8151},
    {4.4880, -0.7122,  0.8976,  0.8151},
    {5.3856, -0.7122,  0.8976,  0.8151},
    {0.0000, -1.5708,  0.9020,  0.9020},
  };

  i = 0; 
  if (argc > 1) i = atoi( argv[argc-1] );
  if ((i>=1)&&(i<=23)) printf("%f  %f  %f  %f\n", patchlist[i-1][0], patchlist[i-1][1], patchlist[i-1][2], patchlist[i-1][3]);
 
  return 0;
}








