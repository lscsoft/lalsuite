/*
*  Copyright (C) 2007 Tania Regimbau
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

#include <stdio.h>
#include <string.h>

#define length 60
                                                                               
main()
 {

  int i, start, stop, N;
  FILE *pf1,*pf2;
  pf1 = fopen ("VirgoP1a.cache","w");
  pf2 = fopen ("LIGOP1a.cache","w");
  start = 782758813;
  stop = 782769553;
  N = ceil ((stop - start)/length) + 1;  
  for (i =0; i < N; i++)
    {
     
      fprintf(pf1,"V2 SIM %d 60 file://localhost/home/tania/lscsoft/src/lalapps/src/stochastic/data/PROJECT1a/simulated_virgo_data/SIM/V2-SIM-%d-60.gwf\n",start,start);

   
     fprintf(pf2,"H SIM %d 60 file://localhost/home/tania/lscsoft/src/lalapps/src/stochastic/data/PROJECT1a/simulated_ligo_noise/SIM/H-SIM-%d-60.gwf\n",start,start);
     start = start+60;
    
    }
  fclose(pf1);fclose(pf2);

 

 }  
