#include <stdio.h>
#include <string.h>

#define length 16
                                                                               
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
