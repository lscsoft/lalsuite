#include <stdio.h>


/*Includes for UserTimeVal()*/
#include <limits.h>
#include <sys/times.h>

/*includes for WallTimeVal()*/
#include <time.h>
#include <sys/time.h>

double UserTimeVal(void);
double WallTimeVal(void);

double UserTimeVal(void) {
  /*returns the user time difference in seconds between calls. 
    It should be used like this:
      UserTimeVal();
      FunctionToBeTimed(ARGS);
      TempTime = UserTime();
      printf("FunctionToBeTimed took %g user seconds\n");
    Please note that this function can only deal with one wrap around.
    Wrap arounds occur about every 35 user minutes.*/

  static clock_t last_usertime = 0;
  struct tms now;
  double val;
  clock_t Check;

  
  Check = times(&now);

  if (Check < 0) {
    printf("There was a problem. Check = %ld\n", Check);
  }
 
  if (now.tms_utime >= last_usertime) {

    val = (now.tms_utime - last_usertime) / (double)CLOCKS_PER_SEC;

  } else if (now.tms_utime < last_usertime) {
 
   printf("Wrapped around. last_usertime = %ld, now.tms_utime = %ld, diff = %ld\n", \
	   last_usertime, now.tms_utime, now.tms_utime + (LONG_MAX - last_usertime));

    val = (now.tms_utime + (LONG_MAX - last_usertime))/(double)CLOCKS_PER_SEC;

  }
  if (val < 0){
    printf("now.tms_utime = %ld, last_usertime = %ld, val = %g\n", \
	   last_usertime, now.tms_utime, val);
  }

  last_usertime = now.tms_utime;
  
  return(val);
}

double WallTimeVal(void) {
  /*returns the user time difference in seconds between calls. 
    It should be used like this:
      WallTimeVal();
      FunctionToBeTimed(ARGS);
      TempTime = WallTime();
      printf("FunctionToBeTimed took %g user seconds\n");
    Please note that this function can only deal with one wrap around.
    Wrap arounds occur about every 35 user minutes.*/

  static clock_t last_walltime = 0;
  static clock_t now_walltime = 0;
  struct tms now;
  double val;
  clock_t Check;
  
  
  Check = times(&now);
 
  if (Check < 0) {
    printf("There was a problem. Check = %ld\n", Check);
  }
  now_walltime = now.tms_utime + now.tms_stime;
 
  if (now_walltime >= last_walltime) {

    val = (now_walltime - last_walltime) / (double)CLOCKS_PER_SEC;
    
  } else if (now_walltime < last_walltime) {
    
    printf("Wrapped around. last_walltime = %ld, now_walltime = %ld, diff = %ld\n", \
           last_walltime, now_walltime, now_walltime + (LONG_MAX - last_walltime));
    
    val = (now_walltime + (LONG_MAX - last_walltime))/(double)CLOCKS_PER_SEC;
    
  }
  if (val < 0){
    printf("now_walltime = %ld, last_walltime = %ld, val = %g\n", \
           last_walltime, now_walltime, val);
  }
  
  last_walltime = now_walltime;
  
  return(val);
}

