#include <stdio.h>


/*Includes for UserTimeVal()*/
#include <limits.h>
#include <sys/times.h>

/*includes for WallTimeVal()*/
#include <sys/time.h>

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
   /*returns the time difference in seconds between calls. 
    It should be used like this:
      WallTimeVal();
      FunctionToBeTimed(ARGS);
      TempTime = WallTime();
      printf("FunctionToBeTimed took %g seconds\n");
   */
  static struct timeval last;
  struct timeval now;
  double val;
  int Check;

  Check = gettimeofday(&now, 0);

  if (Check < 0) {
    printf("There was a problem. Check = %d\n", Check);
  } 

  val = (double) (now.tv_sec - last.tv_sec) + (double)(now.tv_usec - last.tv_usec)/1000000.0;

  last.tv_sec = now.tv_sec;
  last.tv_usec = now.tv_usec;

  return(val);
}
