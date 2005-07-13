
/* B. Krishnan, July 2005 */
/* small modification by M.A. Papa */

#include<stdio.h>
#include<math.h>




int main( int argc, char *argv[]){

  int r;
  double temp1, temp2, temp3, temp4, temp5, temp6; 
  double M11, M12, M21, M22, Y1, Y2, a, b, det;
  double sigma, sigma2, h, C; 
  FILE *fp=NULL;

  fp = fopen( argv[1], "r");

  M11 = 0.0;
  M12 = 0.0;
  M21 = 0.0;
  M22 = 0.0;
  Y1 = 0.0;
  Y2 = 0.0;
  
  do{
    r=fscanf(fp,"%lf%lf%lf%lf%lf\n", &temp1, &sigma, &h, &temp4, &C); 
    if (r==5) {
      sigma2 = sigma * sigma;
      M11 += h/sigma2;
      M12 += 1.0/sigma2;
      M21 += h * h/sigma2;
      M22 += h/sigma2;
      Y1 += C/sigma2;
      Y2 += h * C /sigma2;
    }
  } while (r != EOF);

  fclose(fp);

  det = M11 * M22 - M12 * M21;
  b = (M11 * Y2 - M21 * Y1)/det;
  a = (M22 * Y1 - M12 * Y2)/det;

  fprintf(stdout, "%11.7E %11.7E\n", a*1E-19, (0.95 -b)*1E-19/a);

  return 0;
}

