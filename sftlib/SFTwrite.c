/* $Id$ */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SFTReferenceLibrary.h"

int main() {
  FILE *fp;
  float data[]={1.0, 0.0, 2.0, -1.0, 3.0, -2.0, 4.0, -3.0};
  char comment[]="Hello world";
  int err=0;
  
  fp=fopen("SFT-test", "w");
  if ((err=WriteSFT(fp, 12345, 6789, 60, 1000, sizeof(data)/(2*sizeof(float)), comment, data)))
    fprintf(stderr, "WriteSFT failed with return value %d\n", err);
  fclose(fp);
  return err;
}
