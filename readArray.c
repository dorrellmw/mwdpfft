#include <stdio.h>
#include <stdlib.h>

#include "readArray.h"

int readDoubles(char * path,double * xs) {
  FILE * f = fopen(path,"r");
  if(!f)
    return -1;
  int num=0;
  char buffer[1024];
  char * ptr;
  while(1) {
    ptr = fgets(buffer,1024,f);
    if(ferror(f) || feof(f) || ptr==NULL)
      break;

    xs[num]=atof(buffer);
    num++;
  }
  fclose(f);
  return num;
}
