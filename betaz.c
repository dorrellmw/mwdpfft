#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "betaz.h"

static struct betaz * newBetaZ() {
  struct betaz * bz = malloc(sizeof(struct betaz));
  bz->maxbin = 0;
  bz->dz = 0;
  bz->nframes = 0;
  bz->beta = malloc(sizeof(double));
  bz->beta[0]=0;
  return bz;
}

void * cRealloc(void * old,int oldmax,int newmax,size_t elsize) {
  uintptr_t newbase = (uintptr_t) malloc((2*newmax+1)*elsize);
  uintptr_t new = newbase + newmax*elsize;
  uintptr_t oldbase = ((uintptr_t) old) - oldmax*elsize;
  if(newmax > oldmax) {
    int diff = newmax - oldmax;
    uintptr_t offsetnew = newbase + diff*elsize;
    memcpy((void *) offsetnew,(void *) oldbase,(2*oldmax+1)*elsize);
  } else if (oldmax > newmax) {
    int diff = oldmax - newmax;
    uintptr_t offsetold = oldbase + diff*elsize;
    memcpy((void *) newbase, (void *) offsetold,(2*newmax+1)*elsize);
  } else {
    memcpy((void *) newbase, (void *) oldbase,(2*newmax+1)*elsize);
  }
  free((void *) oldbase);
  return (void *) new;
}

//negative edge is considered in-bounds, positive edge is out-of-bounds
double rewrap(double a, double L) {
  while(a>=L/2)
    a-=L;
  while(a<-L/2)
    a+=L;
  return a;
}

//round toward zero
int rnd(double a) {
  return (a>0)?ceil(a-0.5):floor(a+0.5);
}

struct betaz * runAvgBetaZ(double * bs, float * zs, int natom, double A, double Lz, double dz, struct betaz * bz) {
  if(!bz)
    bz = newBetaZ();
  
  if(bz->nframes==0) {
    bz->dz = dz;
    //floor here because we only want "full" bins
    bz->maxbin = floor((Lz/bz->dz - 1)/2);
    bz->beta = cRealloc(bz->beta, 0, bz->maxbin, sizeof(double));
    for (int i=-bz->maxbin; i<=bz->maxbin; i++) {
      bz->beta[i]=0;
    }
  } else {
    int maxbin = floor((Lz/bz->dz - 1)/2);
    if(maxbin < bz->maxbin) {
      bz->beta = cRealloc(bz->beta, bz->maxbin, maxbin, sizeof(double));
      bz->maxbin = maxbin;
    }
  }

  //Just because I'm lazy.
  int maxbin = bz->maxbin;
  dz = bz->dz;
  double * beta = bz->beta;
  
  double * beta1f = calloc((2*maxbin+1),sizeof(double));
  double * beta1 = &(beta1f[maxbin]);

  for (int i=0; i<natom; i++) {
    int bin = rnd(zs[i]/dz);
    if(abs(bin) <= maxbin)
      beta1[bin]+=bs[i];
  }
  for (int i=-maxbin; i<=maxbin; i++) {
    beta1[i]/=(A*dz);
    beta[i] += (beta1[i]-beta[i])/(bz->nframes + 1);
  }
  bz->nframes++;

  free(beta1f);

  return bz;
}

void freeBetaZ(struct betaz * bz) {
  if(bz->beta) {
    free(bz->beta);
    bz->beta = NULL;
  }
  free(bz);
}
