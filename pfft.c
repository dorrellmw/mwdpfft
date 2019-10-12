#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "pfft.h"

double pi = acos(-1);

double sinc(double x) {
  return sin(x)/x;
}

struct pair {
  int a;
  float b;
};

void zip(int * as, float * bs, struct pair * ps, int n) {
  for(int i=0; i<n; i++) {
    ps[i].a=as[i];
    ps[i].b=bs[i];
  }
}

int comp(const void * a, const void * b) {
  if (((struct pair *)a)->b < ((struct pair *)b)->b) return -1;
  if (((struct pair *)a)->b == ((struct pair *)b)->b) return 0;
  if (((struct pair *)a)->b > ((struct pair *)b)->b) return 1;
  //return ((struct pair *)a)->b - ((struct pair *)b)->b;
}

int nearby1D(struct pair * ps,int natom,double L,int xi, float x,double rc,int * n) {
  struct pair p = {.a = 0, .b = x};
  struct pair * ptr = bsearch(&p,ps,natom,sizeof(struct pair),comp);
  uintptr_t pi = ptr-ps;
  int count = 0;
  int j;
  for(j=pi+1; j<natom; j++) {
    if(ps[j].b - x < rc) {
      n[count]=ps[j].a;
      count++;
    } else break;
  }
  if(j==natom) {
    for(j=0; j<pi; j++) {
      if(ps[j].b + L - x < rc) {
        n[count]=ps[j].a;
        count++;
      } else break;
    }
  }
  //for(j=pi-1; j>=0; j--) {
  //  if(x - ps[j].b < rc) {
  //    n[count]=ps[j].a;
  //    count++;
  //  } else break;
  //}
  //if(j==-1) {
  //  for(j=natom-1; j>pi; j--) {
  //    if(x - ps[j].b + L < rc) {
  //      n[count]=ps[j].a;
  //      count++;
  //    } else break;
  //  }
  //}
  return count;
}

int compi(const void * a, const void * b) {
  return (*(int *)a) - (*(int *)b);
}

int intersect(int * a, int na, int * b, int nb, int * n) {
  qsort(a,na,sizeof(int),compi);
  qsort(b,nb,sizeof(int),compi);
  int count = 0;
  int ai=0;
  int bi=0;
  while(ai<na && bi<nb) {
    if(a[ai]==b[bi]) {
      n[count]=a[ai];
      count++;
      ai++;
      bi++;
      continue;
    } else if (a[ai]>b[bi]) {
      bi++;
    } else {
      ai++;
    }
  }
  return count;
}

//rc = ds*rci
void term1a(struct partsys ps, double rc, struct result res) {
  int ns[ps.natom];
  for(int i=0; i<ps.natom; i++)
    ns[i]=i;
  struct pair * srtx = malloc(ps.natom*sizeof(struct pair));
  struct pair * srty = malloc(ps.natom*sizeof(struct pair));
  zip(ns,ps.xs,srtx,ps.natom);
  zip(ns,ps.ys,srty,ps.natom);
  qsort(srtx,ps.natom,sizeof(struct pair),comp);
  qsort(srty,ps.natom,sizeof(struct pair),comp);
  #pragma omp parallel for
  for(int qi=0; qi<res.nqs; qi++) {
    for(int i=0; i<ps.natom; i++) {
      int nx[ps.natom];
      int nnx = nearby1D(srtx,ps.natom,ps.uc[0],i,ps.xs[i],rc,nx);
      int ny[ps.natom];
      int nny = nearby1D(srty,ps.natom,ps.uc[1],i,ps.ys[i],rc,ny);
      int jis[ps.natom];
      int nj = intersect(nx, nnx, ny, nny, jis);
      if(i%1480==0)
        fprintf(stderr,"%d",i/1480);
      for(int k=0; k<nj; k++) {
        int j = jis[k];
        double dx = ps.xs[j] - ps.xs[i];
        double dy = ps.ys[j] - ps.ys[i];
        if(dx*dx+dy*dy<rc*rc) {
          double dz = ps.zs[j] - ps.zs[i];
          double dr = sqrt(dx*dx + dy*dy + dz*dz);
          res.is[qi]+=2*ps.bs[i]*ps.bs[j]*sinc(res.qs[qi]*dr);
        }
      }
    }
  }
  
  free(srtx);
  free(srty);
}

//rci is the index of the first ring OUTSIDE the cutoff
void term1b(struct partsys ps, struct betaz * bz, double bw, int rci, double ds, int simax, struct result res) {
  #pragma omp parallel for
  for (int qi=0; qi<res.nqs; qi++) {
    for(int i=0; i<ps.natom; i++) {
      if(i%1480==0) fprintf(stderr,"%d",i/1480);
      double base = ps.bs[i]*ds*2*pi*(bz->dz);
      for (int si=rci; si<=simax; si++) {
        double s = ds*0.5+ds*si;
        for (int zi=-bz->maxbin; zi<=bz->maxbin; zi++) {
          double deltaz = zi*(bz->dz) - ps.zs[i];
          //res.is[qi]+=ps.bs[i]*(bz->dz)*ds*2*pi*s*(bz->beta[zi])*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
          res.is[qi]+=base*s*(bz->beta[zi])*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
        //end caps
        {
          //top
          double z = ((bz->maxbin + 0.5)*(bz->dz) + ps.uc[2]/2.0)/2.0;
          double dz = ps.uc[2]/2.0 - (bz->maxbin + 0.5)*(bz->dz);
          double deltaz = z - ps.zs[i];
          res.is[qi]+=ps.bs[i]*dz*ds*2*pi*s*bw*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
          //bottom
          z *= -1;
          deltaz = z - ps.zs[i];
          res.is[qi]+=ps.bs[i]*dz*ds*2*pi*s*bw*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
      }
    }
  }
}

void term234(struct betaz * bz, double * uc, double bw, double ds, int simax, struct result res) {
  for(int qi=0; qi<res.nqs; qi++) {
    for(int si=0; si<=simax; si++) {
      double s = ds*0.5+ds*si;
      for(int zi=-bz->maxbin; zi<=bz->maxbin; zi++) {
        for(int zj=-bz->maxbin; zj<=bz->maxbin; zj++) {
          double deltaz = (zi-zj)*bz->dz;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*bz->dz*bz->dz*ds*(-2*bz->beta[zi]+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
        //end caps (j)
        {
          //top (j)
          double z2 = ((bz->maxbin + 0.5)*(bz->dz) + uc[2]/2.0)/2.0;
          double dz2 = uc[2]/2.0 - (bz->maxbin + 0.5)*(bz->dz);
          double deltaz = zi*bz->dz - z2;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*bz->dz*dz2*ds*(-2*bz->beta[zi]+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
          //bottom (j)
          z2 *= -1;
          deltaz = zi*bz->dz - z2;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*bz->dz*dz2*ds*(-2*bz->beta[zi]+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
      }
      //end caps (i)
      {
        //top (i)
        double z1 = ((bz->maxbin + 0.5)*(bz->dz) + uc[2]/2.0)/2.0;
        double dz1 = uc[2]/2.0 - (bz->maxbin + 0.5)*(bz->dz);
        for(int zj=-bz->maxbin; zj<=bz->maxbin; zj++) {
          double deltaz = z1 - zj*bz->dz;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*dz1*bz->dz*ds*(-2*bw+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
        //end caps (j)
        {
          //top (j)
          double z2 = ((bz->maxbin + 0.5)*(bz->dz) + uc[2]/2.0)/2.0;
          double dz2 = uc[2]/2.0 - (bz->maxbin + 0.5)*(bz->dz);
          double deltaz = z1 - z2;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*dz1*dz2*ds*(-2*bw+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
          //bottom (j)
          z2 *= -1;
          deltaz = z1 - z2;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*dz1*dz2*ds*(-2*bw+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
        //bottom (i)
        z1 *= -1;
        for(int zj=-bz->maxbin; zj<=bz->maxbin; zj++) {
          double deltaz = z1 - zj*bz->dz;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*dz1*bz->dz*ds*(-2*bw+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
        //end caps (j)
        {
          //top (j)
          double z2 = ((bz->maxbin + 0.5)*(bz->dz) + uc[2]/2.0)/2.0;
          double dz2 = uc[2]/2.0 - (bz->maxbin + 0.5)*(bz->dz);
          double deltaz = z1 - z2;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*dz1*dz2*ds*(-2*bw+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
          //bottom (j)
          z2 *= -1;
          deltaz = z1 - z2;
          res.is[qi]+=2*pi*s*bw*uc[0]*uc[1]*dz1*dz2*ds*(-2*bw+bw)*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
        }
      }
    }
  }
}

//smax = (simax+1)*ds
void termHighS(struct partsys ps, double bw, double smax, struct result res) {
  double sum=0;
  for(int i=0; i<ps.natom; i++)
    sum+=ps.bs[i];
  for(int qi=0; qi<res.nqs; qi++) {
    res.is[qi]+=2*pi*(sum - ps.uc[0]*ps.uc[1]*ps.uc[2]*bw)*(sum-ps.uc[0]*ps.uc[1]*ps.uc[2]*bw)*cos(res.qs[qi]*smax)/(res.qs[qi]*res.qs[qi]*ps.uc[0]*ps.uc[1]);
  }
}

void pfft(struct partsys ps, struct betaz * bz, double bw, int rci, double ds, int simax, struct result res) {
  fprintf(stderr,"a");
  term1a(ps, ds*rci, res);
  fprintf(stderr,"b");
  term1b(ps, bz, bw, rci, ds, simax, res);
  fprintf(stderr,"2");
  term234(bz, ps.uc, bw, ds, simax, res);
  fprintf(stderr,"H");
  termHighS(ps, bw, (simax+1)*ds, res);
}
