#include <math.h>
#include <stdio.h>

#include "pfft.h"

double pi = acos(-1);

double sinc(double x) {
  return sin(x)/x;
}

//rc = ds*rci
void term1a(struct partsys ps, double rc, struct result res) {
  for(int qi=0; qi<res.nqs; qi++) {
    for(int i=0; i<ps.natom; i++) {
      #pragma omp parallel for
      for(int j=0; j<ps.natom; j++) {
        double dx = ps.xs[j] - ps.xs[i];
        double dy = ps.ys[j] - ps.ys[i];
        if(dx*dx+dy*dy<rc*rc) {
          double dz = ps.zs[j] - ps.zs[i];
          double dr = sqrt(dx*dx + dy*dy + dz*dz);
          res.is[qi]+=ps.bs[i]*ps.bs[j]*sinc(res.qs[qi]*dr);
        }
      }
    }
  }
}

//rci is the index of the first ring OUTSIDE the cutoff
void term1b(struct partsys ps, struct betaz * bz, double bw, int rci, double ds, int simax, struct result res) {
  for(int i=0; i<ps.natom; i++) {
    for (int qi=0; qi<res.nqs; qi++) {
      for (int si=rci; si<=simax; si++) {
        double s = ds*0.5+ds*si;
        for (int zi=-bz->maxbin; zi<=bz->maxbin; zi++) {
          double deltaz = zi*(bz->dz) - ps.zs[i];
          res.is[qi]+=ps.bs[i]*(bz->dz)*ds*2*pi*s*(bz->beta[zi])*sinc(res.qs[qi]*sqrt(s*s+deltaz*deltaz));
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
