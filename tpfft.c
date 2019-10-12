#include <stdlib.h>
#include <string.h>

#include "dcd.h"
#include "psf.h"
#include "betaz.h"
#include "pfft.h"
#include "readArray.h"

void dumpRes(struct result res) {
  printf("q,I(q)\n");
  for (int i=0; i<res.nqs; i++)
    printf("%lf,%lf\n",res.qs[i],res.is[i]);
}

int main(int argc, char *argv[]) {

  //get scattering lengths
  fprintf(stderr,"Reading structure file.\n");
  struct psf p = readPSF(argv[1]);
  int natom=p.natom;
  double bs[natom];
  for(int i=0; i<natom; i++)
    bs[i]=p.atoms[i].b;
  freePSF(p);

  //prep coords arrays
  fprintf(stderr,"Opening DCD.\n");
  struct dcd * d = openDCD(argv[2]);
  double uc[6];
  float xs[natom];
  float ys[natom];
  float zs[natom];

  struct betaz * bz = NULL;

  //build beta(z)
  fprintf(stderr,"Generating beta(z)\n");
  for (int i=0; i<getNFrames(d); i++) {
    getUnitCell(d,uc);
    getCoords(d,xs,ys,zs);
    bz = runAvgBetaZ(bs,zs,natom, uc[0]*uc[1], uc[2], 0.5, bz);
    nextFrame(d);
  }

  //get solvent NSLD from beta(z)
  double bw=0;
  if(bz->maxbin>30)
    bw = (bz->beta[-bz->maxbin] +
                 bz->beta[-bz->maxbin + 1] +
                 bz->beta[-bz->maxbin + 2] +
                 bz->beta[-bz->maxbin + 3] +
                 bz->beta[ bz->maxbin - 3] +
                 bz->beta[ bz->maxbin - 2] +
                 bz->beta[ bz->maxbin - 1] +
                 bz->beta[ bz->maxbin])/8.0;
  else
    bw = (bz->beta[-bz->maxbin] + bz->beta[ bz->maxbin])/2.0;

  goToFrame(d,0);

  //get qs
  fprintf(stderr,"Reading qs\n");
  double * qs = malloc(1024*1024*sizeof(double));
  int nqs = readDoubles(argv[3],qs);
  while(qs[nqs-1]==0)
    nqs--;
  qs = realloc(qs,nqs*sizeof(double));
  double * is = calloc(nqs,sizeof(double));
  double * tis = calloc(nqs,sizeof(double));

  struct partsys ps = {.natom = natom, .xs = xs, .ys = ys, .zs = zs, .bs = bs, .uc = uc};
  struct result res = {.nqs = nqs, .qs = qs, .is = tis};

  double ds = 1; //1 angstom bin size
  int rci = 30; //30 angstom cutoff
  int simax = 199; //200 angstom limit

  fprintf(stderr,"Running pfft.\n");
  for (int i=0; i<getNFrames(d); i++) {
    getUnitCell(d,uc);
    getCoords(d,xs,ys,zs);
    memset(res.is,0,res.nqs*sizeof(double));
    pfft(ps,bz,bw,rci,ds,simax,res);
    for(int qi=0; qi<res.nqs; qi++) {
      is[qi]+= (res.is[qi]-is[qi])/(i + 1);
    }
  }
  free(tis);

  dumpRes(res);

  closeDCD(d);
  free(qs);
  free(is);
  freeBetaZ(bz);
}
