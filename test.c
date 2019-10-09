#include "stdlib.h"

#include "dcd.h"
#include "psf.h"
#include "betaz.h"

void dumpBetaZ(struct betaz * bz) {
  printf("z,beta(z)\n");
  for (int i=-bz->maxbin; i<=bz->maxbin; i++)
    printf("%lf,%lf\n",0.5*i,bz->beta[i]);
}

int main(int argc, char *argv[]) {

  //get scattering lengths
  struct psf p = readPSF(argv[1]);
  int natom=p.natom;
  double bs[natom];
  for(int i=0; i<natom; i++)
    bs[i]=p.atoms[i].b;
  freePSF(p);

  //get zs
  struct dcd * d = openDCD(argv[2]);
  double uc[6];
  float xs[natom];
  float ys[natom];
  float zs[natom];

  struct betaz * bz = NULL;

  for (int i=0; i<getNFrames(d); i++) {
    getUnitCell(d,uc);
    getCoords(d,xs,ys,zs);
    bz = runAvgBetaZ(bs,zs,natom, uc[0]*uc[1], uc[2], 0.5, bz);
    nextFrame(d);
  }

  closeDCD(d);

  dumpBetaZ(bz);

  free(bz);
}
