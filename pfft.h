#ifndef PFFT_H_
#define PFFT_H_

#include "betaz.h"

struct partsys {
  int natom;
  float * xs;
  float * ys;
  float * zs;
  double * bs;
  double * uc;
};

struct result {
  int nqs;
  double * qs;
  double * is;
};

//void term1a(struct partsys ps, double rc, struct result res);
//void term1b(struct partsys ps, struct betaz * bz, int rci, double ds, int simax, struct result res);
//void term234(struct betaz * bz, double bw, double ds, int simax, struct result res);
//void termHighS(struct partsys ps, double bw, double smax, struct result res);
void pfft(struct partsys ps, struct betaz * bz, double bw, int rci, double ds, int simax, struct result res);

#endif
