struct betaz {
  int maxbin;
  double dz;
  int nframes;
  double * beta;
};

struct betaz * runAvgBetaZ(double * bs, float * zs, int natom, double A, double Lz, double dz, struct betaz * bz);

void freeBetaZ(struct betaz * bz);
