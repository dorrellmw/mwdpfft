#ifndef PSF
#define PSF

struct psfatom {
  char seg[9];
  char resid[9];
  char res[9];
  char name[9];
  char type[7];
  double charge;
  double mass;
  int imove; // 0 if free, 1 if fixed, -1 if lone-pair
  double ech; // CHEQ electronegativity
  double eha; // CHEQ hardness
  double b; // scattering length
};

struct bond {
  int a;
  int b;
};

struct angle {
  int a;
  int b;
  int c;
};

struct dihedral {
  int a;
  int b;
  int c;
  int d;
};

struct psfsig {
  int valid;
  int ext;
  int cmapcheq;
  int xplor;
  int slb;
};

struct psf {
  struct psfsig sig;
  int ntitle;
  int natom;
  int nbond;
  int ntheta;
  int nphi;
  int nimphi;

  char ** titles;
  struct psfatom * atoms;
  struct bond * bonds;
  struct angle * angles;
  struct dihedral * dihedrals;
  struct dihedral * impropers;
};

struct psf readPSF(const char * path);
void freePSF(struct psf p);
#endif
