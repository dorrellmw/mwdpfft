#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "psf.h"

/**
 * Read a count from a PSF section header.
 *
 * Retrieves the count from the header of a section of a PSF. This method does
 * not verify that the header field matches the actual number of elements
 * listed. If the file is unreadable, or the header field cannot be read, a
 * value of -1 is returned instead.
 *
 * @param[in] psf The PSF from which the count should be read.
 * @param[in] bang The PSF section title, starting with an exclamation point.
 * @return The number of elements in the section.
 */
static int readBang(FILE * psf, char * bang) {
  int nbang = -1; // in case bang count is not read successfully.
  char bufferA[1024]; // previous word
  char bufferB[1024]; // current word
  int numread = -1; // for detecting read errors
  while ( (!ferror(psf)) && (!feof(psf)) && numread!=0 ) { // stop looping upon
                                                           // error or eof
    if(!strncmp(bufferB,bang,strlen(bang))) { // watch for bang
      nbang = atoi(bufferA); // read word (number) before bang
      break; // stop looping
    }
    strcpy(bufferA,bufferB); // save current word as previous
    numread = fscanf(psf," %1024s ",bufferB); // load next word
  }
  return nbang;
}

/**
 * Read title lines from a PSF
 *
 * Reads title lines from a PSF file and stores them in the provided psf struct.
 * Also populates the ntitle field of the psf struct.
 *
 * @param[in] psf The PSF from which the title lines will be read.
 * @param[in,out] p The struct into which the titles will be stored.
 * @return The number of title lines read, or -1 if an error occurs.
 */
static int readTitles(FILE * psf, struct psf * p) {
  char buffer[1024]; //For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Read number of title lines from section header
  p->ntitle = readBang(psf, "!NTITLE");
  if(p->ntitle == -1)
    return -1; // Error reading number of title lines from section header

  // Allocate space for title array
  p->titles = malloc(p->ntitle * sizeof(char *));
  for (int i=0; i<p->ntitle; i++) {
    p->titles[i] = malloc(1024 * sizeof(char));
    // Initialize memory to null characters
    for (int j=0; j<1024; j++)
      p->titles[i][j]='\0';
  }

  int n = 0; // track the number of title lines read so far.

  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      break; // Error reading file for next bonds.
    }
    if(buffer[0]=='\n') {
      break; // Reached the end of the titles section.
    }
    if(strstr(buffer,"!")) {
      break; // PSF is missing empty line between sections.
    }
    strcpy(p->titles[n],buffer);
    if(buffer[strlen(buffer)-1]='\n')
      p->titles[n][strlen(buffer)-1]='\0';
    n++;
  } // This loop ends upon interruption by a break statement.
  if(n==p->ntitle)
    return n; // Correct number of title lines found.
  else { // Wrong number of title lines
    // Remove any partial title data from psf structure
    p->ntitle = -1;
    free(p->titles);
    p->titles=NULL;
    return -1;
  }
}

/**
 * Read atom information from a PSF.
 *
 * Reads atom details (including scattering lengths) from a PSF file and stores
 * the contents in the provided psf struct. Also populates the natom field of
 * the psf struct. The format should be a standard PSF, optionally with an
 * appended scattering length column (SLB), with the appropriate first-line
 * signature:
 *
 * PSF
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
 * PSF XPLOR
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)
 * PSF CMAP CHEQ
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
 * PSF CMAP CHEQ XPLOR
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)
 * PSF EXT
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
 * PSF EXT XPLOR
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8)
 * PSF EXT CMAP CHEQ
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
 * PSF EXT CMAP CHEQ XPLOR
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)
 * PSF SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,1X,G14.6)
 * PSF XPLOR SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,1X,G14.6)
 * PSF CMAP CHEQ SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * PSF CMAP CHEQ XPLOR SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * PSF EXT SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,1X,G14.6)
 * PSF EXT XPLOR SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,1X,G14.6)
 * PSF EXT CMAP CHEQ SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * PSF EXT CMAP CHEQ XPLOR SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6,1X,G14.6)
 *
 * The general format is:
 * II,LSEGID,LRESID,LRES,TYPE,IAC,CG,AMASS,IMOVE,[ECH,EHA],[SLD]
 * Formats sourced from charmm src file: source/io/psfres.src
 * The DRUDE extension to the PSF format is NOT supported, simply because the
 * impact of the DRUDE extension on a PSF is not understood by this author.
 *
 * The atom type is represented as a string if the XPLOR extension is used, 
 * otherwise it is represented as an integer. For simplicity, it is always 
 * stored in the psfatom struct as a string.
 *
 * @param[in] psf The PSF from which the atom data will be read.
 * @param[in,out] p The struct into which the atom data will be stored.
 * @return The number of atoms read, or -1 if an error occurs.
 */
static int readAtoms(FILE * psf, struct psf * p) {
  char buffer[1024]; // For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Read number of atoms from section header
  p->natom = readBang(psf, "!NATOM");
  if(p->natom == -1)
    return -1; // Error reading number of atoms from section header

  // Allocate space for atom array
  p->atoms = malloc(p->natom * sizeof(struct psfatom));

  int n = 0; // Track the number of atoms read so far

  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      break; // Error reading file for next atom record.
    }
    if(buffer[0]=='\n') {
      break; // Reached the end of the atoms section.
    }
    if(strstr(buffer,"!")) {
      break; // PSF is missing empty line between sections.
    }

    // Clear atom struct
    memset(p->atoms[n].seg,'\0',9);
    memset(p->atoms[n].resid,'\0',9);
    memset(p->atoms[n].res,'\0',9);
    memset(p->atoms[n].name,'\0',9);
    memset(p->atoms[n].type,'\0',9);
    p->atoms[n].charge=0;
    p->atoms[n].mass=0;
    p->atoms[n].imove=0;
    p->atoms[n].ech=0;
    p->atoms[n].eha=0;
    p->atoms[n].b=0;

    int offset=0; // Keeps track of field alignment
    char substr[15];

    int width=p->sig.ext?10:8; // Width of first field (index)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the atom index substring.
    if(atoi(substr)!=n+1)
      break; // Atom index doesn't match the expected number.
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=p->sig.ext?8:4; // Width of second field (segment name)
    // copy string directly into psfatom struct
    memcpy(p->atoms[n].seg,&buffer[offset],width);
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=p->sig.ext?8:4; // Width of third field (redisue identifier)
    // copy string directly into psfatom struct
    memcpy(p->atoms[n].resid,&buffer[offset],width);
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=p->sig.ext?8:4; // Width of fourth field (redisue name)
    // copy string directly into psfatom struct
    memcpy(p->atoms[n].res,&buffer[offset],width);
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=p->sig.ext?8:4; // Width of fifth field (atom name)
    // copy string directly into psfatom struct
    memcpy(p->atoms[n].name,&buffer[offset],width);
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=(p->sig.ext && p->sig.xplor)?6:4; // Width of sixth field (atom type)
    // copy string directly into psfatom struct
    memcpy(p->atoms[n].type,&buffer[offset],width);
    offset+=width; // advance offset to next field

    offset+=1; // Skip a space

    width=14; // Width of seventh field (atom charge)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the atom charge substring
    p->atoms[n].charge=atof(substr); // write charge value to psfatom struct
    offset+=width; // advance offset to next field

    width=14; // Width of eighth field (atom mass)
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the atom mass substring
    p->atoms[n].mass=atof(substr); // write mass value to psfatom struct
    offset+=width; // advance offset to next field

    width=8; // Width of ninth field ("IMOVE")
    memset(substr,'\0',15); // Clear substring buffer
    memcpy(substr,&buffer[offset],width); // isolate the "IMOVE" substring
    p->atoms[n].imove=atoi(substr); // write "IMOVE" to psfatom struct
    offset+=width; // advance offset to next field

    if(p->sig.cmapcheq) {
      width=14; // Width of first CMAP/CHEQ field (CHEQ electronegativity)
      memset(substr,'\0',15); // Clear substring buffer
      memcpy(substr,&buffer[offset],width); // isolate the "ECH" substring
      p->atoms[n].ech=atof(substr); // write charge value to psfatom struct
      offset+=width; // advance offset to next field

      width=14; // Width of second CMAP/CHEQ field (CHEQ hardness)
      memset(substr,'\0',15); // Clear substring buffer
      memcpy(substr,&buffer[offset],width); // isolate the "EHA" substring
      p->atoms[n].eha=atof(substr); // write mass value to psfatom struct
      offset+=width; // advance offset to next field
    }

    if(p->sig.slb) {
      offset+=1; // Skip a space

      width=14; // Width of last field (atom scattering length "b")
      memset(substr,'\0',15); // Clear substring buffer
      // isolate the scattering length substring
      memcpy(substr,&buffer[offset],width);
      p->atoms[n].b=atof(substr); // write scattering length value to struct
      offset+=width;
    }

    n++; // Advance atom count
  } // This loop ends upon interruption by a break statement.
  if(n==p->natom)
    return n; // Correct number of atoms found.
  else {
    // Remove any partial atom data from psf structure
    p->natom = -1;
    free(p->atoms);
    p->atoms=NULL;
    return -1;
  }
}

/**
 * Read bond data from a PSF
 *
 * Reads bond information from a PSF file and stores it in the provided psf
 * struct. Also populates the nbond field of the psf struct. Bonds are stored
 * as atom pairs, where atoms have been assigned to fields 'a' and 'b' of the
 * bond struct in the order they are listed in the PSF.
 *
 * @param[in] psf The PSF from which the bonds will be read.
 * @param[in,out] p The struct into which the bond data will be stored.
 * @return The number of bonds read, or -1 if an error occurs.
 */
static int readBonds(FILE * psf, struct psf * p) {
  char buffer[1024]; //For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Read number of bonds from section header
  p->nbond = readBang(psf, "!NBOND");
  if(p->nbond == -1)
    return -1; // Error reading number of bonds from section header

  // Allocate space for bond array
  p->bonds = malloc(p->nbond * sizeof(struct bond));

  int n = 0; // track the number of bonds read so far.

  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      break; // Error reading file for next bonds.
    }
    if(buffer[0]=='\n') {
      break; // Reached the end of the bonds section.
    }
    if(strstr(buffer,"!")) {
      break; // PSF is missing empty line between sections.
    }
    int numread = 0;
    numread = sscanf(buffer,"%d %d %d %d %d %d %d %d",&p->bonds[n].a,
      &p->bonds[n].b,&p->bonds[n+1].a,&p->bonds[n+1].b,&p->bonds[n+2].a,
      &p->bonds[n+2].b,&p->bonds[n+3].a,&p->bonds[n+3].b);
    for(int k=0; k<numread/2; k++) {
      p->bonds[n+k].a--;
      p->bonds[n+k].b--;
    }
    n+=numread/2;
  } // This loop ends upon interruption by a break statement.
  if(n==p->nbond)
    return n; // Correct number of bonds found.
  else { // Wrong number of bonds
    // Remove any partial bond data from psf structure
    p->nbond = -1;
    free(p->bonds);
    p->bonds=NULL;
    return -1;
  }
}

/**
 * Read angle data from a PSF
 *
 * Reads angle information from a PSF file and stores it in the provided psf
 * struct. Also populates the ntheta field of the psf struct. Angles are stored
 * as atom triplets, where atoms have been assigned to fields 'a', 'b', and 'c'
 * of the angle struct in the same order they are listed in the PSF.
 *
 * @param[in] psf The PSF from which the angles will be read.
 * @param[in,out] p The struct into which the angle data will be stored.
 * @return The number of angles read, or -1 if an error occurs.
 */
static int readAngles(FILE * psf, struct psf * p) {
  char buffer[1024]; //For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Read number of angles from section header
  p->ntheta = readBang(psf, "!NTHETA");
  if(p->ntheta == -1)
    return -1; // Error reading number of angles from section header

  // Allocate space for angle array
  p->angles = malloc(p->ntheta * sizeof(struct angle));

  int n = 0; // track the number of angles read so far.

  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      break; // Error reading file for next angles.
    }
    if(buffer[0]=='\n') {
      break; // Reached the end of the angles section.
    }
    if(strstr(buffer,"!")) {
      break; // PSF is missing empty line between sections.
    }
    int numread = 0;
    numread = sscanf(buffer,"%d %d %d %d %d %d %d %d %d",
        &p->angles[n].a,&p->angles[n].b,&p->angles[n].c,
        &p->angles[n+1].a,&p->angles[n+1].b,&p->angles[n+1].c,
        &p->angles[n+2].a,&p->angles[n+2].b,&p->angles[n+2].c);
    for(int k=0; k<numread/3; k++) {
      p->angles[n+k].a--;
      p->angles[n+k].b--;
      p->angles[n+k].c--;
    }
    n+=numread/3;
  } // This loop ends upon interruption by a break statement.
  if(n==p->ntheta)
    return n; // Correct number of angles found.
  else { // Wrong number of angles
    // Remove any partial angle data from psf structure
    p->ntheta = -1;
    free(p->angles);
    p->angles=NULL;
    return -1;
  }
}

/**
 * Read dihedral data from a PSF
 *
 * Reads dihedral information from a PSF file and stores it in the provided psf
 * struct. Also populates the nphi field of the psf struct. Dihedrals are stored
 * as atom quadruplets, where atoms have been assigned to fields 'a', 'b', 'c',
 * and 'd' of the dihedral struct in the same order that they are listed in the
 * PSF.
 *
 * @param[in] psf The PSF from which the dihedrals will be read.
 * @param[in,out] p The struct into which the dihedral data will be stored.
 * @return The number of dihedrals read, or -1 if an error occurs.
 */
static int readDihedrals(FILE * psf, struct psf * p) {
  char buffer[1024]; //For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Read number of dihedrals from section header
  p->nphi = readBang(psf, "!NPHI");
  if(p->nphi == -1)
    return -1; // Error reading number of dihedrals from section header

  // Allocate space for dihedral array
  p->dihedrals = malloc(p->nphi * sizeof(struct dihedral));

  int n = 0; // track the number of dihedrals read so far.

  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      break; // Error reading file for next dihedrals.
    }
    if(buffer[0]=='\n') {
      break; // Reached the end of the dihedrals section.
    }
    if(strstr(buffer,"!")) {
      break; // PSF is missing empty line between sections.
    }
    int numread = 0;
    numread = sscanf(buffer,"%d %d %d %d %d %d %d %d",
        &p->dihedrals[n].a,&p->dihedrals[n].b,
        &p->dihedrals[n].c,&p->dihedrals[n].d,
        &p->dihedrals[n+1].a,&p->dihedrals[n+1].b,
        &p->dihedrals[n+1].c,&p->dihedrals[n+1].d);
    for(int k=0; k<numread/4; k++) {
      p->dihedrals[n+k].a--;
      p->dihedrals[n+k].b--;
      p->dihedrals[n+k].c--;
      p->dihedrals[n+k].d--;
    }
    n+=numread/4;
  } // This loop ends upon interruption by a break statement.
  if(n==p->nphi)
    return n; // Correct number of dihedrals found.
  else { // Wrong number of dihedrals
    // Remove any partial dihedral data from psf structure
    p->nphi = -1;
    free(p->dihedrals);
    p->dihedrals=NULL;
    return -1;
  }
}

/**
 * Read improper dihedral data from a PSF
 *
 * Reads impropers from a PSF file and stores them in the provided psf struct.
 * Also populates the nimphi field of the psf struct. Dihedrals are stored as
 * atom quadruplets, where atoms have been assigned to fields 'a', 'b', 'c', and
 * 'd' of the dihedral struct in the same order that they are listed in the PSF.
 *
 * @param[in] psf The PSF from which the impropers will be read.
 * @param[in,out] p The struct into which the impropers will be stored.
 * @return The number of impropers read, or -1 if an error occurs.
 */
static int readImpropers(FILE * psf, struct psf * p) {
  char buffer[1024]; //For storing one line of the file at a time.
  char * ptr; // For detecting read failures

  // Read number of impropers from section header
  p->nimphi = readBang(psf, "!NIMPHI");
  if(p->nimphi == -1)
    return -1; // Error reading number of impropers from section header

  // Allocate space for improper dihedral array
  p->impropers = malloc(p->nimphi * sizeof(struct dihedral));

  int n = 0; // track the number of impropers read so far.

  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      break; // Error reading file for next impropers.
    }
    if(buffer[0]=='\n') {
      break; // Reached the end of the impropers section.
    }
    if(strstr(buffer,"!")) {
      break; // PSF is missing empty line between sections.
    }
    int numread = 0;
    numread = sscanf(buffer,"%d %d %d %d %d %d %d %d",
        &p->impropers[n].a,&p->impropers[n].b,
        &p->impropers[n].c,&p->impropers[n].d,
        &p->impropers[n+1].a,&p->impropers[n+1].b,
        &p->impropers[n+1].c,&p->impropers[n+1].d);
    for(int k=0; k<numread/4; k++) {
      p->impropers[n+k].a--;
      p->impropers[n+k].b--;
      p->impropers[n+k].c--;
      p->impropers[n+k].d--;
    }
    n+=numread/4;
  } // This loop ends upon interruption by a break statement.
  if(n==p->nimphi)
    return n; // Correct number of impropers found.
  else { // Wrong number of impropers
    // Remove any partial improper dihedral data from psf structure
    p->nimphi = -1;
    free(p->impropers);
    p->impropers=NULL;
    return -1;
  }
}

/**
 * Read data from PSF
 *
 * Parses a PSF file and creates a psf struct containing the successfully
 * parsed information, including flags from the PSF file signature and arrays of
 * title lines, atom records, bonds, angles, dihedrals, and improper dihedrals,
 * with associated counts for each array. Atom records are required, all other
 * sections are optional.
 *
 * In the event of a read error, the relevant section will be ignored, but all
 * other data will continue to be read and stored in the psf struct if possible.
 * If a read error occurs during an essential section (such as the PSF file
 * signature or the atom records), the entire PSF is invalidated by returning an
 * empty struct with the 'valid' subfield of the 'sig' field set to 'false'.
 * Before using the returned struct, one should verify that the 'valid'
 * subfield is set to 'true'.
 *
 * @param[in] path The filesystem path to the PSF to be parsed.
 * @return A psf struct containing all of the parsed information.
 */
struct psf readPSF(const char * path) {

  struct psf p = { .sig = { .valid = false,
                            .ext = false,
                            .cmapcheq = false,
                            .xplor = false,
                            .slb = false },
                   .ntitle = -1,
                   .natom = -1,
                   .nbond = -1,
                   .ntheta = -1,
                   .nphi = -1,
                   .nimphi = -1,
                   .atoms = NULL,
                   .bonds = NULL,
                   .angles = NULL,
                   .dihedrals = NULL,
                   .impropers = NULL};

  FILE * psf = fopen(path,"r");
  if(!psf) // Error encountered while opening file.
    return p;

  char buffer[1024]; //For storing one line of the file at a time.

  char * ptr; // For detecting read failures

  ptr = fgets(buffer,1024,psf); //get first line
  if(ferror(psf) || feof(psf) || ptr==NULL) {
    fclose(psf);
    return p; // Error reading first line of file.
  }

  if(strncmp(buffer,"PSF",3)) {
    fclose(psf);
    return p; // File is not a valid PSF.
  } else {
    p.sig.valid=true;
  }

  // Check for each PSF format specifier
  if(strstr(buffer,"EXT"))
    p.sig.ext=true;
  if(strstr(buffer,"CMAP CHEQ"))
    p.sig.cmapcheq=true;
  if(strstr(buffer,"XPLOR"))
    p.sig.xplor=true;
  if(strstr(buffer,"SLB"))
    p.sig.slb=true;

  readTitles(psf, &p);

  // Read atom data
  readAtoms(psf, &p);
  if(p.natom == -1) {
    fclose(psf);
    return p; // Error reading atom information.
  }

  // Read bond data
  readBonds(psf, &p);

  // Read angle data
  readAngles(psf, &p);

  // Read dihedral data
  readDihedrals(psf, &p);

  // Read improper dihedral data
  readImpropers(psf, &p);

  fclose(psf);
  return p;
}

/**
 * Frees the memory allocated for the psf struct
 *
 * Checks to make sure memory has not already been freed, and if not, frees the
 * arrays used to store title lines, atoms, bonds, angles, dihedrals, and
 * improper dihedrals. Once freed, it sets the pointers for each allocation to
 * NULL to prevent double freeing.
 *
 * @param[in] p the psf struct to be freed.
 */
void freePSF(struct psf p) {
  void sfree(void ** ptr) {
    if(*ptr) {
      free(*ptr);
      *ptr = NULL;
    }
  }
  for(int i=0; i<p.ntitle; i++)
    sfree((void **) &p.titles[i]);
  sfree((void **) &p.titles);
  sfree((void **) &p.atoms);
  sfree((void **) &p.bonds);
  sfree((void **) &p.angles);
  sfree((void **) &p.dihedrals);
  sfree((void **) &p.impropers);
}
