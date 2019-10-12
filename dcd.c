#include <stdlib.h>
#include <assert.h>

#include "dcd.h"

struct dcd {
  int is_pdb;
  FILE *hdl;
  uint32_t  nframes;
  uint32_t  natoms;
  long int  offset;
};

/**
 * Opens a DCD file and returns a handle.
 *
 * @param[in] path The path to the DCD file.
 * @return A handle to the DCD file.
 */
struct dcd *openDCD(char *path) {
  struct dcd *d = malloc(sizeof(struct dcd));

  d->hdl=fopen(path,"r");

  if( ! d->hdl )
    return NULL;

  fseek(d->hdl, 8, SEEK_SET);
  assert(1==fread(&(d->nframes), 4, 1, d->hdl));
  fseek(d->hdl, 96, SEEK_SET);
  uint32_t ts;
  assert(1==fread(&ts, 4, 1, d->hdl));
  fseek(d->hdl, 80*((long int) ts) + 8, SEEK_CUR);
  assert(1==fread(&(d->natoms), 4, 1, d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  d->offset = ftell(d->hdl);

  return d;
}

struct dcd *openWritableDCD(char *path) {
  struct dcd *d = malloc(sizeof(struct dcd));

  d->hdl=fopen(path,"r+");

  if( ! d->hdl )
    return NULL;

  fseek(d->hdl, 8, SEEK_SET);
  assert(1==fread(&(d->nframes), 4, 1, d->hdl));
  fseek(d->hdl, 96, SEEK_SET);
  uint32_t ts;
  assert(1==fread(&ts, 4, 1, d->hdl));
  fseek(d->hdl, 80*((long int) ts) + 8, SEEK_CUR);
  assert(1==fread(&(d->natoms), 4, 1, d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  d->offset = ftell(d->hdl);

  return d;
}

/**
 * Closes the DCD file descriptor and frees the associated memory.
 *
 * @param[in] d The dcd handle.
 */
void closeDCD(struct dcd *d) {
  fclose(d->hdl);
  free(d);
}

/**
 * Gets the number of frames in the DCD.
 *
 * @param[in] d The dcd handle.
 * @return The number of frames in the DCD.
 */
uint32_t getNFrames(struct dcd *d) {
  return d->nframes;
}

/**
 * Prepares the DCD handle to read the desired frame.
 *
 * Moves the position indicator within the DCD file descriptor. The specified
 * frame data can then be read using getUnitCell or getCoords.
 *
 * @param[in] d The DCD file handle
 * @param[in] f The zero-indexed frame number.
 */
void goToFrame(struct dcd *d,uint32_t f) {
  fseek(d->hdl,((long int) (d->offset)) + ((long int) f)*(12*((long int) (d->natoms)) + 80),SEEK_SET);
}

/**
 * Prepares the DCD handle to read the next frame.
 *
 * Moves the position indicator within the DCD file descriptor. The next frame
 * data can then be read using getUnitCell or getCoords.
 *
 * @param[in] d The DCD file handle
 */
void nextFrame(struct dcd *d) {
  fseek(d->hdl,(12*((long int) (d->natoms)) + 80),SEEK_CUR);
}

/**
 * Reads the current position of the DCD handle.
 *
 * Gets the number of the current frame (the frame which would be read by a
 * subsequent call to getUnitCell or getCoords).
 *
 * @param[in] d The DCD file handle
 * @return The number of the current frame
 */
uint32_t getFrame(struct dcd *d) {
  return (ftell(d->hdl) - ((long int) (d->offset)))/(12*((long int) (d->natoms)) + 80);
}

/**
 * Reads the unit cell information for the current frame.
 *
 * Gets the unit cell data from the current frame and stores it in the provided
 * array.
 *
 * @param[in] d The DCD file handle
 * @param[out] uc The array into which the unit cell data should be placed.
 */
void getUnitCell(struct dcd *d, double *uc) {
  fseek(d->hdl, 4, SEEK_CUR);
  assert(1==fread(&(uc[0]), 8, 1, d->hdl));
  fseek(d->hdl, 8, SEEK_CUR);
  assert(1==fread(&(uc[1]), 8, 1, d->hdl));
  fseek(d->hdl, 16, SEEK_CUR);
  assert(1==fread(&(uc[2]), 8, 1, d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, -56, SEEK_CUR);
}

/**
 * Reads the coordinate information for the current frame.
 *
 * Gets the coordinate data from the current frame and stores it in the provided
 * arrays.
 *
 * @param[in] d The DCD file handle
 * @param[out] xs The array into which the x-coordinates should be stored.
 * @param[out] ys The array into which the y-coordinates should be stored.
 * @param[out] zs The array into which the z-coordinates should be stored.
 */
void getCoords(struct dcd *d, float *xs, float *ys, float *zs) {
  fseek(d->hdl, 56, SEEK_CUR);
  fseek(d->hdl, 4, SEEK_CUR);
  assert((d->natoms)==fread(xs, 4, (d->natoms), d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, 4, SEEK_CUR);
  assert((d->natoms)==fread(ys, 4, (d->natoms), d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, 4, SEEK_CUR);
  assert((d->natoms)==fread(zs, 4, (d->natoms), d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, (-1)*(12*((long int) (d->natoms)) + 80), SEEK_CUR);
}

void writeCoords(struct dcd *d, float *xs, float *ys, float *zs) {
  fseek(d->hdl, 56, SEEK_CUR);
  fseek(d->hdl, 4, SEEK_CUR);
  assert((d->natoms)==fwrite(xs, 4, (d->natoms), d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, 4, SEEK_CUR);
  assert((d->natoms)==fwrite(ys, 4, (d->natoms), d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, 4, SEEK_CUR);
  assert((d->natoms)==fwrite(zs, 4, (d->natoms), d->hdl));
  fseek(d->hdl, 4, SEEK_CUR);
  fseek(d->hdl, (-1)*(12*((long int) (d->natoms)) + 80), SEEK_CUR);
}
