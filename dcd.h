#ifndef DCD_H_
#define DCD_H_

#include <stdio.h>
#include <stdint.h>

struct dcd;

struct dcd *openDCD(char *);
struct dcd *openWritableDCD(char *);
void closeDCD(struct dcd *);
uint32_t getNFrames(struct dcd *);
void goToFrame(struct dcd *,uint32_t);
void nextFrame(struct dcd *);
uint32_t getFrame(struct dcd *);
void getUnitCell(struct dcd *, double *);
void getCoords(struct dcd *, float *, float *, float *);
void writeCoords(struct dcd *, float *, float *, float *);

#endif
