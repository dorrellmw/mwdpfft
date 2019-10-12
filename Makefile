CFLAGS = -std=c99 -fopenmp
LDLIBS = -lm

all: test tpfft

test: test.c psf.o dcd.o betaz.o

tpfft: tpfft.c psf.o dcd.o betaz.o readArray.o pfft.o

.PHONY: clean
clean:
	-rm -f test tpfft *.o
