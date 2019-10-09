CFLAGS = -std=c99
LDFLAGS = -lm

all: test

test: test.c psf.o dcd.o betaz.o

.PHONY: clean
clean:
	-rm -f test psf.o dcd.o betaz.o
