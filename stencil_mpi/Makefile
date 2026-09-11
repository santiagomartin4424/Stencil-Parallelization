CC      = mpicc
CFLAGS  = -Wall -g -O4 -std=c99
LDLIBS  = -lm -lrt

all: stencil

stencil: stencil.c
	$(CC) $(CFLAGS) stencil.c -o stencil $(LDLIBS)

clean:
	-rm -f stencil

mrproper: clean
	-rm -f *~

