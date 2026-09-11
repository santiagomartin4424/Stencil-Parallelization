CC      = mpicc
# Añadimos -fopenmp aquí
CFLAGS  = -Wall -g -O3 -std=c99 -fopenmp
LDLIBS  = -lm -lrt -fopenmp

all: stencil

stencil: stencil.c
	$(CC) $(CFLAGS) stencil.c -o stencil $(LDLIBS)

clean:
	-rm -f stencil

