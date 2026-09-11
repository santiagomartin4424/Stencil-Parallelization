CC      = gcc
CFLAGS += -Wall -g -O4 -std=c99
LDLIBS += -lm -lrt

all: stencil

clean:
	-rm stencil

mrproper: clean
	-rm *~

archive: stencil.c Makefile
	( cd .. ; tar czf stencil.tar.gz stencil/ )
