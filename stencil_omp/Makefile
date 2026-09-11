CC      = mpicc
CFLAGS += -Wall -g -O3 -std=c99 -fopenmp
LDLIBS += -lm -lrt

TARGET  = stencil_omp
SRC     = stencil_omp.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDLIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
