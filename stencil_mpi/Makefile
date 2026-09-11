CC      = mpicc
CFLAGS += -Wall -g -O3 -std=c99
LDLIBS += -lm -lrt

TARGET  = stencil_mpi
SRC     = stencil_mpi.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDLIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
