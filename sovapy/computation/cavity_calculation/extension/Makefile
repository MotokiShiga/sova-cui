CFLAGS = -O3
SRCS := algorithm.c

all: libalgorithm.so

libalgorithm.so:
	$(CC) -Wall -Wextra -fPIC $(CFLAGS) -o $@ $(SRCS) -shared

clean:
	rm -f *.pyc *.pyo
	rm -f *.o *.so

