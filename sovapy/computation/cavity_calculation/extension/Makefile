CFLAGS = -O3
INCLUDES1 := -I../../../../external_lib/gr/lib/gr
INCLUDES2 := -I../../../../external_lib/gr/lib/gr3
LDFLAGS := -L../../../../external_lib/gr/build
SRCS := algorithm.c gr3_mc.c

all: libalgorithm.so

libalgorithm.so:
	$(CC) -Wall -Wextra -fPIC $(CFLAGS) $(INCLUDES1) $(INCLUDES2) -o $@ $(SRCS) -shared $(LDFLAGS)

clean:
	rm -f *.pyc *.pyo
	rm -f *.o *.so

