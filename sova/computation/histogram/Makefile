CFLAGS = -O3 
UNAME := $(shell uname -s)
INCLUDES := $(shell ./build-flags --include-flags python numpy)
ifeq ($(PYLIBPATH),)
	LIBS := $(shell ./build-flags --lib-flags python numpy)
else
	LIBS = -L$(PYLIBPATH) -lpython2.7
endif


all: histogram.so
	

histogram.so: histogram.o
	$(CC) -shared -o $@ $^ $(LIBS)

histogram.o: histogram.c
	$(CC) -fPIC -c $(CFLAGS) $(INCLUDES) -o $@ $<

clean:
	rm -rf *.o *.so

