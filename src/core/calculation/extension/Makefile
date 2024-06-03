CFLAGS = -O3
INCLUDES := $(shell ./build-flags --include-flags gr)
INCLUDES1 := -I/Users/morita/Data/github/gr/lib/gr
INCLUDES2 := -I/Users/morita/Data/github/gr/lib/gr3
ifeq ($(PYLIBPATH),)
	LDFLAGS := $(shell ./build-flags --lib-flags gr) $(EXTRA_LD_FLAGS)
else
	LDFLAGS := -L$(PYLIBPATH)/python2.7/site-packages/gr3/ -lGR3 $(EXTRA_LD_FLAGS)
endif
LDFLAGS1 := -L/Users/morita/Data/github/gr/build
SUBDIRS=


all: libalgorithm.so subdirs

libalgorithm.so: algorithm.c
	$(CC) -Wall -Wextra -fPIC $(CFLAGS) $(INCLUDES) $(INCLUDES1) $(INCLUDES2) -o $@ $< -shared $(LDFLAGS1)

subdirs:
	for i in $(SUBDIRS); do make -C $$i all; done

clean:
	rm -f *.pyc *.pyo
	rm -f *.o *.so
	for i in $(SUBDIRS); do make -C $$i clean; done

