
CC     = gcc
CFLAGS = -O2 -Wall -Wno-unused-variable -Wno-dangling-else

SRCS   = $(shell find . -type f -name '*.c')
OBJS   = $(patsubst %.c,%.o,$(SRCS))

all: libmathq.dylib mathq

libmathq.dylib: $(OBJS)
	$(CC) $(CFLAGS) -dynamiclib -o $@ $?

libmathq.so: $(OBJS)
	$(CC) $(CFLAGS) -shared -o $@ $?

mathq: MathQ.hs
	ghc $? libmathq.dylib

clean:
	rm -f MathQ.o MathQ.hi libmathq.dylib libmathq.so $(OBJS)
