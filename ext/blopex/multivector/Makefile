include ../../../Config.mak
include ../Makefile.inc
CPPFLAGS = -I../utilities 

all: multivector.o multiserial.o

multivector.o: multivector.c multivector.h interpreter.h
multiserial.o: multiserial.c multiserial.h


clean: FORCE
	rm -f multivector.o multiserial.o

FORCE:
