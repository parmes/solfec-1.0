include ../../../Config.mak
include ../Makefile.inc
CPPFLAGS = -I../utilities 

all: multivector.o multiserial.o pcg_multi.o

multivector.o: multivector.c multivector.h interpreter.h
multiserial.o: multiserial.c multiserial.h
pcg_multi.o: pcg_multi.c pcg_multi.h


clean: FORCE
	rm -f multivector.o multiserial.o

FORCE:
