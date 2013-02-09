#
#  This file is part of MUMPS 4.9.2, built on Thu Nov  5 07:05:08 UTC 2009
#
all: libmpiseq

.PHONY: all libmpiseq clean

include ../../../Config.mak
include ../Makefile.inc

libmpiseq: libmpiseq$(PLAT).a

libmpiseq$(PLAT).a: mpi.o elapse.o
	$(AR) $@ mpi.o elapse.o
	$(RANLIB) $@
.f.o:
	$(FC) $(OPTF) -w -c $*.f
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -w -c $*.c

clean:
	$(RM) *.o *.a
