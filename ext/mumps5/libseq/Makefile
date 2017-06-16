#
#  This file is part of MUMPS 5.1.1, released
#  on Mon Mar 20 14:34:33 UTC 2017
#
libdir = ../../../obj

all: libmpiseq

.PHONY: all libmpiseq clean

include ../../../Config.mak
include ../Makefile.inc

libmpiseq: $(libdir)/libmpiseq$(PLAT).a

$(libdir)/libmpiseq$(PLAT).a: mpi.o mpic.o elapse.o
	$(AR) $@ mpi.o mpic.o elapse.o
	$(RANLIB) $@
.f.o:
	$(FC) $(OPTF) -w -c $*.f
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -c $*.c

clean:
	$(RM) *.o *.a
