include ../Config.mak
include ../Flags.mak

LIBBRICKS = ../obj/libsolfec.a
LIBBRICKSMPI = ../obj/libsolfec-mpi.a
CFLAGS = -std=c99 $(DEBUG) $(PROFILE) $(OPENGL) -I..
LIB = -L../obj -lsolfec -lm $(BLAS) $(LAPACK) $(GLLIB)
ifeq ($(MPI),yes)
  LIBMPI = -L../obj -lsolfec-mpi -lm $(BLAS) $(LAPACK) $(MPILIBS)
endif

TGT = glvtest\
      pbftest\
      hultest\
      cvitest\
      gjktest\
      boxtest\
      listest\

ifeq ($(MPI),yes)

TGTMPI = comtest\

else
TGTMPI = 
endif

all: lib $(TGT) $(TGTMPI)

lib:
	(cd .. && make all)

del:
	rm -f *.dat
	rm -f *.plt
	rm -f *.ps
	rm -f *.eps
	rm -f *.bmp
	rm -f *.tif
	rm -f *.out
	rm -f *.idx
	rm -f *.lab

clean:
	rm -f core *.o
	rm -f *out
	rm -f $(TGT)
	rm -fr *dSYM

../obj/libdebris.a:
	(cd .. && make)

glvtest: glvtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

pbftest: pbftest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

hultest: hultest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

cvitest: cvitest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

gjktest: gjktest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

boxtest: boxtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

listest: listest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

# MPI

comtest: comtest.c $(LIBBRICKSMPI)
	$(MPICC) $(CFLAGS) $(MPITGT) -o $@ $< $(LIBMPI)
