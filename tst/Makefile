include ../Config.mak
include ../Flags.mak

LIBBRICKS = ../obj/libsolfec.a
LIBBRICKSMPI = ../obj/libsolfec-mpi.a
CFLAGS = -std=c99 $(DEBUG) $(PROFILE) $(OPENGL) -I..
LIB = -L../obj -lsolfec -lm $(BLAS) $(LAPACK) $(GLLIB) $(PYTHONLIB)
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
      maptest\
      mtxtest\
      cmptest\

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
	rm -f $(TGT) $(TGTMPI)
	rm -fr *dSYM
	rm -f obj/*

glvtest: obj/glvtest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/glvtest.o: glvtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

pbftest: obj/pbftest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/pbftest.o: pbftest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

hultest: obj/hultest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/hultest.o: hultest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

cvitest: obj/cvitest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/cvitest.o: cvitest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

gjktest: obj/gjktest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/gjktest.o: gjktest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

boxtest: obj/boxtest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/boxtest.o: boxtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

listest: obj/listest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/listest.o: listest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

maptest: obj/maptest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/maptest.o: maptest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

mtxtest: obj/mtxtest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/mtxtest.o: mtxtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

cmptest: obj/cmptest.o $(LIBBRICKS)
	$(CC) $(CFLAGS) -o $@ $< $(LIB)

obj/cmptest.o: cmptest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

# MPI

comtest: obj/comtest.o $(LIBBRICKSMPI)
	$(MPICC) $(CFLAGS) $(MPITGT) -o $@ $< $(LIBMPI)

obj/comtest.o: comtest.c $(LIBBRICKSMPI)
	$(MPICC) $(CFLAGS) $(MPITGT) -c -o $@ $<
