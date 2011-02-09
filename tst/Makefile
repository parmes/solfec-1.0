include ../Config.mak
include ../Flags.mak

LIBBRICKS = ../obj/libsolfec.a
LIBBRICKSMPI = ../obj/libsolfec-mpi.a
CFLAGS = $(STD) $(DEBUG) $(PROFILE) $(OPENGL) $(XDRINC) -I..
LIB = -L../obj -lsolfec -lkrylov -lmetis -ltaucs -lm -lstdc++ $(LAPACK) $(BLAS) $(GLLIB) $(PYTHONLIB) $(XDRLIB) $(CUDALIB)
ifeq ($(MPI),yes)
  LIBMPI = -L../obj -lsolfec-mpi -lkrylov -lmetis -ltaucs -lm -lstdc++ $(LAPACK) $(BLAS) $(MPILIBS) $(XDRLIB) $(CUDALIB)
endif

TGT = glvtest\
      anytest\
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
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/glvtest.o: glvtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

anytest: obj/anytest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/anytest.o: anytest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

pbftest: obj/pbftest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/pbftest.o: pbftest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

hultest: obj/hultest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/hultest.o: hultest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

cvitest: obj/cvitest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/cvitest.o: cvitest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

gjktest: obj/gjktest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/gjktest.o: gjktest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

boxtest: obj/boxtest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/boxtest.o: boxtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

listest: obj/listest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/listest.o: listest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

maptest: obj/maptest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/maptest.o: maptest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

mtxtest: obj/mtxtest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/mtxtest.o: mtxtest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

cmptest: obj/cmptest.o $(LIBBRICKS)
	$(LL) $(CFLAGS) -o $@ $< $(LIB)

obj/cmptest.o: cmptest.c $(LIBBRICKS)
	$(CC) $(CFLAGS) -c -o $@ $<

# MPI

comtest: obj/comtest.o $(LIBBRICKSMPI)
	$(MPILL) $(CFLAGS) $(MPITGT) -o $@ $< $(LIBMPI)

obj/comtest.o: comtest.c $(LIBBRICKSMPI)
	$(MPICC) $(CFLAGS) $(MPITGT) -c -o $@ $<
