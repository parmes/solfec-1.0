include ../Config.mak
include ../Flags.mak

LIBBRICKS = ../obj/libsolfec.a
CFLAGS = -std=c99 $(DEBUG) $(PROFILE) $(OPENGL) -I..
LIB = -L../obj -lsolfec -lm $(BLAS) $(LAPACK) $(GLLIB) $(SPARSELIB)

TGT = glvtest\
      pbftest\
      hultest\
      cvitest\
      gjktest\
      boxtest\
      listest\
      comtest\

all: lib $(TGT)

lib:
	(cd .. && make)

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

comtest: comtest.c $(LIBBRICKS)
	$(MPICC) $(CFLAGS) $(PARALLEL) -o $@ $< $(LIB)
