include ../../Config.mak
include Makefile.inc

DIRS = multivector krylov utilities

HEADERS = \
multivector/interpreter.h \
multivector/multivector.h \
krylov/lobpcg.h \
multivector/temp_multivector.h \
utilities/fortran_matrix.h \
utilities/fortran_interpreter.h \
utilities/fortran_options.h

LIB = libBLOPEX.${AR_LIB_SUFFIX}

##################################################################
# Targets
##################################################################

all: $(LIB) 

objects: 
	@ \
        mkdir -p lib; \
        for i in ${DIRS}; \
        do \
          echo "Making $$i ..."; \
          (cd $$i && $(MAKE)); \
          (cp -fp $$i/***.o ./lib/); \
          echo ""; \
        done

clean: 
	@ \
	rm -fr lib; \
        for i in ${DIRS}; \
        do \
          if [ -d $$i ]; \
          then \
            echo "Cleaning $$i ..."; \
            (cd $$i &&  $(MAKE) $@); \
          fi; \
        done
	

##################################################################
# Rules
##################################################################

$(LIB) : objects 
	@echo  "Building $@ ... "
	cd lib && ${AR} $@ ./*.o 
	cd lib && ${RANLIB} $@
	mv lib/$@ ../../obj/
