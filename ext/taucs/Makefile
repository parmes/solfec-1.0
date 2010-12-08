#**********************************************************************
# TAUCS Makefile
#**********************************************************************

include ../../Config.mak

ifeq ($(DEBUG),yes)
  CFLAGS = -g -I.
else
  CFLAGS = -O2 -I.
endif

AR     = ar -rcu
RANLIB = ranlib

HEADERS =\
 taucs.h\
 taucs_private.h

FILESG =\
 blas_aux.c\
 taucs_ccs_base.c\
 taucs_ccs_factor_llt.c\
 taucs_ccs_io.c\
 taucs_ccs_ooc_llt.c\
 taucs_ccs_ooc_lu.c\
 taucs_ccs_ops.c\
 taucs_ccs_order.c\
 taucs_ccs_solve_llt.c\
 taucs_linsolve.c\
 taucs_logging.c\
 taucs_malloc.c\
 taucs_memory.c\
 taucs_ooc_io.c\
 taucs_recvaidya.c\
 taucs_sn_llt.c\
 taucs_timer.c\
 taucs_vec_base.c

FILESD =\
 taucs_ccs_base.c\
 taucs_ccs_factor_llt.c\
 taucs_ccs_io.c\
 taucs_ccs_ooc_llt.c\
 taucs_ccs_ooc_lu.c\
 taucs_ccs_ops.c\
 taucs_ccs_solve_llt.c\
 taucs_ccs_xxt.c\
 taucs_complex.c\
 taucs_gremban.c\
 taucs_iter.c\
 taucs_logging.c\
 taucs_malloc.c\
 taucs_memory.c\
 taucs_recvaidya.c\
 taucs_sn_llt.c\
 taucs_timer.c\
 taucs_vaidya.c\
 taucs_vec_base.c

OBJSG = ${FILESG:.c=.g}

OBJSD = ${FILESD:.c=.d}

##################################################################
# Targets
##################################################################

all: ../../obj/libtaucs.a

clean:
	rm -f *.g *.d

##################################################################
# Rules
##################################################################

.SUFFIXES: .g .d .c

.c.g:
	$(CC) -DTAUCS_CORE_GENERAL $(CFLAGS) -c -o $@ $<

.c.d:
	$(CC) -DTAUCS_CORE_DOUBLE $(CFLAGS) -c -o $@ $<

../../obj/libtaucs.a: ${OBJSG} ${OBJSD}
	@echo  "Building $@ ... "
	${AR} $@ ${OBJSG} ${OBJSD}
	${RANLIB} $@

${OBJSG}: ${HEADERS}

${OBJSD}: ${HEADERS}
