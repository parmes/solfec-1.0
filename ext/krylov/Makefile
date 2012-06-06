#BHEADER**********************************************************************
# Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of HYPRE.  See file COPYRIGHT for details.
#
# HYPRE is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
#
# $Revision: 2.21 $
#EHEADER**********************************************************************

include ../../Config.mak

ifeq ($(DEBUG),yes)
  CFLAGS = -g -I.
else
  CFLAGS = -O2 -I.
endif

AR     = ar -rcu
RANLIB = ranlib

HEADERS =\
 krylov.h\
 bicgstab.h\
 cgnr.h\
 gmres.h\
 flexgmres.h\
 lgmres.h\
 pcg.h

FILES =\
 krylov.c\
 bicgstab.c\
 cgnr.c\
 gmres.c\
 flexgmres.c\
 lgmres.c\
 pcg.c

OBJS = ${FILES:.c=.o}

##################################################################
# Targets
##################################################################

all: ../../obj/libkrylov.a

clean:
	rm -f *.o lib*

##################################################################
# Rules
##################################################################

../../obj/libkrylov.a: ${OBJS}
	@echo  "Building $@ ... "
	${AR} $@ ${OBJS}
	${RANLIB} $@

${OBJS}: ${HEADERS}
