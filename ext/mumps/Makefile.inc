#
#  This file is part of MUMPS 4.9.2, built on Thu Nov  5 07:05:08 UTC 2009
#
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. If you would like to
# use other orderings, you need to obtain the corresponding package and modify
# the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 5.1 and later) orderings are now available for MUMPS.
#

#SCOTCHDIR  = ${HOME}/scotch_5.1_esmumps
#ISCOTCH    = -I$(SCOTCHDIR)/include 
# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr


#LPORDDIR = $(topdir)/PORD/lib/
#IPORD    = -I$(topdir)/PORD/include/
#LPORD    = -L$(LPORDDIR) -lpord

LMETISDIR = ../../obj
IMETIS    = # Metis doesn't need include files (Fortran interface avail.)

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

LMETIS    = -L$(LMETISDIR) -lmetis
#LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
ORDERINGSF  = -Dmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################
################################################################################

RM = /bin/rm -f
FL = $(FC)
AR = ar vr
RANLIB  = echo

INCSEQ = -Ilibseq
LIBSEQ  =  -Llibseq -lmpiseq

LIBBLAS = $(BLAS)
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
ifeq ($(OS),WIN32)
  OSDEF = -DWITHOUT_PTHREAD
else
  OSDEFS =
endif

#Begin Optimization options
ifeq ($(FC),g95)
  CDEFS   = -DAdd__ $(OSDEFS)
  OPTF   = -O -i4
  OPTL   = -O  
  OPTC   = -O -DMAIN_COMP 
endif
ifeq ($(FC),gfortran)
  CDEFS   = -DAdd_ $(OSDEFS)
  OPTF   = -O -Dintel_ -DALLOW_NON_INIT
  OPTL   = -O  
  OPTC   = -O
endif
ifeq ($(FC),ifc)
  CDEFS   = -DAdd_ $(OSDEFS)
  OPTF    = -O -Dintel_ -DALLOW_NON_INIT -nofor_main
  OPTL    = -O -nofor_main
  OPTC    = -O
endif
ifeq ($(FC),ifort)
  CDEFS   = -DAdd_ $(OSDEFS)
  OPTF    = -O -Dintel_ -DALLOW_NON_INIT -nofor_main
  OPTL    = -O -nofor_main
  OPTC    = -O
endif
#End Optimization options

INC = $(INCSEQ)
LIB = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
