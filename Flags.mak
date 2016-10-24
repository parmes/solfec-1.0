#
# Compilation flags setup
#

ifeq ($(OS),WIN32)
  OS = -DOSTYPE_WIN32
endif

ifeq ($(OS),SOLARIS)
  OS = -DOSTYPE_SOLARIS
endif

ifeq ($(OS),LINUX)
  OS = -DOSTYPE_LINUX
endif

ifeq ($(OS),AIX)
  OS = -DOSTYPE_AIX
endif

ifeq ($(OS),IRIX)
  OS = -DOSTYPE_IRIX
endif

ifeq ($(OS),OSX)
  OS = -DOSTYPE_OSX
endif

ifeq ($(OS),FREEBSD)
  OS = -DOSTYPE_FREEBSD
endif

ifeq ($(POSIX),yes)
  STD = -std=c99 -DPOSIX
else
  STD = -std=c99
endif

ifneq ($(HDF5),yes)
  HDF5INC =
  HDF5LIB = 
  HDF5 = 
else
  HDF5 = -DHDF5 $(HDF5INC)
endif

ifneq ($(ZOLTAN),yes)
  LBINC = -I$(DYNLB)
  LBLIB = -L$(DYNLB) -ldynlb8
  LB = -DDYNLB $(LBINC)
else
  LBINC = $(ZOLTANINC)
  LBLIB = $(ZOLTANLIB)
  LB = -DZOLTAN $(LBINC)
endif

ifneq ($(XDR),yes)
  XDRINC =
  XDRLIB = 
endif

ifeq ($(TIMERS),yes)
  TIMERS = -DTIMERS
else
  TIMERS = 
endif

ifeq ($(OPENGL),yes)
  ifeq ($(VBO),yes)
    OPENGL = -DOPENGL -DVBO $(GLINC)
  else
    OPENGL = -DOPENGL $(GLINC)
  endif
else
  OPENGL =
  GLLIB = 
endif

ifeq ($(LOCAL_BODIES),yes)
  LOCAL_BODIES = -DLOCAL_BODIES
else
  LOCAL_BODIES = 
endif

ifeq ($(DEBUG),yes)
  DBG = yes
  DEBUG =  -W -Wall -Wno-unused-parameter -pedantic -g -DDEBUG
  ifeq ($(PROFILE),yes)
    PROFILE = -p
  else
    PROFILE =
  endif
  ifeq ($(NOTHROW),yes)
    NOTHROW = -DNOTHROW
  else
    NOTHROW =
  endif
  ifeq ($(MEMDEBUG),yes)
    MEMDEBUG = -DMEMDEBUG
  else
    MEMDEBUG =
  endif
  ifeq ($(GEOMDEBUG),yes)
    GEOMDEBUG = -DGEOMDEBUG
  else
    GEOMDEBUG =
  endif
else
  DBG = no
  DEBUG =  -w -pedantic -O3 -funroll-loops
  PROFILE =
  MEMDEBUG =
  GEOMDEBUG =
  NOTHROW =
endif

ifeq ($(MPI),yes)
  ifeq ($(PARDEBUG),yes)
    PARDEBUG = -DPARDEBUG
  else
    PARDEBUG =
  endif
  ifeq ($(PSCTEST),yes)
    PSCTEST = -DPSCTEST
  else
    PSCTEST =
  endif

  MPIFLG = -DMPI $(LB) $(PARDEBUG) $(PSCTEST)
  MPILIBS = $(LBLIB)
endif

ifeq ($(SICONOS),yes)
  SICONOS = -DWITHSICONOS
else
  SICONOS = 
  SICONOSINC =
  SICONOSLIB =
endif

ifdef PARMEC
  PARMECINC = -DPARMEC -I$(PARMEC)
  PARMECLIB = -L$(PARMEC)
else
  PARMECLIB =
endif
