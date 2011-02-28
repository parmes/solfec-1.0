###############################################################################
#                                                                             #
# makefile for TetGen                                                         #
#                                                                             #
# Type "make" to compile TetGen into an executable program (tetgen).          #
# Type "make tetlib" to compile TetGen into a library (libtet.a).             #
# Type "make distclean" to delete all object (*.o) files.                     #
#                                                                             #
###############################################################################

# CXX should be set to the name of your favorite C++ compiler.
# ===========================================================

include ../../Config.mak

# CXXFLAGS is the level of optimiztion, default is -O. One should try
# -O2, -O3 ... to find the best optimization level.
# ===================================================================
#
# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. It should
# always be equal to -O0 (no optimization). Otherwise, TetGen may not
# work properly.

ifeq ($(DEBUG),yes)
  CXXFLAGS = -g
  PREDCXXFLAGS = -g
else
  CXXFLAGS = -O
  PREDCXXFLAGS = -O0
endif

# SWITCHES is a list of switches to compile TetGen.
# =================================================
#
# By default, TetGen uses double precision floating point numbers.  If you
#   prefer single precision, use the -DSINGLE switch. 
#
# The source code of TetGen includes a lot of assertions, which are mainly
#   used for catching bugs at that places.  These assertions somewhat slow
#   down the speed of TetGen.  They can be skipped by define the -DNDEBUG
#   switch.

SWITCHES = -w -DSELF_CHECK

# SWITCHES = -Wall -Wabi -Wctor-dtor-privacy \
#            -Woverloaded-virtual -Wno-pmf-conversions -Wsign-promo \
#            -Wsynth  -Wchar-subscripts -Wconversion -Wsign-compare \
#            -Wcomment  -Wimplicit -Wmissing-braces -Wparentheses \
#            -Wreturn-type -Wswitch -Wswitch-default \
#            -Wswitch-enum -Wtrigraphs -W -DSELF_CHECK

# RM should be set to the name of your favorite rm (file deletion program).

RM = /bin/rm -f

# The action starts here.

tetlib: tetgen.o predicates.o tetsol.o
	ar -rcu ../../obj/libtet.a tetgen.o predicates.o tetsol.o

tetgen:	tetgen.cxx predicates.o
	$(CXX) $(CXXFLAGS) $(SWITCHES) -o tetgen tetgen.cxx predicates.o -lm

tetgen.o: tetgen.cxx tetgen.h
	$(CXX) $(CXXFLAGS) $(SWITCHES) -DTETLIBRARY -c tetgen.cxx

tetsol.o: tetsol.cxx tetsol.h
	$(CXX)  $(CXXFLAGS) $(SWITCHES) -DTETLIBRARY  -c tetsol.cxx

predicates.o: predicates.cxx
	$(CXX) $(PREDCXXFLAGS) -c predicates.cxx

clean:
	$(RM) *.o *.a tetgen *~
