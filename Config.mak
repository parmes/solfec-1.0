#
# Specify C compiler here
#

CC = cc

# 
# Debug or optimized version switch (yes/no)
#

DEBUG = yes
PROFILE = no
MEMDEBUG = no
GEOMDEBUG = no
PARDEBUG = no
NOTHROW = yes

#
# TIMERS (enable/disable detailed solver timings)
#

TIMERS = yes

#
# POSIX
#

POSIX = yes

#
# XDR
#

XDR = no

XDRINC = 
XDRLIB = 

#
# BLAS
#

BLAS = -L/usr/lib -lblas

#
# LAPACK
#

LAPACK = -L/usr/lib -llapack

#
# Python
#

PYTHON = -I/Library/Frameworks/Python.framework/Versions/Current/include/python2.7
PYTHONLIB = -framework Python

#
# OpenGL (yes/no)
#

OPENGL = yes

GLINC =
GLLIB = -framework GLUT -framework OpenGL

#
# VBO  (OPENGL == yes)
#

VBO = yes

#
# MPI (yes/no)
#

MPI = yes
MPICC = mpicc

#
# Zoltan (MPI == yes)
#

ZOLTANINC = -I/Users/tomek/Devel/lib/zoltan/include
ZOLTANLIB = -L/Users/tomek/Devel/lib/zoltan/lib -lzoltan
