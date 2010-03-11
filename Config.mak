#
# Specify C compiler here
#

CC = cc

#
# Specify FORTRAN compiler here
#

FF = g95

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
# POSIX
#

POSIX = yes

#
# BLAS
#

BLAS = -L/usr/lib -lblas

#
# LAPACK
#

LAPACK = -L/usr/lib -llapack

#
# SuperLU
#

SUPERLUINC = -I/Users/tomek/Devel/lib/superlu/include
SUPERLULIB = -L/Users/tomek/Devel/lib/superlu/lib -lsuperlu

#
# Python
#

PYTHON = -I/usr/include/python2.5
PYTHONLIB = -L/usr/lib -lpython2.5

#
# OpenGL (yes/no)
#

OPENGL = yes

GLINC =
GLLIB = -framework GLUT -framework OpenGL

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

#
# HYPRE (MPI == yes)
#

HYPREINC = -I/Users/tomek/Devel/lib/hypre/include
HYPRELIB = -L/Users/tomek/Devel/lib/hypre/lib -lHYPRE
