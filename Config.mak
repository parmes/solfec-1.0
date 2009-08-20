#
# Specify C compiller here
#

CC = cc

# 
# Debug or optimized version switch (yes/no)
#

DEBUG = yes
PROFILE = no
MEMDEBUG = no
GEOMDEBUG = no
NOTHROW = no

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
MPITHREADS = no # use MPI threads
OWNASYNC = on # implement own asynchronous communication using pthreads and blocking communication (used when MPITHREADS == yes)
MPILIB = -L/opt/openmpi/lib -lpthread # paths to MPI libs and pthreads library (used when MPITHREADS == yes)
ZOLTANINC = -I/Users/tomek/Devel/lib/zoltan/include
ZOLTANLIB = -L/Users/tomek/Devel/lib/zoltan/lib -lzoltan
