#
# Operating System (WIN32, SOLARIS, LINUX, AIX, IRIX, OSX, FREEBSD)
#

OS = OSX

#
# Specify C compiler
#

CC = cc

#
# Specify C++ compiler
#

CXX = g++

#
# Specify FORTRAN95 compiler and FORTRAN runtime library
#

FC = gfortran
FCLIB = -lgfortran

# 
# Debug or optimized version switch (yes/no)
#

# general debug tests
DEBUG = yes
# gprof profiling data
PROFILE = no
# switch off memory pools
MEMDEBUG = no
# geometrical debug tests
GEOMDEBUG = no
# light parallel consitency tests
PARDEBUG = yes
# heavy parallel self-consitency tests
PSCTEST = yes
# disable TRY/CATCH
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
# HDF5
#

HDF5 = no

HDF5INC = -I/usr/local/hdf5/include
HDF5LIB = -L/usr/local/hdf5/lib -lhdf5 -lhdf5_hl

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

PYTHON = -I/usr/include/python2.6
PYTHONLIB = -L/usr/lib -lpython2.6

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

ZOLTANINC = -I/usr/local/include
ZOLTANLIB = -L/usr/local/lib -lzoltan

#
# Siconos (yes/no)
#

SICONOS = no
SICONOSINC = -I/usr/local/include/Siconos/Numerics
SICONOSLIB = -L/usr/local/lib -l SiconosNumerics
