ifeq ($(DEBUG),yes)
  CFLAGS = -Wall -g
else
  CFLAGS = -O -I.
endif

cc = $(CC)
AR = ar -rvcu
RANLIB = ranlib
AR_LIB_SUFFIX = a
