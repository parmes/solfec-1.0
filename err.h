/*
 * err.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * basic error handling macros
 */

/* This file is part of Solfec.
 * Solfec is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Solfec is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Solfec. If not, see <http://www.gnu.org/licenses/>. */

#include <setjmp.h>
#include <stdio.h>

#ifndef __err__
#define __err__

typedef struct errstack ERRSTACK;

struct errstack /* error stack */
{
  jmp_buf env;

  ERRSTACK *next;
};

extern ERRSTACK *__errstack__; /* global error stack */

/* error codes */
enum
{
  ERR_OUT_OF_MEMORY = 1,
  ERR_FILE_OPEN,
  ERR_FILE_EMPTY,
  ERR_FILE_FORMAT,
  ERR_NOT_IMPLEMENTED,
  ERR_TMS_INTEGRATE_CONSTANT,
  ERR_PBF_INDEX_FILE_CORRUPTED,
  ERR_PBF_OUTPUT_TIME_DECREASED,
  ERR_MTX_LU_FACTOR,
  ERR_MTX_MATRIX_INVERT,
  ERR_MTX_EIGEN_DECOMP,
  ERR_MTX_KIND,
  ERR_MTX_IFAC,
  ERR_MSH_UNSUPPORTED_ELEMENT,
  ERR_MSH_UNSUPPORTED_FACE,
  ERR_CVX_ZERO_NORMAL,
  ERR_LDY_EIGEN_DECOMP,
  ERR_GLV_WINDOW,
  ERR_BOD_NEW3_SINGULAR_JACOBIAN,
  ERR_BOD_NEW3_NEWTON_DIVERGENCE,
  ERR_BOD_MAX_FREQ_LE0,
  ERR_BOD_KIND,
  ERR_BOD_SCHEME,
  ERR_BOD_MOTION_INVERT,
  ERR_RND_NO_DOMAIN,
  ERR_RND_STACK_OVERFLOW,
  ERR_PCK_UNPACK,
  ERR_DOM_DEPTH,
  ERR_MPI_THREAD_MULTIPLE,
  ERR_ZOLTAN_INIT,
  ERR_ZOLTAN,
  ERR_GAUSS_SEIDEL_DIAGONAL_DIVERGED,
  ERR_GAUSS_SEIDEL_DIVERGED
};

/* get error string */
char* errstring (int error);

/* begin error handling */
#define TRY()\
  {\
    int __code__;\
    ERRSTACK __context__;\
    __context__.next = __errstack__;\
    __errstack__ = &__context__;\
    if (!(__code__ = setjmp (__context__.env)))\
    {\

/* catch an error (read argument) */
#define CATCH(__error__)\
    }\
    else if (__code__ == (__error__))\
    {\
      __errstack__ = __errstack__->next;\

/* catch all errors (write argument) */
#define CATCHANY(__error__)\
    }\
    else if (((__error__) = __code__))\
    {\
      __errstack__ = __errstack__->next;\

/* end error handling */
#define ENDTRY()\
    }\
    else\
    {\
      __errstack__ = __errstack__->next;\
      if (__errstack__) longjmp (__errstack__->env, __code__);\
      else\
      {\
	fprintf (stderr, "%s: %d, Unhandled error: %s\n", __FILE__, __LINE__, errstring (__code__));\
	exit (1);\
      }\
    }\
  }

/* throw an error */
#if NOTHROW
#define THROW(__error__)\
  do\
  {\
    fprintf (stderr, "%s: %d, Unhandled error: %s\n", __FILE__, __LINE__, errstring (__error__));\
    exit (1);\
  }\
  while (0)
#else
#define THROW(__error__)\
  do\
  {\
    if (__errstack__) longjmp (__errstack__->env, __error__);\
    else\
    {\
      fprintf (stderr, "%s: %d, Unhandled error: %s\n", __FILE__, __LINE__, errstring (__error__));\
      exit (1);\
    }\
  }\
  while (0)
#endif

/* general assertion */
#define ASSERT(__test__, __error__) if (! (__test__)) THROW (__error__)

/* memory validity assertion */
#define ERRMEM(__pointer__) if (! (__pointer__)) THROW (ERR_OUT_OF_MEMORY)

/* warning executed only in debug mode */
#if DEBUG
  #define WARNING_DEBUG(Test, ...)\
    if (! (Test)) { fprintf (stderr, "%s: %d => ", __FILE__, __LINE__);\
      fprintf (stderr, __VA_ARGS__);\
      fprintf (stderr, "\n"); }
#else
  #define WARNING_DEBUG(Test, ...)
#endif

/* assertion executed only in debug mode */
#if DEBUG
  #define ASSERT_DEBUG(Test, ...)\
    do {\
    if (! (Test)) { fprintf (stderr, "%s: %d => ", __FILE__, __LINE__);\
      fprintf (stderr, __VA_ARGS__);\
      fprintf (stderr, "\n"); exit (1); } } while (0)
#else
  #define ASSERT_DEBUG(Test, ...)
#endif

/* assertion executed only in debug mode,
 * but still performing the test in reslease mode */
#if DEBUG
  #define ASSERT_DEBUG_EXT(Test, ...)\
    do {\
    if (! (Test)) { fprintf (stderr, "%s: %d => ", __FILE__, __LINE__);\
      fprintf (stderr, __VA_ARGS__);\
      fprintf (stderr, "\n"); exit (1); } } while (0)
#else
  #define ASSERT_DEBUG_EXT(Test, ...) Test
#endif

/* auxiliary numeric file dumps => used for debug plotting */
#if DEBUG
  #define AUXDUMP1(Path, X)\
  {\
    static char *__o__ = "w";\
    FILE *__f__;\
    ASSERT_DEBUG (__f__ = fopen (Path, __o__), "Could not open file");\
    fprintf (__f__, "%.15e\n", X);\
    fclose (__f__);\
    __o__ = "a";\
  }
  #define AUXDUMP2(Path, X, Y)\
  {\
    static char *__o__ = "w";\
    FILE *__f__;\
    ASSERT_DEBUG (__f__ = fopen (Path, __o__), "Could not open file");\
    fprintf (__f__, "%.15e\t%.15e\n", X, Y);\
    fclose (__f__);\
    __o__ = "a";\
  }
  #define AUXDUMP3(Path, X, Y, Z)\
  {\
    static char *__o__ = "w";\
    FILE *__f__;\
    ASSERT_DEBUG (__f__ = fopen (Path, __o__), "Could not open file");\
    fprintf (__f__, "%.15e\t%.15e\t%.15e\n", X, Y, Z);\
    fclose (__f__);\
    __o__ = "a";\
  }
#else
  #define AUXDUMP1(Path, X)
  #define AUXDUMP2(Path, X, Y)
  #define AUXDUMP3(Path, X, Y, Z)
#endif

#endif
