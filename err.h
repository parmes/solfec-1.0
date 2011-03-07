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

#if MPI
#include <mpi.h>
#endif

#ifndef __err__
#define __err__

typedef struct errstack ERRSTACK;

struct errstack /* error stack */
{
  jmp_buf env;

  ERRSTACK *next;
};

extern ERRSTACK *__errstack__; /* global error stack */

extern short WARNINGS_ENABLED; /* warnings flag */

/* error codes */
enum
{
  ERR_OUT_OF_MEMORY = 1,
  ERR_FILE_OPEN,
  ERR_FILE_READ,
  ERR_FILE_WRITE,
  ERR_FILE_CLOSE,
  ERR_FILE_EMPTY,
  ERR_FILE_FORMAT,
  ERR_NOT_IMPLEMENTED,
  ERR_TMS_INTEGRATE_CONSTANT,
  ERR_TMS_TIME_NOT_INCREASED, /* 10 */
  ERR_PBF_INDEX_FILE_CORRUPTED,
  ERR_PBF_OUTPUT_TIME_DECREASED,
  ERR_PBF_WRITE,
  ERR_PBF_READ,
  ERR_MTX_LU_FACTOR,
  ERR_MTX_CHOL_FACTOR,
  ERR_MTX_CHOL_SOLVE,
  ERR_MTX_MATRIX_INVERT,
  ERR_MTX_EIGEN_DECOMP,
  ERR_MTX_KIND, /* 20 */
  ERR_MTX_IFAC,
  ERR_MSH_UNSUPPORTED_ELEMENT,
  ERR_MSH_UNSUPPORTED_FACE,
  ERR_MSH_ELEMENT_WITH_NODE,
  ERR_CVX_ZERO_NORMAL,
  ERR_LDY_EIGEN_DECOMP,
  ERR_GLV_WINDOW,
  ERR_BOD_NEW3_SINGULAR_JACOBIAN,
  ERR_BOD_NEW3_NEWTON_DIVERGENCE,
  ERR_BOD_MAX_FREQ_LE0, /* 30 */
  ERR_BOD_KIND,
  ERR_BOD_SCHEME,
  ERR_BOD_SCHEME_NOT_CONVERGED,
  ERR_BOD_MOTION_INVERT,
  ERR_BOD_ENERGY_CONSERVATION,
  ERR_RND_NO_DOMAIN,
  ERR_RND_STACK_OVERFLOW,
  ERR_PCK_UNPACK,
  ERR_DOM_DEPTH,
  ERR_DOM_TOO_MANY_BODIES, /* 40 */
  ERR_DOM_TOO_MANY_CONSTRAINTS,
  ERR_ZOLTAN_INIT,
  ERR_ZOLTAN,
  ERR_GAUSS_SEIDEL_DIAGONAL_DIVERGED,
  ERR_GAUSS_SEIDEL_DIVERGED,
  ERR_FEM_MASS_NOT_SPD,
  ERR_FEM_COORDS_INVERT,
  ERR_FEM_CUT_VOLUME,
  ERR_FEM_ROT_SINGULAR_JACOBIAN,
  ERR_FEM_ROT_NEWTON_DIVERGENCE, /* 50 */
  ERR_BUG_FOUND
};

/* get error string */
char* errstring (int error);

/* exit function */
#if MPI
#define EXIT(code) MPI_Abort (MPI_COMM_WORLD, code)
#else
#define EXIT(code) exit (code)
#endif

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
	EXIT (1);\
      }\
    }\
  }

/* throw an error */
#if NOTHROW
#define THROW(__error__)\
  do\
  {\
    fprintf (stderr, "%s: %d, Unhandled error: %s\n", __FILE__, __LINE__, errstring (__error__));\
    EXIT (1);\
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
      EXIT (1);\
    }\
  }\
  while (0)
#endif

/* general assertion */
#define ASSERT(__test__, __error__) if (! (__test__)) THROW (__error__)

/* textual assertion */
#define ASSERT_TEXT(Test, ...)\
  do {\
  if (! (Test)) { fprintf (stderr, "%s: %d => ", __FILE__, __LINE__);\
    fprintf (stderr, __VA_ARGS__);\
    fprintf (stderr, "\n"); THROW (ERR_BUG_FOUND); } } while (0)

/* memory validity assertion */
#define ERRMEM(__pointer__) do { if (! (__pointer__)) THROW (ERR_OUT_OF_MEMORY); } while (0)

/* general warning */
#define WARNING(Test, ...)\
  if (! (Test) && WARNINGS_ENABLED) { fprintf (stderr, "WARNING %s: %d => ", __FILE__, __LINE__);\
    fprintf (stderr, __VA_ARGS__);\
    fprintf (stderr, "\n"); }

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
      fprintf (stderr, "\n"); EXIT (1); } } while (0)
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
      fprintf (stderr, "\n"); EXIT (1); } } while (0)
#else
  #define ASSERT_DEBUG_EXT(Test, ...) Test
#endif

/* auxiliary file dump */
#if DEBUG
  #if MPI
  #define FOPEN(Path, File)\
  {\
    static char *__o__ = "w";\
    char __path__ [256];\
    int __rank__;\
    MPI_Comm_rank (MPI_COMM_WORLD, &__rank__);\
    snprintf (__path__, 256, "%s%d", Path, __rank__);\
    ASSERT_DEBUG (File = fopen (__path__, __o__), "File open failed");\
    __o__ = "a";\
  }
  #else
  #define FOPEN(Path, File)\
  {\
    static char *__o__ = "w";\
    ASSERT_DEBUG (File = fopen (Path, __o__), "File open failed");\
    __o__ = "a";\
  }
  #endif
  #define AUXDUMP(Path, ...)\
  {\
    FILE *__f__;\
    FOPEN (Path, __f__);\
    fprintf (__f__, __VA_ARGS__);\
    fclose (__f__);\
  }
#else
  #define AUXDUMP(Path, ...)
#endif

#endif
