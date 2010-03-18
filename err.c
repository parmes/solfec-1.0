/*
 * err.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * basic error handling
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

#include "err.h"

ERRSTACK *__errstack__ = NULL; /* global error context */

/* get error string */
char* errstring (int error)
{
  switch (error)
  {
    case ERR_OUT_OF_MEMORY: return "Out of memory";
    case ERR_FILE_OPEN: return "File open failed";
    case ERR_FILE_EMPTY: return "File empty";
    case ERR_FILE_FORMAT: return "Invalid file format";
    case ERR_NOT_IMPLEMENTED: return "Not implemented";
    case ERR_TMS_INTEGRATE_CONSTANT: return "Cannot integrate constant time series (no limits)";
    case ERR_PBF_INDEX_FILE_CORRUPTED: return "PBF index file corrupted";
    case ERR_PBF_OUTPUT_TIME_DECREASED: return "PBF output time decreased";
    case ERR_MTX_LU_FACTOR: return  "LU factorisation failed";
    case ERR_MTX_MATRIX_INVERT: return  "Matrix inversion failed";
    case ERR_MTX_EIGEN_DECOMP: return  "Eigen decomposition failed";
    case ERR_MTX_KIND: return "Invalid matrix kind";
    case ERR_MTX_IFAC: return "Invalid operation involving a factorised sparse inverse matrix";
    case ERR_MSH_UNSUPPORTED_ELEMENT: return "An element has an unsupported type";
    case ERR_MSH_UNSUPPORTED_FACE: return "An element face has an unsupported type";
    case ERR_CVX_ZERO_NORMAL: return "Zero face normal generated during a CONVEX face update";
    case ERR_LDY_EIGEN_DECOMP: return  "Eigen decomposition of a diagonal W block failed";
    case ERR_GLV_WINDOW: return "Creation of a viewer window failed";
    case ERR_BOD_NEW3_SINGULAR_JACOBIAN: return "NEW3 Singular Jacobian";
    case ERR_BOD_NEW3_NEWTON_DIVERGENCE: return "NEW3 Newton method divergence";
    case ERR_BOD_MAX_FREQ_LE0: return "Maximal eigen frequency of a body not positive";
    case ERR_BOD_KIND: return "Invalid body kind";
    case ERR_BOD_SCHEME: return "Invalid body time integration scheme";
    case ERR_BOD_SCHEME_NOT_CONVERGED: return "Time integration scheme failed to converge";
    case ERR_BOD_MOTION_INVERT: return "Body motion not invertible";
    case ERR_BOD_ENERGY_CONSERVATION: return "Energy conservation failed";
    case ERR_RND_NO_DOMAIN: return "Nothing to visualise";
    case ERR_RND_STACK_OVERFLOW: return "Stack overflow in rendering module";
    case ERR_PCK_UNPACK: return "Trying to unpack more elements then stored";
    case ERR_DOM_DEPTH: return "Unphysical interpenetration has occurred";
    case ERR_DOM_TOO_MANY_BODIES: return "Too many bodies";
    case ERR_DOM_TOO_MANY_CONSTRAINTS: return "Too many constraints";
    case ERR_ZOLTAN_INIT: return "Zoltan initialization failed";
    case ERR_ZOLTAN: return "Zoltan call failed";
    case ERR_GAUSS_SEIDEL_DIAGONAL_DIVERGED: return "GAUSS_SEIDEL_SOLVER failed with error code DIAGONAL_DIVERGED";
    case ERR_GAUSS_SEIDEL_DIVERGED: return "GAUSS_SEIDEL_SOLVER failed with error code DIVERGED";
    case ERR_FEM_MASS_NOT_SPD: return "Mass matrix is not symmetric-positive-definite in FEM module";
    case ERR_FEM_COORDS_INVERT: return "Global to local coordinates convertion failed in FEM module";
    case ERR_FEM_CUT_VOLUME: return "Volume cut of FE mesh does not match the shape volume in FEM module";
    case ERR_UMFPACK_SYMBOLIC: return "UMFPACK symbolic factorization has failed";
    case ERR_UMFPACK_NUMERIC: return "UMFPACK numeric factorization has failed";
    case ERR_UMFPACK_SOLVE: return "UMFPACK solution has failed";
  }

  return "Unknown";
}
