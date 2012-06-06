/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.38 $
 ***********************************************************************EHEADER*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "bicgstab.h"
#include "cgnr.h"
#include "flexgmres.h"
#include "gmres.h"
#include "lgmres.h"
#include "pcg.h"

#ifndef hypre_KRYLOV_HEADER
#define hypre_KRYLOV_HEADER

extern int hypre__global_error;
#define hypre_error_flag            hypre__global_error
#define HYPRE_ERROR_GENERIC         1   /* generic error */
#define HYPRE_ERROR_MEMORY          2   /* unable to allocate memory */
#define HYPRE_ERROR_ARG             4   /* argument error */
/* bits 4-8 are reserved for the index of the argument error */
#define HYPRE_ERROR_CONV          256   /* method did not converge as expected */

void hypre_error_handler(char *filename, int line, int ierr);
#define hypre_error(IERR)  hypre_error_handler(__FILE__, __LINE__, IERR)


#define hypre_CTAllocF(type, count, funcs) ( (type *)(*(funcs->CAlloc))((unsigned int)(count), (unsigned int)sizeof(type)) )
#define hypre_TFreeF( ptr, funcs ) ( (*(funcs->Free))((char *)ptr), ptr = NULL )

#ifndef hypre_max
#define hypre_max(a,b)  (((a)<(b)) ? (b) : (a))
#endif
#ifndef hypre_min
#define hypre_min(a,b)  (((a)<(b)) ? (a) : (b))
#endif

#endif
