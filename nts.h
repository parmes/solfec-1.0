/*
 * nts.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
 */

#include "ldy.h"
#include "mem.h"

#ifndef __nts__
#define __nts__

typedef struct newton NEWTON;

enum ntvariant
{
  NT_NONSMOOTH_HSW          =  0x01,
  NT_NONSMOOTH_HYBRID       =  0x02,
  NT_FIXED_POINT            =  0x04,
  NT_NONSMOOTH_VARIATIONAL  =  0x08,
  NT_SMOOTHED_VARIATIONAL   =  0x10,
  NT_SYMMETRIZE             =  0x20   /* symmetrize linear equations if possible */
};

typedef enum ntvariant NTVARIANT;

struct newton
{
  NTVARIANT variant; /* method variant */

  double epsilon; /* relative accuracy */

  int maxiter; /* iterations bound */

  MEM mapmem; /* memory of DIAB->map and OFFB->map int [3] vectors */

  int length; /* nonmonotone merit function buffer length */
};

/* create solver */
NEWTON* NEWTON_Create (NTVARIANT variant, double epsilon, int maxiter);

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy);

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf);

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt);

#endif
