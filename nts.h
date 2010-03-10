/*
 * nts.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
 */

#include "ldy.h"

#ifndef __nts__
#define __nts__

typedef struct pb NEWTON;

enum ntvariant
{
  NT_NONSMOOTH_HSW,
  NT_NONSMOOTH_HYBRID,
  NT_NONSMOOTH_VARIATIONAL,
  NT_SMOOTHED_VARIATIONAL,
  NT_FIXED_POINT
};

typedef enum ntvariant NTVARIANT;

struct pb
{
  NTVARIANT variant; /* method variant */

  double epsilon; /* relative accuracy */

  int maxiter; /* iterations bound */
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
