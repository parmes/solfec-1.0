/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
 */

#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "bgs.h"
#include "nts.h"
#include "alg.h"
#include "mem.h"
#include "err.h"

/* create solver */
NEWTON* NEWTON_Create (NTVARIANT variant, double epsilon, int maxiter)
{
  return NULL;
}

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy)
{
}

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf)
{
}

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt)
{
}
